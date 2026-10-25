"""Offline queue/receipt tests: python3 -B -m unittest discover -s tests -p test_lab_queue.py."""

import contextlib
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
import unittest
from unittest import mock


SCRIPT = Path(__file__).resolve().parents[1] / "tools/lab_automation/lab.py"
SPEC = importlib.util.spec_from_file_location("lab_queue", SCRIPT)
lab = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(lab)

FAKE_STAGE = r'''
import datetime as dt
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
import uuid

stage, attempt, mode, *extra = sys.argv[1:]
attempt = Path(attempt)
print(stage, attempt.name, flush=True)
if stage == "build" and mode == "exit":
    sys.exit(7)
if stage == "compute":
    with (attempt.parent / "compute.calls").open("a") as calls:
        calls.write(attempt.name + "\n")
    if mode in ("block", "block_once") and (mode == "block" or attempt.name == "attempt-0001"):
        child = subprocess.Popen([sys.executable, "-B", "-c", "import time; time.sleep(90)"])
        def stop(_signal, _frame):
            child.wait(timeout=5)
            sys.exit(143)
        signal.signal(signal.SIGTERM, stop)
        marker = attempt.parent / "running.json"
        temporary = marker.with_suffix(".tmp")
        temporary.write_text(json.dumps({"parent": os.getpid(), "child": child.pid}))
        temporary.replace(marker)
        while True:
            time.sleep(.02)
    if mode == "exit":
        sys.exit(2)
    results = Path(extra[0])
    results.mkdir(parents=True)
    now = dt.datetime.now(dt.timezone.utc).isoformat(timespec="milliseconds").replace("+00:00", "Z")
    manifest = {"status": "completed", "exit_code": 0, "run_id": uuid.uuid4().hex,
                "started_at": "2000-01-01T00:00:00.000Z" if mode == "stale" else now,
                "completed_at": now, "fits": {"requested": 1, "usable": 1, "unusable": 0},
                "config": {"vars": {"output": {"dir": str(results.resolve())}}}}
    (results / "run_manifest.json").write_text(json.dumps(manifest))
    (results / "cf3d.root").write_bytes(b"fake ROOT input for the fake audit")
elif stage in ("report", "audit", "merge"):
    path = Path(extra[0])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"invalid PDF body" if mode == "wrong" else b"%PDF-1.7\n" + stage.encode() + b" body\n%%EOF\n")
'''


class QueueTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="cf-maker-queue-test-")
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)
        self.source = self.base / "checkout"
        self.source.mkdir()
        self.git(self.source, "init", "-q")
        self.git(self.source, "config", "user.name", "Offline test")
        self.git(self.source, "config", "user.email", "test@example.invalid")
        (self.source / "config").mkdir()
        (self.source / "config/one.json").write_text(json.dumps({
            "machine": {"base_input": str(self.base), "base_output": str(self.base)},
            "vars": {"input": {"file": "input.root", "type": "kt"}, "output": {"dir": "unused"}}}))
        (self.source / "fake_stage.py").write_text(FAKE_STAGE)
        self.commit()
        self.path = self.base / "lab.json"
        settings = {"schema_version": 1, "state_dir": str(self.base / "private-state"),
                    "config_files": ["config/one.json"], "upload_lease_seconds": 30,
                    "upload_retry_seconds": 60, "commands": {}}
        extras = {"compute": "{results}", "report": "{results}/presentation.pdf",
                  "audit": "{audit_pdf}", "merge": "{report_pdf}"}
        for stage in ("preflight", "build", "compute", "report", "audit", "merge"):
            command = [sys.executable, "-B", "{source}/fake_stage.py", stage, "{attempt}", "ok"]
            settings["commands"][stage] = command + ([extras[stage]] if stage in extras else [])
        self.path.write_text(json.dumps(settings))
        self.path.chmod(0o600)
        self.config = lab.load_config(self.path)

    def git(self, directory, *args):
        return subprocess.check_output(["git", "-C", str(directory), *args], text=True,
                                       stderr=subprocess.PIPE).strip()

    def commit(self):
        self.git(self.source, "add", ".")
        self.git(self.source, "commit", "-qm", "offline fixture")
        self.sha = self.git(self.source, "rev-parse", "HEAD")

    def quiet(self, function, *args, **kwargs):
        with contextlib.redirect_stdout(io.StringIO()):
            return function(*args, **kwargs)

    def enqueue(self, name="job1"):
        return self.quiet(lab.submit, self.config, self.source, name, self.sha, "100", "owner/repo")

    def run_job(self, name="job1"):
        self.enqueue(name)
        self.quiet(lab.worker, self.config)
        return lab.read_job(self.config, name)

    def mode(self, stage, value):
        self.config["commands"][stage][5] = value
        self.path.write_text(json.dumps(self.config))

    def collect(self, name):
        return self.quiet(lab.collect, self.config, self.base / name)

    def cli_command(self, *args):
        return [sys.executable, "-B", str(SCRIPT), "--config", str(self.path), *args]

    def cli(self, *args):
        return subprocess.run(self.cli_command(*args), capture_output=True, text=True, timeout=15)

    def wait_for(self, predicate, timeout=8):
        end = time.monotonic() + timeout
        while time.monotonic() < end:
            if predicate():
                return
            time.sleep(.02)
        self.fail("Timed out waiting for the worker/process state")

    def stop_worker(self, process):
        if process.poll() is None:
            process.terminate()
        process.wait(timeout=12)

    def blocked_worker(self, mode="block"):
        self.mode("compute", mode)
        self.enqueue()
        process = subprocess.Popen(self.cli_command("worker"), stdout=subprocess.DEVNULL,
                                   stderr=subprocess.DEVNULL)
        self.addCleanup(self.stop_worker, process)
        path = lab.job_dir(self.config, "job1") / "running.json"
        self.wait_for(path.exists)
        return process, json.loads(path.read_text())

    @staticmethod
    def alive(pid):
        try:
            os.kill(pid, 0)
            return True
        except ProcessLookupError:
            return False

    def test_completed_report_receipt_and_duplicate_event_do_not_recompute(self):
        job = self.run_job()
        self.assertEqual((job["status"], job["attempts"], job["computation_exit_code"]), ("uploading", 1, 0))
        receipt = self.collect("delivery")
        bundle = self.base / "delivery"
        manifest = json.loads((bundle / "manifest.json").read_text())
        self.assertEqual(manifest["pdf_sha256"], hashlib.sha256((bundle / "report.pdf").read_bytes()).hexdigest())
        self.assertEqual(manifest["scientific_manifests"], job["scientific_manifests"])
        scientific = Path(job["attempt_dir"]) / "results/one/run_manifest.json"
        self.assertEqual(manifest["scientific_manifests"][0]["manifest_sha256"],
                         hashlib.sha256(scientific.read_bytes()).hexdigest())
        lab.acknowledge(self.config, "job1", receipt["lease_token"], "success", latest_url="https://example.invalid/latest")
        self.assertEqual(lab.read_job(self.config, "job1")["status"], "success")
        self.enqueue()
        self.quiet(lab.worker, self.config)
        self.assertEqual((lab.job_dir(self.config, "job1") / "compute.calls").read_text().splitlines(), ["attempt-0001"])

    def test_stage_failures_and_invalid_outputs_never_become_publishable(self):
        for index, (stage, mode, code) in enumerate((("build", "exit", 7), ("compute", "exit", 2),
                                                    ("compute", "stale", 1), ("report", "wrong", 1))):
            with self.subTest(stage=stage, mode=mode):
                self.mode(stage, mode)
                job = self.run_job("failure" + str(index))
                self.assertEqual((job["status"], job["exit_code"]), ("failed", code))
                self.assertIsNone(self.collect("refused" + str(index)))
                self.assertNotIn("pdf_sha256", job)
                self.mode(stage, "ok")

    def test_duplicate_push_delivery_ids_reuse_job_but_manual_run_is_explicit(self):
        first = self.quiet(lab.submit, self.config, self.source, "push-one", self.sha,
                           "100", "owner/repo", "push")
        repeated = self.quiet(lab.submit, self.config, self.source, "push-two", self.sha,
                              "101", "owner/repo", "push")
        self.assertEqual(first["job_id"], repeated["job_id"])
        self.assertFalse(lab.job_dir(self.config, "push-two").exists())
        self.quiet(lab.submit, self.config, self.source, "manual", self.sha,
                   "102", "owner/repo", "workflow_dispatch")
        self.assertEqual(len(lab.all_jobs(self.config)), 2)
        self.quiet(lab.worker, self.config)
        self.assertTrue(all(job["status"] == "uploading" for job in lab.all_jobs(self.config)))

    def test_durable_marker_recovers_power_loss_before_job_record(self):
        self.enqueue()
        (lab.job_dir(self.config, "job1") / "job.json").unlink()
        self.quiet(lab.worker, self.config)
        job = lab.read_job(self.config, "job1")
        self.assertEqual((job["status"], job["attempts"]), ("uploading", 1))

    def test_child_environment_does_not_inherit_workflow_or_developer_secrets(self):
        with mock.patch.dict(os.environ, {"GITHUB_TOKEN": "secret", "ACTIONS_RUNTIME_TOKEN": "secret",
                                          "CF_DRIVE_REFRESH_TOKEN": "secret", "OTHER_API_KEY": "secret",
                                          "ROOTSYS": "/test/root"}):
            cleaned = lab.clean_environment()
        self.assertNotIn("GITHUB_TOKEN", cleaned)
        self.assertNotIn("ACTIONS_RUNTIME_TOKEN", cleaned)
        self.assertNotIn("CF_DRIVE_REFRESH_TOKEN", cleaned)
        self.assertNotIn("OTHER_API_KEY", cleaned)
        self.assertEqual(cleaned["ROOTSYS"], "/test/root")

    def test_queued_and_running_cancellation_stop_the_whole_process_group(self):
        self.enqueue("queued")
        self.quiet(lab.cancel, self.config, "queued")
        self.quiet(lab.worker, self.config)
        self.assertEqual(lab.read_job(self.config, "queued")["exit_code"], 130)
        self.assertFalse((lab.job_dir(self.config, "queued") / "compute.calls").exists())
        process, pids = self.blocked_worker()
        result = self.cli("cancel", "--job-id", "job1")
        self.assertEqual(result.returncode, 0, result.stderr)
        process.wait(timeout=12)
        self.wait_for(lambda: not any(self.alive(pid) for pid in pids.values()))
        job = lab.read_job(self.config, "job1")
        self.assertEqual((job["status"], job["exit_code"]), ("failed", 130))
        self.assertIsNone(self.collect("cancelled"))

    def test_interrupted_worker_restarts_in_new_attempt_and_retains_previous(self):
        process, pids = self.blocked_worker("block_once")
        process.send_signal(signal.SIGTERM)
        process.wait(timeout=12)
        self.wait_for(lambda: not any(self.alive(pid) for pid in pids.values()))
        first = lab.read_job(self.config, "job1")
        self.assertEqual((first["status"], first["attempts"]), ("queued", 1))
        old_log = Path(first["attempt_dir"]) / "logs/one-compute.stdout.log"
        original = old_log.read_bytes()
        self.quiet(lab.worker, self.config)
        final = lab.read_job(self.config, "job1")
        self.assertEqual((final["status"], final["attempts"]), ("uploading", 2))
        self.assertNotEqual(first["attempt_dir"], final["attempt_dir"])
        self.assertEqual(old_log.read_bytes(), original)
        self.assertEqual((old_log.parent.parent.parent / "compute.calls").read_text().splitlines(),
                         ["attempt-0001", "attempt-0002"])

    def test_concurrent_duplicate_events_export_one_snapshot_with_local_submodule(self):
        child = self.base / "local-module"
        child.mkdir()
        for arguments in (("init", "-q"), ("config", "user.name", "Offline test"),
                          ("config", "user.email", "test@example.invalid")):
            self.git(child, *arguments)
        (child / "tracked.txt").write_text("submodule payload\n")
        self.git(child, "add", ".")
        self.git(child, "commit", "-qm", "module fixture")
        self.git(self.source, "-c", "protocol.file.allow=always", "submodule", "add", "-q", str(child), "vendor/module")
        self.commit()
        (self.source / "untracked-secret.txt").write_text("must stay outside the snapshot")
        command = self.cli_command("submit", "--source", str(self.source), "--job-id", "same-event",
                                   "--sha", self.sha, "--run-id", "100", "--repository", "owner/repo")
        processes = [subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE) for _ in range(2)]
        for process in processes:
            _, stderr = process.communicate(timeout=15)
            self.assertEqual(process.returncode, 0, stderr.decode())
        directory = lab.job_dir(self.config, "same-event")
        self.assertEqual((directory / "source/vendor/module/tracked.txt").read_text(), "submodule payload\n")
        self.assertFalse((directory / "source/.git").exists())
        self.assertFalse((directory / "source/vendor/module/.git").exists())
        self.assertFalse((directory / "source/untracked-secret.txt").exists())
        self.assertFalse((directory / "source.partial").exists())
        self.assertEqual(len(lab.all_jobs(self.config)), 1)
        self.quiet(lab.worker, self.config)
        self.assertEqual((directory / "compute.calls").read_text().splitlines(), ["attempt-0001"])

    def test_lease_busy_expiry_and_stale_receipt(self):
        self.run_job()
        first = self.collect("first")
        self.assertIsNone(self.collect("busy"))
        expires = lab.read_job(self.config, "job1")["lease_until"]
        with mock.patch.object(lab.time, "time", return_value=expires + 1):
            second = self.collect("after-expiry")
        self.assertNotEqual(first["lease_token"], second["lease_token"])
        with self.assertRaisesRegex(ValueError, "stale"):
            lab.acknowledge(self.config, "job1", first["lease_token"], "success")
        lab.acknowledge(self.config, "job1", second["lease_token"], "success")
        lab.acknowledge(self.config, "job1", second["lease_token"], "success")
        self.assertEqual(lab.read_job(self.config, "job1")["status"], "success")

    def test_failed_upload_backs_off_and_can_retry_without_recomputation(self):
        original = self.run_job()
        receipt = self.collect("initial")
        lab.acknowledge(self.config, "job1", receipt["lease_token"], "failed", "simulated upload error")
        failed = lab.read_job(self.config, "job1")
        self.assertEqual((failed["status"], failed["failed_stage"]), ("failed", "upload"))
        self.assertGreater(failed["next_upload_at"], time.time())
        self.assertIsNone(self.collect("too-early"))
        with mock.patch.object(lab.time, "time", return_value=failed["next_upload_at"] + 1):
            automatic = self.collect("automatic-retry")
        lab.acknowledge(self.config, "job1", automatic["lease_token"], "failed", "second error")
        self.quiet(lab.retry_upload, self.config, "job1")
        manual = self.collect("manual-retry")
        lab.acknowledge(self.config, "job1", manual["lease_token"], "success")
        self.quiet(lab.worker, self.config)
        final = lab.read_job(self.config, "job1")
        self.assertEqual((final["attempts"], final["pdf_sha256"], final["pdf_path"]),
                         (1, original["pdf_sha256"], original["pdf_path"]))
        self.assertEqual((lab.job_dir(self.config, "job1") / "compute.calls").read_text().splitlines(), ["attempt-0001"])

    def test_changed_sealed_pdf_is_refused_by_collection_and_retry(self):
        job = self.run_job()
        Path(job["pdf_path"]).write_bytes(b"%PDF-1.7\nchanged sealed content\n%%EOF\n")
        with self.assertRaisesRegex(ValueError, "changed"):
            self.collect("tampered")
        self.assertFalse((self.base / "tampered/report.pdf").exists())
        with self.assertRaisesRegex(ValueError, "checksum"):
            self.quiet(lab.retry_upload, self.config, "job1")

    def test_cli_usage_and_private_state_validation(self):
        result = subprocess.run([sys.executable, "-B", str(SCRIPT), "worker"], capture_output=True)
        self.assertEqual(result.returncode, 2)
        self.assertEqual(self.cli("submit").returncode, 2)
        self.assertEqual(self.cli("ack", "--job-id", "x", "--lease-token", "x", "--result", "unknown").returncode, 2)
        self.assertEqual(Path(self.config["state_dir"]).stat().st_mode & 0o777, 0o700)
        Path(self.config["state_dir"]).chmod(0o755)
        with self.assertRaisesRegex(ValueError, "0700"):
            lab.load_config(self.path)


if __name__ == "__main__":
    unittest.main()
