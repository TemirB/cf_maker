#!/usr/bin/env python3
"""Durable, credential-free computation queue. Python 3.9+, POSIX; no web server."""
from __future__ import annotations

import argparse
import contextlib
import datetime as dt
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import secrets
import shutil
import signal
import subprocess
import sys
import tarfile
import tempfile
import time


INTERRUPTED = False


def utc_now():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="microseconds").replace("+00:00", "Z")


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_json(path, value):
    path = Path(path)
    descriptor, name = tempfile.mkstemp(prefix="." + path.name, dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as output:
            json.dump(value, output, indent=2, ensure_ascii=False, allow_nan=False)
            output.write("\n")
            output.flush()
            os.fsync(output.fileno())
        os.replace(name, path)
        directory = os.open(str(path.parent), os.O_RDONLY)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
    finally:
        if os.path.exists(name):
            os.unlink(name)


def load_config(path):
    path = Path(path).resolve()
    if path.stat().st_mode & 0o022:
        raise ValueError("lab config must not be writable by group/others")
    config = json.loads(path.read_text())
    if config.get("schema_version") != 1:
        raise ValueError("unsupported lab config schema")
    state = Path(config["state_dir"])
    if not state.is_absolute() or "REPLACE_WITH" in str(state):
        raise ValueError("set state_dir to a private absolute path outside the checkout")
    config["state_dir"] = str(state.resolve())
    files = config.get("config_files", [])
    if not files or len({Path(p).stem for p in files}) != len(files):
        raise ValueError("config_files must be nonempty with distinct filenames")
    for name in files:
        if Path(name).is_absolute() or ".." in Path(name).parts:
            raise ValueError("config_files must stay inside the source snapshot")
    for name, default in (("jobs", 2), ("threads", 2), ("max_recovery_attempts", 2),
                          ("upload_lease_seconds", 3600), ("upload_retry_seconds", 900),
                          ("artifact_log_limit_bytes", 10485760)):
        config[name] = int(config.get(name, default))
        if config[name] <= 0:
            raise ValueError(name + " must be positive")
    state.mkdir(parents=True, exist_ok=True, mode=0o700)
    if state.stat().st_mode & 0o077:
        raise ValueError("state_dir must have permissions 0700")
    for name in ("jobs", "queue"):
        (state / name).mkdir(exist_ok=True, mode=0o700)
    return config


def job_dir(config, job_id):
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,119}", job_id):
        raise ValueError("invalid job ID")
    return Path(config["state_dir"]) / "jobs" / job_id


@contextlib.contextmanager
def state_lock(config):
    with (Path(config["state_dir"]) / "state.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        yield


def read_job(config, job_id):
    path = job_dir(config, job_id) / "job.json"
    if not path.exists():
        # The durable queue marker is written first: recover a power loss between
        # that write and job.json without losing an accepted event.
        marker = Path(config["state_dir"]) / "queue" / job_id
        pending = json.loads(marker.read_text())
        if pending.get("job_id") != job_id or not pending.get("commit_sha"):
            raise ValueError("incomplete job record")
        atomic_json(path, pending)
    return json.loads(path.read_text())


def update_job(config, job_id, **changes):
    with state_lock(config):
        job = read_job(config, job_id)
        job.update(changes, updated_at=utc_now())
        atomic_json(job_dir(config, job_id) / "job.json", job)
    return job


def git(source, *args):
    return subprocess.check_output(["git", "-C", str(source), *args], text=True).strip()


def snapshot_repository(source, destination, expected_sha):
    """Export tracked code, recursively including checked-out submodules, without .git."""
    if git(source, "rev-parse", "HEAD") != expected_sha:
        raise ValueError("checkout HEAD does not match the requested SHA")

    def export(repository, target):
        target.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryFile() as archive:
            subprocess.run(["git", "-C", str(repository), "archive", "--format=tar", "HEAD"],
                           stdout=archive, check=True)
            archive.seek(0)
            with tarfile.open(fileobj=archive) as tar:
                for member in tar.getmembers():
                    candidate = (target / member.name).resolve()
                    if candidate != target.resolve() and target.resolve() not in candidate.parents:
                        raise ValueError("unsafe archive path")
                    if member.issym() or member.islnk():
                        link = (candidate.parent / member.linkname).resolve()
                        if link != target.resolve() and target.resolve() not in link.parents:
                            raise ValueError("archive symlink leaves snapshot")
                    if member.isdev() or member.isfifo():
                        raise ValueError("unsupported archive member")
                tar.extractall(target)
        modules = repository / ".gitmodules"
        if modules.is_file():
            result = subprocess.run(["git", "config", "--file", str(modules), "-z", "--get-regexp",
                                     r"^submodule\..*\.path$"], capture_output=True, check=False)
            if result.returncode not in (0, 1):
                raise ValueError("cannot read submodule paths")
            for record in result.stdout.split(b"\0"):
                if not record:
                    continue
                relative = Path(record.split(b"\n", 1)[1].decode())
                if relative.is_absolute() or ".." in relative.parts:
                    raise ValueError("submodule leaves source tree")
                child = repository / relative
                expected = git(repository, "ls-tree", "HEAD", "--", str(relative)).split()[2]
                if not child.is_dir() or git(child, "rev-parse", "HEAD") != expected:
                    raise ValueError("submodule missing or wrong commit: " + str(relative))
                export(child, target / relative)

    export(Path(source).resolve(), Path(destination))


def submit(config, source, job_id, commit_sha, run_id, repository, event_name="workflow_dispatch"):
    if not re.fullmatch(r"[0-9a-f]{40}", commit_sha):
        raise ValueError("commit SHA must be 40 lowercase hexadecimal characters")
    if not re.fullmatch(r"[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+", repository):
        raise ValueError("invalid repository")
    if not re.fullmatch(r"[0-9]+", run_id):
        raise ValueError("run ID must be numeric")
    source = Path(source).resolve()
    state = Path(config["state_dir"])
    if state == source or source in state.parents:
        raise ValueError("state_dir must be outside the checkout")
    directory = job_dir(config, job_id)
    profile_digest = hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest()
    with state_lock(config):
        if (directory / "job.json").exists():
            job = read_job(config, job_id)
            if (job["commit_sha"], job["run_id"], job["repository"]) != (commit_sha, run_id, repository):
                raise ValueError("job ID already belongs to a different event")
            print("Duplicate event: existing job", job_id, job["status"])
            if job["status"] in ("queued", "building", "running"):
                atomic_json(state / "queue" / job_id, job)
            return job
        if event_name == "push":
            for previous in all_jobs(config):
                if (previous.get("event_name") == "push" and previous["repository"] == repository and
                        previous["commit_sha"] == commit_sha and previous.get("profile_digest") == profile_digest):
                    print("Duplicate push: existing job", previous["job_id"], previous["status"])
                    return previous
        directory.mkdir(exist_ok=True, mode=0o700)
        if not (directory / "snapshot.json").exists():
            partial = directory / "source.partial"
            if partial.exists():
                shutil.rmtree(partial)
            snapshot_repository(source, partial, commit_sha)
            if (directory / "source").exists():
                shutil.rmtree(directory / "source")
            partial.rename(directory / "source")
            atomic_json(directory / "snapshot.json", {"commit_sha": commit_sha})
        elif json.loads((directory / "snapshot.json").read_text())["commit_sha"] != commit_sha:
            raise ValueError("interrupted submission has a different snapshot SHA")
        now = utc_now()
        job = dict(schema_version=1, job_id=job_id, commit_sha=commit_sha, run_id=run_id,
                   repository=repository, submitted_at=now, updated_at=now, status="queued",
                   attempts=0, cancel_requested=False, stages=[], exit_code=None,
                   event_name=event_name, profile_digest=profile_digest)
        atomic_json(state / "queue" / job_id, job)
        atomic_json(directory / "job.json", job)
    print("Queued", job_id, commit_sha)
    return job


class StageFailure(Exception):
    def __init__(self, stage, code, reason):
        super().__init__(reason)
        self.stage, self.code = stage, code


class WorkerInterrupted(Exception):
    pass


def clean_environment():
    allowed = {"PATH", "HOME", "USER", "LOGNAME", "SHELL", "LANG", "TMPDIR", "TMP", "TEMP",
               "ROOTSYS", "LD_LIBRARY_PATH", "DYLD_LIBRARY_PATH", "LIBRARY_PATH", "CPATH",
               "CMAKE_PREFIX_PATH", "PKG_CONFIG_PATH", "XDG_CACHE_HOME", "OMP_NUM_THREADS",
               "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"}
    return {key: value for key, value in os.environ.items() if key in allowed or key.startswith("LC_")}


def stop_process(process):
    if process.poll() is not None:
        return
    try:
        os.killpg(process.pid, signal.SIGTERM)
        process.wait(timeout=15)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        process.wait()
    except ProcessLookupError:
        process.wait()


def stage_command(config, name, default, values):
    command = config.get("commands", {}).get(name, default)
    if not isinstance(command, list) or not command or not all(isinstance(arg, str) for arg in command):
        raise ValueError("stage command must be a nonempty argument array")
    return [arg.format_map(values) for arg in command]


def run_stage(config, job_id, attempt, name, phase, command):
    job = read_job(config, job_id)
    if job.get("cancel_requested"):
        raise StageFailure(name, 130, "cancelled by operator")
    if INTERRUPTED:
        raise WorkerInterrupted()
    logs = attempt / "logs"
    logs.mkdir(exist_ok=True)
    prefix = re.sub(r"[^A-Za-z0-9_.-]", "_", name)
    stdout, stderr = logs / (prefix + ".stdout.log"), logs / (prefix + ".stderr.log")
    record = dict(stage=name, status=phase, started_at=utc_now(), command=command,
                  stdout=str(stdout), stderr=str(stderr), exit_code=None)
    update_job(config, job_id, status=phase, current_stage=name, stages=job["stages"] + [record])
    start = time.monotonic()
    wrapped = command
    if config.get("setup_script"):
        wrapped = ["bash", "--noprofile", "--norc", "-c", 'set -e; source "$1"; shift; exec "$@"',
                   "cf-maker-env", config["setup_script"], *command]
    process = None
    code = 127
    interrupted = False
    try:
        with stdout.open("wb") as out, stderr.open("wb") as err:
            process = subprocess.Popen(wrapped, cwd=job_dir(config, job_id) / "source",
                                       env=clean_environment(), stdout=out, stderr=err,
                                       start_new_session=True)
            while True:
                try:
                    code = process.wait(timeout=1)
                    break
                except subprocess.TimeoutExpired:
                    if INTERRUPTED:
                        interrupted = True
                        stop_process(process)
                        code = 143
                        break
                    if read_job(config, job_id).get("cancel_requested"):
                        stop_process(process)
                        code = 130
                        break
    except OSError as error:
        with stderr.open("a") as err:
            err.write(str(error) + "\n")
    finally:
        if process is not None and process.poll() is None:
            stop_process(process)
        record.update(completed_at=utc_now(), exit_code=code,
                      elapsed_seconds=time.monotonic() - start)
        job = read_job(config, job_id)
        update_job(config, job_id, stages=job["stages"][:-1] + [record])
    if interrupted:
        raise WorkerInterrupted()
    if code != 0:
        raise StageFailure(name, code, "stage failed: " + name)


def check_manifest(results, not_before):
    path = Path(results) / "run_manifest.json"
    data = json.loads(path.read_text())
    fits = data.get("fits", {})
    started = dt.datetime.fromisoformat(data.get("started_at", "").replace("Z", "+00:00"))
    earliest = dt.datetime.fromisoformat(not_before.replace("Z", "+00:00"))
    if started.tzinfo is None or earliest.tzinfo is None:
        raise ValueError("manifest timestamps must contain a timezone")
    if (data.get("status") != "completed" or data.get("exit_code") != 0 or
            not isinstance(data.get("run_id"), str) or not data["run_id"] or
            started < earliest - dt.timedelta(milliseconds=1) or not data.get("completed_at") or
            fits.get("requested", 0) <= 0 or fits.get("usable") != fits.get("requested") or
            fits.get("unusable") != 0):
        raise ValueError("analysis manifest does not describe a fresh, fully successful computation")
    declared = data["config"]["vars"]["output"]["dir"]
    if Path(declared).resolve() != Path(results).resolve():
        raise ValueError("manifest refers to another output directory")
    return data


def check_pdf(path, not_before_ns=0):
    path = Path(path)
    stat = path.stat()
    with path.open("rb") as source:
        if source.read(5) != b"%PDF-":
            raise ValueError("not a PDF: " + str(path))
    if stat.st_size < 20 or stat.st_mtime_ns < not_before_ns:
        raise ValueError("empty or stale PDF: " + str(path))
    return path


def execute_job(config, job_id):
    job = read_job(config, job_id)
    directory = job_dir(config, job_id)
    source = directory / "source"
    number = job["attempts"] + 1
    attempt = directory / ("attempt-%04d" % number)
    attempt.mkdir(mode=0o700)
    start_at, start_ns = utc_now(), time.time_ns()
    update_job(config, job_id, attempts=number, attempt_dir=str(attempt),
               started_at=job.get("started_at", start_at), attempt_started_at=start_at,
               status="building", error=None, exit_code=None)
    tools = config.get("tools", {})
    values = dict(source=str(source), attempt=str(attempt), build=str(attempt / "build"),
                  job_id=job_id, jobs=str(config["jobs"]), threads=str(config["threads"]),
                  cxx=config.get("cxx", "g++"), make=tools.get("make", "make"),
                  root=tools.get("root", "root"), latexmk=tools.get("latexmk", "latexmk"),
                  gs=tools.get("gs", "gs"))

    def run(key, label, phase, default):
        command = stage_command(config, key, default, values)
        run_stage(config, job_id, attempt, label, phase, command)

    run("preflight", "report-tools", "building", ["{make}", "check-presentation-tools",
                                                  "BUILD_DIR={build}", "LATEXMK={latexmk}"])
    run("build", "build", "building", ["{make}", "build", "BUILD_DIR={build}",
                                         "JOBS={jobs}", "CXX={cxx}"])
    pdfs, manifests = [], []
    for name in config["config_files"]:
        original = source / name
        effective = json.loads(original.read_text())
        if not effective.get("stages", {}).get("cf3d", True):
            raise ValueError("automatic report requires a new cf3d computation, not resume")
        if config.get("input_root"):
            effective["machine"]["base_input"] = config["input_root"]
        label = original.stem
        results = attempt / "results" / label
        effective["machine"]["base_output"] = str(attempt / "results")
        effective["vars"]["output"]["dir"] = label
        effective["threads"] = config["threads"]
        effective["images"] = {"need": True, "format": "pdf"}
        generated = attempt / (label + ".json")
        atomic_json(generated, effective)
        values.update(config=str(generated), results=str(results), label=label,
                      docs=str(attempt / "docs" / label), audit_pdf=str(results / "chi2_audit.pdf"))
        analysis_started = utc_now()
        run("compute", label + "-compute", "running", ["{build}/main", "{config}"])
        manifest = check_manifest(results, analysis_started)
        manifests.append(dict(config=name, run_id=manifest["run_id"],
                              manifest_sha256=sha256(results / "run_manifest.json")))
        run("report", label + "-presentation", "running", [
            "{make}", "presentation", "CONFIG={config}", "RESULTS={results}", "BUILD_DIR={build}",
            "DOCS_BUILD_DIR={docs}", "PRESENTATION_PLOTS=0", "LATEXMK={latexmk}"])
        # The macro produces all contribution plots and native ROOT/manual diagnostics.
        expression = ("tools/manual_chi2.C(" + json.dumps(str(results / "cf3d.root")) +
                      ",-1,-1,-1,\"\"," + json.dumps(str(results / "chi2_audit.pdf")) + ")")
        values["audit_expression"] = expression
        run("audit", label + "-chi2-audit", "running", ["{root}", "-l", "-b", "-q", "{audit_expression}"])
        if sha256(results / "run_manifest.json") != manifests[-1]["manifest_sha256"]:
            raise ValueError("scientific results changed during report creation")
        pdfs.extend([check_pdf(results / "presentation.pdf", start_ns),
                     check_pdf(results / "chi2_audit.pdf", start_ns)])
    values["report_pdf"] = str(attempt / "report.pdf")
    default_merge = ["{gs}", "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=pdfwrite",
                     "-sOutputFile={report_pdf}", *[str(path) for path in pdfs]]
    run("merge", "merge-report", "running", default_merge)
    pdf = check_pdf(values["report_pdf"], start_ns)
    if read_job(config, job_id).get("cancel_requested"):
        raise StageFailure("seal", 130, "cancelled before publication")
    update_job(config, job_id, status="uploading", current_stage="awaiting-publication",
               completed_at=utc_now(), exit_code=0, computation_exit_code=0,
               pdf_path=str(pdf), pdf_sha256=sha256(pdf), pdf_size=pdf.stat().st_size,
               scientific_manifests=manifests, lease_token=None, lease_until=0,
               next_upload_at=0, elapsed_seconds=(time.time_ns() - start_ns) / 1e9)


def worker(config):
    global INTERRUPTED
    INTERRUPTED = False

    def interrupt(_signal, _frame):
        global INTERRUPTED
        INTERRUPTED = True

    previous = {sig: signal.signal(sig, interrupt) for sig in (signal.SIGINT, signal.SIGTERM)}
    try:
        with (Path(config["state_dir"]) / "worker.lock").open("a") as lock:
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                return
            while not INTERRUPTED:
                with state_lock(config):
                    queued = sorted((Path(config["state_dir"]) / "queue").iterdir(),
                                    key=lambda p: (read_job(config, p.name)["submitted_at"], p.name))
                if not queued:
                    break
                marker = queued[0]
                job_id = marker.name
                job = read_job(config, job_id)
                if job["status"] in ("success", "failed", "uploading"):
                    marker.unlink()
                    continue
                try:
                    if job.get("cancel_requested"):
                        raise StageFailure("queued", 130, "cancelled by operator")
                    if job["attempts"] >= config["max_recovery_attempts"]:
                        raise StageFailure("recovery", 143, "interrupted; recovery attempt limit reached")
                    execute_job(config, job_id)
                except WorkerInterrupted:
                    update_job(config, job_id, status="queued", error="worker interrupted; next start uses a fresh attempt")
                    break
                except StageFailure as error:
                    update_job(config, job_id, status="failed", exit_code=error.code,
                               failed_stage=error.stage, error=str(error), completed_at=utc_now())
                except Exception as error:
                    update_job(config, job_id, status="failed", exit_code=1,
                               failed_stage="validation", error=str(error), completed_at=utc_now())
                marker.unlink(missing_ok=True)
    finally:
        for sig, handler in previous.items():
            signal.signal(sig, handler)


def cancel(config, job_id):
    with state_lock(config):
        job = read_job(config, job_id)
        if job["status"] == "uploading" and job.get("lease_token") and job.get("lease_until", 0) > time.time():
            raise ValueError("publication is already leased: cancel the publisher Actions run as well")
        if job["status"] == "success":
            raise ValueError("job already published; its historical report is retained")
        job.update(cancel_requested=True, updated_at=utc_now())
        if job["status"] in ("queued", "uploading", "failed"):
            job.update(status="failed", exit_code=130, failed_stage="cancelled",
                       error="cancelled by operator", completed_at=utc_now(), next_upload_at=None)
        atomic_json(job_dir(config, job_id) / "job.json", job)
    print("Cancellation requested:", job_id)


def retry_upload(config, job_id):
    with state_lock(config):
        job = read_job(config, job_id)
        if job.get("computation_exit_code") != 0 or not job.get("pdf_path"):
            raise ValueError("job has no sealed successful report; dispatch a new computation")
        if job["status"] == "success":
            print("Already published:", job_id)
            return
        if job.get("lease_token") and job.get("lease_until", 0) > time.time():
            raise ValueError("publication already in progress")
        if sha256(check_pdf(job["pdf_path"])) != job["pdf_sha256"]:
            raise ValueError("sealed PDF checksum changed")
        job.update(status="uploading", exit_code=0, cancel_requested=False,
                   next_upload_at=0, lease_token=None, lease_until=0, error=None, updated_at=utc_now())
        atomic_json(job_dir(config, job_id) / "job.json", job)


def write_github_outputs(path, values):
    if path:
        with Path(path).open("a") as out:
            for key, value in values.items():
                if "\n" in str(value) or "\r" in str(value):
                    raise ValueError("invalid Actions output")
                out.write(key + "=" + str(value) + "\n")


def all_jobs(config):
    return sorted((json.loads(path.read_text()) for path in
                   (Path(config["state_dir"]) / "jobs").glob("*/job.json")),
                  key=lambda j: (j["submitted_at"], j["job_id"]))


def status(config, summary_path=None, job_id=None, full_logs=False):
    lines = ["| Job | Commit | Status | Attempt | Stage | Exit |", "|---|---|---|---:|---|---:|"]
    jobs = [read_job(config, job_id)] if job_id else all_jobs(config)
    for job in jobs[-30:]:
        lines.append("| {job_id} | {commit} | {status} | {attempts} | {stage} | {exit} |".format(
            **job, commit=job["commit_sha"][:12], stage=job.get("current_stage", ""),
            exit=job.get("exit_code", "")))
    text = "\n".join(lines) + "\n"
    print(text)
    if summary_path:
        with Path(summary_path).open("a") as out:
            out.write(text)
    for job in jobs[-5:]:
        if job_id or job["status"] in ("building", "running", "failed"):
            print(job["job_id"], job.get("error") or "")
            for stage in (job["stages"] if full_logs else job["stages"][-1:]):
                for field in ("stdout", "stderr"):
                    path = Path(stage[field])
                    if path.is_file():
                        with path.open("rb") as source:
                            if not full_logs:
                                source.seek(max(0, path.stat().st_size - 8192))
                            print(stage["stage"] + " " + field + ":")
                            guard = secrets.token_hex(16)
                            if os.environ.get("GITHUB_ACTIONS") == "true":
                                print("::stop-commands::" + guard)
                            for chunk in iter(lambda: source.read(65536), b""):
                                print(chunk.decode("utf-8", "replace"), end="")
                            print()
                            if os.environ.get("GITHUB_ACTIONS") == "true":
                                print("::" + guard + "::")


def collect(config, stage_dir, github_output=None, github_summary=None):
    directory = Path(stage_dir)
    if directory.exists() and any(directory.iterdir()):
        raise ValueError("delivery stage directory must be empty")
    selected = None
    now = time.time()
    with state_lock(config):
        for job in all_jobs(config):
            retryable = job["status"] == "failed" and job.get("failed_stage") == "upload"
            if (job.get("cancel_requested") or not (job["status"] == "uploading" or retryable) or
                    job.get("next_upload_at", 0) > now or job.get("lease_until", 0) > now):
                continue
            if job.get("computation_exit_code") != 0:
                continue
            job.update(status="uploading", exit_code=0, lease_token=secrets.token_hex(16),
                       lease_until=now + config["upload_lease_seconds"], updated_at=utc_now())
            atomic_json(job_dir(config, job["job_id"]) / "job.json", job)
            selected = job
            break
    if not selected:
        write_github_outputs(github_output, {"ready": "false"})
        status(config, github_summary)
        return None
    try:
        pdf = check_pdf(selected["pdf_path"])
        if sha256(pdf) != selected["pdf_sha256"] or pdf.stat().st_size != selected["pdf_size"]:
            raise ValueError("sealed PDF changed; publication refused")
        directory.mkdir(parents=True, exist_ok=True, mode=0o700)
        shutil.copyfile(pdf, directory / "report.pdf")
        manifest = {key: selected[key] for key in ("schema_version", "job_id", "commit_sha", "run_id",
                    "repository", "submitted_at", "started_at", "completed_at", "pdf_sha256", "pdf_size")}
        manifest.update(status="uploading", exit_code=0, pdf_path="report.pdf",
                        attempt=selected["attempts"], scientific_manifests=selected["scientific_manifests"])
        atomic_json(directory / "manifest.json", manifest)
        logs = directory / "logs"
        logs.mkdir()
        atomic_json(logs / "stages.json", selected["stages"])
        limit = config["artifact_log_limit_bytes"]
        for stage in selected["stages"]:
            for field in ("stdout", "stderr"):
                source = Path(stage[field])
                if source.is_file():
                    name = "%02d_" % selected["stages"].index(stage) + source.name
                    with source.open("rb") as original, (logs / name).open("wb") as target:
                        if source.stat().st_size > limit:
                            target.write(b"[Artifact contains log tail; complete log retained on lab PC]\n")
                            original.seek(-limit, os.SEEK_END)
                        shutil.copyfileobj(original, target)
    except Exception:
        acknowledge(config, selected["job_id"], selected["lease_token"], "failed", "delivery staging failed")
        raise
    values = dict(ready="true", job_id=selected["job_id"], lease_token=selected["lease_token"],
                  pdf_sha256=selected["pdf_sha256"], commit_sha=selected["commit_sha"],
                  artifact_name="cf-maker-" + selected["job_id"])
    write_github_outputs(github_output, values)
    status(config, github_summary)
    return values


def acknowledge(config, job_id, lease_token, result, error="", latest_url=""):
    with state_lock(config):
        job = read_job(config, job_id)
        if job.get("last_ack_token") == lease_token:
            return
        if not lease_token or job.get("lease_token") != lease_token:
            raise ValueError("stale publication receipt")
        job.update(last_ack_token=lease_token, lease_token=None, lease_until=0, updated_at=utc_now())
        if result == "success":
            job.update(status="success", exit_code=0, published_at=utc_now(), error=None,
                       latest_url=latest_url, next_upload_at=None)
        else:
            failures = job.get("upload_failures", 0) + 1
            delay = min(3600, config["upload_retry_seconds"] * 2 ** min(failures - 1, 8))
            job.update(status="failed", exit_code=1, failed_stage="upload", error=error,
                       upload_failures=failures, next_upload_at=time.time() + delay)
        atomic_json(job_dir(config, job_id) / "job.json", job)


def main():
    os.umask(0o077)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    commands = parser.add_subparsers(dest="operation", required=True)
    enqueue = commands.add_parser("submit")
    for name in ("source", "job-id", "sha", "run-id", "repository"):
        enqueue.add_argument("--" + name, required=True)
    enqueue.add_argument("--event-name", choices=("push", "workflow_dispatch"), default="workflow_dispatch")
    commands.add_parser("worker")
    show = commands.add_parser("status")
    show.add_argument("--github-summary")
    show.add_argument("--job-id")
    show.add_argument("--full-logs", action="store_true")
    for name in ("cancel", "retry-upload"):
        command = commands.add_parser(name)
        command.add_argument("--job-id", required=True)
    delivery = commands.add_parser("collect")
    delivery.add_argument("--stage-dir", required=True)
    delivery.add_argument("--github-output")
    delivery.add_argument("--github-summary")
    receipt = commands.add_parser("ack")
    receipt.add_argument("--job-id", required=True)
    receipt.add_argument("--lease-token", required=True)
    receipt.add_argument("--result", choices=("success", "failed"), required=True)
    receipt.add_argument("--error", default="")
    receipt.add_argument("--latest-url", default="")
    args = parser.parse_args()
    try:
        config = load_config(args.config)
        if args.operation == "submit":
            submit(config, args.source, args.job_id, args.sha, args.run_id, args.repository, args.event_name)
        elif args.operation == "worker":
            worker(config)
        elif args.operation == "status":
            status(config, args.github_summary, args.job_id, args.full_logs)
        elif args.operation == "cancel":
            cancel(config, args.job_id)
        elif args.operation == "retry-upload":
            retry_upload(config, args.job_id)
        elif args.operation == "collect":
            collect(config, args.stage_dir, args.github_output, args.github_summary)
        elif args.operation == "ack":
            acknowledge(config, args.job_id, args.lease_token, args.result, args.error, args.latest_url)
    except Exception as error:
        print("lab automation:", error, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
