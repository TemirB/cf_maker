"""Offline protocol/transaction tests for the hosted Drive publisher.

Run: python3 -m unittest discover -s tests -p test_lab_publish.py
"""

import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import tempfile
import unittest
import urllib.parse


SCRIPT = Path(__file__).resolve().parents[1] / "tools" / "lab_automation" / "publish.py"
SPEC = importlib.util.spec_from_file_location("lab_publish", SCRIPT)
publish = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(publish)


def reply(status=200, value=None, headers=None):
    return publish.Response(status, headers or {}, json.dumps(value or {}).encode("utf-8"))


class FakeDrive:
    """A small stateful fake of actual OAuth, Drive and resumable HTTP requests."""

    def __init__(self):
        self.files = {
            "root_folder": {"id": "root_folder", "name": "LabReports", "mimeType": publish.FOLDER_TYPE,
                            "parents": ["root"], "appProperties": {"cf_role": "root"}},
            "history_folder": {"id": "history_folder", "name": "history", "mimeType": publish.FOLDER_TYPE,
                               "parents": ["root_folder"],
                               "appProperties": {"cf_role": "history_folder"}},
        }
        self.sessions = {}
        self.calls = []
        self.next_id = 0
        self.fail_history = False
        self.drop_ack_role = None
        self.dropped_ack = False
        self.rate_limit_lists = 0
        self.scope = publish.SCOPE
        self.fragment_upload = False
        self.corrupt_checksum_role = None

    def __call__(self, method, url, headers, body):
        self.calls.append((method, url, copy.deepcopy(headers), body))
        parsed = urllib.parse.urlparse(url)
        params = urllib.parse.parse_qs(parsed.query)
        if url == publish.TOKEN_URL:
            return reply(value={"access_token": "DO_NOT_LOG_ACCESS_TOKEN", "scope": self.scope,
                                "expires_in": 3600})
        if headers.get("Authorization") != "Bearer DO_NOT_LOG_ACCESS_TOKEN":
            return reply(401)
        if parsed.path.endswith("/generateIds"):
            self.next_id += 1
            return reply(value={"ids": ["file_{}".format(self.next_id)]})
        if parsed.path == "/drive/v3/files" and method == "GET":
            if self.rate_limit_lists:
                self.rate_limit_lists -= 1
                return reply(429, headers={"Retry-After": "3"})
            query = params["q"][0]
            result = list(self.files.values())
            parent = re.search(r"'([^']+)' in parents", query)
            if parent:
                result = [item for item in result if parent.group(1) in item.get("parents", [])]
            name = re.search(r"name = '([^']+)'", query)
            if name:
                result = [item for item in result if item.get("name") == name.group(1)]
            for key, value in re.findall(r"key='([^']+)' and value='([^']+)'", query):
                result = [item for item in result if item.get("appProperties", {}).get(key) == value]
            mime = re.search(r"mimeType = '([^']+)'", query)
            if mime:
                result = [item for item in result if item.get("mimeType") == mime.group(1)]
            return reply(value={"files": [copy.deepcopy(item) for item in result if not item.get("trashed")]})
        if parsed.path.startswith("/drive/v3/files/") and method == "GET":
            file_id = parsed.path.rsplit("/", 1)[1]
            return reply(value=copy.deepcopy(self.files[file_id])) if file_id in self.files else reply(404)
        if parsed.path == "/drive/v3/files" and method == "POST":
            value = json.loads(body)
            if value["id"] in self.files:
                return reply(409)
            value.setdefault("parents", ["root"])
            self.files[value["id"]] = value
            return reply(value=copy.deepcopy(value))
        if parsed.path.startswith("/upload/drive/v3/files"):
            if "session" not in params:
                value = json.loads(body)
                file_id = value.get("id") if method == "POST" else parsed.path.rsplit("/", 1)[1]
                if self.fail_history and value.get("appProperties", {}).get("cf_role") == "history":
                    return reply(403, {"error": {"errors": [{"reason": "insufficientFilePermissions"}]}})
                if method == "POST" and file_id in self.files:
                    return reply(409)
                session_id = str(len(self.sessions) + 1)
                self.sessions[session_id] = {"id": file_id, "metadata": value, "data": b"",
                                             "total": int(headers["X-Upload-Content-Length"]),
                                             "complete": False}
                return reply(headers={"Location": publish.UPLOAD_API + "?session=" + session_id})
            session = self.sessions[params["session"][0]]
            content_range = headers.get("Content-Range", "")
            if content_range.startswith("bytes */"):
                if session["complete"]:
                    return reply(value=copy.deepcopy(self.files[session["id"]]))
                range_header = {"Range": "bytes=0-{}".format(len(session["data"]) - 1)} if session["data"] else {}
                return reply(308, headers=range_header)
            match = re.fullmatch(r"bytes (\d+)-(\d+)/(\d+)", content_range)
            if not match or int(match.group(1)) != len(session["data"]):
                return reply(400)
            if int(match.group(2)) - int(match.group(1)) + 1 != len(body):
                return reply(400)
            if self.fragment_upload and not session["data"]:
                # Simulate a transport break after the server received only a prefix.
                session["data"] = body[:128 * 1024]
                self.fragment_upload = False
                raise publish.TransportError("simulated network interruption")
            session["data"] += body
            if len(session["data"]) < session["total"]:
                return reply(308, headers={"Range": "bytes=0-{}".format(len(session["data"]) - 1)})
            old = self.files.get(session["id"], {})
            value = dict(old, **session["metadata"])
            value.update(id=session["id"], size=str(len(session["data"])),
                         md5Checksum=hashlib.md5(session["data"]).hexdigest(),
                         sha256Checksum=hashlib.sha256(session["data"]).hexdigest())
            if value.get("appProperties", {}).get("cf_role") == self.corrupt_checksum_role:
                value["md5Checksum"] = "0" * 32
            self.files[session["id"]] = value
            session["complete"] = True
            if not self.dropped_ack and value.get("appProperties", {}).get("cf_role") == self.drop_ack_role:
                self.dropped_ack = True
                raise publish.TransportError("simulated lost final acknowledgement")
            return reply(value=copy.deepcopy(value))
        raise AssertionError("Unexpected HTTP request: {} {}".format(method, parsed.path))

    def role_files(self, role):
        return [item for item in self.files.values() if item.get("appProperties", {}).get("cf_role") == role]

    def latest_starts(self):
        return [call for call in self.calls if call[0] in ("POST", "PATCH") and
                "/upload/" in call[1] and
                json.loads(call[3]).get("appProperties", {}).get("cf_role") == "latest"]


class PublishTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.directory = Path(self.temp.name)
        self.drive = FakeDrive()
        self.sleeps = []
        self.client = publish.DriveClient("CLIENT_ID", "CLIENT_SECRET", "REFRESH_TOKEN",
                                          transport=self.drive, sleep=self.sleeps.append,
                                          attempts=3, chunk_size=256 * 1024)

    def tearDown(self):
        self.temp.cleanup()

    def bundle(self, job="job1", run="10", submitted="2026-10-25T10:00:00Z", data=None):
        data = data or b"%PDF-1.7\nreport body\n%%EOF\n"
        directory = self.directory / job
        directory.mkdir(exist_ok=True)
        (directory / "report.pdf").write_bytes(data)
        manifest = {"schema_version": 1, "job_id": job, "repository": "owner/repo",
                    "commit_sha": "a" * 40, "run_id": run, "submitted_at": submitted,
                    "completed_at": "2026-10-25T20:00:00Z", "exit_code": 0,
                    "pdf_sha256": hashlib.sha256(data).hexdigest(), "pdf_size": len(data),
                    "pdf_path": "report.pdf", "status": "uploading"}
        (directory / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
        return publish.load_bundle(directory)

    def send(self, bundle):
        return publish.publish_bundle(self.client, bundle, "root_folder", "history_folder")

    def quality_bundle(self, quality="incomplete", job="job1", run="10",
                       submitted="2026-10-25T10:00:00Z"):
        bundle = self.bundle(job, run, submitted)
        incomplete = quality == "incomplete"
        bundle.manifest.update(
            report_quality=quality, computation_exit_code=2 if incomplete else 0,
            fit_quality_warnings=[{"stage": "kt-compute", "exit_code": 2,
                                   "error": "6/32 fits are unusable", "requested_fits": 32,
                                   "unusable_fits": 6}] if incomplete else [],
            scientific_manifests=[
                {"config": "config/kt.json", "run_id": "science-kt",
                 "status": "incomplete" if incomplete else "completed",
                 "exit_code": 2 if incomplete else 0},
                {"config": "config/y.json", "run_id": "science-y",
                 "status": "completed", "exit_code": 0}])
        (bundle.pdf.parent / "manifest.json").write_text(json.dumps(bundle.manifest), encoding="utf-8")
        return publish.load_bundle(bundle.pdf.parent)

    def test_explicit_incomplete_report_uploads_history_then_stable_latest_without_hiding_exit_2(self):
        first = self.send(self.bundle())
        bundle = self.quality_bundle(job="job2", run="11", submitted="2026-10-25T11:00:00Z")
        original = copy.deepcopy(bundle.manifest)
        result = self.send(bundle)
        self.assertEqual(result["status"], "published")
        self.assertEqual(result["latest_id"], first["latest_id"])
        self.assertEqual(len(self.drive.role_files("history")), 2)
        self.assertEqual(bundle.manifest, original)
        self.assertEqual(bundle.manifest["exit_code"], 0)
        self.assertEqual(bundle.manifest["computation_exit_code"], 2)
        self.assertEqual(bundle.manifest["scientific_manifests"][0]["exit_code"], 2)
        self.assertEqual(bundle.manifest["scientific_manifests"][0]["status"], "incomplete")
        history_get = next(i for i, call in enumerate(self.drive.calls)
                           if call[0] == "GET" and "/files/" + result["version_id"] + "?" in call[1])
        latest_start = self.drive.calls.index(self.drive.latest_starts()[-1])
        self.assertLess(history_get, latest_start)

    def test_explicit_complete_report_and_legacy_scientific_references_remain_publishable(self):
        complete = self.quality_bundle("complete")
        result = self.send(complete)
        self.assertEqual(result["status"], "published")
        self.assertEqual(complete.manifest["computation_exit_code"], 0)
        legacy = self.bundle("legacy", "11", "2026-10-25T11:00:00Z")
        legacy.manifest["scientific_manifests"] = [{"config": "config/kt.json", "run_id": "old-science",
                                                   "manifest_sha256": "b" * 64}]
        (legacy.pdf.parent / "manifest.json").write_text(json.dumps(legacy.manifest), encoding="utf-8")
        result = self.send(publish.load_bundle(legacy.pdf.parent))
        self.assertEqual(result["status"], "published")

    def test_incomplete_report_requires_explicit_consistent_quality_evidence(self):
        bundle = self.quality_bundle()
        original = copy.deepcopy(bundle.manifest)
        invalid = [
            {"report_quality": "failed"}, {"report_quality": None},
            {"report_quality": "complete"}, {"computation_exit_code": 0},
            {"computation_exit_code": 1}, {"computation_exit_code": 3},
            {"computation_exit_code": True}, {"fit_quality_warnings": []},
            {"fit_quality_warnings": None}, {"fit_quality_warnings": {}},
            {"fit_quality_warnings": [{"stage": "kt-compute", "exit_code": 1, "error": "fatal"}]},
            {"fit_quality_warnings": [{"stage": "kt-compute", "exit_code": True, "error": "fatal"}]},
            {"fit_quality_warnings": [{"stage": "", "exit_code": 2, "error": "unusable"}]},
            {"fit_quality_warnings": [{"stage": 2, "exit_code": 2, "error": "unusable"}]},
            {"fit_quality_warnings": [{"stage": "kt-compute", "exit_code": 2, "error": ""}]},
            {"fit_quality_warnings": [{"stage": "kt-compute", "exit_code": 2, "error": []}]},
            {"fit_quality_warnings": ["kt-compute"]},
            {"scientific_manifests": []}, {"scientific_manifests": None},
            {"scientific_manifests": [{"status": "completed", "exit_code": 0}]},
            {"scientific_manifests": [{"status": "failed", "exit_code": 1}]},
            {"scientific_manifests": [{"status": "incomplete", "exit_code": 0}]},
            {"scientific_manifests": [{"status": "completed", "exit_code": 2}]},
            {"scientific_manifests": [{"status": "completed", "exit_code": False},
                                      {"status": "incomplete", "exit_code": 2}]},
            {"scientific_manifests": [{"status": "incomplete", "exit_code": 2},
                                      {"config": "config/y.json"}]},
            {"exit_code": 2}, {"status": "failed"}, {"pdf_path": "diagnostic.pdf"},
        ]
        manifest_path = bundle.pdf.parent / "manifest.json"
        for change in invalid:
            with self.subTest(change=change):
                manifest_path.write_text(json.dumps(dict(original, **change)), encoding="utf-8")
                with self.assertRaises(publish.PublishError):
                    publish.load_bundle(bundle.pdf.parent)
        for key in ("report_quality", "computation_exit_code", "fit_quality_warnings", "scientific_manifests"):
            with self.subTest(missing=key):
                changed = copy.deepcopy(original)
                del changed[key]
                manifest_path.write_text(json.dumps(changed), encoding="utf-8")
                with self.assertRaises(publish.PublishError):
                    publish.load_bundle(bundle.pdf.parent)
        self.assertFalse(self.drive.calls)

    def test_complete_and_legacy_reports_cannot_conceal_incomplete_scientific_runs(self):
        bundle = self.quality_bundle("complete")
        original = copy.deepcopy(bundle.manifest)
        incomplete = [{"status": "incomplete", "exit_code": 2}]
        manifest_path = bundle.pdf.parent / "manifest.json"
        for change in ({"computation_exit_code": 2}, {"computation_exit_code": 7},
                       {"computation_exit_code": False}, {"fit_quality_warnings": [{}]},
                       {"scientific_manifests": incomplete}):
            with self.subTest(change=change):
                manifest_path.write_text(json.dumps(dict(original, **change)), encoding="utf-8")
                with self.assertRaises(publish.PublishError):
                    publish.load_bundle(bundle.pdf.parent)
        legacy = {key: value for key, value in original.items()
                  if key not in {"report_quality", "computation_exit_code", "fit_quality_warnings"}}
        legacy["scientific_manifests"] = incomplete
        manifest_path.write_text(json.dumps(legacy), encoding="utf-8")
        with self.assertRaises(publish.PublishError):
            publish.load_bundle(bundle.pdf.parent)
        self.assertFalse(self.drive.calls)

    def test_diagnostic_artifact_is_not_a_publishable_bundle(self):
        (self.directory / "diagnostic.pdf").write_bytes(b"%PDF-1.7\ndiagnostic only\n%%EOF\n")
        (self.directory / "diagnostic_manifest.json").write_text(json.dumps({
            "status": "failed", "exit_code": 2, "publishable": False,
            "pdf_path": "diagnostic.pdf"}))
        with self.assertRaisesRegex(publish.PublishError, "manifest"):
            publish.load_bundle(self.directory)
        self.assertFalse(self.drive.calls)

    def test_success_saves_history_before_latest_and_verifies_both(self):
        bundle = self.bundle()
        result = self.send(bundle)
        self.assertEqual(result["status"], "published")
        self.assertEqual(len(self.drive.role_files("history")), 1)
        self.assertEqual(len(self.drive.role_files("latest")), 1)
        self.assertNotEqual(result["version_id"], result["latest_id"])
        self.assertIn(result["latest_id"], result["latest_url"])
        self.assertEqual(self.drive.files[result["latest_id"]]["md5Checksum"], bundle.md5)
        history_get = next(i for i, call in enumerate(self.drive.calls)
                           if call[0] == "GET" and "/files/" + result["version_id"] + "?" in call[1])
        latest_start = self.drive.calls.index(self.drive.latest_starts()[0])
        self.assertLess(history_get, latest_start)

    def test_new_generation_updates_same_latest_id_and_keeps_history(self):
        first = self.send(self.bundle())
        second = self.send(self.bundle("job2", "11", "2026-10-25T11:00:00Z",
                                       b"%PDF-1.7\nnewer report\n%%EOF\n"))
        self.assertEqual(first["latest_id"], second["latest_id"])
        self.assertEqual(len(self.drive.role_files("history")), 2)
        self.assertEqual([call[0] for call in self.drive.latest_starts()], ["POST", "PATCH"])

    def test_failed_history_upload_never_requests_latest(self):
        first = self.send(self.bundle())
        previous = copy.deepcopy(self.drive.files[first["latest_id"]])
        previous_starts = len(self.drive.latest_starts())
        self.drive.fail_history = True
        with self.assertRaises(publish.PublishError):
            self.send(self.bundle("job2", "11", "2026-10-25T11:00:00Z"))
        self.assertEqual(previous, self.drive.files[first["latest_id"]])
        self.assertEqual(previous_starts, len(self.drive.latest_starts()))

    def test_lost_final_acknowledgement_resolves_session_without_false_failure(self):
        for role in ("history", "latest"):
            with self.subTest(role=role):
                self.drive = FakeDrive()
                self.drive.drop_ack_role = role
                self.client.transport = self.drive
                result = self.send(self.bundle())
                self.assertEqual(result["status"], "published")
                self.assertTrue(self.drive.dropped_ack)
                probes = [call for call in self.drive.calls
                          if call[2].get("Content-Range", "").startswith("bytes */")]
                self.assertEqual(len(probes), 1)
                self.assertEqual(len(self.drive.role_files("history")), 1)
                self.assertEqual(len(self.drive.role_files("latest")), 1)

    def test_retry_same_job_does_not_duplicate_archive_or_replace_latest(self):
        bundle = self.bundle()
        first = self.send(bundle)
        starts = len(self.drive.latest_starts())
        second = self.send(bundle)
        self.assertEqual(second["status"], "already_published")
        self.assertEqual(second["version_id"], first["version_id"])
        self.assertEqual(second["latest_id"], first["latest_id"])
        self.assertEqual(len(self.drive.role_files("history")), 1)
        self.assertEqual(len(self.drive.latest_starts()), starts)

    def test_stale_generation_is_archived_without_regressing_latest(self):
        first = self.send(self.bundle("newer", "11", "2026-10-25T11:00:00Z"))
        latest_before = copy.deepcopy(self.drive.files[first["latest_id"]])
        second = self.send(self.bundle("older", "10", "2026-10-25T10:00:00Z"))
        self.assertEqual(second["status"], "archived_stale")
        self.assertEqual(second["latest_id"], first["latest_id"])
        self.assertEqual(latest_before, self.drive.files[first["latest_id"]])
        self.assertEqual(len(self.drive.role_files("history")), 2)

    def test_same_timestamp_run_ids_sort_numerically(self):
        first = self.send(self.bundle("run9", "9"))
        second = self.send(self.bundle("run10", "10"))
        self.assertEqual(second["status"], "published")
        self.assertEqual(first["latest_id"], second["latest_id"])
        self.assertEqual(self.drive.files[second["latest_id"]]["appProperties"]["run_id"], "10")

    def test_interrupted_chunk_resumes_from_server_received_offset(self):
        self.drive.fragment_upload = True
        data = b"%PDF-1.7\n" + b"x" * (600 * 1024) + b"\n%%EOF\n"
        bundle = self.bundle(data=data)
        result = self.send(bundle)
        self.assertEqual(result["status"], "published")
        ranges = [call[2].get("Content-Range", "") for call in self.drive.calls]
        self.assertTrue(any(value.startswith("bytes 131072-") for value in ranges))
        publish.verify_pdf(self.drive.files[result["latest_id"]], bundle)

    def test_rate_limit_honors_retry_after(self):
        self.drive.rate_limit_lists = 1
        self.send(self.bundle())
        self.assertTrue(self.sleeps)
        self.assertGreaterEqual(self.sleeps[0], 3)

    def test_broad_scope_is_rejected_without_uploading(self):
        self.drive.scope = "https://www.googleapis.com/auth/drive"
        with self.assertRaisesRegex(publish.PublishError, "only drive.file"):
            self.send(self.bundle())
        self.assertFalse(self.drive.sessions)

    def test_changed_manifest_checksum_fails_before_any_network(self):
        bundle = self.bundle()
        bundle.pdf.write_bytes(b"%PDF-1.7\nmodified\n")
        with self.assertRaisesRegex(publish.PublishError, "checksum"):
            publish.load_bundle(bundle.pdf.parent)
        self.assertFalse(self.drive.calls)

    def test_failed_source_and_path_traversal_are_rejected(self):
        bundle = self.bundle()
        manifest = bundle.manifest
        for key, value in (("exit_code", 2), ("status", "failed"),
                           ("pdf_path", "../private.pdf"), ("pdf_size", True)):
            with self.subTest(key=key):
                changed = dict(manifest, **{key: value})
                (bundle.pdf.parent / "manifest.json").write_text(json.dumps(changed), encoding="utf-8")
                with self.assertRaises(publish.PublishError):
                    publish.load_bundle(bundle.pdf.parent)

    def test_conflicting_same_job_pdf_is_rejected_before_latest(self):
        self.send(self.bundle())
        starts = len(self.drive.latest_starts())
        different = self.bundle(data=b"%PDF-1.7\ndifferent job output\n%%EOF\n")
        with self.assertRaises(publish.PublishError):
            self.send(different)
        self.assertEqual(len(self.drive.role_files("history")), 1)
        self.assertEqual(starts, len(self.drive.latest_starts()))

    def test_init_reuses_own_folders_and_requires_no_broad_scope(self):
        first = publish.initialize(self.client)
        second = publish.initialize(self.client)
        self.assertEqual(first, second)
        self.assertEqual(first["root_folder_id"], "root_folder")
        self.assertEqual(first["history_folder_id"], "history_folder")
        self.assertEqual(first["scope_verification"], "refresh_response")
        self.assertFalse(self.drive.sessions)

    def test_init_creates_private_folders_with_metadata_requests_only(self):
        self.drive.files.clear()
        output = publish.initialize(self.client)
        root = self.drive.files[output["root_folder_id"]]
        history = self.drive.files[output["history_folder_id"]]
        self.assertEqual(root["name"], "LabReports")
        self.assertEqual(history["parents"], [root["id"]])
        self.assertTrue(all("permissions" not in item for item in self.drive.files.values()))
        self.assertEqual(self.drive.scope, publish.SCOPE)

    def test_duplicate_latest_filename_fails_instead_of_guessing(self):
        first = self.send(self.bundle())
        duplicate = copy.deepcopy(self.drive.files[first["latest_id"]])
        duplicate["id"] = "duplicate_latest"
        self.drive.files[duplicate["id"]] = duplicate
        starts = len(self.drive.latest_starts())
        with self.assertRaisesRegex(publish.PublishError, "Multiple latest.pdf"):
            self.send(self.bundle("job2", "11", "2026-10-25T11:00:00Z"))
        self.assertEqual(starts, len(self.drive.latest_starts()))

    def test_preallocated_id_conflict_is_verified_instead_of_duplicating(self):
        bundle = self.bundle()
        file_id = self.client.allocate_id()
        metadata = {"name": "archive.pdf", "mimeType": "application/pdf",
                    "parents": ["history_folder"], "appProperties": {"cf_role": "history"}}
        first = self.client.upload(bundle, metadata, file_id, create=True)
        second = self.client.upload(bundle, metadata, file_id, create=True)
        self.assertEqual(first["id"], second["id"])
        self.assertEqual(len(self.drive.role_files("history")), 1)
        self.assertEqual(len(self.drive.sessions), 1)

    def test_remote_archive_checksum_mismatch_never_requests_latest(self):
        self.drive.corrupt_checksum_role = "history"
        with self.assertRaisesRegex(publish.PublishError, "checksum"):
            self.send(self.bundle())
        self.assertFalse(self.drive.latest_starts())

    def test_untrusted_resumable_session_host_is_rejected(self):
        response = publish.Response(200, {"Location": "https://example.com/upload/drive/v3/files?session=stolen"}, b"")
        with self.assertRaisesRegex(publish.PublishError, "session URL"):
            publish.DriveClient.session_url(response)


if __name__ == "__main__":
    unittest.main()
