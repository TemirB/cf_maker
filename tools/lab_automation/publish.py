#!/usr/bin/env python3
"""Publish a verified lab report with a drive.file OAuth grant (Python 3.9+).

Run this program on the publisher, never on the laboratory worker.  The caller
must serialize publishers for a Drive root; for Actions use a shared concurrency
group with cancel-in-progress: false.  This program never changes permissions or
deletes files.  The history copy is verified before latest.pdf is touched.
"""

import argparse
import dataclasses
import datetime as dt
import email.utils
import hashlib
import json
import os
from pathlib import Path
import random
import re
import socket
import sys
import time
import urllib.error
import urllib.parse
import urllib.request


SCOPE = "https://www.googleapis.com/auth/drive.file"
API = "https://www.googleapis.com/drive/v3/files"
UPLOAD_API = "https://www.googleapis.com/upload/drive/v3/files"
TOKEN_URL = "https://oauth2.googleapis.com/token"
FOLDER_TYPE = "application/vnd.google-apps.folder"
FIELDS = "id,name,mimeType,parents,trashed,size,md5Checksum,sha256Checksum,appProperties"
ID_PATTERN = re.compile(r"[A-Za-z0-9_-]{1,200}\Z")
RETRYABLE = {429, 500, 502, 503, 504}


class PublishError(Exception):
    """An error whose message is safe to show in workflow logs."""


class TransportError(PublishError):
    pass


class ApiError(PublishError):
    def __init__(self, status, reason=""):
        self.status = status
        self.reason = reason
        super().__init__("Drive request failed (HTTP {}{})".format(
            status, "; " + reason if reason else ""))


@dataclasses.dataclass
class Response:
    status: int
    headers: dict
    body: bytes


def http_transport(method, url, headers, body):
    """Return HTTP errors as responses; never include URL/token in an exception."""
    request = urllib.request.Request(url, data=body, headers=headers, method=method)
    try:
        with urllib.request.urlopen(request, timeout=90) as reply:
            return Response(reply.status, dict(reply.headers), reply.read(2 * 1024 * 1024))
    except urllib.error.HTTPError as error:
        return Response(error.code, dict(error.headers), error.read(2 * 1024 * 1024))
    except (urllib.error.URLError, OSError, socket.timeout):
        raise TransportError("Network request failed; credentials and session URLs are omitted") from None


def decode_json(response):
    try:
        result = json.loads(response.body.decode("utf-8"))
        if not isinstance(result, dict):
            raise ValueError()
        return result
    except (ValueError, UnicodeError):
        raise PublishError("Server returned an invalid JSON response") from None


def header(response, name):
    return next((value for key, value in response.headers.items()
                 if key.lower() == name.lower()), None)


def error_reason(response):
    # Only these documented codes may enter logs, never a free-form server body.
    allowed = {"rateLimitExceeded", "userRateLimitExceeded", "storageQuotaExceeded",
               "insufficientFilePermissions", "notFound", "dailyLimitExceeded"}
    try:
        errors = json.loads(response.body).get("error", {}).get("errors", [])
        return next((item["reason"] for item in errors if item.get("reason") in allowed), "")
    except (ValueError, AttributeError, TypeError, KeyError):
        return ""


def valid_id(value):
    if not isinstance(value, str) or not ID_PATTERN.fullmatch(value):
        raise PublishError("Invalid Drive file/folder ID")
    return value


class DriveClient:
    def __init__(self, client_id, client_secret, refresh_token, transport=http_transport,
                 sleep=time.sleep, attempts=5, chunk_size=8 * 1024 * 1024):
        if not all((client_id, client_secret, refresh_token)):
            raise PublishError("OAuth client ID, client secret and refresh token are required")
        if chunk_size <= 0 or chunk_size % (256 * 1024):
            raise ValueError("Upload chunks must be multiples of 256 KiB")
        self.credentials = {"client_id": client_id, "client_secret": client_secret,
                            "refresh_token": refresh_token, "grant_type": "refresh_token"}
        self.transport = transport
        self.sleep = sleep
        self.attempts = attempts
        self.chunk_size = chunk_size
        self.access_token = None
        self.expires_at = 0
        self.scope_verified = False

    def pause(self, attempt, response=None):
        delay = min(30.0, 2.0 ** attempt) + random.uniform(0, 0.25)
        retry_after = header(response, "Retry-After") if response else None
        if retry_after:
            try:
                requested = float(retry_after)
            except ValueError:
                try:
                    requested = (email.utils.parsedate_to_datetime(retry_after) -
                                 dt.datetime.now(dt.timezone.utc)).total_seconds()
                except (ValueError, TypeError, OverflowError):
                    requested = 0
            # Stop instead of violating a long server-requested delay or hanging a job.
            if requested > 120:
                raise PublishError("Server requested a retry delay over 120 seconds; retry publication later")
            delay = max(delay, requested)
        self.sleep(max(0, delay))

    @staticmethod
    def retryable(response):
        return response.status in RETRYABLE or (
            response.status == 403 and error_reason(response) in
            {"rateLimitExceeded", "userRateLimitExceeded"})

    def refresh(self):
        body = urllib.parse.urlencode(self.credentials).encode("ascii")
        for attempt in range(self.attempts):
            response = None
            try:
                response = self.transport("POST", TOKEN_URL,
                                          {"Content-Type": "application/x-www-form-urlencoded"}, body)
                if response.status == 200:
                    token = decode_json(response)
                    if not isinstance(token.get("access_token"), str) or not token["access_token"]:
                        raise PublishError("OAuth refresh did not return an access token")
                    if "scope" in token:
                        if not isinstance(token["scope"], str) or set(token["scope"].split()) != {SCOPE}:
                            raise PublishError("OAuth grant must contain only drive.file; revoke broader grants and authorize again")
                        self.scope_verified = True
                    self.access_token = token["access_token"]
                    self.expires_at = time.monotonic() + max(0, int(token.get("expires_in", 3600)) - 60)
                    return
                if response.status not in RETRYABLE:
                    raise PublishError("OAuth refresh failed (HTTP {}); check publisher credentials".format(response.status))
            except TransportError:
                pass
            if attempt + 1 == self.attempts:
                raise PublishError("OAuth refresh failed after bounded retries")
            self.pause(attempt, response)

    def raw(self, method, url, body=None, headers=None):
        if not self.access_token or time.monotonic() >= self.expires_at:
            self.refresh()
        headers = dict(headers or {})
        headers["Authorization"] = "Bearer " + self.access_token
        return self.transport(method, url, headers, body)

    def request(self, method, url, body=None, headers=None, accepted=(200,)):
        refreshed = False
        for attempt in range(self.attempts):
            response = None
            try:
                response = self.raw(method, url, body, headers)
                if response.status in accepted:
                    return response
                if response.status == 401 and not refreshed:
                    self.refresh()
                    refreshed = True
                    continue
                if not self.retryable(response):
                    raise ApiError(response.status, error_reason(response))
            except TransportError:
                pass
            if attempt + 1 == self.attempts:
                raise PublishError("Drive request failed after bounded retries")
            self.pause(attempt, response)
        raise PublishError("Drive request could not be authorized")

    def get(self, file_id):
        query = urllib.parse.urlencode({"fields": FIELDS, "supportsAllDrives": "true"})
        return decode_json(self.request("GET", API + "/" + valid_id(file_id) + "?" + query))

    def find(self, query):
        result = []
        page = None
        while True:
            params = {"q": query + " and trashed = false", "spaces": "drive", "pageSize": 100,
                      "fields": "nextPageToken,files({})".format(FIELDS),
                      "supportsAllDrives": "true", "includeItemsFromAllDrives": "true"}
            if page:
                params["pageToken"] = page
            value = decode_json(self.request("GET", API + "?" + urllib.parse.urlencode(params)))
            files = value.get("files")
            if not isinstance(files, list) or not all(isinstance(item, dict) for item in files):
                raise PublishError("Drive list response is missing files")
            result.extend(files)
            page = value.get("nextPageToken")
            if not page:
                return result

    def allocate_id(self):
        response = self.request("GET", API + "/generateIds?count=1&space=drive&type=files")
        value = decode_json(response).get("ids", [])
        if len(value) != 1:
            raise PublishError("Drive did not allocate a file ID")
        return valid_id(value[0])

    def create_folder(self, name, parent=None, role=None):
        file_id = self.allocate_id()
        metadata = {"id": file_id, "name": name, "mimeType": FOLDER_TYPE}
        if parent:
            metadata["parents"] = [valid_id(parent)]
        if role:
            metadata["appProperties"] = {"cf_role": role}
        response = self.request("POST", API + "?fields=" + urllib.parse.quote(FIELDS),
                                json.dumps(metadata).encode("utf-8"),
                                {"Content-Type": "application/json"}, accepted=(200, 201, 409))
        value = self.get(file_id) if response.status == 409 else decode_json(response)
        if value.get("id") != file_id or value.get("name") != name or value.get("mimeType") != FOLDER_TYPE:
            raise PublishError("Created Drive folder could not be verified")
        return value

    @staticmethod
    def session_url(response):
        value = header(response, "Location")
        parsed = urllib.parse.urlparse(value or "")
        if (parsed.scheme != "https" or parsed.hostname not in
                {"www.googleapis.com", "content.googleapis.com"} or parsed.username or
                parsed.password or parsed.port not in {None, 443} or
                not parsed.path.startswith("/upload/drive/v3/files")):
            raise PublishError("Drive returned an invalid resumable session URL")
        return value

    @staticmethod
    def received_bytes(response, total):
        received = header(response, "Range")
        if received is None:
            return 0
        match = re.fullmatch(r"bytes=0-(\d+)", received.strip())
        if not match or int(match.group(1)) >= total:
            raise PublishError("Invalid byte range in resumable upload response")
        return int(match.group(1)) + 1

    def upload(self, bundle, metadata, file_id, create=False):
        """Upload content; resolve a lost final acknowledgement by querying its session."""
        file_id = valid_id(file_id)
        if create:
            metadata = dict(metadata, id=file_id)

        def confirmed(value):
            verify_pdf(value, bundle)
            if (value.get("id") != file_id or
                    any(value.get("appProperties", {}).get(key) != expected
                        for key, expected in metadata.get("appProperties", {}).items()) or
                    ("name" in metadata and value.get("name") != metadata["name"]) or
                    any(parent not in value.get("parents", []) for parent in metadata.get("parents", []))):
                raise PublishError("Uploaded PDF metadata does not match the requested generation")
            return value

        query = urllib.parse.urlencode({"uploadType": "resumable", "supportsAllDrives": "true",
                                        "fields": FIELDS})
        initial_url = UPLOAD_API + ("" if create else "/" + file_id) + "?" + query
        total = bundle.size
        for restart in range(self.attempts):
            response = self.request("POST" if create else "PATCH", initial_url,
                                    json.dumps(metadata).encode("utf-8"),
                                    {"Content-Type": "application/json; charset=UTF-8",
                                     "X-Upload-Content-Type": "application/pdf",
                                     "X-Upload-Content-Length": str(total)},
                                    accepted=(200, 201, 409) if create else (200, 201))
            if response.status == 409:
                return confirmed(self.get(file_id))
            session = self.session_url(response)
            offset = 0
            failures = 0
            with bundle.pdf.open("rb") as source:
                while offset < total:
                    source.seek(offset)
                    chunk = source.read(min(self.chunk_size, total - offset))
                    if not chunk:
                        raise PublishError("Local PDF changed during upload")
                    response = None
                    try:
                        response = self.raw("PUT", session, chunk,
                                            {"Content-Type": "application/pdf",
                                             "Content-Length": str(len(chunk)),
                                             "Content-Range": "bytes {}-{}/{}".format(
                                                 offset, offset + len(chunk) - 1, total)})
                    except TransportError:
                        pass
                    if response and response.status in (200, 201):
                        return confirmed(decode_json(response) if response.body else self.get(file_id))
                    if response and response.status == 308:
                        next_offset = self.received_bytes(response, total)
                        if next_offset > offset:
                            failures = 0
                            offset = next_offset
                            continue
                    elif response and response.status == 401:
                        self.refresh()
                    elif response and response.status == 404:
                        break
                    elif response and not self.retryable(response):
                        raise ApiError(response.status, error_reason(response))
                    failures += 1
                    if failures > self.attempts:
                        break
                    self.pause(failures - 1, response)
                    # Never blindly resend a chunk after a timeout: it may already be committed.
                    probe = self.request("PUT", session, b"",
                                         {"Content-Length": "0",
                                          "Content-Range": "bytes */{}".format(total)},
                                         accepted=(200, 201, 308, 404))
                    if probe.status in (200, 201):
                        return confirmed(decode_json(probe) if probe.body else self.get(file_id))
                    if probe.status == 404:
                        break
                    offset = self.received_bytes(probe, total)
            # A session can expire or its last acknowledgement can be lost.  A matching
            # committed file resolves either case before a new session is opened.
            try:
                return confirmed(self.get(file_id))
            except ApiError as error:
                if error.status != 404:
                    raise
            except PublishError:
                pass
            if restart + 1 < self.attempts:
                self.pause(restart)
        raise PublishError("Upload could not be confirmed after bounded retries; inspect remote state before retrying")


@dataclasses.dataclass
class Bundle:
    manifest: dict
    pdf: Path
    size: int
    sha256: str
    md5: str
    submitted: dt.datetime
    completed: dt.datetime

    @property
    def job_hash(self):
        identity = self.manifest["repository"] + "\0" + self.manifest["job_id"]
        return hashlib.sha256(identity.encode("utf-8")).hexdigest()


def parse_timestamp(value):
    try:
        if not isinstance(value, str) or len(value) > 64:
            raise ValueError()
        parsed = dt.datetime.fromisoformat(value.replace("Z", "+00:00"))
        if parsed.tzinfo is None or parsed.utcoffset() is None:
            raise ValueError()
        return parsed.astimezone(dt.timezone.utc)
    except (ValueError, OverflowError):
        raise PublishError("Manifest timestamps must be ISO 8601 with an explicit time zone") from None


def normalized_timestamp(value):
    return value.isoformat(timespec="microseconds").replace("+00:00", "Z")


def validate_report_quality(manifest):
    """Permit explicit fit-quality warnings without hiding computation failures."""
    quality_fields = {"report_quality", "computation_exit_code", "fit_quality_warnings"}
    legacy = not any(key in manifest for key in quality_fields)
    if not legacy and not all(key in manifest for key in quality_fields):
        raise PublishError("Manifest report quality, computation exit code and warnings must appear together")
    quality = manifest.get("report_quality", "complete")
    computation_exit = manifest.get("computation_exit_code", 0)
    warnings = manifest.get("fit_quality_warnings", [])
    if not isinstance(quality, str) or quality not in {"complete", "incomplete"}:
        raise PublishError("Manifest report_quality must be complete or incomplete")
    if type(computation_exit) is not int or computation_exit not in {0, 2}:
        raise PublishError("Manifest computation_exit_code must be 0 or 2")
    if not isinstance(warnings, list):
        raise PublishError("Manifest fit_quality_warnings must be an array")
    for warning in warnings:
        if (not isinstance(warning, dict) or
                not isinstance(warning.get("stage"), str) or not warning["stage"].strip() or
                type(warning.get("exit_code")) is not int or warning["exit_code"] != 2 or
                not isinstance(warning.get("error"), str) or not warning["error"].strip()):
            raise PublishError("Every fit-quality warning must identify a stage, exit code 2 and error")
    if quality == "incomplete":
        if computation_exit != 2 or not warnings:
            raise PublishError("Incomplete report requires computation exit code 2 and fit-quality warnings")
    elif computation_exit != 0 or warnings:
        raise PublishError("Complete report cannot conceal a computation failure or fit-quality warnings")

    scientific = manifest.get("scientific_manifests")
    if scientific is None and "scientific_manifests" not in manifest and quality == "complete":
        return
    if not isinstance(scientific, list) or quality == "incomplete" and not scientific:
        raise PublishError("Manifest scientific_manifests must list the scientific run statuses")
    incomplete_count = 0
    for record in scientific:
        if not isinstance(record, dict):
            raise PublishError("Manifest scientific run record must be an object")
        # Older successful bundles saved only config/run ID/digest references.
        # This compatibility path never accepts an explicit incomplete status.
        if legacy and "status" not in record and "exit_code" not in record:
            continue
        status, code = record.get("status"), record.get("exit_code")
        if (type(code) is not int or not isinstance(status, str) or
                (status, code) not in {("completed", 0), ("incomplete", 2)}):
            raise PublishError("Scientific run status and exit code must be completed/0 or incomplete/2")
        if status == "incomplete":
            incomplete_count += 1
    if quality == "incomplete" and incomplete_count == 0:
        raise PublishError("Incomplete report must preserve at least one incomplete scientific run")
    if quality == "complete" and incomplete_count:
        raise PublishError("Complete report cannot contain an incomplete scientific run")


def load_bundle(directory):
    directory = Path(directory)
    manifest_path = directory / "manifest.json"
    pdf = directory / "report.pdf"
    if manifest_path.is_symlink() or pdf.is_symlink():
        raise PublishError("Bundle files must be regular files, not symbolic links")
    try:
        if not manifest_path.is_file() or manifest_path.stat().st_size > 1024 * 1024:
            raise PublishError("Bundle manifest is missing or too large")
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if not isinstance(manifest, dict):
            raise ValueError()
    except (OSError, ValueError, UnicodeError):
        raise PublishError("Cannot read a valid bundle manifest") from None
    for key in ("job_id", "repository", "run_id"):
        if (not isinstance(manifest.get(key), str) or not manifest[key].strip() or
                len(manifest[key].encode("utf-8")) > (96 if key == "run_id" else 512)):
            raise PublishError("Manifest has an invalid " + key)
    if not re.fullmatch(r"[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+", manifest["repository"]):
        raise PublishError("Manifest repository must be owner/repository")
    if (type(manifest.get("schema_version")) is not int or manifest["schema_version"] != 1 or
            type(manifest.get("exit_code")) is not int or manifest["exit_code"] != 0 or
            manifest.get("status") != "uploading" or manifest.get("pdf_path") != "report.pdf"):
        raise PublishError("Only finalized schema-1 report.pdf bundles in uploading state may be published")
    validate_report_quality(manifest)
    if (not isinstance(manifest.get("commit_sha"), str) or
            not re.fullmatch(r"[0-9a-fA-F]{40}", manifest["commit_sha"]) or
            not isinstance(manifest.get("pdf_sha256"), str) or
            not re.fullmatch(r"[0-9a-fA-F]{64}", manifest["pdf_sha256"]) or
            type(manifest.get("pdf_size")) is not int or manifest["pdf_size"] <= 0):
        raise PublishError("Manifest has invalid commit/PDF hashes or PDF size")
    submitted = parse_timestamp(manifest.get("submitted_at"))
    completed = parse_timestamp(manifest.get("completed_at"))
    if completed < submitted:
        raise PublishError("Manifest completion precedes submission")
    sha256 = hashlib.sha256()
    md5 = hashlib.md5()
    size = 0
    try:
        if not pdf.is_file():
            raise PublishError("Bundle report.pdf is missing")
        with pdf.open("rb") as source:
            if source.read(5) != b"%PDF-":
                raise PublishError("Bundle report.pdf has no PDF signature")
            source.seek(0)
            for chunk in iter(lambda: source.read(1024 * 1024), b""):
                sha256.update(chunk)
                md5.update(chunk)
                size += len(chunk)
    except OSError:
        raise PublishError("Cannot read bundle report.pdf") from None
    if size != manifest["pdf_size"] or sha256.hexdigest() != manifest["pdf_sha256"].lower():
        raise PublishError("Local PDF size/checksum does not match the manifest")
    return Bundle(manifest, pdf, size, sha256.hexdigest(), md5.hexdigest(), submitted, completed)


def verify_pdf(remote, bundle):
    if (remote.get("trashed") or remote.get("mimeType") != "application/pdf" or
            str(remote.get("size", "")) != str(bundle.size) or
            str(remote.get("md5Checksum", "")).lower() != bundle.md5 or
            (remote.get("sha256Checksum") and str(remote["sha256Checksum"]).lower() != bundle.sha256)):
        raise PublishError("Remote PDF size/checksum/type could not be verified")
    valid_id(remote.get("id"))


def query_literal(value):
    return "'" + value.replace("\\", "\\\\").replace("'", "\\'") + "'"


def require_folder(client, file_id, parent=None):
    value = client.get(valid_id(file_id))
    if (value.get("trashed") or value.get("mimeType") != FOLDER_TYPE or
            (parent and parent not in value.get("parents", []))):
        raise PublishError("Configured Drive folder is missing, trashed, or outside the report root")
    return value


def initialize(client):
    roots = client.find("'root' in parents and mimeType = '{}' and "
                        "appProperties has {{ key='cf_role' and value='root' }}".format(FOLDER_TYPE))
    if len(roots) > 1:
        raise PublishError("Multiple application report roots exist; choose one explicitly")
    root = roots[0] if roots else client.create_folder("LabReports", role="root")
    valid_id(root.get("id"))
    histories = client.find("'{}' in parents and mimeType = '{}' and "
                            "appProperties has {{ key='cf_role' and value='history_folder' }}".format(
                                root["id"], FOLDER_TYPE))
    if len(histories) > 1:
        raise PublishError("Multiple history folders exist inside the report root")
    history = histories[0] if histories else client.create_folder("history", root["id"], "history_folder")
    return {"root_folder_id": root["id"], "history_folder_id": valid_id(history.get("id")),
            "scope": SCOPE,
            "scope_verification": "refresh_response" if client.scope_verified else
            "not_returned_by_google; ensure the original grant requested only drive.file"}


def run_order(submitted, run_id):
    # GitHub numeric run IDs need numeric comparison, e.g. 10 must follow 9.
    return submitted, (1, int(run_id)) if run_id.isdecimal() else (0, run_id)


def publish_bundle(client, bundle, root_id, history_id):
    root_id, history_id = valid_id(root_id), valid_id(history_id)
    if root_id == history_id:
        raise PublishError("Root and history folders must be distinct")
    require_folder(client, root_id)
    require_folder(client, history_id, root_id)
    props = {"cf_role": "history", "job_id": bundle.job_hash, "pdf_sha256": bundle.sha256,
             "commit_sha": bundle.manifest["commit_sha"].lower(),
             "submitted_at": normalized_timestamp(bundle.submitted),
             "run_id": bundle.manifest["run_id"]}
    archives = client.find("'{}' in parents and appProperties has {{ key='job_id' and value='{}' }}".format(
        history_id, bundle.job_hash))
    if len(archives) > 1:
        raise PublishError("Multiple archive files exist for one job; refusing ambiguous publication")
    if archives:
        archive = client.get(valid_id(archives[0].get("id")))
        verify_pdf(archive, bundle)
        saved = archive.get("appProperties", {})
        if (any(saved.get(key) != value for key, value in props.items()) or
                history_id not in archive.get("parents", [])):
            raise PublishError("This job was already archived with a different identity or PDF")
    else:
        safe_job = re.sub(r"[^A-Za-z0-9_-]", "_", bundle.manifest["job_id"])[:48]
        name = "{}_{}_{}.pdf".format(bundle.completed.strftime("%Y%m%dT%H%M%SZ"),
                                     bundle.manifest["commit_sha"][:12], safe_job)
        archive_id = client.allocate_id()
        metadata = {"name": name, "mimeType": "application/pdf", "parents": [history_id],
                    "appProperties": props}
        client.upload(bundle, metadata, archive_id, create=True)
        archive = client.get(archive_id)
        verify_pdf(archive, bundle)
        if archive.get("appProperties") != props or history_id not in archive.get("parents", []):
            raise PublishError("Archive metadata could not be verified; latest was not requested")
    # This is deliberately after successful archive verification.  History errors
    # cannot even initiate a latest upload or change its metadata.
    candidates = client.find("'{}' in parents and name = 'latest.pdf'".format(root_id))
    if len(candidates) > 1:
        raise PublishError("Multiple latest.pdf files exist; refusing to guess a file ID")
    latest = client.get(valid_id(candidates[0].get("id"))) if candidates else None
    latest_props = dict(props, cf_role="latest")
    status = "published"
    if latest:
        previous = latest.get("appProperties", {})
        if (latest.get("trashed") or latest.get("mimeType") != "application/pdf" or previous.get("cf_role") != "latest" or
                root_id not in latest.get("parents", [])):
            raise PublishError("Existing latest.pdf is not a recognized publisher-owned PDF")
        if "submitted_at" not in previous or not isinstance(previous.get("run_id"), str):
            raise PublishError("Existing latest.pdf has no valid publication generation")
        current_order = run_order(parse_timestamp(previous["submitted_at"]), previous["run_id"])
        candidate_order = run_order(bundle.submitted, bundle.manifest["run_id"])
        if current_order > candidate_order:
            status = "archived_stale"
        elif current_order == candidate_order:
            verify_pdf(latest, bundle)
            if any(previous.get(key) != value for key, value in latest_props.items()):
                raise PublishError("One publication generation refers to different jobs/content")
            status = "already_published"
    if status == "published":
        latest_id = latest["id"] if latest else client.allocate_id()
        metadata = {"name": "latest.pdf", "mimeType": "application/pdf", "appProperties": latest_props}
        if not latest:
            metadata["parents"] = [root_id]
        try:
            client.upload(bundle, metadata, latest_id, create=latest is None)
        except PublishError:
            # The last response may have been lost after commit.  Never claim that
            # an error means latest is unchanged; verify the actual remote file.
            try:
                committed = client.get(latest_id)
                verify_pdf(committed, bundle)
                if any(committed.get("appProperties", {}).get(key) != value
                       for key, value in latest_props.items()):
                    raise PublishError("Publication metadata did not commit")
            except PublishError:
                raise PublishError("Latest upload could not be confirmed; history is preserved, inspect remote state before retrying") from None
        latest = client.get(latest_id)
        verify_pdf(latest, bundle)
        if (any(latest.get("appProperties", {}).get(key) != value for key, value in latest_props.items()) or
                root_id not in latest.get("parents", [])):
            raise PublishError("Latest metadata could not be verified; inspect remote state")
    latest_id = valid_id(latest.get("id"))
    return {"status": status, "version_id": valid_id(archive.get("id")), "latest_id": latest_id,
            "latest_url": "https://drive.google.com/file/d/{}/view".format(latest_id)}


def required_env(name):
    value = os.environ.get(name, "")
    if not value:
        raise PublishError("Required publisher environment variable is missing: " + name)
    return value


def emit_outputs(values):
    path = os.environ.get("GITHUB_OUTPUT")
    if path:
        try:
            with open(path, "a", encoding="utf-8") as output:
                for key, value in values.items():
                    if "\n" in str(value) or "\r" in str(value):
                        raise PublishError("Refusing a multiline workflow output")
                    output.write("{}={}\n".format(key, value))
        except OSError:
            raise PublishError("Cannot write GitHub Actions outputs") from None
    print(json.dumps(values, sort_keys=True))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("init", help="Create private app-owned LabReports/history folders")
    publish = subparsers.add_parser("publish", help="Archive a completed report and update stable latest.pdf")
    publish.add_argument("--bundle", required=True, type=Path)
    args = parser.parse_args(argv)
    try:
        # Validate local data before acquiring any access token or making network requests.
        bundle = load_bundle(args.bundle) if args.command == "publish" else None
        client = DriveClient(required_env("CF_DRIVE_CLIENT_ID"), required_env("CF_DRIVE_CLIENT_SECRET"),
                             required_env("CF_DRIVE_REFRESH_TOKEN"))
        if args.command == "init":
            values = initialize(client)
            print("CF_DRIVE_ROOT_FOLDER_ID=" + values["root_folder_id"])
            print("CF_DRIVE_HISTORY_FOLDER_ID=" + values["history_folder_id"])
            emit_outputs(values)
        else:
            values = publish_bundle(client, bundle, required_env("CF_DRIVE_ROOT_FOLDER_ID"),
                                    required_env("CF_DRIVE_HISTORY_FOLDER_ID"))
            emit_outputs(values)
        return 0
    except PublishError as error:
        print("publish: " + str(error), file=sys.stderr)
        return 1
    except Exception as error:
        # Avoid accidentally printing an HTTP library exception containing tokens/URLs.
        print("publish: unexpected {} (details omitted to protect credentials)".format(type(error).__name__),
              file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
