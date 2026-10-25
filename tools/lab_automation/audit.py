#!/usr/bin/env python3
"""Read-only machine/tool/data audit. Network probes are explicit and credential-free."""
import argparse
import datetime
import json
import os
from pathlib import Path
import platform
import re
import shutil
import socket
import subprocess
import urllib.error
import urllib.request


TOOLS = {
    "git": ["--version"], "cmake": ["--version"], "make": ["--version"],
    "g++": ["--version"], "clang++": ["--version"], "root-config": ["--version"],
    "latexmk": ["-v"], "xelatex": ["--version"], "biber": ["--version"],
    "gs": ["--version"], "systemctl": ["--version"], "rclone": ["version"],
}
PROBE_URLS = ["https://github.com", "https://api.github.com",
              "https://results-receiver.actions.githubusercontent.com",
              "https://objects.githubusercontent.com", "https://oauth2.googleapis.com",
              "https://www.googleapis.com/drive/v3/files"]


def probe(url):
    try:
        request = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "cf-maker-audit/1"})
        with urllib.request.urlopen(request, timeout=12) as response:
            code = response.status
        return {"url": url, "reachable": True, "http_status": code}
    except urllib.error.HTTPError as error:
        # 401/403/404 are evidence of an HTTP response, not permission to use the service.
        return {"url": url, "reachable": True, "http_status": error.code,
                "note": "HTTP response does not establish service authorization or VPN policy"}
    except (OSError, urllib.error.URLError) as error:
        return {"url": url, "reachable": False, "error_type": type(error).__name__}


def audit(source, config=None, network=False, policy_confirmed=False):
    source = Path(source).resolve()
    report = {"schema_version": 1, "audited_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
              "scope": "this machine only", "os": platform.system(), "release": platform.release(),
              "architecture": platform.machine(), "hostname": socket.gethostname(),
              "uid": os.getuid() if hasattr(os, "getuid") else None,
              "python": platform.python_version(), "source": str(source),
              "tools": {}, "inputs": [], "historical_durations": [],
              "network_policy_confirmed_by_operator": policy_confirmed,
              "network_probes": [], "unverified": [
                  "laboratory OS/network unless this script actually runs there",
                  "repository ruleset, runner registration and service startup",
                  "authenticated GitHub checkout/artifact upload through the laboratory VPN",
                  "Google OAuth grant and real Drive publication/revocation",
                  "systemd/logout/reboot behavior on the laboratory machine"]}
    setup = config.get("setup_script") if config else None
    for tool, arguments in TOOLS.items():
        command = [tool, *arguments]
        if setup:
            command = ["bash", "--noprofile", "--norc", "-c",
                       'set -e; source "$1"; shift; exec "$@"', "cf-maker-audit", setup, *command]
        try:
            result = subprocess.run(command, text=True, capture_output=True, timeout=15)
            lines = (result.stdout + result.stderr).strip().splitlines()
            report["tools"][tool] = {"found": result.returncode == 0,
                                     "exit_code": result.returncode, "version": lines[:3]}
        except (OSError, subprocess.TimeoutExpired) as error:
            report["tools"][tool] = {"found": False, "error_type": type(error).__name__}
    for name in (config or {}).get("config_files", ["config/kt.json", "config/y.json"]):
        try:
            app = json.loads((source / name).read_text())
            base = (config or {}).get("input_root") or app["machine"]["base_input"]
            path = Path(app["vars"]["input"]["file"])
            if not path.is_absolute():
                path = Path(base) / path
            if not path.is_absolute():
                path = source / path
            path = path.resolve()
            entry = {"config": name, "path": str(path), "exists": path.is_file(),
                     "readable": os.access(path, os.R_OK)}
            if path.is_file():
                entry["size_bytes"] = path.stat().st_size
            report["inputs"].append(entry)
            output = Path(app["machine"]["base_output"]) / app["vars"]["output"]["dir"]
            if not output.is_absolute():
                output = source / output
            log = output / app.get("logging", {}).get("file", "run.log")
            if log.is_file():
                with log.open("rb") as stream:
                    stream.seek(max(0, log.stat().st_size - 1024 * 1024))
                    text = stream.read().decode("utf-8", "replace")
                durations = re.findall(r"(Stage [^\n:]+): done in ([0-9.]+) s", text)
                report["historical_durations"].append({"config": name, "log": str(log),
                                                       "last_stages_seconds": durations[-12:]})
        except (OSError, ValueError, KeyError) as error:
            report["inputs"].append({"config": name, "error_type": type(error).__name__})
    disk_path = Path((config or {}).get("state_dir", source)).resolve()
    while not disk_path.exists():
        disk_path = disk_path.parent
    free = shutil.disk_usage(disk_path).free
    report["free_disk_bytes"] = free
    if network:
        report["network_probes"] = [probe(url) for url in PROBE_URLS]
        report["network_note"] = ("Probes do not prove access to wildcard Actions endpoints, actual "
                                  "runtime blob URLs or authenticated API operations. Verify runner diagnostics.")
    else:
        report["network_note"] = "Not probed; use --network on the laboratory PC after policy approval."
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", default=".")
    parser.add_argument("--config")
    parser.add_argument("--network", action="store_true")
    parser.add_argument("--policy-confirmed", action="store_true")
    parser.add_argument("--output")
    args = parser.parse_args()
    config = json.loads(Path(args.config).read_text()) if args.config else None
    report = audit(args.source, config, args.network, args.policy_confirmed)
    text = json.dumps(report, indent=2, ensure_ascii=False) + "\n"
    if args.output:
        Path(args.output).write_text(text)
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
