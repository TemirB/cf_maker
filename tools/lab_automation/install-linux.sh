#!/usr/bin/env bash
# Run only on an audited Linux lab account. This script never registers a runner.
set -euo pipefail
umask 077

usage() {
  echo 'Usage: install-linux.sh --config /absolute/config.json [--install-dir /absolute/path]'
}

config_path=
install_dir="${XDG_DATA_HOME:-"$HOME/.local/share"}/cf-maker-lab"
while [ "$#" -gt 0 ]; do
  case "$1" in
    --config|--install-dir)
      [ "$#" -ge 2 ] || { usage >&2; exit 2; }
      if [ "$1" = --config ]; then config_path=$2; else install_dir=$2; fi
      shift 2
      ;;
    -h|--help) usage; exit 0 ;;
    *) usage >&2; exit 2 ;;
  esac
done
[ -n "$config_path" ] || { usage >&2; exit 2; }
[ "$(uname -s)" = Linux ] || { echo 'Linux with a systemd user manager is required.' >&2; exit 1; }
[ "$(id -u)" != 0 ] || { echo 'Run as the dedicated, unprivileged lab account, without sudo.' >&2; exit 1; }
command -v python3 >/dev/null
command -v systemctl >/dev/null
command -v loginctl >/dev/null
systemctl --user show-environment >/dev/null
if [ "$(loginctl show-user "$(id -u)" --property=Linger --value)" != yes ]; then
  echo 'User lingering is disabled. Ask the administrator to enable it before installation.' >&2
  exit 1
fi
for unit in cf-maker-worker.path cf-maker-worker.service; do
  unit_state=$(systemctl --user show "$unit" --property=ActiveState --value 2>/dev/null || true)
  case "$unit_state" in
    active|activating|reloading)
      echo 'Stop cf-maker-worker.path and cf-maker-worker.service before replacing installed helpers.' >&2
      exit 1
      ;;
  esac
done

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd -P)
unit_dir="${XDG_CONFIG_HOME:-"$HOME/.config"}/systemd/user"
python_path=$(command -v python3)
python3 - "$config_path" "$install_dir" "$unit_dir" "$script_dir" "$python_path" <<'PY'
import json
import os
from pathlib import Path
import sys

if sys.version_info < (3, 9):
    raise SystemExit("Python 3.9 or newer is required")
config_arg, install_arg, unit_arg, source_arg, python_arg = sys.argv[1:]
source = Path(source_arg).resolve(strict=True)
repository = source.parent.parent

def checked_path(value, label, owned=True):
    path = Path(value)
    if not path.is_absolute():
        raise SystemExit(f"{label} must be an absolute path")
    path = path.resolve()
    # Values are embedded in systemd directives; allow spaces without accepting
    # specifiers, environment expansion, or additional unit-file directives.
    if any(ord(character) < 32 or character in '%$"\\' for character in str(path)):
        raise SystemExit(f"{label} contains characters unsupported by these unit templates")
    if path == repository or repository in path.parents:
        raise SystemExit(f"{label} must be outside the checkout")
    if owned and path.exists() and path.stat().st_uid != os.getuid():
        raise SystemExit(f"{label} must belong to the current user")
    return path

config = checked_path(config_arg, "config")
install = checked_path(install_arg, "install directory")
units = checked_path(unit_arg, "systemd user directory")
python = checked_path(python_arg, "Python executable", owned=False)
if not config.is_file():
    raise SystemExit("config must be an existing JSON file")
data = json.loads(config.read_text(encoding="utf-8"))
if data.get("schema_version") != 1 or not isinstance(data.get("state_dir"), str):
    raise SystemExit("config must have schema_version=1 and an absolute state_dir")
state = checked_path(data["state_dir"], "state directory")
if install == state or install in state.parents or state in install.parents:
    raise SystemExit("install and state directories must be separate")
for directory in (install, units, state, state / "queue", state / "jobs"):
    directory.mkdir(parents=True, exist_ok=True, mode=0o700)
    directory.chmod(0o700)
config.chmod(0o600)

def write_atomic(destination, content, mode):
    if destination.exists() and destination.stat().st_uid != os.getuid():
        raise SystemExit(f"existing file belongs to another user: {destination}")
    temporary = destination.with_name(destination.name + f".tmp-{os.getpid()}")
    with temporary.open("xb") as target:
        target.write(content)
        target.flush()
        os.fsync(target.fileno())
    temporary.chmod(mode)
    os.replace(temporary, destination)

for name in ("lab.py", "audit.py"):
    helper = source / name
    if not helper.is_file():
        raise SystemExit(f"missing helper: {helper}")
    write_atomic(install / name, helper.read_bytes(), 0o700)
replacements = {
    "@INSTALL_DIR@": str(install), "@PYTHON@": str(python),
    "@CONFIG@": str(config), "@QUEUE_DIR@": str(state / "queue"),
}
for name in ("cf-maker-worker.service", "cf-maker-worker.path"):
    text = (source / "systemd" / (name + ".in")).read_text(encoding="utf-8")
    for marker, value in replacements.items():
        text = text.replace(marker, value)
    write_atomic(units / name, text.encode("utf-8"), 0o600)
print(f"Installed queue helpers: {install}")
print(f"State directory: {state}")
print(f"User units: {units}")
PY

systemctl --user daemon-reload
systemctl --user enable cf-maker-worker.service cf-maker-worker.path
systemctl --user --no-block start cf-maker-worker.service cf-maker-worker.path
echo 'Worker and queue watcher enabled for the current user. No runner was registered.'
