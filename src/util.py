"""Shared utility functions for cftk."""

import os
import sys
import time
import re
import json
import shlex
import subprocess
import threading
import uuid
from datetime import datetime, timezone

try:
    import fcntl
except ImportError:  # pragma: no cover - Windows does not provide fcntl.
    fcntl = None


_COMMAND_LOG_ENV = "CFTK_COMMAND_LOG"
_COMMAND_RUN_ID_ENV = "CFTK_COMMAND_RUN_ID"
_COMMAND_LOG_PATH = None
_COMMAND_RUN_ID = None
_COMMAND_LOG_LOCK = threading.Lock()


def disp(msg):
    print(f"@{time.asctime()}\t{msg}", file=sys.stderr)


def configure_command_log(path):
    """Configure the project-scoped JSONL ledger for external commands.

    ``None`` disables recording and is useful for isolated callers/tests.
    The path is also inherited by multiprocessing workers through the
    environment, while the run ID remains stable for the current workflow.
    """
    global _COMMAND_LOG_PATH, _COMMAND_RUN_ID
    if path is None:
        _COMMAND_LOG_PATH = None
        _COMMAND_RUN_ID = None
        os.environ.pop(_COMMAND_LOG_ENV, None)
        os.environ.pop(_COMMAND_RUN_ID_ENV, None)
        return None

    configured_path = os.path.abspath(os.fspath(path))
    os.makedirs(os.path.dirname(configured_path), exist_ok=True)
    if configured_path != _COMMAND_LOG_PATH or _COMMAND_RUN_ID is None:
        _COMMAND_RUN_ID = uuid.uuid4().hex
    _COMMAND_LOG_PATH = configured_path
    os.environ[_COMMAND_LOG_ENV] = _COMMAND_LOG_PATH
    os.environ[_COMMAND_RUN_ID_ENV] = _COMMAND_RUN_ID
    return _COMMAND_LOG_PATH


def _command_log_path():
    return _COMMAND_LOG_PATH or os.environ.get(_COMMAND_LOG_ENV)


def _command_text(command):
    if isinstance(command, str):
        return command
    return shlex.join(str(value) for value in command)


def _write_command_record(record):
    path = _command_log_path()
    if not path:
        return
    line = json.dumps(record, sort_keys=True, separators=(",", ":")) + "\n"
    try:
        with _COMMAND_LOG_LOCK:
            with open(path, "a", encoding="utf-8") as handle:
                if fcntl is not None:
                    fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
                handle.write(line)
                handle.flush()
                if fcntl is not None:
                    fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
    except OSError as exc:
        raise RuntimeError(
            f"[util] ERROR: could not write command provenance {path}: {exc}"
        ) from exc


def recorded_run(command, *args, label="", **kwargs):
    """Run ``subprocess.run`` and record its exact invocation and outcome.

    All positional and keyword arguments are passed through unchanged. This
    keeps captured-output and list-argv callers compatible with their existing
    behavior while extending the same provenance ledger used by
    :func:`run_command`.
    """
    command_id = uuid.uuid4().hex
    command_value = _command_text(command)
    run_id = _COMMAND_RUN_ID or os.environ.get(_COMMAND_RUN_ID_ENV)
    if run_id is None:
        run_id = "unconfigured"
    started = datetime.now(timezone.utc).isoformat()
    run_cwd = kwargs.get("cwd") or os.getcwd()
    common = {
        "command_id": command_id,
        "command": command_value,
        "cwd": os.path.abspath(os.fspath(run_cwd)),
        "label": label,
        "run_id": run_id,
        "shell": bool(kwargs.get("shell", False)),
    }
    _write_command_record({**common, "event": "start", "timestamp": started})
    try:
        completed = subprocess.run(command, *args, **kwargs)
    except BaseException as exc:
        _write_command_record({
            **common,
            "event": "finish",
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "returncode": None,
            "error": f"{type(exc).__name__}: {exc}",
        })
        raise
    _write_command_record({
        **common,
        "event": "finish",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "returncode": completed.returncode,
    })
    return completed


def run_command(cmd, label="", check=True):
    disp(f"CMD [{label}]: {cmd[:120]}" if label else f"CMD: {cmd[:120]}")
    ret = recorded_run(cmd, shell=True, label=label)
    if check and ret.returncode != 0:
        sys.exit(f"[util] ERROR: command failed — {label or cmd[:80]}")
    return ret.returncode


def is_number(s):
    return bool(re.match(r"^-?\d+(?:\.\d+)?$", str(s)))
