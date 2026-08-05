import json
from concurrent.futures import ThreadPoolExecutor
from types import SimpleNamespace

import pytest


@pytest.fixture
def util_module(monkeypatch):
    monkeypatch.syspath_prepend("src")
    import util

    util.configure_command_log(None)
    yield util
    util.configure_command_log(None)


def _records(path):
    return [json.loads(line) for line in path.read_text().splitlines()]


def test_run_command_records_full_command_before_and_after_execution(
    util_module, monkeypatch, tmp_path
):
    log_path = tmp_path / "provenance" / "commands.jsonl"
    util_module.configure_command_log(log_path)
    command = "tool " + ("x" * 180)
    monkeypatch.setattr(
        util_module.subprocess,
        "run",
        lambda *args, **kwargs: SimpleNamespace(returncode=0),
    )

    assert util_module.run_command(command, label="long command") == 0

    records = _records(log_path)
    assert [record["event"] for record in records] == ["start", "finish"]
    assert records[0]["command"] == command
    assert records[1]["command"] == command
    assert records[0]["command_id"] == records[1]["command_id"]
    assert records[1]["returncode"] == 0
    assert records[0]["label"] == "long command"


def test_failed_command_has_completion_record(util_module, monkeypatch, tmp_path):
    log_path = tmp_path / "commands.jsonl"
    util_module.configure_command_log(log_path)
    monkeypatch.setattr(
        util_module.subprocess,
        "run",
        lambda *args, **kwargs: SimpleNamespace(returncode=17),
    )

    with pytest.raises(SystemExit, match="command failed"):
        util_module.run_command("tool --fails", label="failure")

    records = _records(log_path)
    assert records[-1]["event"] == "finish"
    assert records[-1]["returncode"] == 17


def test_list_argv_and_parallel_writes_are_replayable(util_module, monkeypatch, tmp_path):
    log_path = tmp_path / "commands.jsonl"
    util_module.configure_command_log(log_path)
    monkeypatch.setattr(
        util_module.subprocess,
        "run",
        lambda *args, **kwargs: SimpleNamespace(returncode=0),
    )

    def run_one(index):
        return util_module.recorded_run(
            ["tool", "--sample", f"sample {index}"],
            label=f"sample {index}",
        ).returncode

    with ThreadPoolExecutor(max_workers=8) as executor:
        assert list(executor.map(run_one, range(20))) == [0] * 20

    records = _records(log_path)
    assert len(records) == 40
    assert all(record["command"].startswith("tool --sample ") for record in records)
    assert all(record["run_id"] != "unconfigured" for record in records)
