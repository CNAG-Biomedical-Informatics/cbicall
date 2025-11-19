from pathlib import Path

import pytest

from cbicall import dnaseq


def test_submit_cmd_success(monkeypatch, tmp_path):
    called = {}

    def fake_call(cmd, shell):
        called["cmd"] = cmd
        return 0

    monkeypatch.setattr(dnaseq.subprocess, "call", fake_call)
    log_path = tmp_path / "log.txt"

    # Should not raise
    dnaseq.DNAseq._submit_cmd("echo ok", log_path, "RUNOK", debug=False)

    assert "echo ok" in called["cmd"]


def test_submit_cmd_non_zero_raises(monkeypatch, tmp_path):
    def fake_call(cmd, shell):
        return 1

    monkeypatch.setattr(dnaseq.subprocess, "call", fake_call)
    log_path = tmp_path / "log.txt"

    with pytest.raises(RuntimeError) as excinfo:
        dnaseq.DNAseq._submit_cmd("echo fail", log_path, "RUNFAIL", debug=False)

    msg = str(excinfo.value)
    assert "RUNFAIL" in msg
    assert str(log_path) in msg


def test_submit_cmd_exception_raises(monkeypatch, tmp_path):
    def fake_call(cmd, shell):
        raise OSError("boom")

    monkeypatch.setattr(dnaseq.subprocess, "call", fake_call)
    log_path = tmp_path / "log.txt"

    with pytest.raises(RuntimeError) as excinfo:
        dnaseq.DNAseq._submit_cmd("echo fail", log_path, "RUNERR", debug=True)

    msg = str(excinfo.value)
    assert "RUNERR" in msg
    assert str(log_path) in msg

