from pathlib import Path

from cbicall import cli as cli_mod


def test_main_happy_path(monkeypatch, tmp_path):
    # Stub usage() to avoid touching real sys.argv / helpmod logic here
    def fake_usage(version):
        return {
            "threads": 4,
            "paramfile": str(tmp_path / "params.yaml"),
            "debug": 0,
            "verbose": False,
            "noemoji": False,
            "nocolor": False,
        }

    monkeypatch.setattr(cli_mod, "usage", fake_usage)

    # Stub config.read_param_file and set_config_values
    fake_param = {
        "pipeline": "wes",
        "mode": "single",
        "sample": None,
        "sample_map": None,
        "workflow_engine": "bash",
        "gatk_version": "gatk-3.5",
        "cleanup_bam": False,
    }

    def fake_read_param_file(path):
        return fake_param

    def fake_set_config_values(param):
        # Minimal config that main() expects
        return {
            "projectdir": str(tmp_path / "proj"),
            "id": "ID123",
            "bash_wes_single": "/path/to/bash_wes_single.sh",
        }

    monkeypatch.setattr(cli_mod.config_mod, "read_param_file", fake_read_param_file)
    monkeypatch.setattr(cli_mod.config_mod, "set_config_values", fake_set_config_values)

    # Stub write_log so we do not actually write files
    logs = {}

    def fake_write_log(cfg, arg, param):
        logs["cfg"] = cfg
        logs["arg"] = arg
        logs["param"] = param

    monkeypatch.setattr(cli_mod, "write_log", fake_write_log)

    # Fake DNAseq to avoid running real pipeline
    class FakeDNAseq:
        def __init__(self, settings):
            self.settings = settings

        def variant_calling(self):
            # Just record that we were called
            logs["settings"] = self.settings

    monkeypatch.setattr(cli_mod, "DNAseq", FakeDNAseq)

    # Fake GoodBye to avoid randomness
    class FakeGoodBye:
        def say_goodbye(self):
            return "Bye"

    monkeypatch.setattr(cli_mod, "GoodBye", FakeGoodBye)

    rc = cli_mod.main()
    assert rc == 0

    # Sanity checks on what main() passed down
    settings = logs["settings"]
    assert settings["projectdir"] == str(tmp_path / "proj")
    assert settings["pipeline"] == "wes"
    assert settings["mode"] == "single"
    assert settings["workflow_engine"] == "bash"
    assert settings["threads"] == 4
    assert settings["id"] == "ID123"

