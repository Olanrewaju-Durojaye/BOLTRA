import json
from pathlib import Path

from boltzgen_workbench import runner


def test_resume_uses_original_run_with_reuse(tmp_path: Path, monkeypatch):
    output = tmp_path / "boltzgen_output"
    output.mkdir()
    (output / "steps.yaml").write_text("steps: [folding, analysis, filtering]\n", encoding="utf-8")
    (tmp_path / "design_specification.yaml").write_text("entities: []\n", encoding="utf-8")
    (tmp_path / "project_settings.json").write_text(
        json.dumps({"num_designs": 10, "budget": 2, "protocol": "peptide-anything"}), encoding="utf-8"
    )
    captured = {}
    monkeypatch.setattr(runner, "require_boltzgen", lambda: "/usr/bin/boltzgen")
    monkeypatch.setattr(runner, "_run_logged", lambda command, cwd, action, project: captured.setdefault("command", command) or 0)
    runner.resume_pipeline(output, ["folding", "analysis", "filtering"], tmp_path)
    command = captured["command"]
    assert command[1] == "run"
    assert "--reuse" in command
    assert command[command.index("--num_workers") + 1] == "0"
    assert command[command.index("--steps") + 1:] == ["folding", "analysis", "filtering", "--reuse"]


def test_worker_override_is_explicit(tmp_path: Path, monkeypatch):
    captured = {}
    monkeypatch.setenv("BOLTRA_NUM_WORKERS", "3")
    monkeypatch.setattr(runner, "require_boltzgen", lambda: "/usr/bin/boltzgen")
    monkeypatch.setattr(
        runner, "_run_logged",
        lambda command, cwd, action, project: captured.setdefault("command", command) or 0,
    )
    runner.run_pipeline(tmp_path / "design.yaml", tmp_path / "out", "protein-redesign", 10, 2, tmp_path)
    command = captured["command"]
    assert command[command.index("--num_workers") + 1] == "3"


def test_invalid_worker_override_is_rejected(monkeypatch):
    monkeypatch.setenv("BOLTRA_NUM_WORKERS", "many")
    try:
        runner.safe_worker_count()
    except ValueError as exc:
        assert "non-negative whole number" in str(exc)
    else:
        raise AssertionError("invalid worker count was accepted")


def test_low_swap_warning():
    snapshot = runner.ResourceSnapshot(30.0, 20.0, 2.0, 2.0, 100.0)
    warnings = runner.resource_warnings(snapshot, 20, 0)
    assert any("swap" in warning for warning in warnings)
