import json
from pathlib import Path

from boltzgen_workbench.status import inspect_project


def test_status_detects_first_incomplete_stage(tmp_path: Path):
    output = tmp_path / "boltzgen_output"
    output.mkdir()
    (output / "steps.yaml").write_text(
        "steps:\n  - design\n  - inverse_folding\n  - folding\n  - analysis\n  - filtering\n",
        encoding="utf-8",
    )
    (tmp_path / "project_settings.json").write_text(json.dumps({"num_designs": 2}), encoding="utf-8")
    designs = output / "intermediate_designs"
    designs.mkdir()
    for index in range(2):
        (designs / f"design_{index}.cif").write_text("data_test\n", encoding="utf-8")
    status = inspect_project(tmp_path)
    assert not status.complete
    assert status.stages[0].complete
    assert status.remaining_steps[0] == "inverse_folding"


def test_status_accepts_output_directory(tmp_path: Path):
    output = tmp_path / "boltzgen_output"
    output.mkdir()
    (output / "steps.yaml").write_text("steps: [design, filtering]\n", encoding="utf-8")
    status = inspect_project(output)
    assert status.output == output


def test_protein_status_tracks_design_folding(tmp_path: Path):
    output = tmp_path / "boltzgen_output"
    output.mkdir()
    (output / "steps.yaml").write_text(
        "steps: [design, inverse_folding, design_folding, folding, analysis, filtering]\n",
        encoding="utf-8",
    )
    (tmp_path / "project_settings.json").write_text(json.dumps({"num_designs": 2}), encoding="utf-8")
    for directory in (output / "intermediate_designs", output / "intermediate_designs_inverse_folded"):
        directory.mkdir()
        for index in range(2):
            (directory / f"design_{index}.cif").write_text("data_test\n", encoding="utf-8")
    status = inspect_project(tmp_path)
    assert status.remaining_steps[0] == "design_folding"


def test_saved_unstarted_project_is_recognized(tmp_path: Path):
    (tmp_path / "design_specification.yaml").write_text("entities: []\n", encoding="utf-8")
    (tmp_path / "project_settings.json").write_text(json.dumps({
        "num_designs": 10, "protocol": "protein-small_molecule"
    }), encoding="utf-8")
    status = inspect_project(tmp_path)
    assert not status.started
    assert status.output == tmp_path / "boltzgen_output"
    assert [stage.name for stage in status.stages] == [
        "design", "inverse_folding", "design_folding", "folding", "affinity", "analysis", "filtering"
    ]


def test_analysis_and_filtering_count_csv_rows(tmp_path: Path):
    output = tmp_path / "boltzgen_output"
    analysis = output / "intermediate_designs_inverse_folded"
    final = output / "final_ranked_designs"
    analysis.mkdir(parents=True)
    final.mkdir()
    (output / "steps.yaml").write_text("steps: [analysis, filtering]\n", encoding="utf-8")
    (tmp_path / "project_settings.json").write_text(json.dumps({"num_designs": 2}), encoding="utf-8")
    rows = "id,score\na,1\nb,2\n"
    (analysis / "aggregate_metrics_analyze.csv").write_text(rows, encoding="utf-8")
    (final / "all_designs_metrics.csv").write_text(rows, encoding="utf-8")
    (final / "results_overview.pdf").write_bytes(b"%PDF-test")
    status = inspect_project(tmp_path)
    assert status.complete
    assert all("2/2" in stage.evidence for stage in status.stages)
