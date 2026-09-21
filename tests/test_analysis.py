from pathlib import Path

from boltzgen_workbench.analysis import (
    KNOWN_METRICS,
    _label_indices,
    analyze,
    available_metrics,
    compact_design_label,
    load_metrics,
)


FIXTURE = Path(__file__).parent / "fixtures" / "cyclotide_test"


def test_compact_numeric_design_labels():
    assert compact_design_label("design_specification_9", 1) == "9"
    assert compact_design_label("rank1_design_specification_27.cif", 1) == "27"
    assert compact_design_label("binder_alpha", 6) == "6"


def test_large_plot_label_selection_is_capped():
    import pandas as pd

    frame = pd.DataFrame({
        "ptm": [value / 100 for value in range(100)],
        "selection_category": ["neither"] * 40 + ["native_only"] * 20
        + ["custom_only"] * 20 + ["both"] * 20,
    })
    indexes = _label_indices(frame, KNOWN_METRICS["ptm"])
    assert len(indexes) == 30
    assert all(frame.loc[index, "selection_category"] in {"both", "custom_only"} for index in indexes)


def test_redesign_specific_metric_is_recognized():
    import pandas as pd

    names = {metric.name for metric in available_metrics(pd.DataFrame({"design_residue_iptm": [0.8]}))}
    assert "design_residue_iptm" in names


def test_real_boltzgen_fixture_is_readable(tmp_path):
    frame, source = load_metrics(FIXTURE)
    assert source.name == "all_designs_metrics.csv"
    assert len(frame) == 2
    names = {metric.name for metric in available_metrics(frame)}
    assert {"design_ptm", "complex_pde", "design_to_target_ipsae"} <= names
    summary = analyze(
        FIXTURE,
        tmp_path / "analysis",
        KNOWN_METRICS["design_ptm"],
        [(KNOWN_METRICS["complex_pde"], 1.4)],
        0.7,
    )
    assert summary["total_designs"] == 2
    assert summary["passed_native_boltzgen_filters"] == 0
    assert summary["warning"] is not None
    assert (tmp_path / "analysis" / "design_ptm_vs_complex_pde.png").is_file()
    assert (tmp_path / "analysis" / "design_ptm_vs_complex_pde.svg").is_file()
    assert (tmp_path / "analysis" / "design_ptm_vs_complex_pde.pdf").is_file()
    assert (tmp_path / "analysis" / "scientific_report.md").is_file()
    assert (tmp_path / "analysis" / "scientific_report.html").is_file()
    assert summary["thresholds"][0]["statistics"]["count"] == 2
