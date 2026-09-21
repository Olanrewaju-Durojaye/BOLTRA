from pathlib import Path

import yaml

from boltzgen_workbench.scaffolds import (
    copy_scaffolds,
    custom_scaffold_definition,
    discover_scaffolds,
    scaffold_setup_guidance,
)


def test_discovers_and_snapshots_complete_scaffolds(tmp_path: Path):
    source = tmp_path / "example" / "nanobody_scaffolds"
    source.mkdir(parents=True)
    (source / "7eow.cif").write_text("data_7eow\n", encoding="utf-8")
    (source / "7eow.yaml").write_text("path: 7eow.cif\ninclude: all\n", encoding="utf-8")
    scaffolds = discover_scaffolds("nanobody", tmp_path)
    assert [item.name for item in scaffolds] == ["7eow"]
    project = tmp_path / "project"
    (project / "inputs").mkdir(parents=True)
    references = copy_scaffolds(scaffolds, project)
    assert references == ["inputs/scaffolds/7eow.yaml"]
    assert (project / "inputs" / "scaffolds" / "7eow.cif").is_file()


def test_custom_antibody_scaffold_has_two_design_masks():
    scaffold = custom_scaffold_definition("framework.cif", {"H": "26..32,52..57", "L": "24..34"})
    assert scaffold["include"] == [{"chain": {"id": "H"}}, {"chain": {"id": "L"}}]
    assert scaffold["design"][1]["chain"]["res_index"] == "24..34"
    assert len(scaffold["reset_res_index"]) == 2


def test_scaffold_guidance_mentions_environment_and_tutorial():
    guidance = scaffold_setup_guidance()
    assert "BOLTZGEN_ROOT" in guidance
    assert "TUTORIAL.md" in guidance
