from pathlib import Path

import pytest

from boltzgen_workbench.structure import (
    convert_author_ranges,
    read_structure,
    select_boltzgen_ranges,
    write_mapping_csv,
)


def pdb_atom(serial: int, name: str, residue: str, chain: str, number: int, insertion: str = "") -> str:
    return (
        f"ATOM  {serial:5d}  {name:<3s} {residue:>3s} {chain:1s}{number:4d}{insertion:1s}   "
        "  0.000   0.000   0.000  1.00 20.00           C\n"
    )


def test_pdb_author_mapping_handles_gaps_and_insertions(tmp_path: Path):
    target = tmp_path / "target.pdb"
    target.write_text(
        pdb_atom(1, "CA", "ALA", "A", 97)
        + pdb_atom(2, "CA", "GLY", "A", 99)
        + pdb_atom(3, "CA", "SER", "A", 100, "A")
        + pdb_atom(4, "CA", "THR", "A", 100, "B")
        + pdb_atom(5, "CA", "TYR", "B", 7),
        encoding="utf-8",
    )
    chains = read_structure(target)
    assert list(chains) == ["A", "B"]
    assert [r.author_id for r in chains["A"]] == ["97", "99", "100A", "100B"]
    converted, selected = convert_author_ranges("97-100B", chains["A"])
    assert converted == "1..4"
    assert [r.boltzgen_index for r in selected] == [1, 2, 3, 4]


def test_boltzgen_range_validation(tmp_path: Path):
    target = tmp_path / "target.pdb"
    target.write_text(pdb_atom(1, "CA", "ALA", "A", 62), encoding="utf-8")
    residues = read_structure(target)["A"]
    with pytest.raises(ValueError, match="outside 1-1"):
        select_boltzgen_ranges("2", residues)


def test_mapping_csv_records_selected_residues(tmp_path: Path):
    target = tmp_path / "target.pdb"
    target.write_text(
        pdb_atom(1, "CA", "ALA", "A", 62) + pdb_atom(2, "CA", "GLY", "A", 64),
        encoding="utf-8",
    )
    residues = read_structure(target)["A"]
    output = tmp_path / "mapping.csv"
    write_mapping_csv(output, residues, {2})
    text = output.read_text(encoding="utf-8")
    assert "selected_for_design_or_binding" in text
    assert "A,1,62,ALA,,False" in text
    assert "A,2,64,GLY,,True" in text
