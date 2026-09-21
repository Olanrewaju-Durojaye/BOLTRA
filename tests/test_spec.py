from boltzgen_workbench.spec import (
    cyclotide_pattern,
    cyclotide_spec,
    designed_length,
    normalize_ranges,
    protein_binder_spec,
    protein_redesign_spec,
    scaffold_binder_spec,
    small_molecule_binder_spec,
    validate_ligand_structure_file,
)


def test_kalata_b1_pattern():
    pattern = cyclotide_pattern([3, 4, 4, 1, 4, 7])
    assert pattern.sequence == "C3C4C4C1C4C7"
    assert pattern.cysteine_positions == (1, 5, 10, 15, 17, 22)
    assert pattern.bonds == ((1, 15), (5, 17), (10, 22))
    assert pattern.total_length == 29


def test_compact_sequence_length():
    assert designed_length("C3C4C4C1C4C7") == (29, 29)
    assert designed_length("3..5C6C3") == (14, 16)


def test_ranges_are_normalized():
    assert normalize_ranges("90-94, 97, 102-108") == "90..94,97,102..108"


def test_cyclotide_yaml_bonds():
    spec, _ = cyclotide_spec("target.cif", "A", "90-94", [3, 4, 4, 1, 4, 7])
    assert spec["constraints"][0]["bond"]["atom2"] == ["B", 15, "SG"]
    target = spec["entities"][1]["file"]
    assert target["structure_groups"] == "all"
    assert target["binding_types"][0]["chain"]["binding"] == "90..94"


def test_protein_binder_specification():
    spec = protein_binder_spec("target.pdb", "A", "7-12", 100, 150)
    assert spec["entities"][0]["protein"] == {"id": "B", "sequence": "100..150"}
    assert spec["entities"][1]["file"]["binding_types"][0]["chain"]["binding"] == "7..12"


def test_protein_redesign_specification():
    spec = protein_redesign_spec("target.cif", "C", "1-20,25")
    entity = spec["entities"][0]["file"]
    assert entity["include"] == [{"chain": {"id": "C"}}]
    assert entity["design"] == [{"chain": {"id": "C", "res_index": "1..20,25"}}]


def test_small_molecule_ccd_specification():
    spec = small_molecule_binder_spec(140, 180, "ccd", "tsa")
    assert spec == {"entities": [
        {"protein": {"id": "A", "sequence": "140..180"}},
        {"ligand": {"id": "B", "ccd": "TSA"}},
    ]}


def test_small_molecule_smiles_specification():
    spec = small_molecule_binder_spec(150, 200, "smiles", "CC(=O)O")
    assert spec["entities"][1]["ligand"]["smiles"] == "CC(=O)O"


def test_small_molecule_file_specification():
    spec = small_molecule_binder_spec(150, 150, "file", "ligand.pdb")
    assert spec["entities"][1]["file"]["include"] == "all"
    assert spec["entities"][1]["file"]["path"] == "inputs/ligand.pdb"


def test_scaffold_binder_specification():
    spec = scaffold_binder_spec("target.cif", "A", "7-12", [
        "inputs/scaffolds/7eow.yaml", "inputs/scaffolds/7xl0.yaml"
    ])
    assert spec["entities"][0]["file"]["binding_types"][0]["chain"]["binding"] == "7..12"
    assert spec["entities"][1]["file"]["path"] == [
        "inputs/scaffolds/7eow.yaml", "inputs/scaffolds/7xl0.yaml"
    ]


def test_bare_ligand_pdb_is_rejected(tmp_path):
    ligand = tmp_path / "ligand.pdb"
    ligand.write_text(
        "HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00 20.00           C\n",
        encoding="utf-8",
    )
    try:
        validate_ligand_structure_file(ligand)
    except ValueError as exc:
        assert "HETATM-only" in str(exc)
    else:
        raise AssertionError("bare ligand-only PDB was accepted")


def test_complete_pdb_is_accepted_as_ligand_file(tmp_path):
    structure = tmp_path / "complex.pdb"
    structure.write_text(
        "ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00 20.00           C\n"
        "HETATM    2  C1  LIG B   1       1.000   0.000   0.000  1.00 20.00           C\n",
        encoding="utf-8",
    )
    validate_ligand_structure_file(structure)
