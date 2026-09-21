from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import yaml


_RANGE_RE = re.compile(r"^\d+(?:\.\.\d+)?(?:,\d+(?:\.\.\d+)?)*$")


def normalize_ranges(value: str) -> str:
    """Accept human-friendly 5-9,12 and emit BoltzGen 5..9,12 notation."""
    cleaned = re.sub(r"\s+", "", value).replace("-", "..")
    if not cleaned or not _RANGE_RE.fullmatch(cleaned):
        raise ValueError("Use comma-separated residues/ranges, e.g. 90-94,97,102-108")
    for part in cleaned.split(","):
        if ".." in part:
            start, end = map(int, part.split(".."))
            if start > end:
                raise ValueError(f"Range starts after it ends: {part}")
        elif int(part) < 1:
            raise ValueError("Residue indices must start at 1")
    return cleaned


def designed_length(sequence_spec: str) -> tuple[int, int]:
    """Return minimum and maximum length for a BoltzGen compact sequence."""
    tokens = re.findall(r"\d+\.\.\d+|\d+|[A-Z]", sequence_spec.upper())
    if "".join(tokens) != sequence_spec.upper().replace(" ", ""):
        raise ValueError("Invalid compact sequence specification")
    minimum = maximum = 0
    for token in tokens:
        if token.isalpha():
            minimum += 1
            maximum += 1
        elif ".." in token:
            lo, hi = map(int, token.split(".."))
            if lo > hi:
                raise ValueError(f"Invalid variable length: {token}")
            minimum += lo
            maximum += hi
        else:
            minimum += int(token)
            maximum += int(token)
    return minimum, maximum


@dataclass(frozen=True)
class CyclotidePattern:
    sequence: str
    cysteine_positions: tuple[int, ...]
    bonds: tuple[tuple[int, int], ...]
    total_length: int


def cyclotide_pattern(loop_lengths: Iterable[int]) -> CyclotidePattern:
    """Build a six-Cys cyclic cystine-knot pattern.

    Loop lengths are I-II, II-III, III-IV, IV-V, V-VI and VI-I. The final
    loop lies after Cys VI and closes back to Cys I through cyclic=True.
    """
    loops = tuple(int(x) for x in loop_lengths)
    if len(loops) != 6 or any(x < 0 for x in loops):
        raise ValueError("Exactly six non-negative loop lengths are required")
    sequence_parts: list[str] = []
    positions: list[int] = []
    position = 0
    for index in range(6):
        sequence_parts.append("C")
        position += 1
        positions.append(position)
        if loops[index]:
            sequence_parts.append(str(loops[index]))
            position += loops[index]
    bonds = ((positions[0], positions[3]), (positions[1], positions[4]), (positions[2], positions[5]))
    return CyclotidePattern("".join(sequence_parts), tuple(positions), bonds, position)


def target_entity(target_filename: str, chain: str, binding_site: str | None) -> dict:
    file_entity: dict = {
        "path": f"inputs/{target_filename}",
        "include": [{"chain": {"id": chain}}],
        "structure_groups": "all",
    }
    if binding_site:
        file_entity["binding_types"] = [
            {"chain": {"id": chain, "binding": normalize_ranges(binding_site)}}
        ]
    return {"file": file_entity}


def scaffold_binder_spec(
    target_filename: str,
    target_chain: str,
    binding_site: str | None,
    scaffold_references: list[str],
) -> dict:
    """Build a nanobody/antibody specification from one or more scaffold YAMLs."""
    if not scaffold_references:
        raise ValueError("At least one scaffold is required")
    return {
        "entities": [
            target_entity(target_filename, target_chain, binding_site),
            {"file": {"path": list(scaffold_references)}},
        ]
    }


def peptide_spec(
    target_filename: str,
    target_chain: str,
    binding_site: str | None,
    minimum_length: int,
    maximum_length: int,
    cyclic: bool = False,
) -> dict:
    if minimum_length < 1 or maximum_length < minimum_length:
        raise ValueError("Peptide length range is invalid")
    sequence = str(minimum_length) if minimum_length == maximum_length else f"{minimum_length}..{maximum_length}"
    protein: dict = {"id": "B", "sequence": sequence}
    if cyclic:
        protein["cyclic"] = True
    return {"entities": [{"protein": protein}, target_entity(target_filename, target_chain, binding_site)]}


def protein_binder_spec(
    target_filename: str,
    target_chain: str,
    binding_site: str | None,
    minimum_length: int,
    maximum_length: int,
) -> dict:
    """Build a de novo protein-binder specification for protein-anything."""
    if minimum_length < 1 or maximum_length < minimum_length:
        raise ValueError("Protein length range is invalid")
    sequence = str(minimum_length) if minimum_length == maximum_length else f"{minimum_length}..{maximum_length}"
    return {
        "entities": [
            {"protein": {"id": "B", "sequence": sequence}},
            target_entity(target_filename, target_chain, binding_site),
        ]
    }


def small_molecule_binder_spec(
    minimum_length: int,
    maximum_length: int,
    ligand_format: str,
    ligand_value: str,
) -> dict:
    """Build a protein-small_molecule specification.

    ``ligand_value`` is a CCD identifier or SMILES string for those formats,
    and the copied filename for the file format.
    """
    if minimum_length < 1 or maximum_length < minimum_length:
        raise ValueError("Protein length range is invalid")
    sequence = str(minimum_length) if minimum_length == maximum_length else f"{minimum_length}..{maximum_length}"
    protein = {"protein": {"id": "A", "sequence": sequence}}
    format_name = ligand_format.strip().lower()
    if format_name == "ccd":
        ccd = ligand_value.strip().upper()
        if not re.fullmatch(r"[A-Z0-9]{1,5}", ccd):
            raise ValueError("CCD code must contain 1-5 letters or numbers")
        target = {"ligand": {"id": "B", "ccd": ccd}}
    elif format_name == "smiles":
        smiles = ligand_value.strip()
        if not smiles or any(character in smiles for character in "\r\n"):
            raise ValueError("SMILES must be a non-empty single line")
        target = {"ligand": {"id": "B", "smiles": smiles}}
    elif format_name == "file":
        filename = Path(ligand_value).name
        if Path(filename).suffix.lower() not in {".pdb", ".cif", ".mmcif"}:
            raise ValueError("Ligand structure must be PDB or mmCIF")
        target = {"file": {"path": f"inputs/{filename}", "include": "all", "structure_groups": "all"}}
    else:
        raise ValueError("Ligand format must be ccd, smiles, or file")
    return {"entities": [protein, target]}


def validate_ligand_structure_file(path: Path) -> None:
    """Reject bare ligand PDB files that BoltzGen's file parser cannot read.

    A HETATM-only PDB does not contain the polymer/entity information expected by
    BoltzGen 0.3.x. CCD or SMILES is the reproducible route for such ligands.
    Complete PDB/mmCIF structures remain accepted and are ultimately validated by
    ``boltzgen check``.
    """
    suffix = path.suffix.lower()
    if suffix not in {".pdb", ".cif", ".mmcif"}:
        raise ValueError("Ligand structure must be a PDB or mmCIF file")
    if suffix != ".pdb":
        return
    try:
        records = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError as exc:
        raise ValueError(f"Could not read ligand structure: {path}") from exc
    has_atom = any(line.startswith("ATOM  ") for line in records)
    has_hetatm = any(line.startswith("HETATM") for line in records)
    if has_hetatm and not has_atom:
        raise ValueError(
            "This appears to be a bare HETATM-only ligand PDB. BoltzGen 0.3.x may "
            "fail to parse it because polymer/entity records are absent. Use the "
            "ligand's CCD code or SMILES instead, or provide a complete structure "
            "file that passes 'boltzgen check'."
        )
    if not has_atom and not has_hetatm:
        raise ValueError("No ATOM or HETATM coordinates were found in the ligand PDB")


def protein_redesign_spec(
    target_filename: str,
    target_chain: str,
    redesign_residues: str,
) -> dict:
    """Build a structure-guided redesign specification using a design mask."""
    if not redesign_residues:
        raise ValueError("At least one residue must be selected for redesign")
    return {
        "entities": [
            {
                "file": {
                    "path": f"inputs/{target_filename}",
                    "include": [{"chain": {"id": target_chain}}],
                    "structure_groups": "all",
                    "design": [
                        {"chain": {"id": target_chain, "res_index": normalize_ranges(redesign_residues)}}
                    ],
                }
            }
        ]
    }


def cyclotide_spec(
    target_filename: str,
    target_chain: str,
    binding_site: str | None,
    loop_lengths: Iterable[int],
) -> tuple[dict, CyclotidePattern]:
    pattern = cyclotide_pattern(loop_lengths)
    entities = [
        {"protein": {"id": "B", "sequence": pattern.sequence, "cyclic": True}},
        target_entity(target_filename, target_chain, binding_site),
    ]
    constraints = [
        {"bond": {"atom1": ["B", first, "SG"], "atom2": ["B", second, "SG"]}}
        for first, second in pattern.bonds
    ]
    return {"entities": entities, "constraints": constraints}, pattern


def write_spec(specification: dict, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(yaml.safe_dump(specification, sort_keys=False), encoding="utf-8")
