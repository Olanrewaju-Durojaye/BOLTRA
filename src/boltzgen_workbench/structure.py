from __future__ import annotations

import csv
import re
import shlex
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Residue:
    chain: str
    boltzgen_index: int
    author_number: int
    insertion_code: str
    name: str
    label_sequence_id: str = ""

    @property
    def author_id(self) -> str:
        return f"{self.author_number}{self.insertion_code}"


def read_structure(path: Path) -> dict[str, list[Residue]]:
    path = path.expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"Target structure not found: {path}")
    suffix = path.suffix.lower()
    if suffix == ".pdb":
        chains = _read_pdb(path)
    elif suffix in {".cif", ".mmcif"}:
        chains = _read_mmcif(path)
    else:
        raise ValueError("Target must be a PDB or mmCIF structure")
    if not chains:
        raise ValueError("No polymer ATOM residues were found in the target structure")
    return chains


def _read_pdb(path: Path) -> dict[str, list[Residue]]:
    raw: dict[str, list[tuple[int, str, str, str]]] = {}
    seen: dict[str, set[tuple[int, str]]] = {}
    model_seen = False
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            record = line[:6].strip()
            if record == "MODEL":
                if model_seen:
                    break
                model_seen = True
                continue
            if record == "ENDMDL":
                break
            if record != "ATOM":
                continue
            altloc = line[16:17].strip()
            if altloc not in {"", "A", "1"}:
                continue
            chain = line[21:22].strip() or "_"
            number_text = line[22:26].strip()
            if not re.fullmatch(r"-?\d+", number_text):
                continue
            number = int(number_text)
            insertion = line[26:27].strip()
            key = (number, insertion)
            seen.setdefault(chain, set())
            if key in seen[chain]:
                continue
            seen[chain].add(key)
            raw.setdefault(chain, []).append((number, insertion, line[17:20].strip(), ""))
    return _index_residues(raw)


def _read_mmcif(path: Path) -> dict[str, list[Residue]]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    raw: dict[str, list[tuple[int, str, str, str]]] = {}
    seen: dict[str, set[tuple[int, str]]] = {}
    index = 0
    while index < len(lines):
        if lines[index].strip() != "loop_":
            index += 1
            continue
        index += 1
        columns: list[str] = []
        while index < len(lines) and lines[index].lstrip().startswith("_"):
            columns.append(lines[index].strip())
            index += 1
        if not columns or not any(column.startswith("_atom_site.") for column in columns):
            continue
        positions = {column: position for position, column in enumerate(columns)}
        required = ["_atom_site.group_PDB", "_atom_site.auth_seq_id", "_atom_site.label_comp_id"]
        if any(column not in positions for column in required):
            raise ValueError("The mmCIF atom_site table lacks required residue identifiers")
        tokens: list[str] = []
        while index < len(lines):
            stripped = lines[index].strip()
            if not stripped or stripped.startswith("#"):
                index += 1
                if tokens:
                    break
                continue
            if stripped == "loop_" or stripped.startswith("_") or stripped.startswith("data_"):
                break
            tokens.extend(shlex.split(lines[index], posix=True))
            while len(tokens) >= len(columns):
                row, tokens = tokens[: len(columns)], tokens[len(columns) :]
                if row[positions["_atom_site.group_PDB"]].upper() != "ATOM":
                    continue
                model_column = positions.get("_atom_site.pdbx_PDB_model_num")
                if model_column is not None and row[model_column] not in {"1", ".", "?"}:
                    continue
                chain_column = positions.get("_atom_site.auth_asym_id", positions.get("_atom_site.label_asym_id"))
                chain = row[chain_column] if chain_column is not None else "_"
                number_text = row[positions["_atom_site.auth_seq_id"]]
                if not re.fullmatch(r"-?\d+", number_text):
                    continue
                insertion_column = positions.get("_atom_site.pdbx_PDB_ins_code")
                insertion = row[insertion_column] if insertion_column is not None else ""
                if insertion in {".", "?"}:
                    insertion = ""
                number = int(number_text)
                key = (number, insertion)
                seen.setdefault(chain, set())
                if key in seen[chain]:
                    continue
                seen[chain].add(key)
                label_column = positions.get("_atom_site.label_seq_id")
                label = row[label_column] if label_column is not None else ""
                raw.setdefault(chain, []).append(
                    (number, insertion, row[positions["_atom_site.label_comp_id"]], label)
                )
            index += 1
        break
    return _index_residues(raw)


def _index_residues(raw: dict[str, list[tuple[int, str, str, str]]]) -> dict[str, list[Residue]]:
    return {
        chain: [Residue(chain, i, number, insertion, name, label) for i, (number, insertion, name, label) in enumerate(rows, 1)]
        for chain, rows in raw.items()
    }


_AUTHOR_TOKEN = re.compile(r"^(\d+)([A-Za-z]?)$")


def convert_author_ranges(value: str, residues: list[Residue]) -> tuple[str, list[Residue]]:
    lookup = {residue.author_id.upper(): residue for residue in residues}
    selected: list[Residue] = []
    missing: list[str] = []
    for item in (part.strip() for part in value.split(",")):
        if not item:
            raise ValueError("Use comma-separated author residues/ranges, e.g. 97-110,120,150-160")
        if "-" in item:
            start_text, end_text = item.split("-", 1)
            start = _parse_author_id(start_text)
            end = _parse_author_id(end_text)
            start_key = (start[0], _insertion_rank(start[1]))
            end_key = (end[0], _insertion_rank(end[1]))
            if start_key > end_key:
                raise ValueError(f"Range starts after it ends: {item}")
            matches = [r for r in residues if start_key <= (r.author_number, _insertion_rank(r.insertion_code)) <= end_key]
            if not matches:
                missing.append(item)
            selected.extend(matches)
        else:
            number, insertion = _parse_author_id(item)
            residue = lookup.get(f"{number}{insertion}".upper())
            if residue is None:
                missing.append(item)
            else:
                selected.append(residue)
    if missing:
        raise ValueError("These author residues were not found in the selected chain: " + ", ".join(missing))
    unique = {residue.boltzgen_index: residue for residue in selected}
    ordered = [unique[index] for index in sorted(unique)]
    return compress_indices([residue.boltzgen_index for residue in ordered]), ordered


def select_boltzgen_ranges(value: str, residues: list[Residue]) -> tuple[str, list[Residue]]:
    indices: list[int] = []
    for item in (part.strip() for part in value.split(",")):
        match = re.fullmatch(r"(\d+)(?:\s*(?:-|\.\.)\s*(\d+))?", item)
        if not match:
            raise ValueError("Use comma-separated BoltzGen indices/ranges, e.g. 38-51,61,91-101")
        start = int(match.group(1))
        end = int(match.group(2) or start)
        if start < 1 or end < start or end > len(residues):
            raise ValueError(f"BoltzGen range {item} is outside 1-{len(residues)}")
        indices.extend(range(start, end + 1))
    unique = sorted(set(indices))
    return compress_indices(unique), [residues[index - 1] for index in unique]


def _parse_author_id(value: str) -> tuple[int, str]:
    match = _AUTHOR_TOKEN.fullmatch(value.strip())
    if not match:
        raise ValueError(f"Invalid PDB author residue identifier: {value}")
    return int(match.group(1)), match.group(2).upper()


def _insertion_rank(code: str) -> int:
    return 0 if not code else ord(code.upper()) - ord("A") + 1


def compress_indices(indices: list[int]) -> str:
    if not indices:
        raise ValueError("At least one residue is required")
    ranges: list[str] = []
    start = previous = indices[0]
    for index in indices[1:]:
        if index == previous + 1:
            previous = index
            continue
        ranges.append(str(start) if start == previous else f"{start}..{previous}")
        start = previous = index
    ranges.append(str(start) if start == previous else f"{start}..{previous}")
    return ",".join(ranges)


def write_mapping_csv(path: Path, residues: list[Residue], selected: set[int]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["chain", "boltzgen_index", "pdb_author_residue", "residue_name", "label_seq_id", "selected_for_design_or_binding"])
        for residue in residues:
            writer.writerow([
                residue.chain,
                residue.boltzgen_index,
                residue.author_id,
                residue.name,
                residue.label_sequence_id,
                residue.boltzgen_index in selected,
            ])
