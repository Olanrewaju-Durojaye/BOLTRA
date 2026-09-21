from __future__ import annotations

import importlib.metadata
import os
import shutil
from dataclasses import dataclass
from pathlib import Path

import yaml


@dataclass(frozen=True)
class Scaffold:
    name: str
    yaml_path: Path
    structure_paths: tuple[Path, ...]


def _candidate_roots() -> list[Path]:
    roots: list[Path] = []
    configured = os.environ.get("BOLTZGEN_ROOT")
    if configured:
        roots.append(Path(configured).expanduser())
    roots.extend([Path.cwd(), *Path.cwd().parents, Path.home() / "boltzgen-main", Path.home() / "boltzgen"])
    try:
        distribution = importlib.metadata.distribution("boltzgen")
        roots.extend([Path(distribution.locate_file("")), Path(distribution.locate_file("")).parent])
    except importlib.metadata.PackageNotFoundError:
        pass
    unique: list[Path] = []
    for root in roots:
        resolved = root.expanduser().resolve()
        if resolved not in unique:
            unique.append(resolved)
    return unique


def discover_scaffolds(kind: str, extra_root: Path | None = None) -> list[Scaffold]:
    """Discover scaffold YAMLs from the installed package or a BoltzGen checkout."""
    directory_name = {"nanobody": "nanobody_scaffolds", "antibody": "fab_scaffolds"}.get(kind)
    if directory_name is None:
        raise ValueError("Scaffold kind must be nanobody or antibody")
    directories: list[Path] = []
    roots = ([extra_root.expanduser().resolve()] if extra_root is not None else []) + _candidate_roots()
    for root in roots:
        direct = root if root.name == directory_name else root / directory_name
        for candidate in (root / "example" / directory_name, direct,
                          root / "boltzgen" / "example" / directory_name):
            if candidate.is_dir() and candidate not in directories:
                directories.append(candidate)
    # Wheels may place examples or resources below package-specific directories
    # that cannot be inferred reliably from the distribution root alone.
    try:
        distribution = importlib.metadata.distribution("boltzgen")
        for entry in distribution.files or []:
            if directory_name not in entry.parts:
                continue
            located = Path(distribution.locate_file(entry)).resolve()
            candidate = located if located.is_dir() else located.parent
            while candidate.name != directory_name and candidate != candidate.parent:
                candidate = candidate.parent
            if candidate.name == directory_name and candidate.is_dir() and candidate not in directories:
                directories.append(candidate)
    except importlib.metadata.PackageNotFoundError:
        pass
    discovered: dict[str, Scaffold] = {}
    for directory in directories:
        for yaml_path in sorted([*directory.glob("*.yaml"), *directory.glob("*.yml")]):
            try:
                data = yaml.safe_load(yaml_path.read_text(encoding="utf-8")) or {}
                raw_paths = data.get("path")
                raw_paths = raw_paths if isinstance(raw_paths, list) else [raw_paths]
                structures = tuple((yaml_path.parent / value).resolve() for value in raw_paths if isinstance(value, str))
            except (OSError, yaml.YAMLError):
                continue
            if structures and all(path.is_file() for path in structures):
                discovered.setdefault(yaml_path.stem, Scaffold(yaml_path.stem, yaml_path.resolve(), structures))
    return sorted(discovered.values(), key=lambda item: item.name.lower())


def copy_scaffolds(scaffolds: list[Scaffold], project: Path) -> list[str]:
    """Snapshot selected scaffold YAMLs and referenced structures into a project."""
    destination = project / "inputs" / "scaffolds"
    destination.mkdir(exist_ok=True)
    references: list[str] = []
    for scaffold in scaffolds:
        yaml_destination = destination / scaffold.yaml_path.name
        shutil.copy2(scaffold.yaml_path, yaml_destination)
        for structure in scaffold.structure_paths:
            shutil.copy2(structure, destination / structure.name)
        references.append(f"inputs/scaffolds/{yaml_destination.name}")
    return references


def custom_scaffold_definition(structure_filename: str, chains: dict[str, str]) -> dict:
    """Create a fixed-length, structure-conditioned CDR redesign scaffold."""
    include = [{"chain": {"id": chain}} for chain in chains]
    design = [{"chain": {"id": chain, "res_index": residues}} for chain, residues in chains.items()]
    groups = [{"group": {"id": chain, "visibility": 2}} for chain in chains]
    groups.extend({"group": {"id": chain, "res_index": residues, "visibility": 0}}
                  for chain, residues in chains.items())
    return {
        "path": structure_filename,
        "include": include,
        "design": design,
        "structure_groups": groups,
        "reset_res_index": [{"chain": {"id": chain}} for chain in chains],
    }


def scaffold_setup_guidance() -> str:
    """Return actionable setup guidance when BoltzGen examples are not installed."""
    return (
        "Some BoltzGen package installations omit example/nanobody_scaffolds and "
        "example/fab_scaffolds. Download the official BoltzGen source checkout, "
        "then set BOLTZGEN_ROOT to its directory. Detailed commands are provided "
        "in TUTORIAL.md under 'Nanobody and antibody scaffold resources'."
    )
