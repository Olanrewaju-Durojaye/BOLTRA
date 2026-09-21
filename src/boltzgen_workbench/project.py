from __future__ import annotations

import json
import re
import shutil
from pathlib import Path


def safe_project_name(name: str) -> str:
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", name.strip()).strip("._-")
    if not value:
        raise ValueError("Project name must contain a letter or number")
    return value


def create_project(root: Path, name: str, target: Path | None = None) -> tuple[Path, Path | None]:
    if target is not None:
        if not target.is_file():
            raise FileNotFoundError(f"Target structure not found: {target}")
        if target.suffix.lower() not in {".pdb", ".cif", ".mmcif"}:
            raise ValueError("Target must be a PDB or mmCIF structure")
    resolved_root = root.expanduser().resolve()
    if resolved_root.exists() and not resolved_root.is_dir():
        raise ValueError(f"Project location must be a directory, not a file: {resolved_root}")
    project = resolved_root / safe_project_name(name)
    if project.exists():
        raise ValueError(f"Project already exists; choose another name: {project}")
    inputs = project / "inputs"
    inputs.mkdir(parents=True, exist_ok=False)
    for child in ("configuration", "boltzgen_output", "analysis", "selected_designs", "reports"):
        (project / child).mkdir()
    copied = None
    if target is not None:
        copied = inputs / target.name
        shutil.copy2(target.expanduser().resolve(), copied)
    return project, copied


def save_metadata(project: Path, data: dict) -> None:
    (project / "project_settings.json").write_text(
        json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
