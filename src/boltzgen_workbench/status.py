from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import yaml


CANONICAL_STEPS = ["design", "inverse_folding", "design_folding", "folding", "affinity", "analysis", "filtering"]


@dataclass(frozen=True)
class StageStatus:
    name: str
    complete: bool
    evidence: str


@dataclass(frozen=True)
class ProjectStatus:
    project: Path
    output: Path
    expected_designs: int | None
    stages: list[StageStatus]
    problems: list[str]
    started: bool

    @property
    def complete(self) -> bool:
        return bool(self.stages) and all(stage.complete for stage in self.stages)

    @property
    def remaining_steps(self) -> list[str]:
        first = next((i for i, stage in enumerate(self.stages) if not stage.complete), len(self.stages))
        return [stage.name for stage in self.stages[first:]]


def resolve_project(path: Path) -> tuple[Path, Path]:
    resolved = path.expanduser().resolve()
    if (resolved / "boltzgen_output" / "steps.yaml").is_file():
        return resolved, resolved / "boltzgen_output"
    if (resolved / "steps.yaml").is_file():
        return resolved.parent, resolved
    if (resolved / "project_settings.json").is_file() and (resolved / "design_specification.yaml").is_file():
        return resolved, resolved / "boltzgen_output"
    raise FileNotFoundError(
        "Choose a BOLTRA project directory, or a BoltzGen output directory containing steps.yaml"
    )


def _configured_steps(output: Path) -> list[str]:
    try:
        data = yaml.safe_load((output / "steps.yaml").read_text(encoding="utf-8")) or {}
    except (OSError, yaml.YAMLError):
        return CANONICAL_STEPS.copy()
    text = json.dumps(data).lower()
    steps = [step for step in CANONICAL_STEPS if step in text]
    return steps or CANONICAL_STEPS.copy()


def _expected(project: Path) -> int | None:
    settings = project / "project_settings.json"
    if settings.is_file():
        try:
            value = json.loads(settings.read_text(encoding="utf-8")).get("num_designs")
            return int(value) if value is not None else None
        except (ValueError, OSError, json.JSONDecodeError):
            return None
    return None


def _count_cif(directory: Path) -> int:
    return len(list(directory.glob("*.cif"))) if directory.is_dir() else 0


def _csv_rows(path: Path) -> int:
    if not path.is_file():
        return 0
    try:
        return len(pd.read_csv(path))
    except Exception:
        return 0


def _count_affinity(output: Path) -> int:
    files = list(output.rglob("*.npz")) + list(output.rglob("*.csv"))
    return len([path for path in files if "affinity" in str(path.relative_to(output)).lower()])


def inspect_project(path: Path) -> ProjectStatus:
    project, output = resolve_project(path)
    started = (output / "steps.yaml").is_file()
    expected = _expected(project)
    counts = {
        "design": _count_cif(output / "intermediate_designs"),
        "inverse_folding": _count_cif(output / "intermediate_designs_inverse_folded"),
        "design_folding": _count_cif(output / "intermediate_designs_inverse_folded" / "refold_design_cif"),
        "folding": _count_cif(output / "intermediate_designs_inverse_folded" / "refold_cif"),
        "affinity": _count_affinity(output),
        "analysis": _csv_rows(output / "intermediate_designs_inverse_folded" / "aggregate_metrics_analyze.csv"),
        "filtering": _csv_rows(output / "final_ranked_designs" / "all_designs_metrics.csv"),
    }
    stages: list[StageStatus] = []
    if started:
        configured = _configured_steps(output)
    else:
        protocol = ""
        try:
            protocol = json.loads((project / "project_settings.json").read_text(encoding="utf-8")).get("protocol", "")
        except (OSError, json.JSONDecodeError):
            pass
        configured = ["design", "inverse_folding"]
        if protocol in {"protein-anything", "protein-small_molecule"}:
            configured.append("design_folding")
        configured.append("folding")
        if protocol == "protein-small_molecule":
            configured.append("affinity")
        configured.extend(["analysis", "filtering"])
    for step in configured:
        count = counts[step]
        if expected is None:
            complete = count > 0
            target = "at least one artifact"
        else:
            complete = count >= expected
            target = str(expected)
        if step == "filtering":
            overview = (output / "final_ranked_designs" / "results_overview.pdf").is_file()
            complete = complete and overview
            evidence = f"{count}/{target} metric rows; overview PDF: {'yes' if overview else 'no'}"
        else:
            evidence = f"{count}/{target} output artifacts"
        stages.append(StageStatus(step, complete, evidence))

    problems = []
    for candidate in output.rglob("*"):
        if candidate.is_file() and candidate.suffix.lower() in {".cif", ".npz", ".csv", ".yaml"}:
            try:
                if candidate.stat().st_size == 0:
                    problems.append(f"Empty file: {candidate.relative_to(output)}")
            except OSError:
                problems.append(f"Unreadable file: {candidate.relative_to(output)}")
    return ProjectStatus(project, output, expected, stages, problems, started)


def status_text(status: ProjectStatus) -> str:
    lines = [
        f"Project: {status.project}",
        f"BoltzGen output: {status.output}",
        f"Requested designs: {status.expected_designs if status.expected_designs is not None else 'unknown'}",
        f"Run started: {'yes' if status.started else 'no'}",
        "",
        "Pipeline status:",
    ]
    for stage in status.stages:
        marker = "COMPLETE" if stage.complete else "INCOMPLETE"
        lines.append(f"  {stage.name:16} {marker:10} {stage.evidence}")
    if status.problems:
        lines.extend(["", "Potentially damaged artifacts:", *[f"  - {item}" for item in status.problems]])
    if status.complete:
        lines.extend(["", "The configured pipeline appears complete."])
    elif not status.started:
        lines.extend(["", "This validated project is saved but has not been started."])
    else:
        lines.extend(["", "Suggested resume stages: " + " ".join(status.remaining_steps)])
    return "\n".join(lines)
