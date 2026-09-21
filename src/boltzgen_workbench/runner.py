from __future__ import annotations

import shutil
import subprocess
import json
import os
import platform
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path


GIB = 1024 ** 3


@dataclass(frozen=True)
class ResourceSnapshot:
    total_ram_gib: float | None
    available_ram_gib: float | None
    total_swap_gib: float | None
    free_swap_gib: float | None
    free_disk_gib: float


def _meminfo() -> dict[str, int]:
    values: dict[str, int] = {}
    try:
        for line in Path("/proc/meminfo").read_text(encoding="utf-8").splitlines():
            key, raw = line.split(":", 1)
            values[key] = int(raw.strip().split()[0]) * 1024
    except (OSError, ValueError, IndexError):
        pass
    return values


def system_resources(path: Path) -> ResourceSnapshot:
    memory = _meminfo()
    usage = shutil.disk_usage(path)
    gib = lambda value: round(value / GIB, 1) if value is not None else None
    return ResourceSnapshot(
        gib(memory.get("MemTotal")),
        gib(memory.get("MemAvailable")),
        gib(memory.get("SwapTotal")),
        gib(memory.get("SwapFree")),
        round(usage.free / GIB, 1),
    )


def safe_worker_count() -> int:
    raw = os.environ.get("BOLTRA_NUM_WORKERS", "0")
    try:
        workers = int(raw)
    except ValueError as exc:
        raise ValueError("BOLTRA_NUM_WORKERS must be a non-negative whole number") from exc
    if workers < 0:
        raise ValueError("BOLTRA_NUM_WORKERS must be zero or greater")
    return workers


def resource_warnings(snapshot: ResourceSnapshot, num_designs: int, workers: int) -> list[str]:
    warnings: list[str] = []
    if snapshot.total_swap_gib is not None and snapshot.total_swap_gib < 8:
        warnings.append(
            f"Only {snapshot.total_swap_gib:g} GiB swap is configured; 16-32 GiB is recommended "
            "for BoltzGen analysis on workstations with limited RAM."
        )
    if snapshot.available_ram_gib is not None and snapshot.available_ram_gib < 12:
        warnings.append(
            f"Only {snapshot.available_ram_gib:g} GiB RAM is currently available; close browsers "
            "and other memory-heavy applications before running."
        )
    if workers > 0 and snapshot.total_ram_gib is not None and snapshot.total_ram_gib < 64:
        warnings.append(
            f"{workers} workers were requested on a {snapshot.total_ram_gib:g} GiB RAM system; "
            "use BOLTRA_NUM_WORKERS=0 to reduce out-of-memory risk."
        )
    if num_designs >= 100:
        warnings.append(
            f"This run requests {num_designs} designs. Use staged pilot sizes and verify free disk "
            "space and recovery before scaling further."
        )
    if snapshot.free_disk_gib < 25:
        warnings.append(f"Only {snapshot.free_disk_gib:g} GiB free disk space remains.")
    return warnings


def require_boltzgen() -> str:
    executable = shutil.which("boltzgen")
    if not executable:
        raise RuntimeError("BoltzGen was not found. Activate the BoltzGen Conda environment first.")
    return executable


def boltzgen_version() -> str:
    executable = require_boltzgen()
    result = subprocess.run([executable, "--version"], capture_output=True, text=True)
    return (result.stdout or result.stderr).strip() or "unknown"


def gpu_information() -> str | None:
    executable = shutil.which("nvidia-smi")
    if not executable:
        return None
    result = subprocess.run(
        [executable, "--query-gpu=name,driver_version,memory.total", "--format=csv,noheader"],
        capture_output=True,
        text=True,
    )
    return result.stdout.strip() if result.returncode == 0 and result.stdout.strip() else None


def _run_logged(command: list[str], cwd: Path, action: str, project: Path) -> int:
    from . import __version__

    started = datetime.now(timezone.utc)
    record = {
        "action": action,
        "command": command,
        "cwd": str(cwd.resolve()),
        "started_utc": started.isoformat(),
        "boltra_version": __version__,
        "boltzgen_version": boltzgen_version(),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "gpu": gpu_information(),
        "resources": asdict(system_resources(project)),
    }
    interrupted = False
    try:
        code = subprocess.run(command, cwd=cwd).returncode
    except KeyboardInterrupt:
        code = 130
        interrupted = True
    record.update({
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "exit_status": code,
        "status": "completed" if code == 0 else "failed_or_interrupted",
    })
    log = project / "boltra_execution_log.jsonl"
    with log.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record, sort_keys=True) + "\n")
    if interrupted:
        raise KeyboardInterrupt
    return code


def validate(specification: Path, cwd: Path) -> int:
    command = [require_boltzgen(), "check", str(specification)]
    return _run_logged(command, cwd, "validate", cwd)


def run_pipeline(
    specification: Path,
    output: Path,
    protocol: str,
    num_designs: int,
    budget: int,
    cwd: Path,
) -> int:
    if budget < 1 or budget > num_designs:
        raise ValueError("Budget must be between 1 and the number of designs")
    workers = safe_worker_count()
    command = [
        require_boltzgen(), "run", str(specification), "--output", str(output),
        "--protocol", protocol, "--num_designs", str(num_designs),
        "--budget", str(budget), "--devices", "1", "--num_workers", str(workers),
    ]
    return _run_logged(command, cwd, "new_run", cwd)


def resume_pipeline(output: Path, steps: list[str], project: Path) -> int:
    if not (output / "steps.yaml").is_file():
        raise FileNotFoundError(f"BoltzGen steps.yaml not found in {output}")
    settings_path = project / "project_settings.json"
    specification = project / "design_specification.yaml"
    if not settings_path.is_file() or not specification.is_file():
        raise FileNotFoundError(
            "Safe automatic resume requires the original BOLTRA project_settings.json "
            "and design_specification.yaml"
        )
    settings = json.loads(settings_path.read_text(encoding="utf-8"))
    num_designs = int(settings["num_designs"])
    budget = int(settings["budget"])
    protocol = settings.get("protocol", "peptide-anything")
    workers = safe_worker_count()
    command = [require_boltzgen(), "run", str(specification), "--output", str(output),
               "--protocol", protocol, "--num_designs", str(num_designs),
               "--budget", str(budget), "--devices", "1", "--num_workers", str(workers)]
    if steps:
        command.extend(["--steps", *steps])
    command.append("--reuse")
    return _run_logged(command, project, "resume", project)
