from __future__ import annotations

import html
import json
import re
import shutil
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


@dataclass(frozen=True)
class Metric:
    name: str
    label: str
    direction: str
    description: str


KNOWN_METRICS = {
    "ptm": Metric("ptm", "Complex pTM", "higher", "Confidence in the complete predicted complex"),
    "iptm": Metric("iptm", "Complex ipTM", "higher", "Predicted interface confidence for the complex"),
    "design_ptm": Metric("design_ptm", "Design pTM", "higher", "Confidence in the designed binder fold"),
    "design_to_target_iptm": Metric("design_to_target_iptm", "Design-to-target ipTM", "higher", "Predicted interface confidence"),
    "design_iptm": Metric("design_iptm", "Design ipTM", "higher", "Predicted interface confidence involving the design"),
    "design_iiptm": Metric("design_iiptm", "Design iiPTM", "higher", "Interface-focused confidence"),
    "design_residue_iptm": Metric("design_residue_iptm", "Design-residue ipTM", "higher", "Interface confidence across redesigned residues"),
    "design_to_target_ipsae": Metric("design_to_target_ipsae", "Design-to-target ipSAE", "higher", "Interface score derived from aligned errors"),
    "design_ipsae_min": Metric("design_ipsae_min", "Minimum design ipSAE", "higher", "Minimum interface ipSAE for the design"),
    "complex_plddt": Metric("complex_plddt", "Complex pLDDT", "higher", "Local structural confidence"),
    "complex_pde": Metric("complex_pde", "Complex PDE", "lower", "Predicted distance error"),
    "complex_ipde": Metric("complex_ipde", "Complex interface PDE", "lower", "Interface predicted distance error"),
    "min_design_to_target_pae": Metric("min_design_to_target_pae", "Minimum design-to-target PAE", "lower", "Minimum aligned error across the interface"),
    "filter_rmsd": Metric("filter_rmsd", "Complex refolding RMSD", "lower", "Refolded-complex consistency"),
    "filter_rmsd_design": Metric("filter_rmsd_design", "Binder refolding RMSD", "lower", "Refolded-binder consistency"),
    "bb_rmsd": Metric("bb_rmsd", "Complex backbone RMSD", "lower", "Backbone consistency after refolding"),
    "bb_rmsd_design": Metric("bb_rmsd_design", "Binder backbone RMSD", "lower", "Designed-chain backbone consistency"),
    "plip_hbonds_refolded": Metric("plip_hbonds_refolded", "Interface hydrogen bonds", "higher", "PLIP hydrogen-bond count"),
    "plip_saltbridge_refolded": Metric("plip_saltbridge_refolded", "Interface salt bridges", "higher", "PLIP salt-bridge count"),
    "delta_sasa_refolded": Metric("delta_sasa_refolded", "Buried surface area", "higher", "Change in solvent-accessible surface area on binding"),
    "liability_score": Metric("liability_score", "Sequence liability score", "lower", "Predicted sequence-development liabilities"),
    "absolute_score": Metric("absolute_score", "BoltzGen absolute score", "higher", "BoltzGen aggregate ranking score"),
    "quality_score": Metric("quality_score", "BoltzGen quality score", "higher", "BoltzGen aggregate quality ranking"),
    "affinity_pred_value": Metric("affinity_pred_value", "Predicted log10(IC50 µM)", "lower", "Boltz-2 quantitative affinity estimate; lower indicates stronger predicted affinity"),
    "affinity_probability_binary": Metric("affinity_probability_binary", "Binding probability", "higher", "Boltz-2 binary binding probability"),
    "affinity_probability_binary1": Metric("affinity_probability_binary1", "Binding probability", "higher", "Boltz-2 binary binding probability"),
    "affinity_logits_binary": Metric("affinity_logits_binary", "Binding affinity logit", "higher", "Boltz-2 binary binding score before probability conversion"),
}


def locate_metrics(run_directory: Path) -> Path:
    for candidate in (run_directory / "final_ranked_designs" / "all_designs_metrics.csv",
                      run_directory / "intermediate_designs_inverse_folded" / "aggregate_metrics_analyze.csv"):
        if candidate.is_file():
            return candidate
    raise FileNotFoundError("No BoltzGen aggregate metrics CSV was found")


def load_metrics(run_directory: Path) -> tuple[pd.DataFrame, Path]:
    source = locate_metrics(run_directory)
    return pd.read_csv(source), source


def available_metrics(frame: pd.DataFrame) -> list[Metric]:
    return [KNOWN_METRICS[name] for name in frame.select_dtypes(include="number").columns
            if not name.startswith(("rank_", "neg_")) and not name.endswith("_z") and name in KNOWN_METRICS]


def metric_statistics(frame: pd.DataFrame, metric: Metric) -> dict[str, float | int]:
    values = pd.to_numeric(frame[metric.name], errors="coerce").dropna()
    return {"count": int(values.count()), "minimum": float(values.min()),
            "median": float(values.median()), "maximum": float(values.max())}


def apply_threshold(frame: pd.DataFrame, metric: Metric, threshold: float) -> pd.Series:
    values = pd.to_numeric(frame[metric.name], errors="coerce")
    return values.ge(threshold) if metric.direction == "higher" else values.le(threshold)


def _native_pass(frame: pd.DataFrame) -> pd.Series:
    if "pass_filters" not in frame:
        return pd.Series(False, index=frame.index)
    return frame["pass_filters"].astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def _identifier_column(frame: pd.DataFrame) -> str | None:
    return next((name for name in ("id", "design_id", "name", "file_name", "filename") if name in frame), None)


def compact_design_label(identifier: object, fallback: int) -> str:
    """Return a short numerical label while preserving full IDs in data tables."""
    stem = Path(str(identifier)).stem
    match = re.search(r"(\d+)$", stem)
    return match.group(1) if match else str(fallback)


def _label_indices(result: pd.DataFrame, central: Metric, maximum: int = 30) -> list[int]:
    """Label all small sets; prioritize selected candidates in large sets."""
    if len(result) <= 50:
        return list(result.index)
    direction = central.direction == "higher"
    ordered: list[int] = []
    for category in ("both", "custom_only", "native_only"):
        subset = result.loc[result["selection_category"].eq(category)].sort_values(
            central.name, ascending=not direction
        )
        ordered.extend(index for index in subset.index if index not in ordered)
    return ordered[:maximum]


def _find_structure(run_directory: Path, identifier: object) -> Path | None:
    stem = Path(str(identifier)).stem
    candidates = list(run_directory.rglob(f"{stem}.cif")) + list(run_directory.rglob(f"{stem}.pdb"))
    if not candidates:
        candidates = [path for path in run_directory.rglob("*.cif") if stem in path.stem]
    if not candidates:
        return None
    def preference(path: Path) -> tuple[int, int]:
        value = str(path)
        score = 0 if "final_ranked_designs" in value else 1 if "refold_cif" in value else 2
        return score, len(path.parts)
    return sorted(candidates, key=preference)[0]


def _copy_candidates(run_directory: Path, output: Path, result: pd.DataFrame,
                     custom: pd.Series, native: pd.Series) -> dict[str, int]:
    identifier = _identifier_column(result)
    counts = {"selected_custom": 0, "selected_native": 0, "selected_both": 0}
    if identifier is None:
        return counts
    groups = {"selected_custom": custom & ~native, "selected_native": native & ~custom,
              "selected_both": custom & native}
    for group, mask in groups.items():
        destination = output / "candidate_structures" / group
        destination.mkdir(parents=True, exist_ok=True)
        for value in result.loc[mask, identifier]:
            source = _find_structure(run_directory, value)
            if source:
                shutil.copy2(source, destination / source.name)
                counts[group] += 1
    return counts


def _software_version(command: list[str]) -> str:
    try:
        result = subprocess.run(command, capture_output=True, text=True, timeout=10)
        return (result.stdout or result.stderr).strip() or "unknown"
    except (OSError, subprocess.SubprocessError):
        return "unavailable"


def _project_settings(run_directory: Path) -> dict:
    path = run_directory.parent / "project_settings.json"
    try:
        return json.loads(path.read_text(encoding="utf-8")) if path.is_file() else {}
    except (OSError, json.JSONDecodeError):
        return {}


def _write_report(output: Path, summary: dict, result: pd.DataFrame, selected: pd.Series) -> None:
    settings, thresholds = summary["project_settings"], summary["thresholds"]
    lines = ["# BOLTRA post-design analysis report", "", f"Generated: {summary['generated_utc']}", "",
             "## Study and software provenance", "", f"- BOLTRA: {summary['boltra_version']}",
             f"- BoltzGen: {summary['boltzgen_version']}", f"- Metrics source: `{summary['metrics_source']}`",
             f"- Target: {settings.get('target', 'not recorded')}",
             f"- Ligand format: {settings.get('ligand_format', 'not applicable')}",
             f"- Ligand: {settings.get('ligand_display', 'not applicable')}",
             f"- Target chain: {settings.get('target_chain', 'not recorded')}",
             f"- Binding-site input: {settings.get('binding_site_input', 'not recorded')}",
             f"- BoltzGen binding-site indices: {settings.get('binding_site_boltzgen', 'not recorded')}",
             f"- Residues redesigned: {settings.get('residues_to_redesign_input', 'not applicable')}",
             f"- BoltzGen redesign indices: {settings.get('residues_to_redesign_boltzgen', 'not applicable')}",
             f"- Design type: {settings.get('kind', 'not recorded')}",
             f"- Requested designs: {settings.get('num_designs', 'not recorded')}",
             f"- Final budget: {settings.get('budget', 'not recorded')}", "", "## Filtering definition", ""]
    for item in thresholds:
        operator = ">=" if item["direction"] == "higher" else "<="
        stats = item["statistics"]
        lines.append(f"- {item['label']}: {operator} {item['threshold']}; observed min/median/max "
                     f"{stats['minimum']:.4g}/{stats['median']:.4g}/{stats['maximum']:.4g}")
    lines.extend(["", "## Outcomes", "", f"- Total evaluated: {summary['total_designs']}",
                  f"- Passed custom intersection: {summary['selected_by_custom_intersection']}",
                  f"- Passed native BoltzGen filters: {summary['passed_native_boltzgen_filters']}",
                  f"- Passed both: {summary['passed_both']}", "", "### Per-filter counts", "",
                  "| Criterion | Passing designs |", "|---|---:|"])
    for name, count in summary["individual_pass_counts"].items():
        lines.append(f"| {name} | {count} |")
    selected_frame = result.loc[selected]
    lines.extend(["", "## Custom-selected candidates", ""])
    if selected_frame.empty:
        lines.append("No design passed the complete custom intersection.")
    else:
        columns = [name for name in [_identifier_column(result), *[t["name"] for t in thresholds], "pass_filters"]
                   if name and name in result]
        lines.append("| " + " | ".join(columns) + " |")
        lines.append("|" + "|".join("---" for _ in columns) + "|")
        for row in selected_frame[columns].itertuples(index=False, name=None):
            lines.append("| " + " | ".join(str(value) for value in row) + " |")
    lines.extend(["", "## Interpretation note", "",
                  "Custom thresholds and BoltzGen's native `pass_filters` are separate decisions. "
                  "Budget-ranked folders can contain fallback designs that did not pass native filters; "
                  "experimental validation remains necessary.", "", "## Figures", ""])
    for figure in summary["figures"]:
        lines.append(f"- [{figure['caption']}]({figure['png']})")
    markdown = "\n".join(lines) + "\n"
    (output / "scientific_report.md").write_text(markdown, encoding="utf-8")
    body = "<br>\n".join(html.escape(line) for line in lines)
    (output / "scientific_report.html").write_text(
        "<!doctype html><meta charset='utf-8'><title>BOLTRA report</title>"
        "<style>body{font:16px system-ui;max-width:1000px;margin:3rem auto;line-height:1.5}</style>"
        f"<body>{body}</body>\n", encoding="utf-8")


def analyze(run_directory: Path, output_directory: Path, central: Metric,
            supporting: list[tuple[Metric, float]], central_threshold: float) -> dict:
    from . import __version__

    frame, source = load_metrics(run_directory)
    output_directory.mkdir(parents=True, exist_ok=False)
    configured = [(central, central_threshold), *supporting]
    passes = {metric.name: apply_threshold(frame, metric, threshold) for metric, threshold in configured}
    intersection = pd.Series(True, index=frame.index)
    for mask in passes.values():
        intersection &= mask
    native = _native_pass(frame)
    result = frame.copy()
    for name, mask in passes.items():
        result[f"pass_custom_{name}"] = mask
    result["passes_custom_intersection"] = intersection
    result["passes_native_boltzgen_filters"] = native
    result["selection_category"] = "neither"
    result.loc[intersection & ~native, "selection_category"] = "custom_only"
    result.loc[native & ~intersection, "selection_category"] = "native_only"
    result.loc[intersection & native, "selection_category"] = "both"
    result.to_csv(output_directory / "screened_designs.csv", index=False)
    result.loc[intersection].to_csv(output_directory / "selected_candidates.csv", index=False)
    result.loc[native].to_csv(output_directory / "native_candidates.csv", index=False)
    result.loc[intersection & native].to_csv(output_directory / "both_candidates.csv", index=False)

    identifier = _identifier_column(result)
    figures, captions = [], []
    palette = {"neither": "#999999", "custom_only": "#009E73", "native_only": "#E69F00", "both": "#0072B2"}
    markers = {"neither": "x", "custom_only": "o", "native_only": "s", "both": "D"}
    for metric, threshold in supporting:
        fig, ax = plt.subplots(figsize=(7.5, 5.5))
        for category in ("neither", "custom_only", "native_only", "both"):
            mask = result["selection_category"].eq(category)
            if mask.any():
                ax.scatter(result.loc[mask, metric.name], result.loc[mask, central.name],
                           color=palette[category], marker=markers[category], s=42,
                           label=category.replace("_", " "), zorder=3)
        ax.axhline(central_threshold, color="#374151", linestyle="--", linewidth=1.2,
                   label=f"{central.label} threshold = {central_threshold:g}")
        ax.axvline(threshold, color="#6B7280", linestyle=":", linewidth=1.4,
                   label=f"{metric.label} threshold = {threshold:g}")
        if identifier:
            for sequence, index in enumerate(_label_indices(result, central), 1):
                row = result.loc[index]
                ax.annotate(compact_design_label(row[identifier], sequence),
                            (row[metric.name], row[central.name]),
                            xytext=(4, 4), textcoords="offset points", fontsize=6, alpha=.8)
        ax.set(xlabel=metric.label, ylabel=central.label, title=f"{central.label} vs {metric.label}")
        ax.legend(fontsize=7, frameon=False)
        ax.grid(alpha=.15)
        fig.tight_layout()
        stem = f"{central.name}_vs_{metric.name}"
        for suffix in ("png", "svg", "pdf"):
            fig.savefig(output_directory / f"{stem}.{suffix}", dpi=300 if suffix == "png" else None)
        plt.close(fig)
        label_note = ("All designs are labelled numerically." if len(result) <= 50 else
                      "Up to 30 selected designs are labelled numerically, prioritizing joint then custom and native passes.")
        caption = (f"{central.label} versus {metric.label}. Lines mark prespecified custom thresholds; "
                   f"symbols distinguish custom and native BoltzGen filter outcomes. {label_note}")
        figures.append({"png": f"{stem}.png", "svg": f"{stem}.svg", "pdf": f"{stem}.pdf", "caption": caption})
        captions.append(f"**{stem}.** {caption}")
    (output_directory / "figure_captions.md").write_text("\n\n".join(captions) + "\n", encoding="utf-8")

    candidate_counts = _copy_candidates(run_directory, output_directory, result, intersection, native)
    thresholds = [{"name": metric.name, "label": metric.label, "direction": metric.direction,
                   "description": metric.description, "threshold": threshold,
                   "statistics": metric_statistics(frame, metric)} for metric, threshold in configured]
    summary = {"generated_utc": datetime.now(timezone.utc).isoformat(), "boltra_version": __version__,
               "boltzgen_version": _software_version(["boltzgen", "--version"]),
               "metrics_source": str(source), "project_settings": _project_settings(run_directory),
               "total_designs": int(len(frame)), "thresholds": thresholds,
               "individual_pass_counts": {KNOWN_METRICS[name].label: int(mask.sum()) for name, mask in passes.items()},
               "selected_by_custom_intersection": int(intersection.sum()),
               "passed_native_boltzgen_filters": int(native.sum()),
               "passed_both": int((intersection & native).sum()),
               "candidate_structures_copied": candidate_counts, "figures": figures,
               "warning": "BoltzGen budget folders may contain fallback ranked designs" if not native.any() else None}
    (output_directory / "analysis_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    _write_report(output_directory, summary, result, intersection)
    return summary
