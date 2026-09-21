from __future__ import annotations

import sys
import shutil
from pathlib import Path

from .analysis import analyze, available_metrics, load_metrics, metric_statistics
from .project import create_project, safe_project_name, save_metadata
from .runner import (
    resource_warnings,
    resume_pipeline,
    run_pipeline,
    safe_worker_count,
    system_resources,
    validate,
)
from .scaffolds import (
    copy_scaffolds,
    custom_scaffold_definition,
    discover_scaffolds,
    scaffold_setup_guidance,
)
from .spec import (
    cyclotide_spec,
    peptide_spec,
    protein_binder_spec,
    protein_redesign_spec,
    scaffold_binder_spec,
    small_molecule_binder_spec,
    validate_ligand_structure_file,
    write_spec,
)
from .status import inspect_project, status_text
from .structure import convert_author_ranges, read_structure, select_boltzgen_ranges, write_mapping_csv


def ask(prompt: str, default: str | None = None) -> str:
    suffix = f" [{default}]" if default is not None else ""
    value = input(f"{prompt}{suffix}: ").strip()
    return value or (default or "")


def choose(prompt: str, options: list[str]) -> int:
    print(f"\n{prompt}\n")
    for index, option in enumerate(options, 1):
        print(f"  ({index}) {option}")
    while True:
        raw = ask("\nEnter your choice")
        if raw.isdigit() and 1 <= int(raw) <= len(options):
            return int(raw) - 1
        print("Please enter one listed option number.")


def ask_int(prompt: str, default: int, minimum: int = 1, maximum: int | None = None) -> int:
    while True:
        raw = ask(prompt, str(default))
        try:
            value = int(raw)
        except ValueError:
            print("Please enter a whole number.")
            continue
        if value < minimum or (maximum is not None and value > maximum):
            upper = f" and {maximum}" if maximum is not None else ""
            print(f"Enter a value between {minimum}{upper}.")
            continue
        return value


def ask_float(prompt: str) -> float:
    while True:
        try:
            return float(ask(prompt))
        except ValueError:
            print("Please enter a numerical threshold.")


def ask_file(prompt: str) -> Path:
    while True:
        path = Path(ask(prompt)).expanduser().resolve()
        if path.is_file():
            return path
        print(f"File not found: {path}. Include the filename, for example /path/target.pdb.")


def ask_project_path(prompt: str) -> Path:
    while True:
        path = Path(ask(prompt)).expanduser().resolve()
        try:
            inspect_project(path)
            return path
        except FileNotFoundError as exc:
            print(f"Not recognized: {exc}")


def _show_run_preflight(project: Path, num_designs: int) -> None:
    workers = safe_worker_count()
    resources = system_resources(project)
    ram = "unknown" if resources.total_ram_gib is None else f"{resources.total_ram_gib:g} GiB"
    available = "unknown" if resources.available_ram_gib is None else f"{resources.available_ram_gib:g} GiB"
    swap = "unknown" if resources.total_swap_gib is None else f"{resources.total_swap_gib:g} GiB"
    print("\nRun resource preflight:")
    print(f"  RAM: {ram} total; {available} currently available")
    print(f"  Swap: {swap} total")
    print(f"  Disk: {resources.free_disk_gib:g} GiB free at the project location")
    print(f"  BoltzGen data-loader workers: {workers} (safe default; override with BOLTRA_NUM_WORKERS)")
    for warning in resource_warnings(resources, num_designs, workers):
        print(f"  WARNING: {warning}")


def _report_run_result(code: int, project: Path, label: str) -> None:
    print(f"{label} exited with status {code}")
    if code != 0:
        print("Existing project outputs were preserved and can be inspected with menu option 2.")
        return
    try:
        status = inspect_project(project)
    except (FileNotFoundError, ValueError):
        return
    if status.complete:
        print("Post-run audit: every configured stage has the expected output count.")
        return
    print("WARNING: The command exited successfully, but the post-run audit found incomplete output:")
    for stage in status.stages:
        if not stage.complete:
            print(f"  {stage.name}: {stage.evidence}")
    print("Use menu option 2 to resume the incomplete stage. Do not treat this run as complete yet.")


def _run_new_project(project: Path, protocol: str, num_designs: int, budget: int) -> None:
    _show_run_preflight(project, num_designs)
    code = run_pipeline(project / "design_specification.yaml", project / "boltzgen_output",
                        protocol, num_designs, budget, project)
    _report_run_result(code, project, "BoltzGen")


def _residue_selection(chain_residues, prompt: str, blank_allowed: bool):
    while True:
        entered = ask(prompt, "") or None
        if not entered and blank_allowed:
            return None, [], None, None
        if not entered:
            print("At least one residue or range is required.")
            continue
        numbering = choose("Which numbering system did you enter?",
                           ["PDB/mmCIF author residue numbers", "BoltzGen sequential indices"])
        try:
            if numbering == 0:
                converted, selected = convert_author_ranges(entered, chain_residues)
                print(f"\nConverted author residues {entered} to BoltzGen indices {converted}.")
                label = "author"
            else:
                converted, selected = select_boltzgen_ranges(entered, chain_residues)
                print(f"\nValidated BoltzGen binding-site indices: {converted}.")
                label = "boltzgen"
        except ValueError as exc:
            print(f"Residue selection error: {exc}\nPlease enter the residues again.")
            continue
        decision = choose("Use this converted residue selection?", ["Yes", "Enter it again", "Cancel project"])
        if decision == 0:
            return converted, selected, entered, label
        if decision == 2:
            return "CANCEL", [], entered, label


def new_design() -> None:
    kind = choose("What would you like to design?",
                  ["Peptide", "Cyclotide", "De novo protein binder", "Redesign an existing protein",
                   "Protein binder for a small molecule", "Nanobody", "Antibody"])
    if kind == 4:
        new_small_molecule_design()
        return
    if kind in {5, 6}:
        new_scaffold_design("nanobody" if kind == 5 else "antibody")
        return
    target = ask_file("Path to the target PDB/mmCIF file")
    chains = read_structure(target)
    print("\nChains detected in target structure:\n")
    chain_ids = list(chains)
    for index, chain_id in enumerate(chain_ids, 1):
        residues = chains[chain_id]
        print(f"  ({index}) Chain {chain_id}: {len(residues)} resolved residues; "
              f"author range {residues[0].author_id}-{residues[-1].author_id}; BoltzGen range 1-{len(residues)}")
    chain = chain_ids[choose("Select the target chain", [f"Chain {value}" for value in chain_ids])]
    chain_residues = chains[chain]
    if kind == 3:
        scope = choose("Which part of the selected chain should be redesigned?",
                       ["Entire resolved chain", "Specific residues or ranges"])
        if scope == 0:
            site = f"1..{len(chain_residues)}"
            selected = list(chain_residues)
            original_site = site
            numbering = "boltzgen"
            print(f"\nThe complete resolved chain will be redesigned using BoltzGen indices {site}.")
        else:
            site, selected, original_site, numbering = _residue_selection(
                chain_residues, "Residues to redesign", False
            )
    else:
        site, selected, original_site, numbering = _residue_selection(
            chain_residues, "Binding-site residues (blank permits binding anywhere)", True
        )
    if site == "CANCEL":
        print("Project cancelled before any files were created.")
        return

    while True:
        name = safe_project_name(ask("Project name"))
        root = Path(ask("Directory in which to create the project", "~/boltzgen-projects")).expanduser().resolve()
        if root.exists() and not root.is_dir():
            print(f"That is a file, not a project directory: {root}")
        elif (root / name).exists():
            print(f"Project already exists: {root / name}. Choose another name or location.")
        else:
            break
    num_designs = ask_int("Total number of designs", 10, 1, 100000)
    budget = ask_int("Final quality/diversity budget", min(2, num_designs), 1, num_designs)

    project, copied_target = create_project(root, name, target)
    assert copied_target is not None
    write_mapping_csv(project / "residue_mapping.csv", chain_residues,
                      {residue.boltzgen_index for residue in selected})
    if kind == 0:
        minimum = ask_int("Minimum peptide length", 12)
        maximum = ask_int("Maximum peptide length", minimum, minimum)
        cyclic = choose("Select peptide architecture", ["Linear", "Head-to-tail cyclic"]) == 1
        spec = peptide_spec(copied_target.name, chain, site, minimum, maximum, cyclic)
        details = {"kind": "peptide", "length": [minimum, maximum], "cyclic": cyclic}
        protocol = "peptide-anything"
    elif kind == 1:
        print("\nEnter loop lengths I-II, II-III, III-IV, IV-V, V-VI, VI-I.")
        while True:
            try:
                loops = [int(value) for value in ask("Six loop lengths separated by spaces", "3 4 4 1 4 7").split()]
                spec, pattern = cyclotide_spec(copied_target.name, chain, site, loops)
                break
            except ValueError as exc:
                print(f"Invalid loop lengths: {exc}")
        print(f"Generated sequence: {pattern.sequence}")
        print(f"Cysteines: {', '.join(map(str, pattern.cysteine_positions))}")
        print("Disulfides: " + ", ".join(f"{a}-{b}" for a, b in pattern.bonds))
        details = {"kind": "cyclotide", "loops": loops, "sequence": pattern.sequence,
                   "disulfides": pattern.bonds}
        protocol = "peptide-anything"
    elif kind == 2:
        minimum = ask_int("Minimum protein-binder length", 100)
        maximum = ask_int("Maximum protein-binder length", max(150, minimum), minimum)
        spec = protein_binder_spec(copied_target.name, chain, site, minimum, maximum)
        details = {"kind": "protein_binder", "length": [minimum, maximum]}
        protocol = "protein-anything"
    else:
        spec = protein_redesign_spec(copied_target.name, chain, site)
        details = {
            "kind": "protein_redesign",
            "redesign_scope": "entire_resolved_chain" if scope == 0 else "specific_residues",
            "residues_to_redesign_input": original_site,
            "residues_to_redesign_input_numbering": numbering,
            "residues_to_redesign_boltzgen": site,
        }
        protocol = "protein-redesign"

    spec_path = project / "design_specification.yaml"
    write_spec(spec, spec_path)
    binding_input = original_site if kind != 3 else None
    binding_numbering = numbering if kind != 3 else None
    binding_boltzgen = site if kind != 3 else None
    save_metadata(project, {**details, "target": copied_target.name, "target_source": str(target),
                            "target_chain": chain, "binding_site_input": binding_input,
                            "binding_site_input_numbering": binding_numbering,
                            "binding_site_boltzgen": binding_boltzgen, "num_designs": num_designs,
                            "budget": budget, "protocol": protocol})
    print(f"\nSpecification written to {spec_path}")
    if validate(spec_path, project) != 0:
        print("BoltzGen validation failed. The pipeline was not started.")
        return
    if choose("Validation passed. What next?", ["Save without running", "Begin BoltzGen run"]) == 1:
        _run_new_project(project, protocol, num_designs, budget)


def _new_project_location() -> tuple[str, Path]:
    while True:
        name = safe_project_name(ask("Project name"))
        root = Path(ask("Directory in which to create the project", "~/boltzgen-projects")).expanduser().resolve()
        if root.exists() and not root.is_dir():
            print(f"That is a file, not a project directory: {root}")
        elif (root / name).exists():
            print(f"Project already exists: {root / name}. Choose another name or location.")
        else:
            return name, root


def new_small_molecule_design() -> None:
    print("\nDesign a de novo protein binder against one small-molecule target.")
    ligand_choice = choose("How will you specify the small molecule?", [
        "Chemical Component Dictionary (CCD) code",
        "SMILES string",
        "Ligand-only PDB/mmCIF structure file",
    ])
    ligand_file = None
    if ligand_choice == 0:
        ligand_format = "ccd"
        while True:
            ligand_value = ask("Target CCD code").upper()
            try:
                small_molecule_binder_spec(1, 1, ligand_format, ligand_value)
                break
            except ValueError as exc:
                print(f"Invalid CCD code: {exc}")
        ligand_display = ligand_value
    elif ligand_choice == 1:
        ligand_format = "smiles"
        while True:
            ligand_value = ask("Target SMILES string")
            try:
                small_molecule_binder_spec(1, 1, ligand_format, ligand_value)
                break
            except ValueError as exc:
                print(f"Invalid SMILES input: {exc}")
        ligand_display = ligand_value
    else:
        ligand_format = "file"
        while True:
            ligand_file = ask_file("Path to a ligand-only PDB/mmCIF file")
            try:
                validate_ligand_structure_file(ligand_file)
                break
            except ValueError as exc:
                print(f"Ligand structure error: {exc}")
                if choose("How would you like to proceed?", [
                    "Enter another complete PDB/mmCIF file",
                    "Cancel and restart using CCD or SMILES",
                ]) == 1:
                    print("Project cancelled before any files were created.")
                    return
        ligand_value = ligand_file.name
        ligand_display = ligand_file.name
    minimum = ask_int("Minimum protein-binder length", 150)
    maximum = ask_int("Maximum protein-binder length", max(200, minimum), minimum)
    name, root = _new_project_location()
    num_designs = ask_int("Total number of designs", 10, 1, 100000)
    budget = ask_int("Final quality/diversity budget", min(2, num_designs), 1, num_designs)
    project, _ = create_project(root, name)
    if ligand_file is not None:
        shutil.copy2(ligand_file, project / "inputs" / ligand_file.name)
    spec = small_molecule_binder_spec(minimum, maximum, ligand_format, ligand_value)
    spec_path = project / "design_specification.yaml"
    write_spec(spec, spec_path)
    save_metadata(project, {
        "kind": "small_molecule_binder", "ligand_format": ligand_format,
        "ligand_display": ligand_display, "ligand_source": str(ligand_file) if ligand_file else None,
        "length": [minimum, maximum], "num_designs": num_designs, "budget": budget,
        "protocol": "protein-small_molecule",
    })
    print(f"\nSpecification written to {spec_path}")
    if validate(spec_path, project) != 0:
        print("BoltzGen validation failed. The pipeline was not started.")
        return
    if choose("Validation passed. What next?", ["Save without running", "Begin BoltzGen run"]) == 1:
        _run_new_project(project, "protein-small_molecule", num_designs, budget)


def _show_chains(chains: dict) -> list[str]:
    print("\nChains detected in structure:\n")
    chain_ids = list(chains)
    for index, chain_id in enumerate(chain_ids, 1):
        residues = chains[chain_id]
        print(f"  ({index}) Chain {chain_id}: {len(residues)} resolved residues; "
              f"author range {residues[0].author_id}-{residues[-1].author_id}; BoltzGen range 1-{len(residues)}")
    return chain_ids


def _select_scaffold_subset(scaffolds):
    print(f"\nDiscovered {len(scaffolds)} compatible scaffold(s):\n")
    for index, scaffold in enumerate(scaffolds, 1):
        print(f"  ({index}) {scaffold.name}")
    if choose("Which discovered scaffolds should BoltzGen sample?",
              ["All discovered scaffolds", "Select one or more scaffolds"]) == 0:
        return scaffolds
    while True:
        raw = ask("Scaffold numbers separated by spaces")
        try:
            indexes = [int(value) - 1 for value in raw.split()]
            if not indexes or len(indexes) != len(set(indexes)) or any(i < 0 or i >= len(scaffolds) for i in indexes):
                raise ValueError
            return [scaffolds[index] for index in indexes]
        except ValueError:
            print("Enter one or more distinct scaffold numbers separated by spaces.")


def new_scaffold_design(kind: str) -> None:
    title = "nanobody" if kind == "nanobody" else "antibody"
    target = ask_file("Path to the target PDB/mmCIF file")
    target_chains = read_structure(target)
    target_ids = _show_chains(target_chains)
    target_chain = target_ids[choose("Select the target chain", [f"Chain {value}" for value in target_ids])]
    target_residues = target_chains[target_chain]
    site, site_selected, original_site, numbering = _residue_selection(
        target_residues, "Binding-site residues (blank permits binding anywhere)", True
    )
    if site == "CANCEL":
        print("Project cancelled before any files were created.")
        return

    source = choose(f"Select the {title} scaffold source", [
        "BoltzGen template scaffolds discovered from the installation",
        "Custom framework structure with fixed-length CDR redesign",
    ])
    selected_scaffolds = []
    framework = None
    framework_chains: dict[str, str] = {}
    framework_mappings: dict[str, tuple[list, list]] = {}
    if source == 0:
        selected_scaffolds = discover_scaffolds(kind)
        if not selected_scaffolds:
            print("\nNo scaffold directory was found automatically.")
            print(scaffold_setup_guidance())
            print("You may enter the BoltzGen checkout now, or press Ctrl+C and follow TUTORIAL.md.")
            root = Path(ask("Path to the BoltzGen checkout or scaffold directory")).expanduser().resolve()
            selected_scaffolds = discover_scaffolds(kind, root)
        if not selected_scaffolds:
            raise FileNotFoundError(
                f"No complete {title} scaffold YAML/structure pairs were found. "
                "Set BOLTZGEN_ROOT to the BoltzGen checkout and try again."
            )
        selected_scaffolds = _select_scaffold_subset(selected_scaffolds)
    else:
        framework = ask_file(f"Path to the custom {title} framework PDB/mmCIF file")
        framework_structure = read_structure(framework)
        framework_ids = _show_chains(framework_structure)
        if kind == "nanobody":
            chain = framework_ids[choose("Select the nanobody framework chain",
                                         [f"Chain {value}" for value in framework_ids])]
            cdrs, selected, entered, input_numbering = _residue_selection(
                framework_structure[chain], "Nanobody CDR residues to design", False
            )
            if cdrs == "CANCEL":
                print("Project cancelled before any files were created.")
                return
            framework_chains[chain] = cdrs
            framework_mappings[chain] = (framework_structure[chain], selected)
        else:
            heavy = framework_ids[choose("Select the antibody heavy chain",
                                         [f"Chain {value}" for value in framework_ids])]
            remaining = [value for value in framework_ids if value != heavy]
            if not remaining:
                raise ValueError("A custom antibody framework requires distinct heavy and light chains")
            light = remaining[choose("Select the antibody light chain",
                                     [f"Chain {value}" for value in remaining])]
            for label, chain in (("Heavy-chain", heavy), ("Light-chain", light)):
                cdrs, selected, entered, input_numbering = _residue_selection(
                    framework_structure[chain], f"{label} CDR residues to design", False
                )
                if cdrs == "CANCEL":
                    print("Project cancelled before any files were created.")
                    return
                framework_chains[chain] = cdrs
                framework_mappings[chain] = (framework_structure[chain], selected)

    name, root = _new_project_location()
    num_designs = ask_int("Total number of designs", 10, 1, 100000)
    budget = ask_int("Final quality/diversity budget", min(2, num_designs), 1, num_designs)
    project, copied_target = create_project(root, name, target)
    assert copied_target is not None
    write_mapping_csv(project / "target_residue_mapping.csv", target_residues,
                      {residue.boltzgen_index for residue in site_selected})
    if source == 0:
        references = copy_scaffolds(selected_scaffolds, project)
        scaffold_names = [scaffold.name for scaffold in selected_scaffolds]
    else:
        scaffold_dir = project / "inputs" / "scaffolds"
        scaffold_dir.mkdir()
        framework_name = f"custom_{title}_framework{framework.suffix.lower()}"
        shutil.copy2(framework, scaffold_dir / framework_name)
        scaffold_yaml = scaffold_dir / f"custom_{title}.yaml"
        write_spec(custom_scaffold_definition(framework_name, framework_chains), scaffold_yaml)
        references = [f"inputs/scaffolds/{scaffold_yaml.name}"]
        scaffold_names = [scaffold_yaml.stem]
        for chain, (residues, selected) in framework_mappings.items():
            write_mapping_csv(project / f"framework_chain_{chain}_mapping.csv", residues,
                              {residue.boltzgen_index for residue in selected})
    specification = scaffold_binder_spec(copied_target.name, target_chain, site, references)
    spec_path = project / "design_specification.yaml"
    write_spec(specification, spec_path)
    protocol = "nanobody-anything" if kind == "nanobody" else "antibody-anything"
    save_metadata(project, {
        "kind": kind, "target": copied_target.name, "target_source": str(target),
        "target_chain": target_chain, "binding_site_input": original_site,
        "binding_site_input_numbering": numbering, "binding_site_boltzgen": site,
        "scaffold_source": "discovered" if source == 0 else "custom",
        "scaffolds": scaffold_names, "custom_framework_source": str(framework) if framework else None,
        "custom_framework_cdrs_boltzgen": framework_chains or None,
        "num_designs": num_designs, "budget": budget, "protocol": protocol,
    })
    print(f"\nSpecification written to {spec_path}")
    print("Scaffolds snapshotted in: " + str(project / "inputs" / "scaffolds"))
    if validate(spec_path, project) != 0:
        print("BoltzGen validation failed. The pipeline was not started.")
        return
    if choose("Validation passed. What next?", ["Save without running", "Begin BoltzGen run"]) == 1:
        _run_new_project(project, protocol, num_designs, budget)


def project_manager() -> None:
    path = ask_project_path("Path to the BOLTRA project or BoltzGen output directory")
    status = inspect_project(path)
    print("\n" + status_text(status))
    if status.complete:
        return
    if status.problems:
        print("\nRecovery is paused because empty/unreadable artifacts were detected. Inspect them before resuming.")
        return
    if not status.started:
        if choose("Start this saved, validated project?", ["No", "Yes — begin BoltzGen run"]) == 1:
            import json
            settings = json.loads((status.project / "project_settings.json").read_text(encoding="utf-8"))
            _run_new_project(status.project, settings["protocol"], int(settings["num_designs"]),
                             int(settings["budget"]))
        return
    if choose("Resume from the first incomplete stage?",
              ["No", "Yes — use existing compatible artifacts"]) == 1:
        print("\nBOLTRA will preserve existing files. BoltzGen may repeat work from a partially completed stage.")
        if choose("Start recovery now?", ["No", "Yes"]) == 1:
            _show_run_preflight(status.project, status.expected_designs or 0)
            code = resume_pipeline(status.output, status.remaining_steps, status.project)
            _report_run_result(code, status.project, "BoltzGen recovery")


def _analysis_directory(run_dir: Path) -> Path:
    root = run_dir / "workbench_analyses"
    while True:
        name = safe_project_name(ask("Analysis name", "analysis_1"))
        output = root / name
        if not output.exists():
            return output
        print(f"Analysis already exists: {output}. Choose a new name; BOLTRA will not overwrite it.")


def post_design() -> None:
    while True:
        run_dir = Path(ask("Path to a completed BoltzGen output directory")).expanduser().resolve()
        try:
            frame, source = load_metrics(run_dir)
            break
        except FileNotFoundError as exc:
            print(f"Output not ready: {exc}")
    metrics = available_metrics(frame)
    if len(metrics) < 2:
        raise RuntimeError("Fewer than two recognized scientific metrics were found")
    print(f"Loaded {len(frame)} designs from {source}")
    print("\nObserved metric ranges (use these to choose defensible thresholds):\n")
    for metric in metrics:
        stats = metric_statistics(frame, metric)
        print(f"  {metric.label}: min {stats['minimum']:.4g} | median {stats['median']:.4g} | max {stats['maximum']:.4g}")
    central = metrics[choose("Select one central filtering parameter",
                             [f"{m.label} ({m.direction} is better)" for m in metrics])]
    central_threshold = ask_float(f"Threshold for {central.label}")
    choices = [metric for metric in metrics if metric.name != central.name]
    print("\nSupporting parameters (enter one or more numbers separated by spaces):\n")
    for index, metric in enumerate(choices, 1):
        print(f"  ({index}) {metric.label} ({metric.direction} is better)")
    while True:
        raw = ask("Supporting choices")
        try:
            indexes = [int(value) - 1 for value in raw.split()]
            if not indexes or len(indexes) != len(set(indexes)) or any(i < 0 or i >= len(choices) for i in indexes):
                raise ValueError
            break
        except ValueError:
            print("Enter one or more distinct listed numbers separated by spaces.")
    supporting = [(choices[index], ask_float(f"Threshold for {choices[index].label}")) for index in indexes]
    output = _analysis_directory(run_dir)
    summary = analyze(run_dir, output, central, supporting, central_threshold)
    print(f"\nAnalysis written to {output}")
    print(f"Custom intersection: {summary['selected_by_custom_intersection']}/{summary['total_designs']}")
    print(f"Native BoltzGen filters: {summary['passed_native_boltzgen_filters']}/{summary['total_designs']}")
    print(f"Passed both: {summary['passed_both']}/{summary['total_designs']}")
    print("Scientific report: " + str(output / "scientific_report.html"))
    if summary["warning"]:
        print(f"WARNING: {summary['warning']}")


def main() -> None:
    print("\nBOLTRA v1.0.0 — guided BoltzGen design, recovery, and analysis\n")
    option = choose("Select a function", ["Create a guided design or redesign project",
                                           "Open a project: inspect status or resume",
                                           "Perform post-design analysis", "Exit"])
    try:
        if option == 0:
            new_design()
        elif option == 1:
            project_manager()
        elif option == 2:
            post_design()
    except (ValueError, FileNotFoundError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2) from exc
    except KeyboardInterrupt:
        print("\nInterrupted by user. Existing project outputs were preserved.", file=sys.stderr)
        raise SystemExit(130)


if __name__ == "__main__":
    main()
