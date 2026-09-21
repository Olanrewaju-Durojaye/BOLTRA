# Changelog

## 1.0.0

- Declare the first documented, seven-mode public release.
- Add a complete beginner tutorial covering environment setup, model downloads,
  dynamic scaffold resources, all guided workflows, recovery, analysis, and troubleshooting.
- Add citation metadata and release-ready project documentation.
- Reject the common bare HETATM-only ligand PDB form before BoltzGen validation and
  explain when to use a CCD code, SMILES, or a complete structure file.
- Improve nanobody/antibody scaffold setup guidance when example resources are absent.
- Preserve v0.6.1 low-memory execution, artifact auditing, recovery, residue mapping,
  dynamic scaffold discovery, and post-design analysis behavior.

## 0.6.1

- Use zero separate BoltzGen data-loader workers by default to reduce peak host-memory pressure during
  analysis; permit an explicit `BOLTRA_NUM_WORKERS` override for larger systems.
- Report RAM, swap, and free project-disk space before new and resumed runs, with
  actionable warnings for constrained or unusually large runs.
- Audit every configured pipeline stage after a command exits and clearly flag
  partial results such as 19 metric rows for 20 requested designs.
- Fix CSV row counting in project status, which could incorrectly report zero
  analysis and filtering rows in v0.6.0.
- Record the resource snapshot and effective worker setting in the execution log.
- Preserve compatibility with existing v0.1-v0.6 project directories and outputs.

## 0.6.0

- Added nanobody and antibody design through BoltzGen's dedicated protocols.
- Discover complete scaffold YAML/structure pairs dynamically from an installed
  distribution, a BoltzGen checkout, `BOLTZGEN_ROOT`, or a user-supplied location.
- Let users sample every discovered scaffold or a selected subset, and snapshot
  all chosen scaffold files inside the project for reproducibility.
- Added custom fixed-length CDR redesign for one-chain nanobody frameworks and
  two-chain antibody frameworks, including residue-number conversion and mapping CSVs.

## 0.5.0

- Added de novo protein binders for small molecules specified by CCD code, SMILES,
  or a ligand-only PDB/mmCIF structure.
- Added `protein-small_molecule` execution with affinity-stage recovery and recognized
  Boltz-2 binder-probability and quantitative affinity metrics.
- Recognize validated projects saved before their first run and offer to start them.

## 0.4.0

- Add guided de novo protein-binder design using the official `protein-anything` protocol.
- Add structure-guided protein redesign using the official `protein-redesign` design mask.
- Support complete-chain redesign or user-selected residue ranges.
- Apply author-to-BoltzGen residue conversion and auditable mapping to redesign selections.
- Add protein-binder length-range validation and protein-specific project metadata.
- Extend interruption status detection to the protein-binder `design_folding` stage.
- Include redesign selections in generated scientific reports.
- Recognize BoltzGen's redesign-specific `design_residue_iptm` metric.

## 0.3.1

- Replace long design identifiers on figures with their trailing numeric ID.
- For more than 50 designs, label no more than 30 selected candidates.
- Prioritize joint custom/native passes, followed by custom-only and native-only designs.
- Preserve complete design identifiers in every CSV and report for traceability.

## 0.3.0

- Add project inspection and interruption recovery through BoltzGen `execute --reuse`.
- Detect the first incomplete pipeline stage and warn about empty artifacts.
- Record commands, timestamps, exit status, versions, platform, and GPU in JSON Lines.
- Show metric minimum, median, and maximum before threshold entry.
- Report per-filter, custom-intersection, native-filter, and joint pass counts.
- Export labeled, colorblind-friendly plots as PNG, SVG, and PDF with captions.
- Create uniquely named analysis directories without overwriting prior work.
- Collect custom-only, native-only, and jointly selected structures.
- Generate Markdown and HTML scientific reports with settings and provenance.

## 0.2.0

- Adopted the BOLTRA product name while retaining `boltz-workbench` as a command alias.
- Added automatic PDB and mmCIF chain discovery.
- Added explicit author-number versus BoltzGen-index selection.
- Added conversion of author residue ranges, including numbering gaps and insertion codes.
- Added an auditable per-project `residue_mapping.csv`.
- Added friendly errors for file paths entered as project directories and existing project names.
- Added early validation of design count and budget.

## 0.1.0

- Added guided peptide and cyclotide design.
- Added BoltzGen specification validation and execution.
- Added configurable post-design metric screening and plots.
- Re-prompt after common path, number, budget, residue, and menu-entry errors.
