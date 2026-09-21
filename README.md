# BOLTRA v1.0.0

<img width="6000" height="3375" alt="BOLTRA_Main" src="https://github.com/user-attachments/assets/f4045ffe-77ca-4eb4-a413-f47a464b3ac4" />

**BOLTRA - BoltzGen Orchestration Layer for Targeted Design and Results Analysis**
is an independent, menu-driven companion to BoltzGen. It converts common design
questions into validated BoltzGen specifications, manages interrupted projects,
and turns output metrics into traceable candidate-selection reports.

BOLTRA does not replace BoltzGen and does not contain its models. It calls the
locally installed BoltzGen command-line interface and preserves the native
specification and result files.

## Supported workflows

| Mode | BOLTRA guidance |
|---|---|
| Peptide | Linear or head-to-tail cyclic binders with selectable length |
| Cyclotide | Six-loop cystine-knot notation and explicit disulfide topology |
| De novo protein binder | Target chain, binding site, and binder-length range |
| Protein redesign | Whole-chain or selected-residue redesign |
| Small-molecule binder | CCD code, SMILES, or a complete ligand structure file |
| Nanobody | Dynamically discovered BoltzGen VHH templates or a custom framework |
| Antibody | Dynamically discovered BoltzGen Fab templates or a custom framework |

All structure-guided modes support explicit conversion between PDB/mmCIF author
residue numbers and the sequential indices used by BoltzGen.

## Requirements

- Linux workstation
- Python 3.11 or newer; Python 3.12 is recommended
- A working BoltzGen installation and its downloaded model/data artifacts
- NVIDIA GPU and CUDA-capable BoltzGen environment for practical design runs
- Adequate disk space, system RAM, and swap for the requested scale

See [TUTORIAL.md](TUTORIAL.md) for a beginner-oriented installation guide,
hardware checks, scaffold setup, and seven complete worked examples.

## Quick installation

After installing BoltzGen in a dedicated Conda environment, extract this release
and install BOLTRA into the same environment:

```bash
conda activate boltzgen
cd /path/to/boltzgen-workbench
python -m pip install -e .
hash -r
python -c "import boltzgen_workbench; print(boltzgen_workbench.__version__)"
boltra
```

The version check should print `1.0.0`. The historical `boltz-workbench` command
is retained as a compatibility alias.

## Main menu

```text
(1) Create a guided design or redesign project
(2) Open a project: inspect status or resume
(3) Perform post-design analysis
(4) Exit
```

Every new project contains a copied input structure when applicable, a readable
`design_specification.yaml`, metadata, residue-mapping records, and a dedicated
`boltzgen_output` directory after execution begins. BOLTRA validates generated
specifications with `boltzgen check` before offering to run them.

## Recovery and resource safety

BOLTRA inspects expected stage artifacts and resumes from the first incomplete
stage using BoltzGen's reuse mechanism. Existing files are preserved, although a
partially completed stage may be repeated because recovery is stage-aware rather
than a checkpoint inside an individual GPU operation.

Before new and resumed runs, BOLTRA reports available RAM, swap, and project-disk
space. It uses zero separate data-loader workers by default to limit peak host-RAM
use. Advanced users may override this deliberately:

```bash
BOLTRA_NUM_WORKERS=1 boltra
```

After execution, BOLTRA distinguishes generated structures from metric-table
rows. Identifiers such as `00` through `19` represent 20 designs; independently,
an analysis table with 19 rows means that only 19 designs received metric records.

## Post-design analysis

BOLTRA discovers numerical columns from the actual BoltzGen metrics table and
shows each metric's minimum, median, maximum, and optimization direction. A user
chooses one central threshold plus optional supporting thresholds. BOLTRA keeps
custom selection separate from BoltzGen's native `pass_filters` decision.

Each named analysis produces:

- complete screened and selected-candidate CSV files;
- custom-only, native-only, both-pass, and neither classifications;
- PNG, SVG, and PDF scatter plots with scalable numeric labels;
- figure captions and a machine-readable JSON summary;
- Markdown and HTML scientific reports;
- grouped candidate structures when source structures can be located.

Thresholds are project-specific screening choices, not universal biological
cutoffs. Computational prioritization does not replace experimental validation.

## Reproducibility

BOLTRA records commands, timestamps, software versions, platform/GPU information,
effective worker settings, and execution outcomes in
`boltra_execution_log.jsonl`. Scaffold files selected from BoltzGen are copied into
the project so later changes to the installed template collection do not alter an
existing specification.

The accompanying v1.0.0 pilot dataset exercises all seven design modes with 20
generated structures per mode. Its own README documents inputs, counts, analysis
completeness, thresholds, and directory organization.

## Development and tests

```bash
python -m pip install -e '.[test]'
python -m pytest
```

Tests cover specification generation, residue mapping, status/recovery logic,
metric analysis, dynamic scaffold discovery, and ligand-file validation. Model
checkpoints and large structures are not bundled with the software release.

## Scope and independence

BOLTRA was independently implemented against the public BoltzGen CLI, schema,
examples, and output formats. It contains no Tamarind Bio code, assets, or user
interface. BoltzGen is a separate dependency and should be cited alongside BOLTRA
in scientific work.

## Citation and license

Citation metadata are provided in [CITATION.cff](CITATION.cff). Until the BOLTRA
article and archived software DOI are available, cite the GitHub release version
and the corresponding Zenodo record. BOLTRA is distributed under the [MIT
License](LICENSE).
