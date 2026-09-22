# BOLTRA v1.0.0 Tutorial

This tutorial starts with a new Linux workstation and ends with design creation,
safe recovery, and candidate analysis. Commands beginning with `$` are entered in
a terminal; do not type the `$` itself.

## 1. What BOLTRA does

BOLTRA is a guided front end for BoltzGen. It asks plain-language questions,
writes a native BoltzGen YAML specification, validates it, optionally runs the
pipeline, and analyzes the resulting metrics. BoltzGen performs the generative
model calculations; BOLTRA provides orchestration, safeguards, provenance, and
reporting.

The seven modes are peptide, cyclotide, de novo protein binder, protein redesign,
protein binder for a small molecule, nanobody, and antibody.

## 2. Check the workstation

Confirm disk space and an NVIDIA GPU before downloading models:

```bash
df -h ~
nvidia-smi
free -h
```

BoltzGen's model/data download is several gigabytes, while design outputs can grow
substantially with design count and target size. A GPU is strongly recommended.
For a workstation with about 32 GiB RAM, 16–32 GiB swap is a useful safety margin;
swap prevents some abrupt failures but is not a substitute for RAM and is slower.

To inspect existing swap:

```bash
swapon --show
```

If an administrator chooses to add swap, use the operating system's documented
procedure and verify free disk space first. Do not append a second `/etc/fstab`
entry when one already exists.

## 3. Install Conda and create the BoltzGen environment

Install Miniconda or another compatible Conda distribution, then create a clean
Python 3.12 environment:

```bash
conda create -n boltzgen python=3.12 -y
conda activate boltzgen
python --version
```

The shell prompt should begin with `(boltzgen)`. If Conda reports that the
environment does not exist on another computer, create it there; copying BOLTRA's
ZIP file does not copy the Conda environment or BoltzGen models.

## 4. Install and verify BoltzGen

Install BoltzGen in the activated environment:

```bash
python -m pip install boltzgen
boltzgen --version
```

Download all required artifacts:

```bash
boltzgen download all
```

An unauthenticated Hugging Face warning is informational. A Hugging Face token may
improve rate limits but is not required for public artifacts. Verify the CLI:

```bash
boltzgen --help
boltzgen download --help
```

If the environment already contains BoltzGen, do not reinstall it merely to install
BOLTRA. BOLTRA v1.0.0 was developed against BoltzGen 0.3.2; record the exact version
used for every study because later BoltzGen schemas or metrics may change.

## 5. Install BOLTRA

Extract the release into its own directory. Use the exact downloaded filename:

```bash
conda activate boltzgen
cd ~/Downloads
mkdir -p BOLTRA-v1.0.0
unzip BOLTRA-v1.0.0.zip -d BOLTRA-v1.0.0
cd ~/Downloads/BOLTRA-v1.0.0/boltzgen-workbench
python -m pip install -e .
hash -r
python -c "import boltzgen_workbench; print(boltzgen_workbench.__version__)"
boltra
```

If the archive extracts directly without the outer `boltzgen-workbench` directory,
enter the directory that contains `pyproject.toml`. If `boltra` is not found, confirm
that installation and execution occurred in the same active Conda environment:

```bash
which python
which boltzgen
which boltra
```

## 6. Make BoltzGen template scaffolds discoverable

Peptide, cyclotide, protein, redesign, and small-molecule modes need only the
installed package. Nanobody and antibody template modes also need BoltzGen's example
scaffold YAML and structure files. Some package installations omit those examples.

Download an official BoltzGen source snapshot without requiring `git` or `curl`:

```bash
cd ~
python - <<'PY'
from urllib.request import urlretrieve

url = "https://github.com/HannesStark/boltzgen/archive/refs/heads/main.zip"
output = "boltzgen-source.zip"
print("Downloading BoltzGen source...")
urlretrieve(url, output)
print("Saved to:", output)
PY
unzip -q ~/boltzgen-source.zip -d ~
```

Set the checkout location permanently, then load it in the current terminal:

```bash
echo 'export BOLTZGEN_ROOT="$HOME/boltzgen-main"' >> ~/.bashrc
export BOLTZGEN_ROOT="$HOME/boltzgen-main"
```

Confirm that templates exist:

```bash
find "$BOLTZGEN_ROOT/example" -type f \( -name "*.yaml" -o -name "*.yml" \) \
  | grep -Ei 'nanobody_scaffolds|fab_scaffolds' | head
```

For strict reproducibility, archive or identify the exact BoltzGen tag/commit rather
than relying indefinitely on a changing `main` branch. BOLTRA copies every selected
template and its referenced structure into the new project.

## 7. Prepare input structures and residue selections

Use a valid `.pdb`, `.cif`, or `.mmcif` file and enter the complete filename—not a
directory. BOLTRA lists detected chains and their two coordinate systems:

- **Author numbers** are residue identifiers written in the PDB/mmCIF file.
- **BoltzGen indices** are one-based sequential positions among resolved residues in
  the selected chain.

When author numbers are selected, BOLTRA converts them, displays the result for
confirmation, and writes `residue_mapping.csv`. Unresolved residues cannot be
invented or selected.

Accepted selection syntax includes `90-94,97,102-108`. Spaces are optional. A blank
binding site permits binding anywhere when the protocol allows it.

## 8. Start BOLTRA

```bash
conda activate boltzgen
boltra
```

Choose **Create a guided design or redesign project**. A design count controls how
many structures BoltzGen attempts to generate. The final budget controls how many
ranked/diverse structures the native filtering step retains; it cannot exceed the
design count.

The following examples reproduce the input conditions of the v1.0.0 functional
pilot dataset. They are demonstrations, not universal biological defaults.

## 9. Peptide example

- Mode: `Peptide`
- Target: `1T2P.pdb`, chain A
- Binding site, author numbering:
  `90-94,97,102-108,114-120,162-164,172,180-184,187-189,191-195,197,199`
- Designs/budget: `20/20`
- Length: `12–15`
- Architecture: `Head-to-tail cyclic`

BOLTRA maps the binding site to
`29..33,36,41..47,53..59,101..103,111,119..123,126..128,130..134,136,138`,
writes the specification, and runs `boltzgen check` before execution.

## 10. Cyclotide example

Use the same target, chain, binding site, and `20/20` design settings. Enter the six
loop lengths:

```text
3 4 4 1 4 7
```

BOLTRA generates the compact specification `C3C4C4C1C4C7`. Here `C` is a fixed
cysteine and each number is the count of designed residues in the following loop.
The six cysteines occur at positions `1, 5, 10, 15, 17, 22`; the kalata B1-style
cystine-knot topology is `1–15`, `5–17`, and `10–22` (I–IV, II–V, III–VI).
The final loop continues from cysteine VI back to cysteine I because the peptide is
cyclic.

## 11. De novo protein-binder example

Use `1T2P.pdb`, chain A, and the same author-numbered binding site. Set:

- Designs/budget: `20/20`
- Minimum binder length: `100`
- Maximum binder length: `120`

The generated specification uses BoltzGen's protein-binder protocol while retaining
the target chain as context.

## 12. Protein-redesign example

- Mode: `Redesign an existing protein`
- Target: `4J9A.pdb`, chain A
- Choose `Specific residues or ranges`
- Author-numbered redesign selection:
  `11,14,22,32-34,36-42,45,53-58,61-66,84,86-91,95-99,101,105,115,117,119-123`
- Designs/budget: `20/20`

BOLTRA maps that selection to
`10,13,21,31..33,35..41,44,52..57,60..65,83,85..90,94..98,100,104,114,116,118..122`.
To redesign every resolved residue instead, select `Entire resolved chain`.

A BoltzGen warning that designed residues also have a structure group asks you to
confirm the intended backbone-conditioned redesign; it is not automatically an
error.

## 13. Small-molecule-binder example

The pilot used:

- Ligand format: `Chemical Component Dictionary (CCD) code`
- CCD code: `D0G`
- Binder length: fixed at `123`
- Designs/budget: `20/20`

SMILES is also accepted. A structure-file input must be a complete PDB/mmCIF that
BoltzGen can parse. A hand-extracted, HETATM-only PDB often lacks the entity/polymer
metadata required by BoltzGen and can fail with an internal `IndexError`; BOLTRA now
detects the common bare-PDB form before validation and recommends CCD or SMILES.

For a ligand already defined in the PDB Chemical Component Dictionary, the CCD code
is normally the simplest and most reproducible input.

## 14. Nanobody example

Use `1T2P.pdb`, chain A, the common binding site, and `20/20` designs. Choose:

1. `BoltzGen template scaffolds discovered from the installation`.
2. Select the discovered `7eow` scaffold for the pilot configuration.

If no scaffold directory is found, complete section 6 and relaunch BOLTRA. Do not
hard-code a presumed list: BOLTRA discovers complete template pairs from the active
BoltzGen resources so additions or changes are visible.

The alternative custom-framework route asks for a structure, chain, and fixed-length
CDR regions to redesign.

## 15. Antibody example

Use `1T2P.pdb`, chain A, the same binding site, and `20/20` designs. Choose the
discovered `adalimumab.6cr1` Fab scaffold for the pilot configuration.

Custom antibody frameworks require distinct heavy and light chains and allow
fixed-length CDR selections. Always verify chain identities in BOLTRA's detected-chain
summary before proceeding.

## 16. Save now or run now

After validation, BOLTRA offers to save without running or begin the pipeline. A
saved project can be started later through **Open a project: inspect status or
resume**. Point this menu item at either the project directory or its
`boltzgen_output` subdirectory.

Do not enter an input PDB filename when prompted for the directory in which to create
a project. A valid example is `~/Documents/BOLTRA_Projects`, not
`~/Documents/target.pdb`.

## 17. Monitor a run safely

Use another terminal for lightweight monitoring:

```bash
nvidia-smi
free -h
df -h /path/to/project
```

The first large run should be a pilot. Scale gradually; thousands of designs can
require long GPU time, substantial disk space, and a correspondingly expensive
analysis stage. Keep browsers and other memory-heavy applications closed on a
RAM-constrained workstation. BOLTRA's preflight is a warning system, not a guarantee
that arbitrary scale will fit the machine.

## 18. Resume after interruption

Relaunch BOLTRA and choose menu option 2:

```bash
conda activate boltzgen
boltra
```

Enter the project or `boltzgen_output` directory. BOLTRA reports each stage as
complete or incomplete and offers to resume from the first incomplete stage while
preserving compatible files. It may repeat a partially completed stage.

If Linux killed a worker during analysis, reducing workers is safer than repeatedly
resuming with high concurrency. BOLTRA defaults to zero separate loader workers. Only
increase `BOLTRA_NUM_WORKERS` after confirming ample RAM and swap.

## 19. Perform post-design analysis

Choose menu option 3 and enter a completed `boltzgen_output` directory. BOLTRA prints
observed metric ranges, then requests:

1. one central metric and threshold;
2. zero or more supporting metrics and thresholds;
3. a unique analysis name.

Menu entries require option numbers. For example, select the number beside `Complex
ipTM`, then enter `0.13` when the program asks for its threshold. Typing `0.13` at the
menu-choice prompt is not interpreted as a threshold.

Use observed distributions, protocol-specific literature, controls, and study goals
to choose defensible thresholds. Do not treat pilot thresholds as universal.

## 20. Understand selection classes

Every design is assigned to one of four groups:

| Class | Meaning |
|---|---|
| Both | Passes the custom intersection and BoltzGen native filters |
| Custom only | Passes the user's thresholds but not all native filters |
| Native only | Passes native filters but not the custom intersection |
| Neither | Passes neither selection system |

For up to 50 designs, plots label all points with compact trailing numbers. For larger
runs, BOLTRA labels at most 30 prioritized selected points to avoid an unreadable
figure. Full identifiers always remain in the CSV and report.

## 21. Interpret completion counts correctly

Design filenames are zero-based: IDs `00` through `19` are 20 structures. This does
not imply that 19 structures were generated. However, table row counts are literal:
if BOLTRA reports `19/20 metric rows`, one generated design lacks an analysis record.
Inspect or resume the analysis stage rather than inferring completeness from the
largest filename.

## 22. Important outputs

At project level:

```text
design_specification.yaml
project_metadata.yaml
residue_mapping.csv                 # when conversion was used
boltra_execution_log.jsonl
inputs/
boltzgen_output/
```

Within `boltzgen_output`, BoltzGen writes intermediate designs, inverse-folded and
refolded structures, configuration files, metrics, final rankings, and its overview
PDF. Each BOLTRA post-analysis creates a separate directory containing CSV tables,
plots, captions, JSON settings, reports, and candidate structures.

## 23. Troubleshooting

### `EnvironmentNameNotFound: boltzgen`

Create the environment on that computer (section 3) and install BoltzGen there.

### `boltra: command not found`

Activate the correct environment, install from the directory containing
`pyproject.toml`, run `hash -r`, and compare `which python` with `which boltra`.

### Target structure not found

Enter the full filename, such as `~/Documents/Cyclo/1T2P.pdb`, rather than only its
directory.

### `NotADirectoryError` below a `.pdb` path

The project-root question expects a directory. Use a path such as
`~/Documents/BOLTRA_Projects`.

### Residue end is higher than chain length

Author residue numbers were probably supplied as BoltzGen sequential indices. Rerun
the guided flow and choose `PDB/mmCIF author residue numbers` so BOLTRA converts them.

### GitHub asks for a username while cloning a public repository

Use the downloadable source ZIP route in section 6. It avoids local Git credential
or HTTP transport issues.

### Browser crashes or `POOL BROKEN: A worker died`

This commonly indicates host-memory pressure, not interference between BOLTRA and the
browser. Keep the default zero-worker setting, provide swap, close memory-heavy
applications, and scale incrementally. Check the kernel log if permitted:

```bash
journalctl -k --since "30 minutes ago" | grep -Ei 'out of memory|oom|killed process'
```

### No nanobody or antibody scaffold directory

Follow section 6, verify `BOLTZGEN_ROOT`, open a new terminal or export it in the
current one, then relaunch BOLTRA.

## 24. Reproducibility checklist

Before archiving a study, retain:

- input structures and ligand identifiers;
- `design_specification.yaml` and project metadata;
- residue-mapping CSV files;
- scaffold snapshots for nanobody/antibody runs;
- complete `boltzgen_output` directories;
- BOLTRA analysis directories and thresholds;
- `boltra_execution_log.jsonl`;
- BOLTRA, BoltzGen, Python, CUDA/driver, and operating-system versions;
- random seeds or native configuration values when exposed;
- checksums and a description of any incomplete metric records.

Always cite BoltzGen as the underlying design engine, cite the archived BOLTRA release,
and describe BOLTRA as an orchestration and analysis layer rather than a generative
model.

---

## For questions
Please contact the corresponding authors:

Corresponding authors:
    Olanrewaju Ayodeji Durojaye;
    Rachid Daoud

Institution:
    Chemical and Biochemical Sciences, Green Process Engineering,
    University Mohammed VI Polytechnic,
    43150 Ben Guerir, Morocco

Email:
    olanrewaju.ayodeji-durojaye-ext@um6p.ma;
    rachid.daoud@um6p.ma
