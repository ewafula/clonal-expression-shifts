# Clonal expression shifts

Map gene-expression differences onto supplied clone and cellular tree partitions. The primary analysis uses a collapsed clone tree while retaining cell-level expression. An exploratory mode retains eligible cellular branches. Elastic Net, clone-balanced linear models and explicit contrasts describe complementary aspects of the data. The workflow does not infer reliable within-clone ancestry from expression alone. 


## Start here

- [Workflow roadmap](docs/WORKFLOW_ROADMAP.md): the complete analysis process, why each step exists, and how to interpret the evidence. Read this first for a conceptual overview.
- [Detailed reproducible tutorial](docs/TUTORIAL.md): all commands for a public artificial dataset, expected outputs, recovery and adapting your own data.
- [Results guide](docs/RESULTS_GUIDE.md): every result-table column, formulas, scientific limits and reviewed examples.
- [Reproducibility record](docs/REPRODUCIBILITY.md): tested versions, validation and environment-restoration limits.

Run from this directory. Start by checking the R environment:

```bash
Rscript 00_install_packages.R --check_only TRUE
```

If required packages are missing, `Rscript 00_install_packages.R` installs missing dependencies. Installation and later scripts must use the same R library. No new configuration framework is required.

## Repository layout

```text
clonal-expression-shifts/
├── README.md
├── 00_install_packages.R
├── 01_inspect_inputs.R
├── 02_prepare_tree.R
├── 03_prepare_expression.R
├── 04_fit_collapsed_cells.R
├── 05_validate_models.R
├── 06_fit_uncollapsed_cells.R
├── 07_run_gene_panel.R
├── 08_report_results.R
├── 09_make_example_inputs.R
├── docs/
│   ├── WORKFLOW_ROADMAP.md
│   ├── TUTORIAL.md
│   ├── RESULTS_GUIDE.md
│   └── REPRODUCIBILITY.md
├── inputs/                 # local; ignored except .gitkeep
└── results/                # generated; ignored except .gitkeep
```

`09_make_example_inputs.R` supports the tutorial; it generates inputs before Module 1 and is not an extra analysis after reporting. It reads no private data. Follow the tutorial's explicit synthetic paths rather than running Modules 1–3 with their current-data defaults.

## Inputs and ordinary usage

Provide a Seurat RDS with a known normalized layer, a headered TSV with `cell_id` and `clone_id`, and one explicitly rooted Newick cellular tree with matching cell tips. Clone/cell identifiers are essential. Optional sample IDs are descriptive; no sample/site grouping is required. The clone table currently requires these fixed column names; a metadata cell ID can be selected with `--cell_column`, including `rownames`.

Module 3 can use an original-format DEG table (`gene`, `clone_id`, `p_val`, `p_val_adj`, `avg_log2FC`, `pct.1`, `pct.2`) or an explicit gene list. DEG defaults retain the original filtering/ranking logic: adjusted p <0.05, absolute avg_log2FC ≥0.5, pct.1 ≥0.1, initially prepare the top 100. The priority score is a heuristic; upstream selection is not independent validation. Explicit genes permit expression preparation without a DEG table. The tutorial demonstrates that route.

| Script | Purpose | Default output |
|---|---|---|
| 00 | Check/install dependencies | Console and R session information |
| 01 | Inspect layers, metadata, cell overlap and clone counts | results/module01_input_report.txt |
| 02 | Audit monophyly; build both saved branch designs | results/module02.* |
| 03 | Rank candidates and prepare aligned expression | results/module03.* |
| 04 | Fit one collapsed gene; draw its selected-tree view | results/module04.* including .selected_tree.png |
| 05 | Validate the supplied pilot, either tree mode | results/module04_validation.* or results/module06_validation.* (set explicitly) |
| 06 | Fit one exploratory uncollapsed gene | results/module06.* |
| 07 | Run/validate a gene panel with checkpoint/resume | results/module07/ and gene subfolders |
| 08 | Interpret saved panel results without refitting | results/module08/ |
| 09 | Generate entirely artificial tutorial inputs | inputs/tutorial/ |

For optparse-based scripts use `--help` to see arguments. Module 1 also supplies `--help`; installer usage is documented at the top of script 00. Input/output paths can be overridden. Do not edit the original inputs or saved estimates to force alignment.

## Fresh-run fitting and validation sequence

After Modules 1–3 are complete, the current-data default CLDN16 pilot is:

```bash
Rscript 04_fit_collapsed_cells.R
Rscript 05_validate_models.R \
  --pilot results/module04.models.rds \
  --output_prefix results/module04_validation

Rscript 06_fit_uncollapsed_cells.R
Rscript 05_validate_models.R \
  --pilot results/module06.models.rds \
  --output_prefix results/module06_validation
```

Validation filenames identify the fitted analysis: `module04_validation.*` for the collapsed pilot and `module06_validation.*` for the uncollapsed pilot. The script 05 identifies the shared validator, not a separate output group. Always supply the corresponding `--pilot` and `--output_prefix` together. Script 05's defaults match the collapsed pilot.

Script **05 is shared**: `--pilot` selects which fitted gene/tree mode to validate. It is not necessary to repeat a completed collapsed validation to validate uncollapsed models. Separate prefixes preserve both sets of outputs. For other prepared genes/inputs, pass explicit arguments as shown in the tutorial.

## Gene panels and output organization

```bash
Rscript 07_run_gene_panel.R
```

With neither gene argument, this selects up to five different prepared genes: CLDN16 (or the first available usable gene), then strongest remaining clone-mean range, sparsest detection, broadest remaining detection and smallest remaining clone-mean range. The current real-data selection was CLDN16, S100A10, PAMR1, S100A4 and LINC03061. Other datasets select different genes. These are descriptive patterns, not independent null controls.

For chosen genes use one of:

```bash
Rscript 07_run_gene_panel.R --genes CLDN16,PAMR1,S100A4

Rscript 07_run_gene_panel.R \
  --gene_list inputs/panel_genes.txt \
  --output_dir results/module07_custom
```

The text file is headerless, one exact prepared gene per line. Blank lines are ignored, duplicates removed in first-occurrence order, and surrounding whitespace trimmed. `--genes` and `--gene_list` are alternatives. Explicit lists have no ten-gene cap, but every gene must already be finite/variable in the Module 3 output. Larger panels are supported for independent use after development; they are not required for this workflow example.

```text
results/module07/
├── module07.panel.tsv
├── module07.status.tsv
├── module07.summary.tsv
├── module07.checkpoint.rds
├── module07.report.txt
├── CLDN16/
│   ├── CLDN16_collapsed_attempt1.*
│   ├── CLDN16_collapsed_validation_attempt1.*
│   ├── CLDN16_uncollapsed_attempt1.*
│   └── CLDN16_uncollapsed_validation_attempt1.*
└── ...other gene folders...
```

Resume the identical panel command with `--resume`; completed stages are reused and failed/interrupted stages retry under new attempt suffixes. Changed inputs, analysis scripts, settings or gene order require a new output directory. A checkpoint is not a backup of lost outputs. Ordinary gene symbols are preserved in folder names; unsafe path characters become underscores and conflicting folder names are rejected.

## Integrated interpretation

```bash
Rscript 08_report_results.R
```

The main output is `results/module08/module08.report.md`, accompanied by gene/fit/branch/contrast/source TSVs, an overview PNG and session text. It computes interpretations from saved results and verifies matching inputs and source tables. It does not refit models or reload the original Seurat RDS. Incomplete pairs are labelled explicitly. Keep results/module07 alongside module 08 for relative links.

Module 8 supports other compatible workflow datasets:

```bash
Rscript 08_report_results.R \
  --panel_dir results/module07_custom \
  --output_dir results/module08_custom
```

Matching Modules 2–3 prepared inputs and referenced per-gene outputs must remain available. This is not an importer for arbitrary third-party results. Validation covers the reviewed five-gene dataset, an artificial example and consistency/failure tests; it is not validation of every assay, tree or biological sampling design.

## Interpretation and sharing

Clone-balanced apparent training fit, tuning CV, fixed-tree cell-subsampling stability and biological replication are different evidence. Reported coefficient/contrast units follow the prepared response, not automatically DEG log2FC. Working IC uses a different unweighted Gaussian ML objective and does not correct cell dependence, confounding or post-selection. Uncollapsed candidate searches are bounded; no forced consensus or support threshold is applied. Some inputs and derived outputs stay outside the GitHub repository (large files and RDS objects). The s

## Consistent outputs for both tree modes

Each fitting script (04 collapsed or 06 uncollapsed) writes the same eight output types: report, coefficients, contrasts, fit metrics, cell predictions, saved models, comparison PNG and selected-tree PNG. Each comparison has observed expression by clone, observed-versus-training-fitted expression, and the largest descriptive contrasts. The contrast definitions and tree resolution remain specific to the mode.

Shared validator 05 writes the same seven validation outputs for either mode, including stability.png. Module 7 calls these same scripts and requires the complete output sets for every gene/mode. Every completed gene/mode therefore has comparison.png, selected_tree.png and validation stability.png automatically. Module 8 requires and links these figures.

For a clean rerun, archive previous outputs outside the active results path, preserve real inputs, and follow the fresh commands in [TUTORIAL](docs/TUTORIAL.md). Old panel checkpoints are not resumed across changed script versions. The current implementation uses explicit tree_mode metadata in both saved pilots; use matching current scripts 04–08 together.
