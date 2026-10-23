# Reproducible tutorial

For the complete analysis process and the role of every module, start with [WORKFLOW_ROADMAP.md](WORKFLOW_ROADMAP.md).

This tutorial runs the standalone workflow without any private data. Run commands from the `clonal-expression-shifts/` directory in a terminal. Use the existing root scripts; no configuration framework is needed. The artificial example exercises both tree modes, all three analysis methods, shared validation, a five-gene panel and integrated reporting.

For your own data, see the adaptation section below. For column definitions, formulas and worked interpretations, read [RESULTS_GUIDE.md](RESULTS_GUIDE.md). Environment details and limits are in [REPRODUCIBILITY.md](REPRODUCIBILITY.md).

## 1. Check dependencies

```bash
Rscript 00_install_packages.R --check_only TRUE
```

If required packages are missing, install them in your R environment:

```bash
Rscript 00_install_packages.R
```

Then repeat the check. Installation needs a working package repository/network and, where source builds are needed, system build dependencies. The installer installs missing packages, not a frozen historical environment. It does not downgrade already installed packages to the tested versions.

For a custom library, make the library visible to every subsequent R process, not only to the installer. One terminal-session example is:

```bash
mkdir -p local_R_library
export R_LIBS_USER="$PWD/local_R_library"
Rscript 00_install_packages.R --lib local_R_library
Rscript 00_install_packages.R --check_only TRUE
```

Keep machine-specific libraries out of Git. The custom-library setting above is optional; use your site's normal R environment when available. Do not replace a shared server library merely to match this tutorial.

## 2. Generate shareable artificial inputs

```bash
Rscript 09_make_example_inputs.R
```

This creates `inputs/tutorial/` with:

| File | Content |
|---|---|
| `expression.rds` | Small Seurat object: 128 artificial cells, 25 genes, counts and a normalized data layer. |
| `clones.tsv` | Headered `cell_id` / `clone_id` mapping; four clones A–D with 32 cells each. |
| `tree.nwk` | Artificial rooted binary cellular topology with monophyletic clones. |
| `genes.txt` | Five exact gene identifiers, one per line without a header. |
| `README.txt` | Generation seed, scale and interpretation limits. |

The generator reads no previous expression data, clinical records or analysis results. It uses seed 2026. `data` is natural-log `log1p(count / cell library total * 10000)`. All tree topology is constructed; it is not inferred biological ancestry. No branch support values or DEG p-values are fabricated.

The genes illustrate strong Clone A expression (`DEMO-CLONE-A`), sparse Clone C expression (`DEMO-CLONE-C-SPARSE`), a designed subdivision within B (`DEMO-WITHIN-B`), broad expression (`DEMO-BROAD`) and weaker expression (`DEMO-WEAK`). Sampling noise and normalization affect observed values. Broad/weak expression does not make a gene a validated null control. There are 20 additional background features used in library normalization.

Existing inputs are protected. To regenerate separately, use `--output_dir inputs/tutorial_repeat` and change downstream input paths accordingly. Feature names use hyphens because Seurat replaces underscores in feature names.

## 3. Inspect inputs: Module 1

```bash
Rscript 01_inspect_inputs.R \
  --seurat_object inputs/tutorial/expression.rds \
  --concordant_clones inputs/tutorial/clones.tsv \
  --clonal_tree inputs/tutorial/tree.nwk \
  --output results/tutorial/module01_input_report.txt
```

Read the report before proceeding. Expect 128 shared cells, no unmatched clone/tree/assay cells, RNA `counts` and `data` layers, and four clone groups of 32. `CLDN16 present: FALSE` is expected: the generic example uses artificial gene names. Sample/site grouping is not required. Cell barcodes and clone mapping are essential; `sample_id` is optional descriptive metadata.

## 4. Prepare both tree designs: Module 2

```bash
Rscript 02_prepare_tree.R \
  --concordant_clones inputs/tutorial/clones.tsv \
  --clonal_tree inputs/tutorial/tree.nwk \
  --min_descendant_cells 8 \
  --output_prefix results/tutorial/module02
```

The small example uses a descendant-cell threshold of eight so within-clone partitions are eligible. This deliberately differs from the default 25 used for larger trees. Do not silently change scientific filters to recover preferred results on real data.

Expect all four clones monophyletic, five canonical collapsed predictors spanning three independent directions, and 29 canonical uncollapsed predictors with rank 16 including the intercept. Branch names belong to this saved tree; they do not match the private example's IDs. Inspect `.clone_checks.tsv`, `.branches.tsv`, `.report.txt` and the collapsed tree figure. A monophyletic clone permits collapse but does not validate root or ancestral resolution.

## 5. Prepare explicit genes: Module 3

```bash
Rscript 03_prepare_expression.R \
  --seurat_object inputs/tutorial/expression.rds \
  --tree_designs results/tutorial/module02.tree_designs.rds \
  --deg_file inputs/tutorial/not_supplied.tsv \
  --gene_list inputs/tutorial/genes.txt \
  --expression_scale 'Natural-log log1p(count/library total*10000); synthetic tutorial' \
  --output_prefix results/tutorial/module03
```

`not_supplied.tsv` is intentionally absent: an explicit gene list allows preparation without a DEG table. Do not create this file. If a file already exists at that path, choose a different absent filename. This tutorial exercises the explicit-gene route, not the original DEG filtering/ranking route. Real DEG selection is described below.

Expect five ready genes × 128 cells, in exact Module 2 order. Candidate-ranking tables are empty because no DEG table was supplied; no filtered-DEG file is produced in this mode. The expression matrix is sparse and preserves observed zeros. Module 3 performs no imputation or normalization. Its response-scale string records provenance; it does not transform data.

## 6. Fit and validate the collapsed pilot: Modules 4–5

```bash
Rscript 04_fit_collapsed_cells.R \
  --expression results/tutorial/module03.expression.rds \
  --tree_designs results/tutorial/module02.tree_designs.rds \
  --gene DEMO-CLONE-A \
  --output_prefix results/tutorial/module04

Rscript 05_validate_models.R \
  --pilot results/tutorial/module04.models.rds \
  --expression results/tutorial/module03.expression.rds \
  --tree_designs results/tutorial/module02.tree_designs.rds \
  --output_prefix results/tutorial/module04_validation
```

Module 4 fits EN, weighted linear models and equal-clone contrasts. It always creates both `.comparison.png` and `.selected_tree.png`. The comparison panels show observed expression by clone, observed versus fitted expression for EN (otherwise linear or clone means), and the largest available descriptive contrasts. The selected tree maps EN contributions to the clone topology by descendant clone sets, with original branch IDs/coefficients and clone cell counts. No selected EN terms or an omitted EN method is labelled explicitly. Inverse-clone-size weights give equal total weight to each clone; they are balancing weights, not known precision weights. EN uses alpha 0.5, five clone-stratified CV folds and `lambda.1se` by default. CV tunes interpolation among cells in represented clones. Linear basis coefficients, simultaneous penalized coefficients and group contrasts have different reference definitions and should not be called interchangeable fold changes.

Module 5 performs 100 clone-stratified 80% subsamples without replacement and retunes EN each time. Inspect method status counts and per-term denominators. Quantiles are subsampling ranges, not confidence intervals; selection fractions are not p-values. A working IC comparison separately refits ordinary unweighted iid Gaussian ML models on the same gene/cells. It does not score the penalized glmnet loss or treat balancing weights as precision.

## 7. Fit and validate the uncollapsed pilot: Modules 6 and 5

```bash
Rscript 06_fit_uncollapsed_cells.R \
  --expression results/tutorial/module03.expression.rds \
  --tree_designs results/tutorial/module02.tree_designs.rds \
  --gene DEMO-CLONE-A \
  --output_prefix results/tutorial/module06

Rscript 05_validate_models.R \
  --pilot results/tutorial/module06.models.rds \
  --expression results/tutorial/module03.expression.rds \
  --tree_designs results/tutorial/module02.tree_designs.rds \
  --output_prefix results/tutorial/module06_validation
```

The tutorial collapsed validation files use `results/tutorial/module04_validation.*`; uncollapsed validation files use `results/tutorial/module06_validation.*`. The ovarian example uses the same basenames directly under results/. Script 05 is a shared operation, not a standalone output basename. Keep `--pilot` and its `--output_prefix` explicit as shown.

Script **05 is the shared validator**. Its `--pilot` determines which gene and tree mode are validated. This second invocation evaluates Module 6; it does not repeat collapsed validation. Separate prefixes keep both sets of results. The canonical four-step order is **04 → 05 on module04 → 06 → 05 on module06**.

Module6 writes the same eight output types as Module4, including the same comparison-panel structure and a selected-tree PNG. Validator 05 writes the same seven output types in either mode, including its stability PNG.

The cellular view is exploratory. Black background branches show the supplied topology, red/blue show positive/negative EN contributions, and the tip strip identifies clones. Within-clone/sister contrasts require enough cells on both sides; some requests are unavailable by design. A completed contrast iteration does not mean every requested contrast was estimable. The flexible linear basis can lose rank in subsamples; unavailable coefficients must not be replaced with zero.

## 8. Run the five-gene panel: Module 7

```bash
Rscript 07_run_gene_panel.R \
  --expression results/tutorial/module03.expression.rds \
  --tree_designs results/tutorial/module02.tree_designs.rds \
  --gene_list inputs/tutorial/genes.txt \
  --output_dir results/tutorial/module07
```

This independently runs the existing 04/06 fitters and shared validator 05 for each specified gene, with 100 subsamples per mode. It deliberately produces panel-specific fits rather than importing earlier pilot files. Expect ten completed validations and 20 successful child stages, with a directory per gene. Each gene/mode has the same eight fit files plus seven validation files and separate child console logs. Its three figures are comparison.png, selected_tree.png and validation stability.png; no separate plotting command is required. `module07.panel.tsv` identifies the panel and reasons; `.status.tsv` records stages and elapsed time; `.summary.tsv` combines metrics; `.checkpoint.rds` enables resume.

For a short explicit list, use `--genes DEMO-CLONE-A,DEMO-WITHIN-B` instead of `--gene_list`. Do not supply both. A list file is headerless, one exact gene per line. Blank lines and duplicates are removed; case-sensitive IDs must already be prepared. Explicit lists have no ten-gene cap. Larger lists are for later independent analysis, not required to validate this tutorial.

If neither option is given, Module 7 selects up to five distinct prepared genes in order: CLDN16 if available (otherwise the first usable gene), then largest remaining clone-mean range, lowest detection fraction, highest remaining detection fraction and smallest remaining clone-mean range. The gene is removed after each choice. Different datasets select different genes. The real example's five genes are not required for generic operation; details and caveats are in RESULTS_GUIDE.

## 9. Generate the integrated report: Module 8

```bash
Rscript 08_report_results.R \
  --panel_dir results/tutorial/module07 \
  --output_dir results/tutorial/module08
```

Open `results/tutorial/module08/module08.report.md` in a Markdown viewer. Read the panel overview, per-gene clone expression, cross-mode comparisons, EN contributions, contrast support and working IC interpretation. Eight output files include five combined/provenance TSVs, the Markdown report, overview PNG and report-generation session record. No model is refitted. All newly generated table columns are defined in RESULTS_GUIDE.

The report is generic for compatible workflow outputs. Gene/clone names and numerical findings are calculated from each saved run; no CLDN16 conclusion is required. It needs the matching Modules 2–3 inputs and Module 7 outputs with their referenced paths available. It cannot interpret arbitrary third-party model results. Change `--panel_dir` and `--output_dir` for another completed dataset. Keep panel/report folders together so relative image/table links resolve; checkpoint references to prepared inputs must remain available from the working directory. Moving a report alone does not make it self-contained.

The current validation covers the original five-gene dataset plus this independent artificial dataset and failure/consistency tests. This is not evidence that every assay, topology, sample structure or preprocessing pipeline has been validated. Independent biological replication, sample confounding and uncertain real genealogy remain outside what this tutorial establishes.

## 10. Verify completion

Review outputs rather than relying only on exit messages. The following command checks the basic tutorial invariants from the workflow directory:

```bash
Rscript -e 'library(Matrix); t<-readRDS("results/tutorial/module02.tree_designs.rds"); e<-readRDS("results/tutorial/module03.expression.rds"); s<-read.delim("results/tutorial/module08/module08.gene_summary.tsv"); stopifnot(all(t$clone_checks$monophyletic), ncol(e$expression)==128, nrow(e$expression)==5, nrow(s)==10, all(s$status=="complete")); cat("Tutorial invariants passed\n")'
```

Exact branch coefficients, lambda values and tied IC representatives can differ across package/platform versions. Validate identities, dimensions, estimability and numerical consistency, not a screenshot or an exact rank among equivalent branch parameterizations. A designed within-B expression pattern should be inspected against direct within-clone means; topology is known by construction here, unlike biological data.

## Resume, repeat and recover

- **Interrupted Module 7:** repeat the identical panel command with `--resume`. Completed stages with required files present are reused. Failed/interrupted stages use another attempt suffix and preserve old files. Resume requires matching inputs, analysis scripts, gene order and settings.
- **Changed panel/settings/analysis code:** choose a new `--output_dir`, for example `results/tutorial_panel_repeat`; then point Module 8 there. Existing results are not overwritten.
- **Intentional pilot repeat:** use a new `--output_prefix`, for example `results/tutorial/module04_repeat`, and pass its `.models.rds` to validator 05 with another new validation prefix. A finished collapsed validation need not be repeated to validate uncollapsed models.
- **Report text update:** use a new report output directory. No fitting rerun is needed. Reports retain source hashes at generation time; later documentation/code changes do not invalidate the original numerical analysis.
- **Completed files lost:** restore the saved run before resuming; a checkpoint is not a backup of model files.
- **Tree or expression changed:** rerun preparation and downstream analysis under consistent new prefixes. Old hashes intentionally prevent mixing runs. Node/edge numbers from different saved trees are not interchangeable.

## Adapt to your own data

1. Provide a Seurat RDS with a known normalized assay/layer; a headered tab-separated clone table with `cell_id` and `clone_id`; and one explicitly rooted Newick tree whose tips match clone-table IDs. Custom metadata cell IDs can be selected with `--cell_column`; clone-table field names remain fixed. IDs must be unique and nonempty. Sample/site labels are optional descriptive information, not required analysis groups.
2. Run Module 1 with your input paths. Use `--cell_column rownames` if metadata has no cell ID field; use the same convention in Module 3. Resolve overlap problems rather than relabelling cells by position. An object may contain additional cells outside the retained tree/clone analysis.
3. Run Module 2 with your tree/clone paths and reviewed size filters. Non-monophyletic clones cannot be treated as a valid collapsed clone tree. The code does not infer a trustworthy root or interpret all raw node labels as support probabilities.
4. Choose explicit genes (`--genes` or `--gene_list`) or an original-format DEG table. DEG columns must include `gene`, `clone_id`, `p_val`, `p_val_adj`, `avg_log2FC`, `pct.1`, `pct.2`. Default filters are adjusted p <0.05, absolute avg_log2FC ≥0.5 and pct.1 ≥0.1. Ranking is the reused historical heuristic, not a new test. Upstream DEG selection on these same cells is not independent validation of downstream findings.
5. Module 3 reads the exact chosen assay/layer, preserves zeros and checks tree order. It does not automatically normalize counts, join split layers or impute. A `data` layer name alone does not prove the scale. Describe the actual scale with `--expression_scale` when known; do not label all coefficients log2FC.
6. Apply the pilot four-command sequence, then the chosen panel and report commands. Use new prefixes/directories to isolate datasets. Different datasets may have different clone counts, branch rank, eligible contrasts or successful subsamples; honest NA/non-estimable statuses are expected.

For example, real-data DEG selection (replace paths with your own) is:

```bash
Rscript 03_prepare_expression.R \
  --seurat_object inputs/my_expression.rds \
  --tree_designs results/my_module02.tree_designs.rds \
  --deg_file inputs/my_clone_deg.tsv.gz \
  --top_n 100 \
  --output_prefix results/my_module03
```

To run all 100 prepared genes later, supply a text list of those exact genes to Module 7; omitting `--genes`/`--gene_list` selects the small automatic panel, not every prepared gene. Use `--expression results/my_module03.expression.rds`, matching `--tree_designs`, and a new panel output directory. A larger run is not automatically a stronger biological experiment.

## Common problems

| Message/symptom | Action |
|---|---|
| Missing R package | Run 00 in the same library/environment used by later Rscript processes. |
| Output already exists | Use a new prefix/directory; use --resume only for the same Module 7 run. |
| Gene not prepared | Add it during Module 3 preparation in a new consistent run before fitting. |
| Counts or split layers unsupported | Prepare/review normalized data externally; choose its exact layer. The workflow does not guess normalization. |
| Unrooted tree/non-monophyletic clone | Review tree provenance/clone mapping; do not force a biological collapse from expression alone. |
| Clone labels all numeric/unsorted | Treat IDs as labels. Read saved dictionaries rather than assuming numeric order is evolutionary order. |
| Non-estimable resamples | Read individual status messages and term denominators; do not replace estimates with zero. |
| Different AIC and BIC winners | They penalize complexity differently; inspect equivalent spaces and working-likelihood assumptions. |
| Broken report links | Keep referenced result folders/prepared input paths available; regenerate a report in a new directory if paths changed. |
| Module 8 provenance/TSV mismatch | Restore consistent original inputs/results. Do not edit estimates to satisfy a hash check. |

## One clean implementation for fresh runs

Use matching current scripts 04–08 together. Both pilots explicitly save their tree mode. Shared 05 reads that mode; Module 7 uses the same complete fit/validation file contracts; Module 8 verifies and links all three figures.

| Stage | Output suffixes in both modes |
|---|---|
| Fit 04 or 06 | `.report.txt`, `.coefficients.tsv`, `.contrasts.tsv`, `.fit_metrics.tsv`, `.cell_predictions.tsv.gz`, `.models.rds`, `.comparison.png`, `.selected_tree.png` |
| Validate 05 | `.report.txt`, `.stability.tsv`, `.resamples.tsv.gz`, `.iteration_status.tsv`, `.model_comparison.tsv`, `.validation.rds`, `.stability.png` |

The comparison layout is common but the contrasts retain their actual estimands: equal-clone comparisons in collapsed mode, within-clone/sister cell comparisons in uncollapsed mode. Tree drawings likewise show clone versus cellular topology. Uniform output names do not imply identical biological resolution or historical certainty.

### Restart the ovarian-cancer example cleanly

Archive the old results directory before a new full run, and keep the original input files intact. Do not mix old model/checkpoint files with the new implementation. If using the existing ovarian-data filenames/defaults, run these commands from the workflow directory:

```bash
Rscript 00_install_packages.R --check_only TRUE
Rscript 01_inspect_inputs.R
Rscript 02_prepare_tree.R
Rscript 03_prepare_expression.R
Rscript 04_fit_collapsed_cells.R
Rscript 05_validate_models.R \
  --pilot results/module04.models.rds \
  --output_prefix results/module04_validation
Rscript 06_fit_uncollapsed_cells.R
Rscript 05_validate_models.R \
  --pilot results/module06.models.rds \
  --output_prefix results/module06_validation
Rscript 07_run_gene_panel.R
Rscript 08_report_results.R
```

Read each module's report/status before continuing, and stop on failure. These defaults are for the established ovarian input filenames and CLDN16 pilot, not arbitrary new datasets. For different datasets use explicit input/gene paths as described above.

### Restart the artificial tutorial cleanly

Archive old `inputs/tutorial/` and `results/tutorial/` if they exist; preserve all real-data inputs and results. Then run the numbered sections1–10 of this tutorial from the beginning. The generator creates a fresh artificial input set and every figure is produced during its normal fit/validation stage. The synthetic output paths are separate from the ovarian run. `--resume` remains supported only for an interrupted identical Module 7 run with unchanged inputs, scripts, gene order and settings; it is not used to migrate old runs after this cleanup.

Tutorial validation checkpoint October 19 2026: the returned server example completed both pilots and all ten panel gene/mode validations with 100 resamples each; all methods completed successfully and report/provenance checks passed. No corrective rerun is required.
