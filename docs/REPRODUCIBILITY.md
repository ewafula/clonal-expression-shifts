# Environment, reproducibility and validation record

The workflow uses numbered R scripts without a configuration framework or an automatically activated project library. [TUTORIAL.md](TUTORIAL.md) gives the complete public synthetic example. This document distinguishes a reproducible recipe, a pinned version reference and a portable environment restore.

## Tested environment and version reference

The local tutorial was exercised using R 4.6.0 on macOS and 4.6.1 HPC server. Directly relevant versions on  both macOSa and HPC server were Seurat 5.5.1, SeuratObject 5.4.0, Matrix 1.7.6, glmnet 5.0, ape 5.8.1, optparse 1.8.2, dplyr 1.2.1 and readr 2.2.0; consult the exact inventory below if a package version differs from this prose. Matching results across these environments does not establish platform-independent exact coefficients.

```bash
J72MQ09JH7:clonal-expression-shifts wafulae$ Rscript 00_install_packages.R --check_only TRUE

========== R PACKAGE INSTALLATION CHECK ==========
[INFO] CRAN mirror: https://cloud.r-project.org

[INFO] No --lib path provided. Using current .libPaths():
[1] "/Library/Frameworks/R.framework/Versions/4.6/Resources/library"

========== CHECKING REQUIRED PACKAGES ==========
[OK] tidyverse is already installed.
[OK] optparse is already installed.
[OK] Seurat is already installed.
[OK] SeuratObject is already installed.
[OK] glmnet is already installed.
[OK] ape is already installed.
[OK] Matrix is already installed.
[OK] phangorn is already installed.
[OK] igraph is already installed.
[OK] lars is already installed.

[INFO] Skipping optional Seurat extensions.
       To install BPCells, presto, and glmGamPoi, run:
       Rscript 00_install_packages.R --install_seurat_extensions TRUE

========== FINAL PACKAGE CHECK ==========
      package installed
    tidyverse      TRUE
     optparse      TRUE
       Seurat      TRUE
 SeuratObject      TRUE
       glmnet      TRUE
          ape      TRUE
       Matrix      TRUE
     phangorn      TRUE
       igraph      TRUE
         lars      TRUE

[SUCCESS] All required packages are installed and available.

========== LOAD TEST ==========
[LOAD] Testing library(tidyverse)
[LOAD] Testing library(optparse)
[LOAD] Testing library(Seurat)
[LOAD] Testing library(SeuratObject)
[LOAD] Testing library(glmnet)
[LOAD] Testing library(ape)
[LOAD] Testing library(Matrix)
[LOAD] Testing library(phangorn)
[LOAD] Testing library(igraph)
[LOAD] Testing library(lars)
Warning messages:
1: package ‘Seurat’ was built under R version 4.6.1 
2: package ‘sp’ was built under R version 4.6.1 
3: package ‘Matrix’ was built under R version 4.6.1 
4: package ‘igraph’ was built under R version 4.6.1 
[SUCCESS] All required packages loaded successfully.

========== SESSION INFO ==========
R version 4.6.0 (2026-04-24)
Platform: aarch64-apple-darwin23
Running under: macOS Tahoe 26.6

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRblas.0.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] lars_1.3           igraph_2.3.3       phangorn_2.12.1    ape_5.8-1         
 [5] glmnet_5.0         Matrix_1.7-6       Seurat_5.5.1       SeuratObject_5.4.0
 [9] sp_2.2-3           optparse_1.8.2     lubridate_1.9.5    forcats_1.0.1     
[13] stringr_1.6.0      dplyr_1.2.1        purrr_1.2.2        readr_2.2.0       
[17] tidyr_1.3.2        tibble_3.3.1       ggplot2_4.0.3      tidyverse_2.0.0   

loaded via a namespace (and not attached):
  [1] deldir_2.0-4           pbapply_1.7-4          gridExtra_2.3.1       
  [4] rlang_1.3.0            magrittr_2.0.5         RcppAnnoy_0.0.23      
  [7] otel_0.2.0             matrixStats_1.5.0      ggridges_0.5.7        
 [10] compiler_4.6.0         spatstat.geom_3.8-2    png_0.1-9             
 [13] vctrs_0.7.3            reshape2_1.4.5         quadprog_1.5-8        
 [16] shape_1.4.6.1          pkgconfig_2.0.3        fastmap_1.2.0         
 [19] promises_1.5.0         tzdb_0.5.0             jsonlite_2.0.0        
 [22] goftest_1.2-3          later_1.4.8            spatstat.utils_3.2-4  
 [25] irlba_2.3.7            parallel_4.6.0         cluster_2.1.8.3       
 [28] R6_2.6.1               ica_1.0-3              stringi_1.8.7         
 [31] RColorBrewer_1.1-3     spatstat.data_3.1-9    reticulate_1.46.0     
 [34] parallelly_1.48.0      spatstat.univar_3.2-0  lmtest_0.9-40         
 [37] scattermore_1.2        iterators_1.0.14       Rcpp_1.1.2            
 [40] tensor_1.5.1           future.apply_1.20.2    zoo_1.9-0             
 [43] sctransform_0.4.3      timechange_0.4.0       httpuv_1.6.17         
 [46] splines_4.6.0          tidyselect_1.2.1       abind_1.4-8           
 [49] spatstat.random_3.5-1  codetools_0.2-20       miniUI_0.1.2          
 [52] spatstat.explore_3.8-2 listenv_1.0.0          lattice_0.22-9        
 [55] plyr_1.8.9             withr_3.0.3            shiny_1.14.0          
 [58] S7_0.2.2               ROCR_1.0-12            Rtsne_0.17            
 [61] future_1.75.0          fastDummies_1.7.6      survival_3.8-9        
 [64] polyclip_1.10-7        fitdistrplus_1.2-6     pillar_1.11.1         
 [67] KernSmooth_2.23-26     foreach_1.5.2          plotly_4.12.1         
 [70] generics_0.1.4         RcppHNSW_0.7.0         hms_1.1.4             
 [73] scales_1.4.0           globals_0.19.1         xtable_1.8-8          
 [76] glue_1.8.1             tools_4.6.0            data.table_1.18.4     
 [79] RSpectra_0.16-2        RANN_2.6.2             dotCall64_1.2         
 [82] fastmatch_1.1-8        cowplot_1.2.0          grid_4.6.0            
 [85] nlme_3.1-170           patchwork_1.3.2        cli_3.6.6             
 [88] spatstat.sparse_3.2-0  spam_2.11-4            viridisLite_0.4.3     
 [91] uwot_0.2.4             gtable_0.3.6           digest_0.6.39         
 [94] progressr_1.0.0        ggrepel_0.9.8          htmlwidgets_1.6.4     
 [97] farver_2.1.2           htmltools_0.5.9        lifecycle_1.0.5       
[100] httr_1.4.8             mime_0.13              MASS_7.3-66           

========== INSTALLATION CHECK COMPLETE ==========
```

A mismatch signals a different environment; it does not by itself show a scientific error. Review the change and rerun the synthetic checks in that environment. Base-package version matching also depends on the R release. Run the script 01 as show above to if required dependencies are installed. Also read `module08.session.txt` and the `sessionInfo` in earlier RDS outputs to determine the environment that actually generated a particular run, rather than substituting this installed check.

## What is deterministic

- The artificial generator uses a fixed seed (default 2026), fixed clone/cell/gene identities and a constructed rooted topology. Its input provenance is independent of private data.
- Model and panel defaults use seed 12345. Validation records sampling indices, clone counts, folds and tuning choices. Weights and rank checks are explicit.
- Prepared expression records its associated tree content hash. Pilots and validators record input hashes; reports verify compatible inputs, saved responses, linked tables and source hashes.
- Same-version same-platform reruns should reproduce the important numerical results, but timestamps/session records are intentionally run-specific. Identical biological interpretation does not require byte-identical report files.
- Different R versions, package implementations, numerical libraries, character-collation locales or near-tied lambda/IC scores can change numerical details. IC representatives of equivalent prediction spaces can differ without a different fitted space.

## Meaningful validation performed

The artificial tutorial exercises input inspection, clone monophyly/design construction, explicit-gene expression preparation, collapsed fitting/validation, uncollapsed fitting/validation, five-gene panel fitting/validation and report generation. The complete example uses 100 resamples in each validation, not a three-replicate production claim. Check exact cell/clone identities, matrix dimensions, fixed sampling membership, response matching, valid denominator sums, analytic group contrasts, fit-summary consistency and source links. Do not make exact selected-branch counts a universal pass condition: EN selection and nested aliases depend on the saved candidate design. Note that artificial tutorial exercises data is generated from source and requires no inputs or results. `inputs/`, `results/` and optional local package libraries are ignored. 

The previously reviewed five-gene real-data run separately exercises the original DEG selection path, larger cellular topology and non-estimable subsamples. Synthetic tutorial operation does not validate biological root choice, real-tree support, patient independence, expression scale provenance for new datasets or unmeasured confounding. Tutorial genes with broad/weak expression are descriptive examples, not independent null controls, and no power/type-I-error study is claimed.

Resume, changed-settings rejection, failed-stage retry, missing/mismatched provenance, incomplete report handling and overwrite protection have been checked during development. These checks live in temporary development files; users do not need a permanent test framework to execute the workflow. The tutorial includes a compact invariant check and documents expected statuses.

## Scientific interpretation limits

Collapsed analysis remains primary. Uncollapsed expression subdivisions remain exploratory on supplied biological topologies. EN, linear bases and explicit contrasts target different estimands. Training R-squared, cell-stratified tuning CV, cell-subsampling stability and independent biological validation are not interchangeable. Working AIC/BIC assumes independent homoskedastic Gaussian cell errors, scores unweighted ML refits, counts design rank plus one variance parameter, and remains conditional on the evaluated candidate set. It does not correct dependent cells, model-search selection or ancestry uncertainty. No panel-wide confirmatory significance claims are produced.

Module 8 is generic for compatible Modules 2–3/7 outputs; it is not an importer for unrelated analyses results. The artificial dataset tests operation beyond the original dataset, not all possible gene names, clone counts, assay structures or input conventions. Retain matching prepared inputs and original folder paths; copying a Markdown report alone does not make an independent runnable analysis.

## Shared fit, validation and figure contract

Scripts 04/06 always produce the same eight fitting file types, including comparison.png and selected_tree.png. Script 05 always produces the same seven validation file types for either declared tree mode, including stability.png. Script 07 calls those scripts directly, checks every expected file, and preserves per-gene/mode outputs; it does not maintain another plotter. Report 08 requires all three figures for completed pairs and links both mode-specific selected trees. There is no independent plotting/repair execution path.

Both pilots explicitly save tree_mode, so validators/reporting reject missing or inconsistent mode metadata rather than infer a layout. Start a fresh run with matching current scripts 04–08; do not resume checkpoints from earlier script versions. Numerical fitting/tuning/resampling objectives are unchanged by this output cleanup. Tests compare estimates to earlier reviewed results as well as checking files, method-specific empty states, plots and provenance.

The clean implementation was checked with the complete artificial tutorial and the ovarian prepared-input analysis (both pilots and five-gene panels, 100 resamples). Coefficients, contrasts, fit metrics, IC and EN selection summaries agree with earlier reviewed results within numerical tolerance. Some raw linear sign fractions differ between local R v4.6.0 and server R v4.6.1 only where underlying coefficients are at numerical-zero scale; per-sample differences were checked.

## Validation output names

The shared script 05 is run after either fitter. Use `--pilot results/module04.models.rds --output_prefix results/module04_validation` for collapsed validation and `--pilot results/module06.models.rds --output_prefix results/module06_validation` for uncollapsed validation. Tutorial commands use the same basenames under results/tutorial/. Script 05's default prefix now matches its default collapsed pilot. Prefer explicit pairs for provenance clarity. This naming/default change does not modify likelihood, sampling, tuning or fitted results, and a completed explicitly named run need not be repeated. Module 7 already uses gene/mode-specific validation prefixes and is unchanged. Original generated reports retain the source hashes/content of their generation; documentation updates do not retroactively rewrite them.

