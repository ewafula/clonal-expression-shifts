# Integrated clonal-expression report

Generated: 2026-10-18 21:25:48 UTC

Completed gene/mode pairs: 10 of 10 . Incomplete pairs are shown explicitly; this report does not refit models.

## Applicability to other datasets

Module 8 is designed for any dataset processed through compatible versions of this workflow. Gene names, clone IDs, expression patterns, selected branches, stability summaries and IC results are read or computed from the saved run; the current five-gene example is not hard-coded. It requires compatible Modules 2–3 prepared inputs and Module 7 outputs, with referenced files available and matching. It does not directly interpret arbitrary outputs from unrelated software. Use --panel_dir for another completed panel and --output_dir for a separate report. This implementation has been checked on the reviewed five-gene dataset, an independent artificial tutorial dataset and incomplete/mismatched-result tests. Broader validation across independent biological datasets remains outstanding.

## How to interpret this report

Collapsed analysis is primary; uncollapsed cellular analysis is exploratory. Branch selection and fixed-tree cell-subsampling stability do not establish ancestral events, independent biological replication, or causality. No confirmatory p-values, confidence intervals or panel-wide discoveries are supplied.

Response scale: Natural-log log1p(count/library total*10000); synthetic tutorial . Coefficients and contrasts are not automatically DEG log2 fold changes.

Fitting uses inverse-clone-size balancing weights (each clone contributes equal total weight), not known precision weights. EN coefficients are simultaneous penalized contributions; linear coefficients depend on the saved basis; contrasts compare specified groups. They are different estimands. Tuning CV interpolates cells within clones, not held-out patients or clones.

Resampling samples a saved fraction of cells within each clone without replacement and retunes EN; the exact fraction and repetitions appear in each section. Read per-term successful/failed/non-estimable denominators. Quantiles include zero estimates and are descriptive subsampling ranges. Overlap-neighborhood frequencies, where supplied, use Jaccard ≥0.8 including the branch itself; they are not probabilities of a historical region.

Working IC uses ordinary unweighted iid homoskedastic Gaussian cell-error ML refits on the same response/cells, with k = design rank + one residual-variance parameter. This differs from the clone-balanced fit objective. Cell dependence, confounding and post-selection are not corrected. Compare IC only within a gene and evaluated candidate set; never rank genes by raw AIC/BIC. Uncollapsed searches are bounded. No forced agreement between AIC/BIC and EN is expected.

The supplied root and node-label provenance remain unverified. No sample/site adjustment, pseudobulking, or pathway enrichment is added here. Generic small-panel selection is descriptive among prepared candidates, not an independent null-control design.

Prepared expression: 5 genes by 128 cells. Retained clones: 4 . Clone monophyly checks passed for 4 of 4 clones. Original filters/settings remain in the saved Modules 2–3 inputs.

## Panel overview

| gene | selection_reason | fraction_nonzero | clone_mean_range |
| --- | --- | --- | --- |
| DEMO-CLONE-A | user-specified prepared gene | 0.8984375 | 2.9893599 |
| DEMO-CLONE-C-SPARSE | user-specified prepared gene | 0.2109375 | 3.9596430 |
| DEMO-WITHIN-B | user-specified prepared gene | 0.6406250 | 2.3239731 |
| DEMO-BROAD | user-specified prepared gene | 1.0000000 | 0.1286658 |
| DEMO-WEAK | user-specified prepared gene | 0.1484375 | 0.4085494 |

![Training fit overview](module08.overview.png)

Bars compare apparent training fit within each gene. Higher training R-squared in a more flexible model is not evidence of better independent prediction. Genes can have different response scales and variability.

| gene | mode | status | EN_branches | EN_min_frequency | linear_gain | clone_delta_BIC |
| --- | --- | --- | --- | --- | --- | --- |
| DEMO-CLONE-A | collapsed | complete | 2 | 0.30 | 0.000 | 7.661 |
| DEMO-CLONE-A | uncollapsed | complete | 2 | 0.41 | 0.078 | 16.670 |
| DEMO-CLONE-C-SPARSE | collapsed | complete | 2 | 0.60 | 0.000 | 9.350 |
| DEMO-CLONE-C-SPARSE | uncollapsed | complete | 3 | 0.26 | 0.011 | 9.350 |
| DEMO-WITHIN-B | collapsed | complete | 1 | 0.64 | 0.000 | 9.254 |
| DEMO-WITHIN-B | uncollapsed | complete | 1 | 1.00 | 0.217 | 37.165 |
| DEMO-BROAD | collapsed | complete | 0 | NA | 0.000 | 11.416 |
| DEMO-BROAD | uncollapsed | complete | 0 | NA | 0.141 | 14.080 |
| DEMO-WEAK | collapsed | complete | 0 | NA | 0.000 | 13.399 |
| DEMO-WEAK | uncollapsed | complete | 0 | NA | 0.053 | 13.399 |

A positive clone delta BIC means an evaluated model has lower working BIC than full clone means; zero includes equivalent clone-mean parameterizations. It does not establish biological ancestry. The linear gain is the apparent R-squared difference from clone means, not an independent test.

## Comparison between tree modes

DEMO-CLONE-A : uncollapsed minus collapsed apparent R-squared is -0.005 for EN and 0.078 for linear fitting. Full-fit selected EN branches: 2 collapsed versus 2 uncollapsed. These counts use different predictor dictionaries and do not count unique historical events.

DEMO-CLONE-C-SPARSE : uncollapsed minus collapsed apparent R-squared is -0.022 for EN and 0.011 for linear fitting. Full-fit selected EN branches: 2 collapsed versus 3 uncollapsed. These counts use different predictor dictionaries and do not count unique historical events.

DEMO-WITHIN-B : uncollapsed minus collapsed apparent R-squared is 0.177 for EN and 0.217 for linear fitting. Full-fit selected EN branches: 1 collapsed versus 1 uncollapsed. These counts use different predictor dictionaries and do not count unique historical events.

DEMO-BROAD : uncollapsed minus collapsed apparent R-squared is 0.000 for EN and 0.141 for linear fitting. Full-fit selected EN branches: 0 collapsed versus 0 uncollapsed. These counts use different predictor dictionaries and do not count unique historical events.

DEMO-WEAK : uncollapsed minus collapsed apparent R-squared is 0.000 for EN and 0.053 for linear fitting. Full-fit selected EN branches: 0 collapsed versus 0 uncollapsed. These counts use different predictor dictionaries and do not count unique historical events.

Both modes use the same prepared response/cells and clone-balanced weights. Compare observed clone means and explicit contrast definitions; penalized coefficients can redistribute across nested branches and cannot be matched by size alone. IC searches differ between modes, so differing winning models need not indicate a contradiction.

## DEMO-CLONE-A

Selection reason: user-specified prepared gene .

Observed clone means range from 4.271 (Clone C ) to 7.260 (Clone A ). Overall detection is 89.8 percent. These are descriptive expression summaries, not ancestral changes.

| clone | cells | mean | detection_percent |
| --- | --- | --- | --- |
| A | 32 | 7.260 | 100.0 |
| B | 32 | 4.544 | 84.4 |
| C | 32 | 4.271 | 81.2 |
| D | 32 | 4.848 | 93.8 |

### collapsed

Validation settings: 100 subsamples; 80.0 percent per clone without replacement; EN alpha 0.5 and lambda.1se .

Training clone-balanced R-squared: clone means 0.353 ; Elastic Net 0.241 ; linear 0.353 . Linear minus clone means: 0.000 . Greater flexibility can improve training fit without improving independent prediction.

Elastic Net selected 2 branches in the full fit. Successful EN resamples: 100 of 100 . Linear resamples: 100 successful and 0 non-estimable. Contrast-specific usable resamples range from 100 to 100 . Read term-specific denominators; ranges below are not confidence intervals.

Largest full-fit selected EN contributions (up to ten; all branches are in the combined table):

| branch | coefficient | selected_fraction | successful | median | q05 | q95 |
| --- | --- | --- | --- | --- | --- | --- |
| branch_edge_251 | 1.213 | 1.00 | 100 | 0.944 | 0.449 | 1.380 |
| branch_edge_252 | -0.021 | 0.30 | 100 | 0.000 | -0.189 | 0.000 |

Lowest evaluated working BIC: branch_subset_3 = 501.990 ; full clone-mean model delta BIC = 7.661 . Lowest evaluated AIC model: branch_subset_3 . Equivalent prediction spaces can have different branch names. This comparison is within this gene/mode; candidate searches may differ between modes.

Largest absolute descriptive contrasts (not an independent significance ranking):

| contrast | difference | usable | attempted |
| --- | --- | --- | --- |
| {A} vs {B,C,D} | 2.705 | 100 | 100 |
| {B,C} vs {A,D} | -1.647 | 100 | 100 |
| {C} vs {A,B,D} | -1.280 | 100 | 100 |

[Coefficients](../module07/DEMO-CLONE-A/DEMO-CLONE-A_collapsed_attempt1.coefficients.tsv) · [All contrast requests](../module07/DEMO-CLONE-A/DEMO-CLONE-A_collapsed_attempt1.contrasts.tsv) · [Stability](../module07/DEMO-CLONE-A/DEMO-CLONE-A_collapsed_validation_attempt1.stability.tsv) · [Iteration status](../module07/DEMO-CLONE-A/DEMO-CLONE-A_collapsed_validation_attempt1.iteration_status.tsv) · [Working IC](../module07/DEMO-CLONE-A/DEMO-CLONE-A_collapsed_validation_attempt1.model_comparison.tsv)

![DEMO-CLONE-A collapsed expression and training diagnostics](../module07/DEMO-CLONE-A/DEMO-CLONE-A_collapsed_attempt1.comparison.png)

The comparison figure shows observed expression, training fits and descriptive contrasts; it is not independent validation.

![DEMO-CLONE-A collapsed resampling stability](../module07/DEMO-CLONE-A/DEMO-CLONE-A_collapsed_validation_attempt1.stability.png)

The stability figure shows selection frequencies and subsampling ranges. Ranges are not confidence intervals. Full tables retain terms omitted from plots.

[Collapsed selected clone tree](../module07/DEMO-CLONE-A/DEMO-CLONE-A_collapsed_attempt1.selected_tree.png)

Black branches are unselected; red/blue mark positive/negative simultaneous EN contributions. Labels identify original branch terms and coefficients. This topology view does not identify unique historical shifts.

### uncollapsed

Validation settings: 100 subsamples; 80.0 percent per clone without replacement; EN alpha 0.5 and lambda.1se .

Training clone-balanced R-squared: clone means 0.353 ; Elastic Net 0.236 ; linear 0.431 . Linear minus clone means: 0.078 . Greater flexibility can improve training fit without improving independent prediction.

Elastic Net selected 2 branches in the full fit. Successful EN resamples: 100 of 100 . Linear resamples: 100 successful and 0 non-estimable. Contrast-specific usable resamples range from 99 to 100 . Read term-specific denominators; ranges below are not confidence intervals.

Largest full-fit selected EN contributions (up to ten; all branches are in the combined table):

| branch | coefficient | selected_fraction | successful | median | q05 | q95 |
| --- | --- | --- | --- | --- | --- | --- |
| branch_edge_251 | 1.098 | 0.99 | 100 | 0.891 | 0.306 | 1.245 |
| branch_edge_92 | -0.200 | 0.41 | 100 | 0.000 | -1.052 | 0.000 |

Lowest evaluated working BIC: pilot_active_set = 492.982 ; full clone-mean model delta BIC = 16.670 . Lowest evaluated AIC model: pilot_active_set . Equivalent prediction spaces can have different branch names. This comparison is within this gene/mode; candidate searches may differ between modes.

Largest absolute descriptive contrasts (not an independent significance ranking):

| contrast | difference | usable | attempted |
| --- | --- | --- | --- |
| Clone C: branch_edge_123 / within_clone | 1.202 | 100 | 100 |
| Clone C: branch_edge_123 / sister_within_clone | 1.202 | 100 | 100 |
| Clone C: branch_edge_124 / within_clone | -1.202 | 100 | 100 |

All contrast requests: insufficient_cells=32; no_within_clone_split=184; ready=16 . Combined contrast table includes ready contrasts only; unavailable requests remain in the linked original table. No within-clone ancestry is established.

[Coefficients](../module07/DEMO-CLONE-A/DEMO-CLONE-A_uncollapsed_attempt1.coefficients.tsv) · [All contrast requests](../module07/DEMO-CLONE-A/DEMO-CLONE-A_uncollapsed_attempt1.contrasts.tsv) · [Stability](../module07/DEMO-CLONE-A/DEMO-CLONE-A_uncollapsed_validation_attempt1.stability.tsv) · [Iteration status](../module07/DEMO-CLONE-A/DEMO-CLONE-A_uncollapsed_validation_attempt1.iteration_status.tsv) · [Working IC](../module07/DEMO-CLONE-A/DEMO-CLONE-A_uncollapsed_validation_attempt1.model_comparison.tsv)

![DEMO-CLONE-A uncollapsed expression and training diagnostics](../module07/DEMO-CLONE-A/DEMO-CLONE-A_uncollapsed_attempt1.comparison.png)

The comparison figure shows observed expression, training fits and descriptive contrasts; it is not independent validation.

![DEMO-CLONE-A uncollapsed resampling stability](../module07/DEMO-CLONE-A/DEMO-CLONE-A_uncollapsed_validation_attempt1.stability.png)

The stability figure shows selection frequencies and subsampling ranges. Ranges are not confidence intervals. Full tables retain terms omitted from plots.

[Exploratory selected cellular tree](../module07/DEMO-CLONE-A/DEMO-CLONE-A_uncollapsed_attempt1.selected_tree.png)

Black branches show the supplied cellular topology; red/blue mark positive/negative EN contributions. The tip strip identifies clones. Selection does not validate ancestry.

## DEMO-CLONE-C-SPARSE

Selection reason: user-specified prepared gene .

Observed clone means range from 0.000 (Clone A ) to 3.960 (Clone C ). Overall detection is 21.1 percent. These are descriptive expression summaries, not ancestral changes.

| clone | cells | mean | detection_percent |
| --- | --- | --- | --- |
| A | 32 | 0.000 | 0.0 |
| B | 32 | 0.135 | 3.1 |
| C | 32 | 3.960 | 81.2 |
| D | 32 | 0.000 | 0.0 |

### collapsed

Validation settings: 100 subsamples; 80.0 percent per clone without replacement; EN alpha 0.5 and lambda.1se .

Training clone-balanced R-squared: clone means 0.724 ; Elastic Net 0.691 ; linear 0.724 . Linear minus clone means: 0.000 . Greater flexibility can improve training fit without improving independent prediction.

Elastic Net selected 2 branches in the full fit. Successful EN resamples: 100 of 100 . Linear resamples: 100 successful and 0 non-estimable. Contrast-specific usable resamples range from 100 to 100 . Read term-specific denominators; ranges below are not confidence intervals.

Largest full-fit selected EN contributions (up to ten; all branches are in the combined table):

| branch | coefficient | selected_fraction | successful | median | q05 | q95 |
| --- | --- | --- | --- | --- | --- | --- |
| branch_edge_188 | 3.017 | 1.00 | 100 | 2.657 | 2.089 | 3.145 |
| branch_edge_252 | 0.091 | 0.60 | 100 | 0.017 | 0.000 | 0.109 |

Lowest evaluated working BIC: branch_subset_2 = 390.150 ; full clone-mean model delta BIC = 9.350 . Lowest evaluated AIC model: branch_subset_2 . Equivalent prediction spaces can have different branch names. This comparison is within this gene/mode; candidate searches may differ between modes.

Largest absolute descriptive contrasts (not an independent significance ranking):

| contrast | difference | usable | attempted |
| --- | --- | --- | --- |
| {C} vs {A,B,D} | 3.915 | 100 | 100 |
| {B,C} vs {A,D} | 2.047 | 100 | 100 |
| {A} vs {B,C,D} | -1.365 | 100 | 100 |

[Coefficients](../module07/DEMO-CLONE-C-SPARSE/DEMO-CLONE-C-SPARSE_collapsed_attempt1.coefficients.tsv) · [All contrast requests](../module07/DEMO-CLONE-C-SPARSE/DEMO-CLONE-C-SPARSE_collapsed_attempt1.contrasts.tsv) · [Stability](../module07/DEMO-CLONE-C-SPARSE/DEMO-CLONE-C-SPARSE_collapsed_validation_attempt1.stability.tsv) · [Iteration status](../module07/DEMO-CLONE-C-SPARSE/DEMO-CLONE-C-SPARSE_collapsed_validation_attempt1.iteration_status.tsv) · [Working IC](../module07/DEMO-CLONE-C-SPARSE/DEMO-CLONE-C-SPARSE_collapsed_validation_attempt1.model_comparison.tsv)

![DEMO-CLONE-C-SPARSE collapsed expression and training diagnostics](../module07/DEMO-CLONE-C-SPARSE/DEMO-CLONE-C-SPARSE_collapsed_attempt1.comparison.png)

The comparison figure shows observed expression, training fits and descriptive contrasts; it is not independent validation.

![DEMO-CLONE-C-SPARSE collapsed resampling stability](../module07/DEMO-CLONE-C-SPARSE/DEMO-CLONE-C-SPARSE_collapsed_validation_attempt1.stability.png)

The stability figure shows selection frequencies and subsampling ranges. Ranges are not confidence intervals. Full tables retain terms omitted from plots.

[Collapsed selected clone tree](../module07/DEMO-CLONE-C-SPARSE/DEMO-CLONE-C-SPARSE_collapsed_attempt1.selected_tree.png)

Black branches are unselected; red/blue mark positive/negative simultaneous EN contributions. Labels identify original branch terms and coefficients. This topology view does not identify unique historical shifts.

### uncollapsed

Validation settings: 100 subsamples; 80.0 percent per clone without replacement; EN alpha 0.5 and lambda.1se .

Training clone-balanced R-squared: clone means 0.724 ; Elastic Net 0.668 ; linear 0.734 . Linear minus clone means: 0.011 . Greater flexibility can improve training fit without improving independent prediction.

Elastic Net selected 3 branches in the full fit. Successful EN resamples: 100 of 100 . Linear resamples: 100 successful and 0 non-estimable. Contrast-specific usable resamples range from 99 to 100 . Read term-specific denominators; ranges below are not confidence intervals.

Largest full-fit selected EN contributions (up to ten; all branches are in the combined table):

| branch | coefficient | selected_fraction | successful | median | q05 | q95 |
| --- | --- | --- | --- | --- | --- | --- |
| branch_edge_188 | 2.516 | 1.00 | 100 | 2.190 | 1.822 | 2.614 |
| branch_edge_124 | 0.494 | 0.98 | 100 | 0.450 | 0.056 | 0.949 |
| branch_edge_252 | 0.027 | 0.26 | 100 | 0.000 | 0.000 | 0.045 |

Lowest evaluated working BIC: single_20 = 390.150 ; full clone-mean model delta BIC = 9.350 . Lowest evaluated AIC model: ancestor_child_branch_edge_124 . Equivalent prediction spaces can have different branch names. This comparison is within this gene/mode; candidate searches may differ between modes.

Largest absolute descriptive contrasts (not an independent significance ranking):

| contrast | difference | usable | attempted |
| --- | --- | --- | --- |
| Clone C: branch_edge_123 / within_clone | -0.611 | 100 | 100 |
| Clone C: branch_edge_123 / sister_within_clone | -0.611 | 100 | 100 |
| Clone C: branch_edge_124 / within_clone | 0.611 | 100 | 100 |

All contrast requests: insufficient_cells=32; no_within_clone_split=184; ready=16 . Combined contrast table includes ready contrasts only; unavailable requests remain in the linked original table. No within-clone ancestry is established.

[Coefficients](../module07/DEMO-CLONE-C-SPARSE/DEMO-CLONE-C-SPARSE_uncollapsed_attempt1.coefficients.tsv) · [All contrast requests](../module07/DEMO-CLONE-C-SPARSE/DEMO-CLONE-C-SPARSE_uncollapsed_attempt1.contrasts.tsv) · [Stability](../module07/DEMO-CLONE-C-SPARSE/DEMO-CLONE-C-SPARSE_uncollapsed_validation_attempt1.stability.tsv) · [Iteration status](../module07/DEMO-CLONE-C-SPARSE/DEMO-CLONE-C-SPARSE_uncollapsed_validation_attempt1.iteration_status.tsv) · [Working IC](../module07/DEMO-CLONE-C-SPARSE/DEMO-CLONE-C-SPARSE_uncollapsed_validation_attempt1.model_comparison.tsv)

![DEMO-CLONE-C-SPARSE uncollapsed expression and training diagnostics](../module07/DEMO-CLONE-C-SPARSE/DEMO-CLONE-C-SPARSE_uncollapsed_attempt1.comparison.png)

The comparison figure shows observed expression, training fits and descriptive contrasts; it is not independent validation.

![DEMO-CLONE-C-SPARSE uncollapsed resampling stability](../module07/DEMO-CLONE-C-SPARSE/DEMO-CLONE-C-SPARSE_uncollapsed_validation_attempt1.stability.png)

The stability figure shows selection frequencies and subsampling ranges. Ranges are not confidence intervals. Full tables retain terms omitted from plots.

[Exploratory selected cellular tree](../module07/DEMO-CLONE-C-SPARSE/DEMO-CLONE-C-SPARSE_uncollapsed_attempt1.selected_tree.png)

Black branches show the supplied cellular topology; red/blue mark positive/negative EN contributions. The tip strip identifies clones. Selection does not validate ancestry.

## DEMO-WITHIN-B

Selection reason: user-specified prepared gene .

Observed clone means range from 2.561 (Clone C ) to 4.885 (Clone B ). Overall detection is 64.1 percent. These are descriptive expression summaries, not ancestral changes.

| clone | cells | mean | detection_percent |
| --- | --- | --- | --- |
| A | 32 | 2.882 | 62.5 |
| B | 32 | 4.885 | 78.1 |
| C | 32 | 2.561 | 53.1 |
| D | 32 | 2.943 | 62.5 |

### collapsed

Validation settings: 100 subsamples; 80.0 percent per clone without replacement; EN alpha 0.5 and lambda.1se .

Training clone-balanced R-squared: clone means 0.123 ; Elastic Net 0.042 ; linear 0.123 . Linear minus clone means: 0.000 . Greater flexibility can improve training fit without improving independent prediction.

Elastic Net selected 1 branches in the full fit. Successful EN resamples: 100 of 100 . Linear resamples: 100 successful and 0 non-estimable. Contrast-specific usable resamples range from 100 to 100 . Read term-specific denominators; ranges below are not confidence intervals.

Largest full-fit selected EN contributions (up to ten; all branches are in the combined table):

| branch | coefficient | selected_fraction | successful | median | q05 | q95 |
| --- | --- | --- | --- | --- | --- | --- |
| branch_edge_187 | 0.403 | 0.64 | 100 | 0.407 | 0.000 | 0.952 |

Lowest evaluated working BIC: branch_subset_1 = 607.013 ; full clone-mean model delta BIC = 9.254 . Lowest evaluated AIC model: branch_subset_1 . Equivalent prediction spaces can have different branch names. This comparison is within this gene/mode; candidate searches may differ between modes.

Largest absolute descriptive contrasts (not an independent significance ranking):

| contrast | difference | usable | attempted |
| --- | --- | --- | --- |
| {B} vs {A,C,D} | 2.090 | 100 | 100 |
| {C} vs {A,B,D} | -1.009 | 100 | 100 |
| {B,C} vs {A,D} | 0.811 | 100 | 100 |

[Coefficients](../module07/DEMO-WITHIN-B/DEMO-WITHIN-B_collapsed_attempt1.coefficients.tsv) · [All contrast requests](../module07/DEMO-WITHIN-B/DEMO-WITHIN-B_collapsed_attempt1.contrasts.tsv) · [Stability](../module07/DEMO-WITHIN-B/DEMO-WITHIN-B_collapsed_validation_attempt1.stability.tsv) · [Iteration status](../module07/DEMO-WITHIN-B/DEMO-WITHIN-B_collapsed_validation_attempt1.iteration_status.tsv) · [Working IC](../module07/DEMO-WITHIN-B/DEMO-WITHIN-B_collapsed_validation_attempt1.model_comparison.tsv)

![DEMO-WITHIN-B collapsed expression and training diagnostics](../module07/DEMO-WITHIN-B/DEMO-WITHIN-B_collapsed_attempt1.comparison.png)

The comparison figure shows observed expression, training fits and descriptive contrasts; it is not independent validation.

![DEMO-WITHIN-B collapsed resampling stability](../module07/DEMO-WITHIN-B/DEMO-WITHIN-B_collapsed_validation_attempt1.stability.png)

The stability figure shows selection frequencies and subsampling ranges. Ranges are not confidence intervals. Full tables retain terms omitted from plots.

[Collapsed selected clone tree](../module07/DEMO-WITHIN-B/DEMO-WITHIN-B_collapsed_attempt1.selected_tree.png)

Black branches are unselected; red/blue mark positive/negative simultaneous EN contributions. Labels identify original branch terms and coefficients. This topology view does not identify unique historical shifts.

### uncollapsed

Validation settings: 100 subsamples; 80.0 percent per clone without replacement; EN alpha 0.5 and lambda.1se .

Training clone-balanced R-squared: clone means 0.123 ; Elastic Net 0.219 ; linear 0.340 . Linear minus clone means: 0.217 . Greater flexibility can improve training fit without improving independent prediction.

Elastic Net selected 1 branches in the full fit. Successful EN resamples: 100 of 100 . Linear resamples: 100 successful and 0 non-estimable. Contrast-specific usable resamples range from 99 to 100 . Read term-specific denominators; ranges below are not confidence intervals.

Largest full-fit selected EN contributions (up to ten; all branches are in the combined table):

| branch | coefficient | selected_fraction | successful | median | q05 | q95 |
| --- | --- | --- | --- | --- | --- | --- |
| branch_edge_185 | 2.131 | 1.00 | 100 | 1.944 | 1.354 | 2.579 |

Lowest evaluated working BIC: single_17 = 579.103 ; full clone-mean model delta BIC = 37.165 . Lowest evaluated AIC model: single_17 . Equivalent prediction spaces can have different branch names. This comparison is within this gene/mode; candidate searches may differ between modes.

Largest absolute descriptive contrasts (not an independent significance ranking):

| contrast | difference | usable | attempted |
| --- | --- | --- | --- |
| Clone B: branch_edge_185 / within_clone | 4.336 | 100 | 100 |
| Clone B: branch_edge_185 / sister_within_clone | 4.336 | 100 | 100 |
| Clone B: branch_edge_186 / within_clone | -4.336 | 100 | 100 |

All contrast requests: insufficient_cells=32; no_within_clone_split=184; ready=16 . Combined contrast table includes ready contrasts only; unavailable requests remain in the linked original table. No within-clone ancestry is established.

[Coefficients](../module07/DEMO-WITHIN-B/DEMO-WITHIN-B_uncollapsed_attempt1.coefficients.tsv) · [All contrast requests](../module07/DEMO-WITHIN-B/DEMO-WITHIN-B_uncollapsed_attempt1.contrasts.tsv) · [Stability](../module07/DEMO-WITHIN-B/DEMO-WITHIN-B_uncollapsed_validation_attempt1.stability.tsv) · [Iteration status](../module07/DEMO-WITHIN-B/DEMO-WITHIN-B_uncollapsed_validation_attempt1.iteration_status.tsv) · [Working IC](../module07/DEMO-WITHIN-B/DEMO-WITHIN-B_uncollapsed_validation_attempt1.model_comparison.tsv)

![DEMO-WITHIN-B uncollapsed expression and training diagnostics](../module07/DEMO-WITHIN-B/DEMO-WITHIN-B_uncollapsed_attempt1.comparison.png)

The comparison figure shows observed expression, training fits and descriptive contrasts; it is not independent validation.

![DEMO-WITHIN-B uncollapsed resampling stability](../module07/DEMO-WITHIN-B/DEMO-WITHIN-B_uncollapsed_validation_attempt1.stability.png)

The stability figure shows selection frequencies and subsampling ranges. Ranges are not confidence intervals. Full tables retain terms omitted from plots.

[Exploratory selected cellular tree](../module07/DEMO-WITHIN-B/DEMO-WITHIN-B_uncollapsed_attempt1.selected_tree.png)

Black branches show the supplied cellular topology; red/blue mark positive/negative EN contributions. The tip strip identifies clones. Selection does not validate ancestry.

## DEMO-BROAD

Selection reason: user-specified prepared gene .

Observed clone means range from 6.496 (Clone A ) to 6.625 (Clone C ). Overall detection is 100.0 percent. These are descriptive expression summaries, not ancestral changes.

| clone | cells | mean | detection_percent |
| --- | --- | --- | --- |
| A | 32 | 6.496 | 100.0 |
| B | 32 | 6.578 | 100.0 |
| C | 32 | 6.625 | 100.0 |
| D | 32 | 6.566 | 100.0 |

### collapsed

Validation settings: 100 subsamples; 80.0 percent per clone without replacement; EN alpha 0.5 and lambda.1se .

Training clone-balanced R-squared: clone means 0.024 ; Elastic Net 0.000 ; linear 0.024 . Linear minus clone means: 0.000 . Greater flexibility can improve training fit without improving independent prediction.

Elastic Net selected 0 branches in the full fit. Successful EN resamples: 100 of 100 . Linear resamples: 100 successful and 0 non-estimable. Contrast-specific usable resamples range from 100 to 100 . Read term-specific denominators; ranges below are not confidence intervals.

No full-fit EN branch selected (or EN was not fitted); this does not establish absence of expression differences.

Lowest evaluated working BIC: intercept_only = 61.163 ; full clone-mean model delta BIC = 11.416 . Lowest evaluated AIC model: branch_subset_3 . Equivalent prediction spaces can have different branch names. This comparison is within this gene/mode; candidate searches may differ between modes.

Largest absolute descriptive contrasts (not an independent significance ranking):

| contrast | difference | usable | attempted |
| --- | --- | --- | --- |
| {A} vs {B,C,D} | -0.094 | 100 | 100 |
| {C} vs {A,B,D} | 0.078 | 100 | 100 |
| {B,C} vs {A,D} | 0.070 | 100 | 100 |

[Coefficients](../module07/DEMO-BROAD/DEMO-BROAD_collapsed_attempt1.coefficients.tsv) · [All contrast requests](../module07/DEMO-BROAD/DEMO-BROAD_collapsed_attempt1.contrasts.tsv) · [Stability](../module07/DEMO-BROAD/DEMO-BROAD_collapsed_validation_attempt1.stability.tsv) · [Iteration status](../module07/DEMO-BROAD/DEMO-BROAD_collapsed_validation_attempt1.iteration_status.tsv) · [Working IC](../module07/DEMO-BROAD/DEMO-BROAD_collapsed_validation_attempt1.model_comparison.tsv)

![DEMO-BROAD collapsed expression and training diagnostics](../module07/DEMO-BROAD/DEMO-BROAD_collapsed_attempt1.comparison.png)

The comparison figure shows observed expression, training fits and descriptive contrasts; it is not independent validation.

![DEMO-BROAD collapsed resampling stability](../module07/DEMO-BROAD/DEMO-BROAD_collapsed_validation_attempt1.stability.png)

The stability figure shows selection frequencies and subsampling ranges. Ranges are not confidence intervals. Full tables retain terms omitted from plots.

[Collapsed selected clone tree](../module07/DEMO-BROAD/DEMO-BROAD_collapsed_attempt1.selected_tree.png)

Black branches are unselected; red/blue mark positive/negative simultaneous EN contributions. Labels identify original branch terms and coefficients. This topology view does not identify unique historical shifts.

### uncollapsed

Validation settings: 100 subsamples; 80.0 percent per clone without replacement; EN alpha 0.5 and lambda.1se .

Training clone-balanced R-squared: clone means 0.024 ; Elastic Net 0.000 ; linear 0.165 . Linear minus clone means: 0.141 . Greater flexibility can improve training fit without improving independent prediction.

Elastic Net selected 0 branches in the full fit. Successful EN resamples: 100 of 100 . Linear resamples: 100 successful and 0 non-estimable. Contrast-specific usable resamples range from 99 to 100 . Read term-specific denominators; ranges below are not confidence intervals.

No full-fit EN branch selected (or EN was not fitted); this does not establish absence of expression differences.

Lowest evaluated working BIC: single_18 = 58.499 ; full clone-mean model delta BIC = 14.080 . Lowest evaluated AIC model: single_18 . Equivalent prediction spaces can have different branch names. This comparison is within this gene/mode; candidate searches may differ between modes.

Largest absolute descriptive contrasts (not an independent significance ranking):

| contrast | difference | usable | attempted |
| --- | --- | --- | --- |
| Clone B: branch_edge_185 / within_clone | -0.350 | 100 | 100 |
| Clone B: branch_edge_185 / sister_within_clone | -0.350 | 100 | 100 |
| Clone B: branch_edge_186 / within_clone | 0.350 | 100 | 100 |

All contrast requests: insufficient_cells=32; no_within_clone_split=184; ready=16 . Combined contrast table includes ready contrasts only; unavailable requests remain in the linked original table. No within-clone ancestry is established.

[Coefficients](../module07/DEMO-BROAD/DEMO-BROAD_uncollapsed_attempt1.coefficients.tsv) · [All contrast requests](../module07/DEMO-BROAD/DEMO-BROAD_uncollapsed_attempt1.contrasts.tsv) · [Stability](../module07/DEMO-BROAD/DEMO-BROAD_uncollapsed_validation_attempt1.stability.tsv) · [Iteration status](../module07/DEMO-BROAD/DEMO-BROAD_uncollapsed_validation_attempt1.iteration_status.tsv) · [Working IC](../module07/DEMO-BROAD/DEMO-BROAD_uncollapsed_validation_attempt1.model_comparison.tsv)

![DEMO-BROAD uncollapsed expression and training diagnostics](../module07/DEMO-BROAD/DEMO-BROAD_uncollapsed_attempt1.comparison.png)

The comparison figure shows observed expression, training fits and descriptive contrasts; it is not independent validation.

![DEMO-BROAD uncollapsed resampling stability](../module07/DEMO-BROAD/DEMO-BROAD_uncollapsed_validation_attempt1.stability.png)

The stability figure shows selection frequencies and subsampling ranges. Ranges are not confidence intervals. Full tables retain terms omitted from plots.

[Exploratory selected cellular tree](../module07/DEMO-BROAD/DEMO-BROAD_uncollapsed_attempt1.selected_tree.png)

Black branches show the supplied cellular topology; red/blue mark positive/negative EN contributions. The tip strip identifies clones. Selection does not validate ancestry.

## DEMO-WEAK

Selection reason: user-specified prepared gene .

Observed clone means range from 0.415 (Clone C ) to 0.824 (Clone B ). Overall detection is 14.8 percent. These are descriptive expression summaries, not ancestral changes.

| clone | cells | mean | detection_percent |
| --- | --- | --- | --- |
| A | 32 | 0.703 | 15.6 |
| B | 32 | 0.824 | 18.8 |
| C | 32 | 0.415 | 9.4 |
| D | 32 | 0.700 | 15.6 |

### collapsed

Validation settings: 100 subsamples; 80.0 percent per clone without replacement; EN alpha 0.5 and lambda.1se .

Training clone-balanced R-squared: clone means 0.009 ; Elastic Net 0.000 ; linear 0.009 . Linear minus clone means: 0.000 . Greater flexibility can improve training fit without improving independent prediction.

Elastic Net selected 0 branches in the full fit. Successful EN resamples: 100 of 100 . Linear resamples: 100 successful and 0 non-estimable. Contrast-specific usable resamples range from 100 to 100 . Read term-specific denominators; ranges below are not confidence intervals.

No full-fit EN branch selected (or EN was not fitted); this does not establish absence of expression differences.

Lowest evaluated working BIC: intercept_only = 490.530 ; full clone-mean model delta BIC = 13.399 . Lowest evaluated AIC model: intercept_only . Equivalent prediction spaces can have different branch names. This comparison is within this gene/mode; candidate searches may differ between modes.

Largest absolute descriptive contrasts (not an independent significance ranking):

| contrast | difference | usable | attempted |
| --- | --- | --- | --- |
| {C} vs {A,B,D} | -0.327 | 100 | 100 |
| {B} vs {A,C,D} | 0.218 | 100 | 100 |
| {B,C} vs {A,D} | -0.082 | 100 | 100 |

[Coefficients](../module07/DEMO-WEAK/DEMO-WEAK_collapsed_attempt1.coefficients.tsv) · [All contrast requests](../module07/DEMO-WEAK/DEMO-WEAK_collapsed_attempt1.contrasts.tsv) · [Stability](../module07/DEMO-WEAK/DEMO-WEAK_collapsed_validation_attempt1.stability.tsv) · [Iteration status](../module07/DEMO-WEAK/DEMO-WEAK_collapsed_validation_attempt1.iteration_status.tsv) · [Working IC](../module07/DEMO-WEAK/DEMO-WEAK_collapsed_validation_attempt1.model_comparison.tsv)

![DEMO-WEAK collapsed expression and training diagnostics](../module07/DEMO-WEAK/DEMO-WEAK_collapsed_attempt1.comparison.png)

The comparison figure shows observed expression, training fits and descriptive contrasts; it is not independent validation.

![DEMO-WEAK collapsed resampling stability](../module07/DEMO-WEAK/DEMO-WEAK_collapsed_validation_attempt1.stability.png)

The stability figure shows selection frequencies and subsampling ranges. Ranges are not confidence intervals. Full tables retain terms omitted from plots.

[Collapsed selected clone tree](../module07/DEMO-WEAK/DEMO-WEAK_collapsed_attempt1.selected_tree.png)

Black branches are unselected; red/blue mark positive/negative simultaneous EN contributions. Labels identify original branch terms and coefficients. This topology view does not identify unique historical shifts.

### uncollapsed

Validation settings: 100 subsamples; 80.0 percent per clone without replacement; EN alpha 0.5 and lambda.1se .

Training clone-balanced R-squared: clone means 0.009 ; Elastic Net 0.000 ; linear 0.062 . Linear minus clone means: 0.053 . Greater flexibility can improve training fit without improving independent prediction.

Elastic Net selected 0 branches in the full fit. Successful EN resamples: 100 of 100 . Linear resamples: 100 successful and 0 non-estimable. Contrast-specific usable resamples range from 99 to 100 . Read term-specific denominators; ranges below are not confidence intervals.

No full-fit EN branch selected (or EN was not fitted); this does not establish absence of expression differences.

Lowest evaluated working BIC: intercept_only = 490.530 ; full clone-mean model delta BIC = 13.399 . Lowest evaluated AIC model: single_15 . Equivalent prediction spaces can have different branch names. This comparison is within this gene/mode; candidate searches may differ between modes.

Largest absolute descriptive contrasts (not an independent significance ranking):

| contrast | difference | usable | attempted |
| --- | --- | --- | --- |
| Clone B: branch_edge_185 / within_clone | 0.498 | 100 | 100 |
| Clone B: branch_edge_185 / sister_within_clone | 0.498 | 100 | 100 |
| Clone B: branch_edge_186 / within_clone | -0.498 | 100 | 100 |

All contrast requests: insufficient_cells=32; no_within_clone_split=184; ready=16 . Combined contrast table includes ready contrasts only; unavailable requests remain in the linked original table. No within-clone ancestry is established.

[Coefficients](../module07/DEMO-WEAK/DEMO-WEAK_uncollapsed_attempt1.coefficients.tsv) · [All contrast requests](../module07/DEMO-WEAK/DEMO-WEAK_uncollapsed_attempt1.contrasts.tsv) · [Stability](../module07/DEMO-WEAK/DEMO-WEAK_uncollapsed_validation_attempt1.stability.tsv) · [Iteration status](../module07/DEMO-WEAK/DEMO-WEAK_uncollapsed_validation_attempt1.iteration_status.tsv) · [Working IC](../module07/DEMO-WEAK/DEMO-WEAK_uncollapsed_validation_attempt1.model_comparison.tsv)

![DEMO-WEAK uncollapsed expression and training diagnostics](../module07/DEMO-WEAK/DEMO-WEAK_uncollapsed_attempt1.comparison.png)

The comparison figure shows observed expression, training fits and descriptive contrasts; it is not independent validation.

![DEMO-WEAK uncollapsed resampling stability](../module07/DEMO-WEAK/DEMO-WEAK_uncollapsed_validation_attempt1.stability.png)

The stability figure shows selection frequencies and subsampling ranges. Ranges are not confidence intervals. Full tables retain terms omitted from plots.

[Exploratory selected cellular tree](../module07/DEMO-WEAK/DEMO-WEAK_uncollapsed_attempt1.selected_tree.png)

Black branches show the supplied cellular topology; red/blue mark positive/negative EN contributions. The tip strip identifies clones. Selection does not validate ancestry.

## Traceability and reading the tables

Panel source: 
[Module 7 panel](../module07/module07.panel.tsv)
[Module 7 stage history](../module07/module07.status.tsv)

[Gene/mode overview](module08.gene_summary.tsv)
[Training metrics](module08.fit_metrics.tsv)
[All branch stability estimates](module08.branch_stability.tsv)
[Ready contrasts](module08.contrasts.tsv)
[Source hashes](module08.sources.tsv)

[Column definitions and worked examples](../../../docs/RESULTS_GUIDE.md)

Keep this report with the input panel folders so relative figure/table links work. Generated private results remain local and excluded from GitHub. A public tutorial must use shareable synthetic inputs.

