#!/usr/bin/env Rscript
# Candidate filtering/ranking
# Additional work: prepare normalized gene expression in Module 2 cell order.

rank_candidates <- function(deg_df, opt) {
deg_df <- deg_df %>%
  mutate(
    gene = as.character(gene),
    clone_id = as.character(clone_id),
    p_val = as.numeric(p_val),
    avg_log2FC = as.numeric(avg_log2FC),
    pct.1 = as.numeric(pct.1),
    pct.2 = as.numeric(pct.2),
    p_val_adj = as.numeric(p_val_adj),
    abs_avg_log2FC = abs(avg_log2FC),
    pct_detection_difference = pct.1 - pct.2,
    abs_pct_detection_difference = abs(pct_detection_difference),
    neg_log10_padj = -log10(pmax(p_val_adj, .Machine$double.xmin))
  ) %>%
  filter(!is.na(gene), gene != "")

message("  [INFO] Unique genes before filtering: ", n_distinct(deg_df$gene))
message("  [INFO] Clone IDs represented: ", paste(sort(unique(deg_df$clone_id)), collapse = ", "))

# ------------------------------------------------------------
# Step 3: Apply statistical candidate filter
# ------------------------------------------------------------

message("\n[Step 3] Applying statistical DEG candidate filter...")

message("  [FILTER] p_val_adj < ", opt$padj_threshold)
message("  [FILTER] abs(avg_log2FC) >= ", opt$logfc_threshold)
message("  [FILTER] pct.1 >= ", opt$pct1_threshold)

filtered_deg <- deg_df %>%
  filter(
    p_val_adj < opt$padj_threshold,
    abs_avg_log2FC >= opt$logfc_threshold,
    pct.1 >= opt$pct1_threshold
  )

message("  [INFO] Rows passing filters: ", nrow(filtered_deg))
message("  [INFO] Unique genes passing filters: ", n_distinct(filtered_deg$gene))

if (nrow(filtered_deg) == 0) {
  stop(
    "No DEG rows passed the filters. Consider relaxing thresholds or checking the input table."
  )
}

# ------------------------------------------------------------
# Step 4: Create gene-level ranking metrics
# ------------------------------------------------------------

message("\n[Step 4] Ranking candidate genes...")

# Ranking strategy:
#   Genes are prioritized by:
#     1. Strongest clone-specific effect size.
#     2. Strongest adjusted significance.
#     3. Higher detection in focal clone.
#     4. Larger focal-vs-background detection difference.
#     5. Number of clones where the gene passes filters.
#
# Notes:
#   - max_abs_log2FC captures strongest clone-associated shift.
#   - best_padj captures strongest evidence.
#   - max_pct1 protects against very rare/outlier-driven genes.
#   - max_abs_pct_detection_difference favors clone-localized expression.
#   - n_passing_clones helps identify genes repeatedly DE across clones.

candidate_gene_summary <- filtered_deg %>%
  group_by(gene) %>%
  summarise(
    n_passing_deg_rows = n(),
    n_passing_clones = n_distinct(clone_id),
    passing_clones = paste(sort(unique(clone_id)), collapse = ";"),

    best_padj = min(p_val_adj, na.rm = TRUE),
    best_neg_log10_padj = max(neg_log10_padj, na.rm = TRUE),

    max_abs_log2FC = max(abs_avg_log2FC, na.rm = TRUE),
    max_positive_log2FC = suppressWarnings(max(avg_log2FC, na.rm = TRUE)),
    min_negative_log2FC = suppressWarnings(min(avg_log2FC, na.rm = TRUE)),

    max_pct1 = max(pct.1, na.rm = TRUE),
    median_pct1 = median(pct.1, na.rm = TRUE),

    min_pct2 = min(pct.2, na.rm = TRUE),
    median_pct2 = median(pct.2, na.rm = TRUE),

    max_abs_pct_detection_difference =
      max(abs_pct_detection_difference, na.rm = TRUE),

    max_pct_detection_difference =
      max(pct_detection_difference, na.rm = TRUE),

    mean_abs_log2FC = mean(abs_avg_log2FC, na.rm = TRUE),

    top_clone_by_abs_log2FC =
      clone_id[which.max(abs_avg_log2FC)],

    top_clone_avg_log2FC =
      avg_log2FC[which.max(abs_avg_log2FC)],

    top_clone_pct1 =
      pct.1[which.max(abs_avg_log2FC)],

    top_clone_pct2 =
      pct.2[which.max(abs_avg_log2FC)],

    top_clone_padj =
      p_val_adj[which.max(abs_avg_log2FC)],

    .groups = "drop"
  ) %>%
  mutate(
    # Composite score for practical prioritization.
    #
    # This is not a formal statistical test; it is a ranking score to help
    # choose strong candidate genes for the next Elastic Net script.
    #
    # It rewards:
    #   - stronger effect size
    #   - stronger adjusted significance
    #   - higher focal-clone detection
    #   - clone-localized detection difference
    #   - repeated DEG signal across clones
    candidate_priority_score =
      max_abs_log2FC *
      log1p(best_neg_log10_padj) *
      pmax(max_pct1, 0.001) *
      pmax(max_abs_pct_detection_difference, 0.001) *
      log1p(n_passing_clones),

    suggested_for_initial_elastic_net =
      rank(-candidate_priority_score, ties.method = "first") <= opt$top_n
  ) %>%
  arrange(
    desc(candidate_priority_score),
    desc(max_abs_log2FC),
    desc(best_neg_log10_padj),
    desc(max_pct1),
    desc(max_abs_pct_detection_difference)
  ) %>%
  mutate(
    candidate_rank = row_number(),
    filter_padj_threshold = opt$padj_threshold,
    filter_abs_log2FC_threshold = opt$logfc_threshold,
    filter_pct1_threshold = opt$pct1_threshold
  ) %>%
  select(
    candidate_rank,
    gene,
    candidate_priority_score,
    suggested_for_initial_elastic_net,
    n_passing_deg_rows,
    n_passing_clones,
    passing_clones,
    top_clone_by_abs_log2FC,
    top_clone_avg_log2FC,
    top_clone_pct1,
    top_clone_pct2,
    top_clone_padj,
    best_padj,
    best_neg_log10_padj,
    max_abs_log2FC,
    max_positive_log2FC,
    min_negative_log2FC,
    mean_abs_log2FC,
    max_pct1,
    median_pct1,
    min_pct2,
    median_pct2,
    max_abs_pct_detection_difference,
    max_pct_detection_difference,
    filter_padj_threshold,
    filter_abs_log2FC_threshold,
    filter_pct1_threshold
  )

# ------------------------------------------------------------
# Step 5: Add clone-level evidence string per gene
# ------------------------------------------------------------

message("\n[Step 5] Adding clone-level evidence summary...")

gene_evidence <- filtered_deg %>%
  arrange(gene, desc(abs_avg_log2FC), p_val_adj) %>%
  mutate(
    evidence_label = paste0(
      "clone_", clone_id,
      ":log2FC=", signif(avg_log2FC, 4),
      ",pct1=", signif(pct.1, 3),
      ",pct2=", signif(pct.2, 3),
      ",padj=", signif(p_val_adj, 3)
    )
  ) %>%
  group_by(gene) %>%
  summarise(
    clone_level_evidence = paste(evidence_label, collapse = "; "),
    .groups = "drop"
  )

candidate_gene_summary <- candidate_gene_summary %>%
  left_join(gene_evidence, by = "gene")


  list(candidates=candidate_gene_summary, filtered=filtered_deg)
}

for(pkg in c('dplyr','readr','optparse','SeuratObject','Matrix')) if(!requireNamespace(pkg,quietly=TRUE))
  stop('Missing package ',pkg,'. Run Rscript 00_install_packages.R in this environment.')
suppressPackageStartupMessages(library(dplyr))
options_list <- list(
  optparse::make_option('--seurat_object',default='inputs/SPECTRUM-cell-ploidy-prediction.rds',help='Seurat RDS [default %default]'),
  optparse::make_option('--tree_designs',default='results/module02.tree_designs.rds',help='Module 2 saved designs [default %default]'),
  optparse::make_option('--deg_file',default='inputs/all_samples.by_clone.deg.tsv.gz',help='Original Seurat-style DEG TSV [default %default]'),
  optparse::make_option('--assay',default='RNA'),
  optparse::make_option('--layer',default='data'),
  optparse::make_option('--cell_column',default='cell_id',help='Metadata cell ID field, or rownames [default %default]'),
  optparse::make_option('--padj_threshold',type='double',default=0.05),
  optparse::make_option('--logfc_threshold',type='double',default=0.5),
  optparse::make_option('--pct1_threshold',type='double',default=0.10),
  optparse::make_option('--top_n',type='integer',default=100L,help='Prepare the original Top-N flagged genes unless explicit genes supplied [default %default]'),
  optparse::make_option('--genes',default=NULL,type='character',help='Optional comma-separated genes, e.g. CLDN16; overrides Top-N extraction, not DEG ranking'),
  optparse::make_option('--gene_list',default=NULL,type='character',help='Optional plain text file: one gene per line, no header; alternative to --genes'),
  optparse::make_option('--expression_scale',default='Normalized expression; log base and preprocessing provenance unverified',help='Record known normalization/response scale; no transformation is performed'),
  optparse::make_option('--output_prefix',default='results/module03')
)
opt <- optparse::parse_args(optparse::OptionParser(option_list=options_list))
if(!is.null(opt$genes) && !is.null(opt$gene_list)) stop('Choose --genes OR --gene_list, not both.')
if(!is.finite(opt$padj_threshold) || opt$padj_threshold<=0 || opt$padj_threshold>1 ||
   !is.finite(opt$logfc_threshold) || opt$logfc_threshold<0 || !is.finite(opt$pct1_threshold) ||
   opt$pct1_threshold<0 || opt$pct1_threshold>1 || is.na(opt$top_n) || opt$top_n<1L) stop('Invalid DEG thresholds or top_n.')
if(opt$layer=='counts') stop('Use normalized expression for this module, not the counts layer. No normalization is performed here.')
for(p in c(opt$seurat_object,opt$tree_designs,opt$gene_list)) if(!file.exists(p)) stop('Missing input: ',p)
explicit <- NULL
if(!is.null(opt$genes)) explicit <- trimws(strsplit(opt$genes,',',fixed=TRUE)[[1L]])
if(!is.null(opt$gene_list)) explicit <- trimws(readLines(opt$gene_list,warn=FALSE))
if(!is.null(explicit)) {
  explicit <- unique(explicit[nzchar(explicit)])
  if(!length(explicit)) stop('Explicit gene list is empty.')
}
if(!file.exists(opt$deg_file) && is.null(explicit)) stop('DEG file missing; supply a DEG table or an explicit gene list.')
outputs <- paste0(opt$output_prefix,c('.report.txt','.candidate_genes.tsv','.filtered_deg.tsv.gz',
  '.gene_qc.tsv','.clone_expression_summary.tsv','.expression.rds'))
if(any(file.exists(outputs))) stop('Output already exists; use a new --output_prefix. No overwrite performed.')
if(!dir.exists(dirname(opt$output_prefix))) dir.create(dirname(opt$output_prefix),recursive=TRUE)
log <- file(paste0(opt$output_prefix,'.report.txt'),'wt');sink(log,split=TRUE)
exit_code <- 0L
tryCatch({
  cat('MODULE 3 — CANDIDATE GENES AND EXPRESSION\n',format(Sys.time(),tz='UTC',usetz=TRUE),'\n\n')
  cat('Original DEG filters/ranking reused; no fitting, imputation, or normalization.\n')
  if(file.exists(opt$deg_file)) {
    deg <- readr::read_tsv(opt$deg_file,show_col_types=FALSE)
    required <- c('p_val','avg_log2FC','pct.1','pct.2','p_val_adj','clone_id','gene')
    if(length(setdiff(required,names(deg)))) stop('DEG table missing columns: ',paste(setdiff(required,names(deg)),collapse=', '))
    for(n in c('p_val','avg_log2FC','pct.1','pct.2','p_val_adj')) {
      deg[[n]] <- suppressWarnings(as.numeric(deg[[n]]))
      if(any(!is.finite(deg[[n]]))) stop('Nonfinite or nonnumeric DEG column: ',n)
      if(n!='avg_log2FC' && any(deg[[n]]<0 | deg[[n]]>1)) stop('DEG probability/fraction outside [0,1]: ',n)
    }
    if(anyNA(deg$clone_id) || any(!nzchar(trimws(as.character(deg$clone_id))))) stop('Missing DEG clone labels.')
    cat('Input DEG rows:',nrow(deg),'| unique genes:',dplyr::n_distinct(deg$gene),'\n')
    ranked <- rank_candidates(deg,opt)
    candidates <- ranked$candidates
    readr::write_tsv(ranked$filtered,paste0(opt$output_prefix,'.filtered_deg.tsv.gz'))
    cat('Passing rows:',nrow(ranked$filtered),'| candidate genes:',nrow(candidates),'\n')
    cat('Original priority score is a ranking heuristic, not a significance test.\n')
  } else {
    cat('No DEG table available: preparing explicitly supplied genes without DEG selection.\n')
    candidates <- data.frame(gene=character(),candidate_rank=integer())
  }
  requested <- if(!is.null(explicit)) explicit else candidates$gene[candidates$suggested_for_initial_elastic_net]
  if(!length(requested)) stop('No genes requested for extraction.')
  prepared_tree <- readRDS(opt$tree_designs)
  cells <- prepared_tree$cell_metadata
  if(!is.data.frame(cells) || !all(c('cell_id','clone_id') %in% names(cells))) stop('Invalid Module 2 cell metadata.')
  check_ids <- function(ids,label) {
    if(!length(ids) || anyNA(ids) || any(!nzchar(trimws(ids))) || anyDuplicated(ids)) stop(label,' must be nonempty unique IDs.')
  }
  check_ids(cells$cell_id,'Module 2 cell IDs')
  if(anyNA(cells$clone_id) || any(!nzchar(trimws(cells$clone_id)))) stop('Invalid saved clone labels.')
  if(!length(prepared_tree$designs)) stop('No designs in Module 2 file.')
  for(mode in names(prepared_tree$designs)) if(!identical(rownames(prepared_tree$designs[[mode]]$cell_X),cells$cell_id))
    stop('Saved design row order disagrees with Module 2 metadata for ',mode)
  cat('\nLoading Seurat expression once; retaining Module 2 cell order.\n')
  obj <- readRDS(opt$seurat_object)
  if(!inherits(obj,'Seurat')) stop('Expression object is not Seurat.')
  if(!opt$assay %in% names(obj@assays)) stop('Assay absent: ',opt$assay)
  meta <- obj[[]];check_ids(rownames(meta),'Metadata row names')
  if(opt$cell_column=='rownames') ids <- rownames(meta) else if(opt$cell_column %in% names(meta))
    ids <- as.character(meta[[opt$cell_column]]) else if(opt$cell_column=='cell_id') {
      cat('cell_id absent: using metadata row names (original fallback).\n');ids <- rownames(meta)
    } else stop('Requested cell-column absent: ',opt$cell_column)
  check_ids(ids,'Canonical expression cell IDs')
  cells$expression_cell <- rownames(meta)[match(cells$cell_id,ids)]
  if(anyNA(cells$expression_cell)) stop('Module 2 cells missing from Seurat metadata. Resolve alignment; do not silently drop tree cells.')
  aa <- obj[[opt$assay]]
  if(inherits(aa,'Assay5')) {
    available <- SeuratObject::Layers(aa,search=NA)
    if(!opt$layer %in% available) stop('Exact requested layer absent. Available: ',paste(available,collapse=', '),'. Split layers are not joined automatically.')
    mat <- SeuratObject::LayerData(aa,layer=opt$layer)
  } else {
    mat <- if(utils::packageVersion('SeuratObject')>='5.0.0') SeuratObject::GetAssayData(aa,layer=opt$layer)
      else SeuratObject::GetAssayData(aa,slot=opt$layer)
  }
  check_ids(rownames(mat),'Layer gene names');check_ids(colnames(mat),'Layer cell names')
  if(!all(cells$expression_cell %in% colnames(mat))) stop('Selected layer does not contain all saved tree cells. Review split/partial layers before modeling.')
  if(!is.numeric(mat) && !inherits(mat,'dMatrix') && !inherits(mat,'iMatrix')) stop('Selected layer must be numeric.')
  present <- requested %in% rownames(mat)
  qc <- data.frame(gene=requested,requested_order=seq_along(requested),in_assay_layer=present,
    candidate_rank=candidates$candidate_rank[match(requested,candidates$gene)],
    n_cells=nrow(cells),n_nonzero=NA_integer_,mean=NA_real_,median=NA_real_,sd=NA_real_,
    min=NA_real_,max=NA_real_,status=ifelse(present,'pending','missing_from_layer'))
  groups <- split(seq_len(nrow(cells)),cells$clone_id)
  clone_summaries <- list()
  for(i in which(present)) {
    v <- as.numeric(mat[requested[i],cells$expression_cell,drop=TRUE])
    if(any(!is.finite(v))) {qc$status[i] <- 'nonfinite';next}
    qc$n_nonzero[i] <- sum(v!=0);qc$mean[i] <- mean(v);qc$median[i] <- stats::median(v)
    qc$sd[i] <- stats::sd(v);qc$min[i] <- min(v);qc$max[i] <- max(v)
    qc$status[i] <- if(length(v)<2L || is.na(qc$sd[i]) || qc$sd[i]==0) 'constant' else 'ready'
    for(cl in names(groups)) {
      vv <- v[groups[[cl]]]
      clone_summaries[[length(clone_summaries)+1L]] <- data.frame(gene=requested[i],clone_id=cl,
        n_cells=length(vv),mean=mean(vv),median=stats::median(vv),fraction_nonzero=mean(vv!=0))
    }
  }
  # Report missing/nonfinite/constant genes explicitly; they are never replaced by zeros.
  readr::write_tsv(qc,paste0(opt$output_prefix,'.gene_qc.tsv'))
  clone_summary <- if(length(clone_summaries)) dplyr::bind_rows(clone_summaries) else
    data.frame(gene=character(),clone_id=character(),n_cells=integer(),mean=numeric(),median=numeric(),fraction_nonzero=numeric())
  readr::write_tsv(clone_summary,paste0(opt$output_prefix,'.clone_expression_summary.tsv'))
  genes <- qc$gene[qc$status=='ready']
  candidates$in_assay_layer <- candidates$gene %in% rownames(mat)
  candidates$requested_for_expression <- candidates$gene %in% requested
  candidates$prepared_for_modeling <- candidates$gene %in% genes
  readr::write_tsv(candidates,paste0(opt$output_prefix,'.candidate_genes.tsv'))
  cat('\nGene preparation status:\n');print(table(qc$status))
  if(any(qc$status!='ready')) {cat('Excluded from modeling (see gene_qc.tsv):\n');print(qc[qc$status!='ready',c('gene','status')],row.names=FALSE)}
  if(!length(genes)) stop('No variable finite genes remain; QC tables written for review.')
  expression <- mat[genes,cells$expression_cell,drop=FALSE]
  colnames(expression) <- cells$cell_id
  stopifnot(identical(colnames(expression),cells$cell_id))
  outside_clones <- if(exists('ranked')) setdiff(unique(as.character(ranked$filtered$clone_id)),unique(cells$clone_id)) else character()
  if(length(outside_clones)) cat('REVIEW: DEG evidence includes clones outside the saved tree:',paste(outside_clones,collapse=', '),
    '\nRanking is preserved; this mismatch needs review before interpretation.\n')
  saveRDS(list(expression=expression,cell_metadata=cells,gene_qc=qc,candidate_genes=candidates,
    clone_expression_summary=clone_summary,options=opt,tree_designs_path=normalizePath(opt$tree_designs),
    tree_designs_md5=unname(tools::md5sum(opt$tree_designs)),expression_scale=opt$expression_scale,
    sessionInfo=sessionInfo()),paste0(opt$output_prefix,'.expression.rds'))
  cat('\nPrepared expression:',nrow(expression),'genes x',ncol(expression),'cells\n')
  cat('Expression class:',paste(class(expression),collapse='/'),'\nAssay/layer:',opt$assay,'/',opt$layer,'\n')
  cat('Response scale:',opt$expression_scale,'\n')
  cat('Zeros preserved; no imputation or transformation. Regression coefficients are not automatically DEG log2FC.\n')
  cat('Expression uses exact Module 2 cell order for both tree modes.\n')
  cat('\nTop-ranked candidates (up to 10):\n');print(utils::head(candidates[,intersect(c('candidate_rank','gene','candidate_priority_score'),names(candidates)),drop=FALSE],10))
  cat('\nOptions used:\n');print(opt)
  cat('\nFiles to return:\n');cat(paste(outputs[file.exists(outputs)],collapse='\n'),'\n')
  cat('\nEXPRESSION PREPARATION COMPLETE — review before running Module 4.\n')
},error=function(e) {
  cat('\nEXPRESSION PREPARATION FAILED:',conditionMessage(e),'\nReview partial report/QC and terminal error.\n')
  exit_code <<- 1L
})
sink();close(log)
if(exit_code) message('Expression preparation failed; see ',opt$output_prefix,'.report.txt') else message('Report: ', opt$output_prefix,'.report.txt')
quit(status=exit_code)
