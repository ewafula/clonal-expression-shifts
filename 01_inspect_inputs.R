#!/usr/bin/env Rscript
# Module 1: inspect existing inputs and write ONE text report. No analysis/fitting.
# Run from clonal-expression-shifts/; all input paths can be overridden.
defaults <- list(
  seurat_object = 'inputs/SPECTRUM-cell-ploidy-prediction.rds',
  concordant_clones = 'inputs/MERGED.Ovarian-cancer.fasttree.manifold.concordant.clones.tsv.gz',
  clonal_tree = 'inputs/MERGED.Ovarian-cancer.fasttree.manifold.clone.treefile',
  assay = 'RNA', cell_column = 'cell_id', sample_id = NULL,
  output = 'results/module01_input_report.txt'
)
args <- commandArgs(trailingOnly = TRUE)
if ('--help' %in% args) {
  cat('Usage: Rscript 01_inspect_inputs.R [--seurat_object FILE] [--concordant_clones FILE]\n',
      ' [--clonal_tree FILE] [--assay RNA] [--cell_column cell_id|rownames]\n',
      ' [--sample_id COLUMN] [--output results/module01_input_report.txt]\n',
      'Defaults match the existing dataset. Only one text report is written.\n', sep='')
  quit(status = 0)
}
if (length(args) %% 2L) stop('Arguments must be --name value pairs; use --help.')
opt <- defaults
if (length(args)) for (i in seq(1L, length(args), by = 2L)) {
  key <- sub('^--', '', args[i])
  if (!startsWith(args[i], '--') || !key %in% names(defaults)) stop('Unknown option: ', args[i])
  opt[[key]] <- args[i+1L]
}
if (file.exists(opt$output)) stop('Report already exists. Choose a new --output filename: ', opt$output)
if (!dir.exists(dirname(opt$output)) && !dir.create(dirname(opt$output), recursive=TRUE))
  stop('Cannot create output directory: ', dirname(opt$output))
con <- file(opt$output, 'wt')
sink(con, split = TRUE)
exit_code <- 0L
tryCatch({
  cat('MODULE 1 — INPUT INSPECTION\n', format(Sys.time(), tz='UTC', usetz=TRUE), '\n\n')
  cat('Purpose: inspect inputs only; no normalization, filtering, tree collapse, or model fitting.\n\n')
  for (p in c('SeuratObject', 'ape')) {
    if (!requireNamespace(p, quietly=TRUE)) stop('Required package missing: ', p, '. Run Rscript 00_install_packages.R in this R environment.')
    cat(p, 'version:', as.character(utils::packageVersion(p)), '\n')
  }
  for (key in c('seurat_object','concordant_clones','clonal_tree')) {
    if (!file.exists(opt[[key]])) stop('Missing input: ', opt[[key]])
    cat(key, ': ', opt[[key]], '\n', sep='')
  }
  obj <- readRDS(opt$seurat_object)
  if (!inherits(obj, 'Seurat')) stop('Expression RDS must contain a Seurat object.')
  meta <- obj[[]]
  cat('\nEXPRESSION OBJECT\nCells:', ncol(obj), '\n')
  for (a in names(obj@assays)) {
    aa <- obj[[a]]
    layers <- if (inherits(aa, 'Assay5')) SeuratObject::Layers(aa, search=NA) else
      intersect(c('counts','data','scale.data'), methods::slotNames(aa))
    cat('Assay:', a, '| features:', nrow(aa), '| cells:', ncol(aa),
        '| layers/slots:', paste(layers, collapse=', '), '\n')
    if (!inherits(aa, 'Assay5')) cat('  Legacy slots listed above may be empty.\n')
  }
  if (!opt$assay %in% names(obj@assays)) stop('Requested assay absent: ', opt$assay)
  aa <- obj[[opt$assay]]
  if (inherits(aa, 'Assay5')) {
    layers <- SeuratObject::Layers(aa, search=NA)
    if (!'data' %in% layers) cat('REVIEW: no exact data layer; split/absent normalized layers need review before modeling.\n')
    if (!'counts' %in% layers) cat('NOTE: no exact counts layer; inspect layer names above.\n')
  }
  cat('Assay to use:', opt$assay, '\nCLDN16 present:', 'CLDN16' %in% rownames(aa), '\n')
  cat('Layer names do not establish normalization scale or raw-count provenance.\n')
  cat('\nMETADATA COLUMNS\n')
  print(data.frame(column=names(meta), type=vapply(meta, function(x) paste(class(x), collapse='/'), ''),
    missing=vapply(meta, function(x) sum(is.na(x)), 0L)), row.names=FALSE)
  candidates <- names(meta)[grepl('sample|patient|batch|orig.ident', names(meta), ignore.case=TRUE)]
  cat('\nSAMPLE-RELATED METADATA VALUES\n')
  if (!length(candidates)) cat('No obvious sample/patient/batch field names found.\n')
  for (n in candidates) {
    z <- as.character(meta[[n]])
    cat('\n', n, ': ', length(unique(z[!is.na(z)])), ' distinct nonmissing values\n', sep='')
    print(utils::head(sort(table(z, useNA='ifany'), decreasing=TRUE), 50L))
    if (length(unique(z)) > 50L) cat('(Only the 50 most frequent values shown.)\n')
  }
  check_ids <- function(x, label) {
    if (!length(x) || anyNA(x) || any(!nzchar(trimws(x))) || anyDuplicated(x))
      stop(label, ' must contain nonempty unique IDs; duplicates are not silently removed.')
    x
  }
  object_names <- check_ids(colnames(obj), 'Seurat cell names')
  if (!all(object_names %in% rownames(meta))) stop('Seurat cell names do not match metadata row names.')
  meta <- meta[object_names, , drop=FALSE]
  if (opt$cell_column == 'rownames') ids <- rownames(meta) else if (opt$cell_column %in% names(meta))
    ids <- as.character(meta[[opt$cell_column]]) else if (opt$cell_column == 'cell_id') {
      cat('\nNOTE: cell_id absent; using Seurat cell names, as in the original script.\n')
      ids <- object_names
    } else stop('Cell column absent: ', opt$cell_column)
  ids <- check_ids(ids, 'Canonical Seurat cell IDs')
  clones <- utils::read.delim(opt$concordant_clones, colClasses='character', check.names=FALSE)
  if (!all(c('cell_id','clone_id') %in% names(clones))) stop('Clone table needs cell_id and clone_id columns.')
  clones$cell_id <- check_ids(clones$cell_id, 'Clone-table cell IDs')
  if (anyNA(clones$clone_id) || any(!nzchar(trimws(clones$clone_id)))) stop('Missing clone labels.')
  tree <- ape::read.tree(opt$clonal_tree)
  if (!inherits(tree, 'phylo')) stop('Expected exactly one Newick tree.')
  tips <- check_ids(tree$tip.label, 'Tree tips')
  assay_ids <- ids[match(colnames(aa), object_names)]
  if (anyNA(assay_ids)) stop('Assay cells not found in object cell names.')
  common <- Reduce(intersect, list(assay_ids, clones$cell_id, tips))
  cat('\nCELL OVERLAP\nSeurat object:', length(ids), '\nSelected assay:', length(assay_ids),
    '\nClone table:', nrow(clones), '\nTree:', length(tips), '\nShared across selected assay/clone table/tree:', length(common), '\n')
  for (source in c('assay','clone_table','tree')) {
    x <- switch(source, assay=assay_ids, clone_table=clones$cell_id, tree=tips)
    missing <- setdiff(x, common)
    cat(source, 'IDs outside shared set:', length(missing), '\n')
    if (length(missing)) cat('  Examples:', paste(utils::head(missing,10), collapse=', '), '\n')
  }
  if (!length(common)) stop('No shared cells.')
  m <- meta[match(common, ids), , drop=FALSE]
  clone <- clones$clone_id[match(common, clones$cell_id)]
  cat('\nCLONE COUNTS — SHARED CELLS\n'); print(table(clone))
  if (!is.null(opt$sample_id)) {
    if (!opt$sample_id %in% names(m)) stop('Sample ID column absent: ', opt$sample_id)
    sample <- as.character(m[[opt$sample_id]])
    sample[is.na(sample) | !nzchar(trimws(sample))] <- '<MISSING>'
    cat('\nCLONE BY SAMPLE COUNTS — SHARED CELLS\nMetadata field:', opt$sample_id, '\n')
    print(table(sample_id=sample, clone=clone))
  } else cat('\nSample grouping not requested; clone/cell alignment does not require sample IDs.\n')
  cat('\nREVIEW BEFORE NEXT MODULE\n',
    '- Clone IDs and cell IDs are the required grouping/alignment identifiers.\n',
    '- If sample grouping was requested, confirm its metadata meaning; no samples are inferred.\n',
    '- Confirm normalized-expression scale and preprocessing provenance.\n',
    '- Review cell overlap losses above; no input was pruned or overwritten.\n',
    '- This inspection does not validate every layer-specific cell set or every expression value.\n',
    '- Root/support and branch designs belong in Module 2.\n',
    '\nINSPECTION COMPLETE — review before running Module 2.\n', sep='')
}, error=function(e) {
  cat('\nINSPECTION FAILED:', conditionMessage(e), '\nReview this partial report and the terminal error.\n')
  exit_code <<- 1L
})
sink(); close(con)
if (exit_code) message('Inspection failed; see ', opt$output) else message('Report: ', opt$output)
quit(status=exit_code)
