#!/usr/bin/env Rscript
# Prepare reusable topology-only branch designs. No expression models are fitted.

for (pkg in c('optparse','ape','Matrix')) if (!requireNamespace(pkg, quietly=TRUE))
  stop('Missing package ', pkg, '. Run Rscript 00_install_packages.R in this environment.')
options_list <- list(
  optparse::make_option('--concordant_clones', default='inputs/MERGED.Ovarian-cancer.fasttree.manifold.concordant.clones.tsv.gz'),
  optparse::make_option('--clonal_tree', default='inputs/MERGED.Ovarian-cancer.fasttree.manifold.clone.treefile'),
  optparse::make_option('--tree_mode', default='both', help='both, collapsed, or uncollapsed [default %default]'),
  optparse::make_option('--output_prefix', default='results/module02'),
  optparse::make_option('--min_cells_per_clone', type='integer', default=10L),
  optparse::make_option('--min_descendant_cells', type='integer', default=25L),
  optparse::make_option('--max_descendant_fraction', type='double', default=0.90),
  optparse::make_option('--internal_edges_only', default='TRUE', help='Uncollapsed filter only [default %default]')
)
opt <- optparse::parse_args(optparse::OptionParser(option_list=options_list))
if (!opt$tree_mode %in% c('both','collapsed','uncollapsed')) stop('Invalid --tree_mode')
if (!tolower(opt$internal_edges_only) %in% c('true','false')) stop('--internal_edges_only must be TRUE or FALSE')
internal_only <- tolower(opt$internal_edges_only) == 'true'
if (is.na(opt$min_cells_per_clone) || opt$min_cells_per_clone < 1L ||
    is.na(opt$min_descendant_cells) || opt$min_descendant_cells < 2L ||
    !is.finite(opt$max_descendant_fraction) || opt$max_descendant_fraction <= 0 || opt$max_descendant_fraction >= 1)
  stop('Invalid size filters: positive min_cells_per_clone; min_descendant_cells >=2; max_descendant_fraction between 0 and 1.')
for (p in c(opt$concordant_clones,opt$clonal_tree)) if (!file.exists(p)) stop('Input missing: ',p)
want_collapsed <- opt$tree_mode %in% c('both','collapsed')
want_uncollapsed <- opt$tree_mode %in% c('both','uncollapsed')
suffixes <- c('.report.txt','.clone_checks.tsv','.branches.tsv','.tree_designs.rds',
              if(want_collapsed) c('.collapsed_tree.nwk','.collapsed_tree.png'))
outputs <- paste0(opt$output_prefix,suffixes)
if (any(file.exists(outputs))) stop('Output already exists; choose a new --output_prefix. No outputs overwritten.')
if (!dir.exists(dirname(opt$output_prefix))) dir.create(dirname(opt$output_prefix),recursive=TRUE)
write_tsv <- function(d,suffix) utils::write.table(d,paste0(opt$output_prefix,suffix),sep='\t',quote=TRUE,row.names=FALSE,na='NA')
check_ids <- function(x,what) {
  if (!length(x) || anyNA(x) || any(!nzchar(trimws(x))) || anyDuplicated(x)) stop(what,' must be nonempty and unique.')
}
# Adapted original depth-first traversal. Cache children once, rather than per edge.
get_descendant_tip_numbers <- function(tree,node,children_map) {
  n_tip <- length(tree$tip.label); out <- integer(); stack <- as.integer(node)
  while(length(stack)) {
    current <- stack[length(stack)]; stack <- stack[-length(stack)]
    if(current <= n_tip) out <- c(out,current) else stack <- c(stack,children_map[[as.character(current)]])
  }
  sort(unique(out))
}
make_sparse_design <- function(descendants, selected, row_ids, column_ids) {
  dd <- descendants[selected]
  Matrix::sparseMatrix(i=as.integer(unlist(dd,use.names=FALSE)), j=rep(seq_along(dd),lengths(dd)),
    x=rep(1,sum(lengths(dd))),dims=c(length(row_ids),length(selected)),dimnames=list(row_ids,column_ids))
}
# Keep a dictionary for every original edge, including exclusions and orientation.
canonicalize <- function(dictionary,descendants,n_tips) {
  seen <- new.env(parent=emptyenv())
  dictionary$equivalent_to <- NA_character_; dictionary$equivalent_sign <- NA_integer_
  for (i in which(dictionary$included)) {
    key <- paste(descendants[[i]],collapse=',')
    complement <- paste(setdiff(seq_len(n_tips),descendants[[i]]),collapse=',')
    if (exists(key,envir=seen,inherits=FALSE)) {
      dictionary$included[i] <- FALSE; dictionary$reason[i] <- 'duplicate_partition'
      dictionary$equivalent_to[i] <- get(key,envir=seen); dictionary$equivalent_sign[i] <- 1L
    } else if (exists(complement,envir=seen,inherits=FALSE)) {
      dictionary$included[i] <- FALSE; dictionary$reason[i] <- 'complement_with_intercept'
      dictionary$equivalent_to[i] <- get(complement,envir=seen); dictionary$equivalent_sign[i] <- -1L
    } else assign(key,dictionary$branch_id[i],envir=seen)
  }
  dictionary
}
log <- file(paste0(opt$output_prefix,'.report.txt'),'wt');sink(log,split=TRUE)
exit_code <- 0L
tryCatch({
  cat('MODULE 2 — TREE PREPARATION\n',format(Sys.time(),tz='UTC',usetz=TRUE),'\n\n')
  cat('Source: original branch traversal, filtering, and sparse descendant encoding.\n')
  cat('No Seurat RDS loaded; expression-to-tree alignment is rechecked when expression is prepared.\n')
  cat('Using supplied root; no rerooting or branch-support threshold imposed.\n\n')
  clones <- utils::read.delim(opt$concordant_clones,colClasses='character',check.names=FALSE)
  if(!all(c('cell_id','clone_id') %in% names(clones))) stop('Clone table requires cell_id and clone_id.')
  check_ids(clones$cell_id,'Clone-table cell IDs')
  if(anyNA(clones$clone_id) || any(!nzchar(trimws(clones$clone_id)))) stop('Missing clone labels.')
  tree <- ape::read.tree(opt$clonal_tree)
  if(!inherits(tree,'phylo')) stop('Exactly one Newick tree required.')
  check_ids(tree$tip.label,'Tree tips')
  if(!ape::is.rooted(tree)) stop('Tree is structurally unrooted. Supply an explicitly rooted tree; this script does not guess a root.')
  if(length(setdiff(clones$cell_id,tree$tip.label))) stop('Clone-table cells missing from tree; resolve input alignment before preparing designs.')
  clone_counts <- table(clones$clone_id)
  retained_clones <- names(clone_counts)[clone_counts >= opt$min_cells_per_clone]
  removed_clones <- setdiff(names(clone_counts),retained_clones)
  kept <- clones$cell_id[clones$clone_id %in% retained_clones]
  if(length(retained_clones)<2L) stop('Fewer than two clones meet min_cells_per_clone.')
  cat('Input tree tips:',length(tree$tip.label),'\nClone-table cells:',nrow(clones),'\n')
  cat('Tree tips outside retained clone cells:',length(setdiff(tree$tip.label,kept)),'\n')
  cat('Clones excluded by size:',if(length(removed_clones)) paste(removed_clones,collapse=', ') else 'none','\n')
  if(!setequal(tree$tip.label,kept)) tree <- ape::keep.tip(tree,kept)
  tree <- ape::reorder.phylo(tree,order='postorder')
  if(!ape::is.rooted(tree)) stop('Pruning removed the root split; provide/review an explicitly rooted retained tree.')
  cell_metadata <- data.frame(cell_id=tree$tip.label,clone_id=clones$clone_id[match(tree$tip.label,clones$cell_id)])
  n <- length(tree$tip.label); nt <- table(cell_metadata$clone_id)
  cm <- split(tree$edge[,2],tree$edge[,1])
  descendants <- lapply(tree$edge[,2],function(node) get_descendant_tip_numbers(tree,node,cm))
  root <- setdiff(tree$edge[,1],tree$edge[,2])[1L]
  clone_checks <- do.call(rbind,lapply(sort(names(nt)),function(cl) {
    ids <- which(cell_metadata$clone_id==cl)
    node <- if(length(ids)==1L) ids else ape::getMRCA(tree,ids)
    below <- get_descendant_tip_numbers(tree,node,cm)
    support <- if(node>n && !is.null(tree$node.label)) tree$node.label[node-n] else NA_character_
    data.frame(clone_id=cl,n_cells=length(ids),mrca_node=node,n_cells_below_mrca=length(below),
      monophyletic=setequal(ids,below),mrca_label_raw=support)
  }))
  write_tsv(clone_checks,'.clone_checks.tsv');print(clone_checks,row.names=FALSE)
  # All branches retaining complete clones, rather than applying cellular size filters
  # to clone tips. Never force nonmonophyletic clones into a collapsed tip.
  comp <- lapply(descendants,function(ii) table(cell_metadata$clone_id[ii]))
  complete <- vapply(comp,function(tb) all(as.integer(tb)==as.integer(nt[names(tb)])),FALSE)
  stem_labels <- if(is.null(tree$node.label)) rep(NA_character_,nrow(tree$edge)) else
    vapply(tree$edge[,2],function(node) if(node>n) tree$node.label[node-n] else NA_character_,'')
  base <- data.frame(branch_id=paste0('branch_edge_',seq_len(nrow(tree$edge))),
    edge_index=seq_len(nrow(tree$edge)),parent_node=tree$edge[,1],child_node=tree$edge[,2],
    child_label_raw=stem_labels,
    branch_length=if(is.null(tree$edge.length)) NA_real_ else tree$edge.length,
    n_descendant_cells=lengths(descendants),descendant_fraction=lengths(descendants)/n,
    is_terminal_edge=tree$edge[,2]<=n,complete_clones_only=complete,
    descendant_clones=vapply(comp,function(tb) paste(sort(names(tb)),collapse=';'),''),
    splits_within_clone=vapply(comp,function(tb) any(as.integer(tb)<as.integer(nt[names(tb)])),FALSE))
  dictionaries <- list();designs <- list();clone_tree <- NULL
  if(want_collapsed) {
    if(!all(clone_checks$monophyletic)) stop('Nonmonophyletic clone(s): collapsed design not created. Review clone_checks.tsv or explicitly run --tree_mode uncollapsed.')
    dd <- base;dd$mode <- 'collapsed'
    dd$included <- complete & dd$n_descendant_cells<n
    dd$reason <- ifelse(dd$included,'retained',ifelse(!complete,'splits_clone','universal_partition'))
    dd <- canonicalize(dd,descendants,n)
    select <- which(dd$included)
    x <- make_sparse_design(descendants,select,tree$tip.label,dd$branch_id[select])
    clone_order <- sort(names(nt))
    clone_x <- x[match(clone_order,cell_metadata$clone_id),,drop=FALSE];rownames(clone_x) <- clone_order
    # Confirm every cell of a clone actually has its clone's entire predictor row.
    stopifnot(all(x == clone_x[match(cell_metadata$clone_id,clone_order),,drop=FALSE]))
    q <- qr(cbind('(Intercept)'=1,as.matrix(clone_x)))
    basis <- q$pivot[seq_len(q$rank)];basis <- basis[basis>1L]-1L
    designs$collapsed <- list(cell_X=x,clone_X=clone_x,rank_with_intercept=q$rank,
      independent_branch_terms=colnames(clone_x)[basis])
    dictionaries$collapsed <- dd
    # One representative per monophyletic clone preserves inter-clone topology.
    # Erase lengths/labels: representative paths include within-clone distances.
    representatives <- vapply(clone_order,function(cl) sort(cell_metadata$cell_id[cell_metadata$clone_id==cl])[1L],'')
    clone_tree <- ape::keep.tip(tree,representatives)
    clone_tree$tip.label <- clone_order[match(clone_tree$tip.label,representatives)]
    clone_tree$edge.length <- NULL;clone_tree$node.label <- NULL
    ape::write.tree(clone_tree,file=paste0(opt$output_prefix,'.collapsed_tree.nwk'))
    grDevices::png(paste0(opt$output_prefix,'.collapsed_tree.png'),width=1200,height=800,res=150)
    tryCatch({ape::plot.phylo(clone_tree,type='cladogram',use.edge.length=FALSE,
      main='Collapsed clone topology',label.offset=0.03);graphics::mtext('Supplied root; lengths omitted; support not validated',side=1,line=1)},
      finally=grDevices::dev.off())
    cat('\nCollapsed clone-level matrix (binary descendant membership):\n');print(as.matrix(clone_x))
    cat('Collapsed predictors:',ncol(x),'| independent branch directions:',q$rank-1L,'\n')
    cat('Linear-model basis:',paste(colnames(clone_x)[basis],collapse=', '),'\n')
    cat('Other retained branches may be linearly dependent; penalized selection does not establish unique historical shifts.\n')
  }
  if(want_uncollapsed) {
    dd <- base;dd$mode <- 'uncollapsed'
    dd$included <- (!internal_only | !dd$is_terminal_edge) & dd$n_descendant_cells>=opt$min_descendant_cells &
      dd$descendant_fraction<=opt$max_descendant_fraction & dd$n_descendant_cells<n
    dd$reason <- ifelse(dd$included,'retained','original_cellular_edge_filters')
    n_legacy <- sum(dd$included)
    dd <- canonicalize(dd,descendants,n)
    select <- which(dd$included)
    x <- make_sparse_design(descendants,select,tree$tip.label,dd$branch_id[select])
    rank_x <- as.integer(Matrix::rankMatrix(cbind(Matrix::Matrix(1,nrow=n,ncol=1,sparse=TRUE),x),method='qr.R'))
    designs$uncollapsed <- list(cell_X=x,rank_with_intercept=rank_x)
    dictionaries$uncollapsed <- dd
    cat('\nUncollapsed original-filter candidates:',n_legacy,'| after duplicate/complement checks:',ncol(x),
      '| rank with intercept:',rank_x,'\n')
    cat('Uncollapsed within-clone predictors remain exploratory on this supplied topology.\n')
  }
  dictionary <- do.call(rbind,dictionaries);rownames(dictionary) <- NULL
  write_tsv(dictionary,'.branches.tsv')
  saveRDS(list(tree=tree,cell_metadata=cell_metadata,clone_tree=clone_tree,clone_checks=clone_checks,
    branch_dictionary=dictionary,descendant_tip_indices=descendants,designs=designs,options=opt,
    root_provenance='supplied root; biological provenance not verified',
    support_provenance='raw node labels retained; meaning/scale not verified',sessionInfo=sessionInfo()),
    paste0(opt$output_prefix,'.tree_designs.rds'))
  cat('\nRetained tree tips:',n,'| edges:',nrow(tree$edge),'| internal nodes:',tree$Nnode,'\n')
  cat('Root node in saved tree:',root,'\n')
  cat('Raw internal labels:',if(is.null(tree$node.label)) 'absent' else 'present; inspect branches.tsv','\n')
  cat('Node/edge IDs refer ONLY to the saved pruned tree, not to differently ordered input files.\n')
  cat('Branch membership is stable through saved descendant cell sets, not renumbering after resampling.\n')
  cat('No support threshold applied: supported biological resolution requires root/support provenance review.\n')
  cat('\nOptions used:\n');print(opt)
  cat('\nOutput files:\n');cat(paste(outputs,collapse='\n'),'\n')
  cat('\nTREE PREPARATION COMPLETE — review before running Module 3.\n')
},error=function(e) {
  cat('\nTREE PREPARATION FAILED:',conditionMessage(e),'\nReview partial outputs and the terminal error.\n')
  exit_code <<- 1L
})
sink();close(log)
if(exit_code) message('Tree preparation failed; see ',opt$output_prefix,'.report.txt') else message('Report: ', opt$output_prefix,'.report.txt')
quit(status=exit_code)
