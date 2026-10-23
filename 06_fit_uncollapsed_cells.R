#!/usr/bin/env Rscript
# Uncollapsed exploratory pilot. Adapted from Module 4 and original tree-branch fitting.
# No sample adjustment, resampling stability, or confirmatory p-values.

for(pkg in c('optparse','Matrix','glmnet','ape')) if(!requireNamespace(pkg,quietly=TRUE))
  stop('Missing package ',pkg,'. Run Rscript 00_install_packages.R in this R environment.')
options_list <- list(
  optparse::make_option('--expression',default='results/module03.expression.rds',help='Module 3 prepared expression [default %default]'),
  optparse::make_option('--tree_designs',default='results/module02.tree_designs.rds'),
  optparse::make_option('--gene',default='CLDN16',help='One prepared gene for the pilot [default %default]'),
  optparse::make_option('--model',default='all',help='all, elastic_net, linear, or contrasts [default %default]'),
  optparse::make_option('--alpha',type='double',default=0.5),
  optparse::make_option('--nfolds',type='integer',default=5L),
  optparse::make_option('--lambda_choice',default='lambda.1se',help='lambda.1se or lambda.min [default %default]'),
  optparse::make_option('--seed',type='integer',default=12345L),
  optparse::make_option('--min_contrast_cells',type='integer',default=10L,help='Minimum cells on each side within a clone [default %default]'),
  optparse::make_option('--output_prefix',default='results/module06')
)
opt <- optparse::parse_args(optparse::OptionParser(option_list=options_list))
if(is.na(opt$min_contrast_cells) || opt$min_contrast_cells<2L)stop('min_contrast_cells must be >=2.')
if(!opt$model %in% c('all','elastic_net','linear','contrasts')) stop('Invalid --model.')
if(!is.finite(opt$alpha) || opt$alpha<=0 || opt$alpha>1) stop('--alpha must be >0 and <=1 for sparse Elastic Net/lasso.')
if(is.na(opt$nfolds) || opt$nfolds<3L || is.na(opt$seed) || opt$seed<0L) stop('Use nfolds>=3 and a nonnegative integer seed.')
if(!opt$lambda_choice %in% c('lambda.1se','lambda.min')) stop('Invalid --lambda_choice.')
for(p in c(opt$expression,opt$tree_designs)) if(!file.exists(p)) stop('Missing input: ',p)
files <- paste0(opt$output_prefix,c('.report.txt','.coefficients.tsv','.contrasts.tsv','.fit_metrics.tsv',
  '.cell_predictions.tsv.gz','.models.rds','.comparison.png','.selected_tree.png'))
if(any(file.exists(files))) stop('Output exists; select a new --output_prefix. No overwrite performed.')
if(!dir.exists(dirname(opt$output_prefix))) dir.create(dirname(opt$output_prefix),recursive=TRUE)
write_tsv <- function(x,suffix) {
 p <- paste0(opt$output_prefix,suffix)
 con <- if(grepl('\\.gz$',p)) gzfile(p,'wt') else file(p,'wt')
 tryCatch(utils::write.table(x,con,sep='\t',quote=TRUE,row.names=FALSE,na='NA'),finally=close(con))
}
log <- file(paste0(opt$output_prefix,'.report.txt'),'wt');sink(log,split=TRUE)
exit_code <- 0L
tryCatch({
  cat('MODULE 6 — EXPLORATORY UNCOLLAPSED CELL-LEVEL PILOT\n',format(Sys.time(),tz='UTC',usetz=TRUE),'\n\n')
  e <- readRDS(opt$expression);t <- readRDS(opt$tree_designs)
  if(is.null(t$designs$uncollapsed)) stop('Module 2 file has no uncollapsed design.')
  if(is.null(e$tree_designs_md5) || !identical(e$tree_designs_md5,unname(tools::md5sum(opt$tree_designs))))
    stop('Module 3 expression was prepared against a different Module 2 file; use matching inputs.')
  if(!opt$gene %in% rownames(e$expression)) stop('Gene not prepared: ',opt$gene,'. Prepare it with Module 3 --genes first.')
  x <- t$designs$uncollapsed$cell_X;cells <- e$cell_metadata
  if(!identical(colnames(e$expression),rownames(x)) || !identical(cells$cell_id,rownames(x)) ||
     !identical(cells$clone_id,t$cell_metadata$clone_id)) stop('Expression/design cell order or clone labels disagree.')
  clone <- as.character(cells$clone_id)
  if(anyNA(clone) || any(!nzchar(clone)) || anyDuplicated(cells$cell_id)) stop('Invalid cell/clone IDs.')
  y <- as.numeric(e$expression[opt$gene,,drop=TRUE]);n <- length(y)
  if(any(!is.finite(y)) || n<3L || stats::sd(y)==0) stop('Gene response must be finite, variable, with >=3 cells.')
  if(ncol(x)<1L) stop('No uncollapsed branch predictors.')
  if(any(!x@x %in% c(0,1)) || any(Matrix::colSums(x)<=0 | Matrix::colSums(x)>=n)) stop('Uncollapsed design must contain varying binary predictors.')
  levels_clone <- sort(unique(clone));sizes <- table(factor(clone,levels=levels_clone))
  if(length(levels_clone)<2L) stop('At least two clones required.')
  weights <- 1/as.numeric(sizes[match(clone,levels_clone)]);weights <- weights/mean(weights)
  clone_means <- vapply(levels_clone,function(cl) mean(y[clone==cl]),0.0)
  names(clone_means) <- levels_clone
  reference <- unname(clone_means[clone]);intercept <- stats::weighted.mean(y,weights)
  # Same clone-balanced loss for baselines, Elastic Net, and linear fitting.
  metric <- function(label,pred) {
    weighted_sse <- sum(weights*(y-pred)^2);weighted_sst <- sum(weights*(y-intercept)^2)
    data.frame(model=label,n_cells=n,n_clones=length(levels_clone),
      weighted_rmse=sqrt(weighted_sse/sum(weights)),weighted_mae=sum(weights*abs(y-pred))/sum(weights),
      weighted_r_squared_apparent=1-weighted_sse/weighted_sst,
      rmse=sqrt(mean((y-pred)^2)),mae=mean(abs(y-pred)),
      r_squared_apparent=1-sum((y-pred)^2)/sum((y-mean(y))^2),
      evaluation='training fit; not independent validation')
  }
  predictions <- data.frame(cell_id=cells$cell_id,clone_id=clone,observed=y,weight=weights,
    fold_id=NA_integer_,intercept_only=intercept,clone_means_reference=reference)
  metrics <- list(metric('intercept_only',rep(intercept,n)),metric('clone_means_reference',reference))
  coefficient_rows <- list();models <- list();cv <- NULL;linear_fit <- NULL
  # A deterministic full-rank basis spans the same uncollapsed prediction space.
  q_all <- qr(cbind(1,as.matrix(x)))
  idx <- q_all$pivot[seq_len(q_all$rank)];idx <- idx[idx>1L]-1L
  basis <- colnames(x)[idx]
  if(n<=q_all$rank)stop('Uncollapsed design is saturated; restrict candidate branches.')
  cat('Gene:',opt$gene,'| cells:',n,'| clones:',length(levels_clone),'\n')
  cat('Response scale:',e$expression_scale,'\n')
  cat('Clone-balanced weights: inverse clone cell count, rescaled to mean 1.\n')
  cat('Each clone contributes equal total weight. These are balancing, not known precision weights.\n')
  cat('Retained branch predictors:',ncol(x),'| independent directions:',q_all$rank-1L,'\n')
  if(opt$model %in% c('all','elastic_net')) {
    if(ncol(x)<2L) stop('glmnet requires >=2 predictor columns. Choose --model linear or contrasts.')
    k <- min(opt$nfolds,min(sizes))
    if(k<3L) stop('Insufficient cells per clone for >=3 stratified folds. Choose linear/contrasts or review the data.')
    set.seed(opt$seed)
    fold <- integer(n)
    for(cl in levels_clone) {
      ii <- which(clone==cl);fold[ii] <- sample(rep(seq_len(k),length.out=length(ii)))
    }
    predictions$fold_id <- fold
    # Original cv.glmnet Gaussian fit; new weights and saved clone-stratified folds.
    cv <- glmnet::cv.glmnet(x=x,y=y,family='gaussian',alpha=opt$alpha,weights=weights,
      foldid=fold,nfolds=k,standardize=TRUE,intercept=TRUE,keep=TRUE)
    chosen_lambda <- cv[[opt$lambda_choice]]
    beta <- as.matrix(stats::coef(cv,s=opt$lambda_choice))
    pred <- as.numeric(stats::predict(cv,newx=x,s=opt$lambda_choice))
    coefficient_rows$elastic_net <- data.frame(model='elastic_net',term=rownames(beta),
      coefficient=as.numeric(beta[,1]),selected=as.numeric(beta[,1])!=0,
      interpretation='penalized simultaneous branch contribution; response units')
    predictions$elastic_net <- pred;predictions$residual_elastic_net <- y-pred
    mm <- metric('elastic_net',pred);metrics[[length(metrics)+1L]] <- mm
    models$elastic_net <- cv
    cat('\nElastic Net: alpha=',opt$alpha,' | ',opt$lambda_choice,'=',chosen_lambda,
      ' | selected branches=',sum(beta[rownames(beta)!='(Intercept)',1]!=0),'\n',sep='')
    cat('lambda.min:',cv$lambda.min,'| lambda.1se:',cv$lambda.1se,'| folds:',k,'\n')
    cat('CV is clone-stratified cell interpolation used for tuning; not held-out-clone/patient validation.\n')
    cat('Fold IDs and tuning out-of-fold predictions are saved, not presented as unbiased final test performance.\n')
  }
  if(opt$model %in% c('all','linear')) {
    dat <- data.frame(y=y,as.matrix(x[,basis,drop=FALSE]),check.names=FALSE)
    dat$balancing_weight <- weights
    linear_fit <- stats::lm(stats::reformulate(basis,response='y'),data=dat,weights=balancing_weight)
    if(anyNA(stats::coef(linear_fit))) stop('Unexpected aliased linear coefficients.')
    pred <- as.numeric(stats::fitted(linear_fit))
    coefficient_rows$linear <- data.frame(model='linear',term=names(stats::coef(linear_fit)),
      coefficient=as.numeric(stats::coef(linear_fit)),selected=NA,
      interpretation='unpenalized coefficient in saved full-rank basis; not variable selection')
    predictions$linear <- pred;predictions$residual_linear <- y-pred
    metrics[[length(metrics)+1L]] <- metric('linear',pred);models$linear <- linear_fit
    cat('\nLinear basis:',paste(basis,collapse=', '),'\n')
    cat('Linear basis spans the eligible cellular partitions; lower training error can reflect extra flexibility. Basis coefficients do not identify unique shifts.\n')
  }
  # Within-clone contrasts compare actual cells, not clone-average differences.
  contrast_rows <- list();membership <- list();ready_membership <- list()
  dd <- t$branch_dictionary[t$branch_dictionary$mode=='uncollapsed' & t$branch_dictionary$included,]
  if(opt$model %in% c('all','contrasts')) {
    for(j in seq_len(nrow(dd))) {
      inside_all <- t$descendant_tip_indices[[dd$edge_index[j]]]
      parent <- dd$parent_node[j];child <- dd$child_node[j]
      sibling_nodes <- t$tree$edge[t$tree$edge[,1]==parent & t$tree$edge[,2]!=child,2]
      sibling_edges <- match(sibling_nodes,t$tree$edge[,2])
      sisters <- sort(unique(unlist(t$descendant_tip_indices[sibling_edges],use.names=FALSE)))
      for(cl in levels_clone)for(type in c('within_clone','sister_within_clone')) {
        inside <- intersect(inside_all,which(clone==cl))
        outside <- if(type=='within_clone')setdiff(which(clone==cl),inside_all) else intersect(sisters,which(clone==cl))
        id <- paste(dd$branch_id[j],cl,type,sep='__')
        ni <- length(inside);no <- length(outside)
        state <- if(!ni || !no)'no_within_clone_split' else if(min(ni,no)<opt$min_contrast_cells)'insufficient_cells' else 'ready'
        mean_in <- if(ni)mean(y[inside]) else NA_real_;mean_out <- if(no)mean(y[outside]) else NA_real_
        contrast_rows[[length(contrast_rows)+1L]] <- data.frame(contrast_id=id,branch_id=dd$branch_id[j],
          clone_id=cl,contrast_type=type,n_inside=ni,n_outside=no,status=state,
          mean_inside=mean_in,mean_outside=mean_out,
          difference=if(state=='ready')mean_in-mean_out else NA_real_,
          fraction_nonzero_inside=if(ni)mean(y[inside]!=0) else NA_real_,
          fraction_nonzero_outside=if(no)mean(y[outside]!=0) else NA_real_,
          contrast=paste0('Clone ',cl,': ',dd$branch_id[j],' / ',type),
          interpretation='descriptive within-clone cell contrast; not a validated ancestral shift')
        if(state=='ready')ready_membership[[id]] <- list(inside=inside,outside=outside)
      }
    }
  }
  contrasts <- if(length(contrast_rows))do.call(rbind,contrast_rows) else data.frame(
    contrast_id=character(),branch_id=character(),clone_id=character(),contrast_type=character(),
    n_inside=integer(),n_outside=integer(),status=character(),mean_inside=numeric(),mean_outside=numeric(),
    difference=numeric(),fraction_nonzero_inside=numeric(),fraction_nonzero_outside=numeric(),contrast=character(),interpretation=character())
  ready_contrasts <- contrasts[contrasts$status=='ready',]
  if(opt$model %in% c('all','contrasts')) {
    cat('\nWithin-clone contrast availability:\n');print(table(contrasts$status))
    cat('Sister comparisons use other children of the immediate cellular parent, restricted to the same clone.\n')
    cat('No sample/site adjustment or confirmatory p-values; small or unsplit groups have no estimated difference.\n')
    if(nrow(ready_contrasts))models$contrasts<-list(table=ready_contrasts)
    else cat('No within-clone contrasts meet support thresholds; no contrast resampling model will be available.\n')
  }
  coefficients <- if(length(coefficient_rows)) do.call(rbind,coefficient_rows) else
    data.frame(model=character(),term=character(),coefficient=numeric(),selected=logical(),interpretation=character())
  rownames(coefficients) <- NULL
  dd <- t$branch_dictionary[t$branch_dictionary$mode=='uncollapsed',]
  coefficients$descendant_clones <- dd$descendant_clones[match(coefficients$term,dd$branch_id)]
  for(field in c('n_descendant_cells','complete_clones_only','splits_within_clone'))coefficients[[field]]<-dd[[field]][match(coefficients$term,dd$branch_id)]
  chosen<-coefficients[coefficients$model=='elastic_net' & coefficients$term!='(Intercept)' & coefficients$selected,]
  if(nrow(chosen)){cat('\nSelected EN branch classes (complete clone partitions versus within-clone splits):\n');print(table(complete_clones_only=chosen$complete_clones_only));print(chosen[,c('term','coefficient','n_descendant_cells','descendant_clones','splits_within_clone')],row.names=FALSE)}
  metrics <- do.call(rbind,metrics)
  write_tsv(coefficients,'.coefficients.tsv');write_tsv(contrasts,'.contrasts.tsv')
  write_tsv(metrics,'.fit_metrics.tsv');write_tsv(predictions,'.cell_predictions.tsv.gz')
  saveRDS(list(gene=opt$gene,models=models,coefficients=coefficients,contrasts=ready_contrasts,all_contrasts=contrasts,
    fit_metrics=metrics,cell_predictions=predictions,expression_scale=e$expression_scale,tree_mode='uncollapsed',
    contrast_membership=ready_membership,contrast_type='within_clone',
    independent_branch_terms=basis,options=opt,tree_designs_md5=e$tree_designs_md5,
    expression_md5=unname(tools::md5sum(opt$expression)),
    inference_status='exploratory on uncertain within-clone topology; no historical localization established',
    sessionInfo=sessionInfo()),paste0(opt$output_prefix,'.models.rds'))
  # One compact diagnostic figure, all panels describe the same cells and response.
  grDevices::png(paste0(opt$output_prefix,'.comparison.png'),width=1800,height=750,res=150)
  tryCatch({
    graphics::par(mfrow=c(1,3),mar=c(5,4,4,1))
    graphics::boxplot(y~factor(clone,levels=levels_clone),xlab='Clone',ylab='Prepared expression units',
      main=paste(opt$gene,'observed expression'),outline=FALSE,col='grey90')
    graphics::points(seq_along(levels_clone),clone_means,pch=19,col='firebrick')
    fitted <- if('elastic_net' %in% names(predictions))predictions$elastic_net else if('linear' %in% names(predictions))predictions$linear else reference
    fitted_label <- if('elastic_net' %in% names(predictions))'Elastic Net' else if('linear' %in% names(predictions))'Linear' else 'Clone means'
    graphics::plot(fitted,y,pch=16,cex=.3,col=grDevices::adjustcolor('steelblue',alpha.f=.3),
      xlab=paste(fitted_label,'training fitted expression'),ylab='Observed expression',main='Training fit only')
    graphics::abline(a=0,b=1,lty=2)
    if(nrow(ready_contrasts)) {
      show<-utils::head(ready_contrasts[order(-abs(ready_contrasts$difference)),],10)
      graphics::par(mar=c(5,11,4,1))
      labels<-paste(show$clone_id,show$branch_id,ifelse(show$contrast_type=='within_clone','rest','sister'),sep=' / ')
      graphics::barplot(show$difference,names.arg=labels,horiz=TRUE,las=1,cex.names=.55,
        col=ifelse(show$difference<0,'steelblue','tan3'),main='Largest descriptive contrasts',xlab='Within-clone mean difference')
      graphics::abline(v=0)
    } else {graphics::plot.new();graphics::text(.5,.5,'No eligible within-clone contrasts')}
  },finally=grDevices::dev.off())
  # Adapted original selected-branch tree view; no within-clone support implied.
  edge_col<-rep('black',nrow(t$tree$edge));edge_width<-rep(1,nrow(t$tree$edge))
  selected<-coefficients[coefficients$model=='elastic_net' & coefficients$term!='(Intercept)' & coefficients$selected,]
  selected<-selected[!is.na(selected$term),]
  for(i in seq_len(nrow(selected))) {
    edge<-dd$edge_index[match(selected$term[i],dd$branch_id)]
    edge_col[edge]<-if(selected$coefficient[i]<0)'steelblue' else 'firebrick';edge_width[edge]<-3
  }
  grDevices::png(paste0(opt$output_prefix,'.selected_tree.png'),
    width=2000,height=2200,res=180,pointsize=14)
  tryCatch({
    palette<-grDevices::hcl.colors(length(levels_clone),'Dark 3')
    # Reserve a header panel so legends never overlap the cellular topology.
    header_rows<-ceiling(length(levels_clone)/4)
    graphics::layout(matrix(c(1,2),ncol=1),heights=c(2+header_rows,15))
    graphics::par(mar=c(0,1,1,1),fg='black',col.axis='black',col.lab='black')
    graphics::plot.new()
    graphics::text(.5,.88,paste(opt$gene,'exploratory uncollapsed selection'),
      cex=1.25,font=2)
    graphics::legend(x=.5,y=.70,xjust=.5,yjust=1,
      legend=paste('Clone',levels_clone),col=palette,lwd=7,
      ncol=min(4,length(levels_clone)),cex=1,bty='n')
    graphics::legend(x=.5,y=.18,xjust=.5,yjust=.5,
      legend=c('Cellular topology','Positive EN contribution','Negative EN contribution'),
      col=c('black','firebrick','steelblue'),lwd=c(1,3,3),ncol=3,cex=.85,bty='n')
    graphics::par(mar=c(4,1,0,2))
    ape::plot.phylo(t$tree,show.tip.label=FALSE,use.edge.length=FALSE,
      edge.color=edge_col,edge.width=edge_width,
      x.lim=c(0,max(ape::node.depth(t$tree))+5))
    # A thick strip beside the tips shows each cell's clone, including small clones.
    coords<-get('last_plot.phylo',envir=ape::.PlotPhyloEnv)
    tip_x<-max(coords$xx[seq_len(length(clone))])+1
    graphics::segments(tip_x,coords$yy[seq_along(clone)]-.5,
      tip_x,coords$yy[seq_along(clone)]+.5,
      col=palette[match(clone,levels_clone)],lwd=8,lend='butt')
    note<-if(!'elastic_net' %in% names(models))'Elastic Net was not fitted; topology only.' else if(!nrow(selected))'No Elastic Net branches selected.' else 'Color denotes simultaneous fitted EN contributions.'
    graphics::mtext(note,side=1,line=1,cex=.85,col='black')
    graphics::mtext('Exploratory cellular topology; within-clone ancestry remains uncertain',
      side=1,line=2.3,cex=.85,col='black')
  },finally=grDevices::dev.off())
  cat('\nTraining fit metrics:\n');print(metrics,row.names=FALSE)
  cat('\nLIMITATIONS\n',
    '- Training fit and tuning CV do not establish independent biological replication.\n',
    '- Cells may be dependent; no naive cell-level p-values or confidence intervals are supplied.\n',
    '- No sample/site covariates are fitted in this pilot; confounding is not resolved by balancing.\n',
    '- Response units are not automatically log2FC. Root/support limitations remain.\n',
    '- Use the shared Module 5 validator with this pilot for stability and bounded working-likelihood comparisons.\n',sep='')
  cat('\nOptions used:\n');print(opt)
  cat('\nFiles to return:\n');cat(paste(files,collapse='\n'),'\n')
  cat('\nUNCOLLAPSED FIT COMPLETE — run 05_validate_models.R on the saved pilot for validation.\n')
},error=function(e) {
  cat('\nUNCOLLAPSED PILOT FAILED:',conditionMessage(e),'\nReview partial outputs and terminal error.\n')
  exit_code <<- 1L
})
sink();close(log)
if(exit_code)message('Uncollapsed pilot failed; see ',opt$output_prefix,'.report.txt') else message('Report: ', opt$output_prefix,'.report.txt')
quit(status=exit_code)
