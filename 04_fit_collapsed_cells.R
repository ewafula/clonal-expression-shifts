#!/usr/bin/env Rscript
# Collapsed-tree cell-level pilot. Reuses original Gaussian cv.glmnet fitting.
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
  optparse::make_option('--output_prefix',default='results/module04')
)
opt <- optparse::parse_args(optparse::OptionParser(option_list=options_list))
# Match edges by descendant clone sets, never by cellular or collapsed node numbers.
plot_collapsed_selection <- function(t,coefficients,gene,has_en,path) {
  tree<-t$clone_tree
  if(is.null(tree))stop('Saved Module 2 clone tree missing.')
  dd<-t$branch_dictionary[t$branch_dictionary$mode=='collapsed',]
  cc<-coefficients[coefficients$model=='elastic_net' & coefficients$term!='(Intercept)' & coefficients$selected,]
  children<-split(tree$edge[,2],tree$edge[,1]);nt<-length(tree$tip.label)
  below<-function(node){if(node<=nt)return(tree$tip.label[node]);unlist(lapply(children[[as.character(node)]],below),use.names=FALSE)}
  key<-function(x)paste(sort(x),collapse=';')
  keys<-vapply(tree$edge[,2],function(n)key(below(n)),'')
  edge<-match(vapply(strsplit(dd$descendant_clones[match(cc$term,dd$branch_id)],';',fixed=TRUE),key,''),keys)
  if(anyNA(edge)||anyDuplicated(edge))stop('Selected branches cannot be mapped uniquely to the collapsed tree.')
  colors<-rep('black',nrow(tree$edge));widths<-rep(1.5,nrow(tree$edge))
  colors[edge]<-ifelse(cc$coefficient<0,'steelblue','firebrick');widths[edge]<-4
  sizes<-table(as.character(t$cell_metadata$clone_id))
  labels<-paste0('Clone ',tree$tip.label,' (n=',sizes[tree$tip.label],')')
  grDevices::png(path,width=2000,height=max(1100,100*nt),res=180,pointsize=14)
  tryCatch({
    graphics::layout(matrix(c(1,2),ncol=1),heights=c(2,10))
    graphics::par(mar=c(0,1,1,1));graphics::plot.new()
    graphics::text(.5,.8,paste(gene,'collapsed selected branches'),font=2,cex=1.3)
    graphics::legend(.5,.35,xjust=.5,legend=c('Unselected branch','Positive EN contribution','Negative EN contribution'),col=c('black','firebrick','steelblue'),lwd=c(1.5,4,4),ncol=3,bty='n',cex=.85)
    graphics::par(mar=c(4,1,0,1))
    tree$tip.label<-labels
    ape::plot.phylo(tree,use.edge.length=FALSE,show.tip.label=TRUE,edge.color=colors,edge.width=widths,cex=1,font=1,
      x.lim=c(0,max(ape::node.depth(tree))+2))
    if(length(edge)&&length(edge)<=12)ape::edgelabels(paste0(cc$term,'\n',formatC(cc$coefficient,format='f',digits=3)),edge=edge,frame='none',cex=.8,col=colors[edge],adj=c(.5,-.3))
    note<-if(!has_en)'Elastic Net was not fitted; topology only.' else if(!length(edge))'No Elastic Net branches selected.' else if(length(edge)>12)'Color denotes fitted EN contribution; branch labels omitted for legibility (see coefficients table).' else 'Color denotes fitted EN contribution; labels give branch ID and coefficient.'
    graphics::mtext(note,side=1,line=1,cex=.8)
    graphics::mtext('Topology only; response units are not automatically log2FC. Selected branches are not unique historical shifts.',side=1,line=2.3,cex=.75)
  },finally=grDevices::dev.off())
}
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
  cat('MODULE 4 — COLLAPSED CELL-LEVEL PILOT\n',format(Sys.time(),tz='UTC',usetz=TRUE),'\n\n')
  e <- readRDS(opt$expression);t <- readRDS(opt$tree_designs)
  if(is.null(t$designs$collapsed)) stop('Module 2 file has no collapsed design.')
  if(is.null(e$tree_designs_md5) || !identical(e$tree_designs_md5,unname(tools::md5sum(opt$tree_designs))))
    stop('Module 3 expression was prepared against a different Module 2 file; use matching inputs.')
  if(!opt$gene %in% rownames(e$expression)) stop('Gene not prepared: ',opt$gene,'. Prepare it with Module 3 --genes first.')
  x <- t$designs$collapsed$cell_X;cells <- e$cell_metadata
  if(!identical(colnames(e$expression),rownames(x)) || !identical(cells$cell_id,rownames(x)) ||
     !identical(cells$clone_id,t$cell_metadata$clone_id)) stop('Expression/design cell order or clone labels disagree.')
  clone <- as.character(cells$clone_id)
  if(anyNA(clone) || any(!nzchar(clone)) || anyDuplicated(cells$cell_id)) stop('Invalid cell/clone IDs.')
  y <- as.numeric(e$expression[opt$gene,,drop=TRUE]);n <- length(y)
  if(any(!is.finite(y)) || n<3L || stats::sd(y)==0) stop('Gene response must be finite, variable, with >=3 cells.')
  if(ncol(x)<1L) stop('No collapsed branch predictors.')
  if(any(!x@x %in% c(0,1)) || any(Matrix::colSums(x)<=0 | Matrix::colSums(x)>=n)) stop('Collapsed design must contain varying binary predictors.')
  levels_clone <- sort(unique(clone));sizes <- table(factor(clone,levels=levels_clone))
  if(length(levels_clone)<2L) stop('At least two clones required.')
  clone_x <- x[match(levels_clone,clone),,drop=FALSE]
  if(!all(x==clone_x[match(clone,levels_clone),,drop=FALSE])) stop('Collapsed predictors split cells within a clone.')
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
  basis <- t$designs$collapsed$independent_branch_terms
  if(is.null(basis) || !all(basis %in% colnames(x))) stop('Saved independent linear basis is missing/invalid.')
  q_all <- qr(cbind(1,as.matrix(clone_x)))
  q_basis <- qr(cbind(1,as.matrix(clone_x[,basis,drop=FALSE])))
  if(q_basis$rank!=length(basis)+1L || q_basis$rank!=q_all$rank) stop('Saved basis does not span collapsed design.')
  if(q_basis$rank!=length(levels_clone)) stop('Collapsed design does not span all clone means; review tree design.')
  if(n<=q_basis$rank) stop('No residual degrees of freedom for a cell-level pilot.')
  cat('Gene:',opt$gene,'| cells:',n,'| clones:',length(levels_clone),'\n')
  cat('Response scale:',e$expression_scale,'\n')
  cat('Clone-balanced weights: inverse clone cell count, rescaled to mean 1.\n')
  cat('Each clone contributes equal total weight. These are balancing, not known precision weights.\n')
  cat('Retained branch predictors:',ncol(x),'| independent directions:',q_basis$rank-1L,'\n')
  if(opt$model %in% c('all','elastic_net')) {
    if(ncol(x)<2L) stop('glmnet requires >=2 predictor columns. For this tiny design choose --model linear or contrasts.')
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
    if(!isTRUE(all.equal(pred,reference,tolerance=1e-8))) stop('Linear basis predictions disagree with direct clone means.')
    coefficient_rows$linear <- data.frame(model='linear',term=names(stats::coef(linear_fit)),
      coefficient=as.numeric(stats::coef(linear_fit)),selected=NA,
      interpretation='unpenalized coefficient in saved full-rank basis; not variable selection')
    predictions$linear <- pred;predictions$residual_linear <- y-pred
    metrics[[length(metrics)+1L]] <- metric('linear',pred);models$linear <- linear_fit
    cat('\nLinear basis:',paste(basis,collapse=', '),'\n')
    cat('Basis predictions equal the full clone-mean reference. Basis choice does not identify unique historical shifts.\n')
  }
  contrasts <- data.frame(branch_id=character(),contrast=character(),descendant_clones=character(),
    comparison_clones=character(),inside_equal_clone_mean=numeric(),outside_equal_clone_mean=numeric(),
    difference=numeric(),clone_weights=character(),equivalent_to=character(),
    interpretation=character())
  if(opt$model %in% c('all','contrasts')) {
    dd <- t$branch_dictionary
    dd <- dd[dd$mode=='collapsed' & dd$complete_clones_only & dd$n_descendant_cells<n,]
    rows <- lapply(seq_len(nrow(dd)),function(i) {
      inside <- sort(unique(clone[t$descendant_tip_indices[[dd$edge_index[i]]]]))
      outside <- setdiff(levels_clone,inside)
      if(!length(inside) || !length(outside)) stop('Invalid contrast partition.')
      cw <- setNames(rep(0,length(levels_clone)),levels_clone)
      cw[inside] <- 1/length(inside);cw[outside] <- -1/length(outside)
      data.frame(branch_id=dd$branch_id[i],
        contrast=paste0('{',paste(inside,collapse=','),'} vs {',paste(outside,collapse=','),'}'),
        descendant_clones=paste(inside,collapse=';'),comparison_clones=paste(outside,collapse=';'),
        inside_equal_clone_mean=mean(clone_means[inside]),outside_equal_clone_mean=mean(clone_means[outside]),
        difference=sum(cw*clone_means[names(cw)]),
        clone_weights=paste(names(cw),format(cw,digits=8),sep=':',collapse=';'),
        equivalent_to=dd$equivalent_to[i],
        interpretation='descriptive equal-clone mean contrast; no confirmatory p-value')
    })
    contrasts <- do.call(rbind,rows)
    cat('\nExplicit contrasts (inside minus outside; equal weight per clone on each side):\n')
    print(contrasts[,c('contrast','difference')],row.names=FALSE)
    cat('Complementary root contrasts are reported in both biological orientations and are redundant sign reversals.\n')
    cat('Contrasts differ from simultaneous Elastic Net coefficients and linear basis coefficients.\n')
    models$contrasts <- list(clone_means=clone_means,table=contrasts)
  }
  coefficients <- if(length(coefficient_rows)) do.call(rbind,coefficient_rows) else
    data.frame(model=character(),term=character(),coefficient=numeric(),selected=logical(),interpretation=character())
  rownames(coefficients) <- NULL
  dd <- t$branch_dictionary[t$branch_dictionary$mode=='collapsed',]
  coefficients$descendant_clones <- dd$descendant_clones[match(coefficients$term,dd$branch_id)]
  metrics <- do.call(rbind,metrics)
  write_tsv(coefficients,'.coefficients.tsv');write_tsv(contrasts,'.contrasts.tsv')
  write_tsv(metrics,'.fit_metrics.tsv');write_tsv(predictions,'.cell_predictions.tsv.gz')
  saveRDS(list(gene=opt$gene,models=models,coefficients=coefficients,contrasts=contrasts,
    fit_metrics=metrics,cell_predictions=predictions,expression_scale=e$expression_scale,
    tree_mode='collapsed',contrast_type='equal_clone',
    independent_branch_terms=basis,options=opt,tree_designs_md5=e$tree_designs_md5,
    expression_md5=unname(tools::md5sum(opt$expression)),
    inference_status='descriptive pilot; cell dependence and selection uncertainty not resolved',
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
    if(nrow(contrasts)) {
      show<-utils::head(contrasts[order(-abs(contrasts$difference)),],10)
      graphics::par(mar=c(5,11,4,1))
      graphics::barplot(show$difference,names.arg=show$contrast,horiz=TRUE,las=1,
        col=ifelse(show$difference<0,'steelblue','tan3'),cex.names=.65,
        xlab='Equal-clone mean difference',main='Largest descriptive contrasts')
      graphics::abline(v=0)
    } else {
      graphics::plot.new();graphics::text(.5,.5,'Contrasts not requested')
    }
  },finally=grDevices::dev.off())
  plot_collapsed_selection(t,coefficients,opt$gene,'elastic_net' %in% names(models),paste0(opt$output_prefix,'.selected_tree.png'))
  cat('\nTraining fit metrics:\n');print(metrics,row.names=FALSE)
  cat('\nLIMITATIONS\n',
    '- Training fit and tuning CV do not establish independent biological replication.\n',
    '- Cells may be dependent; no naive cell-level p-values or confidence intervals are supplied.\n',
    '- No sample/site covariates are fitted in this pilot; confounding is not resolved by balancing.\n',
    '- Response units are not automatically log2FC. Root/support limitations remain.\n',
    '- Stability assessment and compatible-likelihood AIC/BIC comparisons belong to Module 5.\n',sep='')
  cat('\nOptions used:\n');print(opt)
  cat('\nFiles to return:\n');cat(paste(files,collapse='\n'),'\n')
  cat('\nCOLLAPSED FIT COMPLETE — run 05_validate_models.R on the saved pilot for validation.\n')
},error=function(e) {
  cat('\nCOLLAPSED PILOT FAILED:',conditionMessage(e),'\nReview partial outputs and terminal error.\n')
  exit_code <<- 1L
})
sink();close(log)
if(exit_code)message('Collapsed Pilot failed; see ',opt$output_prefix,'.report.txt') else message('Report: ', opt$output_prefix,'.report.txt')
quit(status=exit_code)
