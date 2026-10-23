#!/usr/bin/env Rscript
# Clone-stratified subsampling and explicitly declared working-likelihood IC.
# No naive p-values, bootstrap confidence intervals, or independent-patient claims.

for(pkg in c('optparse','Matrix','glmnet')) if(!requireNamespace(pkg,quietly=TRUE))
  stop('Missing package ',pkg,'. Run Rscript 00_install_packages.R in this R environment.')
options_list <- list(
  optparse::make_option('--pilot',default='results/module04.models.rds'),
  optparse::make_option('--expression',default='results/module03.expression.rds'),
  optparse::make_option('--tree_designs',default='results/module02.tree_designs.rds'),
  optparse::make_option('--reps',type='integer',default=100L,help='Number of subsamples [default %default]'),
  optparse::make_option('--fraction',type='double',default=0.8,help='Fraction sampled without replacement within each clone [default %default]'),
  optparse::make_option('--seed',type='integer',default=12345L),
  optparse::make_option('--skip_ic',action='store_true',default=FALSE,help='Skip working iid Gaussian ML comparisons'),
  optparse::make_option('--output_prefix',default='results/module04_validation',help='Output prefix; explicitly match the supplied pilot [default %default]')
)
opt <- optparse::parse_args(optparse::OptionParser(option_list=options_list))
if(is.na(opt$reps) || opt$reps<1L || !is.finite(opt$fraction) || opt$fraction<=0 || opt$fraction>=1 ||
   is.na(opt$seed) || opt$seed<0 || as.double(opt$seed)+opt$reps>.Machine$integer.max)
  stop('Use reps>=1, fraction strictly between 0 and 1, and a valid nonnegative seed.')
for(p in c(opt$pilot,opt$expression,opt$tree_designs))if(!file.exists(p))stop('Missing input: ',p)
outputs <- paste0(opt$output_prefix,c('.report.txt','.stability.tsv','.resamples.tsv.gz',
  '.iteration_status.tsv','.model_comparison.tsv','.validation.rds','.stability.png'))
if(any(file.exists(outputs)))stop('Output exists; choose a new --output_prefix. No overwrite performed.')
if(!dir.exists(dirname(opt$output_prefix)))dir.create(dirname(opt$output_prefix),recursive=TRUE)
write_tsv <- function(d,suffix) {
 p <- paste0(opt$output_prefix,suffix);con <- if(grepl('\\.gz$',p))gzfile(p,'wt') else file(p,'wt')
 tryCatch(utils::write.table(d,con,sep='\t',quote=TRUE,row.names=FALSE,na='NA'),finally=close(con))
}
# Return NA rather than inventing a coefficient when a resample cannot be fit.
not_estimable <- function(message) structure(list(message=message),class=c('not_estimable','error','condition'))
log <- file(paste0(opt$output_prefix,'.report.txt'),'wt');sink(log,split=TRUE)
exit_code <- 0L
tryCatch({
  cat('SHARED MODEL VALIDATOR (MODULES 5–6) — STABILITY AND WORKING-LIKELIHOOD COMPARISONS\n',format(Sys.time(),tz='UTC',usetz=TRUE),'\n\n')
  a <- readRDS(opt$pilot);e <- readRDS(opt$expression);t <- readRDS(opt$tree_designs)
  if(!identical(a$expression_md5,unname(tools::md5sum(opt$expression))) ||
     !identical(a$tree_designs_md5,unname(tools::md5sum(opt$tree_designs))) ||
     !identical(e$tree_designs_md5,a$tree_designs_md5))stop('Pilot/expression/tree files do not match; use the original compatible inputs.')
  tree_mode<-a$tree_mode
  if(length(tree_mode)!=1L || !tree_mode %in% c('collapsed','uncollapsed'))stop('Pilot must declare collapsed or uncollapsed tree_mode; use a current Module 4 or 6 fit.')
  x <- t$designs[[tree_mode]]$cell_X;cells <- e$cell_metadata;gene <- a$gene
  if(is.null(x) || !gene %in% rownames(e$expression))stop('Requested design/prepared pilot gene absent.')
  if(!identical(cells$cell_id,rownames(x)) || !identical(colnames(e$expression),rownames(x)) ||
     !identical(cells$clone_id,t$cell_metadata$clone_id) || !identical(a$cell_predictions$cell_id,cells$cell_id))stop('Cell order/clone alignment mismatch.')
  y <- as.numeric(e$expression[gene,]);clone <- as.character(cells$clone_id);n <- length(y)
  if(any(!is.finite(y)) || stats::sd(y)==0)stop('Pilot response is nonfinite/constant.')
  if(!isTRUE(all.equal(y,a$cell_predictions$observed,tolerance=0)))stop('Pilot response does not match prepared expression.')
  levels_clone <- sort(unique(clone));groups <- split(seq_len(n),clone)
  basis <- a$independent_branch_terms
  if(!all(basis %in% colnames(x)))stop('Saved linear basis absent from design.')
  engines <- intersect(c('elastic_net','linear','contrasts'),names(a$models))
  if(!length(engines))stop('No recognized fitted methods in pilot.')
  definitions <- list()
  for(method in engines) {
    if(method=='contrasts') {
      definitions[[method]] <- data.frame(model=method,term=if(is.null(a$contrast_membership))a$contrasts$branch_id else a$contrasts$contrast_id,full_value=a$contrasts$difference,
        effect_type=if(is.null(a$contrast_membership))'equal-clone contrast' else 'within-clone cell contrast',
        label=if(is.null(a$contrast_membership))a$contrasts$contrast else paste(a$contrasts$clone_id,a$contrasts$branch_id,ifelse(a$contrasts$contrast_type=='within_clone','rest','sister'),sep=' / '))
    } else {
      cc <- a$coefficients[a$coefficients$model==method & a$coefficients$term!='(Intercept)',]
      definitions[[method]] <- data.frame(model=method,term=cc$term,full_value=cc$coefficient,
        effect_type=if(method=='elastic_net')'penalized branch coefficient' else 'linear basis coefficient',label=cc$term)
    }
  }
  definition <- do.call(rbind,definitions);rownames(definition) <- NULL
  if(!is.null(a$contrast_membership) && 'contrasts' %in% engines && !identical(names(a$contrast_membership),definitions$contrasts$term))stop('Within-clone membership order mismatch.')
  if('elastic_net' %in% engines && !identical(definitions$elastic_net$term,colnames(x)))stop('Pilot branch coefficient order differs from design.')
  # Derive contrast membership from saved edge descendants, not delimited clone names.
  contrast_groups <- if('contrasts' %in% engines && is.null(a$contrast_membership))lapply(a$contrasts$branch_id,function(id) {
    dd <- t$branch_dictionary[t$branch_dictionary$mode=='collapsed' & t$branch_dictionary$branch_id==id,]
    if(nrow(dd)!=1L)stop('Contrast branch missing from dictionary.')
    inside <- sort(unique(clone[t$descendant_tip_indices[[dd$edge_index]]]))
    list(inside=inside,outside=setdiff(levels_clone,inside))
  }) else list()
  cat('Tree mode:',tree_mode,'\n')
  cat('Gene:',gene,'| methods:',paste(engines,collapse=', '),'| cells:',n,'\n')
  cat('Subsamples:',opt$reps,'| fraction per clone:',opt$fraction,'| sampling WITHOUT replacement\n')
  cat('Elastic Net is retuned inside every resample using the pilot alpha/lambda rule.\n')
  cat('Weights recomputed by retained clone sizes; branch identities fixed to saved Module 2 descendants.\n')
  cat('No sample/site stratification or between-patient replication is implied.\n\n')
  values <- list();status_rows <- list();samples <- vector('list',opt$reps)
  for(b in seq_len(opt$reps)) {
    set.seed(opt$seed+b)
    ii <- sort(unlist(lapply(groups,function(g)g[sample.int(length(g),max(1L,floor(opt$fraction*length(g))))]),use.names=FALSE))
    yy <- y[ii];xx <- x[ii,,drop=FALSE];cl <- clone[ii];sz <- table(factor(cl,levels=levels_clone))
    w <- 1/as.numeric(sz[match(cl,levels_clone)]);w <- w/mean(w)
    samples[[b]] <- list(seed=opt$seed+b,cell_indices=ii,clone_sizes=sz,fold_id=NULL)
    for(method in engines) {
      warnings <- character();lambda <- NA_real_;k <- NA_integer_
      result <- tryCatch(withCallingHandlers({
        if(method=='elastic_net') {
          k <- min(a$options$nfolds,min(sz))
          if(k<3L || ncol(xx)<2L || stats::sd(yy)==0)stop(not_estimable('Insufficient cells/folds/predictors or constant resampled response.'))
          fold <- integer(length(ii))
          for(c in levels_clone) {
            jj <- which(cl==c);fold[jj] <- sample(rep(seq_len(k),length.out=length(jj)))
          }
          samples[[b]]$fold_id <- fold
          fit <- glmnet::cv.glmnet(x=xx,y=yy,weights=w,foldid=fold,nfolds=k,
            family='gaussian',alpha=a$options$alpha,standardize=TRUE,intercept=TRUE)
          if(fit$glmnet.fit$jerr!=0)stop('glmnet reported nonzero convergence status.')
          lambda <- fit[[a$options$lambda_choice]]
          z <- as.matrix(stats::coef(fit,s=a$options$lambda_choice))
          val <- z[definitions[[method]]$term,1]
        } else if(method=='linear') {
          z <- cbind(1,as.matrix(xx[,basis,drop=FALSE]))
          if(nrow(z)<=ncol(z) || qr(z)$rank<ncol(z))stop(not_estimable('Linear basis has no residual degrees of freedom or is not estimable.'))
          fit <- stats::lm.wfit(z,yy,w=w)
          val <- fit$coefficients[-1L]
        } else if(!is.null(a$contrast_membership)) {
          val<-vapply(a$contrast_membership,function(g) {
            inside<-which(ii %in% g$inside);outside<-which(ii %in% g$outside)
            if(min(length(inside),length(outside))<a$options$min_contrast_cells)return(NA_real_)
            mean(yy[inside])-mean(yy[outside])
          },0.0)
          if(all(is.na(val)))stop(not_estimable('No within-clone contrast meets resampled support thresholds.'))
        } else {
          means <- vapply(levels_clone,function(c)mean(yy[cl==c]),0.0);names(means)<-levels_clone
          val <- vapply(contrast_groups,function(g)mean(means[g$inside])-mean(means[g$outside]),0.0)
        }
        if(any(!is.finite(val)) && !(method=='contrasts' && !is.null(a$contrast_membership) && all(is.finite(val)|is.na(val))))stop('Nonfinite resample estimates.')
        list(value=as.numeric(val))
      },warning=function(z){warnings <<- c(warnings,conditionMessage(z));invokeRestart('muffleWarning')}),
      error=function(z)z)
      ok <- !inherits(result,'error')
      state <- if(ok)'success' else if(inherits(result,'not_estimable'))'not_estimable' else 'failed'
      status_rows[[length(status_rows)+1L]] <- data.frame(iteration=b,seed=opt$seed+b,model=method,
        n_cells=length(ii),status=state,lambda=lambda,nfolds=k,
        message=if(ok)'' else conditionMessage(result),warnings=paste(unique(warnings),collapse=' | '))
      if(ok) {
        good<-is.finite(result$value)
        values[[length(values)+1L]] <- data.frame(iteration=b,model=method,term=definitions[[method]]$term[good],value=result$value[good])
      }
    }
    if(b==1L || b%%10L==0L || b==opt$reps)cat('Completed subsample',b,'of',opt$reps,'\n')
  }
  statuses <- do.call(rbind,status_rows)
  draws <- if(length(values))do.call(rbind,values) else data.frame(iteration=integer(),model=character(),term=character(),value=numeric())
  summary_rows <- lapply(seq_len(nrow(definition)),function(j) {
    d<-definition[j,];st<-statuses[statuses$model==d$model,]
    v<-draws$value[draws$model==d$model & draws$term==d$term]
    nv<-length(v);nsel<-sum(v!=0);np<-sum(v>0);nn<-sum(v<0)
    data.frame(d,attempted=opt$reps,successful=nv,failed=sum(st$status=='failed'),not_estimable=sum(st$status=='not_estimable')+sum(st$status=='success')-nv,
      full_selected=if(d$model=='elastic_net')d$full_value!=0 else NA,
      selection_frequency=if(d$model=='elastic_net' && nv>0)nsel/nv else NA_real_,
      positive_frequency=if(nv)np/nv else NA_real_,negative_frequency=if(nv)nn/nv else NA_real_,
      sign_consistency_among_nonzero=if(nsel)max(np,nn)/nsel else NA_real_,
      agreement_with_full_sign_among_nonzero=if(nsel && d$full_value!=0)sum(sign(v[v!=0])==sign(d$full_value))/nsel else NA_real_,
      median=if(nv)stats::median(v) else NA_real_,
      q05=if(nv)unname(stats::quantile(v,.05)) else NA_real_,q95=if(nv)unname(stats::quantile(v,.95)) else NA_real_)
  })
  stability <- do.call(rbind,summary_rows)
  stability$overlap_neighborhood_frequency<-NA_real_;stability$n_overlap_neighbors<-NA_integer_
  if(tree_mode=='uncollapsed' && 'elastic_net' %in% engines) {
    inter<-as.matrix(Matrix::crossprod(x));sizes<-as.numeric(Matrix::colSums(x))
    similarity<-inter/(outer(sizes,sizes,'+')-inter)
    for(i in which(stability$model=='elastic_net')) {
      j<-match(stability$term[i],colnames(x));neighbors<-colnames(x)[similarity[j,]>=.8]
      valid_iterations<-statuses$iteration[statuses$model=='elastic_net' & statuses$status=='success']
      stability$n_overlap_neighbors[i]<-length(neighbors)
      if(length(valid_iterations))stability$overlap_neighborhood_frequency[i]<-mean(vapply(valid_iterations,function(b)
        any(draws$value[draws$iteration==b & draws$model=='elastic_net' & draws$term %in% neighbors]!=0),FALSE))
    }
    cat('Overlap neighborhood uses Jaccard >=0.8 on fixed descendant cell sets, including the branch itself. It is a sensitivity summary, not a distinct historical region.\n')
  }
  write_tsv(stability,'.stability.tsv');write_tsv(draws,'.resamples.tsv.gz');write_tsv(statuses,'.iteration_status.tsv')
  cat('\nIteration outcomes:\n');print(table(statuses$model,statuses$status))
  cat('Frequencies use successful estimates only; failed/non-estimable fits remain explicit and never become zero coefficients.\n')
  cat('Non-Elastic-Net methods are not variable selection: their selection_frequency is NA.\n')
  cat('q05/q95 are subsampling ranges including zero estimates, NOT confidence intervals or p-values.\n')

  # Separate ordinary unweighted ML refits: balancing weights are NOT precisions.
  ic <- data.frame(model_id=character(),source=character(),terms=character(),n=integer(),rank=integer(),k=integer(),
    residual_df=integer(),logLik=numeric(),AIC=numeric(),BIC=numeric(),AICc=numeric(),
    delta_AIC=numeric(),delta_BIC=numeric(),equivalence_group=character(),status=character())
  if(!opt$skip_ic) {
    cat('\nINFORMATION CRITERIA: ordinary unweighted iid Gaussian ML refits on the SAME gene/cells.\n')
    cat('This working likelihood assumes independent homoskedastic cell errors, which are not established.\n')
    cat('It targets ordinary cell-weighted fit, not the clone-balanced objective. No penalized glmnet objective is scored.\n')
    p<-ncol(x);sets<-list(intercept_only=character());sources<-'baseline'
    if(p<=10L) {
      for(s in seq_len(p))for(v in utils::combn(colnames(x),s,simplify=FALSE)) {
        sets[[paste0('branch_subset_',length(sets))]]<-v;sources<-c(sources,'prespecified_candidate_dictionary_search')
      }
      cat('Enumerating',2^p,'branch subsets including the intercept-only model; dependent subsets are excluded.\n')
    } else {
      for(v in colnames(x)){sets[[paste0('single_',length(sets))]]<-v;sources<-c(sources,'single_branch')}
      cat('More than 10 predictors: bounded comparison of null, single branches, pilot active set and linear basis; no exhaustive search.\n')
    }
    if('elastic_net' %in% engines) {
      cc<-a$coefficients
      sets$pilot_active_set<-cc$term[cc$model=='elastic_net' & cc$term!='(Intercept)' & cc$selected]
      sources<-c(sources,'post_selection_unpenalized_active_set_refit')
    }
    if(p>10L && 'elastic_net' %in% engines && length(sets$pilot_active_set)) {
      active<-sets$pilot_active_set
      counts<-as.numeric(Matrix::colSums(x));intersection<-as.matrix(Matrix::crossprod(x))
      for(term in active) {
        sets[[paste0('active_minus_',term)]]<-setdiff(active,term);sources<-c(sources,'post_selection_drop_one_refit')
        j<-match(term,colnames(x));anc<-which(counts>counts[j] & intersection[j,]==counts[j])
        if(length(anc)) {
          parent<-anc[which.min(counts[anc])]
          sets[[paste0('ancestor_child_',term)]]<-c(colnames(x)[parent],term)
          sources<-c(sources,'selected_child_with_nearest_candidate_ancestor_refit')
        }
      }
    }
    sets$linear_basis<-basis;sources<-c(sources,'unpenalized_full_rank_basis_refit')
    matrices<-lapply(sets,function(v)cbind('(Intercept)'=1,as.matrix(x[,v,drop=FALSE])))
    matrices$clone_means_reference<-stats::model.matrix(~factor(clone))
    sets$clone_means_reference<-character();sources<-c(sources,'unpenalized_clone_means_refit')
    # Distinct joint clone/branch row patterns preserve each candidate prediction space.
    space_rows<-which(!duplicated(data.frame(clone=clone,as.matrix(x),check.names=FALSE)))
    space_keys<-character()
    rows<-lapply(seq_along(matrices),function(i) {
      z<-matrices[[i]];q<-qr(z);r<-q$rank;k<-r+1L
      row<-data.frame(model_id=names(matrices)[i],source=sources[i],
        terms=if(names(matrices)[i]=='clone_means_reference')'full clone-mean parameterization' else paste(sets[[i]],collapse=';'),
        n=n,rank=r,k=k,residual_df=n-r,logLik=NA_real_,AIC=NA_real_,BIC=NA_real_,AICc=NA_real_,
        delta_AIC=NA_real_,delta_BIC=NA_real_,equivalence_group=NA_character_,status='working_iid_cell_gaussian_ML')
      if(r<ncol(z)){row$status<-'rank_deficient';return(row)}
      if(n<=r){row$status<-'no_residual_degrees_of_freedom';return(row)}
      fit<-stats::lm.fit(z,y);rss<-sum(fit$residuals^2)
      if(!is.finite(rss) || rss<=.Machine$double.eps*max(1,sum(y^2))) {
        row$status<-'zero_or_numerically_zero_residual_variance';return(row)
      }
      row$logLik<- -n/2*(log(2*pi)+1+log(rss/n))
      row$AIC<- -2*row$logLik+2*k;row$BIC<- -2*row$logLik+k*log(n)
      row$AICc<-if(n>k+1L)row$AIC+2*k*(k+1)/(n-k-1) else NA_real_
      # Same prediction space can have different branch explanations. Group them.
      cq<-qr(z[space_rows,,drop=FALSE]);qq<-qr.Q(cq)[,seq_len(r),drop=FALSE]
      projection<-tcrossprod(qq);projection[abs(projection)<1e-10]<-0
      key<-paste(format(round(projection,10),scientific=FALSE,trim=TRUE),collapse='|')
      group<-match(key,space_keys)
      if(is.na(group)){space_keys<<-c(space_keys,key);group<-length(space_keys)}
      row$equivalence_group<-paste0('space_',group);row
    })
    ic<-do.call(rbind,rows)
    valid<-is.finite(ic$AIC)&is.finite(ic$BIC)
    if(any(valid)) {
      ic$delta_AIC[valid]<-ic$AIC[valid]-min(ic$AIC[valid]);ic$delta_BIC[valid]<-ic$BIC[valid]-min(ic$BIC[valid])
      best<-ic[valid,];best<-best[order(best$BIC),];best<-best[!duplicated(best$equivalence_group),]
      cat('\nLowest working BIC scores by distinct prediction space (up to five):\n')
      print(utils::head(best[,c('model_id','equivalence_group','rank','AIC','BIC','delta_BIC')],5),row.names=FALSE)
    } else cat('No valid likelihood comparisons; see explicit status column.\n')
    cat('k includes design rank (intercept included) plus ONE residual variance parameter. n is the cell count under this declared likelihood.\n')
    cat('Equivalent spaces do not identify unique historical branches. Search/post-selection and cell dependence are not corrected by IC.\n')
  } else cat('\nIC skipped by request.\n')
  write_tsv(ic,'.model_comparison.tsv')
  saveRDS(list(gene=gene,tree_mode=tree_mode,stability=stability,resamples=draws,iteration_status=statuses,
    samples=samples,model_comparison=ic,options=opt,pilot_options=a$options,
    response_scale=e$expression_scale,comparison_likelihood=if(opt$skip_ic)'skipped' else 'ordinary unweighted iid cell Gaussian ML; exploratory working assumption',
    pilot_md5=unname(tools::md5sum(opt$pilot)),expression_md5=a$expression_md5,tree_designs_md5=a$tree_designs_md5,
    sessionInfo=sessionInfo()),paste0(opt$output_prefix,'.validation.rds'))
  grDevices::png(paste0(opt$output_prefix,'.stability.png'),width=1700,height=800,res=150)
  tryCatch({
    graphics::par(mfrow=c(1,2),mar=c(7,5,4,1))
    en<-stability[stability$model=='elastic_net',]
    if(nrow(en)>15L)en<-utils::head(en[order(-en$selection_frequency,-abs(en$full_value),na.last=TRUE),],15)
    if(nrow(en) && any(is.finite(en$selection_frequency))) {
      graphics::barplot(en$selection_frequency,names.arg=en$term,las=2,ylim=c(0,1),col='steelblue',
        ylab='Selection fraction among successful fits',main=paste(gene,'EN stability (up to 15 branches)'),cex.names=.7)
    } else {graphics::plot.new();graphics::text(.5,.5,'No usable Elastic Net selection frequencies')}
    contrasts<-stability[stability$model=='contrasts' & is.finite(stability$median),]
    if(!nrow(contrasts))contrasts<-stability[stability$model=='linear' & is.finite(stability$median),]
    if(nrow(contrasts)>12L)contrasts<-utils::head(contrasts[order(-abs(contrasts$full_value)),],12)
    if(nrow(contrasts)) {
      graphics::par(mar=c(7,9,4,1))
      graphics::plot(contrasts$median,seq_len(nrow(contrasts)),type='n',xlim=range(c(contrasts$q05,contrasts$q95)),
        ylim=c(.5,nrow(contrasts)+.5),yaxt='n',pch=19,xlab='Prepared response units',ylab='',
        main='Largest effects: median and 5–95% range')
      graphics::axis(2,at=seq_len(nrow(contrasts)),labels=contrasts$label,las=1,cex.axis=.7)
      graphics::segments(contrasts$q05,seq_len(nrow(contrasts)),contrasts$q95,seq_len(nrow(contrasts)),lwd=2,col='steelblue')
      for(q in c('q05','q95'))graphics::segments(contrasts[[q]],seq_len(nrow(contrasts))-.08,contrasts[[q]],seq_len(nrow(contrasts))+.08,col='steelblue')
      graphics::points(contrasts$median,seq_len(nrow(contrasts)),pch=19,cex=.6)
      graphics::abline(v=0,lty=2);graphics::mtext('Descriptive ranges, not confidence intervals',side=1,line=4,cex=.8)
    } else {graphics::plot.new();graphics::text(.5,.5,'No usable contrast/linear resampling summaries')}
  },finally=grDevices::dev.off())
  cat('\nLIMITATIONS: cell-sampling/tuning stability on one fixed tree; no new samples/patients, topology support, or confounding adjustment.\n')
  cat('Selection frequencies are not p-values; no automatic strong-support cutoff is imposed.\n')
  cat('Information criteria cannot be compared across genes, datasets, transformations, or balancing/precision assumptions.\n')
  cat('\nOptions used:\n');print(opt)
  cat('\nFiles to return:\n');cat(paste(outputs,collapse='\n'),'\n')
  cat('\nVALIDATION COMPLETE — review method success counts and working-model assumptions before running gene panel Module 7.\n')
},error=function(z){cat('\nVALIDATION FAILED:',conditionMessage(z),'\nReview partial outputs and terminal error.\n');exit_code<<-1L})
sink();close(log)
if(exit_code)message('Validation failed; see ',opt$output_prefix,'.report.txt') else message('Report: ', opt$output_prefix,'.report.txt')
quit(status=exit_code)
