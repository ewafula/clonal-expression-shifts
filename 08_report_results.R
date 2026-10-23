#!/usr/bin/env Rscript
# Interpret saved Module 7 results; no fitting or new statistical tests.
for(pkg in c('optparse','Matrix'))if(!requireNamespace(pkg,quietly=TRUE))stop('Missing package ',pkg,'. Run 00_install_packages.R.')

suppressPackageStartupMessages(library(Matrix))
opt<-optparse::parse_args(optparse::OptionParser(option_list=list(
  optparse::make_option('--panel_dir',default='results/module07'),
  optparse::make_option('--output_dir',default='results/module08')
)))
checkpoint<-file.path(opt$panel_dir,'module07.checkpoint.rds')
if(!file.exists(checkpoint))stop('Missing panel checkpoint: ',checkpoint)
state<-readRDS(checkpoint);panel<-state$signature$panel;status<-state$status
inputs<-names(state$signature$input_md5)
if(length(inputs)!=2||!all(file.exists(inputs))||!identical(unname(tools::md5sum(inputs)),unname(state$signature$input_md5)))stop('Prepared expression/tree inputs missing or changed since Module 7.')
e<-readRDS(inputs[1]);t<-readRDS(inputs[2])
if(!identical(e$tree_designs_md5,unname(tools::md5sum(inputs[2]))))stop('Expression/tree provenance mismatch.')
outnames<-c('module08.report.md','module08.gene_summary.tsv','module08.fit_metrics.tsv','module08.branch_stability.tsv','module08.contrasts.tsv','module08.sources.tsv','module08.overview.png','module08.session.txt')
if(any(file.exists(file.path(opt$output_dir,outnames))))stop('Report outputs exist; use a new --output_dir.')
if(!dir.exists(opt$output_dir))dir.create(opt$output_dir,recursive=TRUE)
fmt<-function(x,digits=3){x[is.finite(x)&abs(x)<.5*10^(-digits)]<-0;ifelse(is.na(x),'NA',formatC(x,format='f',digits=digits))}
escape<-function(x)gsub('|',' / ',gsub('[\r\n]',' ',as.character(x)),fixed=TRUE)
mdtable<-function(x)c(paste0('| ',paste(names(x),collapse=' | '),' |'),paste0('| ',paste(rep('---',ncol(x)),collapse=' | '),' |'),apply(x,1,function(r)paste0('| ',paste(escape(r),collapse=' | '),' |')),'')
relative<-function(path){a<-strsplit(normalizePath(path,mustWork=TRUE),'/',fixed=TRUE)[[1]];b<-strsplit(normalizePath(opt$output_dir,mustWork=TRUE),'/',fixed=TRUE)[[1]];n<-0L;while(n<min(length(a),length(b))&&a[n+1]==b[n+1])n<-n+1L;parts<-c(rep('..',length(b)-n),if(n<length(a))a[(n+1):length(a)]);paste(vapply(parts,utils::URLencode,'',reserved=TRUE),collapse='/')}
link<-function(label,path)if(file.exists(path))paste0('[',label,'](',relative(path),')') else paste0(label,' (unavailable)')
sources<-data.frame(path=character(),md5=character(),role=character())
record<-function(paths,role){paths<-unique(paths[file.exists(paths)]);if(length(paths))sources<<-rbind(sources,data.frame(path=paths,md5=unname(tools::md5sum(paths)),role=role))}
record(c(checkpoint,inputs),'checkpoint/prepared input')
record(c('08_report_results.R','docs/RESULTS_GUIDE.md',file.path(opt$panel_dir,c('module07.panel.tsv','module07.status.tsv','module07.summary.tsv','module07.report.txt'))),'report code/dictionary/panel tables')
latest<-status[!duplicated(status[c('gene','tree_mode','stage')],fromLast=TRUE),]
summaries<-list();fits<-list();branches<-list();contrasts<-list();sections<-character();scales<-character()
for(i in seq_len(nrow(panel))){
 g<-panel$gene[i];y<-as.numeric(e$expression[g,]);clone<-as.character(e$cell_metadata$clone_id)
 if(any(!is.finite(y)))stop('Nonfinite prepared response: ',g)
 means<-tapply(y,clone,mean);detection<-tapply(y!=0,clone,mean)
 sections<-c(sections,paste0('## ',g),'',paste('Selection reason:',panel$selection_reason[i],'.'),'',
   paste('Observed clone means range from',fmt(min(means)),'(Clone',names(which.min(means)),') to',fmt(max(means)),'(Clone',names(which.max(means)),'). Overall detection is',fmt(100*mean(y!=0),1),'percent. These are descriptive expression summaries, not ancestral changes.'),'',
   mdtable(data.frame(clone=names(means),cells=as.integer(table(clone)[names(means)]),mean=fmt(means),detection_percent=fmt(100*detection,1))))
 for(mode in state$signature$modes){
  rows<-latest[latest$gene==g&latest$tree_mode==mode,];frow<-rows[rows$stage=='fit',];vrow<-rows[rows$stage=='validation',]
  complete<-nrow(frow)==1&&nrow(vrow)==1&&frow$status=='success'&&vrow$status=='success'
  sections<-c(sections,paste0('### ',mode), '')
  if(!complete){sections<-c(sections,'Incomplete fit or validation. No numerical interpretation generated; inspect Module 7 stage logs.','');summaries[[length(summaries)+1]]<-data.frame(gene=g,tree_mode=mode,status='incomplete',n_cells=length(y),n_clones=length(means),en_selected=NA_integer_,en_selection_min=NA_real_,en_selection_max=NA_real_,en_success=NA_integer_,linear_success=NA_integer_,linear_not_estimable=NA_integer_,contrast_success_min=NA_integer_,contrast_success_max=NA_integer_,clone_R2=NA_real_,en_R2=NA_real_,linear_R2=NA_real_,linear_gain_over_clone=NA_real_,best_BIC_model=NA_character_,best_BIC=NA_real_,clone_delta_BIC=NA_real_);next}
  fp<-frow$output_prefix;vp<-vrow$output_prefix
  required<-c(paste0(fp,c('.models.rds','.coefficients.tsv','.contrasts.tsv','.fit_metrics.tsv','.comparison.png','.selected_tree.png')),paste0(vp,c('.validation.rds','.stability.tsv','.iteration_status.tsv','.model_comparison.tsv','.stability.png')))
  if(!all(file.exists(required)))stop('Completed pair has missing files: ',g,' / ',mode)
  a<-readRDS(required[1]);v<-readRDS(paste0(vp,'.validation.rds'))
  if(!identical(v$pilot_md5,unname(tools::md5sum(required[1])))||!identical(v$expression_md5,unname(state$signature$input_md5[1]))||!identical(v$tree_designs_md5,unname(state$signature$input_md5[2]))||a$gene!=g||v$gene!=g||!identical(a$tree_mode,mode)||v$tree_mode!=mode||!identical(a$cell_predictions$cell_id,e$cell_metadata$cell_id)||!isTRUE(all.equal(y,a$cell_predictions$observed,tolerance=0)))stop('Saved pilot/validation provenance or alignment mismatch: ',g,' / ',mode)
  if(!identical(a$expression_md5,v$expression_md5)||!identical(a$tree_designs_md5,v$tree_designs_md5))stop('Pilot input hash mismatch.')
  for(pair in list(c(paste0(fp,'.coefficients.tsv'),'coefficients','a'),c(paste0(fp,'.contrasts.tsv'),if(mode=='uncollapsed')'all_contrasts' else 'contrasts','a'),c(paste0(fp,'.fit_metrics.tsv'),'fit_metrics','a'),c(paste0(vp,'.stability.tsv'),'stability','v'),c(paste0(vp,'.iteration_status.tsv'),'iteration_status','v'),c(paste0(vp,'.model_comparison.tsv'),'model_comparison','v'))){
    disk<-utils::read.delim(pair[1],check.names=FALSE,colClasses='character',na.strings='NA')
    saved<-if(pair[3]=='a')a[[pair[2]]] else v[[pair[2]]]
    if(!identical(names(disk),names(saved))||nrow(disk)!=nrow(saved))stop('Linked TSV schema differs from saved result: ',pair[1])
    for(field in names(saved)){
      x<-saved[[field]];d<-disk[[field]]
      equal<-if(is.numeric(x))isTRUE(all.equal(as.numeric(d),as.numeric(x),tolerance=1e-10)) else identical(as.character(d),as.character(x))
      if(!equal)stop('Linked TSV differs from saved result: ',pair[1],' / ',field)
    }
  }
  record(required,paste(g,mode));record(paste0(fp,c('.comparison.png','.selected_tree.png')),paste(g,mode,'figure'));record(paste0(vp,'.stability.png'),paste(g,mode,'figure'))
  scales<-unique(c(scales,a$expression_scale))
  m<-a$fit_metrics;s<-v$stability;it<-v$iteration_status;ic<-v$model_comparison
  if(any(s$successful+s$failed+s$not_estimable!=s$attempted))stop('Invalid stability denominators: ',g)
  branch<-s[s$model!='contrasts',];dd<-t$branch_dictionary[t$branch_dictionary$mode==mode,]
  branch$descendant_clones<-dd$descendant_clones[match(branch$term,dd$branch_id)]
  branch$splits_within_clone<-dd$splits_within_clone[match(branch$term,dd$branch_id)]
  branches[[length(branches)+1]]<-cbind(gene=g,tree_mode=mode,branch)
  ct<-a$contrasts;ctid<-if(mode=='collapsed')ct$branch_id else ct$contrast_id
  cs<-s[s$model=='contrasts',];idx<-match(ctid,cs$term)
  contrasts[[length(contrasts)+1]]<-data.frame(gene=g,tree_mode=mode,contrast_id=ctid,branch_id=ct$branch_id,contrast=ct$contrast,clone_id=if(mode=='collapsed')NA_character_ else ct$clone_id,contrast_type=if(mode=='collapsed')'equal_clone' else ct$contrast_type,difference=ct$difference,successful=cs$successful[idx],attempted=cs$attempted[idx],failed=cs$failed[idx],not_estimable=cs$not_estimable[idx],median=cs$median[idx],q05=cs$q05[idx],q95=cs$q95[idx])
  fits[[length(fits)+1]]<-cbind(gene=g,tree_mode=mode,m)
  en<-s[s$model=='elastic_net'&s$full_selected,];cc<-a$coefficients[a$coefficients$model=='elastic_net'&a$coefficients$term!='(Intercept)'&a$coefficients$selected,]
  good<-ic$status=='working_iid_cell_gaussian_ML'&is.finite(ic$BIC);z<-ic[good,];best<-if(nrow(z))z[which.min(z$BIC),] else NULL
  ref<-z[z$model_id=='clone_means_reference',];getr<-function(model){x<-m$weighted_r_squared_apparent[m$model==model];if(length(x))x else NA_real_}
  count<-function(model,stat){if(!any(it$model==model))NA_integer_ else sum(it$model==model&it$status==stat)}
  enit<-it[it$model=='elastic_net',];ecount<-count('elastic_net','success');lcount<-count('linear','success');lnon<-count('linear','not_estimable')
  cmin<-if(nrow(cs))min(cs$successful) else NA_integer_;cmax<-if(nrow(cs))max(cs$successful) else NA_integer_
  bname<-if(is.null(best))NA_character_ else best$model_id;bval<-if(is.null(best))NA_real_ else best$BIC;delta<-if(nrow(ref)&&!is.null(best))ref$BIC-bval else NA_real_
  summaries[[length(summaries)+1]]<-data.frame(gene=g,tree_mode=mode,status='complete',n_cells=length(y),n_clones=length(means),en_selected=nrow(cc),en_selection_min=if(nrow(en))min(en$selection_frequency) else NA_real_,en_selection_max=if(nrow(en))max(en$selection_frequency) else NA_real_,en_success=ecount,linear_success=lcount,linear_not_estimable=lnon,contrast_success_min=cmin,contrast_success_max=cmax,clone_R2=getr('clone_means_reference'),en_R2=getr('elastic_net'),linear_R2=getr('linear'),linear_gain_over_clone=getr('linear')-getr('clone_means_reference'),best_BIC_model=bname,best_BIC=bval,clone_delta_BIC=delta)
  sections<-c(sections,paste('Validation settings:',v$options$reps,'subsamples;',fmt(100*v$options$fraction,1),'percent per clone without replacement; EN alpha',a$options$alpha,'and',a$options$lambda_choice,'.'),'',paste('Training clone-balanced R-squared: clone means',fmt(getr('clone_means_reference')),'; Elastic Net',fmt(getr('elastic_net')),'; linear',fmt(getr('linear')),'. Linear minus clone means:',fmt(getr('linear')-getr('clone_means_reference')),'. Greater flexibility can improve training fit without improving independent prediction.'),'',
   paste('Elastic Net selected',nrow(cc),'branches in the full fit. Successful EN resamples:',ecount,'of',nrow(enit),'. Linear resamples:',lcount,'successful and',lnon,'non-estimable. Contrast-specific usable resamples range from',cmin,'to',cmax,'. Read term-specific denominators; ranges below are not confidence intervals.'),'')
  if(nrow(en)){
   show<-en[order(-abs(en$full_value)),];show<-head(show,10)
   sections<-c(sections,'Largest full-fit selected EN contributions (up to ten; all branches are in the combined table):','',mdtable(data.frame(branch=show$term,coefficient=fmt(show$full_value),selected_fraction=fmt(show$selection_frequency,2),successful=show$successful,median=fmt(show$median),q05=fmt(show$q05),q95=fmt(show$q95))))
  }else sections<-c(sections,'No full-fit EN branch selected (or EN was not fitted); this does not establish absence of expression differences.','')
  if(!is.null(best))sections<-c(sections,paste('Lowest evaluated working BIC:',bname,'=',fmt(bval),'; full clone-mean model delta BIC =',fmt(delta),'. Lowest evaluated AIC model:',z$model_id[which.min(z$AIC)],'. Equivalent prediction spaces can have different branch names. This comparison is within this gene/mode; candidate searches may differ between modes.'),'') else sections<-c(sections,'No valid working IC comparison available.','')
  if(nrow(ct)){
   show<-head(ct[order(-abs(ct$difference)),],3);ix<-match(if(mode=='collapsed')show$branch_id else show$contrast_id,cs$term)
   sections<-c(sections,'Largest absolute descriptive contrasts (not an independent significance ranking):','',mdtable(data.frame(contrast=show$contrast,difference=fmt(show$difference),usable=cs$successful[ix],attempted=cs$attempted[ix])))
  }
  if(mode=='uncollapsed'){
   availability<-table(a$all_contrasts$status)
   sections<-c(sections,paste('All contrast requests:',paste(paste(names(availability),as.integer(availability),sep='='),collapse='; '),'. Combined contrast table includes ready contrasts only; unavailable requests remain in the linked original table. No within-clone ancestry is established.'),'')
  }
  sections<-c(sections,paste(link('Coefficients',paste0(fp,'.coefficients.tsv')),link('All contrast requests',paste0(fp,'.contrasts.tsv')),link('Stability',paste0(vp,'.stability.tsv')),link('Iteration status',paste0(vp,'.iteration_status.tsv')),link('Working IC',paste0(vp,'.model_comparison.tsv')),sep=' · '),'',
    paste0('![',g,' ',mode,' expression and training diagnostics](',relative(paste0(fp,'.comparison.png')),')'),'','The comparison figure shows observed expression, training fits and descriptive contrasts; it is not independent validation.','',
    paste0('![',g,' ',mode,' resampling stability](',relative(paste0(vp,'.stability.png')),')'),'','The stability figure shows selection frequencies and subsampling ranges. Ranges are not confidence intervals. Full tables retain terms omitted from plots.','')
  if(mode=='collapsed')sections<-c(sections,link('Collapsed selected clone tree',paste0(fp,'.selected_tree.png')),'','Black branches are unselected; red/blue mark positive/negative simultaneous EN contributions. Labels identify original branch terms and coefficients. This topology view does not identify unique historical shifts.','')
  if(mode=='uncollapsed')sections<-c(sections,link('Exploratory selected cellular tree',paste0(fp,'.selected_tree.png')),'','Black branches show the supplied cellular topology; red/blue mark positive/negative EN contributions. The tip strip identifies clones. Selection does not validate ancestry.','')
 }
}
summary<-do.call(rbind,summaries);bind<-function(x)if(length(x))do.call(rbind,x) else data.frame(gene=character(),tree_mode=character())
write_table<-function(x,name)utils::write.table(x,file.path(opt$output_dir,name),sep='\t',quote=TRUE,row.names=FALSE,na='NA')
write_table(summary,'module08.gene_summary.tsv');write_table(bind(fits),'module08.fit_metrics.tsv');write_table(bind(branches),'module08.branch_stability.tsv');write_table(bind(contrasts),'module08.contrasts.tsv');write_table(unique(sources),'module08.sources.tsv')
# Descriptive cross-mode plot, with no significance or support threshold.
grDevices::png(file.path(opt$output_dir,'module08.overview.png'),width=1800,height=max(1200,180*nrow(panel)),res=150)
tryCatch({graphics::par(mar=c(5,10,5,2));ids<-paste(summary$gene,summary$tree_mode,sep=' / ');xx<-rbind(summary$clone_R2,summary$en_R2,summary$linear_R2);colnames(xx)<-ids
 if(any(is.finite(xx))){graphics::barplot(xx,horiz=TRUE,beside=TRUE,las=1,col=c('grey45','steelblue','darkorange3'),cex.names=.7,xlab='Clone-balanced training R-squared (apparent)',main='Training fit across genes and tree modes',ylim=c(0,nrow(summary)*4+3));graphics::legend('top',legend=c('Clone means','Elastic Net','Linear'),fill=c('grey45','steelblue','darkorange3'),horiz=TRUE,cex=.8,bty='n')} else {graphics::plot.new();graphics::text(.5,.5,'No completed gene/mode pairs')}
},finally=grDevices::dev.off())
complete<-summary$status=='complete'
cross<-c('## Comparison between tree modes','')
for(g in panel$gene){
  z<-summary[summary$gene==g&summary$status=='complete',]
  if(all(c('collapsed','uncollapsed') %in% z$tree_mode)){
    a<-z[z$tree_mode=='collapsed',];b<-z[z$tree_mode=='uncollapsed',]
    cross<-c(cross,paste(g,': uncollapsed minus collapsed apparent R-squared is',fmt(b$en_R2-a$en_R2),'for EN and',fmt(b$linear_R2-a$linear_R2),'for linear fitting. Full-fit selected EN branches:',a$en_selected,'collapsed versus',b$en_selected,'uncollapsed. These counts use different predictor dictionaries and do not count unique historical events.'),'')
  }
}
cross<-c(cross,'Both modes use the same prepared response/cells and clone-balanced weights. Compare observed clone means and explicit contrast definitions; penalized coefficients can redistribute across nested branches and cannot be matched by size alone. IC searches differ between modes, so differing winning models need not indicate a contradiction.','')
intro<-c('# Integrated clonal-expression report','',paste('Generated:',format(Sys.time(),tz='UTC',usetz=TRUE)), '',paste('Completed gene/mode pairs:',sum(complete),'of',nrow(summary),'. Incomplete pairs are shown explicitly; this report does not refit models.'),'',
 '## Applicability to other datasets','',
 'Module 8 is designed for any dataset processed through compatible versions of this workflow. Gene names, clone IDs, expression patterns, selected branches, stability summaries and IC results are read or computed from the saved run; the current five-gene example is not hard-coded. It requires compatible Modules 2–3 prepared inputs and Module 7 outputs, with referenced files available and matching. It does not directly interpret arbitrary outputs from unrelated software. Use --panel_dir for another completed panel and --output_dir for a separate report. This implementation has been checked on the reviewed five-gene dataset, an independent artificial tutorial dataset and incomplete/mismatched-result tests. Broader validation across independent biological datasets remains outstanding.','',
 '## How to interpret this report','',
 'Collapsed analysis is primary; uncollapsed cellular analysis is exploratory. Branch selection and fixed-tree cell-subsampling stability do not establish ancestral events, independent biological replication, or causality. No confirmatory p-values, confidence intervals or panel-wide discoveries are supplied.','',
 paste('Response scale:',paste(scales,collapse='; '),'. Coefficients and contrasts are not automatically DEG log2 fold changes.'),'',
 'Fitting uses inverse-clone-size balancing weights (each clone contributes equal total weight), not known precision weights. EN coefficients are simultaneous penalized contributions; linear coefficients depend on the saved basis; contrasts compare specified groups. They are different estimands. Tuning CV interpolates cells within clones, not held-out patients or clones.','',
 'Resampling samples a saved fraction of cells within each clone without replacement and retunes EN; the exact fraction and repetitions appear in each section. Read per-term successful/failed/non-estimable denominators. Quantiles include zero estimates and are descriptive subsampling ranges. Overlap-neighborhood frequencies, where supplied, use Jaccard ≥0.8 including the branch itself; they are not probabilities of a historical region.','',
 'Working IC uses ordinary unweighted iid homoskedastic Gaussian cell-error ML refits on the same response/cells, with k = design rank + one residual-variance parameter. This differs from the clone-balanced fit objective. Cell dependence, confounding and post-selection are not corrected. Compare IC only within a gene and evaluated candidate set; never rank genes by raw AIC/BIC. Uncollapsed searches are bounded. No forced agreement between AIC/BIC and EN is expected.','',
 'The supplied root and node-label provenance remain unverified. No sample/site adjustment, pseudobulking, or pathway enrichment is added here. Generic small-panel selection is descriptive among prepared candidates, not an independent null-control design.','',
 paste('Prepared expression:',nrow(e$expression),'genes by',ncol(e$expression),'cells. Retained clones:',length(unique(e$cell_metadata$clone_id)),'. Clone monophyly checks passed for',sum(t$clone_checks$monophyletic),'of',nrow(t$clone_checks),'clones. Original filters/settings remain in the saved Modules 2–3 inputs.'),'',
 '## Panel overview','',mdtable(panel[,c('gene','selection_reason','fraction_nonzero','clone_mean_range')]),
 '![Training fit overview](module08.overview.png)','',
 'Bars compare apparent training fit within each gene. Higher training R-squared in a more flexible model is not evidence of better independent prediction. Genes can have different response scales and variability.','',
 mdtable(data.frame(gene=summary$gene,mode=summary$tree_mode,status=summary$status,EN_branches=summary$en_selected,EN_min_frequency=fmt(summary$en_selection_min,2),linear_gain=fmt(summary$linear_gain_over_clone),clone_delta_BIC=fmt(summary$clone_delta_BIC))),
 'A positive clone delta BIC means an evaluated model has lower working BIC than full clone means; zero includes equivalent clone-mean parameterizations. It does not establish biological ancestry. The linear gain is the apparent R-squared difference from clone means, not an independent test.','')
tail<-c('## Traceability and reading the tables','',
 'Panel source: ',link('Module 7 panel',file.path(opt$panel_dir,'module07.panel.tsv')),link('Module 7 stage history',file.path(opt$panel_dir,'module07.status.tsv')),'',
 link('Gene/mode overview',file.path(opt$output_dir,'module08.gene_summary.tsv')),link('Training metrics',file.path(opt$output_dir,'module08.fit_metrics.tsv')),link('All branch stability estimates',file.path(opt$output_dir,'module08.branch_stability.tsv')),link('Ready contrasts',file.path(opt$output_dir,'module08.contrasts.tsv')),link('Source hashes',file.path(opt$output_dir,'module08.sources.tsv')),'',
 link('Column definitions and worked examples','docs/RESULTS_GUIDE.md'),'','Keep this report with the input panel folders so relative figure/table links work.','')
writeLines(c(intro,cross,sections,tail),file.path(opt$output_dir,'module08.report.md'),useBytes=TRUE)
capture.output(sessionInfo(),file=file.path(opt$output_dir,'module08.session.txt'))
cat('REPORT COMPLETE:',file.path(opt$output_dir,'module08.report.md'),'\n',sum(complete),'of',nrow(summary),'pairs complete. No models refitted.\n')
