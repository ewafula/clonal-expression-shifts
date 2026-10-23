#!/usr/bin/env Rscript
# Small descriptive panel; reuse the single-gene scripts unchanged.

for(pkg in c('optparse','Matrix')) if(!requireNamespace(pkg,quietly=TRUE))
  stop('Missing package ',pkg,'. Run Rscript 00_install_packages.R.')
suppressPackageStartupMessages(library(Matrix))
opt<-optparse::parse_args(optparse::OptionParser(option_list=list(
  optparse::make_option('--expression',default='results/module03.expression.rds'),
  optparse::make_option('--tree_designs',default='results/module02.tree_designs.rds'),
  optparse::make_option('--genes',default=NULL,help='Optional comma-separated prepared genes; default descriptive panel of up to five'),
  optparse::make_option('--gene_list',default=NULL,help='Headerless text file: one prepared gene per line; alternative to --genes'),
  optparse::make_option('--tree_mode',default='both',help='both, collapsed, or uncollapsed'),
  optparse::make_option('--reps',type='integer',default=100L),
  optparse::make_option('--seed',type='integer',default=12345L),
  optparse::make_option('--output_dir',default='results/module07',help='Panel output folder; each gene gets its own subfolder [default %default]'),
  optparse::make_option('--resume',action='store_true',default=FALSE)
)))
if(!is.null(opt$genes)&&!is.null(opt$gene_list))stop('Use either --genes or --gene_list, not both.')
if(!opt$tree_mode %in% c('both','collapsed','uncollapsed'))stop('Invalid --tree_mode.')
if(is.na(opt$reps)||opt$reps<1||is.na(opt$seed)||opt$seed<0||as.double(opt$seed)+opt$reps>.Machine$integer.max)stop('Invalid reps/seed.')
scripts<-c(collapsed='04_fit_collapsed_cells.R',uncollapsed='06_fit_uncollapsed_cells.R',validation='05_validate_models.R')
for(f in c(opt$expression,opt$tree_designs,scripts))if(!file.exists(f))stop('Missing file: ',f,'. Run from the workflow directory.')
e<-readRDS(opt$expression);t<-readRDS(opt$tree_designs)
if(!identical(e$tree_designs_md5,unname(tools::md5sum(opt$tree_designs))))stop('Expression/tree inputs do not match.')
if(!identical(colnames(e$expression),e$cell_metadata$cell_id))stop('Expression/metadata order differs.')
clone<-as.character(e$cell_metadata$clone_id);genes<-rownames(e$expression)
# Describe the available prepared genes; these are not independent negative controls.
profile<-do.call(rbind,lapply(genes,function(g){y<-as.numeric(e$expression[g,]);means<-tapply(y,clone,mean)
  data.frame(gene=g,mean=mean(y),fraction_nonzero=mean(y!=0),sd=sd(y),clone_mean_range=diff(range(means)))}))
valid<-with(profile,is.finite(sd)&sd>0&is.finite(clone_mean_range))
profile<-profile[valid,];panel<-profile[FALSE,];panel$selection_reason<-character()
add<-function(g,reason){if(length(g)&&!is.na(g)&&!g %in% panel$gene){row<-profile[match(g,profile$gene),,drop=FALSE];row$selection_reason<-reason;panel<<-rbind(panel,row)}}
if(!is.null(opt$genes)||!is.null(opt$gene_list)){
  if(!is.null(opt$gene_list)){
    if(!file.exists(opt$gene_list))stop('Missing gene-list file: ',opt$gene_list)
    requested<-trimws(readLines(opt$gene_list,warn=FALSE,encoding='UTF-8'))
    requested<-sub('^\ufeff','',requested)
    requested<-requested[nzchar(requested)]
    if(any(grepl('[[:space:],]',requested)))stop('--gene_list requires exactly one gene per line, without a header or other columns.')
  }else requested<-trimws(strsplit(opt$genes,',',fixed=TRUE)[[1]])
  requested<-unique(requested)
  if(!length(requested)||any(!nzchar(requested)))stop('Supply a nonempty gene list.')
  missing<-setdiff(requested,profile$gene)
  if(length(missing))stop('Genes not prepared or not finite/variable: ',paste(missing,collapse=', '),'. Prepare them with Module 3 first.')
  for(g in requested)add(g,'user-specified prepared gene')
}else{
  add(if('CLDN16' %in% profile$gene)'CLDN16' else profile$gene[1],'pilot anchor (CLDN16 when prepared; otherwise first available)')
  rules<-list(c('clone_mean_range','decreasing','largest clone-mean range'),c('fraction_nonzero','increasing','lowest detection fraction'),c('fraction_nonzero','decreasing','highest detection fraction'),c('clone_mean_range','increasing','smallest clone-mean range among prepared genes'))
  for(rule in rules){available<-profile[!profile$gene %in% panel$gene,];if(nrow(available)){idx<-order(if(rule[2]=='decreasing')-available[[rule[1]]] else available[[rule[1]]],available$gene);add(available$gene[idx[1]],rule[3])}}
}
if(!nrow(panel))stop('No usable prepared genes.')
panel$panel_index<-seq_len(nrow(panel))
# Preserve ordinary gene symbols; replace unsafe path characters without changing gene IDs.
folder<-gsub('[^A-Za-z0-9._-]','_',panel$gene)
if(any(folder %in% c('', '.', '..')) || anyDuplicated(tolower(folder)))
  stop('Gene symbols produce empty or conflicting folder names; choose an unambiguous panel.')
panel$gene_directory<-file.path(opt$output_dir,folder)
opt$output_prefix<-file.path(opt$output_dir,'module07')
modes<-if(opt$tree_mode=='both')c('collapsed','uncollapsed') else opt$tree_mode
for(mode in modes)if(is.null(t$designs[[mode]]$cell_X)||!identical(rownames(t$designs[[mode]]$cell_X),e$cell_metadata$cell_id))stop('Missing or misaligned design: ',mode)
state_path<-paste0(opt$output_prefix,'.checkpoint.rds')
summary_path<-paste0(opt$output_prefix,'.summary.tsv');status_path<-paste0(opt$output_prefix,'.status.tsv')
panel_path<-paste0(opt$output_prefix,'.panel.tsv');report_path<-paste0(opt$output_prefix,'.report.txt')
signature<-list(layout_version=2L,input_md5=tools::md5sum(c(opt$expression,opt$tree_designs)),script_md5=tools::md5sum(scripts),panel=panel,modes=modes,reps=opt$reps,seed=opt$seed)
if(opt$resume){
  if(!file.exists(state_path))stop('No checkpoint to resume.')
  state<-readRDS(state_path)
  if(!identical(state$signature,signature))stop('Inputs, scripts, panel or options changed; use a new output prefix.')
}else{
  if(any(file.exists(c(state_path,summary_path,status_path,panel_path,report_path))))stop('Panel outputs exist; use --resume or a new prefix.')
  state<-list(signature=signature,status=data.frame(gene=character(),tree_mode=character(),stage=character(),attempt=integer(),status=character(),elapsed_seconds=numeric(),exit_code=integer(),output_prefix=character(),log_file=character(),message=character()),sessionInfo=sessionInfo())
}
if(!dir.exists(dirname(opt$output_prefix)))dir.create(dirname(opt$output_prefix),recursive=TRUE)
write_table<-function(x,path)utils::write.table(x,path,sep='\t',quote=TRUE,row.names=FALSE,na='NA')
save_state<-function(){saveRDS(state,state_path);write_table(state$status,status_path);write_table(panel,panel_path)}
save_state()
run_stage<-function(g,mode,stage,script,args,base,suffixes){
  previous<-state$status[state$status$gene==g & state$status$tree_mode==mode & state$status$stage==stage,]
  if(nrow(previous)&&tail(previous$status,1)=='success'){
    prefix<-tail(previous$output_prefix,1)
    if(!all(file.exists(paste0(prefix,suffixes))))stop('Completed outputs missing: ',prefix,'. Restore them or use a new panel prefix.')
    return(prefix)
  }
  attempt<-if(nrow(previous))max(previous$attempt)+1L else 1L
  prefix<-paste0(base,'_attempt',attempt)
  # Interrupted runs may leave outputs without a recorded status; preserve them.
  while(any(file.exists(paste0(prefix,c(suffixes,'.console.txt'))))){attempt<-attempt+1L;prefix<-paste0(base,'_attempt',attempt)}
  if(!dir.exists(dirname(prefix)))dir.create(dirname(prefix),recursive=TRUE)
  logfile<-paste0(prefix,'.console.txt');started<-proc.time()[['elapsed']]
  cat(g,mode,stage,'->',prefix,'\n')
  message<-''
  code<-tryCatch(system2(file.path(R.home('bin'),'Rscript'),vapply(c(script,args,'--output_prefix',prefix),shQuote,''),stdout=logfile,stderr=logfile),error=function(err){message<<-conditionMessage(err);1L})
  good<-identical(as.integer(code),0L)&&all(file.exists(paste0(prefix,suffixes)))
  if(!good&&!nzchar(message))message<-'Child script failed or expected outputs missing; inspect console log.'
  state$status<<-rbind(state$status,data.frame(gene=g,tree_mode=mode,stage=stage,attempt=attempt,status=if(good)'success' else 'failed',elapsed_seconds=proc.time()[['elapsed']]-started,exit_code=as.integer(code),output_prefix=prefix,log_file=logfile,message=message))
  save_state()
  if(good)prefix else NA_character_
}
fit_suffix<-c('.models.rds','.report.txt','.coefficients.tsv','.contrasts.tsv','.fit_metrics.tsv','.cell_predictions.tsv.gz','.comparison.png','.selected_tree.png')
val_suffix<-c('.validation.rds','.report.txt','.stability.tsv','.resamples.tsv.gz','.iteration_status.tsv','.model_comparison.tsv','.stability.png')
rows<-list()
for(i in seq_len(nrow(panel)))for(mode in modes){
  g<-panel$gene[i];base<-file.path(panel$gene_directory[i],paste0(folder[i],'_',mode))
  fit<-run_stage(g,mode,'fit',scripts[[mode]],c('--expression',opt$expression,'--tree_designs',opt$tree_designs,'--gene',g,'--seed',opt$seed),base, fit_suffix)
  if(is.na(fit))next
  val<-run_stage(g,mode,'validation',scripts[['validation']],c('--pilot',paste0(fit,'.models.rds'),'--expression',opt$expression,'--tree_designs',opt$tree_designs,'--reps',opt$reps,'--seed',opt$seed),paste0(base,'_validation'),val_suffix)
  if(is.na(val))next
  a<-readRDS(paste0(fit,'.models.rds'));v<-readRDS(paste0(val,'.validation.rds'))
  if(v$pilot_md5!=unname(tools::md5sum(paste0(fit,'.models.rds'))))stop('Validation/pilot hash mismatch: ',val)
  metrics<-a$fit_metrics;s<-v$stability;ic<-v$model_comparison
  good<-ic$status=='working_iid_cell_gaussian_ML' & is.finite(ic$BIC)
  best<-if(any(good))ic[which(good)[which.min(ic$BIC[good])],] else NULL
  en<-s[s$model=='elastic_net',];cc<-a$coefficients
  for(j in seq_len(nrow(metrics))){method<-metrics$model[j];ss<-v$iteration_status[v$iteration_status$model==method,]
    rows[[length(rows)+1L]]<-data.frame(gene=g,tree_mode=mode,model=method,n_cells=metrics$n_cells[j],n_clones=metrics$n_clones[j],weighted_rmse=metrics$weighted_rmse[j],weighted_r_squared_apparent=metrics$weighted_r_squared_apparent[j],n_selected_branches=if(method=='elastic_net')sum(cc$model==method & cc$term!='(Intercept)' & cc$selected) else NA_integer_,resamples_successful=if(nrow(ss))sum(ss$status=='success') else NA_integer_,resamples_not_estimable=if(nrow(ss))sum(ss$status=='not_estimable') else NA_integer_,resamples_failed=if(nrow(ss))sum(ss$status=='failed') else NA_integer_,full_selected_min_frequency=if(method=='elastic_net'&&any(en$full_selected))min(en$selection_frequency[en$full_selected]) else NA_real_,best_working_BIC_model=if(is.null(best))NA_character_ else best$model_id,best_working_BIC=if(is.null(best))NA_real_ else best$BIC,fit_prefix=fit,validation_prefix=val)
  }
  write_table(do.call(rbind,rows),summary_path)
}
if(!length(rows))write_table(data.frame(gene=character(),tree_mode=character(),model=character()),summary_path)
latest<-state$status[!duplicated(state$status[c('gene','tree_mode','stage')],fromLast=TRUE),]
complete<-sum(latest$stage=='validation' & latest$status=='success')
writeLines(c('MODULE 7 — GENE PANEL',paste('Completed gene/tree-mode validations:',complete,'of',nrow(panel)*length(modes)),
 'Reuses scripts 04/06 and shared validator 05; default fitting/tuning and 80% clone-stratified sampling retained.',
 'Panel patterns are relative to the prepared genes. DEG-selected genes are not independent negative controls.',
 'Training metrics are apparent fit; stability is conditional on fixed topology. No pooled p-values or significance claims.',
 'IC comparisons are within each gene/mode only. Do not rank genes by AIC/BIC; bounded searches may differ between modes.',
 'Elapsed time is per child process; output bytes below are file storage, not peak RAM. Peak RAM is not measured.',
 paste('Child output storage bytes:',sum(file.info(unique(unlist(lapply(seq_len(nrow(state$status)),function(i)Sys.glob(paste0(state$status$output_prefix[i],'.*'))))))$size,na.rm=TRUE)),
 'Review panel/status/summary tables and per-gene outputs before larger batches. Resume with the same command plus --resume.',
 'Failed or interrupted stages retry under a new attempt suffix; previous files are preserved.'),report_path)
save_state();cat('Panel finished:',complete,'of',nrow(panel)*length(modes),'validations. Review',report_path,'\n')
if(complete<nrow(panel)*length(modes))quit(status=1L)
