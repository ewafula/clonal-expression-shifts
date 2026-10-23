#!/usr/bin/env Rscript
# Public, entirely artificial example. Does not read patient sample data or prior results.

for(pkg in c('optparse','SeuratObject','Matrix','ape'))if(!requireNamespace(pkg,quietly=TRUE))stop('Missing ',pkg,'. Run 00_install_packages.R.')
suppressPackageStartupMessages(library(Matrix))
opt<-optparse::parse_args(optparse::OptionParser(option_list=list(
 optparse::make_option('--output_dir',default='inputs/tutorial'),
 optparse::make_option('--seed',type='integer',default=2026L)
)))
if(is.na(opt$seed)||opt$seed<0L)stop('Use a nonnegative integer seed.')
files<-file.path(opt$output_dir,c('expression.rds','clones.tsv','tree.nwk','genes.txt','README.txt'))
if(any(file.exists(files)))stop('Example inputs exist; use a different --output_dir. Nothing overwritten.')
if(!dir.exists(opt$output_dir))dir.create(opt$output_dir,recursive=TRUE)
set.seed(opt$seed)
clone<-rep(c('A','B','C','D'),each=32L);cells<-paste0('synthetic_',clone,'_',sprintf('%03d',sequence(rep(32L,4))))
genes<-c('DEMO-CLONE-A','DEMO-CLONE-C-SPARSE','DEMO-WITHIN-B','DEMO-BROAD','DEMO-WEAK')
features<-c(genes,paste0('BACKGROUND-',seq_len(20)))
counts<-matrix(rpois(length(features)*length(cells),lambda=5),nrow=length(features),dimnames=list(features,cells))
counts[1,]<-rpois(length(cells),ifelse(clone=='A',18,2))
counts[2,]<-rpois(length(cells),ifelse(clone=='C',1.4,.03))
counts[3,]<-rpois(length(cells),ifelse(clone=='B'&seq_along(cells)%%32 %in% 1:16,15,1))
counts[4,]<-rpois(length(cells),9)
counts[5,]<-rpois(length(cells),.25)
meta<-data.frame(cell_id=cells,clone_id=clone,row.names=cells)
obj<-SeuratObject::CreateSeuratObject(counts=as(counts,'dgCMatrix'),meta.data=meta,project='SyntheticTutorial',min.cells=0,min.features=0)
normalized<-log1p(t(t(counts)/colSums(counts))*10000)
obj<-SeuratObject::SetAssayData(obj,assay='RNA',layer='data',new.data=as(normalized,'dgCMatrix'))
saveRDS(obj,files[1]);utils::write.table(meta[,c('cell_id','clone_id')],files[2],sep='\t',quote=FALSE,row.names=FALSE)
balanced<-function(tips){if(length(tips)==1)return(tips);n<-length(tips)%/%2;paste0('(',balanced(tips[seq_len(n)]),',',balanced(tips[(n+1):length(tips)]),')')}
clades<-vapply(split(cells,clone),balanced,'')
newick<-paste0('((',clades['A'],',(',clades['B'],',',clades['C'],')),',clades['D'],');')
tree<-ape::read.tree(text=newick);stopifnot(ape::is.rooted(tree),setequal(tree$tip.label,cells))
ape::write.tree(tree,file=files[3]);writeLines(genes,files[4])
writeLines(c('ENTIRELY SYNTHETIC: 128 artificial cells, four monophyletic clones, five example genes and 20 background features.',
 paste('Seed:',opt$seed),'Normalization: natural-log log1p(count / per-cell library total * 10000).',
 'Artificial binary topology; no inferred ancestry, branch lengths, support values, patient identifiers or real observations.',
 'Genes illustrate clone differences, sparse detection, a designed within-B subdivision, broad expression and weaker expression.',
 'DEMO-BROAD and DEMO-WEAK are not independently validated null controls. Library normalization introduces compositional effects.',
 'No DEG table or synthetic p-values are supplied. Module 3 uses the explicit gene-list route.',
 'Use docs/TUTORIAL.md for the exact commands and verification guidance.'),files[5])
cat('Synthetic inputs written to',opt$output_dir,'\nNo private input was read.\n')
