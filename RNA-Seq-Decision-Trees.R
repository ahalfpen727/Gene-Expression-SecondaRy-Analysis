#building classification trees with RNA-Seq data
# load libraries
library(rpart);library(rpart.plot);library(partykit);library(multtest)
library(cummeRbund);library(nortest);library(stats);library(stats4)
library(readr);library(tidyr);library(magrittr);library(tidyverse)
library(qvalue);library(RColorBrewer);library(VennDiagram)
library(ggplot2);library(gplots);library(gridExtra)
library(limma);library(multtest);library(dplyr)
library(edgeR);library(tree);library(plyr)
library(ALL);data(golub);library(rsample)
library(GenomicFeatures);library(genefilter)
library(caret);library(tree)
scrs=file.path("/media/drew/easystore/umb_triley/urine1/Sample-Library-Preparation/MichiganUrineSpecimensAUASIscores.csv")
AUASI<-read_csv(scrs,trim_ws = T,col_names = TRUE)
AUASI.df<-AUASI %>%
   mutate(Sample.Num = paste("Sample", Samples, sep="_"))
#lanes<-read_csv(file="/media/drew/easystore/umb_triley/urine1/Sample-Library-Preparation/lane-and-sample-numbers.csv", col_names = T,trim_ws = T)
pool<-read_delim(file="/media/drew/easystore/umb_triley/urine1/Sample-Library-Preparation/Pool-A-and-Pool-B-barcode-summary.csv",
                 trim_ws = T,  delim=",")
lanes<-read_delim(file="/media/drew/easystore/umb_triley/urine1/UrineScores.csv",trim_ws = T,  delim=",")
# sequence data
genome.fa="/media/drew/easystore/ReferenceGenomes/UCSC_hg38/genome.fa"
refgtf="/media/drew/easystore/ReferenceGenomes/UCSC_hg38/genes.gtf"
inDir="/media/drew/easystore/umb_triley/urine1/cuffdiff_results_hg38_default/LUTS-over-CTRL"
gtfDir="/media/drew/easystore/umb_triley/urine1/cuffcompare_results_hg38_gtf_guided"
refDir="/media/drew/easystore/ReferenceGenomes/UCSC_hg38"
reffile <- file.path(refDir,"genes.gtf")
gtffile <- file.path(gtfDir,"cuffcmp.combined.gtf")
cuffcmp="/media/drew/easystore/umb_triley/urine1/cuffcompare_results_hg38_gtf_guided/cuffcmp.combined.gtf"
mergedgtf <- readGFF(cuffcmp)
hg38.genes.gtf<-as.data.frame(mergedgtf)
# novel gtf
mergedgtf <- readGFF(cuffcmp)
novelhg38.genes.gtf<-as.data.frame(mergedgtf)
novelmerged<-novelhg38.genes.gtf[which(novelhg38.genes.gtf["class_code"] != "="),]
novelmerged<-as.data.frame(mergedgtf)
novel.hg38.granges<-makeGRangesFromDataFrame(novelmerged, keep.extra.columns=TRUE)
novelhg38.granges<-makeGRangesFromDataFrame(novelhg38.genes.gtf, keep.extra.columns=TRUE)
ntxdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs=character())
seqlevels(ntxdb)<-seqlevels0(ntxdb)
nseq.txdb<-seqlevels(ntxdb)
ntbg <- transcriptsBy(ntxdb,by="gene")
#reference gtf
ref.gtf <- readGFF(refgtf)
refhg38.genes.gtf<-as.data.frame(ref.gtf)
ref.hg38.granges<-makeGRangesFromDataFrame(refhg38.genes.gtf, keep.extra.columns=TRUE)
rtxdb <- makeTxDbFromGFF(reffile, format="gtf", circ_seqs=character())
seqlevels(rtxdb)<-seqlevels0(rtxdb)
rseq.txdb<-seqlevels(rtxdb)
rtbg <- transcriptsBy(rtxdb,by="gene")
#reference
rtbg
#novel
ntbg
# gene expr
cuff<-cummeRbund::readCufflinks(dir=inDir,genome=genome.fa,gtfFile=refgtf, rebuild=F)
cuff
lanes
pool
Lane.factor <- as.factor(lanes$`Lanes_1:4_or_5:8`)
replicate.info<-cummeRbund::replicates(cuff)
cuffdir<-dirname(replicate.info$file)
sample_number<-basename(cuffdir)
replicates.info<-replicate.info[-1]
replicates.info$sample_number<-basename(cuffdir)
replicates.info$lanes<-c("L1T4","L1T4","L1T4","L5T8","L5T8","L5T8","L1T4","L5T8","L5T8","L1T4","L1T4","L1T4","L5T8","L5T8","L1T4","L1T4","L5T8","L5T8")
lane.info<-replicates.info$lanes
lanes <- as.factor(lane.info)
replicates.info$sample_number <- gsub("_out","",replicates.info$sample_number)
replicates.info<-replicates.info %>%
   mutate(batch.effects = paste(sample_name, lanes, sep="_"))
batcheffects <- as.factor(replicates.info$batch.effects)
replicates.info

grp.factor<-c(rep(0,9),rep(1,9))
grp.fac<-factor(grp.factor, levels=0:1, labels=c(under, over))
groups<-replicates.info$sample_name
samples<-replicates.info$rep_name
under=groups[1]
over=groups[((length(groups)/2)+1)]
conditions <- as.factor(groups)
#sampleCondition<-c(rep(under,length(groups)/2),rep(over,length(groups)/2))cuff_annotation_data<-featureNames(cuff@genes)
# gene expr data
cuff_annotation_data<-featureNames(cuff@genes)
gene_exp.diff<-diffData(cummeRbund::genes(cuff))
g.cnt.df<-repCountMatrix(cummeRbund::genes(cuff))
g.count.df=as.data.frame(g.cnt.df)
g.cnt.ma=as.matrix(g.cnt.df)
# filter NAs
g.cnt.df<-subset(g.cnt.df,!is.na(row.names(g.cnt.df)))
# factors for conditions
under.group<-grep(pattern=under, colnames(g.cnt.df))
over.group<-grep(pattern=over, colnames(g.cnt.df))
g.o.cnt.df<-g.cnt.df[,over.group]
g.u.cnt.df<-g.cnt.df[,under.group]
# set real number ceiling and floor
ma<-max(gene_exp.diff$log2_fold_change[is.finite(gene_exp.diff$log2_fold_change)])
mi<-min(gene_exp.diff$log2_fold_change[is.finite(gene_exp.diff$log2_fold_change)])
# set Inf and -Inf with ceiling and floor
gene_exp.diff$log2_fold_change<-replace(gene_exp.diff$log2_fold_change,
                                         gene_exp.diff$log2_fold_change == "Inf", ma)
gene_exp.diff$log2_fold_change<-replace(gene_exp.diff$log2_fold_change,
                                         gene_exp.diff$log2_fold_change == "-Inf", mi)
# sig gene expr data
sig_gene_exp.diff<-subset(gene_exp.diff, gene_exp.diff$significant=="yes")
mySigGenes<-getSig(cuff,x=over,y=under,alpha=0.05,level='genes')
sigGenes<-getGenes(cuff, mySigGenes)
# select directionally similar genes for each group
sig.h.gene_exp.diff<-subset(sig_gene_exp.diff,
                             sig_gene_exp.diff$log2_fold_change > 0
                             & sig_gene_exp.diff$q_value < 0.05)
sig.l.gene_exp.diff<-subset(sig_gene_exp.diff,
                             sig_gene_exp.diff$log2_fold_change < 0
                             & sig_gene_exp.diff$q_value < 0.05)

s.g.h.rep.matrix<-g.cnt.df[which(row.names(g.cnt.df) %in% sig.h.gene_exp.diff$gene_id),]
s.g.l.rep.matrix<-g.cnt.df[which(row.names(g.cnt.df) %in% sig.l.gene_exp.diff$gene_id),]
sig.gene.matrix<-g.cnt.df[which(row.names(g.cnt.df) %in% sig_gene_exp.diff$gene_id),]

sig.cnt.df<-g.cnt.df[which(row.names(g.cnt.df) %in% sig_gene_exp.diff$gene_id),]

s.g.h.rep.matrix<-s.g.h.rep.matrix[,over.group]
s.g.l.rep.matrix<-s.g.l.rep.matrix[,under.group]

# design matrix \
#design.sample <- model.matrix(~0 + rep_name, data=replicates.info)
design <- model.matrix(~batcheffects-1, data=replicates.info)
row.names(design) <- samples
# contrast matrix
contr.matrix <- makeContrasts(CTRL_L1T4 - LUTS_L1T4, CTRL_L5T8 - LUTS_L5T8, CTRL_L1T4 - CTRL_L5T8, LUTS_L1T4 - LUTS_L5T8,
                              levels = c("CTRL_L1T4", "CTRL_L5T8", "LUTS_L1T4", "LUTS_L5T8"))
design
dgel<- DGEList(counts=g.cnt.df, group=factor(lanes))
dge.norm <- calcNormFactors(dgel)
log2.cpm <- voom(dgel,design,plot=FALSE)
log2.cpm$design
fit.lm <- lmFit(log2.cpm,design)
fit.bayes <- eBayes(fit.lm)
f.bayes.dt <- decideTests(fit.bayes)
summary(f.bayes.dt)
# Fold-change thresholding
tfit <- treat(fit.bayes, lfc = 2)
# Q-Q plot of moderated t-statistics
qqt(fit.bayes$t[,2],df=fit.bayes$df.residual+fit.bayes$df.prior)
dgeObj <- DGEList(counts=g.cnt.df, group=groups)
disp <- estimateDisp(dgeObj,design)
plotBCV(disp, main="CV^2 as a Function of Counts/Million Mapped Reads")
# fisher exact test
et<-exactTest(disp)
diff_exp.et.df<-et$table
topTags(et)
o.sig <- subset(diff_exp.et.df, (PValue < 0.01))
dim(o.sig)
# Classic Approach after filtering
logcounts <- cpm(dgeObj,log=TRUE,prior.count=1)
str(logcounts)
## Identify genes with at least 1 cpm in at least 9 samples
keep.exprs <- rowSums(logcounts > 1) >= length(groups)/2
table(keep.exprs)
# Subset the rows of countdata to keep the more highly expressed genes
g.f.cnt.ma <- g.cnt.df[keep.exprs,]
dgeObj=DGEList(g.f.cnt.ma,group=groups)
f.disp <- estimateDisp(dgeObj,design)
plotMeanVar(f.disp)
plotBCV(f.disp, main="CV^2 as a Function of Counts/Million Mapped Reads")
levels(f.disp$samples$group)
# fisher exact test
et<-exactTest(f.disp)
diff_exp.et.df<-et$table
topTags(et)
o.sig.filt <- subset(diff_exp.et.df, (PValue < 0.05))
dim(o.sig.filt)
dge.norm <- calcNormFactors(dgeObj)
plotMDS(dge.norm, method="bcv", col=as.numeric(dge.norm$samples$group, main="Distance based on CV^2"))
plotMDS(dge.norm, method="logFC", col=as.numeric(dge.norm$samples$group, main="Distance based on logFC"))
d <- estimateGLMCommonDisp(dgeObj, design, verbose=TRUE)
d$common.dispersion
dgeObj.disp <- estimateDisp(dgeObj, design)
plotBCV(dgeObj.disp)

glmfit <- glmFit(dgeObj.disp, design, dispersion=d$common.dispersion)
glmfit
glm.lrt <- glmLRT(glmfit, contrast=contr.matrix)
# overdispersion
dispers<-deviance(glm.lrt)/df.residual(glm.lrt)
over.disp<-dispers[which(dispers > 1)]
length(over.disp)
not.over.disp<-dispers[which(dispers < 1)]
length(not.over.disp)
# Fix overdisperion
q.fit <- glmFit(dgeObj.disp, design, dispersion=d$common.dispersion, family="quasipoisson")
head(coef(q.fit))
### Adjusted p-values from limma
g.f.cnt.fit.lm = lmFit(g.cnt.df,design)
g.f.cnt.fit.lm
g.f.cnt.ebayes.lm = eBayes(g.f.cnt.fit.lm)
g.f.cnt.lm.pvals = topTable(g.f.cnt.ebayes.lm,number=dim(g.cnt.ma)[1])$P.Value<0.01
dt<-decideTests(g.f.cnt.ebayes.lm,p.value = 0.01,lfc = 4,adjust.method = "fdr",)
dt
summary(dt)
replicates.info

# so that we calculate distances between patients
kclu=kmeans(t(sig.cnt.df),centers=2)
#kmclu=cluster::pam(t(sig.cnt.df),k=4) #  cluster using k-medoids
table(kclu$cluster)

goldata <- as.data.frame(t(sig.cnt.df),col.names=replicates.info$sample_name)
gol.fac <- factor(colnames(sig.cnt.df))
gol.rp  <-  rpart(gol.fac~., data=goldata, method="class")
predictedclass <- predict(gol.rp, type="class")
table(predictedclass, gol.fac)
summary(gol.rp)
predict(gol.rp,type="class")
predict(gol.rp, type="prob")
plot(gol.rp, branch=0,margin=0.1)
text(gol.rp, digits=3, use.n=TRUE)











n<-10 ; sigma <- 0.5
fac <- factor(c(rep(1,n),rep(2,n),rep(3,n)))
levels(fac) <- c("ALL1","ALL2","AML")
geneA <- c(rnorm(20,0,sigma),rnorm(10,2,sigma))
geneB <- c(rnorm(10,0,sigma),rnorm(20,2,sigma))
geneC <- c(rnorm(30,1,sigma))
dat <- data.frame(fac,geneA,geneB,geneC)
rp <- rpart(fac ~ geneA + geneB + geneC, method="class",data=dat)

plot(rp, branch=0,margin=0.1)
text(rp, digits=3, use.n=TRUE)




