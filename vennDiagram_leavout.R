# Analysis of cross validation leaveout results
###################################################################################
install.packages("Vennerable", repos="http://R-Forge.R-project.org")
biocM(c("gageData", "DLBCL","Rcy3","globaltest","KEGG.db","Rsamtools", "ALL","All"))

### GTF_ONLY (original prediction) significantly differentialy expressed features
library(Vennerable); library(VennDiagram); library(reshape);
library(R.utils); library(grid); library(cummeRbund);library(MASS);
library(limma);require(stats)
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/vennDia.R")
http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/R_BioCondManual.html
#R_graphics_venn
cuff<-readCufflinks(gtfFile="~/genomes/GRCh38_gtf_only/chromosome_genome.fa.gtf/cuffcmp.combined.gtf",genome="~/genomes/GRCh38_gtf_only/chromosome_genome.fa.gtf/genome.fa",rebuild=F)
cuff
over="LUTS"
under="CTRL"

mySigGenes<-getSig(cuff,x=over,y=under,alpha=.05,level='genes')
sigGenes<-getGenes(cuff, mySigGenes)
sig_genes_exp.diff<-diffData(sigGenes)
sig_genes_annot<-sigGenes@annotation
sig_genes_exp.diff<-cbind2(sig_genes_annot["gene_short_name"], sig_genes_exp.diff)
sig_genes_gtf<-sig_genes_exp.diff$gene_short_name

mySigCDS<-getSig(cuff,x=over,y=under,alpha=.05,level='CDS')
sigCDS<-getGenes(cuff, mySigCDS)
sig_cds_exp.diff<-diffData(sigCDS)
sig_cds_annot<-sigCDS@annotation
sig_cds_exp.diff<-cbind2(sig_cds_annot["gene_short_name"], sig_cds_exp.diff)
sig_cds_gtf<-sig_cds_exp.diff$gene_short_name

mySigIsos<-getSig(cuff,x=over,y=under,alpha=.05,level='isoforms')
sigIsos<-getGenes(cuff, mySigIsos)
sig_isos_exp.diff<-diffData(sigIsos)
sig_isos_annot<-sigIsos@annotation
sig_isos_exp.diff<-cbind2(sig_isos_annot["gene_short_name"], sig_isos_exp.diff)
sig_isos_gtf<-sig_isos_exp.diff$gene_short_name

mySigTSS<-getSig(cuff,x=over,y=under,alpha=.05,level='TSS')
sigTSS<-getGenes(cuff, mySigTSS)
sig_tss_exp.diff<-diffData(sigTSS)
sig_tss_annot<-sigTSS@annotation
sig_tss_exp.diff<-cbind2(sig_tss_annot["gene_short_name"], sig_tss_exp.diff)
sig_tss_gtf<-sig_tss_exp.diff$gene_short_name

# formating gtf_only data for comparison with leaveout analysis

gtfISOS <- vector("list", length(sig_isos_gtf))
for (i in 1:length(gtfISOS)) gtfISOS[[i]] <- 18
names(gtfISOS)<-sig_isos_gtf
gtfISO<-as.data.frame(gtfISOS, col.names=as.character("count"))
names(gtfISO) = row.names=as.character(sig_isos_gtf)
gtf_sigISOS<-t(gtfISO)

gtfCDS <- vector("list", length(sig_cds_gtf))
for (i in 1:length(gtfCDS)) gtfCDS[[i]] <- 18
names(gtfCDS)<-sig_cds_gtf
gtfCDS<-as.data.frame(gtfCDS, col.names=as.character("count"))
names(gtfCDS) = row.names=as.character(sig_cds_gtf)
gtf_sigCDS<-t(gtfCDS)

gtfGENES <- vector("list", length(sig_genes_gtf))
for (i in 1:length(gtfGENES)) gtfGENES[[i]] <- 18
names(gtfGENES)<-sig_genes_gtf
gtfGENES<-as.data.frame(gtfGENES, col.names=as.character("count"))
names(gtfGENES) = row.names=as.character(sig_genes_gtf)
gtf_sigGENES<-t(gtfGENES)

gtfTSS <- vector("list", length(sig_tss_gtf))
for (i in 1:length(gtfTSS)) gtfTSS[[i]] <- 18
names(gtfTSS)<-sig_tss_gtf
gtfTSS<-as.data.frame(gtfTSS, col.names=as.character("count"))
names(gtfTSS) = row.names=as.character(sig_tss_gtf)
gtf_sigTSS<-t(gtfTSS)

################################################################
# Leaveout_data significantly differentially expressed features
# 4 levels [genes, isoforms, cds, tss]
#################################################################

sigCDS<-read.table("sigcds_counts", col.names=c("counts", "names"))
sigTSS<-read.table("sigtss_counts", col.names=c("counts","names"))
sigISOS<-read.table("sigiso_counts", col.names=c("counts","names"))
sigGENES<-read.table("siggene_counts", col.names=c("counts","names"))

################################################################

topsigCDS<-order(sigCDS[,"counts"], decreasing=T)
sigCDS_ordered<-as.data.frame(sigCDS[,"names"][topsigCDS])
dim(sigCDS_ordered)
head(sigCDS_ordered)
topsigCDS<-as.list(sigCDS_ordered)

topsigGENES<-order(sigGENES[,"counts"], decreasing=T)
sigGENES_ordered<-as.data.frame(sigGENES[,"names"][topsigGENES])
dim(sigGENES_ordered)
head(sigGENES_ordered)
topsigGENES<-as.list(sigGENES_ordered)

topsigISOS<-order(sigISOS[,"counts"], decreasing=T)
sigISOS_ordered<-as.data.frame(sigISOS[,"names"][topsigISOS])
dim(sigISOS_ordered)
head(sigISOS_ordered)
topsigISOS<-as.list(sigISOS_ordered)

topsigTSS<-order(sigTSS[,"counts"], decreasing=T)
sigTSS_ordered<-as.data.frame(sigTSS[,"names"][topsigTSS])
dim(sigTSS_ordered)
head(sigTSS_ordered)
topsigTSS<-as.list(sigTSS_ordered)
###############################################################
pdf("leave-out_analsysis.pdf")
plot(density(sigISOS$count),ylab=c("Frequency"),xlab=c("trials"), main="Isoform distribution")
plot(sigISOS, xlab="# of trials which identified isoform signifcance",
	 ylab="number of significant isoforms",col=c("green", "yellow"),
	 main="Isoform significance consitency")
ISOS_sigcount<-factor(sigISOS[,"counts"],  labels=c("18","17","16","15","14","13","12","11","10","9","8","7","6","5","4","3","2","1"))
ISOS_names<-factor(sigISOS[,"names"])
plot(ISOS_sigcount, col=c("red", "green"),
	 xlab="# of trials which identified the isoforms as signifcant",
	 ylab="# of Isoforms identified as signifcant",
	 main="Isoforms significance consitency")

plot(density(sigCDS$count),ylab=c("Frequency"),xlab=c("trials"), main="CDS distribution")
plot(sigCDS, xlab="# of trials which identified CDS as signifcant",
	 ylab="number of significant CDS",col=c("green", "yellow"),
	 main="CDS significance consitency")
CDS_sigcount<-factor(sigCDS[,"counts"], labels=c("18","17","16","15","14","13","12","11","10","9","8","7","6","5","4","3","2","1"))
CDS_names<-factor(sigCDS[,"names"])
plot(CDS_sigcount, col=c("blue", "green"),
	 xlab="# of trials which identified the CDS as signifcant",
	 ylab="total number of CDS identified as significant",
	 main="CDS significance consitency")

plot(density(sigTSS$count),ylab=c("Frequency"),xlab="trials", main="TSS distribution")
plot(sigTSS, xlab="# of trials which identified tss signifcance",
	 ylab="total # of TSS_groups identified as significant",col=c("green", "red"),
	 main="TSS significance consitency")
TSS_sigcount<-factor(sigTSS[,"counts"] ,  labels=c("18","17","16","15","14","13","12","11","10","9","8","7","6","5","4","3","2","1"))
TSS_names<-factor(sigTSS[,"names"])
plot(TSS_sigcount, col=c("blue", "purple"),
	 xlab="# of TSS_group identified as signifcant",
	 ylab="total # of trials which identified the TSS_groups as significant",
	 main="TSS_group significance consitency")

plot(density(sigGENES$count),ylab=c("Frequency"),xlab="trials", main="Genes distribution")
plot(sigGENES, xlab="# of trials which identified the genes as signifcant",
	 ylab="total # of Genes identified as significant",col=c("darkgreen", "purple"),
	 main="Gene significance consitency")
GENES_sigcount<-factor(sigGENES[,"counts"], labels=c("18","17","16","15","14","13","12","11","10","9","8","7","6","5","4","3","2","1"))
GENES_names<-factor(sigGENES[,"names"])
plot(GENES_sigcount, col=c("orange", "green"),
	 xlab="# of trials which identified the genes as signifcant",
	 ylab="total # Genes identified as significant",
	 main="Genes significance consitency")

# VennDiagram
#################################################################

source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/vennDia.R")
#source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/R_BioCondManual.html#R_graphics_venn")
TSSs<-sigTSS[,"names"]
length(TSSs)
length(gtf_sigTSS)

ISOSs<-sigISOS[,"names"]
length(ISOSs)
length(gtf_sigISOS)

CDSs<-sigCDS[,"names"]
length(CDSs)
length(gtf_sigCDS)

GENESs<-sigGENES[,"names"]
length(GENESs)
length(gtf_sigGENES)

venGENE<-venn.diagram(list(A=GENESs, B=gtf_sigGENES), filename=NULL,
					  fill=rainbow(2), category.names=c("leaveout Genes", "gtf_only Genes"))
grid.newpage()
grid.draw(venGENE)
x<-grid.draw(venGENE)

leaveout_sigs <- venn.diagram(list(A=GENESs, B=CDSs, C=TSSs, D=ISOSs),
				   filename=NULL , fill=rainbow(4), category.names=c("sigGENES","sigCDS","sigTSS","sigISOS"))
grid.newpage()
grid.draw(leaveout_sigs)

Genes_tss<-GENESs[GENESs %in% TSSs]
Genes_tss_iso<-Genes_tss[Genes_tss %in% ISOSs]
leaveout_sigfeat<-Genes_tss_iso[Genes_tss_iso %in% CDSs]
length(leaveout_sigfeat)

leaveout_cts <- venn.diagram(list(A=topsigGENES[,"names"], B=topsigISOS[,"names"], C=topsigTSS[,"names"], D=topsigCDS[,"names"]),
				   filename=NULL , fill=rainbow(4),
				   category.names=c("sigGENES","sigCDS","sigTSS","sigISOS"))
grid.newpage()
grid.draw(leaveout_cts)

vennD=Venn(SetNames = c("sigGENES","sigCDS","sigTSS","sigISOS")) #

v1 <- venn.diagram(list(A=GENESs, B=CDSs, C=TSSs, D=ISOSs),
				   filename=NULL , fill=rainbow(4), category.names=)


vennD=Venn(SetNames = c("sigGENES","sigCDS","sigTSS","sigISOS")) #
v1 <- venn.diagram(list(A=GENESs, B=CDSs, C=TSSs, D=ISOSs),
				   filename=NULL , fill=rainbow(4), category.names=c("sigGENES","sigCDS","sigTSS","sigISOS"))
grid.newpage()
grid.draw(v1)

features.list<-list(A = CDSs, B = TSSs, C = ISOSs, D = GENESs)
venn.diagram(features.list,, filename="leavout_venn.pdf")
venn.diagram(features.list)
vennDiagram(features.list, include ="both", names=vennD)
## data
TSS <- data.frame(names = sample(topsigTSS), A = 1)
GENE <- data.frame(names = sample(topsigGENES), B = 1)
CDS <- data.frame(names = sample(topsigCDS), C = 1)
ISO <- data.frame(names = sample(topsigISOS), D = 1)

abcdgenes<-as.list(a=c(CDS),b=c(GENE), c=c(ISO), d=c(TSS))
colnames(B234) <- c("B2_genes","B3_genes","B4_genes")
VCB234 <- vennCounts(vennD, include="both")
vennDiagram(VCB234, circle.col=c("blue","red","green"), lwd=3)


# grab Thomas Girke's Venn Diagram functions

## Define venndiagram function
##  http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/R_BioCondManual.html#R_graphics_venn
vennDiagram(vennD, circle.col=c("blue","red","green"), lwd=3)

table(ISO)
table(GENE)
table(CDS)
table(TSS)
table(sort(CDS),sort(GENE), sort(ISO), sort(TSS))
sort(CDSs, decreasing=T)
B234<-as.data.frame(cbind(B2cell,B3cell, B4cell))
colnames(B234) <- c("B2_genes","B3_genes","B4_genes")
VCB234 <- vennCounts(vennD, include="both")
vennDiagram(VCB234, circle.col=c("blue","red","green"), lwd=3)



dev.off()

