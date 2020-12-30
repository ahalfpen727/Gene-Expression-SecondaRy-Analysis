## ----setup, message=FALSE, echo = FALSE---------------------------------------
library(limma);library(edgeR);library(EGSEAdata);library(EGSEA)
#biocM(c("EGSEAdata","EGSEA"))
## ----collectionlist-----------------------------------------------------------
egsea.data("human")
info = egsea.data("human", returnInfo = TRUE)
names(info)
info$msigdb$info$collections

## ----indexing-----------------------------------------------------------------
gs.annots = buildIdx(entrezIDs=v$genes$ENTREZID, species="human",
                     msigdb.gsets=c("c2", "c5"), go.part = TRUE)
names(gs.annots)


url = "http://bioinf.wehi.edu.au/EGSEA/mam.rnaseq.rdata"
utils::download.file(url, destfile="mam.rnaseq.rdata", mode="wb")
load("mam.rnaseq.rdata")
names(mam.rnaseq.data)
dim(mam.rnaseq.data)
x = calcNormFactors(mam.rnaseq.data, method = "TMM")
design = model.matrix(~0+x$samples$group+x$samples$lane)
colnames(design) = gsub("x\\$samples\\$group", "", colnames(design))
colnames(design) = gsub("x\\$samples\\$lane", "", colnames(design))
head(design)
contr.matrix = makeContrasts(
         BasalvsLP = Basal-LP,
         BasalvsML = Basal - ML,
         LPvsML = LP - ML,
         levels = colnames(design))
head(contr.matrix)

## ----setup2, eval=TRUE--------------------------------------------------------
v = voom(x, design, plot=FALSE)
names(v)


## ----exploresets--------------------------------------------------------------
class(gs.annots$c2)
summary(gs.annots$c2)
show(gs.annots$c2)
s = getSetByName(gs.annots$c2, "SMID_BREAST_CANCER_LUMINAL_A_DN")
class(s)
names(s)
names(s$SMID_BREAST_CANCER_LUMINAL_A_DN)

## ----indexclass---------------------------------------------------------------
slotNames(gs.annots$c2)

## ----symbolmap----------------------------------------------------------------
colnames(v$genes)
symbolsMap = v$genes[, c(1, 2)]
colnames(symbolsMap) = c("FeatureID", "Symbols")
symbolsMap[, "Symbols"] = as.character(symbolsMap[, "Symbols"])

## ----base---------------------------------------------------------------------
egsea.base()

## ----selectbasemethods--------------------------------------------------------
baseMethods = egsea.base()[-2]
baseMethods

## ----combine------------------------------------------------------------------
egsea.combine()

## ----sort---------------------------------------------------------------------
egsea.sort()

## ----egseatest----------------------------------------------------------------
gsa = egsea(voom.results=v, contrasts=contr.matrix,
         gs.annots=gs.annots, symbolsMap=symbolsMap,
         baseGSEAs=baseMethods, sort.by="med.rank",
         num.threads = 16, report = FALSE)

## ----showegsea----------------------------------------------------------------
show(gsa)

## ----summariseegsea-----------------------------------------------------------
summary(gsa)

## ----topsets------------------------------------------------------------------
topSets(gsa, gs.label="c2", contrast = "comparison", names.only=TRUE)

## ----topsetslim---------------------------------------------------------------
t = topSets(gsa, contrast = "comparison",
             names.only=FALSE, number = Inf, verbose = FALSE)
t[grep("LIM_", rownames(t)), c("p.adj", "Rank", "med.rank", "vote.rank")]

## ----topsets2-----------------------------------------------------------------
topSets(gsa, gs.label="kegg", contrast="BasalvsLP", sort.by="med.rank")
topSets(gsa, gs.label="kegg", contrast="comparison", sort.by="med.rank")

## ----heatmaps-----------------------------------------------------------------
plotHeatmap(gsa, gene.set="LIM_MAMMARY_STEM_CELL_UP", gs.label="c2",
         contrast = "comparison", file.name = "hm_cmp_LIM_MAMMARY_STEM_CELL_UP", format="png")
plotHeatmap(gsa, gene.set="LIM_MAMMARY_STEM_CELL_DN", gs.label="c2",
         contrast = "comparison", file.name = "hm_cmp_LIM_MAMMARY_STEM_CELL_DN", format="png")

## ----pathwayplot1, eval=FALSE-------------------------------------------------
#  plotPathway(gsa, gene.set = "Vascular smooth muscle contraction",
#               contrast = "BasalvsLP", gs.label = "kegg",
#               file.name = "Vascular_smooth_muscle_contraction")

## ----pathwayplot2, eval=FALSE-------------------------------------------------
#  plotPathway(gsa, gene.set = "Vascular smooth muscle contraction",
#               contrast = "comparison", gs.label = "kegg",
#               file.name = "Vascular_smooth_muscle_contraction_cmp")

## ----mdsplot------------------------------------------------------------------
plotMethods(gsa, gs.label = "c2", contrast = "BasalvsLP",
         file.name = "mds_c2_BasalvsLP", format="png")
plotMethods(gsa, gs.label = "c5BP", contrast = "BasalvsLP",
         file.name = "mds_c5_BasalvsLP", format="png")

## ----keggsummaryplot1---------------------------------------------------------
plotSummary(gsa, gs.label = 3, contrast = 3,
         file.name = "summary_kegg_LPvsML", format="png")

## ----c2summaryplot2-----------------------------------------------------------
plotSummary(gsa, gs.label = 1, contrast = 3,
         file.name = "summary_c2_LPvsML",
         x.axis = "med.rank", format="png")

## ----c2summaryplot3-----------------------------------------------------------
plotSummary(gsa, gs.label = 1, contrast = 3,
         file.name = "summary_sig_c2_LPvsML",
         x.axis = "med.rank", x.cutoff=300, format="png")

## ----summaryplotkegg1and2-----------------------------------------------------
plotSummary(gsa, gs.label = "kegg", contrast = c(1,2),
         file.name = "summary_kegg_1vs2", format="png")

## ----gographs-----------------------------------------------------------------
plotGOGraph(gsa, gs.label="c5BP", contrast = 1, file.name="BasalvsLP-c5BP-top-", format="png")
plotGOGraph(gsa, gs.label="c5CC", contrast = 1, file.name="BasalvsLP-c5CC-top-", format="png")

## ----summarybarplot-----------------------------------------------------------
plotBars(gsa, gs.label = "c2", contrast="comparison", file.name="comparison-c2-bars", format="png")

## ----summaryheatmap-----------------------------------------------------------
plotSummaryHeatmap(gsa, gs.label="c2", hm.vals = "avg.logfc.dir",
         file.name="summary_heatmaps_c2", format="png")
plotSummaryHeatmap(gsa, gs.label="kegg", hm.vals = "avg.logfc.dir",
         file.name="summary_heatmaps_kegg", format="png")

## ----toptable-----------------------------------------------------------------
t = limmaTopTable(gsa, contrast=1)
head(t)

## ----htmlreport, warning=FALSE, eval=FALSE------------------------------------
#  generateReport(gsa, number = 20, report.dir="./mam-rnaseq-egsea-report")

## ----mareadidats, eval=TRUE---------------------------------------------------
library(limma)
url = "http://bioinf.wehi.edu.au/EGSEA/arraydata.zip"
utils::download.file(url, destfile="arraydata.zip", mode="wb")
utils::unzip("arraydata.zip", exdir = ".")
targets = read.delim("targets.txt", header=TRUE, sep=" ")
data = read.idat(as.character(targets$File),
                   bgxfile="GPL6887_MouseWG-6_V2_0_R0_11278593_A.bgx",
                   annotation=c("Entrez_Gene_ID","Symbol", "Chromosome"))
data$other$Detection = detectionPValues(data)
data$targets = targets
colnames(data) = targets$Sample

## ----manormalize, eval=TRUE---------------------------------------------------
data = neqc(data)

## ----mafilter, eval=TRUE------------------------------------------------------
table(targets$Celltype)
keep.exprs = rowSums(data$other$Detection<0.05)>=5
table(keep.exprs)
data = data[keep.exprs,]
dim(data)
head(data$genes)
sum(is.na(data$genes$Entrez_Gene_ID))
data1 = data[!is.na(data$genes$Entrez_Gene_ID), ]
dim(data1)
ord = order(lmFit(data1)$Amean, decreasing=TRUE)
ids2keep = data1$genes$Array_Address_Id[ord][!duplicated(data1$genes$Entrez_Gene_ID[ord])]
data1 = data1[match(ids2keep, data1$genes$Array_Address_Id),]
dim(data1)

expr = data1$E
group = as.factor(data1$targets$Celltype)
probe.annot = data1$genes[, 2:4]
head(probe.annot)

## ----malinearmodel, eval=TRUE-------------------------------------------------
head(data1$targets)
experiment = as.character(data1$targets$Experiment)
design = model.matrix(~0 + group + experiment)
colnames(design) = gsub("group", "", colnames(design))
design
contr.matrix = makeContrasts(
         BasalvsLP = Basal-LP,
         BasalvsML = Basal-ML,
         LPvsML = LP-ML,
         levels = colnames(design))
contr.matrix

## ----maindex, eval=FALSE------------------------------------------------------
#  library(EGSEA)
#  library(EGSEAdata)
#  gs.annots = buildIdx(entrezIDs=probe.annot[, 2],
#               species="mouse",
#               msigdb.gsets=c("c2", "c5"), go.part = TRUE)
#  names(gs.annots)

## ----maegsea, eval=FALSE------------------------------------------------------
#  baseMethods = egsea.base()[-2]
#  baseMethods
#
#  gsam = egsea.ma(expr=expr, group=group,
#   		  	 probe.annot = probe.annot,
#   		  	 design = design,
#           contrasts=contr.matrix,
#           gs.annots=gs.annots,
#           baseGSEAs=baseMethods, sort.by="med.rank",
#           num.threads = 8, report = FALSE)

## ----mareport, warning=FALSE, eval=FALSE--------------------------------------
#  generateReport(gsam, number = 20, report.dir="./mam-ma-egsea-report")

## ----matopsets, eval=FALSE----------------------------------------------------
#  topSets(gsam, gs.label="c2", contrast="comparison", names.only=TRUE, number=5)

## ----softwareinfo-------------------------------------------------------------
sessionInfo()

