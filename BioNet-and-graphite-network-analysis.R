### R code from vignette source 'graphite.Rnw'
library(graph)
library(graphite)
humanReactome <- pathways("hsapiens", "reactome")
names(humanReactome)[1:10]
p <- humanReactome[["ABC-family proteins mediated transport"]]
p


###################################################
### code chunk number 4: base2
###################################################
p <- humanReactome[[21]]
pathwayTitle(p)


###################################################
### code chunk number 5: base3
###################################################
head(nodes(p))


###################################################
### code chunk number 6: base4
###################################################
head(edges(p))


###################################################
### code chunk number 7: base5
###################################################
head(nodes(p), which = "mixed")
head(edges(p), which = "mixed")


###################################################
### code chunk number 8: base6
###################################################
pathwayDatabases()


###################################################
### code chunk number 9: graph1
###################################################
g <- pathwayGraph(p)
g


###################################################
### code chunk number 10: graph2
###################################################
edgeData(g)[1]


###################################################
### code chunk number 11: graph3
###################################################
g <- pathwayGraph(p, which = "mixed")
g


###################################################
### code chunk number 12: graph4
###################################################
edgeData(g)[1]


###################################################
### code chunk number 13: ident1
###################################################
pSymbol <- convertIdentifiers(p, "SYMBOL")
pSymbol
head(nodes(pSymbol))


###################################################
### code chunk number 14: ident2
###################################################
reactomeSymbol <- convertIdentifiers(humanReactome[1:5], "SYMBOL")


###################################################
### code chunk number 15: spia1
###################################################
library(SPIA)
data(colorectalcancer)

library(hgu133plus2.db)
top$ENTREZ <- mapIds(hgu133plus2.db,
                     top$ID, "ENTREZID", "PROBEID", multiVals = "first")
top <- top[!is.na(top$ENTREZ) & !duplicated(top$ENTREZ), ]
top$ENTREZ <- paste("ENTREZID", top$ENTREZ, sep = ":")
tg1 <- top[top$adj.P.Val < 0.05, ]

DE_Colorectal = tg1$logFC
names(DE_Colorectal) <- tg1$ENTREZ
ALL_Colorectal <- top$ENTREZ

biocarta <- pathways("hsapiens", "biocarta")[1:10]
biocarta <- convertIdentifiers(biocarta, "ENTREZID")
prepareSPIA(biocarta, "biocartaEx")
res <- runSPIA(de = DE_Colorectal, all = ALL_Colorectal, "biocartaEx")
res[1:5,]


###################################################
### code chunk number 16: tg1
###################################################
library(topologyGSA)
data(examples)
colnames(y1) <- paste("SYMBOL", colnames(y1), sep = ":")
colnames(y2) <- paste("SYMBOL", colnames(y2), sep = ":")

kegg <- pathways("hsapiens", "kegg")
p <- convertIdentifiers(kegg[["Fc epsilon RI signaling pathway"]], "SYMBOL")
runTopologyGSA(p, "var", y1, y2, 0.05)


###################################################
### code chunk number 17: tg2
###################################################
library(ALL)
library(a4Preproc)
library(clipper)

data(ALL)
pheno <- as(phenoData(ALL), "data.frame")
samples <- unlist(lapply(c("NEG", "BCR/ABL"), function(t) {
  which(grepl("^B\\d*", pheno$BT) & (pheno$mol.biol == t))[1:10]
}))
classes <- c(rep(1,10), rep(2,10))

expr <- exprs(ALL)[,samples]
rownames(expr) <- paste("ENTREZID",
                        featureData(addGeneInfo(ALL))$ENTREZID,
                        sep = ":")

k <- as.list(pathways("hsapiens", "kegg"))
selected <- k[c("Chronic myeloid leukemia",
                "Bladder cancer",
                "Cytosolic DNA-sensing pathway")]

clipped <- runClipper(selected, expr, classes, "mean", pathThr = 0.1)
resClip <- do.call(rbind, clipped$results)
resClip[, c("startIdx", "endIdx", "maxIdx", "lenght",
           "maxScore", "aScore", "involvedGenes")]


###################################################
### code chunk number 18: bldp
###################################################
edges <- data.frame(src_type = "ENTREZID", src="672",
                    dest_type = "ENTREZID", dest="7157",
                    direction="undirected", type="binding")
pathway <- buildPathway("#1", "example", "hsapiens", "database",
                        proteinEdges = edges)


###################################################
### code chunk number 19: bldp2
###################################################
edges <- data.frame(src_type = "ENTREZID", src="672",
                    dest_type = "ENTREZID", dest="7157",
                    direction="undirected", type="binding")
edgemix <- data.frame(src_type = "CHEBI", src="77750",
                    dest_type = "ENTREZID", dest="7157",
                    direction="undirected", type="biochemicalReaction")
edgemet <- data.frame(src_type = "CHEBI", src="15351",
                    dest_type = "CHEBI", dest="77750",
                    direction="undirected", type="biochemicalReaction")
pathway <- buildPathway("#1", "example", "hsapiens", "database",
                        proteinEdges = edges,
                        mixedEdges = edgemix,
                        metaboliteEdges = edgemet)


###################################################
### code chunk number 20: para1 (eval = FALSE)
###################################################
## options(Ncpus = 6)


###################################################
### code chunk number 21: para1 (eval = FALSE)
###################################################
## original <- pathways("hsapiens", "reactome")
##
## # Do (will exploit parallelism)
## converted <- convertIdentifiers(original, "SYMBOL")
##
## # Don't (no parallelism here)
## converted <- lapply(original, convertIdentifiers, "SYMBOL")


