#####################################################################
# library(STRINGdb); KEGG_+_GO_analysis post-featurecounts.R
#####################################################################
source("https://bioconductor.org/biocLite.R")
# biocLite("limma");biocLite("Rsubread"); biocLite("STRINGdb")
browseVignettes("STRINGdb")
#STRINGdb$help("get_graph")
# get_interactions(string_ids)
###### returns the interactions in between the input proteins
# get_neighbors(string_ids)
###### Get the neighborhoods of a protein (or of a vector of proteins).
# get_subnetwork(string_ids)
###### returns a subgraph from the given input proteins

#####################################################################
# GO Pathway DOSE korg
# Sinorhizobium meliloti strain 1021
# [TAX:382] Taxonomy ID: 266834
# GO:0003674 - MF ,GO:0005575 - CC , GO:0008150 - BP
# 000006965* sme  Sinorhizobium meliloti 1021 --> NC_003047 11474104,11481430
#cuff<-readCufflinks(dbFile = "cuffData.db", genome = "1021_genome.fa", rebuild=T)

#####################################################################
### LibraryLoading
#####################################################################

library(GO.db);library(GOstats);library(topGO)
library(gage);library(pathview)
library(KEGG.db);library(KEGGgraph)
library(clusterProfiler);library(DOSE)
library(ReactomePA);library(STRINGdb)
library(igraph);library(biomaRt)
library(keggorthology)
#library("EnrichmentBrowser"); vignette("EnrichmentBrowser")
# library(org.Hs.eg.db); library(keggorthology);library(Path2PPI)

#####################################################################
# Standard GO analysis from edgeR
#####################################################################
library(STRINGdb);library(Rsubread)
library(limma); library(edgeR)
library(biomaRt)
# library(biomaRt) functions to create a genetable from a gff3
Gff2GeneTable("1021_genome.gff3")
load("geneTable.rda")
edb<-geneTable$GeneID
head(geneTable)

#########################################################################
# Grab DE tables from each comparison made in the featurecounts script
#########################################################################

head(A.vs.AB.DE)
dim(A.vs.AB.DE)

head(A.vs.wt1021.DE)
dim(A.vs.wt1021.DE)

head(AB.vs.wt1021.DE)
dim(AB.vs.wt1021.DE)

head(AB.vs.wt1021B.DE)
dim(AB.vs.wt1021B.DE)

head(wt1021.over.wt1021B.DE)
dim(wt1021B.vs.wt1021.DE)

#####################################################################
# library(STRINGdb); KEGG_IDS
#####################################################################

 browseVignettes("STRINGdb") ; #STRINGdb$help("get_graph")
## get_interactions(string_ids)   # returns the interactions in between the input proteins
## get_neighbors(string_ids)      # Get the neighborhoods of a protein (or of a vector of proteins).
## get_subnetwork(string_ids)     # returns a subgraph from the given input proteins
#000006965* sme  Sinorhizobium meliloti 1021 --> NC_003047 11474104,11481430

###########################################################################################
## Query STRINGdb database for species, get KEGGids, GOids, and STRINGids
###########################################################################################

sme1021 <- search_kegg_organism('Sinorhizobium meliloti 1021', by='scientific_name')
dim(sme1021); head(sme1021)
sme1021$kegg_code
Smeliloti <- search_kegg_organism('Sinorhizobium meliloti', by='scientific_name')
Smeliloti$scientific_name

smelil<-search_kegg_organism(sme1021$kegg_code, by='kegg_code')
dim(smelil);head(smelil)

species.all<-get_STRING_species(version="10", species_name=NULL)
colnames(species.all)
sm1021<-grep(pattern='Sinorhizobium meliloti', species.all$official_name, ignore.case = T)
taxa.info<-species.all[sm1021,]
taxa.info
taxID<-taxa.info$species_id
taxID
string.db.sme1021 <- STRINGdb$new(version="10", species=taxID)
string.db.sme1021

sme.kegg1021<-search_kegg_organism('sme', by='kegg_code')
sme.kegg.org1021<- search_kegg_organism('Sinorhizobium meliloti 1021', by='scientific_name')
dim(sme.kegg.org1021)
head(sme.kegg.org1021)
sme.pwys <- download.kegg.pathways("sme")
kegg.gs <- get.kegg.genesets("sme")
head(sme.pwys)
head(kegg.gs)

library(gage)
data(gse16873)
sme.kegg.sets<-kegg.gsets(species = "sme",id.type = "kegg")
sme.kegg.sets


###########################################################################################
## Write the KEGGids, the genes involved, and the human readable pathway names to a file
###########################################################################################

keggfile<-file.path("Sme1021.kegg.genesets.txt", "w")
keggfile<-file("Sme1021.kegg.genesets.txt", "a")
KEGGid<-names(kegg.gs)
x=0
for (keggpath in kegg.gs){
  x<-c(x + 1)
  kegg.df <-c(x,KEGGid[x],keggpath)
  write(kegg.df, file=keggfile, append=T)
}

###########################################################################################
## For each comparison DE table, map the gene symbols to the KEGGids/STRINGids
###########################################################################################


A.vs.wt1021.DE_mapped <- string.db.sme1021$map( A.vs.wt1021.DE, "GeneSymbol", removeUnmappedRows = TRUE )
write.table(A.vs.wt1021.DE_mapped, file="A.vs.wt1021.KEGG.difftable")
head(A.vs.wt1021.DE_mapped)
dim(A.vs.wt1021.DE_mapped)

AB.vs.wt1021B.DE_mapped <- string.db.sme1021$map( AB.vs.wt1021B.DE, "GeneSymbol", removeUnmappedRows = TRUE )
write.table(AB.vs.wt1021B.DE_mapped, file="AB.vs.wt1021B.KEGG.difftable")
head(AB.vs.wt1021B.DE_mapped)

AB.vs.wt1021.DE_mapped <- string.db.sme1021$map( AB.vs.wt1021.DE, "GeneSymbol", removeUnmappedRows = TRUE )
write.table(AB.vs.wt1021.DE_mapped, file="AB.vs.wt1021.KEGG.difftable")
head(AB.vs.wt1021.DE_mapped)

A.vs.AB.DE_mapped <- string.db.sme1021$map( A.vs.AB.DE, "GeneSymbol", removeUnmappedRows = TRUE )
write.table(A.vs.AB.DE_mapped, file="A.vs.AB.KEGG.difftable")
head(A.vs.AB.DE_mapped)

wt1021.vs.wt1021B.DE_mapped <- string.db.sme1021$map( wt1021.vs.wt1021B.DE, "GeneSymbol", removeUnmappedRows = TRUE )
write.table(wt1021.vs.wt1021B.DE_mapped, file="wt1021.vs.wt1021B.KEGG.difftable")
head(wt1021.vs.wt1021B.DE_mapped)

#####################################################################
# enrichment
###################################################

AB.vs.wt1021B.DE.df<-as.data.frame(cbind(gene=AB.vs.wt1021B.DE_mapped$GeneSymbol,
                              pvalue=AB.vs.wt1021B.DE_mapped$Pval,
                              logFC=AB.vs.wt1021B.DE_mapped$LogFoldChange), stringsAsFactors=F)
dim(AB.vs.wt1021B.DE.df)
head(AB.vs.wt1021B.DE.df)

AB.vs.wt1021B.intersected<-string.db.sme1021$map(AB.vs.wt1021B.DE.df, "gene", removeUnmappedRows=T)
head(AB.vs.wt1021B.intersected)
class(AB.vs.wt1021B.intersected)
string.db.sme1021$plot_network(AB.vs.wt1021B.intersected$STRING_id[1:400],)

AB.vs.wt1021B.subnets<-string.db.sme1021$get_subnetwork(AB.vs.wt1021B.intersected)
AB.vs.wt1021B.subnets

AB.vs.wt1021B.mapped_sig<-as.data.frame(cbind(genes=c(AB.vs.wt1021B.intersected$gene[AB.vs.wt1021B.intersected$pvalue < 0.05]),
                                              pvalue=c(AB.vs.wt1021B.intersected$pvalue[AB.vs.wt1021B.intersected$pvalue < 0.05]),
                                        logFC=c(AB.vs.wt1021B.intersected$logFC[AB.vs.wt1021B.intersected$pvalue < 0.05]),
                                        STRING_id=c(AB.vs.wt1021B.intersected$STRING_id[AB.vs.wt1021B.intersected$pvalue < 0.05])),
                                        stringsAsFactors=F, row.names=F)
head(AB.vs.wt1021B.mapped_sig)

AB.vs.wt1021B.DE.pv.fc.STRING<-as.data.frame(cbind(gene=AB.vs.wt1021B.DE_mapped$GeneSymbol,
                              pvalue=AB.vs.wt1021B.DE_mapped$Pval,
                              logFC=AB.vs.wt1021B.DE_mapped$LogFoldChange,
                              STRING_id=AB.vs.wt1021B.DE_mapped$STRING_id), stringsAsFactors=F, row.names=F, col.names=T)
head(AB.vs.wt1021B.DE.pv.fc.STRING)

# post payload information to the STRING server
AB.vs.wt1021B_pval01 <- string.db.sme1021$post_payload(AB.vs.wt1021B.mapped_sig$STRING_id,
                                         colors=AB.vs.wt1021B.mapped_sig["pvalue"]$color )

# display a STRING network png with the "halo"

string.db.sme1021$plot_network( AB.vs.wt1021B.DE.pv.fc.STRING$STRING_id[1:50],
                                payload_id=AB.vs.wt1021B_pval01,
                                required_score=AB.vs.wt1021B.DE.pv.fc.STRING$logFC[1:50])

# plot the enrichment for the best 100 genes
ab.wt1021.top100<-string.db.sme1021$plot_ppi_enrichment( AB.vs.wt1021B.intersected$STRING_id[1:500], quiet=TRUE )

#####################################################################
# enrichment A.vs.wt1021
###################################################

A.vs.wt1021.DE.df<-as.data.frame(cbind(gene=A.vs.wt1021.DE_mapped$GeneSymbol,
                              pvalue=A.vs.wt1021.DE_mapped$Pval,
                              logFC=A.vs.wt1021.DE_mapped$LogFoldChange), stringsAsFactors=F)
dim(A.vs.wt1021.DE.df)
head(A.vs.wt1021.DE.df)

A.vs.wt1021.intersected<-string.db.sme1021$map(A.vs.wt1021.DE.df, "gene", removeUnmappedRows=T)
head(A.vs.wt1021.intersected)
class(A.vs.wt1021.intersected)
string.db.sme1021$plot_network(A.vs.wt1021.intersected$STRING_id[1:400],)

A.vs.wt1021.subnets<-string.db.sme1021$get_subnetwork(A.vs.wt1021.intersected)
A.vs.wt1021.subnets

A.vs.wt1021.mapped_sig<-as.data.frame(cbind(genes=c(A.vs.wt1021.intersected$gene[A.vs.wt1021.intersected$pvalue < 0.05]),
                                              pvalue=c(A.vs.wt1021.intersected$pvalue[A.vs.wt1021.intersected$pvalue < 0.05]),
                                        logFC=c(A.vs.wt1021.intersected$logFC[A.vs.wt1021.intersected$pvalue < 0.05]),
                                        STRING_id=c(A.vs.wt1021.intersected$STRING_id[A.vs.wt1021.intersected$pvalue < 0.05])),
                                        stringsAsFactors=F, row.names=F)
head(A.vs.wt1021.mapped_sig)

A.vs.wt1021.DE.pv.fc.STRING<-as.data.frame(cbind(gene=A.vs.wt1021.DE_mapped$GeneSymbol,
                              pvalue=A.vs.wt1021.DE_mapped$Pval,
                              logFC=A.vs.wt1021.DE_mapped$LogFoldChange,
                              STRING_id=A.vs.wt1021.DE_mapped$STRING_id), stringsAsFactors=F, row.names=F, col.names=T)
head(A.vs.wt1021.DE.pv.fc.STRING)
# post payload information to the STRING server
A.vs.wt1021_pval01 <- string.db.sme1021$post_payload(A.vs.wt1021.mapped_sig$STRING_id,
                                         colors=A.vs.wt1021.mapped_sig["pvalue"]$color )

# display a STRING network png with the "halo"

string.db.sme1021$plot_network( A.vs.wt1021.DE.pv.fc.STRING$STRING_id[1:50],
                                payload_id=A.vs.wt1021_pval01,
                                required_score=A.vs.wt1021.DE.pv.fc.STRING$logFC[1:50])

# plot the enrichment for the best 100 genes
ab.wt1021.top100<-string.db.sme1021$plot_ppi_enrichment( A.vs.wt1021.intersected$STRING_id[1:500], quiet=TRUE )

#####################################################################
# enrichment
###################################################

AB.vs.wt1021.DE.df<-as.data.frame(cbind(gene=AB.vs.wt1021.DE_mapped$GeneSymbol,
                              pvalue=AB.vs.wt1021.DE_mapped$Pval,
                              logFC=AB.vs.wt1021.DE_mapped$LogFoldChange), stringsAsFactors=F)
dim(AB.vs.wt1021.DE.df)
head(AB.vs.wt1021.DE.df)

AB.vs.wt1021.intersected<-string.db.sme1021$map(AB.vs.wt1021.DE.df, "gene", removeUnmappedRows=T)
head(AB.vs.wt1021.intersected)
class(AB.vs.wt1021.intersected)
string.db.sme1021$plot_network(AB.vs.wt1021.intersected$STRING_id[1:400],)

AB.vs.wt1021.subnets<-string.db.sme1021$get_subnetwork(AB.vs.wt1021.intersected)
AB.vs.wt1021.subnets

AB.vs.wt1021.mapped_sig<-as.data.frame(cbind(genes=c(AB.vs.wt1021.intersected$gene[AB.vs.wt1021.intersected$pvalue < 0.05]),
                                              pvalue=c(AB.vs.wt1021.intersected$pvalue[AB.vs.wt1021.intersected$pvalue < 0.05]),
                                        logFC=c(AB.vs.wt1021.intersected$logFC[AB.vs.wt1021.intersected$pvalue < 0.05]),
                                        STRING_id=c(AB.vs.wt1021.intersected$STRING_id[AB.vs.wt1021.intersected$pvalue < 0.05])),
                                        stringsAsFactors=F, row.names=F)
head(AB.vs.wt1021.mapped_sig)

AB.vs.wt1021.DE.pv.fc.STRING<-as.data.frame(cbind(gene=AB.vs.wt1021.DE_mapped$GeneSymbol,
                              pvalue=AB.vs.wt1021.DE_mapped$Pval,
                              logFC=AB.vs.wt1021.DE_mapped$LogFoldChange,
                              STRING_id=AB.vs.wt1021.DE_mapped$STRING_id), stringsAsFactors=F, row.names=F, col.names=T)
head(AB.vs.wt1021.DE.pv.fc.STRING)
# post payload information to the STRING server
AB.vs.wt1021_pval01 <- string.db.sme1021$post_payload(AB.vs.wt1021.mapped_sig$STRING_id,
                                         colors=AB.vs.wt1021.mapped_sig["pvalue"]$color )

# display a STRING network png with the "halo"

string.db.sme1021$plot_network( AB.vs.wt1021.DE.pv.fc.STRING$STRING_id[1:50],
                                payload_id=AB.vs.wt1021_pval01,
                                required_score=AB.vs.wt1021.DE.pv.fc.STRING$logFC[1:50])

# plot the enrichment for the best 100 genes
ab.wt1021.top100<-string.db.sme1021$plot_ppi_enrichment( AB.vs.wt1021.intersected$STRING_id[1:500], quiet=TRUE )

#############################################################################
#####################################################################
# enrichment
###################################################

A.vs.AB.DE.df<-as.data.frame(cbind(gene=A.vs.AB.DE_mapped$GeneSymbol,
                              pvalue=A.vs.AB.DE_mapped$Pval,
                              logFC=A.vs.AB.DE_mapped$LogFoldChange), stringsAsFactors=F)
dim(A.vs.AB.DE.df)
head(A.vs.AB.DE.df)

A.vs.AB.intersected<-string.db.sme1021$map(A.vs.AB.DE.df, "gene", removeUnmappedRows=T)
head(A.vs.AB.intersected)
class(A.vs.AB.intersected)
string.db.sme1021$plot_network(A.vs.AB.intersected$STRING_id[1:400],)

A.vs.AB.subnets<-string.db.sme1021$get_subnetwork(A.vs.AB.intersected)
A.vs.AB.subnets

A.vs.AB.mapped_sig<-as.data.frame(cbind(genes=c(A.vs.AB.intersected$gene[A.vs.AB.intersected$pvalue < 0.05]),
                                              pvalue=c(A.vs.AB.intersected$pvalue[A.vs.AB.intersected$pvalue < 0.05]),
                                        logFC=c(A.vs.AB.intersected$logFC[A.vs.AB.intersected$pvalue < 0.05]),
                                        STRING_id=c(A.vs.AB.intersected$STRING_id[A.vs.AB.intersected$pvalue < 0.05])),
                                        stringsAsFactors=F, row.names=F)
head(A.vs.AB.mapped_sig)

A.vs.AB.DE.pv.fc.STRING<-as.data.frame(cbind(gene=A.vs.AB.DE_mapped$GeneSymbol,
                              pvalue=A.vs.AB.DE_mapped$Pval,
                              logFC=A.vs.AB.DE_mapped$LogFoldChange,
                              STRING_id=A.vs.AB.DE_mapped$STRING_id), stringsAsFactors=F, row.names=F, col.names=T)
head(A.vs.AB.DE.pv.fc.STRING)

# post payload information to the STRING server
A.vs.AB_pval01 <- string.db.sme1021$post_payload(A.vs.AB.mapped_sig$STRING_id,
                                         colors=A.vs.AB.mapped_sig["pvalue"]$color )

# display a STRING network png with the "halo"

string.db.sme1021$plot_network( A.vs.AB.DE.pv.fc.STRING$STRING_id[1:50],
                                payload_id=A.vs.AB_pval01,
                                required_score=A.vs.AB.DE.pv.fc.STRING$logFC[1:50])

# plot the enrichment for the best 100 genes
ab.wt1021.top100<-string.db.sme1021$plot_ppi_enrichment( A.vs.AB.intersected$STRING_id[1:500], quiet=TRUE )

#############################################################################
#####################################################################
# enrichment
###################################################

wt1021.vs.wt1021B.DE.df<-as.data.frame(cbind(gene=wt1021.vs.wt1021B.DE_mapped$GeneSymbol,
                              pvalue=wt1021.vs.wt1021B.DE_mapped$Pval,
                              logFC=wt1021.vs.wt1021B.DE_mapped$LogFoldChange), stringsAsFactors=F)
dim(wt1021.vs.wt1021B.DE.df)
head(wt1021.vs.wt1021B.DE.df)

wt1021.vs.wt1021B.intersected<-string.db.sme1021$map(wt1021.vs.wt1021B.DE.df, "gene", removeUnmappedRows=T)
head(wt1021.vs.wt1021B.intersected)
class(wt1021.vs.wt1021B.intersected)
string.db.sme1021$plot_network(wt1021.vs.wt1021B.intersected$STRING_id[1:400],)

wt1021.vs.wt1021B.subnets<-string.db.sme1021$get_subnetwork(wt1021.vs.wt1021B.intersected)
wt1021.vs.wt1021B.subnets

wt1021.vs.wt1021B.mapped_sig<-as.data.frame(cbind(genes=c(wt1021.vs.wt1021B.intersected$gene[wt1021.vs.wt1021B.intersected$pvalue < 0.05]),
                                              pvalue=c(wt1021.vs.wt1021B.intersected$pvalue[wt1021.vs.wt1021B.intersected$pvalue < 0.05]),
                                        logFC=c(wt1021.vs.wt1021B.intersected$logFC[wt1021.vs.wt1021B.intersected$pvalue < 0.05]),
                                        STRING_id=c(wt1021.vs.wt1021B.intersected$STRING_id[wt1021.vs.wt1021B.intersected$pvalue < 0.05])),
                                        stringsAsFactors=F, row.names=F)
head(wt1021.vs.wt1021B.mapped_sig)

wt1021.vs.wt1021B.DE.pv.fc.STRING<-as.data.frame(cbind(gene=wt1021.vs.wt1021B.DE_mapped$GeneSymbol,
                              pvalue=wt1021.vs.wt1021B.DE_mapped$Pval,
                              logFC=wt1021.vs.wt1021B.DE_mapped$LogFoldChange,
                              STRING_id=wt1021.vs.wt1021B.DE_mapped$STRING_id), stringsAsFactors=F, row.names=F, col.names=T)
head(wt1021.vs.wt1021B.DE.pv.fc.STRING)
# post payload information to the STRING server
wt1021.vs.wt1021B_pval01 <- string.db.sme1021$post_payload(wt1021.vs.wt1021B.mapped_sig$STRING_id,
                                         colors=wt1021.vs.wt1021B.mapped_sig["pvalue"]$color )

# display a STRING network png with the "halo"

string.db.sme1021$plot_network( wt1021.vs.wt1021B.DE.pv.fc.STRING$STRING_id[1:50],
                                payload_id=wt1021.vs.wt1021B_pval01,
                                required_score=wt1021.vs.wt1021B.DE.pv.fc.STRING$logFC[1:50])

# plot the enrichment for the best 100 genes
wt1021.wt1021b.top100<-string.db.sme1021$plot_ppi_enrichment( wt1021.vs.wt1021B.intersected$STRING_id[1:500], quiet=TRUE )


## ----eval = FALSE--------------------------------------------------------
wt1021.vs.wt1021B.kegg.rich <- enrichKEGG(gene = wt1021.vs.wt1021B.DE_mapped$GeneSymbol,organism='sme',pvalueCutoff = 0.05)
head(wt1021.vs.wt1021B.kegg.rich)
dim(wt1021.vs.wt1021B.kegg.rich)
wt1021.vs.wt1021B.mkegg.rich <- enrichMKEGG(gene = wt1021.vs.wt1021B.DE_mapped$GeneSymbol,organism='sme',pvalueCutoff = 0.05)
head(wt1021.vs.wt1021B.mkegg.rich)
dim(wt1021.vs.wt1021B.mkegg.rich)

barplot(wt1021.vs.wt1021B.mkegg.rich, drop=TRUE, showCategory=12)
barplot(wt1021.vs.wt1021B.kegg.rich,drop=T, showCategory=12)
dotplot(wt1021.vs.wt1021B.mkegg.rich)
dotplot(wt1021.vs.wt1021B.kegg)
cnetplot(wt1021.vs.wt1021B.mkegg.rich, categorySize="pvalue") # ,wt1021.vs.wt1021B.kegg
enrichMap(wt1021.vs.wt1021B.mkegg.rich)
cnetplot(wt1021.vs.wt1021B.kegg.rich, categorySize="pvalue") # ,wt1021.vs.wt1021B.kegg
enrichMap(wt1021.vs.wt1021B.kegg.rich)
cnetplot(wt1021.vs.wt1021B.kegg.rich,categorySize="pvalue", foldChange=,wt1021.vs.wt1021B.DE_mapped$LogFoldChange, )
cnetplot(wt1021.vs.wt1021B.kegg.rich, categorySize="pvalue") # ,wt1021.vs.wt1021B.kegg
enrichMap(wt1021.vs.wt1021B.kegg.rich)


## ----fig.height=12, fig.width=8------------------------------------------
AB.vs.wt1021.mkegg <- enrichMKEGG(gene = AB.vs.wt1021.DE_mapped$GeneSymbol,
                    organism = 'sme')
AB.vs.wt1021.kegg <- enrichKEGG(gene = AB.vs.wt1021.DE_mapped$GeneSymbol,
                    organism = 'sme')

barplot(AB.vs.wt1021.mkegg, drop=TRUE, showCategory=12)
barplot(AB.vs.wt1021.kegg, showCategory=8)
dotplot(AB.vs.wt1021.mkegg)
dotplot(AB.vs.wt1021.kegg)
cnetplot(AB.vs.wt1021.mkegg, categorySize="pvalue") # ,wt1021.vs.wt1021B.kegg
enrichMap(AB.vs.wt1021.mkegg)
cnetplot(AB.vs.wt1021.kegg, categorySize="pvalue",foldChange=AB.vs.wt1021.DE_mapped$LogFoldChange)
enrichMap(AB.vs.wt1021.kegg)
# GO analysis adjusting for gene length bias
# (assuming that y$genes$Length contains gene lengths)
library(EnrichmentBrowser)
go.abund <- goa(A.vs.wt1021.DE, geneid = "GeneID", trend = T)
go.abund

go.len <- goanna(A.vs.wt1021.DE, geneid = "GeneID", trend = "Length")

topGO(go.len, sort = "Qval")
topGO(go.len, sort = "Pval")
#Avswt1021.up_reg_genes<-topGO(go.abund, sort = "Qval")
#tAvswt1021.down_reg_genes<-topGO(go.abund, sort = "Qval")


## Default usage with a list of gene sets:

go.de <- goana(list(DE1 = EG.DE1, DE2 = EG.DE2, DE3 = EG.DE3))
topGO(go.de, sort = "DE1")
topGO(go.de, sort = "DE2")
topGO(go.abund, ontology = "BP")
topGO(go.de, ontology = "CC", sort = "DE3")
topGO(go.de, ontology = "MF", sort = "DE3")

## Standard KEGG analysis

AB.vs.wt1021.DE.kegg <- kegga(AB.vs.wt1021.DE$GeneSymbol, species.KEGG="sme") # equivalent to previous
AB.vs.wt1021.DE.kegg

barplot(AB.vs.wt1021.DE.kegg$DE,drop=T, showCategory=8)
## ------------------------------------------------------------------------
dotplot(AB.vs.wt1021.DE.kegg)
## ----fig.cap="enrichment map of enrichment result", fig.align="center", fig.height=16, fig.width=16, eval=FALSE----
 enrichMap(mkk)
## ## categorySize can be scaled by 'pvalue' or 'geneNum'
 cnetplot(AB.vs.wt1021.DE_mapped$STRING_id, categorySize="pvalue", foldChange=AB.vs.wt1021.DE_mapped$LogFoldChange)

## ----fig.height=12, fig.width=8------------------------------------------
plotGOgraph(ego)

## ----fig.cap="plotting gsea result", fig.align="center", fig.height=6, fig.width=8----
gseaplot(AB.vs.wt1021.DE.kegg, geneSetID = "sme")

head(ggo)
AB.vs.wt1021.DE_gse<-as.data.frame(AB.vs.wt1021.DE_mapped$Pval, row.names=c(AB.vs.wt1021.DE_mapped$GeneSymbol), stringsAsFactors=F)
names(AB.vs.wt1021.DE_sorted.gse)<-AB.vs.wt1021.DE_mapped$GeneSymbol
AB.vs.wt1021.DE_sorted.gse<-sort(AB.vs.wt1021.DE.kegg$P.DE, decreasing=T)
kk2 <- gseKEGG(geneList     = AB.vs.wt1021.DE_sorted.gse,
               organism     = 'sme',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)
 ##############################################################################
### DOSE/ClusterProfile - KEGG and Gene Ontology analysis
###############################################################################
wt1021.vs.wt1021B.mkegg.rich <- enrichMKEGG(gene = wt1021.vs.wt1021B.DE_mapped$GeneSymbol,organism = 'sme')
wt1021.vs.wt1021B.kegg.rich <- enrichKEGG(gene = wt1021.vs.wt1021B.DE_mapped$GeneSymbol,organism='sme',pvalueCutoff = 0.05)
head(wt1021.vs.wt1021B.kegg.rich)
dim(wt1021.vs.wt1021B.kegg.rich)
wt1021.vs.wt1021B.DE.kegg.rich <- enrichMKEGG(gene = wt1021.vs.wt1021B.DE_mapped$GeneSymbol,organism='sme')
head(wt1021.vs.wt1021B.DE.kegg.rich)
dim(wt1021.vs.wt1021B.DE.kegg.rich)


AB.vs.wt1021.kegg.rich <- enrichMKEGG(gene = AB.vs.wt1021.DE_mapped$GeneSymbol,organism='sme',pvalueCutoff = 0.05)
head(AB.vs.wt1021.kegg.rich)
dim(AB.vs.wt1021.kegg.rich)
AB.vs.wt1021B.kegg.rich <- enrichKEGG(gene = AB.vs.wt1021B.DE_mapped$GeneSymbol,organism='sme',pvalueCutoff = 0.05)
head(AB.vs.wt1021B.kegg.rich)
dim(AB.vs.wt1021B.kegg.rich)
AB.vs.wt1021B.kegg.rich
A.vs.wt1021.kegg.rich <- enrichKEGG(gene = A.vs.wt1021.DE_mapped$GeneSymbol,organism='sme',pvalueCutoff = 0.05)
head(A.vs.wt1021.kegg.rich)
dim(A.vs.wt1021.kegg.rich)
A.vs.wt1021.kegg.rich
## ----KEGG Download------------------------------------------------------------
sme.kegg.code<-search_kegg_organism('sme', by='kegg_code')
sme.kegg.code
go.abund
Smeliloti.kegg <- search_kegg_organism('Sinorhizobium meliloti 1021', by='scientific_name')
dim(Smeliloti.kegg)
head(Smeliloti.kegg)
sme.1021.kegg<-download_KEGG(species="sme", keggType = "KEGG", keyType = "kegg")
sme.1021.kegg$KEGGPATHID2EXTID[1:10,1]
sme.1021.kegg$KEGGPATHID2EXTID[1:10,2]
length(sme.1021.kegg)
names(sme.1021.kegg)
bitr_kegg
sme.1021.kegg
gene.df <- bitr(names(AB.vs.wt1021.DE_mapped$GeneSymbol), fromType = "SYMBOL",toType = c("ENTREZID", "KEGG"),OrgDb=sme.1021.kegg)


#######################################################################
##
#######################################################################

A.vs.wt1021.DE_mapped <- string.db.sme1021$map( A.vs.wt1021.DE, "GeneSymbol", removeUnmappedRows = TRUE )
write.table(A.vs.wt1021.DE_mapped, file="A.vs.wt1021.KEGG.difftable")
head(A.vs.wt1021.DE_mapped)

AB.vs.wt1021B.DE_mapped <- string.db.sme1021$map( AB.vs.wt1021B.DE, "GeneSymbol", removeUnmappedRows = TRUE )
write.table(AB.vs.wt1021B.DE_mapped, file="AB.vs.wt1021B.KEGG.difftable")
head(AB.vs.wt1021B.DE_mapped)

AB.vs.wt1021.DE_mapped <- string.db.sme1021$map( AB.vs.wt1021.DE, "GeneSymbol", removeUnmappedRows = TRUE )
write.table(AB.vs.wt1021.DE_mapped, file="AB.vs.wt1021.KEGG.difftable")
head(AB.vs.wt1021.DE_mapped)

A.vs.AB.DE_mapped <- string.db.sme1021$map( A.vs.AB.DE, "GeneSymbol", removeUnmappedRows = TRUE )
write.table(A.vs.AB.DE_mapped, file="A.vs.AB.KEGG.difftable")
head(A.vs.AB.DE_mapped)

wt1021B.vs.wt1021.DE_mapped <- string.db.sme1021$map( wt1021B.vs.wt1021.DE, "GeneSymbol", removeUnmappedRows = TRUE )
write.table(wt1021B.vs.wt1021.DE_mapped, file="wt1021B.vs.wt1021.KEGG.difftable")
head(tw1021B.vs.wt1021.DE_mapped)

# GO analysis adjusting for gene length bias
# (assuming that y$genes$Length contains gene lengths)
library(EnrichmentBrowser)
go.abund <- goa(A.vs.wt1021.DE, geneid = "GeneID", trend = T)
go.abund

go.len <- goanna(A.vs.wt1021.DE, geneid = "GeneID", trend = "Length")

topGO(go.len, sort = "Qval")
topGO(go.len, sort = "Pval")
#Avswt1021.up_reg_genes<-topGO(go.abund, sort = "Qval")
#tAvswt1021.down_reg_genes<-topGO(go.abund, sort = "Qval")


## Default usage with a list of gene sets:

go.de <- goana(list(DE1 = EG.DE1, DE2 = EG.DE2, DE3 = EG.DE3))
topGO(go.de, sort = "DE1")
topGO(go.de, sort = "DE2")
topGO(go.abund, ontology = "BP")
topGO(go.de, ontology = "CC", sort = "DE3")
topGO(go.de, ontology = "MF", sort = "DE3")

## Standard KEGG analysis

AB.vs.wt1021.DE.kegg <- kegga(AB.vs.wt1021.DE$GeneSymbol, species.KEGG="sme") # equivalent to previous
AB.vs.wt1021.DE.kegg

barplot(AB.vs.wt1021.DE.kegg$DE,drop=T, showCategory=8)
## ------------------------------------------------------------------------
dotplot(AB.vs.wt1021.DE.kegg)
## ----fig.cap="enrichment map of enrichment result", fig.align="center", fig.height=16, fig.width=16, eval=FALSE----
 enrichMap(mkk)
## ## categorySize can be scaled by 'pvalue' or 'geneNum'
 cnetplot(AB.vs.wt1021.DE_mapped$STRING_id, categorySize="pvalue", foldChange=AB.vs.wt1021.DE_mapped$LogFoldChange)

## ----fig.height=12, fig.width=8------------------------------------------
plotGOgraph(ego)

## ----fig.cap="plotting gsea result", fig.align="center", fig.height=6, fig.width=8----
gseaplot(AB.vs.wt1021.DE.kegg, geneSetID = "sme")

head(ggo)
AB.vs.wt1021.DE_gse<-as.data.frame(AB.vs.wt1021.DE_mapped$Pval, row.names=c(AB.vs.wt1021.DE_mapped$GeneSymbol), stringsAsFactors=F)
names(AB.vs.wt1021.DE_sorted.gse)<-AB.vs.wt1021.DE_mapped$GeneSymbol
AB.vs.wt1021.DE_sorted.gse<-sort(AB.vs.wt1021.DE.kegg$P.DE, decreasing=T)
kk2 <- gseKEGG(geneList     = AB.vs.wt1021.DE_sorted.gse,
               organism     = 'sme',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)
## ------------------------------------------------------------------------
ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)


## ------------------------------------------------------------------------
kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kk)

## ------------------------------------------------------------------------
kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)

## ----eval = FALSE--------------------------------------------------------
## mkk <- enrichMKEGG(gene = gene,
##                    organism = 'hsa')

## ----eval=FALSE----------------------------------------------------------
## mkk2 <- gseMKEGG(geneList = geneList,
##                  species = 'hsa')

## ----eval=FALSE----------------------------------------------------------
## david <- enrichDAVID(gene = gene,
##                      idType = "ENTREZ_GENE_ID",
##                      listType = "Gene",
##                      annotation = "KEGG_PATHWAY",
##                      david.user = "clusterProfiler@hku.hk")

## ----fig.height=5, fig.width=9-------------------------------------------
barplot(ggo, drop=TRUE, showCategory=12)

## ----fig.height=5, fig.width=8-------------------------------------------
barplot(mkk, showCategory=8)
## ------------------------------------------------------------------------
dotplot(mkk)
## ----fig.cap="enrichment map of enrichment result", fig.align="center", fig.height=16, fig.width=16, eval=FALSE----
 enrichMap(mkk)
## ## categorySize can be scaled by 'pvalue' or 'geneNum'
 cnetplot(mkk, categorySize="pvalue", foldChange=geneList)

## ----fig.height=12, fig.width=8------------------------------------------
plotGOgraph(ego)

## ----fig.cap="plotting gsea result", fig.align="center", fig.height=6, fig.width=8----
gseaplot(kk, geneSetID = "sme")


## ----eval=FALSE----------------------------------------------------------
 library("pathview")
hsa04110 <- pathview(gene.data  = geneList,
                     species    = "sme",
                     limit      = list(gene=max(abs(geneList)), cpd=1))

## ------------------------------------------------------------------------
data(gcSample)
lapply(gcSample, head)

## ------------------------------------------------------------------------
ck <- compareCluster(geneCluster = gcSample, fun = "enrichKEGG")
head(as.data.frame(ck))

## ------------------------------------------------------------------------
mydf <- data.frame(Entrez=names(geneList), FC=geneList)
mydf <- mydf[abs(mydf$FC) > 1,]
mydf$group <- "upregulated"
mydf$group[mydf$FC < 0] <- "downregulated"
mydf$othergroup <- "A"
mydf$othergroup[abs(mydf$FC) > 2] <- "B"

formula_res <- compareCluster(Entrez~group+othergroup, data=mydf, fun="enrichKEGG")

head(as.data.frame(formula_res))

## ----fig.height=7, fig.width=9-------------------------------------------
dotplot(ck)

## ----fig.height=6, fig.width=10------------------------------------------
dotplot(formula_res)
dotplot(formula_res, x=~group) + ggplot2::facet_grid(~othergroup)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
# [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"     "ENZYME"
# [8] "EVIDENCE"     "EVIDENCEALL"  "GENENAME"     "GO"           "GOALL"        "IPI"          "MAP"
#[15] "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"         "PROSITE"
#[22] "REFSEQ"       "SYMBOL"       "UCSCKG"       "UNIGENE"      "UNIPROT"
eg = bitr(glist, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(eg)
uniprot_ids <- bitr(glist, fromType="SYMBOL", toType=c("UNIPROT"), OrgDb="org.Hs.eg.db")
head(uniprot_ids)
refseq_ids <- bitr(glist, fromType="SYMBOL", toType=c("REFSEQ"), OrgDb="org.Hs.eg.db")
head(refseq_ids)
go_ids <- bitr(glist, fromType="SYMBOL", toType=c("UCSCKG"), OrgDb="org.Hs.eg.db")
head(go_ids)
go_ids <- bitr(glist, fromType="SYMBOL", toType=c("GOALL"), OrgDb="org.Hs.eg.db")
head(go_ids)
#eg2np <- bitr_kegg(glist, fromType='ncbi-geneid', toType='kegg', organism='hsa')
#bitr_kegg("Z5100", fromType="kegg", toType='ncbi-proteinid', organism='ece')
#bitr_kegg("Z5100", fromType="kegg", toType='uniprot', organism='ece')
library(DOSE)
na.omit(genelist)
gene <- names(glist)
gene.df <- bitr(glist, fromType = "SYMBOL",toType = c("ENTREZID", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
str(gene.df)
entrezgenes<-gene.df[,"ENTREZID"]
ggo <- groupGO(gene=entrezgenes, OrgDb=org.Hs.eg.db, ont="CC",
               level    = 3,readable = TRUE)
head(ggo)
kk <- enrichKEGG(gene = entrezgenes,organism='hsa',pvalueCutoff = 0.05)
head(kk)

gene.df <- bitr(AB.vs.wt1021.DE$GeneSymbol, fromType = "SYMBOL",toType = c("ENTREZID", "KEGG"),OrgDb=Org.Hs.egOMIM2EG@datacache)

S.me1021 <- enrichKEGG(gene = geneList,organism='sme',pvalueCutoff = 0.05)
head(S.me1021)

#############################################################################


AB.vs.wt1021.npid <- bitr_kegg(AB.vs.wt1021.DE$GeneSymbol, fromType='kegg', toType='ncbi-proteinid', organism='sme',drop=T)
head(AB.vs.wt1021.npid)
dim(AB.vs.wt1021.npid)

AB.vs.wt1021.geneid <- bitr_kegg(AB.vs.wt1021.DE$GeneSymbol, fromType='kegg', toType='ncbi-geneid', organism='sme')
dim(AB.vs.wt1021.geneid)

## ------------------------------------------------------------------------
ego <- enrichGO(gene=entrezgenes,
                universe=names(geneList),
                OrgDb= org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)

## ----eval=FALSE----------------------------------------------------------
## ego2 <- enrichGO(gene         = gene.df$ENSEMBL,
##                 OrgDb         = org.Hs.eg.db,
## 		keytype       = 'ENSEMBL',
##                 ont           = "CC",
##                 pAdjustMethod = "BH",
##                 pvalueCutoff  = 0.01,
##                 qvalueCutoff  = 0.05)

## ----eval=FALSE----------------------------------------------------------
## ego2 <- setReadable(ego2, OrgDb = org.Hs.eg.db)

## ----eval=FALSE----------------------------------------------------------
## ego3 <- gseGO(geneList     = geneList,
##               OrgDb        = org.Hs.eg.db,
##               ont          = "CC",
##               nPerm        = 1000,
##               minGSSize    = 100,
##               maxGSSize    = 500,
##               pvalueCutoff = 0.05,
##               verbose      = FALSE)



## ------------------------------------------------------------------------
barcodeplot(AB.vs.wt1021B.DE[,8], index = AB.vs.wt1021B.DE[,7],index2 = AB.vs.wt1021B.DE[,8], col.bars = "dodgerblue",alpha=.01,
            labels = "LogFoldChange",xlab="FoldChange")

barcodeplot(A.vs.wt1021.DE[,8], index = A.vs.wt1021.DE[,7],index2 = A.vs.wt1021.DE[,8], col.bars = "dodgerblue",alpha=.01,
            labels = "LogFoldChange",xlab="FoldChange")

barcodeplot(wt1021B.vs.wt1021.DE[,8], index = wt1021B.vs.wt1021.DE[,7],index2 = wt1021B.vs.wt1021.DE[,8], col.bars = "dodgerblue",alpha=.01,
            labels = "LogFoldChange",xlab="FoldChange")

barcodeplot(A.vs.AB.DE[,8], index = A.vs.AB.DE[,7],index2 = A.vs.AB.DE[,8], col.bars = "dodgerblue",alpha=.01,
            labels = "LogFoldChange",xlab="FoldChange")

barcodeplot(AB.vs.wt1021.DE[,8], index = AB.vs.wt1021.DE[,7],index2 = AB.vs.wt1021.DE[,8], col.bars = "dodgerblue",alpha=.01,
            labels = "LogFoldChange",xlab="FoldChange")

## ------------------------------------------------------------------------
kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)

## ----eval = FALSE--------------------------------------------------------
sme.genes <- enrichMKEGG(gene = geneList,
                         organism = 'sme')

## ----eval=FALSE----------------------------------------------------------
## mkk2 <- gseMKEGG(geneList = geneList,
##                  species = 'hsa')

## ----eval=FALSE----------------------------------------------------------
## david <- enrichDAVID(gene = gene,
##                      idType = "ENTREZ_GENE_ID",
##                      listType = "Gene",
##                      annotation = "KEGG_PATHWAY",
##                      david.user = "clusterProfiler@hku.hk")

## ------------------------------------------------------------------------
gmtfile <- system.file("extdata", "c5.cc.v5.0.entrez.gmt", package="clusterProfiler")
c5 <- read.gmt(gmtfile)

egmt <- enricher(gene, TERM2GENE=c5)
head(egmt)

egmt2 <- GSEA(geneList, TERM2GENE=c5, verbose=FALSE)
head(egmt2)

## ----fig.height=5, fig.width=9-------------------------------------------
barplot(mkk, drop=TRUE, showCategory=12)

## ----fig.height=5, fig.width=8-------------------------------------------
barplot(ego, showCategory=8)

## ------------------------------------------------------------------------
dotplot(mkk)

## ----
#fig.cap="enrichment map of enrichment result", fig.align="center", fig.height=16, fig.width=16, eval=FALSE----
# enrichMap(ego)

## ----fig.height=14, fig.width=14, eval=FALSE-----------------------------
## ## categorySize can be scaled by 'pvalue' or 'geneNum'
## cnetplot(ego, categorySize="pvalue", foldChange=geneList)

## ----fig.height=12, fig.width=8------------------------------------------
plotGOgraph(ego)

## ----fig.cap="plotting gsea result", fig.align="center", fig.height=6, fig.width=8----
gseaplot(kk2, geneSetID = "hsa04145")

## ----eval=FALSE----------------------------------------------------------
## browseKEGG(kk, 'hsa04110')

## ----eval=FALSE----------------------------------------------------------
## library("pathview")
## hsa04110 <- pathview(gene.data  = geneList,
##                      pathway.id = "hsa04110",
##                      species    = "hsa",
##                      limit      = list(gene=max(abs(geneList)), cpd=1))

## ------------------------------------------------------------------------
data(gcSample)
lapply(gcSample, head)

## ------------------------------------------------------------------------
ck <- compareCluster(geneCluster = gcSample, fun = "enrichKEGG")
head(as.data.frame(ck))

## ------------------------------------------------------------------------
mydf <- data.frame(Entrez=names(geneList), FC=geneList)
mydf <- mydf[abs(mydf$FC) > 1,]
mydf$group <- "upregulated"
mydf$group[mydf$FC < 0] <- "downregulated"
mydf$othergroup <- "A"
mydf$othergroup[abs(mydf$FC) > 2] <- "B"

formula_res <- compareCluster(Entrez~group+othergroup, data=mydf, fun="enrichKEGG")

head(as.data.frame(formula_res))

## ----fig.height=7, fig.width=9-------------------------------------------
dotplot(ck)

## ----fig.height=6, fig.width=10------------------------------------------
dotplot(formula_res)
dotplot(formula_res, x=~group) + ggplot2::facet_grid(~othergroup)

library(package = affyLib, character.only = TRUE)
## the distribution of the adjusted p-values
hist(geneList, 100)

## how many differentially expressed genes are:
sum(topDiffGenes(geneList))

## build the topGOdata class
GOdata <- new("topGOdata",ontology = "BP",
              allGenes = geneList,geneSel = topDiffGenes,
              annot = annFUN.db,affylib = affyLib)

## display the GOdata object
GOdata

 ##########################################################
## Examples on how to use the methods
##########################################################

 ## description of the experiment
description(GOdata)

## obtain the genes that will be used in the analysis
a <- genes(GOdata)
str(a)
numGenes(GOdata)

## obtain the score (p-value) of the genes
selGenes <- names(geneList)[sample(1:length(geneList), 10)]
gs <- geneScore(GOdata, whichGenes = selGenes)
print(gs)

## if we want an unnamed vector containing all the feasible genes
gs <- geneScore(GOdata, use.names = FALSE)
str(gs)

 ## the list of significant genes
 sg <- sigGenes(GOdata)
str(sg)
numSigGenes(GOdata)

                ## to update the gene list
.geneList <- geneScore(GOdata, use.names = TRUE)
GOdata ## more available genes
GOdata <- updateGenes(GOdata, .geneList, topDiffGenes)
GOdata ## the available genes are now the feasible genes

                ## the available GO terms (all the nodes in the graph)
go <- usedGO(GOdata)
length(go)

                ## to list the genes annotated to a set of specified GO terms
sel.terms <- sample(go, 10)
ann.genes <- genesInTerm(GOdata, sel.terms)
str(ann.genes)

                ## the score for these genes
 ann.score <- scoresInTerm(GOdata, sel.terms)
str(ann.score)

## to see the number of annotated genes
num.ann.genes <- countGenesInTerm(GOdata)
str(num.ann.genes)

## to summarise the statistics
termStat(GOdata, sel.terms)

