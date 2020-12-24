source("http://bioconductor.org/biocLite.R")

options(error=traceback) # causes a traceback to appear if there is an error, and the traceback has the line number, prefixed by #
## options(echo=TRUE)

sessionInfo() # Keep track of versions used, especially packages


####################################################################################
# Evaluate command line arguments ##################################################
####################################################################################
# dir = "~/Documents/UMB/Riley_Lab/HiSeq_data/jose/"
## args = commandArgs(trailingOnly = T)
args = commandArgs(TRUE)

for (i in 1:length(args)) {
    eval(parse(text=args[[i]]))
}

####################################################################################
# Functions ########################################################################
####################################################################################

head2 = function(object)
{
    if (is.null(object)) {
        return("NULL")
    }
    return(head(object))
}

####################################################################################
####################################################################################
####################################################################################


## cuffdiff.table = read.xlsx2(paste(inDir, inFile, sep = ""), startRow = 2)
tablePathName = paste( "gene_exp.diff") #, sep = "gene_exp.diff")
writeLines(c("\ntablePathName:", tablePathName))

cuffdiff.table = read.table(tablePathName, header=TRUE)

## writeLines("\ncuffdiff.table:", head2(cuffdiff.table))
writeLines("\ncuffdiff.table:")
cuffdiff.table[0:10, ]

## biocLite("pathview")
## biocLite("gage")

require(gage)
require(pathview)

 cuffdiff.table = read.delim(file = "gene_exp.diff", sep = "\t")

## Cufflinks is a very popular RNA-Seq analysis tool.  It is developed by the same group as TopHat.  It is
## implemented independent of Bioconductor.  However, we can read its differential expression analysis
## results into R easily.  The result file is named *_exp_diff
## Notice that the gene symbols need to be converted to Entrez Gene IDs, which are used in
## KEGG pathways (for many research species) and GO gene sets.

# notice the column name special character changes. The column used to be
# cuffdiff.table$log2.fold_change. for older versions of Cufflinks.
# cuff.fc = cuffdiff.table$log2.FPKMy.FPKMx
cuff.fc = cuffdiff.table$log2.fold_change.

writeLines(c("\ncuff.fc: ", head2(cuff.fc)))

gnames = cuffdiff.table$gene
sel = gnames != "-"
## writeLines("sel:")
## sel[0:10]

gnames = as.character(gnames[sel])
## writeLines("gnames:")
## gnames[0:10]

cuff.fc = cuff.fc[sel]
## cuff.fc[0:10]

names(cuff.fc) = gnames
## gnames[0:10]

# Notice that the gene sybmols need to be converted to Entrez Gene IDs, which are used in
# KEGG pathways (many research species) and GO gene sets.
gnames.eg = pathview::id2eg(gnames, category = "symbol")
sel2 = gnames.eg[,2]>""

cuff.fc = cuff.fc[sel2]
names(cuff.fc) = gnames.eg[sel2,2]
writeLines(c("\nrange(cuff.fc): ", range(cuff.fc)))
head(cuff.fc)
# remove the -Inf and Inf values, which block the downstream analysis
cuff.fc[cuff.fc > 10] = 10
cuff.fc[cuff.fc < -10] = -10

exp.fc = cuff.fc

head(exp.fc)

## Next, we use GAGE for pathway analysis, and Pathview for visualization.  Notice that
## this step (the same code) is identical for DESeq, edgeR, Limma and Cufflinks workflows.

# In pathview function, please be conscious of the argument gene.idtype. Default is gene.idtype="entrez",i.e. Entrez
# Gene, which are the primary KEGG gene ID for many common model organisms. For other species, gene.idtype
# should be set to "KEGG" since KEGG uses other types of gene IDs.

data(kegg.gs)
data("paths.hsa")
writeLines(c("\nlapply(kegg.gs[1:3], head2, 3):", lapply(kegg.gs[1:3], head2,3)))
writeLines(c("\nlength(kegg.gs): ", length(kegg.gs)))
writeLines(c("\nexp.fc:", head2(exp.fc)))
writeLines(c("\nstr(exp.fc):", str(exp.fc)))

fc.kegg.p <- gage(exp.fc, gsets = kegg.gs, ref = NULL, samp = NULL)

writeLines(c("\nfc.kegg.p$greater:", head2(fc.kegg.p$greater)))

sel <- fc.kegg.p$greater[, "p.val"] < 0.25 & !is.na(fc.kegg.p$greater[, "p.val"])
sel

path.ids <- rownames(fc.kegg.p$greater)[sel]
writeLines(c("\npath.ids:", head2(path.ids)))
writeLines(c("\nlength(path.ids):", length(path.ids)))

sel.l <- fc.kegg.p$less[, "q.val"] < 0.25 & !is.na(fc.kegg.p$less[,"q.val"])

path.ids.l <- rownames(fc.kegg.p$less)[sel.l]
writeLines(c("\npath.ids.l:", head2(path.ids.l)))
writeLines(c("\nlength(path.ids.l):", length(path.ids.l)))

path.ids2 <- substr(c(path.ids, path.ids.l), 1, 8)
writeLines(c("\npath.ids2:", head2(path.ids2)))
writeLines(c("\nlength(path.ids2):", length(path.ids2)))


# view first 3 pathways as demo
# species = "hsa" is for human
top10 = max(20, length(path.ids2))
pv.out.list <- sapply(path.ids2[1:top10], function(pid) pathview(gene.data =  exp.fc, pathway.id = pid,
                                                             species = "hsa", out.suffix = "cuff"))

pv.out.list


## Here we used GAGE to infer the significant pathways.  But we are not limited to these pathways.  We can
## use Pathview to visualize RNA-Seq data (exp.fc here) on all interesting pathways directly.  We may also do GO
## and other types of gene set analysis as described in native workflow above.

