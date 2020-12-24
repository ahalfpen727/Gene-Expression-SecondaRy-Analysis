#!/bin/sh
# sequence data
export lanes="/media/drew/easystore/umb_triley/urine1/Sample-Library-Preparation/lane-and-sample-numbers.csv"
export genomeR="/media/drew/easystore/umb_triley/ReferenceGenomes/UCSC_hg38/genome.fa"
export refgtf="/media/drew/easystore/umb_triley/ReferenceGenomes/UCSC_hg38/genes.gtf"
export inDir="/media/drew/easystore/umb_triley/urine1/cuffdiff_results_hg38_default/LUTS-over-CTRL"
export gtfDir="/media/drew/easystore/umb_triley/urine1/cuffcompare_results_hg38_gtf_guided"
export cuffcmp="/media/drew/easystore/umb_triley/urine1/cuffcompare_results_hg38_gtf_guided/cuffcmp.combined.gtf"
export outDir="/media/drew/easystore/umb_triley/urine1/gene_expRession_analysis"
export countma="/media/drew/easystore/umb_triley/urine1/gene_expRession_analysis/count_matrix.csv"
export Difftable="/media/drew/easystore/umb_triley/urine1/gene_expRession_analysis/DiffTable.csv"
export Rplots="/media/drew/easystore/umb_triley/urine1/gene_expRession_analysis/Rplots.pdf"
export over="LUTS"
export under="CTRL"
#----CUMMERBUND PARAMETERS--------------------------------------------------------------------------------------------------------------------------
CUMMERBUND_OUTPUT=${CUMMERBUND}${OUTPUT_LABEL}
CUMMERBUND_GENE_LIST='CXCL12,CXCR4,TGFB,MAPK,IL8,CCL5'
CUMMERBUND_INPUT_GTF="${CUFFCOMPARE}${OUTPUT_LABEL}/cuffcmp.combined.gtf"
#CUMMERBUND_INPUT_GTF="${CUFFMERGE}${OUTPUT_LABEL1}/merged.gtf"
#CUMMERBUND_INPUT_GTF="${REFERENCE_GTF}"
R_GENOME="hg38"
FPKM_MATRIX="fpkmMatrix.csv"

export FPKM_MATRIX
export GENOME_PATH
export R_GENOME
export CUMMERBUND_INPUT_GTF
export GENE_LIST
#---------------------------------------------------------------------------------------------------------------------------------------------------


#----BIONET PARAMETERS------------------------------------------------------------------------------------------------------------------------------
 BIONET_OUTPUT=${BIONET}${OUTPUT_LABEL}
 BIONET_GENE_ALPHA_VALUE=0.5
 BIONET_SUBNET_ALPHA_VALUE=0.1
 BIONET_NUMBER_NETWORKS=20

 export BIONET_GENE_ALPHA_VALUE
 export BIONET_SUBNET_ALPHA_VALUE
 export BIONET_NUMBER_NETWORKS

#----GO PARAMETERS----------------------------------------------------------------------------------------------------------------------------------
GO_OUTPUT=${GO}${OUTPUT_LABEL}
GO_DIFF_EXPRESSED_ALPHA_VALUE=0.1   # Alpha value for classifying genes as differentially expressed
GO_HYPER_GEO_ALPHA_VALUE=0.01   #Alpha value for enrichment of differentially expressed genes is a category

export GO_DIFF_EXPRESSED_ALPHA_VALUE
export GO_HYPER_GEO_ALPHA_VALUE
#----PATHVIEW PARAMETERS--------------------------------------------------------------------------\
--------------------------------------------------
PATHVIEW_OUTPUT=${PATHVIEW}${OUTPUT_LABEL}
PATHVIEW_SPECIES="human"   # pathview common.name ("human" or "mouse", etc.)
PATHVIEW_ALPHA_DOWN=0.2  # Alpha value for genes that are under expressed in the case
PATHVIEW_ALPHA_UP=0.2   # Alpha value for genes that are over expressed in the case
PATHVIEW_NUMBER_PATHWAYS_DOWN=30
PATHVIEW_NUMBER_PATHWAYS_UP=30
PATHVIEW_PATHWAYS_LIST=""

export PATHVIEW_SPECIES
export PATHVIEW_ALPHA_DOWN
export PATHVIEW_ALPHA_UP
export PATHVIEW_NUMBER_PATHWAYS_DOWN
export PATHVIEW_NUMBER_PATHWAYS_UP
export PATHVIEW_PATHWAYS_LIST
# OUTPUT
if [ ! -d ${outDir} ]; then
    mkdir ${outDir}
fi

# Plots
if [ ! -f ${Rplots} ]; then
    touch ${Rplots}
fi
