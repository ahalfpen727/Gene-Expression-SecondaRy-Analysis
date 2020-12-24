#!bin/bash
 
if [ ! -d ${OUTPUT} ]; then
    mkdir ${OUTPUT}
fi

R CMD BATCH --no-restore --no-save '--args diffDir=\"${diffDir}\" inDir=\"${INPUT}/${diffDir}\" outDir=\"${OUTPUT}/${diffDir}\"  under=\"${under}\" over=\"${over}\" Rplots=\"${OUTPUT}/${diffDir}/${R_PLOTS}\" refgtf=\"${CUMMERBUND_INPUT_GTF}\" genomeR=\"${R_GENOME}\" FPKMmatrix=\"${OUTPUT}/${diffDir}/${FPKM_MATRIX}\" geneList=\"${GENES_OF_INTEREST}\"' ${1}.R ${OUTPUT}/${diffDir}/${1}.Rout"
