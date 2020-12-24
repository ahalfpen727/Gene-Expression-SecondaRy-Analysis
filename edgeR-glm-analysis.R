#1
genecounts <- read.table("18H18L.txt", row.names = 1)

colnames(genecounts) <- c(paste("18H",1:10,sep = "_"), paste("18L",1:10,sep="_"))

condition <-c(rep("H", 10), rep("L", 10))

dge <- DGEList(counts = genecounts, group = condition)

countsPerMillion <- cpm(dge)

countCheck <- countsPerMillion > 1

keep <- which(rowSums(countCheck)>=2)

dge <- dge[keep,]

dge <- calcNormFactors(dge)

design <- model.matrix(~ condition + 0, data = dge$samples)

colnames(design) = c("High", "Low")

disp <- estimateGLMCommonDisp(dge, design)

disp <- estimateGLMTrendedDisp(disp, design)

disp <- estimateGLMTagwiseDisp(disp, design)

fit <- glmFit(disp, design)

lrt <- glmLRT(fit)

res <- as.data.frame(lrt$table)

res <- cbind(res, FDR = p.adjust(res$PValue, method = "BH"))

fc1 <- res[res[,1] >= 1 | res[,1] <= -1,]

fc1p <- fc1[fc1$PValue <= 0.05,]

fc1q <- fc1[fc1$FDR <= 0.05,]


#2 
genecounts <- read.table("18H18L.txt", row.names = 1)

colnames(genecounts) <- c(paste("18H",1:10,sep = "_"),paste("18L",1:10,sep="_"))

condition <-c(rep("H", 10), rep("L", 10))

dge <- DGEList(counts = genecounts, group = condition)

countsPerMillion <- cpm(dge)

countCheck <- countsPerMillion > 0.5

keep <- which(rowSums(countCheck)>=2)

dge <- dge[keep,]

dge <- calcNormFactors(dge)

design <- model.matrix(~ condition + 0, data = dge$samples)

colnames(design) = c("High", "Low")

disp <- estimateGLMCommonDisp(dge, design)

disp <- estimateGLMTrendedDisp(disp, design)

disp <- estimateGLMTagwiseDisp(disp, design)

fit <- glmFit(disp, design)

highlow <- makeContrasts(High-Low,levels = design)

lrt <- glmLRT(fit, contrast =highlow)

res <- as.data.frame(lrt$table)

res <- cbind(res, FDR = p.adjust(res$PValue, method = "BH"))

fc1 <- res[res[,1] >= 1 | res[,1] <= -1,]

fc1p <- fc1[fc1$PValue <= 0.05,]

fc1q <- fc1[fc1$FDR <= 0.05,]

write.table(fc1p,"fc1p_18H18L_1.txt", quote = F, sep = "\t", row.names = T, col.names = T)



#3
genecounts <- read.table("18H18L.txt", row.names = 1)

colnames(genecounts) <- c(paste("18H",1:10,sep = "_"),paste("18L",1:10,sep="_"))

condition <-c(rep("H", 10), rep("L", 10))

dge <- DGEList(counts = genecounts, group = condition)

countsPerMillion <- cpm(dge)

countCheck <- countsPerMillion > 0.5

keep <- which(rowSums(countCheck)>=2)

dge <- dge[keep,]

dge <- calcNormFactors(dge)

design <- model.matrix(~ condition + 0, data = dge$samples)

colnames(design) = c("High", "Low")

disp <- estimateGLMCommonDisp(dge, design)

disp <- estimateGLMTrendedDisp(disp, design)

disp <- estimateGLMTagwiseDisp(disp, design)

fit <- glmFit(disp, design)

lrt <- glmLRT(fit)

res <- as.data.frame(lrt$table)

res <- cbind(res, FDR = p.adjust(res$PValue, method = "BH"))

fc1 <- res[res[,1] >= 1 | res[,1] <= -1,]

fc1p <- fc1[fc1$PValue <= 0.05,]

fc1q <- fc1[fc1$FDR <= 0.05,]

#2-1
library(gplots)

data1 <- read.table("fc1p_18H18L_1.txt", header = T)

data2 <- read.table("fc1p_18H30H_1.txt", header = T)

data3 <- read.table("fc1p_18L30L_1.txt", header = T)

data4 <- read.table("fc1p_30H30L_1.txt", header = T)

test1 <- data1[,1]

test2 <- data2[,1]

test3 <- data3[,1]

test4 <- data4[,1]

venn(list(A = test1, B = test2, C = test3, D = test4))

A group = 105, B group = 1778, C group = 1665, D group = 47
