#building the classification tree
#install if necessary
install.packages('rpart');install.packages('rpart.plot');install.packages('multtest')
install.packages('rattle');install.packages(golub);install.packages("tree")
biocM(c("car", "VennDiagram","nortest"))
# load libraries
library(rpart);library(rpart.plot);library(multtest);library(tree)
library(cummeRbund);library(outliers);library(nortest);library(stats)
library(knitr);library(car); library(VennDiagram);library(stats4)
library(genefilter);library(qvalue);library(RColorBrewer)
library(ggplot2);library(gplots);library(gridExtra)
library(limma);library(multtest);library(edgeR)
library(genefilter);library(qvalue);library(dplyr)
library(caret);library(rsample);data(golub)

# plotting regression trees
cuff<- readCufflinks(dir='/media/drew/easystore/umb_triley/urine1/cuffdiff_results_hg38_default/LUTS-over-CTRL/',
                     genome="hg38",gtfFile='/media/drew/easystore/ReferenceGenomes/UCSC_hg38/genes.gtf',
                     rebuild=F)
cuff
replicates.info<-cummeRbund::replicates(cuff)
groups<-replicates.info$sample_name
samples<-replicates.info$rep_name
conditions<-as.factor(groups)
replicates.info
conditions
samples
groups

over=groups[((length(groups)/2)+1)]
under=groups[1]
over;under

# Simple design matrix and contrast matrix
design <- model.matrix(~0 + groups, data=replicates.info)
colnames(design) <- levels(conditions)
row.names(design) <- samples
design

contr.matrix <- makeContrasts(lutsVSctrl = CTRL - LUTS,
                              levels = colnames(design))
contr.matrix
# gene expr data
genes_exp.df<-diffData(cummeRbund::genes(cuff))
sig.genes_exp.df<-subset(genes_exp.df, genes_exp.df$significant=="yes")

#all_exp.df<-diffTable(cummeRbund::genes(cuff))
g.cnt.matrix<-repCountMatrix(cummeRbund::genes(cuff))
# factors for conditions
under.group<-grep(pattern=under, colnames(g.cnt.matrix))
over.group<-grep(pattern=over, colnames(g.cnt.matrix))

# Lane bias correction design matrix and contrast matrix
Lane.factor <- as.factor(lanes$Lane)
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
groups<-replicates.info$sample_name
samples<-replicates.info$rep_name
under=groups[1]
over=groups[((length(groups)/2)+1)]
conditions <- as.factor(groups)
sampleCondition<-c(rep(under,length(groups)/2),rep(over,length(groups)/2))
# design matrix \
#design.sample <- model.matrix(~0 + rep_name, data=replicates.info)
design <- model.matrix(~0 + batcheffects, data=replicates.info)
row.names(design) <- samples
# contrast matrix
contr.matrix <- makeContrasts(CTRL_L1T4 - CTRL_L5T8, CTRL_L5T8 - LUTS_L1T4, LUTS_L1T4 - LUTS_L5T8,
                              levels = c("CTRL_L1T4", "CTRL_L5T8", "LUTS_L1T4", "LUTS_L5T8"))
# gene expr data
cuff_annotation_data<-featureNames(cuff@genes)
gene_exp.diff<-diffData(cummeRbund::genes(cuff))
# reformat all_exp
sig_genes.df<-subset(gene_exp.diff, gene_exp.diff$significant=="yes")
g.cnt.df<-repCountMatrix(cummeRbund::genes(cuff))
# gene expr data
g.cnt.ma<-repCountMatrix(cummeRbund::genes(cuff))
# factors for conditions
under.group<-grep(pattern=under, colnames(g.cnt.ma))
over.group<-grep(pattern=over, colnames(g.cnt.ma))

g.count.df=as.data.frame(g.cnt.df)
g.count.ma=as.matrix(g.cnt.ma)
inds <- which(!is.na(g.count.df$GOid) & !is.na(g.count.df$EntrezID))
# factors for conditions
g.count.df<-g.count.df[inds,]
g.cnt.df<-g.cnt.df[inds,]
g.cnt.ma<-g.cnt.ma[inds,]

under.group<-grep(pattern=under, colnames(g.cnt.df))
over.group<-grep(pattern=over, colnames(g.cnt.df))

g.o.cnt.df<-g.cnt.df[,over.group]
g.u.cnt.df<-g.cnt.df[,under.group]

tree1 <- tree(Species ~ Sepal.Width + Petal.Width, data = iris)
summary(tree1)

plot(tree1)
text(tree1)

plot(iris$Petal.Width,iris$Sepal.Width,pch=19,col=as.numeric(iris$Species))
partition.tree(tree1,label="Species",add=TRUE)
legend(1.75,4.5,legend=unique(iris$Species),col=unique(as.numeric(iris$Species)),pch=19)
graph <- qplot(Petal.Width, Sepal.Width, data=iris, colour=Species, size=I(4))
graph + geom_hline(aes(yintercept=2.65)) + geom_vline(aes(xintercept=0.8)) + geom_vline(aes(xintercept=1.75)) + geom_vline(aes(xintercept=1.35))

tree1 <- tree(Species ~ Sepal.Width + Sepal.Length + Petal.Length + Petal.Width, data = iris)
summary(tree1)


plot(tree1)
text(tree1)

rpart <- rpart(Species ~ ., data=iris, method="class",)
rpart
# plot decision tree
fancyRpartPlot(rpart, main="Iris")

set.seed(123)
n<-10 ; sigma <- 0.5
fac <- factor(c(rep(1,n),rep(2,n),rep(3,n)))
levels(fac) <- c("ALL1","ALL2","AML")
geneA <- c(rnorm(10,0,sigma),rnorm(10,2,sigma),rnorm(10,4,sigma))
dat <- data.frame(fac,geneA)library(rpart)
rp  <-  rpart(fac ~ geneA, method="class",data=dat)plot(rp, branch=0,margin=0.1)
text(rp, digits=3, use.n=TRUE)

geneA <- c(rnorm(20,0,sigma),rnorm(10,2,sigma))
geneB <- c(rnorm(10,0,sigma),rnorm(20,2,sigma))
geneC <- c(rnorm(30,1,sigma))
dat <- data.frame(fac,geneA,geneB,geneC)
rp  <-  rpart(fac ~ geneA  + geneB + geneC, method="class",data=dat)

gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
gol.rp  <-  rpart(gol.fac ~ golub[1042,] , method="class")
predictedclass <- predict(gol.rp, type="class")
table(predictedclass, gol.fac)

summary(gol.rp)summary(gol.rp)
rpart(formula = gol.fac ~ golub[1042, ], method = "class")
predict(gol.rp,type="class")
predict(gol.rp, type="prob")

row.names(golub)<- paste("gene", 1:3051, sep = "")
goldata <- data.frame(t(golub[1:3051,]))
gol.fac <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
gol.rp  <-  rpart(gol.fac~., data=goldata, method="class", cp=0.001)
plot(gol.rp, branch=0,margin=0.1);
text(gol.rp, digits=3, use.n=TRUE)
golub.gnames[896,]

library("hgu95av2.db");library(ALL);data(ALL)
ALLB123 <- ALL[,ALL$BT %in% c("B1","B2","B3")]
pano <- apply(exprs(ALLB123), 1, function(x) anova(lm(x ~ ALLB123$BT))$Pr[1])
names <- featureNames(ALL)[pano<0.000001]
symb <- mget(names, env = hgu95av2SYMBOL)
ALLBTnames <- ALLB123[names, ]
probedat <- as.matrix(exprs(ALLBTnames))
row.names(probedat)<-unlist(symb)
diagnosed <- factor(ALLBTnames$BT)
tr <- rpart(factor(ALLBTnames$BT) ~ ., data = data.frame(t(probedat)))
plot(tr, branch=0,margin=0.1);  text(tr, digits=3, use.n=TRUE)
rpartpred <- predict(tr, type="class")
table(rpartpred,diagnosed)
predicted.class <- predict(tr, type="class")
predicted.probabilities <- predict(tr, type="prob")
out <- data.frame(predicted.probabilities,predicted.class,diagnosis=factor(ALLBTnames$BT))
print(out,digits=2)

# Training and validation
i <- sample(1:78, 39, replace = FALSE)
noti <- setdiff(1:78,i)
df <- data.frame(Y = factor(ALLBTnames$BT), X =t(probedat))
rpart.est <- rpart(Y ~ ., data = df, subset=i)
rpart.pred.t <- predict(rpart.est, df[i,], type="class")
table(rpart.pred.t,factor(ALLBTnames$BT[i]))

rpart.pred.v <- predict(rpart.est,df[noti,], type="class")
table(rpart.pred.v,factor(ALLBTnames$BT[noti]))

#SVM
library(e1071)
df <- data.frame(Y = factor(ALLBTnames$BT), X =t(probedat))
Y <- factor(ALLBTnames$BT);X <- t(probedat)
svmest <- svm(X, Y, data=df, type = "C-classification", kernel = "linear")
svmpred <- predict(svmest, X, probability=TRUE)
table(svmpred, factor(ALLBTnames$BT))

Yt <- factor(ALLBTnames$BT)[i]
Yv <- factor(ALLBTnames$BT)[noti]
X <- t(probedat); Xt <- X[i,]; Xv <- X[noti,]
svmest <- svm(Xt, Yt, type = "C-classification", kernel = "linear")
svmpredt <- predict(svmest, Xt, probability=TRUE)
table(svmpredt, Yt)
svmpredv <- predict(svmest, Xv, probability=TRUE)
table(svmpredv, Yv)


#NeuralNets
Y <- factor(ALLBTnames$BT)
X <- t(probedat)
> library(nnet)> df <- data.frame(Y = Y, X = X[, sample(ncol(X), 20)])> nnest <- nnet(Y ~ .,data = df, size = 5, maxit = 500, decay = 0.01,+ MaxNWts = 5000)> pred <- predict(nnest, type = "class")> table(pred, \
                                                                                                                                                                                                            Y) # prints confusion ma


> nnest.t <- nnet(Y ~ ., data = df,subset=i, size = 5,decay = 0.01,+ maxit=500)> prednnt <- predict(nnest.t, df[i,], type = "class")> table(prednnt,Ytrain=Y[i])

prednnv <- predict(nnest.t, df[noti,], type = "class")
table(prednnv, Yval= Y[noti])

