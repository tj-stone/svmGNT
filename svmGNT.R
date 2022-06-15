#GNT SVM 3 Group Classifier build process:
library(minfi)
library(ConsensusClusterPlus)
library(Rtsne)
library(ggplot2)
library(caret)
library(doParallel)
library(MLmetrics)

#Read in raw training (450k) and test (EPIC) data and combine into one array
baseDir <- getwd()
targets <- read.metharray.sheet(baseDir, pattern = "SVMGNT_cohort.csv")
RGSetTraining <- read.metharray.exp(targets = targets[1:42,], force = TRUE)
RGSetTest <- read.metharray.exp(targets = targets[43:76,], force = TRUE)
RGSet <- combineArrays(RGSetTraining, RGSetTest)

#After reading in data to RGSet, check proportion of failed probes, identify samples with >10% failure
RGSet$Sample_Plate[1:42] <- "Training"
RGSet$Sample_Plate[43:76] <- "Test"
detP <- detectionP(RGSet)
failed <- detP > 0.01
failedProbes <- colMeans(failed)
sum(failedProbes > 0.1)
badArrays <- paste(names(failedProbes[failedProbes > 0.1]),
                  failedProbes[failedProbes > 0.1], 
                  sep="-")
subLabels <- names(failedProbes[failedProbes > 0.1])
RGSet$Sample_Plate[colnames(RGSet) %in% subLabels] <- "Suboptimal"

#Normalise dataset with SWAN normalisation.
MSet <- preprocessFunnorm(RGSet)

#Get annotation & remove probes on X and Y chromosomes
A <- getAnnotation(MSet)
remove <- rownames(A)[A$chr == "chrX" | A$chr == "chrY"]
MSet2 <- MSet[!rownames(MSet) %in% remove,]

#Filter SNP probes, mismatches, & MAF 0.05 from Pidsley et al.
filter_probes <- read.table("850Kremove_probes.txt")
MSet3 <- MSet2[!rownames(MSet2) %in% filter_probes$V1]

#Get beta values and split cohort into training and test subcohorts
betaAll <- na.omit(getBeta(MSet3))
betaTrain <- betaAll[,1:42]
betaTest <- betaAll[,43:76]

#Calculate median absolute deviation and select 10,000 most variable probes in training cohort
madsTraining <- apply(betaTrain, 1, mad)
top10kTraining <- betaTrain[rev(order(madsTraining))[1:10000],]
colnames(top10kTraining) <- MSet3$Sample_Name[1:42]

#Extract corresponding 10k training probes from all (training + test) cohort data
clusterProbes <- rownames(top10kTraining)[1:10000]
top10kAll <- betaAll[rownames(betaAll) %in% clusterProbes,]

#Consensus clustering to derive putative class membership for test cases alongside training cohort
ConsensusAll <- ConsensusClusterPlus(top10kAll, 
                                     maxK = 6, 
                                     reps = 10000, 
                                     clusterAlg = "hc", 
                                     title = "WardD2 Combined Consensus K1-6", 
                                     innerLinkage = "ward.D2", 
                                     finalLinkage = "ward.D2")

#t-SNE to visualise consensus clustering classification of training and test cohorts
set.seed(14)
tsne <- Rtsne(t(top10kAll),
              check_duplicates = FALSE,
              pca = FALSE,
              perplexity = 20,
              theta = 0.5,
              dims = 2)

#Plot t-SNE result, coloured by consensus classification w/ highlights for bad arrays
tiff("GNT_all_clustering_TSNE.tiff", units = "in", width = 5.5, height = 4, res = 600)
qplot(tsne$Y[,1], 
      tsne$Y[,2],
      color = as.factor(ConsensusAll[[3]]$consensusClass)) +
      labs(x = "TSNE-1",
      y = "TSNE-2") + 
      scale_colour_manual(values = c("indianred",
                                     "springgreen3", 
                                     "steelblue2")) +
      theme_classic() +
      theme(legend.position = "none") +
      geom_point(size = 3,
                 shape = factor(MSet3$Sample_Plate,
                                levels = c("Training", "Test", "Suboptimal"))) +
      theme(text=element_text(size = 14))
dev.off()

#Append training probes with target class of training cohort
top10kTraining <- rbind(top10kTraining, ConsensusAll[[3]]$consensusClass[1:42])
row.names(top10kTraining)[10001] <- "class"
top10kTraining <- as.data.frame(t(top10kTraining))
top10kTraining$class[top10kTraining$class == 1] <- "Group1"
top10kTraining$class[top10kTraining$class == 2] <- "Control"
top10kTraining$class[top10kTraining$class == 3] <- "Group2"

#Check for predictors with zero or near zero variance
nzv <- nearZeroVar(top10kTraining[,1:10000], saveMetrics = TRUE)
summary(c(as.factor(nzv$nzv), as.factor(nzv$nzv)))

#Check for abundance of highly correlated predictors
trainingCor <-  cor(top10kTraining[,1:10000])
summary(trainingCor[upper.tri(trainingCor)])
trainingHighCor <- findCorrelation(trainingCor, cutoff = .95)
length(trainingHighCor)

#Set training parameters
#Extract estimated sigma using kernlab sigest (mirrors caret sigma estimation while allowing custom Cost scaling)
models <- getModelInfo("svmRadial",
                       regex=FALSE)[[1]]
preProcValues <- preProcess(top10kTraining,
                            method = c("center", "scale"))
processData <- predict(preProcValues, top10kTraining)
x <- model.matrix(class~.,
                  data = processData)[,-1]
y <- processData$class
set.seed(1)
sigmas <- kernlab::sigest(as.matrix(x),
                          na.action = na.omit,
                          scaled = TRUE)

#Expand parameter search grid using estimated mean sigma alongside Cost on log10 scale
svmGrid <- expand.grid(sigma = c(mean(sigmas)[-2]), C = 10^(-3:1))

#Specify training control parameters with leave-one-out crossval, prediction probabilities, and mult-class performance metrics.
#Also upsample all training groups to accomodate for differences in original group sizes.
trctrl <- trainControl(method = "LOOCV",
                       classProbs = T,
                       summaryFunction = multiClassSummary,
                       sampling = "up")

#Enable parallel processing to reduce training time
cl <- makePSOCKcluster(4)
registerDoParallel(cl)

#Train support vector machine model using radial bias function kernel, data normalised/scaled.
svmGNT <- train(class~.,
                data = top10kTraining, 
                method = "svmRadial",
                trControl = trctrl,
                tuneGrid = svmGrid,
                preProcess = c("center", "scale"))
stopCluster(cl)

#Plot accuracy of the trained model
plot(svmGNT)

#Assess all training + test cohort data using the prediction model, gather classification and probability scores
preds <- predict.train(svmGNT,
                       t(top10kAll),
                       type = "prob")

predictions <- colnames(preds)[apply(preds, 1, which.max)]
SVMScores <- apply(preds, 1, max)
classifications <- data.frame(predictions, SVMScores, row.names = colnames(top10kAll))
colnames(classifications) <- c("Predicted Class", "Score")
write.csv(classifications, file = 'svmGNT_classifications.csv')

#Visualise classification against previously generated t-SNE
tiff("GNT_svm_predictions_TSNE.tiff", units = "in", width = 5.5, height = 4, res = 600)
qplot(tsne$Y[,1], 
      tsne$Y[,2],
      color=as.factor(classifications$`Predicted Class`)) +
      labs(x = "TSNE-1",
           y = "TSNE-2") + 
      scale_colour_manual(values = c("springgreen3",
                                     "indianred", 
                                     "steelblue2")) +
      theme_classic() +
      theme(legend.position = "none") +
      geom_point(size = 3,
                 shape = factor(MSet3$Sample_Plate,
                 levels = c("Training", "Test", "Suboptimal"))) +
      theme(text=element_text(size = 14))
dev.off()
