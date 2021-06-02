#sourcePath = "C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Butterfly/"

#source(paste0(sourcePath,"ButterflyHelpers.R"))

#localDataPath = "C:/Work/R/ButterflyQuant/"
#cachedDataPath = paste0(localDataPath,"savedData/")
dataPath = "E:/Butterfly/"
#figure_data_path = "C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Butterfly/FigureData/"
#figure_path = "Z:/projects/Butterfly/figures/"

library(ggplot2)
library("ggpubr")
library(BUSpaRse)
library(textTinyR)


#read data prediction:
########################################################################

fullMatrix = read_count_output(paste0(dataPath,"pbmcv2/bus_output/count/"), "output", tcc=F)
dsMatrix = read_count_output(paste0(dataPath,"pbmcv2/bus_output/count_ds/"), "output", tcc=F)
predMatrix = read_count_output(paste0(dataPath,"pbmcv2/bus_output/count_ds_pred/"), "output", tcc=F)



#sum up the counts of all cells
full = sparse_Sums(fullMatrix, T)
ds = sparse_Sums(dsMatrix, T)
pred = sparse_Sums(predMatrix, T)

#check that the genes are on the same rows (they should be)
fullNames = row.names(fullMatrix)
dsNames = row.names(dsMatrix)
predNames = row.names(predMatrix)
sum(fullNames != dsNames) #should be 0
length(fullNames) == length(dsNames)#should be TRUE
sum(fullNames != predNames)#should be 0
length(fullNames) == length(predNames)#should be TRUE
#OK!

#now TPM + log transform
tpmFull = full*10^6/sum(full)
sum(tpmFull) #should be 10^6
fullL = log2(tpmFull + 1)

tpmDs = ds*10^6/sum(ds)
sum(tpmDs) #should be 10^6
dsL = log2(tpmDs + 1)

tpmPred = pred*10^6/sum(pred)
sum(tpmPred) #should be 10^6
predL = log2(tpmPred + 1)

#now compare prediction and ds vs the full dataset

x = dsL

yDS = dsL - fullL
plot(x,yDS, main="Downsampled", xlab="log2(CPM + 1)", ylab="Log2 FC vs full sample", ylim = c(-2.5,1.5))
lines(c(0,14), c(0,0), col="red", lw=2)

yPred = predL - fullL
plot(x,yPred, main="Downsampled and Predicted", xlab="log2(CPM + 1)", ylab="Log2 FC vs full sample", ylim = c(-2.5,1.5))
lines(c(0,14), c(0,0), col="red", lw=2)

#The predicted is expected to have less spread

