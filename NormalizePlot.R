#These functions can be implemented to normalize NanoString binding data (DxData), assess sample species based on probe binding ratios (ProbeRatio), and visualise correlation of binding profiles between two sample sets (e.g. test and reference) (DxCorrelation).

# Required R packages: compositions, reshape2, ggplot2. R 3.6 or higher.
library(compositions)
library(reshape2)
library(ggplot2)

# DxData: Normalize NanoString output with positive and negative control probe values. Data can also be normalized via subtraction of the average blank sample binding values. Input data should be provided as a matrix, with columns corresponding to samples and rows corresponding to probes. Vectors containing the corresponding column sample names for test samples, reference samples, and blank samples should also be provided. Please see example data for further details. If you are missing either test or reference samples, please provide the same vector for both test and reference arguments. Please provide unique sample names. This function will return a csv file of normalized binding values per probe for each sample.
#---------------------------------------------------------------

DxData <- function(InputData, ReferenceID, TestID, BlankID){
  #InputData should be a matrix of probe binding values, all other arguments should be sample name vectors. All sample names provided should be unique.
  DataVec <- unique(c(ReferenceID, TestID))
  if(length(BlankID > 1)){
    InputData$DxBlank <- meanRow(InputData[,colnames(InputData) %in% BlankID])
    BlankID <- c("DxBlank")
  }
  Norm <- InputData[14:34,] - InputData[14:34, colnames(InputData) %in% BlankID]
  Norm <- rbind(InputData[1:13,], Norm)
  Norm[Norm < 0] <- 0
  Norm <- as.matrix(Norm[,colnames(Norm) %in% DataVec])
  Geo <- as.vector(geometricmeanCol(Norm[2:7,]))
  MeanPos <- mean(Norm[2:7,])
  NormFactor <- as.vector(MeanPos/Geo)
  Norm <- rbind(Norm, NormFactor)
  PosData <- matrix(nrow=nrow(Norm[14:34, ]),ncol=0)
  temp <- matrix(nrow=nrow(Norm[14:34, ]),ncol=0)
  for (i in 1:ncol(Norm)){
    temp <- Norm[14:34, i]*Norm[35,i]
    PosData <- cbind(PosData, temp)
  }
  colnames(PosData) <- colnames(Norm)
  GeoMean <- geometricmean(Norm[8:13, ])
  NormData <- pmax(PosData[,]-GeoMean, 0)
  write.csv(NormData, "./NormalizedData.csv")
}

# ProbeRatio: Calculate ratio between strongest binding Candida species-specific probe and strongest binding pan-Candida probe per sample. Input should be normalized binding values for all candida probes, without positive/negative control probe binding values. This outputs a true/false statement per sample, based on whether that sample meets the probe binding ratio criteria indicating that is likely a Candida species (true) or is unlikely to be a Candida species (false). Ratio's of normalized binding < 10 and > 0.01 are used to set this criteria.
#--------------------------------------------------------------

ProbeRatio <- function(NormalizedData){
  #NormalizedData should be a matrix of normalized Candida probe binding values.
  Norm <- NormalizedData
  CandidaCall <- c()
  CandidaRatio <- c()
  for (i in 1:ncol(Norm)){
    species <- max(Norm[1:19,i])
    pan <- max(Norm[20:21,i])
    ratio <- species/pan
    Call <- ratio>0.01 & ratio<10
    CandidaCall <- c(CandidaCall, Call)
    CandidaRatio <- c(CandidaRatio, ratio)
  }
  return(CandidaCall)
}

# DxCorrelation: Calculate correlation and visualise the correlation of binding profiles between two sample sets (e.g. test and reference). Input should be normalized binding values for all Candida probes, without positive/negative control probe binding values. Vectors containing the corresponding column sample names for test samples and reference samples are required. If you are missing either test or reference samples, please provide the same vector for both test and reference arguments. Optional: a vector specifying color values can be provided. ColVec will color correlation points according to reference sample used. Output will be a csv file of correlation values between samples, and a plot visualising these correlation values per sample.
#-------------------------------------------------------------

DxCorrelation <- function(NormalizedData, ReferenceID, TestID, ColVec=NULL){
  #NormalizedData should be a matrix of normalized Candida probe binding values, all other arguments should be sample name or color vectors.
  Norm <- NormalizedData[1:(nrow(NormalizedData)-2), ]
  Cor1 <- cor(Norm)
  write.csv(Cor1, "./CorrelationValues.csv")
  CorPlot <- Cor1[rownames(Cor1) %in% ReferenceID, colnames(Cor1) %in% TestID]
  CorPlot <- melt(CorPlot)
  if(!is.null(ColVec)){
    p <- ggplot(CorPlot, aes(x=Var2, y=value, color = Var1)) + geom_point() + theme_classic() + xlab("Sample") + ylab("Pearson Correlation") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + scale_color_manual(values = ColVec)
  }
  else{
    p <- ggplot(CorPlot, aes(x=Var2, y=value, color = Var1)) + geom_point() + theme_classic() + xlab("Sample") + ylab("Pearson Correlation") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  }
  return(p)
}