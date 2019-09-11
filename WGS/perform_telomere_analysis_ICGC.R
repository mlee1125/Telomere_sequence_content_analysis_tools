#!/usr/bin/env Rscript

require(dplyr)
args = commandArgs(trailingOnly=TRUE)
seqErrCSV = "./telo_seq_error.csv"

normalit <- function(m){(m/sum(m))}

# test if there are two areguments: if not, return an error
if (length(args)==2) {
  seqError <- data.frame(read.csv(seqErrCSV))
  inputCSV <- data.frame(read.csv(args[1]))
  inputCSV[,2] <- factor(inputCSV[,2])
  inputCSV <- inputCSV[inputCSV[,"Sample_Type"] == 2,]
  
  strandSums <- aggregate(inputCSV[,"Count"], by=as.list(inputCSV[,c(1:2,4)]), FUN=sum)
  strandProp <- data.frame(strandSums %>% group_by(Sample_ID, Sample_Type) %>% mutate(x = normalit(x)))
  
  outputCSV <- data.frame(inputCSV %>% group_by(Sample_ID, Sample_Type, Strand) %>% mutate(Count = normalit(Count)))
  outputCSV[outputCSV[,4]=="C" & outputCSV[,"Repeat"] != "TTAGGG", "Count"] - seqError[seqError[,1]=="C",3] -> outputCSV[outputCSV[,4]=="C" & outputCSV[,"Repeat"] != "TTAGGG", "Count"]
  outputCSV[outputCSV[,4]=="G" & outputCSV[,"Repeat"] != "TTAGGG", "Count"] - seqError[seqError[,1]=="G",3] -> outputCSV[outputCSV[,4]=="G" & outputCSV[,"Repeat"] != "TTAGGG", "Count"]
  
  for (x in unique(strandSums[strandSums[,"x"] < 1000, "Sample_ID"])){
    outputCSV <- outputCSV[outputCSV[,"Sample_ID"] != x,]
    strandProp <- strandProp[strandProp[,"Sample_ID"] != x,]
    strandSums <- strandSums[strandSums[,"Sample_ID"] != x,]
  }
  for (x in unique(outputCSV[is.na(outputCSV[,"Count"]), "Sample_ID"])){
    outputCSV <- outputCSV[outputCSV[,"Sample_ID"] != x,]
    strandProp <- strandProp[strandProp[,"Sample_ID"] != x,]
    strandSums <- strandSums[strandSums[,"Sample_ID"] != x,]
  }
  outputCSV[outputCSV[,"Count"] < 0, "Count"] <- 0
  outputCSV[order(outputCSV[,"Repeat"], outputCSV[,"Strand"], outputCSV[,"Sample_ID"]), "Count"] <- outputCSV[order(outputCSV[,"Repeat"], outputCSV[,"Strand"], outputCSV[,"Sample_ID"]), "Count"] * strandProp[order(strandProp[,"Strand"], strandProp[,"Sample_ID"]), "x"]
  outputCSV <- data.frame(aggregate(outputCSV[,"Count"], by=as.list(outputCSV[,c(1:2,3)]), FUN=sum))
  
  outputCSV2 <- data.frame(TTAGGG=outputCSV[outputCSV[,"Repeat"] == "TTAGGG", "x"])
  for (x in c("ATAGGG", "CTAGGG", "GTAGGG", "TAAGGG", "TCAGGG", "TGAGGG", "TTCGGG", "TTGGGG","TTTGGG", "TTAAGG", "TTACGG", "TTATGG", "TTAGAG", "TTAGCG", "TTAGTG", "TTAGGA", "TTAGGC", "TTAGGT")){
    outputCSV2 <- data.frame(outputCSV2, x=outputCSV[outputCSV[,"Repeat"] == x, "x"])
  }
  colnames(outputCSV2) <- c("TTAGGG","ATAGGG", "CTAGGG", "GTAGGG", "TAAGGG", "TCAGGG", "TGAGGG", "TTCGGG", "TTGGGG","TTTGGG", "TTAAGG", "TTACGG", "TTATGG", "TTAGAG", "TTAGCG", "TTAGTG", "TTAGGA", "TTAGGC", "TTAGGT")
  rownames(outputCSV2) <- outputCSV[outputCSV[,"Repeat"] == "TTAGGG", "Sample_ID"]
  outputCSV2[,"TTAGGG"] <- (1 - apply(outputCSV2[,-1], 1, sum))
  write.csv(outputCSV2,file = args[2])
} else {
  stop("Must supply input file and output file.n", call.=FALSE)
}

