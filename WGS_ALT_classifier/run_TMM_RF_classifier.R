require(ggplot2)
require(randomForest)
require(e1071)
require(ROCR)
require(stringr)
require(RColorBrewer)
require(reshape2)

# Data required to be the same format and column names as the training_data.csv table

input_data <- read.csv('./TCGA_telomere_counts_-_RF_input.csv', header=TRUE, sep=',', row.names = 1)

results <- predict(RF, input_data[,c('TTAGGG','ATAGGG','CTAGGG','GTAGGG','TAAGGG','TCAGGG','TGAGGG','TTCGGG','TTGGGG','TTTGGG','TTAAGG','TTACGG','TTATGG','TTAGAG','TTAGCG','TTAGTG','TTAGGA','TTAGGC','TTAGGT','TL_TN')], type = "prob")

write.csv(results, file = "./output_TMM_predictions.csv")
