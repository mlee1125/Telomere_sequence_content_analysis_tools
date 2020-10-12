require(ggplot2)
require(randomForest)
require(e1071)
require(ROCR)
require(stringr)
require(RColorBrewer)
require(reshape2)

training_data <- read.csv('./training_data.csv', header=TRUE, sep=',', row.names = 1)

RF <- randomForest(TMM ~ ., data = training_data[,-15], replace = TRUE, mtry = 5, ntree = 500, proximity=TRUE, localImp=TRUE)
print(RF)
print(RF$importance)

plot(RF)
MDSplot(RF, training_data$TMM)

pred <- prediction(RF$votes[,2], training_data[,"TMM"])
perf <- performance(pred, "tpr", "fpr")
plot(perf, colorize=T)

ggplot(data.frame(TPR=perf@y.values[[1]], FPR=perf@x.values[[1]]), aes(FPR,TPR, colour='x'))+
  scale_color_manual(values='blue')+
  geom_step(size=1.5) +
  theme(
    panel.background = element_rect(colour='white', fill='white'), 
    panel.grid.major = element_line(colour='#525252'),
    axis.text = element_text(size=rel(2)),
    axis.title = element_text(size=rel(2))
  ) +
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))


write.csv(RF$votes, file = "training_data_OOB_testing_votes.csv")
