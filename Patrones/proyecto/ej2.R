
pm10 = read.table('pm10.txt', header=TRUE)
pm10$highpm10 <- as.factor(pm10$highpm10)

train_pct <- 0.75
n_train <- floor(train_pct * nrow(pm10))
rand_ind <- sample(nrow(pm10))

train <- pm10[head(rand_ind, n=n_train),]
test <- pm10[tail(rand_ind, n=nrow(pm10) - n_train),]

getCorrectPreds <- function(model, testData){
  testData.pred <- predict(model, newdata=testData, type='response')
  correct = 0
  for (i in 1:length(testData.pred)){
    if(round(testData.pred[i]) == testData$highpm10[i]) correct = correct + 1
  }
  return( correct/length(testData.pred) )
}

f1 <- 'highpm10 ~ cars+temp2m+winddirection+time'
f2 <- 'highpm10 ~ cars+temp2m'
f3 <- 'highpm10 ~ cars*temp2m'
pm10.fit1 = glm(as.formula(f1), data=train, family='binomial')
# summary(pm10.fit1)
# coef(pm10.fit1)
correct.percent = getCorrectPreds(pm10.fit1, test)
print(paste(c("Correctly predicted :", correct.percent), collapse=" "))

pm10.fit2 = glm(as.formula(f2), data=train, family='binomial')
# summary(pm10.fit2)
# coef(pm10.fit2)
correct.percent = getCorrectPreds(pm10.fit2, test)
print(paste(c("Correctly predicted :", correct.percent), collapse=" "))

pm10.fit3 = glm(as.formula(f3), data=train, family='binomial')
# summary(pm10.fit3)
# coef(pm10.fit3)
correct.percent = getCorrectPreds(pm10.fit3, test)
print(paste(c("Correctly predicted :", correct.percent), collapse=" "))

# Neural Nets
library('nnet')
neuronal1 <- nnet(as.formula(f1), size=10, data=train)
res <- predict(neuronal1, test)
correct = 0
for (i in 1:length(res)){
  if(round(res[i]) == test$highpm10[i]) correct = correct + 1
}
correct.percent <- correct/length(res)
print(paste(c("Correctly predicted :", correct.percent), collapse=" "))

neuronal2 <- nnet(as.formula(f2), size=10, data=train)
res <- predict(neuronal2, test)
correct = 0
for (i in 1:length(res)){
  if(round(res[i]) == test$highpm10[i]) correct = correct + 1
}
correct.percent <- correct/length(res)
print(paste(c("Correctly predicted :", correct.percent), collapse=" "))

neuronal3 <- nnet(as.formula(f3), size=10, data=train)
res <- predict(neuronal3, test)
correct = 0
for (i in 1:length(res)){
  if(round(res[i]) == test$highpm10[i]) correct = correct + 1
}
correct.percent <- correct/length(res)
print(paste(c("Correctly predicted :", correct.percent), collapse=" "))

