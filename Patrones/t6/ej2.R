library("tm")
library("SnowballC")
library("e1071")

alldir=DirSource('./texts', encoding = "UTF-8", recursive=TRUE)
news <- Corpus(alldir, readerControl=list(reader=readPlain,language="en"))
dtm <- DocumentTermMatrix(news, control=list(removePunctuation=TRUE, 
                        tolower=TRUE, stopwords=c("english"), stripWhitespace=TRUE, 
                        minWordLength=2))
d <- DocumentTermMatrix(news, list(dictionary=findFreqTerms(dtm, 50)))

classvec <- vector()
# código para categorías de 
# https://github.com/chenmiao/Big_Data_Analytics_Web_Text
for (filedir in alldir$filelist) {
  classlabel=basename(dirname(filedir))
  classvec=c(classvec,classlabel)
}
summary(classvec)
classvec <- factor(classvec)

#separar conjunto de pruebas
rand_idxs <- sample(nrow(d))
n_test = floor(nrow(d) * 0.75)

train <- d[head(rand_idxs, n=n_test), ]
test <- d[tail(rand_idxs, n=(nrow(d)-n_test)), ]

train_class <- classvec[head(rand_idxs, n=n_test)]
test_class <-classvec[tail(rand_idxs, n=(nrow(d)-n_test))]

svm_model <- svm(train, train_class, kernel="linear")

pred <- predict(svm_model, test)
table(pred,test_class)
