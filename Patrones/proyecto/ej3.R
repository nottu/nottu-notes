library("ada")
library("rpart")

spambase <- read.table("spambase.data", sep = ",")
names(spambase)[length(names(spambase))] <- "Y" #cambia nombre del ultimo valor
spambase$Y <- factor(spambase$Y)

train_pct <- 0.75
n_train   <- floor(train_pct * nrow(spambase))
rand_ind  <- sample(nrow(spambase))

train <- spambase[head(rand_ind, n=n_train),]
test  <- spambase[tail(rand_ind, n=nrow(spambase) - n_train),]


# Forest Stuff
trainForest <- function(data, n_trees){
  n_pred <- trunc(sqrt(ncol(data)))

  forest <- vector(mode="list", length=n_trees)
  for (i in 1:n_trees) {
    p <- names(train)[sample(ncol(data) - 1, n_pred)] #no usamos la columna Y
    f <- paste("Y ~ ", paste(p, collapse = "+"))
    train_set <- train[sample(nrow(train), replace = T), ]
    fit <- rpart(as.formula(f), method='class', data=train_set)
    forest[[i]] <- fit
  }
  return(forest)
}

testForest <- function(forest, test, n_trees){
  res <- rep(0, nrow(test))
  for (i in 1:nrow(test)) {
    val <- 0
    for (j in 1:n_trees) {
        if(predict(forest[[j]], newdata = test[i,], type = "class") == 1) val <- val +1 
    }
    if(val >= n_trees/2) res[i] = 1
  }
  return(res)
}

forest <- trainForest(test, 51)
res <- testForest(forest, test, 51)

table(res, test$Y)

#boost stuff
control <- rpart.control(cp = -1, maxdepth = 10, maxcompete = 1, xval = 0)
boost <- ada(Y ~ ., data = train, type = "discrete", control = control, iter = 50)
res_boost <- predict(boost, newdata =test)
table(res_boost, test$Y)


reglog = glm(Y ~., data=train, family='binomial')
reglog.pred <- round(predict(reglog, newdata=test, type='response'))

table(reglog.pred, test$Y)