library("rpart")

data <- read.csv("datos.data")
names(data)[length(names(data))] <- "Y" #cambia nombre del ultimo valor
data$Y <- factor(data$Y)

train_pct <- 0.75
n_train <- floor(train_pct * nrow(data))
rand_ind <- sample(nrow(data))

train <- data[head(rand_ind, n=n_train),]
test <- data[tail(rand_ind, n=nrow(data) - n_train),]

#entrenamos arboles
n_trees <- 51
n_pred <- trunc(sqrt(ncol(data)))

forest <- vector(mode="list", length=n_trees)
for (i in 1:n_trees) {
  p <- names(train)[sample(ncol(data) - 1, n_pred)] #no usamos la columna Y
  f <- paste("Y ~ ", paste(p, collapse = "+"))
  train_set <- train[sample(nrow(train), replace = T), ]
  fit <- rpart(as.formula(f), method='class', data=train_set)
  forest[[i]] <- fit
}

res <- rep(0, nrow(test))
for (i in 1:nrow(test)) {
  val <- 0
  for (j in 1:n_trees) {
      if(predict(forest[[j]], newdata = test[i,], type = "class") == 1) val <- val +1 
  }
  if(val >= n_trees/2) res[i] = 1
}

table(res, test$Y)


#boosting
control <- rpart.control(cp = -1, maxdepth = 10, maxcompete = 1, xval = 0)
boost <- ada(Y ~ ., data = train, type = "discrete", control = control, iter = 50)
res_boost <- predict(boost, newdata =test)
table(res_boost, test$Y)

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

testForest <- function(forest, test){
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