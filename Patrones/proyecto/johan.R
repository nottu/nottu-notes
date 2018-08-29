library("gmapsdistance")
key1 <- 'AIzaSyCrHhon55lhJLIOhHRgOmlYJR7ZQEJWK-U'
key2 <- 'AIzaSyCHCNx_vob2bQEXEDzqWyZyOYBdwGkZY6o'
# key1 <- 'AIzaSyAV05A6COKknh5SYbRs3Z8utJ7arCwYYng'

# keys <- c(key1, key2)
key <- set.api.key(key2)
mins <- 60
hour <- 60 * mins
sleepTime <- 5 * mins # in seconds

# returns time since epoch
getCurrentTime <- function(){
  return ( as.numeric(Sys.time()) );
}
# returnes time since param
getTimeSince <- function(time){ 
  return (as.numeric(Sys.time()) - time);
}

places_north <- list(americas = '19.5188850+-99.1723660', 
                     soumaya  = '19.4395855+-99.2036083',
                     arena    = '19.4402330+-99.2032650',
                     cosmopol = '19.6044867+-99.2298526',
                     aragon   = '19.4508515+-99.1143098')

places_south <- list(frida      = '19.3552163+-99.1612616',
                    anahuacalli = '19.3259224+-99.1533190',
                    anahuac     = '19.3060979+-99.0995032',
                    sixflags    = '19.2984631+-99.2135507',
                    olmedo      = '19.2668719+-99.1277232')
n_places <- length(places_north)

getData <- function(){
  current_time = getCurrentTime()
  nts <- gmapsdistance(unlist(places_north), unlist(places_south), mode="driving", key=key)
  stn <- gmapsdistance(unlist(places_south), unlist(places_north), mode="driving", key=key)

  vec <- c()
  for (i in 1: n_places){
    for (j in 1: n_places){
      vec[(i - 1) * n_places + j] = nts$Time[,-1][i,j]
    }
  }
  vec[26] = current_time
  vec[27] = -1 #category
  df <- data.frame(t(matrix(vec)))
  vec <- c()
  for (i in 1: n_places){
    for (j in 1: n_places){
      vec[(i - 1) * n_places + j] = stn$Time[,-1][j,i] #data is transposed...
    }
  }
  vec[26] = current_time
  vec[27] = 1 #category
  mtx <- t(matrix(vec))
  df <- rbind(df, data.frame( mtx ))
  return(df)
}

getSamples <-function(n_samples){
  ctime <- getCurrentTime()
  df <- getData()
  print("data sampled")
  getTimeSince(ctime)
  if(n_samples == 1) return(df)
  for (i in 1:n_samples) {
    print("new sample")
    # key <- set.api.key(keys[(1 + i %% 2)]) #alternate between keys...
    Sys.sleep(sleepTime - getTimeSince(ctime))
    ctime <- getCurrentTime()
    df <- rbind(df, getData())
  }
  return(df)
}
ctime = getCurrentTime()
print("Sleeping now")
Sys.sleep((hour * 1) - getTimeSince(ctime))
print("Awake, working.... ")
df <-getSamples( (4 * 60)/10 ) #sample for 4hrs every 1 minuts
write.csv(df, file = "north_south_data.csv") #save data in case of disaster...

# PCA
names(df)[length(names(df))] <- "Y" #cambia nombre del ultimo valor
df.dat = df[,-length(df)]
p <-prcomp((df.dat), scale = TRUE)
plot(p$x[,1], p$x[,2], col=ddf$Y+2) #plot with colors for categories

# svm classifier
svm_model <- svm(df.dat, df$Y, kernel="linear")
pred <- predict(svm_model, df.dat)
table(pred > 0, df$Y)

svm_model_r <- svm(df.dat, df$Y)