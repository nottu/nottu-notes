library("gdata")
require(kohonen)
# Load Data
data = read.xls('data.xlsx', sheet=1)
data2016 = data[data[, 'year'] == 2016,][-1][-1][-1]
rownames(data2016) <- data[data[, 'year'] == 2016,][,1]
for(i in 1:ncol(data2016)){
  data2016[is.na(data2016[,i]), i] <- mean(data2016[,i], na.rm = TRUE)
}
data2016 <- data2016[, colSums(is.na(data2016)) != nrow(data2016)]
# PCA
p <- prcomp(data2016, scale=TRUE) #remove country, country, year
# plot(p$x[,1], p$x[,2], ylab="PC2", xlab='PC1')
biplot(p, cex=0.5)
px = p$x[order(p$x[,1], decreasing=TRUE),]
par(oma=c(2,2,2,2))
plot(head(px[,1],10), type='l', xaxt='n', ann=FALSE)
axis(1,at=1:10, 
  labels=head(rownames(px), 10), las=2, cex=0.2)

data2016s <- scale(data2016)
# SOM
som1 <- som(data2016s, grid = somgrid(xdim = 11, ydim=11, topo="hexagonal"))
# plot(som1, type="dist.neighbours")
# plot(som1, type="codes")
plot(som1, type = "mapping", main = "Mapping Type SOM", labels=rownames(data2016s),data=data2016s, cex=0.5, font=1)
# ISOMAP
# dis <- vegdist(data2016s)

# ord <- isomap(dis, k = 3)
# pl <- plot(ord, main="isomap k=3", pch=rownames(data2016s))
dis <- dist(data2016s, method = "euclidean")
ord <- isomap(dis, k = 3)
pl <- plot(ord, main="isomap k=3", pch='.')
a <- data.frame(pl$sites)
text(a[,1], a[,2], labels=rownames(data2016s), cex=0.5)
