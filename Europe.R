rm(list=ls())
nfa <- read.csv("NFA 2018.csv", header=TRUE, sep=",")
nfa_con <- nfa[which(nfa$record=="EFConsPerCap"),]
nfa14 <- nfa_con[which(nfa_con$year==2014),]
ef <- nfa14[which(nfa14$UN_region=="Europe"),]
ef <- ef[,-c(3,5,6)]
colnames(ef) <- c("Country", "ISO_code","Subregion", "crop", "grazing", "forest", "fishing", 
                  "built_up", "carbon", "total", "GDPPC", "population")
rownames(ef) <- c(1:nrow(ef))


#################### Descriptive
v <- c(4:9,12)

#par(cex= 1.5, cex.axis=1.5, cex.lab=1.5, cex.main=1.5)

# Basic visualization
panel.boxplot <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  bx <- boxplot.stats(x, do.conf = FALSE)
  arrows(bx$stats[1], 0.5, bx$stats[5], 0.5, angle=90, length=0.1, code=3)
  rect(bx$stats[2], 0.3, bx$stats[4], 0.7, ...)
  segments(bx$stats[3], 0.3, y1=0.7, lwd=3, lend=1)
  if(n <- length(bx$out))
    points(bx$out, rep(0.5, n), pch=1, ...)
  box()
}
pairs(ef[,v], diag.panel=panel.boxplot, panel=function(x,y,...){
  points(x, y, ...)
  abline(lm(y ~ x), col="blue")
})

summary(ef)

round(cor(ef[,v]), 4)

# Outliers detection
outlier <- order(ef$total)[38]
panel.boxplot2 <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  bx <- boxplot.stats(x, do.conf = FALSE)
  arrows(bx$stats[1], 0.5, bx$stats[5], 0.5, angle=90, length=0.1, code=3)
  rect(bx$stats[2], 0.3, bx$stats[4], 0.7, ...)
  segments(bx$stats[3], 0.3, y1=0.7, lwd=3, lend=1)
  if (n <- length(bx$out)){
    cols <- rep("black", length(bx$out))
    cols[which(bx$out == x[outlier])] <- "red"
    points(bx$out, rep(0.5, n), pch=1, col=cols, ...)}
  box()
}
cols <- rep("black", nrow(ef))
cols[outlier] <- "red"
pairs(ef[,v], diag.panel=panel.boxplot2, panel=function(x,y,...){
  points(x, y, col=cols, ...)
  abline(lm(y ~ x), col="blue")
})

pairs(ef[-outlier,v], diag.panel=panel.boxplot, panel=function(x,y,...){
  points(x, y, ...)
  abline(lm(y ~ x), col="blue")
})

region <- c("Southern", "Western","Northern","Eastern")
cformat <- c("black", "red", "forestgreen", "blue")
sformat <- c(1:4)
#sformat <- c("S", "W", "N", "E")
pt.fn <- function(x, y,...) {
  for (i in 1:4) {
    subregion <- paste(region[i], "Europe")  
    vec <- which(ef$Subregion==subregion)
    points(x[vec], y[vec], col=cformat[i], pch=sformat[i], cex=1.5)
  }
}
legend.fn <- function(position, ...) {
  legend(position, legend=region, col=cformat, pch=sformat, ...)
}


########## PCA
# pca summary
pca <- princomp(ef[,v], cor=TRUE)
summary(pca, loadings=TRUE)

# scree plot
eigvals <- pca$sdev^2
k <- length(eigvals)
plot(1:k,eigvals,type="b",xlab="i",ylab=expression(lambda[i]), main="Scree plot", cex=1.5)

# pca scores
library(MVA)
r <- range(pca$scores[,1:4])
pairs(pca$scores[,1:4],xlim=r,ylim=r,
      panel=function(x,y,...){
        text(x,y,ef$ISO_code,cex=1,col="red")
        bvbox(cbind(x,y),add=TRUE)
      })

plot(pca$scores[,1:2], type="n", xlab="PC1", ylab="PC2")
pt.fn(pca$scores[,1], pca$scores[,2])
legend.fn("topright")

# pca loadings
plot(pca$loadings[,1:2], type="n", xlab="PC1", ylab="PC2", cex=1.5)
text(pca$loadings[,1], pca$loadings[,2], labels= colnames(ef)[v], cex=1.5)

# biplot
biplot(pca, xlabs=ef$ISO_code, xlab="PC1", ylab="PC2", cex=1.5)

#################### Multiple lm
par(mfrow=c(3,3))
for(i in v){
  plot(ef$GDPPC ~ ef[,i],ylab="per capita GDP",xlab=colnames(ef)[i])
  abline(lm(ef$GDPPC ~ ef[,i]),col="red")
}

summary(lm(ef$GDPPC ~ pca$scores))

mfrow(par=c(1,2))
plot(ef$GDPPC ~ pca$scores[,1],type="n", xlab=paste("PC1"),
     ylab="per capita GDP")
abline(lm(ef$GDPPC ~ pca$scores[,1]),col="red")
text(pca$scores[,1], ef$GDPPC, labels=ef$ISO_code)


par(mfrow=c(1,1))

################### CCA
x1 <- scale(ef$crop)
x2 <- scale(ef$grazing+ef$fishing)
y1 <- scale(ef$forest+ef$carbon)
y2 <- scale(ef$built_up)
round(cor(cbind(x1,x2,y1,y2)), 4)
cormat <- cor(cbind(x1,x2,y1,y2))
r11 <- cormat[1:2,1:2]
r12 <- cormat[1:2,3:4]
r22 <- cormat[3:4,3:4]
r21 <- t(r12)
e1 <- solve(r11) %*% r12 %*% solve(r22) %*% r21
e2 <- solve(r22) %*% r21 %*% solve(r11) %*% r12
eigen(e1)
eigen(e2)

x <- cbind(x1,x2)
y <- cbind(y1,y2)
u <- x %*% eigen(e1)$vectors
cor(u,x)
v <- y %*% eigen(e2)$vectors
cor(v,y)

cc <- sqrt(eigen(e1)$values)
library(CCP)
p.asym(cc, 38, 2, 2, tstat = "Pillai")

plot(u[,1],v[,1],xlab=expression(U[1]),ylab=expression(V[1]), type="n")
text(u[,1],v[,1], labels=ef$ISO_code)
abline(0,0, lty="dashed", col="red")
abline(v=0, lty="dashed", col="red")
