####UNCHUNKED#######
library(ape)
library(plotrix)
setwd("/Volumes/KINGSTON/matK_Sim_Files/STEP_5_DATE_TREES")
shit <- c(6, 11, 14, 15, 17, 31)
net.div.bd <- c()
eps <- c()
for(i in 1:50){
  if(i %in% shit){
    next
  }
  else if(i == 34){
    next
  }
  else{
    filename = paste("tree.dated.", i, sep="")
    bd.tree <- read.tree(filename)
    net.div.bd <- c(net.div.bd, birthdeath(bd.tree)$par[2])
    eps <- c(eps, birthdeath(bd.tree)$par[1])
  }
}

print(net.div.bd)
print(eps)
m_div <- mean(net.div.bd)
s_div <- sd(net.div.bd)
m_eps <- mean(eps)
s_eps <- sd(eps)

#GET TRUE ABLINE MEAN FOR DIV & EPS BY RUNNING BIRTHDEATH FUNCTION ON dated.matk.tre FROM STEP 1
bd.matk <- read.tree("true.tre")
matk.div.bd <- birthdeath(bd.matk)$par[2]
matk.eps <- birthdeath(bd.matk)$par[1]

#DIVERSIFICATION CURVE
hist(net.div.bd, prob = TRUE, col = 'blue', ylim = c(0,65))
d = density(net.div.bd)
lines(d, col = "green", lwd = 3)
abline(v = m_div, col = "red", lwd = 3, lty = "dotted") #Unchunked
abline(v = 0.1009396, col = "orange", lwd = 3, lty = "dotted") #Chunked
ablineclip(b = 0, h = 55, x1 = 0.0979002, x2 = 0.1039790, lty = 1, col = "orange")
points(x = 0.1009396, y = 55, type = "p", pch = 19, col = "black")
abline(v = matk.div.bd, col = "purple", lwd = 3, lty = "dotted") #True tree
t.test(net.div.bd, conf.level = .95)
ablineclip(b = 0, h = 50, x1 = 0.09828769, x2 = 0.10764976, lty = 1, col = "red")
points(x = 0.1029687, y = 50, type = "p", pch = 19, col = "black")

#EXTINCTION CURVE
hist(eps, prob = TRUE, col = 'red', ylim = c(0, 8))
dens <- density(eps)
lines(dens, col = "orange", lwd = 3)
abline(v = m_eps, col = "blue", lwd = 3, lty = "dotted")
abline(v = matk.eps, col = "purple", lwd = 3, lty = "dotted")
t.test(eps, conf.level = .95)
ablineclip(b = 0, h = 6, x1 = 0.09963284, x2 = 0.28454725, lty = 1, col = "blue")
points(x = 0.19209, y = 6, type = "p", pch = 19, col = "black")
