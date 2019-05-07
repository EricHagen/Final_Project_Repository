####CHUNKED#######

library(ape)
library(plyr)
library(psych)
setwd("~/Desktop/STEP_5_DATE_TREES_CHUNKED")
shit <- c(13, 14, 16, 22, 30, 31, 40, 45, 48)
net.div.bd <- c()
eps <-c() 
for(i in 1:50){
  if(i %in% shit){
    next
  }
  else if(i == 34){
    next
  }
  else{
    filename = paste("tree.dated.", i, sep="")
    bd.trees <- read.tree(filename)
    net.div.bd <- c(net.div.bd, birthdeath(bd.trees)$par[2])
    eps <- c(eps, birthdeath(bd.trees)$par[1])
  }
}

#print(net.div.bd)
#print(eps)

m_div <- mean(net.div.bd)
m_div_df <- data.frame(m_div)
s_div <- sd(net.div.bd)
m_eps <- mean(eps)
s_eps <- sd(eps)
st_err_netdivbd <- (s_div / sqrt(length(net.div.bd)))
sterrnet <- data.frame(st_err_netdivbd)
st_err_eps <- (s_eps / sqrt(length(eps)))

#GET TRUE ABLINE MEAN FOR DIV & EPS BY RUNNING BIRTHDEATH FUNCTION ON dated.matk.tre FROM STEP 1
bd.matk <- read.tree("true.tre")
matk.div.bd <- birthdeath(bd.matk)$par[2]
matk.eps <- birthdeath(bd.matk)$par[1]

#WRITE THE ERROR.BAR FUNCTION (from Dave Tang: https://davetang.org/muse/2014/06/25/plotting-error-bars-with-r/)
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

#DIVERSIFICATION CURVE
barcenters <- hist(net.div.bd, prob = TRUE, col = 'blue', ylim = c(0,95))
d = density(net.div.bd)
error.bar(net.div.bd, m_div_df, sterrnet)
#segments(net.div.bd, m_div - st_err_netdivbd * 2, net.div.bd, m_div + st_err_netdivbd * 2, lwd = 1.5)
#arrows(net.div.bd, m_div - st_err_netdivbd * 2, net.div.bd, m_div + st_err_netdivbd * 2, lwd = 1.5, code = 3, length = 0.05)
lines(d, col = "orange", lwd = 3)
t.test(net.div.bd, conf.level = .95)
abline(v = m_div, col = "red", lwd = 3, lty = "dotted")
abline(v = matk.div.bd, col = "purple", lwd = 3, lty = "dotted")

#EXTINCTION CURVE
hist(eps, prob = TRUE, col = 'red', ylim = c(0, 8))
dens <- density(eps)
lines(dens, col = "orange", lwd = 3)
abline(v = m_eps, col = "blue", lwd = 3, lty = "dotted")
abline(v = matk.eps, col = "purple", lwd = 3, lty = "dotted")
