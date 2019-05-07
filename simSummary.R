#
# n number of simulations
# file.name the name of the file to be used as an output 
# analysis.type = c("all", "root", "clades", "diversification")
# true.tree

library(ape)
library(paleotree)


SummarizeCalibrationTest <- function(true.tree, true.data, n=100, file.name="all.calibrations") {
    
    root.age <- c()
    net.div <- c()
    extinct.frac <- c()
    transition.edge.lengths <- c()
    woody.edge.lengths <- c()
    herb.edge.lengths <- c()
    
    for(tree.index in 1:100){
        #        if(tree.index != 42 & tree.index != 82){
            print(tree.index)
            tree.tmp <- read.tree(paste("tree.dated.", tree.index, sep=""))
            root.age <- c(root.age, max(branching.times(tree.tmp)))
            net.div <- c(net.div,  birthdeath(tree.tmp)$par[2])
            extinct.frac <- c(extinct.frac,  birthdeath(tree.tmp)$par[1])
            edge.length.mat <- EdgeLengthByState(tree.tmp, true.data)
            transition.edge.lengths <- rbind(transition.edge.lengths, as.numeric(as.character(edge.length.mat[which(edge.length.mat[,4] != edge.length.mat[,5]),3])))
            woody.edge.lengths <- rbind(woody.edge.lengths, as.numeric(as.character(edge.length.mat[which(edge.length.mat[,4]==1 & edge.length.mat[,5]==1),3])))
            herb.edge.lengths <- rbind(herb.edge.lengths, as.numeric(as.character(edge.length.mat[which(edge.length.mat[,4]==2 & edge.length.mat[,5]==2),3])))
            # }else{
            #root.age <- c(root.age, NA)
            #net.div <- c(net.div,  NA)
            #extinct.frac <- c(extinct.frac, NA)
            #edge.length.mat <- NA
            #transition.edge.lengths <- rbind(transition.edge.lengths, NA)
            #woody.edge.lengths <- rbind(woody.edge.lengths, NA)
            #herb.edge.lengths <- rbind(herb.edge.lengths, NA)
            #}
    }
    save(root.age, net.div, extinct.frac, transition.edge.lengths, woody.edge.lengths, herb.edge.lengths, file=file.name)
}




EdgeLengthByState <- function(phy, data){
    
    data <- data.frame(data[,2], data[,2], row.names=data[,1])
    data <- data[phy$tip.label,]
    
    tip.states <- data[,2]
    node.states <- phy$node.label
    total.state.vector <- c(tip.states, node.states)
    
    edge.mat <- cbind(phy$edge, phy$edge.length)
    begin.state <- c()
    end.state <- c()
    for(node.index in 1:dim(edge.mat)[1]){
        begin.state <- c(begin.state, total.state.vector[edge.mat[node.index,1]])
        end.state <- c(end.state, total.state.vector[edge.mat[node.index,2]])
    }
    
    edge.mat <- cbind(edge.mat, begin.state, end.state)
    final.mat <- data.frame(anc=edge.mat[,1], desc=edge.mat[,2], edge.length=edge.mat[,3], begin.state=edge.mat[,4], end.state=edge.mat[,5])
    return(final.mat)
}



#SummarizeRatedTrees <- function(file.name){
#    branch.var <- c()
#    sum.branch <- c()
#    max.length <- c()
#    for(index in 1:100){
#        tree.tmp <- read.tree(paste("camp.rated.", index, ".tre", sep=""))
#        branch.var <- c(branch.var, sd(tree.tmp$edge.length))
#        sum.branch <- c(sum.branch, sum(tree.tmp$edge.length))
#        max.length <- c(max.length, max(dateNodes(tree.tmp)))
#    }
#   return(
#}


starting.point.tree <- function(phy, yule=FALSE) {
    p.yule <- c(yule(phy)$lambda, 0)
    if(yule){
        p.yule
    }else{
        suppressWarnings(c(birthdeath(phy)$para[2] / (1-birthdeath(phy)$para[1]), ((birthdeath(phy)$para[1] * birthdeath(phy)$para[2]) / (1-birthdeath(phy)$para[1]))))
    }
}


starting.point.generator <- function(phy, k, samp.freq.tree, q.div=5, yule=FALSE) {
    pars.bd <- suppressWarnings(starting.point.tree(phy, yule))
    #Rescale parameters to account for sampling, if necessary, using Stadler 2013:
    pars.bd[1] = pars.bd[1] / samp.freq.tree
    pars.bd[2] = pars.bd[2] - (pars.bd[1]*samp.freq.tree) * (1 - 1/samp.freq.tree)
    r <- if  ( pars.bd[1] > pars.bd[2] )
    (pars.bd[1] - pars.bd[2]) else pars.bd[1]
    p <- rep(c(pars.bd, r / q.div), c(k, k, k * (k-1)))
    names(p) <- NULL
    p
}


true.tree <- read.tree("camp.matk.habit.recon.asterwood.carrotwood.tre")
true.data <- read.delim("camp.matk.habit.trait.asterwood.carrotwood.txt")
#This generates the summary files -- commented out for now.
#SummarizeCalibrationTest(true.tree, true.data, file.name="all3xSHIFTS")

true.root.age <- max(branching.times(true.tree))
true.net.div <- birthdeath(true.tree)$par[2]
true.extinct.frac <- birthdeath(true.tree)$par[1]
true.edge.mat <- EdgeLengthByState(true.tree, true.data)
true.transition.edge.lengths <- as.numeric(as.character(true.edge.mat[which(true.edge.mat[,4] != true.edge.mat[,5]),3]))
true.woody.edge.lengths <- as.numeric(as.character(true.edge.mat[which(true.edge.mat[,4]==1 & true.edge.mat[,5]==1),3]))
true.herb.edge.lengths <- as.numeric(as.character(true.edge.mat[which(true.edge.mat[,4]==2 & true.edge.mat[,5]==2),3]))


######################################################################################################################################
######################################################################################################################################
### Summarize results
######################################################################################################################################
######################################################################################################################################

load("allNOSHIFTS")
pd.transition.resid <- c()
pd.woody.resid <- c()
pd.herb.resid <- c()
for(index in 1:100){
    pd.transition.resid <- rbind(pd.transition.resid, 100*((transition.edge.lengths[index,] - true.transition.edge.lengths)/true.transition.edge.lengths))
    pd.woody.resid <- rbind(pd.woody.resid, 100*((woody.edge.lengths[index,] - true.woody.edge.lengths)/true.woody.edge.lengths))
    pd.herb.resid <- rbind(pd.herb.resid, 100*((herb.edge.lengths[index,] - true.herb.edge.lengths)/true.herb.edge.lengths))
}
mean.pd.transition.resid <- apply(pd.transition.resid, 1, mean)
mean.pd.woody.resid <- apply(pd.woody.resid, 1, mean)
mean.pd.herb.resid <- apply(pd.herb.resid, 1, mean)

#Table1 -- sum stats
sum.stats <- c()
sum.stats <- rbind(sum.stats, c("all, no var", mean(mean.pd.woody.resid), mean(mean.pd.herb.resid), mean(root.age), mean(net.div), mean(extinct.frac)))

#regress
regress.root.age <- c()
fit1 <- lm(root.age~mean.pd.woody.resid)
aa1 <- summary(fit1)
fit2 <- lm(root.age~mean.pd.herb.resid)
aa2 <- summary(fit2)
regress.root.age <- rbind(regress.root.age, c("all, no var", aa1$coefficients[2], aa1$r.squared, aa2$coefficients[2], aa2$r.squared))

regress.ef <- c()
fit1 <- lm(extinct.frac~mean.pd.woody.resid)
aa1 <- summary(fit1)
fit2 <- lm(extinct.frac~mean.pd.herb.resid)
aa2 <- summary(fit2)
regress.ef <- rbind(regress.ef , c("all, no var", aa1$coefficients[2], aa1$r.squared, aa2$coefficients[2], aa2$r.squared))

regress.netdiv <- c()
fit1 <- lm(net.div~mean.pd.woody.resid)
aa1 <- summary(fit1)
fit2 <- lm(net.div~mean.pd.herb.resid)
aa2 <- summary(fit2)
regress.netdiv <- rbind(regress.netdiv, c("all, no var", aa1$coefficients[2], aa1$r.squared, aa2$coefficients[2], aa2$r.squared))


load("all3xSHIFTS")

pd.transition.resid <- c()
pd.woody.resid <- c()
pd.herb.resid <- c()
for(index in 1:100){
    pd.transition.resid <- rbind(pd.transition.resid, 100*((transition.edge.lengths[index,] - true.transition.edge.lengths)/true.transition.edge.lengths))
    pd.woody.resid <- rbind(pd.woody.resid, 100*((woody.edge.lengths[index,] - true.woody.edge.lengths)/true.woody.edge.lengths))
    pd.herb.resid <- rbind(pd.herb.resid, 100*((herb.edge.lengths[index,] - true.herb.edge.lengths)/true.herb.edge.lengths))
}
mean.pd.transition.resid2 <- apply(pd.transition.resid, 1, mean)
mean.pd.woody.resid2 <- apply(pd.woody.resid, 1, mean)
mean.pd.herb.resid2 <- apply(pd.herb.resid, 1, mean)

sum.stats <- rbind(sum.stats, c("all, 3x", mean(mean.pd.woody.resid), mean(mean.pd.herb.resid), mean(root.age), mean(net.div), mean(extinct.frac)))

fit1 <- lm(root.age~mean.pd.woody.resid)
aa1 <- summary(fit1)
fit2 <- lm(root.age~mean.pd.herb.resid)
aa2 <- summary(fit2)
regress.root.age <- rbind(regress.root.age , c("all, 3x", aa1$coefficients[2], aa1$r.squared, aa2$coefficients[2], aa2$r.squared))

fit1 <- lm(extinct.frac~mean.pd.woody.resid)
aa1 <- summary(fit1)
fit2 <- lm(extinct.frac~mean.pd.herb.resid)
aa2 <- summary(fit2)
regress.ef <- rbind(regress.ef , c("all, 3x", aa1$coefficients[2], aa1$r.squared, aa2$coefficients[2], aa2$r.squared))

fit1 <- lm(net.div~mean.pd.woody.resid)
aa1 <- summary(fit1)
fit2 <- lm(net.div~mean.pd.herb.resid)
aa2 <- summary(fit2)
regress.netdiv <- rbind(regress.netdiv, c("all, 3x", aa1$coefficients[2], aa1$r.squared, aa2$coefficients[2], aa2$r.squared))


load("postKTNOSHIFTS")

pd.transition.resid <- c()
pd.woody.resid <- c()
pd.herb.resid <- c()
for(index in 1:100){
    pd.transition.resid <- rbind(pd.transition.resid, 100*((transition.edge.lengths[index,] - true.transition.edge.lengths)/true.transition.edge.lengths))
    pd.woody.resid <- rbind(pd.woody.resid, 100*((woody.edge.lengths[index,] - true.woody.edge.lengths)/true.woody.edge.lengths))
    pd.herb.resid <- rbind(pd.herb.resid, 100*((herb.edge.lengths[index,] - true.herb.edge.lengths)/true.herb.edge.lengths))
}
mean.pd.transition.resid <- apply(pd.transition.resid, 1, mean)
mean.pd.woody.resid <- apply(c(pd.woody.resid,pd.herb.resid), 1, mean)
mean.pd.herb.resid <- apply(pd.herb.resid, 1, mean)

sum.stats <- rbind(sum.stats, c("postKT, no var", mean(mean.pd.woody.resid), mean(mean.pd.herb.resid), mean(root.age), mean(net.div), mean(extinct.frac)))

#regress
fit1 <- lm(root.age~mean.pd.woody.resid)
aa1 <- summary(fit1)
fit2 <- lm(root.age~mean.pd.herb.resid)
aa2 <- summary(fit2)
regress.root.age <- rbind(regress.root.age , c("postKT, no var", aa1$coefficients[2], aa1$r.squared, aa2$coefficients[2], aa2$r.squared))

fit1 <- lm(extinct.frac~mean.pd.woody.resid)
aa1 <- summary(fit1)
fit2 <- lm(extinct.frac~mean.pd.herb.resid)
aa2 <- summary(fit2)
regress.ef <- rbind(regress.ef , c("postKT, no var", aa1$coefficients[2], aa1$r.squared, aa2$coefficients[2], aa2$r.squared))

fit1 <- lm(net.div~mean.pd.woody.resid)
aa1 <- summary(fit1)
fit2 <- lm(net.div~mean.pd.herb.resid)
aa2 <- summary(fit2)
regress.netdiv <- rbind(regress.netdiv, c("postKT, no var", aa1$coefficients[2], aa1$r.squared, aa2$coefficients[2], aa2$r.squared))


load("postKT3xSHIFTS")

pd.transition.resid <- c()
pd.woody.resid <- c()
pd.herb.resid <- c()
for(index in 1:100){
    pd.transition.resid <- rbind(pd.transition.resid, 100*((transition.edge.lengths[index,] - true.transition.edge.lengths)/true.transition.edge.lengths))
    pd.woody.resid <- rbind(pd.woody.resid, 100*((woody.edge.lengths[index,] - true.woody.edge.lengths)/true.woody.edge.lengths))
    pd.herb.resid <- rbind(pd.herb.resid, 100*((herb.edge.lengths[index,] - true.herb.edge.lengths)/true.herb.edge.lengths))
}
mean.pd.transition.resid <- apply(pd.transition.resid, 1, mean)
mean.pd.woody.resid <- apply(pd.woody.resid, 1, mean)
mean.pd.herb.resid <- apply(pd.herb.resid, 1, mean)

sum.stats <- rbind(sum.stats, c("postKT, 3x", mean(mean.pd.woody.resid), mean(mean.pd.herb.resid), mean(root.age), mean(net.div), mean(extinct.frac)))

fit1 <- lm(root.age~mean.pd.woody.resid)
aa1 <- summary(fit1)
fit2 <- lm(root.age~mean.pd.herb.resid)
aa2 <- summary(fit2)
regress.root.age <- rbind(regress.root.age , c("postKT, 3x", aa1$coefficients[2], aa1$r.squared, aa2$coefficients[2], aa2$r.squared))

fit1 <- lm(extinct.frac~mean.pd.woody.resid)
aa1 <- summary(fit1)
fit2 <- lm(extinct.frac~mean.pd.herb.resid)
aa2 <- summary(fit2)
regress.ef <- rbind(regress.ef , c("postKT, 3x", aa1$coefficients[2], aa1$r.squared, aa2$coefficients[2], aa2$r.squared))

fit1 <- lm(net.div~mean.pd.woody.resid)
aa1 <- summary(fit1)
fit2 <- lm(net.div~mean.pd.herb.resid)
aa2 <- summary(fit2)
regress.netdiv <- rbind(regress.netdiv, c("postKT, 3x", aa1$coefficients[2], aa1$r.squared, aa2$coefficients[2], aa2$r.squared))


load("postKT20NOSHIFTS")

pd.transition.resid <- c()
pd.woody.resid <- c()
pd.herb.resid <- c()
for(index in 1:100){
    pd.transition.resid <- rbind(pd.transition.resid, 100*((transition.edge.lengths[index,] - true.transition.edge.lengths)/true.transition.edge.lengths))
    pd.woody.resid <- rbind(pd.woody.resid, 100*((woody.edge.lengths[index,] - true.woody.edge.lengths)/true.woody.edge.lengths))
    pd.herb.resid <- rbind(pd.herb.resid, 100*((herb.edge.lengths[index,] - true.herb.edge.lengths)/true.herb.edge.lengths))
}
mean.pd.transition.resid <- apply(pd.transition.resid, 1, mean)
mean.pd.woody.resid <- apply(pd.woody.resid, 1, mean)
mean.pd.herb.resid <- apply(pd.herb.resid, 1, mean)

sum.stats <- rbind(sum.stats, c("postKT20, no var", mean(mean.pd.woody.resid), mean(mean.pd.herb.resid), mean(root.age), mean(net.div), mean(extinct.frac)))

#regress
fit1 <- lm(root.age~mean.pd.woody.resid)
aa1 <- summary(fit1)
fit2 <- lm(root.age~mean.pd.herb.resid)
aa2 <- summary(fit2)
regress.root.age <- rbind(regress.root.age , c("postKT20, no var", aa1$coefficients[2], aa1$r.squared, aa2$coefficients[2], aa2$r.squared))

fit1 <- lm(extinct.frac~mean.pd.woody.resid)
aa1 <- summary(fit1)
fit2 <- lm(extinct.frac~mean.pd.herb.resid)
aa2 <- summary(fit2)
regress.ef <- rbind(regress.ef , c("postKT20, no var", aa1$coefficients[2], aa1$r.squared, aa2$coefficients[2], aa2$r.squared))

fit1 <- lm(net.div~mean.pd.woody.resid)
aa1 <- summary(fit1)
fit2 <- lm(net.div~mean.pd.herb.resid)
aa2 <- summary(fit2)
regress.netdiv <- rbind(regress.netdiv, c("postKT20, no var", aa1$coefficients[2], aa1$r.squared, aa2$coefficients[2], aa2$r.squared))


load("postKT203xSHIFTS")

pd.transition.resid <- c()
pd.woody.resid <- c()
pd.herb.resid <- c()
for(index in 1:100){
    pd.transition.resid <- rbind(pd.transition.resid, 100*((transition.edge.lengths[index,] - true.transition.edge.lengths)/true.transition.edge.lengths))
    pd.woody.resid <- rbind(pd.woody.resid, 100*((woody.edge.lengths[index,] - true.woody.edge.lengths)/true.woody.edge.lengths))
    pd.herb.resid <- rbind(pd.herb.resid, 100*((herb.edge.lengths[index,] - true.herb.edge.lengths)/true.herb.edge.lengths))
}
mean.pd.transition.resid <- apply(pd.transition.resid, 1, mean)
mean.pd.woody.resid <- apply(pd.woody.resid, 1, mean)
mean.pd.herb.resid <- apply(pd.herb.resid, 1, mean)

sum.stats <- rbind(sum.stats, c("postKT20, 3x", mean(mean.pd.woody.resid), mean(mean.pd.herb.resid), mean(root.age), mean(net.div), mean(extinct.frac)))

fit1 <- lm(root.age~mean.pd.woody.resid)
aa1 <- summary(fit1)
fit2 <- lm(root.age~mean.pd.herb.resid)
aa2 <- summary(fit2)
regress.root.age <- rbind(regress.root.age , c("postKT20, 3x", aa1$coefficients[2], aa1$r.squared, aa2$coefficients[2], aa2$r.squared))

fit1 <- lm(extinct.frac~mean.pd.woody.resid)
aa1 <- summary(fit1)
fit2 <- lm(extinct.frac~mean.pd.herb.resid)
aa2 <- summary(fit2)
regress.ef <- rbind(regress.ef , c("postKT20, 3x", aa1$coefficients[2], aa1$r.squared, aa2$coefficients[2], aa2$r.squared))

fit1 <- lm(net.div~mean.pd.woody.resid)
aa1 <- summary(fit1)
fit2 <- lm(net.div~mean.pd.herb.resid)
aa2 <- summary(fit2)
regress.netdiv <- rbind(regress.netdiv, c("postKT20, 3x", aa1$coefficients[2], aa1$r.squared, aa2$coefficients[2], aa2$r.squared))

colnames(sum.stats) <- c("cal scheme", "woody %diff", "herb. %diff", "root age", "net.div", "ef")
colnames(regress.root.age) <- c("cal scheme", "woody slope", "woody R-sq", "herb slope", "herb R-sq")
colnames(regress.ef) <- c("cal scheme", "woody slope", "woody R-sq", "herb slope", "herb R-sq")
colnames(regress.netdiv) <- c("cal scheme", "woody slope", "woody R-sq", "herb slope", "herb R-sq")

sum.stats.df <- as.data.frame(sum.stats)
regress.root.age.df <- as.data.frame(regress.root.age)
regress.ef.df <- as.data.frame(regress.ef)
regress.netdiv.df <- as.data.frame(regress.netdiv)




######################################################################################################################################
######################################################################################################################################
### PLOT CODE 1: 3x Summary
######################################################################################################################################
######################################################################################################################################

#### Plot 1 <- Net Div., Extinct. Frac., Root Age top row -- bias by branch type second row.
library(viridis)

pdf("Fig3_basic_summary.pdf", height=7.5, width=7.5)
par(mfcol=c(3,3), mar=c(4,4,0.5,0.5), oma=c(1.5,2,1,1))

#Panel A -- Root Age
load("all3xSHIFTS")
h <- hist(root.age, plot=FALSE, breaks=seq(75,155,8))
h$counts = h$counts/sum(h$counts)
plot(h, axes=FALSE, xlab="", ylab="", ylim=c(0,.5), xlim=c(80,160), col="gray", main="")
par(tck=.01)
axis(2, at = seq(0,.5, by = .1), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(80,160, by = 20), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
text(80,.5, "All cals.", pos=4)
#title(xlab="Root age (MA)", line=2.5)
title(ylab="Frequency", line=2.5)
abline(v=true.root.age, lty=2)
mtext("(A)",side=3, line=0, adj=0)


load("postKT3xSHIFTS")
h <- hist(root.age, plot=FALSE, breaks=seq(75,155,8))
h$counts = h$counts/sum(h$counts)
plot(h, axes=FALSE, xlab="", ylab="", ylim=c(0,.5), xlim=c(80,160), col="gray", main="")
par(tck=.01)
axis(2, at = seq(0,.5, by = .1), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(80,160, by = 20), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
text(80,.5, "Post KP cals.", pos=4)
#title(xlab="Root age (MA)", line=2.5)
title(ylab="Frequency", line=2.5)
abline(v=true.root.age, lty=2)
mtext("(D)",side=3, line=0, adj=0)

load("postKT203xSHIFTS")
h <- hist(root.age, plot=FALSE, breaks=seq(75,155,8))
h$counts = h$counts/sum(h$counts)
plot(h, axes=FALSE, xlab="", ylab="", ylim=c(0,.5), xlim=c(80,160), main="", col="gray")
par(tck=.01)
axis(2, at = seq(0,.5, by = .1), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(80,160, by = 20), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
text(80,.5, "20 random post-KP cals", pos=4)
title(xlab="Root age (MA)", line=2)
title(ylab="Frequency", line=2.5)
abline(v=true.root.age, lty=2)
mtext("(G)",side=3, line=0, adj=0)



load("all3xSHIFTS")
h <- hist(net.div, plot=FALSE, breaks=seq(0,.25,.025))
h$counts = h$counts/sum(h$counts)
plot(h, axes=FALSE, xlab="", ylab="", ylim=c(0,.5), xlim=c(0,.250), main="", col="gray")
par(tck=.01)
axis(2, at = seq(0,.5, by = .1), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(0,.25, by = .025), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
#title(xlab="Root age (MA)", line=2.5)
#title(ylab="Frequency", line=2)
abline(v=true.net.div, lty=2)
mtext("(B)",side=3, line=0, adj=0)

load("postKT3xSHIFTS")
h <- hist(net.div, plot=FALSE, breaks=seq(0,.25,.025))
h$counts = h$counts/sum(h$counts)
plot(h, axes=FALSE, xlab="", ylab="", ylim=c(0,.5), xlim=c(0,.250), main="", col="gray")
par(tck=.01)
axis(2, at = seq(0,.5, by = .1), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(0,.25, by = .025), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
#title(xlab="Root age (MA)", line=2.5)
#title(ylab="Frequency", line=2)
abline(v=true.net.div, lty=2)
mtext("(E)",side=3, line=0, adj=0)

load("postKT203xSHIFTS")
h <- hist(net.div, plot=FALSE, breaks=seq(0,.25,.025))
h$counts = h$counts/sum(h$counts)
plot(h, axes=FALSE, xlab="", ylab="", ylim=c(0,.5), xlim=c(0,.250), main="", col="gray")
par(tck=.01)
axis(2, at = seq(0,.5, by = .1), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(0,.25, by = .025), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(xlab=expression(Net~diversification~rate~(r)), line=2)
#title(ylab="Frequency", line=2)
abline(v=true.net.div, lty=2)
mtext("(H)",side=3, line=0, adj=0)



load("all3xSHIFTS")
h <- hist(extinct.frac, plot=FALSE, breaks=seq(0,1,.1))
h$counts = h$counts/sum(h$counts)
plot(h, axes=FALSE, xlab="", ylab="", ylim=c(0,.5), xlim=c(0,1), main="", col="gray")
par(tck=.01)
axis(2, at = seq(0,.5, by = .1), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(0,1, by = .2), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
#title(xlab="Root age (MA)", line=2.5)
#title(ylab="Frequency", line=2)
abline(v=true.extinct.frac, lty=2)
mtext("(C)",side=3, line=0, adj=0)


load("postKT3xSHIFTS")
h <- hist(extinct.frac, plot=FALSE, breaks=seq(0,1,.1))
h$counts = h$counts/sum(h$counts)
plot(h, axes=FALSE, xlab="", ylab="", ylim=c(0,.5), xlim=c(0,1), main="", col="gray")
par(tck=.01)
axis(2, at = seq(0,.5, by = .1), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(0,1, by = .2), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
#title(xlab="Root age (MA)", line=2.5)
#title(ylab="Frequency", line=2)
abline(v=true.extinct.frac, lty=2)
mtext("(F)",side=3, line=0, adj=0)

load("postKT203xSHIFTS")
h <- hist(extinct.frac, plot=FALSE, breaks=seq(0,1,.1))
h$counts = h$counts/sum(h$counts)
plot(h, axes=FALSE, xlab="", ylab="", ylim=c(0,.5), xlim=c(0,1), main="", col="gray")
par(tck=.01)
axis(2, at = seq(0,.5, by = .1), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(0,1, by = .2), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(xlab=expression(Extinction~fraction~(epsilon)), line=2)
#title(ylab="Frequency", line=2)
abline(v=true.extinct.frac, lty=2)
mtext("(I)",side=3, line=0, adj=0)

dev.off()




######################################################################################################################################
######################################################################################################################################
### PLOT CODE 2: By branch
######################################################################################################################################
######################################################################################################################################

pdf("Fig4_basic_summary.pdf", height=7.5, width=7.5)
par(mfcol=c(3,3), mar=c(4,4,0.5,0.5), oma=c(1.5,2,1,1))

#Panel A -- Root Age
load("all3xSHIFTS")
pd.transition.resid <- c()
pd.woody.resid <- c()
pd.herb.resid <- c()
for(index in 1:100){
    pd.transition.resid <- rbind(pd.transition.resid, 100*((transition.edge.lengths[index,] - true.transition.edge.lengths)/true.transition.edge.lengths))
    pd.woody.resid <- rbind(pd.woody.resid, 100*((woody.edge.lengths[index,] - true.woody.edge.lengths)/true.woody.edge.lengths))
    pd.herb.resid <- rbind(pd.herb.resid, 100*((herb.edge.lengths[index,] - true.herb.edge.lengths)/true.herb.edge.lengths))
}
mean.pd.transition.resid <- apply(pd.transition.resid, 2, mean)
mean.pd.woody.resid <- apply(pd.woody.resid, 2, mean)
mean.pd.herb.resid <- apply(pd.herb.resid, 2, mean)
plot(true.transition.edge.lengths, mean.pd.transition.resid, axes=FALSE, xlab="", ylab="", ylim=c(-500,500), xlim=c(0,100), main="", col=viridis(11)[4], pch=19)
par(tck=.01)
axis(2, at = seq(-500,500, by = 250), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(0,100, by = 20), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
#title(xlab="True edge length", line=2)
title(ylab="Percent difference from true length", line=2.5)
text(10,500, "All cals.", pos=4)
abline(h=0, lty=2)
mtext("(A)",side=3, line=0, adj=0)

load("postKT3xSHIFTS")
pd.transition.resid <- c()
pd.woody.resid <- c()
pd.herb.resid <- c()
for(index in 1:100){
    pd.transition.resid <- rbind(pd.transition.resid, 100*((transition.edge.lengths[index,] - true.transition.edge.lengths)/true.transition.edge.lengths))
    pd.woody.resid <- rbind(pd.woody.resid, 100*((woody.edge.lengths[index,] - true.woody.edge.lengths)/true.woody.edge.lengths))
    pd.herb.resid <- rbind(pd.herb.resid, 100*((herb.edge.lengths[index,] - true.herb.edge.lengths)/true.herb.edge.lengths))
}
mean.pd.transition.resid <- apply(pd.transition.resid, 2, mean)
mean.pd.woody.resid <- apply(pd.woody.resid, 2, mean)
mean.pd.herb.resid <- apply(pd.herb.resid, 2, mean)

plot(true.transition.edge.lengths, mean.pd.transition.resid, axes=FALSE, xlab="", ylab="", ylim=c(-500,500), xlim=c(0,100), main="", col=viridis(11)[4], pch=19)
par(tck=.01)
axis(2, at = seq(-500,500, by = 250), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(0,100, by = 20), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
#title(xlab="True edge length", line=2)
title(ylab="Percent difference from true length", line=2.5)
text(10,500, "Post KP cals.", pos=4)
#text(30,500, "Transition edges", cex=0.75)
abline(h=0, lty=2)
mtext("(D)",side=3, line=0, adj=0)

load("postKT203xSHIFTS")
pd.transition.resid <- c()
pd.woody.resid <- c()
pd.herb.resid <- c()
for(index in 1:100){
    pd.transition.resid <- rbind(pd.transition.resid, 100*((transition.edge.lengths[index,] - true.transition.edge.lengths)/true.transition.edge.lengths))
    pd.woody.resid <- rbind(pd.woody.resid, 100*((woody.edge.lengths[index,] - true.woody.edge.lengths)/true.woody.edge.lengths))
    pd.herb.resid <- rbind(pd.herb.resid, 100*((herb.edge.lengths[index,] - true.herb.edge.lengths)/true.herb.edge.lengths))
}
mean.pd.transition.resid <- apply(pd.transition.resid, 2, mean)
mean.pd.woody.resid <- apply(pd.woody.resid, 2, mean)
mean.pd.herb.resid <- apply(pd.herb.resid, 2, mean)
plot(true.transition.edge.lengths, mean.pd.transition.resid, axes=FALSE, xlab="", ylab="", ylim=c(-500,500), xlim=c(0,100), main="", col=viridis(11)[4], pch=19)
par(tck=.01)
axis(2, at = seq(-500,500, by = 250), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(0,100, by = 20), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(xlab="True edge length of transition edges", line=2)
title(ylab="Percent difference from true length", line=2.5)
text(10,500, "20 random post-KP cals", pos=4)
#text(30,500, "Transition edges", cex=0.75)
abline(h=0, lty=2)
mtext("(G)",side=3, line=0, adj=0)


load("all3xSHIFTS")
pd.transition.resid <- c()
pd.woody.resid <- c()
pd.herb.resid <- c()
for(index in 1:100){
    pd.transition.resid <- rbind(pd.transition.resid, 100*((transition.edge.lengths[index,] - true.transition.edge.lengths)/true.transition.edge.lengths))
    pd.woody.resid <- rbind(pd.woody.resid, 100*((woody.edge.lengths[index,] - true.woody.edge.lengths)/true.woody.edge.lengths))
    pd.herb.resid <- rbind(pd.herb.resid, 100*((herb.edge.lengths[index,] - true.herb.edge.lengths)/true.herb.edge.lengths))
}
mean.pd.transition.resid <- apply(pd.transition.resid, 2, mean)
mean.pd.woody.resid <- apply(pd.woody.resid, 2, mean)
mean.pd.herb.resid <- apply(pd.herb.resid, 2, mean)
plot(true.woody.edge.lengths, mean.pd.woody.resid, axes=FALSE, xlab="", ylab="", ylim=c(-500,500), xlim=c(0,100), main="", col=viridis(11)[1], pch=19)
ll <- loess.smooth(true.woody.edge.lengths, mean.pd.woody.resid)
lines(ll, col="red", lty=3, lwd=2)
par(tck=.01)
axis(2, at = seq(-500,500, by = 250), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(0,100, by = 20), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
#title(xlab="True edge length", line=2.5)
#title(ylab="Percent difference", line=2)
abline(h=0, lty=2)
mtext("(B)",side=3, line=0, adj=0)

load("postKT3xSHIFTS")
pd.transition.resid <- c()
pd.woody.resid <- c()
pd.herb.resid <- c()
for(index in 1:100){
    pd.transition.resid <- rbind(pd.transition.resid, 100*((transition.edge.lengths[index,] - true.transition.edge.lengths)/true.transition.edge.lengths))
    pd.woody.resid <- rbind(pd.woody.resid, 100*((woody.edge.lengths[index,] - true.woody.edge.lengths)/true.woody.edge.lengths))
    pd.herb.resid <- rbind(pd.herb.resid, 100*((herb.edge.lengths[index,] - true.herb.edge.lengths)/true.herb.edge.lengths))
}
mean.pd.transition.resid <- apply(pd.transition.resid, 2, mean)
mean.pd.woody.resid <- apply(pd.woody.resid, 2, mean)
mean.pd.herb.resid <- apply(pd.herb.resid, 2, mean)
plot(true.woody.edge.lengths, mean.pd.woody.resid, axes=FALSE, xlab="", ylab="", ylim=c(-500,500), xlim=c(0,100), main="", col=viridis(11)[1], pch=19)
ll <- loess.smooth(true.woody.edge.lengths, mean.pd.woody.resid)
lines(ll, col="red", lty=3, lwd=2)
par(tck=.01)
axis(2, at = seq(-500,500, by = 250), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(0,100, by = 20), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
#title(xlab="True edge length", line=2.5)
#title(ylab="Percent difference", line=2)
#text(30,500, "'woody' edges", cex=0.75)
abline(h=0, lty=2)
mtext("(E)",side=3, line=0, adj=0)

load("postKT203xSHIFTS")
pd.transition.resid <- c()
pd.woody.resid <- c()
pd.herb.resid <- c()
for(index in 1:100){
    pd.transition.resid <- rbind(pd.transition.resid, 100*((transition.edge.lengths[index,] - true.transition.edge.lengths)/true.transition.edge.lengths))
    pd.woody.resid <- rbind(pd.woody.resid, 100*((woody.edge.lengths[index,] - true.woody.edge.lengths)/true.woody.edge.lengths))
    pd.herb.resid <- rbind(pd.herb.resid, 100*((herb.edge.lengths[index,] - true.herb.edge.lengths)/true.herb.edge.lengths))
}
mean.pd.transition.resid <- apply(pd.transition.resid, 2, mean)
mean.pd.woody.resid <- apply(pd.woody.resid, 2, mean)
mean.pd.herb.resid <- apply(pd.herb.resid, 2, mean)
plot(true.woody.edge.lengths, mean.pd.woody.resid, axes=FALSE, xlab="", ylab="", ylim=c(-500,500), xlim=c(0,100), main="", col=viridis(11)[1], pch=19)
ll <- loess.smooth(true.woody.edge.lengths, mean.pd.woody.resid)
lines(ll, col="red", lty=3, lwd=2)
par(tck=.01)
axis(2, at = seq(-500,500, by = 250), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(0,100, by = 20), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(xlab="True edge length of 'woody' edges", line=2)
#title(ylab="Percent difference", line=2)
#text(30,500, "'woody' edges", cex=0.75)
abline(h=0, lty=2)
mtext("(H)",side=3, line=0, adj=0)


load("all3xSHIFTS")
pd.transition.resid <- c()
pd.woody.resid <- c()
pd.herb.resid <- c()
for(index in 1:100){
    pd.transition.resid <- rbind(pd.transition.resid, 100*((transition.edge.lengths[index,] - true.transition.edge.lengths)/true.transition.edge.lengths))
    pd.woody.resid <- rbind(pd.woody.resid, 100*((woody.edge.lengths[index,] - true.woody.edge.lengths)/true.woody.edge.lengths))
    pd.herb.resid <- rbind(pd.herb.resid, 100*((herb.edge.lengths[index,] - true.herb.edge.lengths)/true.herb.edge.lengths))
}
mean.pd.transition.resid <- apply(pd.transition.resid, 2, mean)
mean.pd.woody.resid <- apply(pd.woody.resid, 2, mean)
mean.pd.herb.resid <- apply(pd.herb.resid, 2, mean)

plot(true.herb.edge.lengths, mean.pd.herb.resid, axes=FALSE, xlab="", ylab="", ylim=c(-500,500), xlim=c(0,100), main="", col=viridis(11)[7], pch=19)
ll <- loess.smooth(true.herb.edge.lengths, mean.pd.herb.resid)
lines(ll, col="red", lty=3, lwd=2)
par(tck=.01)
axis(2, at = seq(-500,500, by = 250), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(0,100, by = 20), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
#title(xlab="True edge length", line=2)
#title(ylab="Percent difference", line=2)
abline(h=0, lty=2)
mtext("(C)",side=3, line=0, adj=0)

load("postKT3xSHIFTS")
pd.transition.resid <- c()
pd.woody.resid <- c()
pd.herb.resid <- c()
for(index in 1:100){
    pd.transition.resid <- rbind(pd.transition.resid, 100*((transition.edge.lengths[index,] - true.transition.edge.lengths)/true.transition.edge.lengths))
    pd.woody.resid <- rbind(pd.woody.resid, 100*((woody.edge.lengths[index,] - true.woody.edge.lengths)/true.woody.edge.lengths))
    pd.herb.resid <- rbind(pd.herb.resid, 100*((herb.edge.lengths[index,] - true.herb.edge.lengths)/true.herb.edge.lengths))
}
mean.pd.transition.resid <- apply(pd.transition.resid, 2, mean)
mean.pd.woody.resid <- apply(pd.woody.resid, 2, mean)
mean.pd.herb.resid <- apply(pd.herb.resid, 2, mean)
plot(true.herb.edge.lengths, mean.pd.herb.resid, axes=FALSE, xlab="", ylab="", ylim=c(-500,500), xlim=c(0,100), main="", col=viridis(11)[7], pch=19)
ll <- loess.smooth(true.herb.edge.lengths, mean.pd.herb.resid)
lines(ll, col="red", lty=3, lwd=2)
par(tck=.01)
axis(2, at = seq(-500,500, by = 250), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(0,100, by = 20), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
#title(xlab="True edge length", line=2)
#title(ylab="Percent difference", line=2)
#text(30,500, "'herb' edges", cex=0.75)
abline(h=0, lty=2)
mtext("(F)",side=3, line=0, adj=0)

load("postKT203xSHIFTS")
pd.transition.resid <- c()
pd.woody.resid <- c()
pd.herb.resid <- c()
for(index in 1:100){
    pd.transition.resid <- rbind(pd.transition.resid, 100*((transition.edge.lengths[index,] - true.transition.edge.lengths)/true.transition.edge.lengths))
    pd.woody.resid <- rbind(pd.woody.resid, 100*((woody.edge.lengths[index,] - true.woody.edge.lengths)/true.woody.edge.lengths))
    pd.herb.resid <- rbind(pd.herb.resid, 100*((herb.edge.lengths[index,] - true.herb.edge.lengths)/true.herb.edge.lengths))
}
mean.pd.transition.resid <- apply(pd.transition.resid, 2, mean)
mean.pd.woody.resid <- apply(pd.woody.resid, 2, mean)
mean.pd.herb.resid <- apply(pd.herb.resid, 2, mean)
plot(true.herb.edge.lengths, mean.pd.herb.resid, axes=FALSE, xlab="", ylab="", ylim=c(-500,500), xlim=c(0,100), main="", col=viridis(11)[7], pch=19)
ll <- loess.smooth(true.herb.edge.lengths, mean.pd.herb.resid)
lines(ll, col="red", lty=3, lwd=2)
par(tck=.01)
axis(2, at = seq(-500,500, by = 250), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(0,100, by = 20), las =1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
title(xlab="True edge length of 'herb' edges", line=2)
#title(ylab="Percent difference", line=2)
#text(30,500, "'herb' edges", cex=0.75)
abline(h=0, lty=2)
mtext("(I)",side=3, line=0, adj=0)

dev.off()

