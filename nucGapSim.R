
######################################################################################################################################
######################################################################################################################################
### Simulator for alignments with indels
######################################################################################################################################
######################################################################################################################################
library(ape)
library(seqinr)
library(expm)
#written by Jeremy M. Beaulieu

SingleSiteUpPass <- function(phy, Q, root.value){
	# Randomly choose starting state at root using the root.values as the probability:
	root.value = sample.int(dim(Q)[2], 1, FALSE, prob=root.value)
	# Reorder the phy:
    ntips <- length(phy$tip.label)
    N <- dim(phy$edge)[1]
    ROOT <- ntips + 1 #perhaps use an accessor to get the root node id
	# Generate vector that contains the simulated states:
    sim.data.site <- integer(ntips + phy$Nnode)
    sim.data.site[ROOT] <- as.integer(root.value)
    anc <- phy$edge[, 1]
    des <- phy$edge[, 2]
    edge.length <- phy$edge.length
	for (i in N:1) {
		p <- expm(Q * edge.length[i], method="Ward77")[sim.data.site[anc[i]], ]
		sim.data.site[des[i]] <- sample.int(dim(Q)[2], size = 1, FALSE, prob = p)
	}
	sim.data.site <- sim.data.site[1:ntips]
	return(sim.data.site)
}


NucIndelSimulator <- function(phy, gtr.pars, shape.par, ncats, base.freqs, nsites, nuc.model, insert.rate, delete.rate, block.size=3){
    if(nuc.model == "JC") {
        base.freqs=base.freqs
        nuc.mutation.rates <- CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
    }
    if(nuc.model == "GTR") {
        base.freqs=base.freqs
        nuc.mutation.rates <- CreateNucleotideMutationMatrix(gtr.pars, model=nuc.model, base.freqs=base.freqs)
    }
    if(nuc.model == "UNREST") {
        nuc.mutation.rates <- CreateNucleotideMutationMatrix(gtr.pars, model=nuc.model)
    }

    Q_mat <- nuc.mutation.rates
    diag(Q_mat) = 0
    diag(Q_mat) <- -rowSums(Q_mat)
    scale.factor <- -sum(diag(Q_mat) * base.freqs, na.rm=TRUE)
    Q_mat_scaled = Q_mat * (1/scale.factor)

    # Perform simulation by looping over desired number of sites. The optimal aa for any given site is based on the user input vector of optimal AA:
    sim.nuc.data <- matrix(0, Ntip(phy), nsites)
    if(block.size==1){
        indel.data <- matrix(0, Ntip(phy), nsites)
    }
    if(block.size==3){
        indel.data <- matrix(0, Ntip(phy), nsites/3)
    }
    phy <- reorder(phy, "postorder")
    for(site in 1:nsites){
        Q_tmp <- Q_mat_scaled * rgamma(1, 1.08, 1.08)
        sim.nuc.data[,site] = SingleSiteUpPass(phy=phy, Q=Q_tmp, root.value=base.freqs)
        if(block.size==1){
            indel.data[,site] <- IndelSimulator(phy=phy, insert.rate=insert.rate, delete.rate=delete.rate)
        }
        if(block.size==3){
            if((site%%3)==0){
                indel.data[,site/3] <- IndelSimulator(phy=phy, insert.rate=insert.rate, delete.rate=delete.rate)
            }
        }
    }
    nuc.names <- n2s(0:3)
    # Finally, translate this information into a matrix of nucleotides -- this format allows for write.dna() to write a fasta formatted file:
    nucleotide.data <- matrix(NA, Ntip(phy), nsites)
    if(block.size == 1){
        for(nuc.sequence in seq(1, nsites, by=1)){
            nucleotide.data[which(indel.data[,nuc.sequence]==1),nuc.sequence] = nuc.names[sim.nuc.data[which(indel.data[,nuc.sequence]==1),nuc.sequence]]
            nucleotide.data[which(indel.data[,nuc.sequence]==2),nuc.sequence] = nuc.names[sim.nuc.data[which(indel.data[,nuc.sequence]==2),nuc.sequence]]
        }
    }
    if(block.size == 3){
        for(nuc.sequence in seq(1, nsites, by=1)){
            if((nuc.sequence%%3)==0){
                nucleotide.data[which(indel.data[,nuc.sequence/3]==1),(nuc.sequence-2):nuc.sequence] = nuc.names[sim.nuc.data[which(indel.data[,nuc.sequence/3]==1),(nuc.sequence-2):nuc.sequence]]
                nucleotide.data[which(indel.data[,nuc.sequence/3]==2),(nuc.sequence-2):nuc.sequence] = nuc.names[sim.nuc.data[which(indel.data[,nuc.sequence/3]==2),(nuc.sequence-2):nuc.sequence]]
            }
        }
    }
    nucleotide.data[is.na(nucleotide.data)] = "-"
    rownames(nucleotide.data) <- phy$tip.label
    # Done.
    return(nucleotide.data)
}


IndelSimulator <- function(phy, insert.rate, delete.rate){
    root.value <- sample.int(3, 1, FALSE, prob=c(1/2, 1/2, 0))
    root.p <- rep(0, 3)
    root.p[root.value] = 1
    sim.indel.data <- matrix(0, nrow=Ntip(phy), 1)
    Q_indelmat = matrix(c(0,0,0,insert.rate,0,0,0,delete.rate,0), 3, 3)
    diag(Q_indelmat) = -rowSums(Q_indelmat)
    scale.factor <- -sum(diag(Q_indelmat))
    Q_indelmat <- Q_indelmat * (1/scale.factor)
    sim.indel.data[,1] = SingleSiteUpPass(phy, Q=Q_indelmat, root.value=root.p)
    rownames(sim.indel.data) <- phy$tip.label
	# Done.
	return(sim.indel.data)
}


DiscreteGamma <- function (shape, ncats){
    quantiles <- qgamma((1:(ncats - 1))/ncats, shape = shape, rate = shape)
    return(diff(c(0, pgamma(quantiles * shape, shape + 1), 1)) * ncats)
}


CreateNucleotideMutationMatrix <- function(rates, model="JC", base.freqs=NULL) {
    if(model == "JC") {
        nuc.mutation.rates <- matrix(data=rates[1], nrow=4, ncol=4)
        rownames(nuc.mutation.rates) <- n2s(0:3)
        colnames(nuc.mutation.rates) <- n2s(0:3)
        diag(nuc.mutation.rates) <- 0
        diag(nuc.mutation.rates) <- -rowSums(nuc.mutation.rates)
        if(!is.null(base.freqs)){
            diag(nuc.mutation.rates) = 0
            nuc.mutation.rates = t(nuc.mutation.rates * base.freqs)
            diag(nuc.mutation.rates) = -rowSums(nuc.mutation.rates)
        }
        return(nuc.mutation.rates)
    }
    if(model == "GTR") {
        index <- matrix(NA, 4, 4)
        np <- 5
        sel <- col(index) < row(index)
        sel[4,3] = FALSE
        index[sel] <- 1:np
        index <- t(index)
        index[sel] <- 1:np
        nuc.mutation.rates <- matrix(0, nrow=4, ncol=4)
        nuc.mutation.rates<-matrix(rates[index], dim(index))
        rownames(nuc.mutation.rates) <- n2s(0:3)
        colnames(nuc.mutation.rates) <- n2s(0:3)
        nuc.mutation.rates[4,3] = nuc.mutation.rates[3,4] = 1
        diag(nuc.mutation.rates) <- 0
        diag(nuc.mutation.rates) <- -rowSums(nuc.mutation.rates)
        if(!is.null(base.freqs)){
            diag(nuc.mutation.rates) = 0
            nuc.mutation.rates = t(nuc.mutation.rates * base.freqs)
            diag(nuc.mutation.rates) = -rowSums(nuc.mutation.rates)
        }
        return(nuc.mutation.rates)
    }
    if(model == "UNREST") {
        index <- matrix(NA, 4, 4)
        np <- 11
        index[col(index) != row(index)] <- 1:np
        nuc.mutation.rates <- matrix(0, nrow=4, ncol=4)
        nuc.mutation.rates<-matrix(rates[index], dim(index))
        rownames(nuc.mutation.rates) <- n2s(0:3)
        colnames(nuc.mutation.rates) <- n2s(0:3)
        nuc.mutation.rates[4,3] = 1
        diag(nuc.mutation.rates) <- 0
        diag(nuc.mutation.rates) <- -rowSums(nuc.mutation.rates)
        return(nuc.mutation.rates)
    }
}



