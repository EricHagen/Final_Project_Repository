#RateEvolve for simulations in divergence time analysis

#written by Jeremy M. Beaulieu

library(ape)
library(phytools)

rateEvolve<-function(phy, data, type="ucln", nsims){
	for(j in 1:nsims){
		temp<-rateSim(phy,data, type=type)
		write.tree(temp,file="camp.matk.rated.start0025.rate003.3xshift.trees",append=TRUE)
	}
}


rateSim<-function(phy, data, type=c("autocor", "ucln"), rate=0.03, starting.val=0.0025, discrete.shift.scalar=3){
	
	data<-data.frame(data[,2], data[,2], row.names=data[,1])
	data<-data[phy$tip.label,]
	
	#Edge matrix formatting taken from OUwie
	int.states<-factor(phy$node.label)
	phy$node.label=as.numeric(int.states)
	tip.states<-factor(data[,1])
	data[,1]<-as.numeric(tip.states)
	k<-length(levels(int.states))
	n=max(phy$edge[,1])
	ntips=length(phy$tip.label)
		
	#Obtain root state and internal node labels
	root.state<-phy$node.label[1]
	int.state<-phy$node.label[-1]
	#New tree matrix to be used for subsetting regimes
	edges=cbind(c(1:(n-1)),phy$edge, nodeHeights(phy))
	edges=edges[sort.list(edges[,3]),]
	
	mm<-c(data[,1],int.state)

	regime <- matrix(0,nrow=length(mm),ncol=length(unique(mm)))
	#Generates an indicator matrix from the regime vector
	for (i in 1:length(mm)) {
		regime[i,mm[i]] <- 1 
	}
    if(type=="autocor"){
        #Finishes the edges matrix
        edges=cbind(edges,regime)
        #Resort the edge matrix so that it looks like the original matrix order
        edges=edges[sort.list(edges[,1]),]
        
        oldregime=root.state
        nodevar=rep(0,max(edges[,3]))
        
        for(i in 1:length(edges[,1])){
            anc = edges[i,2]
            oldtime=edges[i,4]
            newtime=edges[i,5]
            if(anc%in%edges[,3]){
                start = which(edges[,3]==anc)
                oldregime = which(edges[start,6:(k+5)]==1)
                prev.rate = nodevar[start]
            }else{
                #For the root:
                oldregime = oldregime
                prev.rate = starting.val
            }
            newregime=which(edges[i,6:(k+5)]==1)
            if(oldregime == newregime){
                newtime = newtime-oldtime
                if(length(nodevar[i-1])==0){
                    if(newregime == 1){
                        nodevar[i] = exp(rnorm(1, mean=log(prev.rate), sd=sqrt(newtime) * rate))
                    }else{
                        nodevar[i] = exp(rnorm(1, mean=log(prev.rate), sd=sqrt(newtime) * (rate * discrete.shift.scalar)))
                    }
                }else{
                    if(newregime == 1){
                        nodevar[i] = exp(rnorm(1, mean=log(prev.rate), sd=sqrt(newtime) * rate))
                    }else{
                        nodevar[i] = exp(rnorm(1, mean=log(prev.rate), sd=sqrt(newtime) * (rate * discrete.shift.scalar)))
                    }
                }
            }else{
                newtime = (newtime-oldtime)/2
                if(oldregime == 1){
                    epoch1 = exp(rnorm(1, mean=log(prev.rate), sd=sqrt(newtime) * rate))
                }else{
                    epoch1 = exp(rnorm(1, mean=log(prev.rate), sd=sqrt(newtime) * (rate * discrete.shift.scalar)))
                }
                if(newregime==1){
                    epoch2 = exp(rnorm(1, mean=log(epoch1), sd=sqrt(newtime) * rate))
                }
                if(newregime==2){
                    epoch2 = exp(rnorm(1, mean=log(epoch1), sd=sqrt(newtime) * (rate * discrete.shift.scalar)))
                }
                nodevar[i] <- epoch2
            }
            oldregime=newregime
        }
        phy$edge.length <- nodevar[-length(nodevar)] * phy$edge.length
    }
    
    if(type=="ucln"){
        #Finishes the edges matrix
        edges=cbind(edges,regime)
        #Resort the edge matrix so that it looks like the original matrix order
        edges=edges[sort.list(edges[,1]),]
        
        #3x differences based on Smith and Donoghue 2008:
        rate1<-rlnorm(10000, log(5e-4), .75)
        rate2<-rlnorm(10000, log(0.002), .75)
        rate3<-rlnorm(10000, log(0.006), .75)
        
        #	rate2<-rate1
        rates<-cbind(rate1,rate2,rate3)
        oldregime=root.state
        nodevar=rep(0,max(edges[,3]))
        
        for(i in 1:length(edges[,1])){
            anc = edges[i,2]
            oldtime=edges[i,4]
            newtime=edges[i,5]
            if(anc%in%edges[,3]){
                start=which(edges[,3]==anc)
                oldregime=which(edges[start,6:(k+5)]==1)
            }
            else{
                #For the root:
                oldregime=oldregime
            }
            newregime=which(edges[i,6:(k+5)]==1)
            if(oldregime==newregime){
                newtime=newtime-oldtime
                nodevar[i]=newtime*sample(rates[,oldregime],1)
            }
            else{
                newtime=(newtime-oldtime)/2
                epoch1=newtime*sample(rates[,oldregime],1)
                epoch2=newtime*sample(rates[,newregime],1)
                nodevar[i]<-epoch1+epoch2 
            }
            oldregime=newregime
        }
        phy$edge.length <- nodevar
    }
	phy
}


#tree<-read.tree("camp.matk.habit.recon.asterwood.carrotwood.tre")
#trait<-read.delim("camp.matk.habit.trait.asterwood.carrotwood.txt")
#rateEvolve(tree,trait, type="autocor", nsims=100)


