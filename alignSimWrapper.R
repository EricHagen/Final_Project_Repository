
setwd("STEP_3_SIMULATE_DATA/")
source("nucGapSim.R")
source("cleanSeqs.R")

#A->C; A->G; A->T; C->G; C->T; G->T

MakeChunky <- function(seq.mat, min.block.size=25){
  for(i in 1:(dim(seq.mat)[1])){
    count = 0
    for(j in 1:(dim(seq.mat)[2])){
      if(seq.mat[i,j] == TRUE){
        count = count + 1
      }else{
        if(count < min.block.size){
          seq.mat[i,(j-count):j] = FALSE
        }
        count = 0
      }
    }
  }
  seq.mat <- as.data.frame(seq.mat)
  return(seq.mat)
}

alignSimWrap <- function(start.rep, end.rep, rate.phy.set, nsites){
    #This is for providing real missing data chunks to the data:
    major.seq <- read.alignment("matk.matched.pruned.nametrunc.legit.stop.removed.cleaned.brunFix.fasta", format="fasta")
    major.seq.mat <- as.matrix.alignment(major.seq)
    missing.set <- major.seq.mat == "-"

    #Now let us smooth this out so that it is missing chunks not just the simple single gaps or smallish gaps that may represent real indels:
    final.missing.set <- MakeChunky(missing.set)
    final.missing.set <- final.missing.set[,1:1515]
    ##
    
    system("mkdir FINISHED_ALIGNMENTS/")
    
    #Alright. Let us go to work...
    spscores.mat <- matrix(0, 10, 4)
    for(i in start.rep:end.rep){
        print(i)
        final.missing.set <- final.missing.set[rate.phy.set[[i]]$tip.label,]
        gtr.pars <- c(1.35, 1.61, 0.254, 0.979, 1.56, 1.0)
        base.freqs <- c(0.297, 0.174, 0.165, 0.363)
        shape.par <- 1.08
        indel.rate <- mean(gtr.pars) * 0.1
        
        #Part 1: simulate sequence evolution according known parameters but also allow for gaps.
        sim.fasta <- NucIndelSimulator(rate.phy.set[[i]], gtr.pars=gtr.pars, shape.par=shape.par, ncats=4, base.freqs=base.freqs, nsites=nsites, nuc.model="GTR", insert.rate=indel.rate, delete.rate=indel.rate)
        write.dna(sim.fasta, file=paste("simRep", i, ".fasta", sep=""), colw=1000000, format="fasta")
        
        sim.fasta.chunks <- sim.fasta
        sim.fasta.chunks[final.missing.set==TRUE] = "-"
        write.dna(sim.fasta.chunks, file=paste("simRepChunked", i, ".fasta", sep=""), colw=1000000, format="fasta")
        
        #Part 2: Align using MAFFT:
        system(paste("mafft --auto", paste("simRep", i, ".fasta", sep=""), ">", paste("simRep", i, ".mafft.aln", sep="")))
        system(paste("mafft --auto", paste("simRepChunked", i, ".fasta", sep=""), ">", paste("simRepChunked", i, ".mafft.aln", sep="")))
        ##############
        
        system(paste("mv", paste("simRep", i, ".mafft.aln", sep=""), "FINISHED_ALIGNMENTS/", sep=" "))
        system(paste("mv", paste("simRepChunked", i, ".mafft.aln", sep=""), "FINISHED_ALIGNMENTS/", sep=" "))

        #Part 3: Clean and format for RAxML:
        seq <- read.alignment(paste("FINISHED_ALIGNMENTS/simRep", i, ".mafft.aln", sep=""), format="fasta")
        seq.clean <- CleanSeqs(seq, cutoff=.25)
        write.dna(seq.clean, file=paste("FINISHED_ALIGNMENTS/simRep", i, ".mafft.clean.aln", sep=""), format="fasta", colw=10000000)
        write.dna(seq.clean, file=paste("FINISHED_ALIGNMENTS/simRep", i, ".mafft.nex", sep=""), format="sequential", colw=10000000)
        
        #Part 4: Run RAxML on MAFFT only alignment:
        system(paste("../standard-RAxML-master/raxmlHPC-PTHREADS-SSE3 -T 24 -f a -p 12345 -s ", paste("FINISHED_ALIGNMENTS/simRep", i, ".mafft.nex", sep=""), "-x 12345 -#10 -m GTRCAT -n FINISHED_ALIGNMENTS/simMAFFT.", i, sep=""))
        ##############
        
        #Part 5: Clean and format for RAxML:
        seq <- read.alignment(paste("FINISHED_ALIGNMENTS/simRepChunked", i, ".mafft.aln", sep=""), format="fasta")
        seq.clean <- CleanSeqs(seq, cutoff=.25)
        write.dna(seq.clean, file=paste("FINISHED_ALIGNMENTS/simRepChunked", i, ".mafft.clean.aln", sep=""), format="fasta", colw=10000000)
        write.dna(seq.clean, file=paste("FINISHED_ALIGNMENTS/simRepChunked", i, ".mafft.nex", sep=""), format="sequential", colw=10000000)
        
        #Part 6: Run RAxML on MAFFT only alignment:
        system(paste("../standard-RAxML-master/raxmlHPC-PTHREADS-SSE3 -T 24 -f a -p 12345 -s ", paste("FINISHED_ALIGNMENTS/simRepChunked", i, ".mafft.nex", sep=""), "-x 12345 -#10 -m GTRCAT -n FINISHED_ALIGNMENTS/simMAFFT.chunked.", i, sep=""))
        ##############
    }
}


BlockBreaker <- function(list, fasta, filename){
    ##there is a lot of converting between lists and matrix formats to make this thing work
    as.matrix(fasta)->list.of.death
    n<-length(list.of.death)
    quip<-matrix(nrow=length(list.of.death[,1]), ncol=length(list.of.death[1,]))
    ##need to this to track names and omit NA values later
    namer<-matrix(nrow=length(list.of.death[,1]), ncol=length(list.of.death[1,]))
    names(fasta)->masterlist
    as.matrix(masterlist)->uber.list
    #the loop
    for(i in 1:n){
        #first get the name of the taxon at row i from the fasta file to see if it will be retained
        uber.list[i,]->love.me
        ##if this name is in the taxon list, get the name and sequence data and store into two separate matrices
        if(love.me%in%list==TRUE)
        {
            love.me->namer[i,]
            list.of.death[i,]->quip[i,]
            as.matrix(quip)->quip
        }
        
    }
    ###get rid of all NA values
    quip[rowSums(is.na(quip))<ncol(quip),]->quip3
    na.omit(namer)->namer2
    ##write fasta file, need to make changing outputs automatic
    write.fasta(quip3, names=namer2, nbchar = 60, file.out=filename, open = "w")
    #test for missing names with list%in%namer2
    #return(namer2)
}


#setwd("~/Beaulieu-lab-repo/matk_Sim_Files/STEP_2_RATING")
#setwd("~/Beaulieu-lab-repo/matk_Sim_Files/STEP_3_SIMULATE_DATA")
rate.phy.set <- read.tree("camp.matk.rated.start0025.rate003.3xshift.trees")
alignSimWrap(start.rep=1, end.rep=5, rate.phy.set=rate.phy.set, nsites=1515)




