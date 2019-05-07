##treePL batch script

#written by Jeremy M. Beaulieu
library(geiger)
library(parallel)

GetFullHashCalibrations <- function(reference.tree, target.tree, tax.table){
    #Do congruify on tree i:
	congruify.output.full = congruify.phylo(reference=reference.tree, target=target.tree, taxonomy=tax.table)
    congruify.calibrations = congruify.output.full$calibrations
    n.nodes <- Nnode(congruify.output.full$target)
    n.tips <- Ntip(congruify.output.full$target)
    node.labels <- c(rep(1,n.tips), congruify.output.full$target$node.label)
    node.state <- c()
    for(node.index in 1:dim(congruify.calibrations)[1]){
        tmp <- getMRCA(congruify.output.full$target, tip=c(congruify.calibrations$taxonA[node.index],congruify.calibrations$taxonB[node.index]))
        node.state <- c(node.state, node.labels[tmp])
    }
    names_mat <- names(congruify.calibrations)
    congruify.calibrations <- cbind(congruify.calibrations, node.state)
    names(congruify.calibrations) <- c(names_mat, "habit")
    return(congruify.calibrations)
}


ChooseCalibrationSet <- function(congruify.output, num.cals=NULL, type=c("any", "woody.only", "herb.only")){
    if(type == "any"){
        if(is.null(num.cals)){
            congruify.output.red <- congruify.output
        }else{
            rows.to.choose <- sample(1:dim(congruify.output)[1], num.cals)
            congruify.output.red <- congruify.output[rows.to.choose,]
        }
    }
    if(type == "woody.only"){
        if(is.null(num.cals)){
            congruify.output.red <- congruify.output[which(congruify.output$habit==1),]
        }else{
            congruify.output <- congruify.output[which(congruify.output$habit==1),]
            rows.to.choose <- sample(1:dim(congruify.output)[1], num.cals)
            congruify.output.red <- congruify.output[rows.to.choose,]
        }
    }
    if(type == "herb.only"){
        if(is.null(num.cals)){
            congruify.output.red <- congruify.output[which(congruify.output$habit==2),]
        }else{
            congruify.output <- congruify.output[which(congruify.output$habit==2),]
            rows.to.choose <- sample(1:dim(congruify.output)[1], num.cals)
            congruify.output.red <- congruify.output[rows.to.choose,]
        }
    }
    return(congruify.output.red)
}


DoPrimeAnalysis <- function(calibration.set, target.tree, index){
    start.index <- end.index <- index
    ###The overall loop is looping over trees you are referencing:
    for (tree.index in start.index:end.index) {
       write.tree(target.tree[[tree.index]], paste("camp.rated.", index, ".tre", sep=""))
        cat("treefile =", paste("camp.rated.", index, ".tre", sep=""),"\nnumsites = 1515\n", file=paste("tree", tree.index, ".cppr8s", sep=""),sep="")
        
        ####### CALIBRATION BIT IS ADDED HERE #######
        for(mrca.index in sequence(dim(calibration.set)[1])){
        #Now organize the MRCA names and the two taxa that specify them from congruify output and into our treePL configurator:
        	cat("mrca =",calibration.set[mrca.index,1],calibration.set[mrca.index,4], calibration.set[mrca.index,5],"\n", file=paste("tree", tree.index, ".cppr8s", sep=""), sep=" ", append=TRUE)
        }
        #Now organize our minimum age estimates them from congruify output and into our treePL configurator:
        for(minage.index in sequence(dim(calibration.set)[1])){
        	cat("min =",calibration.set[minage.index,1],calibration.set[minage.index,3], "\n", file=paste("tree",tree.index, ".cppr8s", sep=""), sep=" ", append=TRUE)
        }
        #Now organize our maximum age estimates them from congruify output and into our treePL configurator:
        for(maxage.index in sequence(dim(calibration.set)[1])){
        	cat("max =",calibration.set[maxage.index,1],calibration.set[maxage.index,2], "\n", file=paste("tree",tree.index, ".cppr8s", sep=""), sep=" ", append=TRUE)
        }
        #############################################
        
        #Now some other bits. Note here you might want to alter the range of the smoothing parameters:
        cat("outfile=tree.dated.",tree.index,"\nprime\ncv\nrandomcv\nthorough\ncvstart=1000\ncvstop=0.0000001\nplsimaniter=200000\ncvsimaniter=50000\nlogpen\n", file=paste("tree",tree.index, ".cppr8s", sep=""), sep="", append=TRUE)

        #This now runs treePL based on the file you just made -- at this point, it is going to take a while -- so go get a beer and chillax.
        system(paste("treePL ", paste("tree", tree.index, ".cppr8s", sep=""), " > prime.out.", tree.index, sep=""))
    }
}


DoFullAnalysis <- function(calibration.set, target.tree, index){
    start.index <- end.index <- index
    ###The overall loop is looping over trees you are referencing:
    for (tree.index in start.index:end.index) {
        write.tree(target.tree[[tree.index]], paste("camp.rated.", index, ".tre", sep=""))
        cat("treefile =", paste("camp.rated.", index, ".tre", sep=""),"\nnumsites = 1515\n", file=paste("tree", tree.index, ".cppr8s", sep=""),sep="")
        
        ####### CALIBRATION BIT IS ADDED HERE #######
        for(mrca.index in sequence(dim(calibration.set)[1])){
            #Now organize the MRCA names and the two taxa that specify them from congruify output and into our treePL configurator:
            cat("mrca =",calibration.set[mrca.index,1],calibration.set[mrca.index,4], calibration.set[mrca.index,5],"\n", file=paste("tree", tree.index, ".cppr8s", sep=""), sep=" ", append=TRUE)
        }
        #Now organize our minimum age estimates them from congruify output and into our treePL configurator:
        for(minage.index in sequence(dim(calibration.set)[1])){
            cat("min =",calibration.set[minage.index,1],calibration.set[minage.index,3], "\n", file=paste("tree",tree.index, ".cppr8s", sep=""), sep=" ", append=TRUE)
        }
        #Now organize our maximum age estimates them from congruify output and into our treePL configurator:
        for(maxage.index in sequence(dim(calibration.set)[1])){
            cat("max =",calibration.set[maxage.index,1],calibration.set[maxage.index,2], "\n", file=paste("tree",tree.index, ".cppr8s", sep=""), sep=" ", append=TRUE)
        }
        #############################################
        
        ####### OTHER BITS HERE #######
        cat("outfile=tree.dated.",tree.index,"\ncv\nrandomcv\nthorough\ncvstart=0.00001\ncvstop=0.0000001\nplsimaniter=200000\ncvsimaniter=50000\nlogpen\nnthreads=2\n", file=paste("tree",tree.index, ".cppr8s", sep=""), sep="", append=TRUE)
        cat("cvoutfile=cv.",tree.index,"\n", file=paste("tree",tree.index, ".cppr8s", sep=""), sep="", append=TRUE)
        
        ##### ADD THE PRIMING OPTIONS HERE #######
        prime.output <- readLines(paste("prime.out.", tree.index, sep=""))
        prime.options.start <- which(prime.output %in% "PLACE THE LINES BELOW IN THE CONFIGURATION FILE")
        prime.options <- prime.output[(prime.options.start+1):length(prime.output)]
        for(prime.index in 1:length(prime.options)){
            cat(prime.options[prime.index], "\n", file=paste("tree", tree.index, ".cppr8s", sep=""), sep="", append=TRUE)
        }
        ###############################
        
        #This now runs treePL based on the file you just made -- at this point, it is going to take a while -- so go get a beer and chillax.
        system(paste("treePL", paste("tree", tree.index, ".cppr8s", sep="")))
    }
}


######################################################################################################################################
######################################################################################################################################
### Running Simulation
######################################################################################################################################
######################################################################################################################################

reference.tree.set = read.tree("Camp.mcc.plus100.tre")
reference.tree <- reference.tree.set[[1]]
#target.tree <- read.tree("camp.matk.rated.start0025.rate003.noshift.trees")

###################
#target.tree <- read.tree()
###################

tax.table <- read.csv("Campanulid.matK.reference.table.csv", as.is=TRUE, row=1)

#Analysis -- All calibrations, uniform "prior" at root.
tot.cals <- GetFullHashCalibrations(reference.tree, target.tree[[1]], tax.table)
subset.cals <- tot.cals[tot.cals$MaxAge<65,]
subset.cals <- rbind(tot.cals[1,], subset.cals)
subset.cals$MaxAge[1] = 150.0
subset.cals$MinAge[1] = 80.0

DoParallel <- function(index){
    DoPrimeAnalysis(subset.cals, target.tree, index)
    DoFullAnalysis(subset.cals, target.tree, index)
}

mclapply(1:50, DoParallel, mc.cores=20)

