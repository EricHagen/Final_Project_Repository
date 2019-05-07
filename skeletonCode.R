

### Skeleton code for the pipeline of crazy.

## Step 1 True -- Get True tree ##
# I already did this for you. It's basically the matK tree of campanulids, dated and such.
library(ape)
setwd("STEP_1_TRUE")
tree <- read.tree(file = "tree.dated.matK.tre")

## Step 2 Rating -- In this step you are going to take the tree in step 1 and turn into 
#a molecular rated tree that contains realistic heterogeneity among taxa. So, the input 
#are two files: 1) the same tree in step 1 but with growth habit designated at a node 
#(1=woody, 2=herbaceous), 2) a trait file that designates the woody and herbaceous 
#states for each tip (now, why would I do that?!!?!? -> Do woody & herbaceous states 
#affect diversification rates differently? Can we estimate diversification rates using
#trees at all @Rabosky? SpEcIeS sElEcTiOn???). I am applying an autocorrelated model of
#heterogeneity. There is a baseline rate of 0.03 which "evolves" on the tree. However, 
#when there is a shift to herbaceous the rate is increased 3x based on Smith and Donoghue 
#(2008). But you really should consult my paper on the specifics on this.

setwd("STEP_2_RATING") #Just use the rated trees that I published, but this is how 
#you'd run that code:
source("STEP_2_RATING/rateEvolve.R")
tree_woodyherbaceous <- read.tree("STEP_2_RATING/camp.matk.habit.recon.asterwood.carrotwood.tre")
trait <- read.delim("STEP_2_RATING/camp.matk.habit.trait.asterwood.carrotwood.txt")
rateEvolve(tree_woodyherbaceous, trait, type="autocor", nsims=10)


## Step 3 Simulate and align data, infer tree -- In this step, you will use the rated 
#trees in step 2 and then simulate nucleotide data that matches the empirical data as 
#best we can. It will then take these alignments and use RAxML to infer the tree resulting 
#from each. Note that in this step we will doing TWO versions of the same alignment --
#1) where it is a full matrix, 2) where we chunk it up so that it's messy and looks like 
#data you'd get from GenBank (why would you do that?!!?!? -> to control for whether 
#missing data is driving the answer). The RAxML lines are 53 and 63. BE SURE TO KNOW 
#WHAT THESE LINES ARE DOING. The output for this will be the trees for step 4.

setwd("STEP_3_SIMULATE_DATA/")
source("alignSimWrapper.R")
rate.phy.set <- read.tree("../STEP_2_RATING/camp.matk.rated.start0025.rate003.3xshift.trees")
alignSimWrap(start.rep=1, end.rep=1, rate.phy.set=rate.phy.set, nsites=1515)

##Step 4 Dating the tree -- Let's wait until you finish step 3 before moving forward.
## MORE SOON.
