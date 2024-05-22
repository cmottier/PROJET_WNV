library(here) # To deal with paths : project must be root, and data must be in a folder named "data".
library(ape)  # R package for phylogenetic trees

################################################################################
## Tree
################################################################################
## Import the tree
tree <- read.tree(here("data", "WNV_Pybus2012_MCC.newick"))
str(tree)
n_nodes <- tree$Nnode
n_tips <- length(tree$tip.label)

################################################################################
## Data
################################################################################
## Import the data
dat <- read.table(here("data", "WNV_lat_long.txt"), header = TRUE)
head(dat)

## Match
# Match the data with the tips of the tree, so that they are in the same order
match(tree$tip.label, dat$traits)
dat <- dat[match(tree$tip.label, dat$traits), ]
match(tree$tip.label, dat$traits)

## Plot tree and data
lat <- dat$lat
long <- dat$long
names(lat) <- names(long) <- dat$traits
plot(tree, show.tip.label = FALSE, x.lim = 25)
phydataplot(lat, tree, offset = 3, scaling = 0.1)
mtext("Latitude",1,at=15,line=3)
phydataplot(long, tree, offset = 17, scaling = 0.05)
mtext("Longitude",1,at=23,line=3)
axisPhylo(backward = F)
# axisPhylo(backward = F, root.time=1998)
