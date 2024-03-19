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
# See help(read.tree) for the structure of the tree in R

## Plot the tree
plot(tree)
axisPhylo()

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
plot(tree, show.tip.label = FALSE, x.lim = 25)
phydataplot(dat$lat, tree, offset = 1, scaling = 0.1)
phydataplot(dat$long, tree, offset = 15, scaling = 0.05)

# Same plot, without the warnings : vectors need to be named.
lat <- dat$lat
long <- dat$long
names(lat) <- names(long) <- dat$traits
plot(tree, show.tip.label = FALSE, x.lim = 25)
phydataplot(lat, tree, offset = 1, scaling = 0.1)
phydataplot(long, tree, offset = 15, scaling = 0.05)

################################################################################
## treedata
################################################################################
library(treeio) # "tree input/output"
# This is to export the data to a format that can be read by EvoLaps
# It exports the data as an "extended newick"
# First create a table with tips and node reference numbers
rec_table <- list(node = seq_len(n_tips + n_nodes))

# # This is a very dumb ancestral reconstruction, just for the example : all ancestral nodes are at the same location
# dumb_ancestral <- c(mean(lat), mean(long))
# # Allocate locations to the data table (first tips, then nodes)
# rec_table[["location1"]] <- c(dat$lat, rep(dumb_ancestral[1], n_nodes))
# rec_table[["location2"]] <- c(dat$long, rep(dumb_ancestral[2], n_nodes))

# Estimateur MB

MB_ancestral <- mut
# # Allocate locations to the data table (first tips, then nodes)
# rec_table[["location1"]] <- c(dat$lat, rep(MB_ancestral[1], n_nodes))
# rec_table[["location2"]] <- c(dat$long, rep(MB_ancestral[2], n_nodes))

# Estimateur des noeuds intermédiaires
X <- Estimateur
rec_table[["location1"]] <- c(dat$lat, MB_ancestral[1], Estimateur[1:102])
rec_table[["location2"]] <- c(dat$long, MB_ancestral[2], Estimateur[103:204])

# Estimateur MB avec dérive

MB_ancestral_d <- theta[1,]
# Estimateur des noeuds intermédiaires
rec_table[["location1"]] <- c(dat$lat, MB_ancestral_d[1], Estimateur_d[1:102])
rec_table[["location2"]] <- c(dat$long, MB_ancestral_d[2], Estimateur_d[103:204])

# Format the data to export it
rec_table <- as_tibble(rec_table)
tree_tibble <- as_tibble(tree)
tree_data <- full_join(tree_tibble, rec_table, by = 'node')
tree_data <- as.treedata(tree_data)
# Write the extended newick. The resulting file can be read by Evolaps.
write.beast(tree_data, file = here("results", "tree_MB_derive.tree"), tree.name = "TREE_MB_derive")
