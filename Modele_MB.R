################################################################################
# Modèle brownien simple
# Y ~ MN(Un*mut, Cn, R)
# Y = Un*mut+E où E ~ MN(0,Cn,R)
################################################################################

n <- 104
m <- 102

################################################################################
# Estimations
################################################################################

# Construction de Cn, basée sur les feuilles (temps d'évolution commune)
Cn <- vcv(tree)

# Estimateurs de mu et R
Un <- matrix(rep(1,n),ncol=1)
Y <- as.matrix(dat[,2:3])

mut_hat <- solve(t(Un)%*%solve(Cn)%*%Un)%*%t(Un)%*%solve(Cn)%*%Y
R_hat <- 1/(n-1)*t(Y-Un%*%mut_hat)%*%solve(Cn)%*%(Y-Un%*%mut_hat)

# Construction de Cmn, entre nœuds internes et feuilles

Cmn <- matrix(nrow=m, ncol=n) 
for (i in 1:m) {
  for (j in (1:n)) {
    ancetre <- getMRCA(tree,c(i+n+1,j))
    tij <- node.depth.edgelength(tree)[ancetre]
    Cmn[i,j] <- tij
  }
}

# Estimateur des positions aux nœuds internes
Um <- matrix(rep(1,m),ncol=1)
Z_hat <- Um%*%mut_hat + Cmn%*%solve(Cn)%*%(Y-Un%*%mut_hat)


################################################################################
# Représentation Evolaps
# This is to export the data to a format that can be read by EvoLaps
# It exports the data as an "extended newick"
################################################################################

library(ellipse)
library(treeio)

# First create a table with tips and node reference numbers
rec_table <- list(node = seq_len(n_tips + n_nodes))

# Positions des nœuds estimées
rec_table[["location1"]] <- c(dat$lat, mut_hat[1], Z_hat[,1])
rec_table[["location2"]] <- c(dat$long, mut_hat[2], Z_hat[,2])

# Ellipses de confiance
# mut_hat = AY
A <- solve(t(Un)%*%solve(Cn)%*%Un)%*%t(Un)%*%solve(Cn)
# Z_hat = QY
Q <- Um%*%A + Cmn%*%solve(Cn)%*%(diag(n)-Un%*%A)

# Ellipse de confiance d'un nœud interne (j=1 latitude, j=2 longitude)
level <- 0.80
ell <- function(i,j) {
  if (i==105) {
    res = ellipse(kronecker(R_hat,A%*%Cn%*%t(A)), center = mut_hat, level = level)
  }
  else {
    res = ellipse(kronecker(R_hat,(Q%*%Cn%*%t(Q))[i,i]), center = Z_hat[i,], level = level)
  }
  return(res[,j])
}

# Ajout des ellipses de confiance dans rec_table
trait_name_lat <- paste0("location1_", level * 100, "%HPD")
trait_name_long <- paste0("location2_", level * 100, "%HPD")
rec_table[[paste0(trait_name_lat, "_", 1)]] <- c(rep(NA, n), lapply(c(105,(1:m)), function(i) ell(i,1)))
rec_table[[paste0(trait_name_long, "_", 1)]] <- c(rep(NA, n), lapply(c(105,(1:m)), function(i) ell(i,2)))

# Format the data to export it
rec_table <- as_tibble(rec_table)
tree_tibble <- as_tibble(tree)
tree_data <- full_join(tree_tibble, rec_table, by = 'node')
tree_data <- as.treedata(tree_data)
# Write the extended newick. The resulting file can be read by Evolaps.
write.beast(tree_data, file = here("results", "tree_MB_ellipses.tree"), tree.name = "TREE_MB_ellipses")


################################################################################
# AIC et BIC
################################################################################

l_hat <- -n*log(2*pi)-n/2*log(det(R_hat))-log(det(Cn))-1/2*sum(diag(solve(R_hat)%*%t(Y-Un%*%mut_hat)%*%solve(Cn)%*%(Y-Un%*%mut_hat)))
AIC <- -2*l_hat + 2*5
BIC <- -2*l_hat + log(104)*5
AIC
BIC


