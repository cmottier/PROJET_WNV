################################################################################
# Modèle brownien avec dérive b
# Y ~ MN(Un*mu^t+T*b^t, Cn, R) 
# Y = X*theta+E où X =[Un|T], theta^t = [mu|b] et E ~ MN(0,Cn,R)
################################################################################

n <- 104
m <- 102

################################################################################
# Estimations
################################################################################

# Construction de C, basée sur les feuilles (temps d'évolution commune)
Cn <- vcv(tree)

# Matrice T des temps d'évolution des feuilles
T <- matrix(nrow=n, ncol=1)
for (i in 1:n) {
  T[i] <- C[i,i]
}

# Estimateurs de theta (ie mu et b) et R
Un <- matrix(rep(1,n),ncol=1)
X <- cbind(Un,T)
Y <- as.matrix(dat[,2:3])

theta_hat_d <- solve(t(X)%*%solve(Cn)%*%X)%*%t(X)%*%solve(Cn)%*%Y
R_hat_d <- 1/(n-2)*t(Y)%*%solve(Cn)%*%(Y-X%*%theta_hat_d)

# Construction de Cmn, entre nœuds internes et feuilles
Cmn <- matrix(nrow=m, ncol=n) 
for (i in 1:m) {
  for (j in (1:n)) {
    ancetre <- getMRCA(tree,c(i+n+1,j))
    tij <- node.depth.edgelength(tree)[ancetre]
    Cmn[i,j] <- tij
  }
}

# Matrice Tm des temps d'évolution des nœuds internes
Tm <- matrix(nrow=m, ncol=1)
for (i in 1:m) {
  Tm[i] <- Cnoeud[i+n,i+n]
}

# Matrice Xm des prédicteurs pour les nœuds internes
Um <- matrix(rep(1,m),ncol=1)
Xm <- cbind(Um,Tm)

# Estimateur des positions des nœuds internes
Z_hat_d <- Xm%*%theta_hat_d + Cmn%*%solve(Cn)%*%(Y-X%*%theta_hat_d)


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
rec_table[["location1"]] <- c(dat$lat, theta_hat_d[1,1], Z_hat_d[,1])
rec_table[["location2"]] <- c(dat$long, theta_hat_d[1,2], Z_hat_d[,2])

# Ellipses de confiance
# theta_hat = AY
A_d <- solve(t(X)%*%solve(Cn)%*%X)%*%t(X)%*%solve(Cn)
# Z_hat = QY
Q_d <- Xm%*%A_d + Cmn%*%solve(Cn)%*%(diag(n)-X%*%A_d)

# Ellipse de confiance d'un nœud interne (j=1 latitude, j=2 longitude)
level <- 0.80
ell_d <- function(i,j) {
  if (i==105) {
    res = ellipse(kronecker(R_hat_d,(A_d%*%Cn%*%t(A_d))[1,1]), center = theta_hat_d[1,], level = level)
  }
  else {
    res = ellipse(kronecker(R_hat_d,(Q_d%*%Cn%*%t(Q_d))[i,i]), center = Z_hat_d[i,], level = level)
  }
  return(res[,j])
}

# Ajout des ellipses de confiance dans rec_table
trait_name_lat <- paste0("location1_", level * 100, "%HPD")
trait_name_long <- paste0("location2_", level * 100, "%HPD")
rec_table[[paste0(trait_name_lat, "_", 1)]] <- c(rep(NA, n), lapply(c(105,(1:m)), function(i) ell_d(i,1)))
rec_table[[paste0(trait_name_long, "_", 1)]] <- c(rep(NA, n), lapply(c(105,(1:m)), function(i) ell_d(i,2)))

# Format the data to export it
rec_table <- as_tibble(rec_table)
tree_tibble <- as_tibble(tree)
tree_data <- full_join(tree_tibble, rec_table, by = 'node')
tree_data <- as.treedata(tree_data)
# Write the extended newick. The resulting file can be read by Evolaps.
write.beast(tree_data, file = here("results", "tree_MB_d_ellipses.tree"), tree.name = "TREE_MB_d_ellipses")


################################################################################
# AIC, BIC et test de rapport de vraisemblance
################################################################################

l_hat_d <- -n*log(2*pi)-n/2*log(det(R_hat_d))-log(det(Cn))-1/2*sum(diag(solve(R_hat_d)%*%t(Y-X%*%theta_hat_d)%*%solve(Cn)%*%(Y-X%*%theta_hat_d)))
AIC_d <- -2*l_hat_d + 2*7
BIC_d <- -2*l_hat_d + log(104)*7
AIC_d
BIC_d

# Test de rapport de vraisemblance : (H_0) dérive nulle
lambda <- -2*(l_hat-l_hat_d)
lambda
lambda < qchisq(0.95, 2) 
lambda < qchisq(0.975, 2)
