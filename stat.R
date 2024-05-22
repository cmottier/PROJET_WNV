##################################################################
## import des données (les mêmes que dans WNV_DATA)
##################################################################
library(ape)
library(here)

tree <- read.tree(here("data", "WNV_Pybus2012_MCC.newick"))
dat <- read.table(here("data", "WNV_lat_long.txt"), header = TRUE)
dat <- dat[match(tree$tip.label, dat$traits), ]

C <- vcv(tree)
invC <- solve(C)
Un <- matrix(rep(1,104),ncol=1)
Y <- as.matrix(dat[,2:3])

mut <- solve(t(Un)%*%invC%*%Un)%*%t(Un)%*%invC%*%Y
R <- 1/103*t(Y-Un%*%mut)%*%invC%*%(Y-Un%*%mut)

##################################################################
##################################################################



### première idée: simuler un brownien, fitter suivant le modèle brownien et faire un boxplot des erreurs

## paramètres du brownien simulé

# variance
sigma <- matrix(R,nrow = 2) #on utilise la matrice R estimée

# moyenne (ancetre)
mu = c(1,1)

# dérive
b <- c(1,2)

## paramètres pour l'estimation
C <- vcv(tree)
invC<-solve(C)
Un <- matrix(rep(1,104),ncol = 1)

## simulations sans dérive:
sim <- mvSIM(tree = tree, nsim = 10000, model = "BM1", param = list(sigma = sigma, theta = mu))

err <- NULL
err_mu1 <- NULL
err_mu2 <- NULL
err_norm <- NULL
mu1 <- NULL
mu2 <- NULL

for(i in 1:10000){
  
  y <- sim[[i]] # simulation n°i
  mu_est <- solve(t(Un)%*%invC%*%Un)%*%t(Un)%*%invC%*%y #estimateur de mu (ancêtre)
  mu1[i] <- mu_est[1]
  mu2[i] <- mu_est[2]
  err_mu1[i] <- (mu[1] - mu_est[1]) # écart sur la première composante
  err_mu2[i] <- (mu[2] - mu_est[2]) # écart sur la deuxième composante
  err[i]<- sqrt(sum((mu - mu_est)^2)) # 'erreur globale"
  err_norm[i] <- sqrt((sum((mu-mu_est)^2))) / sqrt(sum(mu^2)) #erreur normalisée en norme 2
  
}

## boxplot de l'erreur normalisée
boxplot(err_norm)


#################################################
## boxplots sur des arbres de tailles différentes
#################################################

### arbre à 10 feuilles
tree10 <- rtree(n=10,equiprob = FALSE) # simulation d'un arbre à 10 feuilles
MB10 <- mvSIM(tree = tree10, nsim = 10000, model = "BM1", param = list(sigma = sigma, theta = mu)) #simulation du brownien à partir de l'abre simulé

C10 <- vcv(tree10)
U10 <- matrix(rep(1,10),ncol = 1)
invC10 <- solve(C10)
mu10 <- list()
mu1_10 <- NULL
mu2_10 <- NULL
err_norm10 <- NULL
for(i in 1:10000){
  mu10[[i]] <- solve(t(U10)%*%invC10%*%U10)%*%t(U10)%*%invC10%*%MB10[[i]]
  mu1_10[i] <- mu10[[i]][1]
  mu2_10[i] <- mu10[[i]][2]
  err_norm10[i] <-sqrt(sum((mu - mu10[[i]])^2)) / sqrt(sum(mu^2))
}



### arbre à 50 feuilles
tree50 <- rtree(n = 50, equiprob = FALSE)
MB50 <- mvSIM(tree = tree50, nsim = 10000, model = "BM1", param = list(sigma = sigma, theta = mu))
C50 <- vcv(tree50)
U50 <- matrix(rep(1,50),ncol = 1)
invC50 <- solve(C50)
mu50 <- list()
mu1_50 <- NULL
mu2_50 <- NULL
err_norm50 <- NULL
for(i in 1:10000){
  mu50[[i]] <- solve(t(U50)%*%invC50%*%U50)%*%t(U50)%*%invC50%*%MB50[[i]]
  mu1_50[i] <- mu50[[i]][1]
  mu2_50[i] <- mu50[[i]][2]
  err_norm50[i] <-sqrt(sum((mu - mu50[[i]])^2)) / sqrt(sum(mu^2))
}


### arbre à 100 feuilles
tree100 <- rtree(n = 100, equiprob = FALSE)
MB100 <- mvSIM(tree = tree100, nsim = 10000, model = "BM1", param = list(sigma = sigma, theta = mu))
C100 <- vcv(tree100)
U100 <- matrix(rep(1,100),ncol = 1)
invC100 <- solve(C100)
mu100 <- list()
mu1_100 <- NULL
mu2_100 <- NULL
err_norm100 <- NULL
for(i in 1:10000){
  mu100[[i]] <- solve(t(U100)%*%invC100%*%U100)%*%t(U100)%*%invC100%*%MB100[[i]]
  mu1_100[i] <- mu100[[i]][1]
  mu2_100[i] <- mu100[[i]][2]
  err_norm100[i] <-sqrt(sum((mu - mu100[[i]])^2)) / sqrt(sum(mu^2))
}


### arbre à 150 feuilles
tree150 <- rtree(n = 150, equiprob = FALSE)
MB150 <- mvSIM(tree = tree150, nsim = 10000, model = "BM1", param = list(sigma = sigma, theta = mu))
C150 <- vcv(tree150)
U150 <- matrix(rep(1,150),ncol = 1)
invC150 <- solve(C150)
mu150 <- list()
mu1_150 <- NULL
mu2_150 <- NULL
err_norm150 <- NULL
for(i in 1:10000){
  mu150[[i]] <- solve(t(U150)%*%invC150%*%U150)%*%t(U150)%*%invC150%*%MB150[[i]]
  mu1_150[i] <- mu150[[i]][1]
  mu2_150[i] <- mu150[[i]][2]
  err_norm150[i] <-sqrt(sum((mu - mu150[[i]])^2)) / sqrt(sum(mu^2))
}


### arbre à 200 feuilles
tree200 <- rtree(n = 200, equiprob = FALSE)
MB200 <- mvSIM(tree = tree200, nsim = 10000, model = "BM1", param = list(sigma = sigma, theta = mu))
C200 <- vcv(tree200)
U200 <- matrix(rep(1,200),ncol = 1)
invC200 <- solve(C200)
mu200 <- list()
mu1_200 <- NULL
mu2_200 <- NULL
err_norm200 <- NULL
for(i in 1:10000){
  mu200[[i]] <- solve(t(U200)%*%invC200%*%U200)%*%t(U200)%*%invC200%*%MB200[[i]]
  mu1_200[i] <- mu200[[i]][1]
  mu2_200[i] <- mu200[[i]][2]
  err_norm200[i] <-sqrt(sum((mu - mu200[[i]])^2)) / sqrt(sum(mu^2))
}

data_error <- data.frame(cbind(err_norm10,err_norm50,err_norm100,err_norm150,err_norm200))

## boxplot des erreurs normalisées par nombre de feuilles
boxplot(data_error,main = "Erreurs normalisées par nombre de feuilles")
