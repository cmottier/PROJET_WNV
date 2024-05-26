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

# fonction pour simuler un arbre et un/des brownien(s) à partir de cet arbre
sim_error <- function(n,sigma = R,mu = c(1,1), precision = 10000){
  set.seed(674)
  tree <- rtree(n = n, equiprob = FALSE) #simuler l'arbre
  tree$edge.length <- tree$edge.length/max(vcv(tree)) #renormaliser la variance
  
  MB_sim <- mvSIM(tree = tree, nsim = precision, model = "BM1", param = list(sigma = sigma, theta = mu))
  C <- vcv(tree)
  Un <- matrix(rep(1,n),ncol = 1)
  invC <- solve(C)
  
  mu_est <- NULL
  R <- list()
  mu1 <- NULL
  mu2 <- NULL
  err_mu1 <- NULL
  err_mu2 <- NULL
  for(i in 1:precision){
    mu_est <- solve(t(Un)%*%invC%*%Un)%*%t(Un)%*%invC%*%MB_sim[[i]]
    mu1[i] <- mu_est[1]
    mu2[i] <- mu_est[2]
    err_mu1[i] <- (mu1[i] - mu[1]) / mu[1]
    err_mu2[i] <- (mu2[i] - mu[2]) / mu[2]
    R[[i]] <- 1/(n-1)*t(MB_sim[[i]]-Un%*%mu_est)%*%invC%*%(MB_sim[[i]]-Un%*%mu_est)
  }
  

  est_moy <- c(mean(mu1),mean(mu2))
  
  res <- list(est_moy,mu1,mu2 ,err_mu1 ,err_mu2,R)
  
  return(res)
  
}

# simuler des arbres, MB sur 10, 50, 100, 200, 250 feuilles

sim10 <- sim_error(n = 10)
sim50 <- sim_error(n = 50)
sim100 <- sim_error(n = 100)
sim200 <- sim_error(n = 200) 
sim250 <- sim_error(n = 250) 

data_err1 <- data.frame(cbind(sim10[[4]] , sim50[[4]] , sim100[[4]], sim200[[4]], sim250[[4]]))
colnames(data_err1) <- c("10","50","100","200","250")
boxplot(data_err1)


## observer l'évolution de l'erreur en fonction du nombre de feuilles
# peut prendre du temps à compiler
err1 <- NULL
for(i in 10:350){
  err1[i] <- mean(abs(sim_error(n=i, precision = 200)[[4]]))
}

plot(x = 10:350, y = err1[10:350], type = "l")


### pareil pour le MB avec dérive
sim_errord <- function(n,sigma = R,mu = c(1,1), precision = 10000, b =c(1,2)){
  set.seed(674)
  tree <- rtree(n = n, equiprob = FALSE)
  tree$edge.length <- tree$edge.length/max(vcv(tree))
  MB_sim <- mvSIM(tree = tree, nsim = precision, model = "BM1", param = list(sigma = sigma, theta = mu, trend = b))
  
  C <- vcv(tree)
  Un <- matrix(rep(1,n),ncol = 1)
  T<-matrix(diag(C),ncol = 1)
  invC <- solve(C)
  X <- cbind(Un,T)
  
  theta_est <- NULL
  R <- list()
  mu1 <- NULL
  mu2 <- NULL
  err_mu1 <- NULL
  err_mu2 <- NULL
  
  for(i in 1:precision){
    theta_est <- solve(t(X)%*%invC%*%X)%*%t(X)%*%invC%*%MB_sim[[i]]
    mu1[i] <- theta_est[1,1]
    mu2[i] <- theta_est[1,2]
    err_mu1[i] <- (mu1[i] - mu[1]) / mu[1]
    err_mu2[i] <- (mu2[i] - mu[2]) / mu[2]
    R[[i]] <- 1/(n-2)*t(MB_sim[[i]])%*%invC%*%(MB_sim[[i]]-X%*%theta_est)
  }
  
  
  est_moy <- c(mean(mu1),mean(mu2))
  
  res <- list(est_moy,mu1,mu2 ,err_mu1 ,err_mu2,R)
  
  return(res)
  
}

sim10d <- sim_errord(n=10)
sim50d <- sim_errord(n = 50)
sim100d <- sim_errord(n = 100)
sim200d <- sim_errord(n = 200)
sim250d <- sim_errord(n = 250)

data_err1d <- data.frame(cbind(sim10d[[4]], sim50d[[4]], sim100d[[4]], sim200d[[4]], sim250d[[4]]))
colnames(data_err1d) <- c("10","50","100","200","250")
boxplot(data_err1d)

#################################################################################################
### stats sur notre arbre
#################################################################################################

simWNV <- function(tree = tree_WNV,nsim = 1000,sigma = R,mu = c(1,1)){
  set.seed(674)
  MBsim <- mvSIM(tree = tree, nsim = nsim, model = "BM1", param = list(sigma = sigma, theta = mu))
  mu_est <- NULL
  mu1 <- NULL
  mu2 <- NULL
  err1 <- NULL
  err2 <- NULL
  C<-vcv(tree)
  invC <- solve(C)
  Un <- matrix(rep(1,length(tree$tip.label)),ncol = 1)
  for(i in 1:nsim){
    mu_est  <- solve(t(Un)%*%invC%*%Un)%*%t(Un)%*%invC%*%MBsim[[i]]
    mu1[i] <- mu_est[1]
    mu2[i] <- mu_est[2]
    err1[i] <- (mu_est[1] - mu[1]) / mu[1]
    err2[i] <- (mu_est[2] - mu[2]) / mu[2]
  }
  
  mu_moy <- c(mean(mu1),mean(mu2))
  
  res <- list(err1 = err1,err2 = err2, moy = mu_moy)
  
  return(res)
}

sim <- simWNV(nsim = 1000)

boxplot(sim$err1)
plot(sim$err1,type = "l", xlab = "Nombre de simulations", ylab = "Erreur")
plot(sim$err2,type = "l", xlab = "Nombre de simulations", ylab = "Erreur")

#stat sur les erreurs
summary(sim$err1)
summary(sim$err2)

sim$moy #estimation moyenne

### pareil avec le MB avec dérive

simWNVd <- function(tree = tree_WNV,nsim = 1000,sigma = R,mu = c(1,1), b = c(1,2)){
  set.seed(674)
  MBsim <- mvSIM(tree = tree, nsim = nsim, model = "BM1", param = list(sigma = sigma, theta = mu, trend = b))
  theta_est <- NULL
  mu1 <- NULL
  mu2 <- NULL
  err1 <- NULL
  err2 <- NULL
  C<-vcv(tree)
  invC <- solve(C)
  T <- matrix(diag(C),ncol=1)
  Un <- matrix(rep(1,length(tree$tip.label)),ncol = 1)
  X <- cbind(Un,T)
  for(i in 1:nsim){
    theta_est  <- solve(t(X)%*%invC%*%X)%*%t(X)%*%invC%*%MBsim[[i]]
    mu1[i] <- theta_est[1,1]
    mu2[i] <- theta_est[1,2]
    err1[i] <- (mu1[i] - mu[1]) / mu[1]
    err2[i] <- (mu2[i] - mu[2]) / mu[2]
  }
  
  mu_moy <- c(mean(mu1),mean(mu2))
  
  res <- list(err1 = err1, err2 = err2, moy = mu_moy,theta = theta_est)
}

sim_d <- simWNVd(n =1000)

sim_d$moy # moins bonne estimation ?

plot(sim_d$err1, type = "l", xlab = "Nombre de simulation", ylab = "Erreur")
plot(sim_d$err2, type = "l", xlab = "Nombre de simulation", ylab = "Erreur")

summary(sim_d$err1)
summary(sim_d$err2)
# l'erreur semble avoir une plus grande variance