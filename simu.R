set.seed(14)

### EXEMPLE DE MVMORPH

# Generating a random tree with 50 species
tree_test<-pbtree(n=50)

tree_test

# Setting the regime states of tip species
sta<-as.vector(c(rep("Forest",20),rep("Savannah",30))); names(sta)<-tree_test$tip.label

# Making the simmap tree with mapped states
tree_test<-make.simmap(tree_test, sta , model="ER", nsim=1)

tree_test

# Number of simulated datasets
nsim<-1

# Rates matrices for the "Forest" and the "Savannah" regimes

# Note: use lists for multiple rates (matrices or scalars)
sigma<-list(Forest=matrix(c(2,0.5,0.5,1),2), Savannah=matrix(c(5,3,3,4),2))

sigma

# ancestral states for each traits
theta<-c(0,0)

# Simulate the "BMM" model
simul_1<-mvSIM(tree_test, nsim=nsim, model="BMM", param=list(sigma=sigma, theta=theta))
View(simul_1)
head(simul_1)

# fit the BMM model on simulated data
fit <- mvBM(tree_test, simul_1)

# simulate 100 datasets from the fitted object
simul_2 <- simulate(fit, tree=tree_test, nsim=100)

# parametric bootstrap; e.g.:
bootstrap <- lapply(simul_2, function(x) mvBM(tree_test, x, echo=F, diagnostic=F))

# retrieve results; e.g. for the log-likelihood
log_distribution <- sapply(bootstrap, logLik)
hist(log_distribution, main="Log-likelihood distribution")
abline(v=fit$LogLik, lty=2, lwd=5, col="red")

#tree_test



##### TEST #####

## matrice de covariance R
Sigma <- matrix(c(1,-1,-1,2),nrow = 2)
Sigma

## simulation d'un brownien (équivalent à dat pour nous)
sim <- mvSIM(tree = tree, nsim = 1, model = "BM1",param = list(sigma=Sigma, theta=theta))


### appliquer notre modèle à cette situation

n <- 104

# matrice des temps d'évolution partagés
C_test <- vcv(tree)

# inverse
C_test_inv <- solve(C_test)

# vecteur 1n
un <- matrix(rep(1,n), ncol = 1)

# matrice Y (sans dérive)
Y<-sim

# vecteur temps
Temps<-diag(C_test)

# vecteur dérive (déterministe)
b = c(2,3)

# simulation avec dérive
sim_derive <- mvSIM(tree = tree, nsim = 1, model = "BM1",param = list(sigma=Sigma, theta=theta)) + t(b%*%t(Temps))

# matrice Y (avec dérive)
Y_derive <- as.matrix(sim_derive)

X <- matrix(cbind(un,Temps),ncol = 2)
X

# estimation modèle sans dérive
mu_test <- solve(t(un)%*%C_test_inv%*%un)%*%t(un)%*%C_test_inv%*%Y
mu_test
R_test <- 1/103*t(Y-un%*%mu_test)%*%C_test_inv%*%(Y-un%*%mu_test)
R_test

# calcul d'erreur
sum((theta - mu_test)^2)
(theta - mu_test)

simu <- mvSIM(tree = tree, nsim = 10000, model = "BM1",param = list(sigma=Sigma, theta=theta))
err<-NULL
for(i in 1:10000){
  Y <- simu[[i]]
  mu <- solve(t(un)%*%C_test_inv%*%un)%*%t(un)%*%C_test_inv%*%Y
  #R <- 1/103*t(Y-un%*%mu_test)%*%C_test_inv%*%(Y-un%*%mu_test)
  err[i]<-sum((theta - mu)^2)
}
mean(err)

# estimation modèle avec dérive
theta_hat <- solve(t(X)%*%C_test_inv%*%X)%*%t(X)%*%C_test_inv%*%Y_derive
R_hat_d <- 1/102*t(Y_derive)%*%C_test_inv%*%(Y_derive-X%*%theta_hat)
theta_hat
R_hat_d

# calcul d'erreur

mu_derive <- theta_hat[1,]
sum(mu_derive^2)



dim(theta_hat)
dim(Y_derive)
dim(C_test_inv)
dim(X%*%theta)

fit <- mvBM(tree,data = dat2)

dat2<-dat[,c(2,3)]
rownames(dat2)  <- dat[,1]
View(dat2)

summary(fit)

