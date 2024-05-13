################################################################################
# Modèle brownien avec dérive b
# Y ~ MN(Un*mu^t+T*b^t, C, R) 
# Y = X*theta+E où X =[Un|t], theta^t = [mu|b] et E ~ MN(0,C,R)
################################################################################

# Construction de C, basée sur les feuilles (temps d'évolution commune)

C <- vcv(tree)

# Matrice T des temps des feuilles

T <- matrix(nrow=104, ncol=1)
for (i in 1:104) {
  T[i] <- C[i,i]
}

# Estimateurs de theta (ie mu et b) et R

invC <- solve(C)
Un <- matrix(rep(1,104),ncol=1)
X <- cbind(Un,T)
Y <- as.matrix(dat[,2:3])

theta <- solve(t(X)%*%invC%*%X)%*%t(X)%*%invC%*%Y
R_d <- 1/102*t(Y)%*%invC%*%(Y-X%*%theta)


# ################################################################################
# # Estimateur des valeurs des ancêtres avec E[Z|Y]=mZ + SigmaZY*SigmaY^(-1)*(Y-mY)
# 
# # Construction de Cnoeud, basée sur tous les noeuds
# 
# f <- function(x) {
#   if (x<105) return(x)
#   else return(x+1)
# }
# 
# Cnoeud <- matrix(nrow=206, ncol=206) # On enlève la racine 105
# for (i in 1:206) {
#   Cnoeud[i,i] <- node.depth.edgelength(tree)[f(i)]
# }
# for (i in 1:205) {
#   for (j in ((i+1):206)) {
#     ancetre <- getMRCA(tree,c(f(i),f(j)))
#     tij <- node.depth.edgelength(tree)[ancetre]
#     Cnoeud[i,j] <- tij
#     Cnoeud[j,i] <- tij
#   }
# }
# 
# # Matrice Tm des temps des noeuds internes
# 
# Tm <- matrix(nrow=102, ncol=1)
# for (i in 1:102) {
#   Tm[i] <- Cnoeud[i+104,i+104]
# }
# 
# # Matrice de variance-covariance de (X,Y)
# RC_d <- kronecker(R_d,Cnoeud, FUN = "*")
# 
# # Mettre sous forme vec
# Um <- matrix(rep(1,102),ncol=1)
# vecY <- matrix(Y, ncol=1, byrow = F)
# mY_d <- matrix(Un%*%theta[1,]+T%*%theta[2,], ncol=1, byrow =F)
# mX_d <- matrix(Um%*%theta[1,]+Tm%*%theta[2,], ncol=1, byrow=F)
# 
# # Extraction des matrices SigmaXY et SigmaY
# SigmaXY_d <- rbind(cbind(RC_d[105:206,1:104],RC_d[105:206,207:310]),cbind(RC_d[311:412,1:104],RC_d[311:412,207:310]))
# SigmaY_d <- rbind(cbind(RC_d[1:104,1:104],RC_d[1:104,207:310]),cbind(RC_d[207:310,1:104],RC_d[207:310,207:310]))
# 
# # Estimateur
# Estimateur_d <- mX_d + SigmaXY_d%*%solve(SigmaY_d)%*%(vecY-mY_d)

################################################################################
# Estimateur des valeurs des ancêtres avec E[Z|Y]=mZ + I2@(CmnCn^(-1))*(Y-mY)

# Construction de Cmn entre nœuds internes et feuilles

Cmn <- matrix(nrow=102, ncol=104) 
for (i in 1:102) {
  for (j in (1:104)) {
    ancetre <- getMRCA(tree,c(i+105,j))
    tij <- node.depth.edgelength(tree)[ancetre]
    Cmn[i,j] <- tij
  }
}

# Matrice Tm des temps des noeuds internes

Tm <- matrix(nrow=102, ncol=1)
for (i in 1:102) {
  Tm[i] <- Cnoeud[i+104,i+104]
}

# Mettre sous forme vec
vecY <- matrix(Y, ncol=1, byrow = F)
mY_d <- matrix(X%*%theta, ncol=1, byrow =F)
Um <- matrix(rep(1,102),ncol=1)
mZ_d <- matrix(cbind(Um,Tm)%*%theta, ncol=1, byrow=F)


# Estimateur
Estimateur_d <- mZ_d + kronecker(diag(2),(Cmn%*%solve(C)), FUN = "*")%*%(vecY-mY_d)



################################################################################
# Calcul de l'AIC ?

L_hat_d <- -103 * log(2*pi) - 103/2 * log(det(R_d)) - 1/2*sum(diag(solve(R_d)%*%t(Y-X%*%theta)%*%(Y-X%*%theta)))
AIC_d <- -2*L_hat_d + 2*8
# 2933
