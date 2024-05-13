################################################################################
# Modèle brownien simple
# Y ~ MN(Un*mut, C, R)
# Y = Un*mut+E où E ~ MN(0,C,R)
################################################################################

# Construction de C, basée sur les feuilles (temps d'évolution commune)

C <- vcv(tree)

# C <- matrix(nrow=104, ncol=104)
# for (i in 1:104) {
#   C[i,i] <- node.depth.edgelength(tree)[i]
# }
# for (i in 1:103) {
#   for (j in ((i+1):104)) {
#     ancetre <- getMRCA(tree,c(i,j))
#     tij <- node.depth.edgelength(tree)[ancetre]
#     C[i,j] <- tij
#     C[j,i] <- tij
#   }
# }

# Estimateurs de mu et R

invC <- solve(C)
Un <- matrix(rep(1,104),ncol=1)
Y <- as.matrix(dat[,2:3])

mut <- solve(t(Un)%*%invC%*%Un)%*%t(Un)%*%invC%*%Y
R <- 1/103*t(Y-Un%*%mut)%*%invC%*%(Y-Un%*%mut)


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
# # Matrice de variance-covariance de (Z,Y)
# RC <- kronecker(R,Cnoeud, FUN = "*")
# 
# # Mettre sous forme vec
# vecY <- matrix(Y, ncol=1, byrow = F)
# mY <- matrix(Un%*%mut, ncol=1, byrow =F)
# Um <- matrix(rep(1,102),ncol=1)
# mZ <- matrix(Um%*%mut, ncol=1, byrow=F)
# 
# # Extraction des matrices SigmaXY et SigmaY
# SigmaZY <- rbind(cbind(RC[105:206,1:104],RC[105:206,207:310]),cbind(RC[311:412,1:104],RC[311:412,207:310]))
# SigmaY <- rbind(cbind(RC[1:104,1:104],RC[1:104,207:310]),cbind(RC[207:310,1:104],RC[207:310,207:310]))
# 
# # Estimateur
# Estimateur <- mZ + SigmaZY%*%solve(SigmaY)%*%(vecY-mY)


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

# Mettre sous forme vec
vecY <- matrix(Y, ncol=1, byrow = F)
mY <- matrix(Un%*%mut, ncol=1, byrow =F)
Um <- matrix(rep(1,102),ncol=1)
mZ <- matrix(Um%*%mut, ncol=1, byrow=F)

# Estimateur
Estimateur <- mZ + kronecker(diag(2),(Cmn%*%solve(C)), FUN = "*")%*%(vecY-mY)

################################################################################
# Calcul de l'AIC ?

L_hat <- -103 * log(2*pi) - 103/2 * log(det(R)) - 1/2*sum(diag(solve(R)%*%t(Y-Un%*%mut)%*%(Y-Un%*%mut)))
AIC <- -2*L_hat + 2*6
# 4405


