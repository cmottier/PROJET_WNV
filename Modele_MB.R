################################################################################
# Modèle brownien simple
# Y ~ MN(Un*mut, C, R)
# Y = Un*mut+E où E ~ MN(0, C,R)
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

# Construction de Cnoeud, basée sur tous les noeuds

f <- function(x) {
  if (x<105) return(x)
  else return(x+1)
}

Cnoeud <- matrix(nrow=206, ncol=206) # On enlève la racine 105
for (i in 1:206) {
  Cnoeud[i,i] <- node.depth.edgelength(tree)[f(i)]
  }
for (i in 1:205) {
  for (j in ((i+1):206)) {
    ancetre <- getMRCA(tree,c(f(i),f(j)))
    tij <- node.depth.edgelength(tree)[ancetre]
    Cnoeud[i,j] <- tij
    Cnoeud[j,i] <- tij
  }
}

################################################################################
# Estimateur des valeurs des ancêtres avec E[X|Y]=mX + SigmaXY*SigmaY^(-1)*(Y-mY)

# Matrice de variance-covariance de (X,Y)
RC <- kronecker(R,Cnoeud, FUN = "*")

# Mettre sous forme vec
vecY <- matrix(Y, ncol=1, byrow = F)
mY <- matrix(Un%*%mut, ncol=1, byrow =F)
Um <- matrix(rep(1,102),ncol=1)
mX <- matrix(Um%*%mut, ncol=1, byrow=F)

# Extraction des matrices SigmaXY et SigmaY
SigmaXY <- rbind(cbind(RC[105:206,1:104],RC[105:206,207:310]),cbind(RC[311:412,1:104],RC[311:412,207:310]))
SigmaY <- rbind(cbind(RC[1:104,1:104],RC[1:104,207:310]),cbind(RC[207:310,1:104],RC[207:310,207:310]))

# Estimateur
Estimateur <- mX + SigmaXY%*%solve(SigmaY)%*%(vecY-mY)



