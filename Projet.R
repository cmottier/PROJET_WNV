# # Donne l'ancêtre direct et le poids de l'arête associée
# 
# ancetre_direct <- function(noeud) {
#   res = NULL
#   aretes <- tree$edge
#   poids <- tree$edge.length
#   nb_aretes <- nrow(tree$edge)
#   for (i in 1:nb_aretes) {
#     if (aretes[i,2]==noeud) {return(matrix(c(aretes[i,1],poids[i])))}
#   }
#   return(res)
# }
# 
# # Donne tous les ancêtres et les poids des arêtes associées (poids fils-ancêtre affiché sous l'ancêtre)
# 
# ancetres <- function(noeud) {
#   res = NULL
#   if (is.null(ancetre_direct(noeud))) {return(res)}
#   else {
#     a <- ancetre_direct(noeud)
#     res = cbind(res, cbind(a, ancetres(a[1])))
#   }
#   return(res)
# }
# 
# # Donne le poids de l'arête commune à deux noeuds
# 
# temps_commun <- function(noeud1,noeud2) {
#   if (noeud1 == noeud2) {
#     return(sum(ancetres(noeud1)[2,]))}
#   else {
#     anc1 <- cbind(matrix(c(noeud1,0), nrow=2),ancetres(noeud1))
#     anc2 <- cbind(matrix(c(noeud2,0), nrow=2),ancetres(noeud2))
#     res = 0
#     l1 <- ncol(anc1)
#     l2 <- ncol(anc2)
#     while (anc1[1,l1-1]==anc2[1,l2-1]) {
#       res = res +anc1[2,l1]
#       l1 <- l1-1
#       l2 <- l2-1
#       if (l1==1 | l2==1) {return(res)}
#     }
#     return(res)}
# }
# 
# # Crée la matrice C
# 
# C <- matrix(nrow=104, ncol=104)
# for (i in 1:104) {
#   for (j in i:104) {
#     C[i,j] <- temps_commun(i,j)
#     C[j,i] <- temps_commun(j,i)
#   }
# }
#
# # Estimateurs de mu et R
# 
# invC <- solve(C)
# Un <- matrix(rep(1,104),ncol=1)
# Y <- as.matrix(dat[,2:3])
# 
# mut <- solve(t(Un)%*%invC%*%Un)%*%t(Un)%*%invC%*%Y
# R <- 1/103*t(Y-Un%*%mut)%*%invC%*%(Y-Un%*%mut)


# Construction de C, basée sur les feuilles (temps d'évolution commune)

C <- matrix(nrow=104, ncol=104)
for (i in 1:104) {
  C[i,i] <- node.depth.edgelength(tree)[i]
}
for (i in 1:103) {
  for (j in ((i+1):104)) {
    ancetre <- getMRCA(tree,c(i,j))
    tij <- node.depth.edgelength(tree)[ancetre]
    C[i,j] <- tij
    C[j,i] <- tij
  }
}

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

Um <- matrix(rep(1,102),ncol=1)

RC <- kronecker(R,Cnoeud, FUN = "*")
vecY <- matrix(Y, ncol=1, byrow = F)
mY <- matrix(Un%*%mut, ncol=1, byrow =F)
mX <- matrix(Um%*%mut, ncol=1, byrow=F)
SigmaXY <- rbind(cbind(RC[105:206,1:104],RC[105:206,207:310]),cbind(RC[311:412,1:104],RC[311:412,207:310]))
SigmaY <- rbind(cbind(RC[1:104,1:104],RC[1:104,207:310]),cbind(RC[207:310,1:104],RC[207:310,207:310]))
Estimateur <- mX + SigmaXY%*%solve(SigmaY)%*%(vecY-mY)



