################################################################################
# Simulation d’un mouvement brownien en deux dimensions
################################################################################
set.seed(1235)

### discrétisation du temps
pas_temps = 1/50000
### Simulation des accroissements
B_acc1 = rnorm(50000,sd=sqrt(pas_temps))
B_acc2 = rnorm(50000,sd=sqrt(pas_temps))
### Simulation d’une trajectoire
B_sim1 = c(0,cumsum(B_acc1))
B_sim2 = c(0,cumsum(B_acc2))
plot(B_sim1, B_sim2 ,type="l", xlab="B_t1", ylab="B_t2", col="darkblue")
points(0,0, col="red", pch=19, cex=1)

