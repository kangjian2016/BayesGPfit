library(lattice)
grids = GP.generate.grids(d=2L,num_grids=200L)
Psi_mat = GP.eigen.funcs.fast.orth(grids)
A = crossprod(Psi_mat)
range(A[upper.tri(A,diag=FALSE)])
image(A)

fig = list()
a = 62
for(i in a+1:4){
  fig[[i-a]] = levelplot(Psi_mat[,i]~grids[,1]+grids[,2])
}
plot(fig[[1]],split=c(1,1,2,2),more=TRUE)
plot(fig[[2]],split=c(1,2,2,2),more=TRUE)
plot(fig[[3]],split=c(2,1,2,2),more=TRUE)
plot(fig[[4]],split=c(2,2,2,2))

Psi_mat = GP.eigen.funcs.fast(grids)
fig = list()
for(i in a+1:4){
  fig[[i-a]] = levelplot(Psi_mat[,i]~grids[,1]+grids[,2])
}
plot(fig[[1]],split=c(1,1,2,2),more=TRUE)
plot(fig[[2]],split=c(1,2,2,2),more=TRUE)
plot(fig[[3]],split=c(2,1,2,2),more=TRUE)
plot(fig[[4]],split=c(2,2,2,2))
