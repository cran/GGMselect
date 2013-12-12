library(GGMselect)
itest=5
# +++++++++++++++++++++++++++++++++++++++++++++++++
p=30
n=30
eta=0.17
dmax=4
iG = 7
iS = 9
set.seed(iG)
Gr <- simulateGraph(p,eta)
set.seed(iS*(pi/3.1415)**iG)
X <- rmvnorm(n, mean=rep(0,p), sigma=Gr$C)

K=2.5

# ptm <- proc.time()

mon5GRest <- selectQE(X, dmax, K, verbose=TRUE)

# cat("Elapsed time ",  proc.time()-ptm,"\n")
cat ("End of test ", itest, "\n")


