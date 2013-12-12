library(GGMselect)
data2 = read.table("GaussianData.txt", header=TRUE)
data=as.matrix(data2)
set.seed(4)
#GRest <- selectFast(data, family="C01")# Matrix singular.
 print(GRest <- selectFast(data, family="C01", min.ev=0.01)) #OK

print(GRest <- selectFast(data, family="LA")) # OK, mais voir le"VOIR SH" dans loop.c
#selectQE(data) # Matrix singular.
print(GRest <-selectQE(data, min.ev=0.01)) # OK
#GRest <- selectFast(data, family="EW") # all results are Na
#GRest <- selectFast(data, family="EW", beta=400) # Matrix singular.
# GRest <- selectFast(data, family="EW", beta=4000, tau=0.1) # Matrix singular.
#GRest <- selectFast(data, family="EW", beta=4000, tau=0.1, h=0.0001) # Matrix singular.

print(GRest <- selectFast(data, family="EW", beta=4000,min.ev=0.01)) # OK


