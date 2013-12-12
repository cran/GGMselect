library(GGMselect)
attach("G3simone.Rdata")
attach("genes.lm.Rdata")

print(selectMyFam(genes.lm, MyFamily = G3simone, K=10))
