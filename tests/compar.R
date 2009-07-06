# +++++++++++++++++++++++++++++++++++++++++++++++++
# Fonctions pour comparer les resultats de selectFast par
# rapport aux originaux (ceux calculés avec le code tout en R)

# +++++++++++++++++++++++++++++++++++++++++++++++++
prmess <- function(a, nom) {
  if (!is.logical(a)) {
  cat("Error on component", nom, ":\n")
  print(a)
  return(FALSE)
  }
else
  return(TRUE)
}
# ------------------------------------------------

comparOR <- function(res, or)
  {
ret <- rep(TRUE,3)
lesnoms <- c("LA","EW", "C01", "C01.LA", "C01.LA.EW")
lesnomsor <- c("ET","EWOU", "C01", "C01.ET", "C01.ET.EWOU")
for (inom in 1:length(lesnoms)) {
  nom <- lesnoms[inom]
  nomor <- lesnomsor[inom]
  a <- all.equal(res[[nom]]$crit.min, or[[nomor]]$crit.min)
  ret[[1]] <- prmess(a, nom)
  dimnames(or[[nomor]]$Neighb) <-
   list(sub("p=","",dimnames(or[[nomor]]$Neighb)[[1]]),
   sub("p=","",dimnames(or[[nomor]]$Neighb)[[2]]),
   sub("p=","",dimnames(or[[nomor]]$Neighb)[[3]]))
  dimnames(or[[nomor]]$G) <-
   list(sub("p=","",dimnames(or[[nomor]]$G)[[1]]),
   sub("p=","",dimnames(or[[nomor]]$G)[[2]]),
   sub("p=","",dimnames(or[[nomor]]$G)[[3]]))

  a <- all.equal(res[[nom]]$Neighb, or[[nomor]]$Neighb[,,1])
  ret[[2]] <- prmess(a, nom)
  a <- all.equal(res[[nom]]$G, or[[nomor]]$G[,,1])
  ret[[3]] <- prmess(a, nom)
  } # fin inom
if (!is.logical(all.equal(ret, rep(TRUE,3))))
  return(FALSE)
else
  return(TRUE)
}
