# +++++++++++++++++++++++++++++++++++++++
# Supprime les labels des résulats originaux
# (ceux obtenus avec la version tout en R)
# pour pouvoir comparer avec ceux obtenus avec selectQE
# +++++++++++++++++++++++++++++++++++++++
supLab <- function (or) {
    dimnames(or) <-
   list(sub("p= ","",dimnames(or)[[1]]),
   sub("p= ","",dimnames(or)[[2]]),
   sub("p= ","",dimnames(or)[[3]]))
    return(or)
  }
