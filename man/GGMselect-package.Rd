\name{GGMselect-package}
\alias{GGMselect-package}
\alias{GGMselect}
\docType{package}
\title{Gaussian Graphs Models selection}
\description{GGMselect is a package dedicated to graph estimation in Gaussian
Graphical Models. The main functions return the adjacency matrix of an
undirected graph estimated from a data matrix.

This package is developed in the Applied Mathematics and Informatics
(\url{http://www.jouy.inra.fr/mia_eng/}) Lab of INRA - Jouy-en-Josas, France.

To cite GGMselect, please use \code{citation("GGMselect")}.}
\details{\tabular{ll}{
Package: \tab GGMselect\cr
Version: \tab 0.1-0\cr
Date: \tab 2009-07-06\cr
Depends: \tab mvtnorm, lars, mvbutils (>= 2.0), gtools\cr
Suggests: \tab network, glasso\cr
License:\tab GPL (>= 3)\cr
URL: \tab \url{http://www.inra.fr/miaj/public/logiciels/GGMselect}\cr
}}
\author{Annie Bouvier, Christophe Giraud, Sylvie Huet, Verzelen Nicolas

  Maintainer: Annie Bouvier <Annie.Bouvier@jouy.inra.fr>
  
  More details are available on \url{../doc/Notice.pdf}}
\references{Please use \code{citation("GGMselect")}.}
\seealso{ \code{\link{selectFast}}, \code{\link{selectQE}},
 \code{\link{selectMyFam}},\code{\link{convertGraph}},
  \code{\link{simulateGraph}}, \code{\link{penalty}}}
\examples{
p=30
n=30
# simulate graph
eta=0.11
Gr <- simulateGraph(p,eta)
# simulate data
X <- rmvnorm(n, mean=rep(0,p), sigma=Gr$C)
# estimate graph
GRest <- selectFast(X)

# plot result
library(network)
par(mfrow=c(1,2))
gV <- network(Gr$G)
plot(gV,jitter=TRUE, usearrows = FALSE, label=1:p,displaylabels=TRUE)
g <- network(GRest$EW$G)
plot(g, jitter=TRUE, usearrows = FALSE,  label=1:p,displaylabels=TRUE)
}
\keyword{ package }
\keyword{ multivariate }
\keyword{ graphs }
\keyword{ models }
