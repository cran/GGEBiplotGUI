\name{GGEBiplot}
\alias{GGEBiplot}
\title{
GGEBiplotGUI: Interactive GGE Biplots in R
}
\description{
The GGEBiplotGUI package
provides a graphical user interface for the construction of, interaction with, and manipulation
of GGE biplots in R. Some of the functions of the package are: (i) ranking the cultivars
based on their performance in any given environment, (ii) ranking the environments based
on the relative performance of any given cultivar, (iii) comparing the performance of any
pair of cultivars in diferent environments, (iv) identifying the best cultivar in each environment,
(v) grouping the environments based on the best cultivars, (vi) evaluating the
cultivars based on both average yield and stability and (vii) evaluating the environments
based on both discriminating ability and representativeness. Three-dimensional biplots
are incorporated via the rgl package.
}
\usage{
GGEBiplot(Data)
}
\arguments{
  \item{Data}{ A data.frame or matrix}
}
\references{
Gabriel, K. R. (1971) The biplot graphical display of matrices with application to principal component analysis. Biometrika 58, 453-467.
Galindo, M. P. (1986) Una alternativa de representacion simultanea: HJ-Biplot. Questiio, 10(1):13-23, 1986.
Yan W, Hunt LA, Sheng Q, Szlavnics Z (2000). Cultivar evaluation and mega-environment investigation based on GGE biplot. Crop Sci, 40, 597-605.
Yan W, Kang M (2003). GGE Biplot Analysis: A Graphical Tool for Breeders, Geneticists, and Agronomists. CRC Press, Boca Raton, FL, USA.
}
\author{
Elisa Frutos Bernal <efb@usal.es>
}
\examples{
data(Ontario)
GGEBiplot(Data = Ontario)
}
