% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{dose2domDev}
\alias{dose2domDev}
\title{Converts a dosage matrix into a matrix of centered "classical"-coded dominance deviations.}
\usage{
dose2domDev(M)
}
\arguments{
\item{M}{dosage matrix. Assumes SNPs in M coded 0, 1, 2 (requires rounding dosages to integers). M is Nind x Mrow, numeric matrix, with row/columanes to indicate SNP/ind ID}
}
\value{

}
\description{
Converts a dosage matrix into a matrix of centered dominance deviations.
This sets up the "classical" (aka "Statistical") partition of additive dominance in terms of breeding values and dom. deviations.
See function \code{\link{dose2domDevGenotypic}} to set-up the "genotypic" (aka "biological") partition in terms of genotypic effects.
See Vitezica et al. 2013.
Also see Varona et al. 2018. https://doi.org/10.3389/fgene.2018.00078
}
\seealso{
Other helper: 
\code{\link{backsolveSNPeff}()},
\code{\link{centerDosage}()},
\code{\link{crosses2predict}()},
\code{\link{dose2domDevGenotypic}()},
\code{\link{effectsArray2list}()},
\code{\link{genmap2recombfreq}()},
\code{\link{genoVarCovarMatFunc}()},
\code{\link{getAF}()},
\code{\link{getMAF}()},
\code{\link{getPropHom}()},
\code{\link{intensity}()},
\code{\link{kinship}()},
\code{\link{maf_filter}()},
\code{\link{quadform}()},
\code{\link{remove_invariant}()}
}
\concept{helper}
