% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predCrossVar.R
\name{predCrossMeans}
\alias{predCrossMeans}
\title{Predict cross means}
\usage{
predCrossMeans(
  CrossesToPredict,
  predType,
  AddEffectList,
  DomEffectList = NULL,
  doseMat,
  ncores = 1,
  ...
)
}
\arguments{
\item{CrossesToPredict}{data.frame or tibble, col/colnames: sireID, damID.
sireID and damID must both be in the haploMat.}

\item{predType}{string, "BV" or "TGV". "BV" predicts cross mean breeding
values as the mean GEBV of parents. "TGV" predicts the cross total genetic
value. Warning: prediction of meanTGV with F-M Eqn. 14.6 appropriate only
using a+d partition not allele sub. + dom. dev.; genotypic NOT classical in
terms used by Vitezica et al. 2013. For that reason,
\code{\link{predCrossMeans}} has a "predType" not a "modelType" argument
predType="TGV" uses Falconer-MacKay Eqn. 14.6 and takes add and dom
effects. predType="BV" input should be allele subst. effs, computes
mid-parent GEBV there is no equivalent to predicting the dominance variance
for the mean thus the difference from the predCrossVars() function. NOTICE:
NOT SAME as predType argument used in \code{\link{predCrossVars}}, sorry.}

\item{AddEffectList}{list of ADDITIVE effect matrices, one matrix per trait,
Each element of the list is named with a string identifying the trait and
the colnames of each matrix are labelled with snpIDs.}

\item{DomEffectList}{list of DOMINANCE effect matrices, one matrix per trait,
Each element of the list is named with a string identifying the trait and
the colnames of each matrix are labelled with snpIDs.}

\item{doseMat}{dosage matrix. required only for modelType=="DirDom". Assumes
SNPs coded 0, 1, 2. Nind rows x Nsnp cols, numeric matrix, with rownames
and colnames to indicate SNP/ind ID}

\item{ncores}{number of cores, parallelizes across \code{CrossesToPredict},
in multi-trait cases, process traits for each family in serial within each
worker.}

\item{...}{}
}
\value{
tibble, each row contains predictions for a single cross. Columns:
\itemize{
 \item \code{"Trait"}:
 \item \code{"sireID"}:
 \item \code{"damID"}:
 \item \code{"sireGEBV"}: genomic estimated breeding value (GEBV) of the male parent of the cross
 \item \code{"damGEBV"}: genomic estimated breeding value (GEBV) of the female parent of the cross
 \item \code{"predOf"}: "MeanBV" or "MeanTGV"
 \item \code{"predMean"}: The predicted mean value for the cross
}
}
\description{
Function to predict the mean performances of the offspring of crosses. Takes
a list of crosses to predict, marker effects, parental allele dosage matrix
as input. Predicts potentially over multiple crosses and multiple traits.
With \code{predType="BV"} predicts the mid-parent of crosses by computing
parental GEBV. With \code{predType="TGV"} predicts the mean total merit of
cross offspring using a Falconer-MacKay Eqn. 14.6 and takes user supplied
additive and dominance effects as input. The additive-dominance effects
should be partitioned according to the "genotypic" marker codings (see
Vitezica et al. 2013. GENETICS).
}
\seealso{
Other predCrossVar: 
\code{\link{calcCrossLD}()},
\code{\link{calcGameticLD}()},
\code{\link{predCrossVars}()}
}
\concept{predCrossVar}
