#' Compute the per-individual proportion homozygous
#'
#' Compute the per-individual proportion homozygous.
#' For example, to use as a predictor for a directional dominance model.
#'
#' @param M dosage matrix. Assumes SNPs in M coded 0, 1, 2 (requires rounding dosages to integers). M is Nind x Mrow, numeric matrix, with row/columanes to indicate SNP/ind ID
#'
#' @return vector of per-individual proportion homozygous
#' @export
#' @family helper
getPropHom<-function(M){
  W<-M; W[which(W==2)]<-0;
  # f = 1 − h/N,
  # where N is the number of SNPs
  f<-1-(rowSums(W)/ncol(W))
  return(f)
}

#' Centers dosage matrix
#'
#' Centers dosage matrix, e.g. for use in whole-genome regressions like rrBLUP
#'
#' @param M dosage matrix. Assumes SNPs in M coded 0, 1, 2 (requires rounding dosages to integers). M is Nind x Mrow, numeric matrix, with row/columanes to indicate SNP/ind ID
#'
#'
#' @return
#' @export
#' @family helper
centerDosage<-function(M){
  freq <- colMeans(M,na.rm=T)/2
  P <- matrix(rep(freq,nrow(M)),byrow=T,ncol=ncol(M))
  Z <- M-2*P
  return(Z)
}

#' Converts a dosage matrix into a matrix of centered "classical"-coded dominance deviations.
#'
#' Converts a dosage matrix into a matrix of centered dominance deviations.
#' This sets up the "classical" (aka "Statistical") partition of additive dominance in terms of breeding values and dom. deviations.
#' See function \code{\link{dose2domDevGenotypic}} to set-up the "genotypic" (aka "biological") partition in terms of genotypic effects.
#' See Vitezica et al. 2013.
#' Also see Varona et al. 2018. https://doi.org/10.3389/fgene.2018.00078
#'
#' @param M dosage matrix. Assumes SNPs in M coded 0, 1, 2 (requires rounding dosages to integers). M is Nind x Mrow, numeric matrix, with row/columanes to indicate SNP/ind ID
#'
#' @return
#' @export
#' @family helper
dose2domDev<-function(M){
  freq <- colMeans(M,na.rm=T)/2
  P <- matrix(rep(freq,nrow(M)),byrow=T,ncol=ncol(M))
  W<-M;
  W[which(W==1)]<-2*P[which(W==1)];
  W[which(W==2)]<-(4*P[which(W==2)]-2);
  W <- W-2*(P^2)
  return(W)
}


#' Converts a dosage matrix into a matrix of centered "genotypic"-coded dominance deviations.
#'
#' This sets up the "genotypic" (aka "biological") partition of additive dominance in terms of their genotypic effects instead of in terms of breeding values or dominance deviations.
#' See function \code{\link{dose2domDev}} to set-up the "statistical" (aka "classical") partition in terms of genotypic effects.
#' See Vitezica et al. 2013.
#' Also see Varona et al. 2018. https://doi.org/10.3389/fgene.2018.00078
#'
#' @param M dosage matrix. Assumes SNPs in M coded 0, 1, 2 (requires rounding dosages to integers). M is Nind x Mrow, numeric matrix, with row/columanes to indicate SNP/ind ID
#'
#' @return
#' @export
#' @family helper
dose2domDevGenotypic<-function(M){
  freq <- colMeans(M,na.rm=T)/2
  P <- matrix(rep(freq,nrow(M)),byrow=T,ncol=ncol(M))
  W<-M; W[which(W==2)]<-0;
  W <- W-(2*P*(1-P))
  return(W)
}

#' Compute allele frequencies
#'
#' get a vector of allele frequences from a dosage matrix
#'
#' @param M dosage matrix. Assumes SNPs in M coded 0, 1, 2 (requires rounding dosages to integers). M is Nind x Mrow, numeric matrix, with row/columanes to indicate SNP/ind ID
#'
#' @return vector of allele frequencies, names = SNP IDs if in cols of M
#' @export
#' @family helper
getAF<-function(M){ colMeans(M,na.rm=T)/2 }

#' Compute minor allele frequencies
#'
#' get a vector of \emph{minor} allele frequences from a dosage matrix
#'
#' @param M dosage matrix. Assumes SNPs in M coded 0, 1, 2 (requires rounding dosages to integers). M is Nind x Mrow, numeric matrix, with row/columanes to indicate SNP/ind ID
#'
#' @return vector of \emph{minor} allele frequencies, names = SNP IDs if in cols of M
#' @export
#' @family helper
getMAF<-function(M){
  freq<-colMeans(M, na.rm=T)/2; maf<-freq;
  maf[which(maf > 0.5)]<-1-maf[which(maf > 0.5)]
  return(maf) }

#' Filter a dosage matrix by a minor allele frequency
#'
#' filter a dosage matrix by minor allele frequence.
#' get a vector of allele frequences from a dosage matrix
#'
#' @param M dosage matrix. Assumes SNPs in M coded 0, 1, 2 (requires rounding dosages to integers). M is Nind x Mrow, numeric matrix, with row/columanes to indicate SNP/ind ID
#' @param thresh threshold value. Columns of M with maf<thresh will be removed
#' @return dosage matrix potentially with columns removed
#' @export
#' @family helper
maf_filter<-function(M,thresh){
  freq<-colMeans(M, na.rm=T)/2; maf<-freq;
  maf[which(maf > 0.5)]<-1-maf[which(maf > 0.5)]
  snps1<-M[,which(maf>thresh)];
  return(snps1) }

#' Remove invariant SNPs from dosage matrix
#'
#' filter a dosage matrix, removing invariant markers. Removes e.g. cases where MAF=0.5 but all dosages == 1 (het).
#'
#' @param M dosage matrix. Assumes SNPs in M coded 0, 1, 2 (requires rounding dosages to integers). M is Nind x Mrow, numeric matrix, with row/columanes to indicate SNP/ind ID
#' @param thresh threshold value. Columns of M with maf<thresh will be removed
#' @return dosage matrix potentially with columns removed
#' @export
#' @family helper
remove_invariant<-function(M){
  snps1<-M[ ,apply(M, 2, var) != 0]
  return(snps1) }

#' Compute ridge-regression equivalent marker effects from GBLUP results ("Backsolve SNP effects")
#'
#' From the GBLUP solutions and a centered SNP matrix backsolve SNP effects
#'
#' @param Z Centered marker matrix (dominance deviations must also be centered)
#' @param g The solutions (blups, i.e. GEBVs) from the GBLUP model
#'
#' @return matrix of SNP effects matching RR-BLUP / SNP-BLUP
#' @export
#' @family helper
backsolveSNPeff<-function(Z,g){
  # New version of this function attempts to be
  # robust to singularities that sometimes arisen
  # should return "NA" if all else fails

  # debug: # Z=M; g=ga;
  ZZt<-tcrossprod(Z);
  # setting tol=rcond(ZZt) below seems to avoid compute-singularities
  ## better than adding a small value to the diag of ZZt
  bslv<-function(Z,ZZt,g){
    return(crossprod(Z,solve(ZZt,tol = rcond(ZZt)))%*%g) }
  possibly_bslv<-possibly(bslv,NA_real_)
  bslv_out<-possibly_bslv(Z,ZZt,g)

  # last ditch attempt
  if(!"matrix" %in% class(bslv_out)){
    if(is.na(bslv_out)){
      # try adding small diag to ZZt
      diag(ZZt)<-diag(ZZt)+1e-8
      bslv_out<-possibly_bslv(Z,ZZt,g) } }
  # last LAST ditch attempt
  if(!"matrix" %in% class(bslv_out)){
    if(is.na(bslv_out)){
      # try adding slightly less small diag to ZZt
      diag(ZZt)<-diag(ZZt)+1e-6
      bslv_out<-possibly_bslv(Z,ZZt,g) } }
  return(bslv_out)
}

#' Compute the LD matrix from dosages
#'
#' Compute the \eqn{p_{SNP} \times p_{SNP}} variance-covariance matrix of SNP dosages.
#' This is an estimator of the LD between loci within a given population.
#'
#' @param Z column-centered matrix of SNP dosages. Assumes SNPs in Z were originally coded 0, 1, 2 were column centered.
#'
#' @return
#'  NOTE: this matrix is going to be big in practice.
#'  The \eqn{p_{SNP} \times p_{SNP}} variance-covariance matrix of SNP dosages.
#'  may be worth computing in an R session using multi-threaded BLAS
#' @export
#' @family helper
genoVarCovarMatFunc<-function(Z){
  SigmaM<-1/nrow(Z)*t(Z)%*%Z
  return(SigmaM)
}


#' Create a matrix of pairwise recombination frequencies from a genetic map
#'
#' Compute the pairwise recombination frequencies between all loci from genetic map positions.
#'
#' @param m vector of centiMorgan-scale genetic positions. names(m) correspond to a SNP_ID. Since m potentially contains all chromosomes, sets recomb. freq. b/t chrom. to 0.5
#' @param nChr number of chromosomes
#'
#' @details names(m) must be formatted as "chr"_"id" with "chr" being integer. For example: 2_QTL1 for a locus on chr. 2.
#' May be worth computing in an R session using multi-threaded BLAS.
#' @return potentially really large matrix of pairwise recombination frequencies between loci
#' @export
#' @family helper
genmap2recombfreq<-function(m,nChr){
  d<-as.matrix(dist(m,upper=T,diag = T,method='manhattan'))
  c1<-0.5*(1-exp(-2*(d/100)))
  # Since m contains all chromosomes, set recomb. freq. b/t chrom. to 0.5
  for(i in 1:nChr){
    c1[grepl(paste0("^",i,"_"),rownames(c1)),!grepl(paste0("^",i,"_"),colnames(c1))]<-0.5
    c1[!grepl(paste0("^",i,"_"),rownames(c1)),grepl(paste0("^",i,"_"),colnames(c1))]<-0.5
  }
  return(c1)
}

#' Standardized selection intensity
#'
#' Compute the standardized selection intensity from the proportion selection.
#'
#' @param propSel proportion selection
#'
#' @return
#' @export
#' @family helper
intensity<-function(propSel){ dnorm(qnorm(1-propSel))/propSel } # standardized selection intensity

#' Construct several types of kinship matrix
#'
#' Function to create additive and dominance genomic relationship matrices from biallelic dosages.
#'
#'
#' @param M dosage matrix. Assumes SNPs in M coded 0, 1, 2 (requires rounding dosages to integers). M is Nind x Mrow, numeric matrix, with row/columanes to indicate SNP/ind ID.
#' @param type string, "add" or "domGenotypic" and "domClassic".
#' \enumerate{
#'    \item \code{type="add"} gives same as \code{\link[rrBLUP]{A.mat}}, i.e. Van Raden, Method 1.
#'    \item \code{type="domClassic"} and \code{type="domGenotypic"} give the classical and genotypic parameterization according to Vitezica et al. 2013. Genetics. Both dominance matrices, in combo with the additive matrix, predict total merit identically. Difference in partition of add-dom, meaning of variance components and genomic predictions.
#' @return square symmetric genomic relationship matrix
#' @export
#' @family helper
kinship<-function(M,type){
  M<-round(M)
  freq <- colMeans(M,na.rm=T)/2
  P <- matrix(rep(freq,nrow(M)),byrow=T,ncol=ncol(M))
  if(type=="add"){
    Z <- M-2*P
    varD<-sum(2*freq*(1-freq))
    K <- tcrossprod(Z)/ varD
    return(K)
  }
  if(type=="domGenotypic"){
    W<-M; W[which(W==2)]<-0;
    W <- W-(2*P*(1-P))
    varD<-sum((2*freq*(1-freq))*(1-(2*freq*(1-freq))))
    D <- tcrossprod(W) / varD
    return(D)
  }
  if(type=="domClassic"){
    W<-M;
    W[which(W==1)]<-2*P[which(W==1)];
    W[which(W==2)]<-(4*P[which(W==2)]-2);
    W <- W-2*(P^2)
    varD<-sum((2*freq*(1-freq))^2)
    D <- tcrossprod(W) / varD
    return(D)
  }
}

#' Make a data.frame of pairwise crosses between a set parents
#'
#' Make a data.frame of all pairwise matings given a vector of parent IDs.
#' Include selfs. No reciprocal crosses, i.e. use as male == use as female.
#' Diagonal and upper-triangle of mating matrix.
#'
#' @param parents
#'
#' @return tibble, two columns, sireID and damID, all pairwise crosses (see details).
#' @export
#' @family helper
crosses2predict<-function(parents){
  CrossesToPredict<-matrix(NA,nrow=length(parents),ncol=length(parents))
  CrossesToPredict[upper.tri(CrossesToPredict,diag = T)]<-1
  rownames(CrossesToPredict)<-colnames(CrossesToPredict)<-parents
  CrossesToPredict<-CrossesToPredict %>%
    as.data.frame %>%
    tibble::rownames_to_column(var = "sireID") %>%
    tidyr::pivot_longer(cols = (-sireID), names_to = "damID", values_to = "keep") %>%
    dplyr::filter(keep==1) %>%
    dplyr::select(-keep)
  return(CrossesToPredict)
}

#' Quadratic matrix form
#'
#' does \eqn{x \%\ast\% D \%\ast\% t(x)} a bit faster
#'
#' @param D square symmetric matrix
#' @param x vector 1
#' @param y vector 2
#'
#' @return square symmetrix matrix
#' @export
#' @family helper
quadform<-function(D,x,y){ return(as.numeric(colSums(x*(D%*%y)))) }


#' Converts an array of posterior samples of multi-trait marker effects to a named list (one for each trait).
#'
#' Converts an array of posterior samples of multi-trait marker effects to a named list (one for each trait).
#'
#' @param effectsArray According to BGLR documentation: 3D array, with dim=c(nRow,p,traits), where nRow number of MCMC samples saved, p is the number of predictors and traits is the number of traits. \code{\link[BGLR]{Multitrait}}writes a binary file to disk when saveEffects=TRUE is specified. It can be read to R with \code{\link[BGLR]{readBinMatMultitrait}}.
#' @param snpIDs character vector with labels for the predictors (SNPs), numeric should work too, but untested.
#' @param traits character vector to label the traits.
#' @param nIter number of iterations used for MCMC; used internally only to exclude burn-in samples from computation
#' @param burnIn burnIn for MCMC; used internally only to exclude burn-in samples from computation
#' @param thin thin for MCMC; used internally only to exclude burn-in samples from computation
#' @family helper
#'
#' @return list of matrices, one matrix per trait, each matrix has \code{nrow((nIter-burnIn)/thin)} and \code{ncol(length(snpIDs))}. Each element of the list is named with a string identifying the trait and the colnames of each matrix are labelled with snpIDs.
#' @export
#' @family helper
effectsArray2list<-function(effectsArray, snpIDs, traits, nIter, burnIn, thin){
  # Discard burnIn
  effectsArray<-effectsArray[-c(1:burnIn/thin),,]
  # Add dimnames
  dimnames(effectsArray)[[2]]<-snpIDs
  dimnames(effectsArray)[[3]]<-traits
  # 3D arrays of effects to lists-of-matrices (by trait)
  effectsArray<-purrr::array_branch(effectsArray,3)
  return(effectsArray)
}
