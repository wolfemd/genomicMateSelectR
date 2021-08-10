#' Predict cross variances
#'
#' Function to predict the variances (and trait-trait co-variances) expected in the offspring of a cross. Takes a list of crosses to predict, marker effects, parental haplotype matrix and recombination frequency matrix as input. Predicts potentially over multiple crosses and multiple traits. When multiple traits are supplied, trait-trait co-variances are also predicted.
#'
#' @param CrossesToPredict data.frame or tibble, col/colnames: sireID, damID. sireID and damID must both be in the haploMat.
#' @param modelType string, "A" or "AD", should additive only or additive and dominance variances be predicted?
#' @param AddEffectList list of ADDITIVE effect matrices, one matrix per trait, Each element of the list is named with a string identifying the trait and the colnames of each matrix are labelled with snpIDs.
#' @param DomEffectList list of DOMINANCE effect matrices, one matrix per trait, Each element of the list is named with a string identifying the trait and the colnames of each matrix are labelled with snpIDs.
#' @param predType, string, default \code{predType="VPM"}, "VPM" or "PMV". Choose option "VPM" if you have REML marker effect estimates (or posterior-means from MCMC) one set of marker effect estimates per trait. Variance of posterior means is faster but the alternative predType=="PMV" is expected to be less biassed. PMV requires user to supply a (probably LARGE) variance-covariance matrix of effects estimates.
#' @param haploMat matrix of phased haplotypes, 2 rows per sample, cols = loci, {0,1}, rownames assumed to contain GIDs with a suffix, separated by "_" to distinguish haplotypes
#' @param recombFreqMat a square symmetric matrix with values = (1-2*c1), where c1=matrix of expected recomb. frequencies. The choice to do 1-2c1 outside the function was made for computation efficiency; every operation on a big matrix takes time.
#' @param ncores number of cores, parallelizes across \code{CrossesToPredict}, in multi-trait cases, process traits for each family in serial within each worker.
#' @param nBLASthreads number of cores for each worker to use for multi-thread BLAS
#' @param ...
#'
#' @return tibble, each row contains predictions for a single cross. Columns:
#' \itemize{
#'  \item \code{"Nsegsnps"}: give the number of segregating sites in the cross and thus those used for variance predictions
#'  \item \code{"ComputeTime"}: time in minutes
#'  \item \code{"predVars"}: list-column, each element is a tibble, containing all predicted variances and covariances in a long-format.
#' }
#' @export
#' @family predCrossVar
predCrossVars<-function(CrossesToPredict,modelType,
                       AddEffectList,DomEffectList=NULL,
                       predType="VPM",
                       haploMat,recombFreqMat,
                       ncores=1,nBLASthreads=NULL,...){
  starttime<-proc.time()[3]
  # Center posterior distribution of effects
  ## on posterior mean across MCMC samples
  AddEffectList<-purrr::map(AddEffectList,~scale(.,center = T, scale = F));
  if(modelType=="AD"){
    DomEffectList<-purrr::map(DomEffectList,~scale(.,center = T, scale = F)) }

  ## Extract the posterior mean effects vectors
  ## If predType="VPM" this just recovers the original effects
  postMeanAddEffects<-purrr::map(AddEffectList,~attr(.,which = "scaled:center"))
  if(modelType=="AD"){
    postMeanDomEffects<-purrr::map(DomEffectList,~attr(.,which = "scaled:center"))
  } else {
    postMeanDomEffects<-NULL
    }
  parents<-union(CrossesToPredict$sireID,
                 CrossesToPredict$damID)
  haploMat<-haploMat[sort(c(paste0(parents,"_HapA"),paste0(parents,"_HapB"))),]

  if(predType=="VPM"){
    AddEffectList<-NULL;
    if(predType=="VPM" & modelType=="AD"){ DomEffectList<-NULL; } }

  # Set-up a loop over the crosses
  require(furrr); plan(multisession, workers = ncores)
  options(future.globals.maxSize=+Inf); options(future.rng.onMisuse="ignore")
  crossespredicted<-CrossesToPredict %>%
    dplyr::mutate(predVars=future_pmap(.,
                                predOneCross,
                                modelType=modelType,
                                haploMat=haploMat,
                                recombFreqMat=recombFreqMat,
                                predType=predType,
                                postMeanAddEffects=postMeanAddEffects,
                                postMeanDomEffects=postMeanDomEffects,
                                AddEffectList=AddEffectList,
                                DomEffectList=DomEffectList,
                                nBLASthreads=nBLASthreads))
  plan(sequential)

  crossespredicted %<>%
    unnest(predVars)

  totcomputetime<-proc.time()[3]-starttime
  print(paste0("Done predicting fam vars. ",
               "Took ",round((totcomputetime)/60,2),
               " mins for ",nrow(crossespredicted)," crosses"))
  return(crossespredicted)
}
# INTERNAL FUNCTION - predict all variance-covariance components for one cross
#' Title
#'
#' @param sireID
#' @param damID
#' @param modelType
#' @param haploMat
#' @param recombFreqMat
#' @param predType
#' @param postMeanAddEffects
#' @param postMeanDomEffects
#' @param AddEffectList
#' @param DomEffectList
#' @param nBLASthreads
#' @param ...
#'
#' @return
predOneCross<-function(sireID,damID,modelType,
                       haploMat,recombFreqMat,
                       predType,
                       postMeanAddEffects,
                       postMeanDomEffects=NULL,
                       AddEffectList=NULL,DomEffectList=NULL,
                       nBLASthreads,...){
  starttime<-proc.time()[3]

  if(!is.null(nBLASthreads)) { RhpcBLASctl::blas_set_num_threads(nBLASthreads) }

  # Before predicting variances
  # check for and remove SNPs that
  # won't segregate, i.e. are fixed in parents
  ### hopes to save time / mem
  x<-colSums(rbind(haploMat[grep(paste0("^",sireID,"_"),rownames(haploMat)),],
                   haploMat[grep(paste0("^",damID,"_"),rownames(haploMat)),]))
  segsnps2keep<-names(x[x>0 & x<4])

  if(length(segsnps2keep)>0){
    recombFreqMat<-recombFreqMat[segsnps2keep,segsnps2keep,drop=F]
    haploMat<-haploMat[,segsnps2keep,drop=F]
    postMeanAddEffects<-purrr::map(postMeanAddEffects,~.[segsnps2keep])
    if(modelType=="AD"){
      postMeanDomEffects<-purrr::map(postMeanDomEffects,~.[segsnps2keep]) }

    # calc cross LD matrix
    progenyLD<-calcCrossLD(sireID,damID,recombFreqMat,haploMat)
    rm(recombFreqMat,haploMat); gc()

    # Set-up loop over variance and covarance parameters
    ## Trait variances to-be-predicted
    traits<-names(postMeanAddEffects)
    varcovars<-tibble::tibble(Trait1=traits,
                              Trait2=traits)
    ## If multiple traits
    if(length(traits)>1){
      ## Add covariances to-be-predicted
      varcovars<-dplyr::bind_rows(varcovars, # trait variances
                                  combn(traits,2,simplify = T) %>% # covariances
                                    t(.) %>% #
                                    `colnames<-`(.,c("Trait1","Trait2")) %>%
                                    tibble::as_tibble(.)) }
    varcovars<-varcovars %>%
      dplyr::mutate(predVars=purrr::pmap(.,
                                         predOneCrossVar,
                                         modelType=modelType,
                                         progenyLD=progenyLD,
                                         # haploMat=haploMat,
                                         # recombFreqMat=recombFreqMat,
                                         predType=predType,
                                         postMeanAddEffects=postMeanAddEffects,
                                         postMeanDomEffects=postMeanDomEffects,
                                         AddEffectList=AddEffectList,
                                         DomEffectList=DomEffectList)) %>%
      unnest(predVars)
    computetime<-proc.time()[3]-starttime
    out_thiscross<-tibble::tibble(Nsegsnps=length(segsnps2keep),
                          ComputeTime=computetime,
                          predVars=list(varcovars))
  } else {
    computetime<-proc.time()[3]-starttime
    out_thiscross<-tibble::tibble(Nsegsnps=length(segsnps2keep),
                          ComputeTime=computetime,
                          predVars=list()) }
  return(out_thiscross)
}
# INTERNAL FUNCTION - predict one variance-covariance component for one cross
#' Title
#'
#' @param Trait1
#' @param Trait2
#' @param progenyLD
#' @param modelType
#' @param predType
#' @param postMeanAddEffects
#' @param postMeanDomEffects
#' @param AddEffectList
#' @param DomEffectList
#' @param segsnps2keep
#' @param ...
#'
#' @return
predOneCrossVar<-function(Trait1,Trait2,progenyLD,modelType,
                          #haploMat,recombFreqMat,
                          predType,
                          postMeanAddEffects,
                          postMeanDomEffects=NULL,
                          AddEffectList=NULL,
                          DomEffectList=NULL,
                          segsnps2keep=NULL,...){

  if(predType=="PMV"){
    # Posterior Sample Variance-Covariance Matrix of Marker Effects
    postVarCovarOfAddEffects<-(1/(nrow(AddEffectList[[Trait1]])-1))*crossprod(AddEffectList[[Trait1]],AddEffectList[[Trait2]])
    postVarCovarOfAddEffects<-postVarCovarOfAddEffects[segsnps2keep,
                                                       segsnps2keep,drop=F]
    if(modelType=="AD"){
      postVarCovarOfDomEffects<-(1/(nrow(DomEffectList[[Trait1]])-1))*crossprod(DomEffectList[[Trait1]],DomEffectList[[Trait2]])
      postVarCovarOfDomEffects<-postVarCovarOfDomEffects[segsnps2keep,
                                                         segsnps2keep,drop=F]
      }
  } else {
    postVarCovarOfAddEffects<-NULL;
    postVarCovarOfDomEffects<-NULL;
  }

  ## Predicts variance for one cross and
  ### one variance or covariance paramater (Trait1-Trait2)

  ## Predict cross additive variance
  #### posterior mean (VPM)
  predVarA_vpm<-genomicMateSelectR::quadform(D=progenyLD,
                                             x=postMeanAddEffects[[Trait1]],
                                             y=postMeanAddEffects[[Trait2]])

  #### posterior mean (co)variance (PMV)
  if(predType=="PMV"){ predVarA_pmv<-predVarA_vpm+sum(diag(progenyLD%*%postVarCovarOfAddEffects)) }

  if(modelType=="AD"){
    ## Predict cross dominance variance
    #### VPM
    progenyLDsq<-progenyLD*progenyLD
    predVarD_vpm<-genomicMateSelectR::quadform(D=progenyLDsq,
                                               x=postMeanDomEffects[[Trait1]],
                                               y=postMeanDomEffects[[Trait2]])
    #### PMV
    if(predType=="PMV"){ predVarD_pmv<-predVarD_vpm+sum(diag(progenyLDsq%*%postVarCovarOfDomEffects)) }
  }

  if(modelType=="A"){ rm(progenyLD); gc() }
  if(modelType=="AD"){ rm(progenyLD,progenyLDsq); gc() }

  # Tidy the results
  out<-tibble::tibble(predOf="VarA",predVar=predVarA_vpm)
  ### VarA
  if(predType=="PMV"){
    out<-out %>%
      rename(VPM=predVar) %>%
      mutate(PMV=predVarA_pmv) }
  ### If modelType=="AD" - VarD
  if(modelType=="AD"){
    if(predType=="VPM"){
      out %<>% dplyr::bind_rows(tibble::tibble(predOf="VarD",predVar=predVarD_vpm)) }
    if(predType=="PMV"){
      out %<>% dplyr::bind_rows(tibble::tibble(predOf="VarD",
                                 VPM=predVarD_vpm,
                                 PMV=predVarD_pmv)) }
  }
  return(out)
}


#' Calculate the gametic LD matrix for a single parent
#'
#' Function to compute gametic LD matrix for a single parent. Uses a matrix of recombination frequencies and the supplied haplotypes of the parent as in Bijma et al. 2020 and Lehermeier et al. 2017b (and others).
#'
#' @param parentGID string, the GID of an individual. Needs to correspond to renames in haploMat
#' @param recombFreqMat a square symmetric matrix with values = (1-2*c1), where c1=matrix of expected recomb. frequencies. The choice to do 1-2c1 outside the function was made for computation efficiency; every operation on a big matrix takes time.
#' @param haploMat matrix of phased haplotypes, 2 rows per sample, cols = loci, {0,1}, rownames assumed to contain GIDs with a suffix, separated by "_" to distinguish haplotypes
#'
#' @details Columns of haploMat and row/cols of recombFreqMat should be in same order. May be worth computing in an R session using multi-threaded BLAS.
#' @return Potentially really large matrix representing the LD between loci in gametes
#' @family predCrossVar
calcGameticLD<-function(parentGID,recombFreqMat,haploMat){
  X<-haploMat[paste0(parentGID,c("_HapA","_HapB")),,drop=F]
  p<-colMeans(X);
  D<-recombFreqMat*((0.5*t(X)%*%X)-p%*%t(p));
  return(D) }

#' Calculate the cross LD matrix based on both parents
#'
#' Wraps the \code{\link{calcGameticLD}} function. CrossLD = sireLD + damLD.
#'
#' @param sireID string, Male parent genotype ID.
#' @param damID string, Female parent genotype ID.
#' @param recombFreqMat a square symmetric matrix with values = (1-2*c1), where c1=matrix of expected recomb. frequencies. The choice to do 1-2c1 outside the function was made for computation efficiency; every operation on a big matrix takes time.
#' @param haploMat matrix of phased haplotypes, 2 rows per sample, cols = loci, {0,1}, rownames assumed to contain GIDs with a suffix, separated by "_" to distinguish haplotypes
#'
#' @details Columns of haploMat and row/cols of recombFreqMat should be in same order. May be worth computing in an R session using multi-threaded BLAS.
#' @return Potentially really large matrix representing the LD between loci in the cross
#' @export
#' @family predCrossVar
calcCrossLD<-function(sireID,damID,recombFreqMat,haploMat){
  return(calcGameticLD(sireID,recombFreqMat,haploMat)+
           calcGameticLD(damID,recombFreqMat,haploMat)) }

#' Predict cross means
#'
#' Function to predict the mean performances of the offspring of crosses. Takes
#' a list of crosses to predict, marker effects, parental allele dosage matrix
#' as input. Predicts potentially over multiple crosses and multiple traits.
#' With \code{predType="BV"} predicts the mid-parent of crosses by computing
#' parental GEBV. With \code{predType="TGV"} predicts the mean total merit of
#' cross offspring using a Falconer-MacKay Eqn. 14.6 and takes user supplied
#' additive and dominance effects as input. The additive-dominance effects
#' should be partitioned according to the "genotypic" marker codings (see
#' Vitezica et al. 2013. GENETICS).
#'
#'
#' @param CrossesToPredict data.frame or tibble, col/colnames: sireID, damID.
#'   sireID and damID must both be in the haploMat.
#' @param predType string, "BV" or "TGV". "BV" predicts cross mean breeding
#'   values as the mean GEBV of parents. "TGV" predicts the cross total genetic
#'   value. Warning: prediction of meanTGV with F-M Eqn. 14.6 appropriate only
#'   using a+d partition not allele sub. + dom. dev.; genotypic NOT classical in
#'   terms used by Vitezica et al. 2013. For that reason,
#'   \code{\link{predCrossMeans}} has a "predType" not a "modelType" argument
#'   predType="TGV" uses Falconer-MacKay Eqn. 14.6 and takes add and dom
#'   effects. predType="BV" input should be allele subst. effs, computes
#'   mid-parent GEBV there is no equivalent to predicting the dominance variance
#'   for the mean thus the difference from the predCrossVars() function. NOTICE:
#'   NOT SAME as predType argument used in \code{\link{predCrossVars}}, sorry.
#' @param AddEffectList list of ADDITIVE effect matrices, one matrix per trait,
#'   Each element of the list is named with a string identifying the trait and
#'   the colnames of each matrix are labelled with snpIDs.
#' @param DomEffectList list of DOMINANCE effect matrices, one matrix per trait,
#'   Each element of the list is named with a string identifying the trait and
#'   the colnames of each matrix are labelled with snpIDs.
#' @param doseMat dosage matrix. required only for modelType=="DirDom". Assumes
#'   SNPs coded 0, 1, 2. Nind rows x Nsnp cols, numeric matrix, with rownames
#'   and colnames to indicate SNP/ind ID
#' @param ncores number of cores, parallelizes across \code{CrossesToPredict},
#'   in multi-trait cases, process traits for each family in serial within each
#'   worker.
#' @param ...
#'
#' @return tibble, each row contains predictions for a single cross. Columns:
#' \itemize{
#'  \item \code{"Trait"}:
#'  \item \code{"sireID"}:
#'  \item \code{"damID"}:
#'  \item \code{"sireGEBV"}: genomic estimated breeding value (GEBV) of the male parent of the cross
#'  \item \code{"damGEBV"}: genomic estimated breeding value (GEBV) of the female parent of the cross
#'  \item \code{"predOf"}: "MeanBV" or "MeanTGV"
#'  \item \code{"predMean"}: The predicted mean value for the cross
#' }
#' @export
#' @family predCrossVar
predCrossMeans<-function(CrossesToPredict,predType,
                         AddEffectList,DomEffectList=NULL,
                         doseMat,
                         ncores=1,
                         ...){
  #nBLASthreads=NULL,
  means<-tibble::tibble(Trait=names(AddEffectList))
  parents<-CrossesToPredict %$% union(sireID,damID)
  doseMat<-doseMat[parents,colnames(AddEffectList[[1]])]

  require(furrr); plan(multisession, workers = ncores)
  options(future.globals.maxSize=+Inf); options(future.rng.onMisuse="ignore")

  if(predType=="BV"){
    means<-means %>%
      mutate(predictedMeans=future_map(Trait,function(Trait,
                                                      #nBLASthreads=nBLASthreads,
                                                      ...){

        #if(!is.null(nBLASthreads)) { RhpcBLASctl::blas_set_num_threads(nBLASthreads) }

        parentGEBVs<-tcrossprod(doseMat,AddEffectList[[Trait]])

        predmeans<-CrossesToPredict %>%
          dplyr::left_join(tibble::tibble(sireID=rownames(parentGEBVs),
                           sireGEBV=as.numeric(parentGEBVs))) %>%
          dplyr::left_join(tibble::tibble(damID=rownames(parentGEBVs),
                           damGEBV=as.numeric(parentGEBVs))) %>%
          dplyr::mutate(predOf="MeanBV",
                 predMean=(sireGEBV+damGEBV)/2)
        return(predmeans) }))
    plan(sequential)
  }

  if(predType=="TGV"){
    plan(multisession, workers = ncores)
    means<-means %>%
      dplyr::mutate(predictedMeans=future_map(Trait,function(Trait,
                                                      #BLASthreads=nBLASthreads,
                                                      ...){

        #if(!is.null(nBLASthreads)) { RhpcBLASctl::blas_set_num_threads(nBLASthreads) }

        predmeans<-CrossesToPredict %>%
          dplyr::mutate(predOf="MeanTGV",
                 predMean=map2_dbl(sireID,damID,function(sireID,damID,...){
                   # Eqn 14.6 from Falconer+MacKay
                   p1<-doseMat[sireID,]/2
                   p2<-doseMat[damID,]/2
                   q<-1-p1
                   y<-p1-p2
                   g<-AddEffectList[[Trait]]*(p1-q-y) + DomEffectList[[Trait]]*((2*p1*q)+y*(p1-q))
                   meanG<-sum(g)
                   return(meanG)}))
        return(predmeans) }))
    plan(sequential)
  }
  means<-means %>%
    unnest(predictedMeans)
  return(means)
}
