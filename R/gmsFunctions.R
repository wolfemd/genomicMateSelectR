#' Run k-fold cross-validation
#'
#' Assess the accuracy of predicted previously unobserved genotypes
#' (individuals) based on the available training data. Runs k-fold
#' cross-validation for potentially multiple traits and optionally computing
#' prediction accuracy on user-specified selection index. Three models are
#' enabled: additive-only ("A"), additive-plus-dominance ("AD") and a
#' directional-dominance model that incorporates a genome-wide homozygosity
#' effect ("DirDom"). The union of all genotypes scored for all traits is broken
#' into k-folds a user specified number of times. Subsequently each train-test
#' pair is predicted for each trait and accuracies are computed.
#'
#' @param blups nested data.frame with list-column "TrainingData" containing
#'   BLUPs. Each element of "TrainingData" list, is data.frame with de-regressed
#'   BLUPs, BLUPs and weights (WT) for training and test.
#' @param modelType string, "A", "AD", "DirDom". modelType="A": additive-only,
#'   GEBVS modelType="AD": the "classic" add-dom model, GEBVS+GEDDs = GETGVs
#'   modelType="DirDom": the "genotypic" add-dom model with prop. homozygous fit
#'   as a fixed-effect, to estimate a genome-wide inbreeding effect. obtains
#'   add-dom effects, computes allele sub effects (\eqn{\alpha = a + d(q-p)})
#'   incorporates into GEBV and GETGV. "DirDom" requires dosages
#' @param selInd logical, TRUE/FALSE, selection index accuracy estimates,
#'   requires input weights via \code{SIwts}
#' @param SIwts required if \code{selInd=FALSE}, named vector of selection index
#'   weights, names match the "Trait" variable in \code{blups}
#' @param grms list of GRMs where each element is named either A, D, or, AD.
#'   Matrices supplied must match required by A, AD and ADE models. For ADE
#'   grms=list(A=A,D=D)
#' @param dosages dosage matrix. required only for modelType=="DirDom". Assumes
#'   SNPs coded 0, 1, 2. Nind rows x Nsnp cols, numeric matrix, with rownames
#'   and colnames to indicate SNP/ind ID
#' @param nrepeats number of repeats
#' @param nfolds number of folds
#' @param ncores number of cores, parallelizes across repeat-folds
#' @param nBLASthreads number of cores for each worker to use for multi-thread
#'   BLAS
#' @param gid string variable name used for genotype ID's/ in e.g. \code{blups}
#'   (default="GID")
#' @param seed numeric, use seed to achieve reproducibile train-test folds.
#' @param ...
#'
#' @return Returns tidy results in a tibble with accuracy estimates for each rep-fold in a list-column "accuracyEstOut".
#'
#' @export
#' @family CrossVal
runCrossVal<-function(blups,
                      modelType,
                      selInd,SIwts = NULL,
                      grms,dosages=NULL,
                      nrepeats,nfolds,
                      ncores=1,nBLASthreads=NULL,
                      gid="GID",seed=NULL,...){

  # SET-UP CROSS-VALIDATION TRAINING-TEST FOLDS

  # same train-test folds across traits
  gids<-blups %>%
    unnest(TrainingData) %>%
    distinct(!!sym(gid)) %>%
    .[[gid]]

  # whether or not user inputs a master seed
  # generate and store in output
  # seeds to make each replicate reproducible.
  if(!is.null(seed)){
    set.seed(seed);
    seeds<-sample(1:1e6,replace = F,size = nrepeats)
  } else {
    seeds-sample(1:1e6,replace = F,size = nrepeats) }

  # Set-up replicated cross-validation folds
  # splitting by clone (if clone in training dataset, it can't be in testing)
  require(rsample)
  cvsamples<-tibble(repeats=1:nrepeats,
                    seeds=seeds,
                    splits=map(seeds,function(seeds,...){
                      set.seed(seeds);
                      cvfolds<-vfold_cv(tibble(GID=gids),v = nfolds)
                      return(cvfolds)})) %>%
    unnest(splits)

  # FIT GENOMIC PREDICTION MODELS

  ## The fitModels() internal function
  ## Now runs _across_ traits and, if requested,
  ## computes the selection index accuracy
  ## runCrossVal() now parallelizes over repeat-folds using ncores
  ## Traits are handled in serial by each parallel worker
  ## nBLASthreads controls the number of additional cores each worker
  ## uses to speed matrix computations
  ## Internal function
  ## fits prediction model and calcs. accuracy for each train-test split

  fitModels<-function(splits,gids,
                      modelType,
                      blups,selInd,SIwts,
                      gid="GID",
                      grms,dosages=NULL,
                      nBLASthreads){
    # internal testing of fitModels() inputs - one rep-fold
    # splits<-cvsamples$splits[[1]]
    # rm(splits)

    if(!is.null(nBLASthreads)) { RhpcBLASctl::blas_set_num_threads(nBLASthreads) }
    # workers in plan(multisession) need this call internal to the function, it seems.

    # subset kinship (and if modelType=="DirDom" also the dosages) matrices
    ### only lines with BLUPs for cross-validation
    A<-grms[["A"]][gids,gids]
    if(modelType %in% c("AD","DirDom")){ D<-grms[["D"]][gids,gids]  }
    if(modelType=="DirDom"){
      snps<-dosages[gids,]
      Mmat<-centerDosage(snps);
      Dmat<-dose2domDevGenotypic(snps)
      f<-getPropHom(snps);
      p<-getAF(snps)
      q<-1-p
      rm(snps);
    }

    predictOneTrait<-possibly(function(TrainingData,splits,gid,
                                       modelType,A,D=NULL,
                                       Mmat=NULL,Dmat=NULL,f=NULL,p=NULL,q=NULL){
      #TrainingData<-blups$TrainingData[[1]]

      trainingdata<-TrainingData %>%
        dplyr::rename(GID=!!sym(gid)) %>%
        filter(GID %in% training(splits)[[gid]],
               GID %in% rownames(A))

      trainingdata[[paste0(gid,"a")]]<-factor(trainingdata[["GID"]],
                                              levels=rownames(A))
      if(modelType %in% c("AD")){
        trainingdata[[paste0(gid,"d")]]<-trainingdata[[paste0(gid,"a")]] }
      if(modelType %in% c("DirDom")){
        trainingdata[[paste0(gid,"d_star")]]<-trainingdata[[paste0(gid,"a")]] }

      # Set-up random model statements
      randFormula<-paste0("~vs(",gid,"a,Gu=A)")
      if(modelType %in% c("AD")){
        randFormula<-paste0(randFormula,"+vs(",gid,"d,Gu=D)") }
      if(modelType=="DirDom"){
        randFormula<-paste0(randFormula,"+vs(",gid,"d_star,Gu=D)")
        trainingdata %<>%
          left_join(tibble(GID=names(f),f=as.numeric(f)))
      }

      # Fixed model statements
      fixedFormula<-ifelse(modelType=="DirDom",
                           "drgBLUP ~1+f","drgBLUP ~1")
      # Fit genomic prediction model
      require(sommer)
      fit <- sommer::mmer(fixed = as.formula(fixedFormula),
                          random = as.formula(randFormula),
                          weights = WT,
                          data=trainingdata,
                          date.warning = F,
                          getPEV = FALSE)

      # reduce memory footprint
      rm(A); if(modelType %in% c("AD","DirDom")){ rm(D); gc() }

      print(paste0("GBLUP model complete - one trait"))
      if(modelType=="DirDom"){

        # Backsolve SNP effects
        # Compute allele sub effects
        ## Every model has an additive random term
        ga<-as.matrix(fit$U[[paste0("u:",gid,"a")]]$drgBLUP,ncol=1)

        # model DirDom is a different add-dom partition,
        ### add effects are not allele sub effects and gblups are not GEBV
        addsnpeff<-backsolveSNPeff(Z=Mmat,g=ga)
        ### dom effects are called d*, gd_star or domstar
        ### because of the genome-wide homoz. term included in model
        gd_star<-as.matrix(fit$U[[paste0("u:",gid,"d_star")]]$drgBLUP,ncol=1)
        domstar_snpeff<-backsolveSNPeff(Z=Dmat,g=gd_star)
        ### b = the estimate (BLUE) for the genome-wide homoz. effect
        b<-fit$Beta[fit$Beta$Effect=="f","Estimate"]
        ### calc. domsnpeff including the genome-wide homoz. effect
        ### divide the b effect up by number of SNPs and _subtract_ from domstar
        domsnpeff<-domstar_snpeff-(b/length(domstar_snpeff))

        ### allele substitution effects using a+d(q-p) where d=d*-b/p
        allelesubsnpeff<-addsnpeff+(domsnpeff*(q-p))
      }

      # Gather the GBLUPs
      if(modelType %in% c("A","AD")){
        gblups<-tibble(GID=as.character(names(fit$U[[paste0("u:",gid,"a")]]$drgBLUP)),
                       GEBV=as.numeric(fit$U[[paste0("u:",gid,"a")]]$drgBLUP)) }
      if(modelType=="AD"){
        gblups %<>% # compute GEDD (genomic-estimated dominance deviation)
          mutate(GEDD=as.numeric(fit$U[[paste0("u:",gid,"d")]]$drgBLUP),
                 # compute GETGV
                 GETGV=rowSums(.[,grepl("GE",colnames(.))])) }

      if(modelType=="DirDom"){
        # re-calc the GBLUP, GEdomval using dom. effects where d=d*-b/p
        ge_domval<-Dmat%*%domsnpeff
        # calc. the GEBV using allele sub. effects where alpha=a+d(p-q), and d=d*-b/p
        gebv<-Mmat%*%allelesubsnpeff
        # Tidy tibble of GBLUPs
        gblups<-tibble(GID=as.character(names(fit$U[[paste0("u:",gid,"a")]]$drgBLUP)),
                       GEadd=as.numeric(fit$U[[paste0("u:",gid,"a")]]$drgBLUP),
                       GEdom_star=as.numeric(fit$U[[paste0("u:",gid,"d_star")]]$drgBLUP)) %>%
          left_join(tibble(GID=rownames(ge_domval),GEdomval=as.numeric(ge_domval))) %>%
          left_join(tibble(GID=rownames(gebv),GEBV=as.numeric(gebv))) %>%
          # GETGV from GEadd + GEdomval
          mutate(GETGV=GEadd+GEdomval)

        # free up the memory footprint
        rm(ga,addsnpeff,gd_star,domstar_snpeff,b,domsnpeff,allelesubsnpeff,
           ge_domval,gebv,Dmat,Mmat,fit); gc()
        print(paste0("Backsolving SNP effects for DirDom model compete - one trait"))
      }

      # this is to remove conflicts with dplyr function select() downstream
      detach("package:sommer",unload = T); detach("package:MASS",unload = T)

      # Calculate accuracy for each trait
      ## Convert predicted gblups to a long-format
      gblups %<>%
        dplyr::select(GID,any_of(c("GEBV","GETGV"))) %>%
        pivot_longer(any_of(c("GEBV","GETGV")),
                     names_to = "predOf",
                     values_to = "GBLUP")

      ## Grab the test set BLUPs as validation data
      validationData<-TrainingData %>%
        dplyr::rename(GID=!!sym(gid)) %>%
        dplyr::select(GID,BLUP) %>%
        filter(GID %in% testing(splits)[[gid]])

      # Measure accuracy in test set
      ## cor(GEBV,BLUP)
      ## cor(GETGV,BLUP)
      accuracy<-gblups %>%
        left_join(validationData) %>%
        nest(predVSobs=c(GID,GBLUP,BLUP)) %>%
        mutate(Accuracy=map_dbl(predVSobs,~cor(.$GBLUP,.$BLUP, use = 'complete.obs')))
      return(accuracy)
    },
    otherwise = NA)

    # Predict for one trait to each trait's training dataset
    if(modelType=="A"){
      predictions<-blups %>%
        mutate(modelOut=map(TrainingData,~predictOneTrait(TrainingData=.,
                                                          splits=splits,gid=gid,
                                                          modelType=modelType,
                                                          A=A))) }
    if(modelType=="AD"){
      predictions<-blups %>%
        mutate(modelOut=map(TrainingData,~predictOneTrait(TrainingData=.,
                                                          splits=splits,gid=gid,
                                                          modelType=modelType,
                                                          A=A,D=D))) }
    if(modelType=="DirDom"){
      predictions<-blups %>%
        mutate(modelOut=map(TrainingData,~predictOneTrait(TrainingData=.,
                                                          splits=splits,gid=gid,
                                                          modelType=modelType,
                                                          A=A,D=D,
                                                          Mmat=Mmat,Dmat=Dmat,
                                                          f=f,p=p,q=q))) }

    rm(A); if(modelType=="AD"){ rm(D) };
    if(modelType=="DirDom"){ rm(D,Mmat,Dmat,f,p,q) }

    print(paste0("Genomic predictions done for all traits in one repeat-fold"))

    predictions %<>%
      select(-TrainingData) %>%
      unnest(modelOut,keep_empty = F)

    if(selInd){
      # calc. SELIND and SELIND accuracy

      gblups<-predictions %>%
        select(-Accuracy) %>%
        unnest(predVSobs) %>%
        select(-BLUP) %>%
        pivot_wider(values_from = "GBLUP",
                    names_from = "Trait")

      if(all(names(SIwts) %in% colnames(gblups))){
        gblups %<>%
          mutate(GBLUP=as.numeric((gblups %>%
                                     select(names(SIwts)) %>%
                                     as.matrix(.))%*%SIwts)) %>%
          select(predOf,GID,GBLUP)


        validationData<-blups %>%
          unnest(TrainingData) %>%
          select(Trait,GID,BLUP) %>%
          pivot_wider(names_from = "Trait", values_from = "BLUP")
        validationData %<>%
          mutate(BLUP=as.numeric((validationData %>%
                                    select(names(SIwts)) %>%
                                    as.matrix(.))%*%SIwts)) %>%
          select(GID,BLUP)

        predictions %<>%
          bind_rows(gblups %>%
                      left_join(validationData) %>%
                      nest(predVSobs=c(GID,GBLUP,BLUP)) %>%
                      mutate(Trait="SELIND") %>%
                      relocate(Trait,.before = 1) %>%
                      mutate(Accuracy=map_dbl(predVSobs,~cor(.$GBLUP,.$BLUP,
                                                             use = 'na.or.complete'))))
      }
    }

    predictions %<>%
      mutate(NcompleteTestPairs=map_dbl(predVSobs,function(predVSobs){
        if(!is.null(predVSobs)){
          out<-na.omit(.) %>% nrow(.) } else { out<-NA }
        return(out) }))

    return(predictions)

  }

  require(furrr); plan(multisession, workers = ncores)
  options(future.globals.maxSize=+Inf); options(future.rng.onMisuse="ignore")
  # quick test
  #cvsamples %<>% slice(1:2)
  cvsamples %<>%
    mutate(accuracyEstOut=future_map(splits,
                                     ~fitModels(splits=.,
                                                modelType=modelType,
                                                blups=blups,
                                                selInd=selInd,SIwts=SIwts,
                                                gid=gid,grms=grms,dosages=dosages,
                                                nBLASthreads=nBLASthreads)))
  plan(sequential)
  return(cvsamples)
}


#' Run genomic predictions
#'
#' Run GBLUP model using \code{\link[sommer]{mmer}}, potentially on multiple
#' traits. Returns genomic BLUPs (GEBV and GETGV). If requested, returns
#' backsolved marker effects (equivalent to ridge regression / SNP-BLUP).
#' Three models are
#' enabled: additive-only ("A"), additive-plus-dominance ("AD") and a
#' directional-dominance model that incorporates a genome-wide homozygosity
#' effect ("DirDom"). Inbreeding effect is included in output GEBV/GETGV
#' predictions *after* backsolving SNP effects. If requested, returns
#' GEBV/GETGV computed for a selection index using \code{selInd=TRUE}
#' and supplying \code{SIwts}.
#'
#' @param modelType string, "A", "AD", "DirDom".
#' modelType="A": additive-only, GEBVS
#' modelType="AD": the "classic" add-dom model, GEBVS+GEDDs = GETGVs
#' modelType="DirDom":
#' the "genotypic" add-dom model with prop. homozygous
#' fit as a fixed-effect, to estimate a genome-wide inbreeding effect.
#' obtains add-dom effects, computes
#' allele sub effects (\eqn{\alpha = a + d(q-p)})
#' incorporates into GEBV and GETGV
#' @param selInd logical, TRUE/FALSE, selection index accuracy estimates,
#' requires input weights via \code{SIwts}
#' @param SIwts required if \code{selInd=FALSE}, named vector of selection
#' index weights, names match the "Trait" variable in \code{blups}
#' @param getMarkEffs T/F return marker effects, backsolved from GBLUP
#' @param returnPEV T/F return PEVs in GBLUP
#' @param blups nested data.frame with list-column "TrainingData" containing
#'   BLUPs. Each element of "TrainingData" list, is data.frame with de-regressed
#'   BLUPs, BLUPs and weights (WT) for training and test.
#' @param grms list of genomic relation matrices (GRMs, aka kinship matrices).
#' Any genotypes in the GRMs get predicted with, or without phenotypes.
#' Each element is named either A or D. Matrices supplied must match
#' required by A, AD and DirDom models. e.g. grms=list(A=A,D=D).
#' @param dosages dosage matrix. required only for modelType=="DirDom".
#' Assumes SNPs coded 0, 1, 2. Nind rows x Nsnp
#' cols, numeric matrix, with rownames and colnames to indicate SNP/ind ID
#' @param gid string variable name used for genotype ID's in e.g. \code{blups} (default="GID")
#' @param ncores number of cores
#' @param nBLASthreads number of cores for each worker to use for multi-thread BLAS
#'
#' @return tibble, one row, two list columns (basically a named two-element
#' list of lists): \code{gblups[[1]]} and \code{genomicPredOut[[1]]}.
#' code{gblups[[1]]}: tibble of predicted GEBV/GETGV, all traits and potentially
#' SELIND genomic BLUPs along the columns.
#'
#' \code{genomicPredOut[[1]]} is a tibble that contains
#' some combination of lists-columns:
#' \itemize{
#'  \item gblups
#'  \item varcomps,
#'  \item fixeffs,
#'  \item allelesubsnpeff,
#'  \item addsnpeff,
#'  \item domstar_snpeff,
#'  \item domsnpeff
#' }
#'
#' @export
#' @family prediction_functions
runGenomicPredictions<-function(modelType,
                                selInd,SIwts = NULL,
                                getMarkEffs=FALSE,
                                returnPEV=FALSE,
                                blups,grms,dosages=NULL,gid="GID",
                                ncores=1,
                                nBLASthreads=NULL){

  fitModel<-function(trainingdata,modelType,getMarkEffs,returnPEV,
                     gid="GID",grms,dosages,
                     nBLASthreads,...){
    if(!is.null(nBLASthreads)) { RhpcBLASctl::blas_set_num_threads(nBLASthreads) }
    # workers in plan(multisession) need this call internal to the function, it seems.

    A<-grms[["A"]]
    if(modelType %in% c("AD","DirDom")){ D<-grms[["D"]] }

    trainingdata %<>%
      dplyr::rename(gid=!!sym(gid)) %>%
      filter(gid %in% rownames(A))

    trainingdata[[paste0(gid,"a")]]<-factor(trainingdata[["gid"]],
                                            levels=rownames(A))
    if(modelType %in% c("AD")){
      trainingdata[[paste0(gid,"d")]]<-trainingdata[[paste0(gid,"a")]] }
    if(modelType %in% c("DirDom")){
      trainingdata[[paste0(gid,"d_star")]]<-trainingdata[[paste0(gid,"a")]] }

    # Set-up random model statements
    randFormula<-paste0("~vs(",gid,"a,Gu=A)")
    if(modelType %in% c("AD")){
      randFormula<-paste0(randFormula,"+vs(",gid,"d,Gu=D)") }
    if(modelType=="DirDom"){
      randFormula<-paste0(randFormula,"+vs(",gid,"d_star,Gu=D)")
      f<-getPropHom(dosages)
      trainingdata %<>% mutate(f=f[trainingdata$gid]) }

    # Fixed model statements
    fixedFormula<-ifelse(modelType=="DirDom",
                        "drgBLUP ~1+f","drgBLUP ~1")
    # Fit genomic prediction model
    require(sommer)
    fit <- sommer::mmer(fixed = as.formula(fixedFormula),
                        random = as.formula(randFormula),
                        weights = WT,
                        data=trainingdata,
                        date.warning = F,
                        getPEV = returnPEV)

    # Shrink the memory footprint
    rm(grms,A); if(model %in% c("DirDom","AD")){ rm(D) };

    if(getMarkEffs==TRUE | modelType=="DirDom"){

      # Backsolve SNP effects
      # Compute allele sub effects
      ## Every model has an additive random term
      ga<-as.matrix(fit$U[[paste0("u:",gid,"a")]]$drgBLUP,ncol=1)
      M<-centerDosage(dosages)

      if(modelType %in% c("A","AD")){
        # models A and AD give add effects corresponding to allele sub effects
        allelesubsnpeff<-backsolveSNPeff(Z=M,g=ga) }

      if(modelType %in% c("AD")){
        # model AD the dom effects are dominance deviation effects
        gd<-as.matrix(fit$U[[paste0("u:",gid,"d")]]$drgBLUP,ncol=1)
        domdevsnpeff<-backsolveSNPeff(Z=dose2domDev(dosages),g=gd) }

      if(modelType %in% c("DirDom")){
        # model DirDom is a different add-dom partition,
        ### add effects are not allele sub effects and gblups are not GEBV
        addsnpeff<-backsolveSNPeff(Z=M,g=ga)
        ### dom effects are called d*, gd_star or domstar
        ### because of the genome-wide homoz. term included in model
        gd_star<-as.matrix(fit$U[[paste0("u:",gid,"d_star")]]$drgBLUP,ncol=1)
        domdevMat_genotypic<-dose2domDevGenotypic(dosages)
        domstar_snpeff<-backsolveSNPeff(Z=domdevMat_genotypic,g=gd_star)
        ### b = the estimate (BLUE) for the genome-wide homoz. effect
        b<-fit$Beta[fit$Beta$Effect=="f","Estimate"]
        ### calc. domsnpeff including the genome-wide homoz. effect
        ### divide the b effect up by number of SNPs and _subtract_ from domstar
        domsnpeff<-domstar_snpeff-(b/length(domstar_snpeff))

        ### allele substitution effects using a+d(q-p) where d=d*-b/p
        p<-getAF(dosages)
        q<-1-p
        allelesubsnpeff<-addsnpeff+(domsnpeff*(q-p))
      } }

    # Gather the GBLUPs
    if(modelType %in% c("A","AD")){
      gblups<-tibble(GID=as.character(names(fit$U[[paste0("u:",gid,"a")]]$drgBLUP)),
                     GEBV=as.numeric(fit$U[[paste0("u:",gid,"a")]]$drgBLUP))
      if(returnPEV){
        pev_bv<-diag((fit$PevU[[paste0("u:",gid,"a")]]$drgBLUP))
        gblups %<>% left_join(tibble(GID=names(pev_bv),PEVbv=pev_bv))
      }
    }

    if(modelType=="AD"){
      gblups %<>% # compute GEDD (genomic-estimated dominance deviation)
        mutate(GEDD=as.numeric(fit$U[[paste0("u:",gid,"d")]]$drgBLUP),
               # compute GETGV
               GETGV=rowSums(.[,grepl("GE",colnames(.))]))
      if(returnPEV){
        pev_dd<-diag((fit$PevU[[paste0("u:",gid,"d")]]$drgBLUP))
        gblups %<>% left_join(tibble(GID=names(pev_dd),PEVdd=pev_dd))
      }
    }
    if(modelType=="DirDom"){
      # re-calc the GBLUP, GEdomval using dom. effects where d=d*-b/p
      ge_domval<-domdevMat_genotypic%*%domsnpeff
      # calc. the GEBV using allele sub. effects where alpha=a+d(p-q), and d=d*-b/p
      gebv<-M%*%allelesubsnpeff
      # Tidy tibble of GBLUPs
      gblups<-tibble(GID=as.character(names(fit$U[[paste0("u:",gid,"a")]]$drgBLUP)),
                     GEadd=as.numeric(fit$U[[paste0("u:",gid,"a")]]$drgBLUP),
                     GEdom_star=as.numeric(fit$U[[paste0("u:",gid,"d_star")]]$drgBLUP)) %>%
        left_join(tibble(GID=rownames(ge_domval),GEdomval=as.numeric(ge_domval))) %>%
        left_join(tibble(GID=rownames(gebv),GEBV=as.numeric(gebv))) %>%
        # GETGV from GEadd + GEdomval
        mutate(GETGV=GEadd+GEdomval)
      if(returnPEV){
        pev_a<-diag((fit$PevU[[paste0("u:",gid,"a")]]$drgBLUP))
        pev_dstar<-diag((fit$PevU[[paste0("u:",gid,"d_star")]]$drgBLUP))
        gblups %<>%
          left_join(tibble(GID=names(pev_a),PEVd_star=pev_a)) %>%
          left_join(tibble(GID=names(pev_dstar),PEVd_star=pev_dstar))
      }
    }

    # Extract variance components
    varcomps<-summary(fit)$varcomp

    # Exract fixed effects
    # for modelType="DirDom", contains estimate of genome-wide homoz. effect
    fixeffs<-summary(fit)$betas

    # Shrink the memory footprint again
    rm(fit); gc()

    # Collect results
    results<-tibble(gblups=list(gblups),
                    varcomps=list(varcomps),
                    fixeffs=list(fixeffs))

    if(getMarkEffs==TRUE){
      # Add snpeffects to output
      results %<>% mutate(allelesubsnpeff=list(allelesubsnpeff))
      if(modelType=="AD"){ results %<>% mutate(domdevsnpeff=list(domdevsnpeff)) }
      if(modelType=="DirDom"){
        results %<>% mutate(addsnpeff=list(addsnpeff),
                            domstar_snpeff=list(domstar_snpeff),
                            domsnpeff=list(domsnpeff)) } }
    # this is to remove conflicts with dplyr function select() downstream
    detach("package:sommer",unload = T); detach("package:MASS",unload = T)
    # return results
    return(results)
  }

  #require(furrr); require(future.callr); plan(callr, workers = ncores)
  require(furrr); plan(multisession, workers = ncores)
  options(future.globals.maxSize=+Inf); options(future.rng.onMisuse="ignore")
  predictions<-blups %>%
    mutate(genomicPredOut=future_map(TrainingData,
                                     ~fitModel(trainingdata=.,
                                               modelType=modelType,
                                               getMarkEffs=getMarkEffs,
                                               returnPEV=returnPEV,
                                               gid=gid,
                                               grms=grms,
                                               dosages=dosages,
                                               nBLASthreads=nBLASthreads)))
  plan(sequential)

  predictions %<>%
    select(-TrainingData) %>%
    unnest(genomicPredOut)
  # tidy GBLUP output for e.g. breeders / selections
  gblups<-predictions %>%
    select(Trait,gblups) %>%
    unnest(gblups) %>%
    select(!!sym(gid),Trait,any_of(c("GEBV","GETGV"))) %>%
    pivot_longer(any_of(c("GEBV","GETGV")),
                 values_to = "GBLUP",
                 names_to = "predOf") %>%
    pivot_wider(names_from = "Trait",
                values_from = "GBLUP")
  if(selInd){
    # calc. SELIND and add to tidy output
    gblups %<>%
      mutate(SELIND=as.numeric((gblups %>%
                                  select(names(SIwts)) %>%
                                  as.matrix(.))%*%SIwts)) %>%
      relocate(SELIND, .after = predOf)
  }
  predictions<-tibble(gblups=list(gblups),
                      genomicPredOut=list(predictions))
  return(predictions)
}

#' Predict crosses
#'
#' Predict potentially for multiple traits, the means, variances and
#' trait-trait covariances in a set ofuser supplied crosses.l
#' If requested, computed the selection index means and variances.
#' Computes the usefulness criteria \eqn{UC_{parent}} and \eqn{UC_{variety}}
#' potentially with a user supplied standardized selection intensity value
#' \code{stdSelInt}. Output enables easy ranking of potential crosses.
#' This function takes the matrices of snpeffects output
#' (\code{genomicPredOut[[1]]}) from the \code{\link{runGenomicPredictions}}
#' function (when \code{getMarkEffs=TRUE}).
#' This is a wrapper function around \code{\link{predCrossVars} and
#' \link{predCrossMeans}}.
#'
#' @param modelType string, A, AD or DirDom. A and AD representing model with
#   Additive-only, Add. plus Dominance, respectively. **NEW**:
#   modelType="DirDom" includes a genome-wide homozygosity effect as in Xiang
#   et al. 2016, uses a different dominance GRM and will probably be a little
#   slower.#' @param stdSelInt
#' @param selInd logical, TRUE/FALSE, selection index accuracy estimates,
#' requires input weights via \code{SIwts}
#' @param CrossesToPredict data.frame or tibble, col/colnames: sireID, damID. sireID and damID must both be in the haploMat.
#' @param snpeffs the element \code{genomicPredOut[[1]]} of the output of
#' \code{\link{runGenomicPredictions}}.
#' @param dosages dosage matrix. required only for modelType=="DirDom".
#' Assumes SNPs coded 0, 1, 2. Nind rows x Nsnp
#' cols, numeric matrix, with rownames and colnames to indicate SNP/ind ID
#' @param haploMat matrix of phased haplotypes, 2 rows per sample, cols = loci, {0,1}, rownames assumed to contain GIDs with a suffix, separated by "_" to distinguish haplotypes
#' @param recombFreqMat a square symmetric matrix with values = (1-2*c1), where c1=matrix of expected recomb. frequencies. The choice to do 1-2c1 outside the function was made for computation efficiency; every operation on a big matrix takes time.
#' @param ncores number of cores
#' @param nBLASthreads number of cores for each worker to use for multi-thread BLAS
#' @param predTheMeans default: TRUE, t/f whether to predict cross means
#' @param predTheVars default: TRUE, t/f whether to predict cross vars
#'
#' @return tibble, one row, two list columns (basically a named two-element
#' list of lists): \code{tidyPreds[[1]]} and \code{rawPreds[[1]]}.
#' code{tidyPreds[[1]]}: tidy output, fewer details. sireID, damID, Nsegsnps, predOf,Trait, predMean, predVar, predSD, predUsefulnesstibble of predicted GEBV/GETGV, all traits and potentially
#' SELIND genomic BLUPs along the columns. \code{rawPreds[[1]]}: more detailed output, list of 2 ("predMeans" tibble and "predVars" tibble).
#'
#' @export
#' @family prediction_functions
predictCrosses<-function(modelType,
                         stdSelInt = 2.67,
                         selInd,SIwts = NULL,
                         CrossesToPredict,
                         snpeffs,dosages,
                         haploMat,recombFreqMat,
                         ncores=1,nBLASthreads=NULL,
                         predTheMeans=TRUE,
                         predTheVars=TRUE){
  ## Format SNP effect matrices ~~~~~~~~~~~~~~~~~~~~~~~~

  AlleleSubEffectList<-snpeffs$allelesubsnpeff %>%
    `names<-`(.,snpeffs$Trait) %>%
    map(.,~t(.))

  if(modelType=="AD"){
    # DomDevEffects for model "AD" to predict VarTGV = VarBV + VarDD
    DomDevEffectList=snpeffs$domdevsnpeff %>%
      `names<-`(.,snpeffs$Trait) %>%
      map(.,~t(.)) }

  if(modelType=="DirDom"){
    # AddEffectList + DomEffectList --> VarTGV; AlleleSubEffectList --> VarBV;
    AddEffectList<-snpeffs$addsnpeff %>%
      `names<-`(.,snpeffs$Trait) %>%
      map(.,~t(.))
    DomEffectList<-snpeffs$domsnpeff %>%
      `names<-`(.,snpeffs$Trait) %>%
      map(.,~t(.)) }
  # store raw mean and var preds
  if(predTheVars){
  ## Predict cross variances ~~~~~~~~~~~~~~~~~~~~~~~~
  print("Predicting cross variance parameters")
  if(modelType=="A"){
    predictedvars<-predCrossVars(CrossesToPredict=CrossesToPredict,
                                 AddEffectList=AlleleSubEffectList,
                                 modelType="A",
                                 haploMat=haploMat,
                                 recombFreqMat=recombFreqMat,
                                 ncores=ncores,nBLASthreads=nBLASthreads) %>%
      unnest(predVars) %>%
      mutate(predOf="VarBV") }
  if(modelType=="AD"){
    predictedvars<-predCrossVars(CrossesToPredict=CrossesToPredict,
                                 AddEffectList=AlleleSubEffectList,
                                 DomEffectList=DomDevEffectList,
                                 modelType="AD",
                                 haploMat=haploMat,
                                 recombFreqMat=recombFreqMat,
                                 ncores=ncores,nBLASthreads=nBLASthreads)
    predictedvars %<>%
      unnest(predVars) %>%
      mutate(predOf=ifelse(predOf=="VarA","VarBV","VarDD"))
  }
  if(modelType=="DirDom"){
    predictedvarTGV<-predCrossVars(CrossesToPredict=CrossesToPredict,
                                   AddEffectList=AddEffectList,
                                   DomEffectList=DomEffectList,
                                   modelType="AD", # no "DirDom" model in predCrossVars() nor is it needed
                                   haploMat=haploMat,
                                   recombFreqMat=recombFreqMat,
                                   ncores=ncores,nBLASthreads=nBLASthreads)
    predictedvarBV<-predCrossVars(CrossesToPredict=CrossesToPredict,
                                  AddEffectList=AlleleSubEffectList,
                                  DomEffectList=NULL,
                                  modelType="A", # no "DirDom" model in predCrossVars() nor is it needed
                                  haploMat=haploMat,
                                  recombFreqMat=recombFreqMat,
                                  ncores=ncores,nBLASthreads=nBLASthreads)
    predictedvars<-predictedvarBV %>%
      unnest(predVars) %>%
      mutate(predOf="VarBV") %>%
      bind_rows(predictedvarTGV %>%
                  unnest(predVars)) }
  }

  if(predTheMeans){
    ## Predict cross means ~~~~~~~~~~~~~~~~~~~~~~~~
  print("Predicting cross means")
  ### predict MeanBVs
  predictedmeans<-predCrossMeans(AddEffectList=AlleleSubEffectList,
                                 CrossesToPredict=CrossesToPredict,
                                 doseMat=dosages,
                                 ncores=ncores,
                                 nBLASthreads=nBLASthreads,
                                 predType="BV")
  if(modelType=="AD"){
    ### DO NOT predict MeanTGV ~but~ duplicate MeanBV as MeanTGV prediction
    ### there IS predVarTGV for this model, output predUC-TGV (i.e. UC_variety)
    predictedmeans %<>%
      bind_rows(predictedmeans %>% mutate(predOf="TGV")) }

  if(modelType=="DirDom"){
    ### predict MeanTGVs
    ####  Prediction of MeanTGV is only available for the DirDom model
    #### or a model with "genotypic" additive-dominance SNP effects
    #### As implemented, modelType="AD" is the "classical" partition (BVs+ DomDevs)
    predictedmeans %<>%
      bind_rows(predCrossMeans(AddEffectList=AddEffectList,
                               DomEffectList=DomEffectList,
                               CrossesToPredict=CrossesToPredict,
                               doseMat=dosages,
                               ncores=ncores,
                               nBLASthreads=nBLASthreads,
                               predType="TGV"))
  }

  }

  ## SIMPLIFIED, TIDY, CROSS-WISE OUTPUT ~~~~~~~~~~~~~~~~~~~~~~~~
  rawPreds<-list()
  if(predTheMeans){ rawPreds[["predMeans"]]<-list(predictedmeans) }
  if(predTheVars){ rawPreds[["predVars"]]<-list(predictedvars) }

  if(predTheMeans){
  ## tidy pred. means ~~~~~~
  predictedmeans %<>%
    mutate(predOf=gsub("Mean","",predOf),
           Trait2=Trait) %>% # to match with variance pred. output
    rename(Trait1=Trait) %>% # to match with variance pred. output
    select(sireID,damID,predOf,Trait1,Trait2,predMean)
}
  if(predTheVars){
    ## tidy pred. vars ~~~~~~
  predictedvars %<>%
    select(sireID,damID,Nsegsnps,predOf,Trait1,Trait2,predVar) %>%
    mutate(predOf=gsub("Var","",predOf))
  if(modelType=="AD"){
    predictedvars %<>%
      filter(predOf=="BV") %>%
      bind_rows(predictedvars %>%
                  pivot_wider(names_from = "predOf",
                              values_from = "predVar",
                              names_prefix = "predVar") %>%
                  mutate(predVar=predVarBV+predVarDD,
                         predOf="TGV") %>%
                  select(-predVarBV,-predVarDD))
  }
  if(modelType=="DirDom"){
    predictedvars %<>%
      filter(predOf=="BV") %>%
      bind_rows(predictedvars %>%
                  filter(predOf!="BV") %>%
                  pivot_wider(names_from = "predOf",
                              values_from = "predVar",
                              names_prefix = "predVar") %>%
                  mutate(predVar=predVarA+predVarD,
                         predOf="TGV") %>%
                  select(-predVarA,-predVarD))
  }
}
  ## SELECTION INDEX MEANS AND VARIANCES ~~~~~~~~~~~~~~~~~~~~~~~~
  #### Compute and add to tidy output, if requested
  if(selInd){
    print("Computing SELECTION INDEX means and variances.")
    if(predTheMeans){
    traits<-unique(predictedmeans$Trait1)
    ## Compute Mean SELIND
    predictedmeans %<>%
      select(-Trait2) %>%
      spread(Trait1,predMean) %>%
      select(sireID,damID,predOf,all_of(traits)) %>%
      mutate(SELIND=as.numeric((predictedmeans %>%
                                  select(-Trait2) %>%
                                  spread(Trait1,predMean) %>%
                                  select(all_of(names(SIwts))) %>%
                                  as.matrix(.))%*%SIwts)) %>%
      relocate(SELIND, .after = predOf) %>%
      pivot_longer(cols = c(SELIND,all_of(traits)),
                   names_to = "Trait1",
                   values_to = "predMean") %>%
      mutate(Trait2=Trait1) %>%
      select(sireID,damID,predOf,Trait1,Trait2,predMean)
    }

    if(predTheVars){
      ## Compute Var SELIND
    require(furrr); plan(multisession, workers = ncores)
    options(future.globals.maxSize=+Inf); options(future.rng.onMisuse="ignore")

    predictedvars %<>%
      nest(predVars=c(Trait1,Trait2,predVar)) %>%
      ## loop over each rep-fold-predOf-sireIDxdamID
      mutate(predVars=future_map(predVars,function(predVars,...){

        gmat<-predVars %>%
          pivot_wider(names_from = "Trait2",
                      values_from = "predVar") %>%
          column_to_rownames(var = "Trait1") %>%
          as.matrix
        gmat[lower.tri(gmat)]<-t(gmat)[lower.tri(gmat)]
        gmat %<>% .[names(SIwts),names(SIwts)]
        predSelIndVar<-SIwts%*%gmat%*%SIwts
        ## add sel index predictions to component trait
        ## var-covar predictions

        predVars<-tibble(Trait1="SELIND",Trait2="SELIND",
                         predVar=as.numeric(predSelIndVar)) %>%
          bind_rows(predVars)
        return(predVars) })) %>%
      unnest(predVars)
    plan(sequential)
  }
  }
  if(predTheMeans & predTheVars){
  ## USEFULNESS CRITERIA ~~~~~~~~~~~~~~~~~~~~~~~~
  tidyPreds<-predictedvars %>%
    inner_join(predictedmeans) %>%
    rename(Trait=Trait1) %>%
    select(sireID,damID,Nsegsnps,predOf,Trait,predMean,predVar) %>%
    mutate(predSD=sqrt(predVar),
           predUsefulness=predMean+(stdSelInt*predSD)) }
  if(predTheMeans & !predTheVars){
    tidyPreds<-predictedvars %>%
      rename(Trait=Trait1) %>%
      select(sireID,damID,Nsegsnps,predOf,Trait,predVar)
  }
  if(!predTheMeans & predTheVars){
    tidyPreds<-predictedmeans %>%
      select(sireID,damID,Nsegsnps,predOf,Trait,predMean) %>%
      mutate(predSD=sqrt(predVar))
  }
  predictions<-tibble(tidyPreds=list(tidyPreds),
                      rawPreds=list(rawPreds))
  return(predictions)
}

