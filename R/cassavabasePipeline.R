#' Read database downloads into R
#'
#' @param phenotypeFile CSV file of field trial phenotype data downloaded from e.g. cassavabase
#' @param metadataFile CSV file of field trial metadata downloaded from e.g. cassavabase
#'
#' @details \strong{NOTICE:} This function is part of a family of functions (\code{"cassavabase_pheno_pipeline"}) developed as part of the NextGen Cassava Breeding Project genomic selection pipeline.
#' For some examples of their useage:
#' \itemize{
#'  \item \href{https://wolfemd.github.io/IITA_2021GS/}{IITA_2021GS (\strong{\emph{should be first use of this package}})}
#'  \item \href{https://wolfemd.github.io/NRCRI_2021GS/}{NRCRI_2021GS}
#'  \item \href{https://wolfemd.github.io/TARI_2020GS/}{TARI_2020GS}
#'  \item \href{https://wolfemd.github.io/IITA_2020GS/}{IITA_2020GS}
#' }
#' @family cassavabase_pheno_pipeline
#' @return merged data.frame phenos + metadata
#' @export
readDBdata<-function(phenotypeFile,metadataFile=NULL){
  indata<-read.csv(phenotypeFile,
                   na.strings = c("#VALUE!",NA,".",""," ","-","\""),
                   stringsAsFactors = F)
  if(!is.null(metadataFile)){
    meta<-read.csv(metadataFile,
                   na.strings = c("#VALUE!",NA,".",""," ","-","\""),
                   stringsAsFactors = F) %>%
      rename(programName=breedingProgramName,
             programDescription=breedingProgramDescription,
             programDbId=breedingProgramDbId)
    indata<-left_join(indata,meta) }
  indata %<>%
    filter(observationLevel=="plot")
  return(indata) }

#' Add TrialType classifier
#'
#'
#' @param indata data.frame, database data, e.g. read by \code{\link{readDBdata}}
#'
#' @details \strong{NOTICE:} This function is part of a family of functions (\code{"cassavabase_pheno_pipeline"}) developed as part of the NextGen Cassava Breeding Project genomic selection pipeline.
#' For some examples of their useage:
#' \itemize{
#'  \item \href{https://wolfemd.github.io/IITA_2021GS/}{IITA_2021GS (\strong{\emph{should be first use of this package}})}
#'  \item \href{https://wolfemd.github.io/NRCRI_2021GS/}{NRCRI_2021GS}
#'  \item \href{https://wolfemd.github.io/TARI_2020GS/}{TARI_2020GS}
#'  \item \href{https://wolfemd.github.io/IITA_2020GS/}{IITA_2020GS}
#' }
#' @family cassavabase_pheno_pipeline
#' @return same as "indata", adds column "TrialType"
#' @export
makeTrialTypeVar<-function(indata){
  # So far, this function is not very general
  # Handles IITA and NRCRI trial names as of September 2020.
  # Can customize this or add lines to grab TrialTypes for each breeding program
  if(indata$programName=="IITA"){
    outdata<-indata %>%
      mutate(TrialType=ifelse(grepl("CE|clonal|13NEXTgenC1",studyName,ignore.case = T),"CET",NA),
             TrialType=ifelse(grepl("EC",studyName,ignore.case = T),"ExpCET",TrialType),
             TrialType=ifelse(grepl("PYT",studyName,ignore.case = T),"PYT",TrialType),
             TrialType=ifelse(grepl("AYT",studyName,ignore.case = T),"AYT",TrialType),
             TrialType=ifelse(grepl("UYT",studyName,ignore.case = T),"UYT",TrialType),
             TrialType=ifelse(grepl("geneticgain|gg|genetic gain",studyName,ignore.case = T),"GeneticGain",TrialType),
             TrialType=ifelse(grepl("Cassava",studyName,ignore.case = T) & grepl("/",studyName),"GeneticGain",TrialType),
             TrialType=ifelse((grepl("clonal evaluation trial",!grepl("genetic gain",studyDescription,ignore.case = T),
                                     ignore.case = T)),"CET",TrialType),
             TrialType=ifelse(grepl("preliminary yield trial",studyDescription,ignore.case = T),"PYT",TrialType),
             TrialType=ifelse(grepl("Crossingblock|\\.CB\\.|cross",studyName) & is.na(TrialType),"CrossingBlock",TrialType),
             TrialType=ifelse(grepl("NCRP",studyName) & is.na(TrialType),"NCRP",TrialType),
             TrialType=ifelse(grepl("conservation",studyName) & is.na(TrialType),"Conservation",TrialType),
             TrialType=ifelse(grepl("seedling|\\.SN",studyName),"SN",TrialType)) }
  if(indata$programName=="NRCRI"){
    outdata<-indata %>%
      mutate(TrialType=ifelse(grepl("TP1",studyName,ignore.case = T),"TP1",NA),
             TrialType=ifelse(grepl("TP2",studyName,ignore.case = T),"TP2",TrialType),
             TrialType=ifelse(grepl("C1a",studyName,ignore.case = T),"C1a",TrialType),
             TrialType=ifelse(grepl("C1b",studyName,ignore.case = T),"C1b",TrialType),
             TrialType=ifelse(grepl("C2a",studyName,ignore.case = T),"C2a",TrialType),
             TrialType=ifelse(grepl("C2b",studyName,ignore.case = T),"C2b",TrialType),
             TrialType=ifelse(grepl("NCRP",studyName) & is.na(TrialType),"NCRP",TrialType),
             TrialType=ifelse(grepl("15nextgen60gs-cbUM|crossnblk|crossingblock",studyName,ignore.case = T) &
                                !grepl("CET",studyName),
                              "CrossingBlock",TrialType),
             TrialType=ifelse(grepl("seedling",studyName,ignore.case = T),NA,TrialType)) }

  if(indata$programName=="TARI"){
    outdata<-indata %>%
      mutate(TrialType=ifelse(grepl("Advanced Yield|AYT", trialType, ignore.case = T),"AYT",NA),
             TrialType=ifelse(grepl("Clonal|CET", trialType, ignore.case = T),"CET",TrialType),
             TrialType=ifelse(grepl("Preliminary|PYT", trialType, ignore.case = T),"PYT",TrialType),
             TrialType=ifelse(grepl("Regional", trialType, ignore.case = T),"RegionalTrial",TrialType),
             TrialType=ifelse(grepl("Uniform|UYT", trialType, ignore.case = T),"UYT",TrialType),
             TrialType=ifelse(grepl("Variety Release", trialType, ignore.case = T),"VarietyTrial",TrialType),
             TrialType=ifelse(grepl("CROSSING", trialType, ignore.case = T),"CrossingBlock",TrialType),
             TrialType=ifelse(grepl("GWAS", trialType, ignore.case = T),"GWAS",TrialType)) }

  return(outdata) }

#' Rename columns and remove everything unecessary
#'
#' @param traitabbrevs data.frame with 2 cols (TraitAbbrev and TraitName). TraitName should match exactly to cassava ontology names
#' @param indata data.frame read from cassavabase download
#' @param customColsToKeep char. vec. of any custom cols you added and want to keep
#'
#' @details \strong{NOTICE:} This function is part of a family of functions (\code{"cassavabase_pheno_pipeline"}) developed as part of the NextGen Cassava Breeding Project genomic selection pipeline.
#' For some examples of their useage:
#' \itemize{
#'  \item \href{https://wolfemd.github.io/IITA_2021GS/}{IITA_2021GS (\strong{\emph{should be first use of this package}})}
#'  \item \href{https://wolfemd.github.io/NRCRI_2021GS/}{NRCRI_2021GS}
#'  \item \href{https://wolfemd.github.io/TARI_2020GS/}{TARI_2020GS}
#'  \item \href{https://wolfemd.github.io/IITA_2020GS/}{IITA_2020GS}
#' }
#' @family cassavabase_pheno_pipeline
#' @return
#' @export
renameAndSelectCols<-function(traitabbrevs,indata,
                              customColsToKeep=NULL){
  outdata<-indata %>%
    select(any_of(c("studyYear","programName","locationName","studyName","studyDesign",
                    "plotWidth","plotLength","fieldSize","plantingDate","harvestDate",
                    "germplasmName","observationUnitDbId",
                    "replicate","blockNumber","plotNumber","rowNumber","colNumber","entryType",
                    "trialType","plantsPerPlot","numberBlocks","numberReps")),
           any_of(customColsToKeep),
           any_of(traitabbrevs$TraitName)) %>% ungroup() %>%
    mutate(across(any_of(traitabbrevs$TraitName), as.numeric)) %>% ungroup() %>%
    pivot_longer(cols = any_of(traitabbrevs$TraitName),
                 names_to = "TraitName",
                 values_to = "Value") %>%
    left_join(.,traitabbrevs) %>%
    select(-TraitName) %>%
    pivot_wider(names_from = TraitAbbrev,
                values_from = "Value")
  return(outdata) }

#' Nest data.frame by trials
#'
#' Also, Create some explicitly nested variables including loc and year to nest with the trial data
#'
#' @param indata
#'
#' @details \strong{NOTICE:} This function is part of a family of functions (\code{"cassavabase_pheno_pipeline"}) developed as part of the NextGen Cassava Breeding Project genomic selection pipeline.
#' For some examples of their useage:
#' \itemize{
#'  \item \href{https://wolfemd.github.io/IITA_2021GS/}{IITA_2021GS (\strong{\emph{should be first use of this package}})}
#'  \item \href{https://wolfemd.github.io/NRCRI_2021GS/}{NRCRI_2021GS}
#'  \item \href{https://wolfemd.github.io/TARI_2020GS/}{TARI_2020GS}
#'  \item \href{https://wolfemd.github.io/IITA_2020GS/}{IITA_2020GS}
#' }
#' @family cassavabase_pheno_pipeline
#' @return
#' @export
nestByTrials<-function(indata){
  nested_indata<-indata %>%
    # Create some explicitly nested variables including loc and year to nest with the trial data
    mutate(yearInLoc=paste0(programName,"_",locationName,"_",studyYear),
           trialInLocYr=paste0(yearInLoc,"_",studyName),
           repInTrial=paste0(trialInLocYr,"_",replicate),
           blockInRep=paste0(repInTrial,"_",blockNumber)) %>%
    nest(TrialData=-c(programName,locationName,studyYear,TrialType,studyName))
  return(nested_indata)
}

#' "Detect" experimental designs
#'
#' After cleaning up cassavabase trial data, the next step is to run some \emph{ad hoc} code that checks the experimental design of each trial.
#' If you are absolutely certain of the usage of the design variables in your dataset, you might not need this step.
#'
#' Returns the input data.frame with two new columns indicating TRUE/FALSE, for each trial (location-year-studyName),
#' whether the trial \code{replicate / repInTrial} columns indicates \strong{\code{CompleteBlocks}} and
#' if the \code{blockInRep / blockNumber} contains \strong{\code{IncompleteBlocks}}.
#'
#' Examples of reasons to do the step below:
#' \itemize{
#'  \item Some trials appear to be complete blocked designs and the blockNumber is used instead of replicate, which is what most use.
#'  \item Some complete block designs have nested, incomplete sub-blocks, others simply copy the "replicate" variable into the "blockNumber variable"
#'  \item  Some trials have only incomplete blocks \emph{but} the incomplete block info might be in the replicate \emph{and/or} the blockNumber column
#' }
#'
#' One reason it might be important to get this right is that the variance among complete blocks might not be the same among incomplete blocks. If we treat a mixture of complete and incomplete blocks as part of the same random-effect (replicated-within-trial), we assume they have the same variance.
#'
#' @param indata
#'
#' @details \strong{NOTICE:} This function is part of a family of functions (\code{"cassavabase_pheno_pipeline"}) developed as part of the NextGen Cassava Breeding Project genomic selection pipeline.
#' For some examples of their useage:
#' \itemize{
#'  \item \href{https://wolfemd.github.io/IITA_2021GS/}{IITA_2021GS (\strong{\emph{should be first use of this package}})}
#'  \item \href{https://wolfemd.github.io/NRCRI_2021GS/}{NRCRI_2021GS}
#'  \item \href{https://wolfemd.github.io/TARI_2020GS/}{TARI_2020GS}
#'  \item \href{https://wolfemd.github.io/IITA_2020GS/}{IITA_2020GS}
#' }
#' @family cassavabase_pheno_pipeline
#' @return Returns the input data.frame with two new columns indicating TRUE/FALSE, for each trial (location-year-studyName),
#' whether the trial \code{replicate / repInTrial} columns indicates \strong{\code{CompleteBlocks}} and
#' if the \code{blockInRep / blockNumber} contains \strong{\code{IncompleteBlocks}}.
#' @export
detectExptDesigns<-function(indata){
  # nestByTrials
  nestedDBdata<-indata %>%
    # Create some explicitly nested variables including loc and year to nest with the trial data
    mutate(yearInLoc=paste0(programName,"_",locationName,"_",studyYear),
           trialInLocYr=paste0(yearInLoc,"_",studyName),
           repInTrial=paste0(trialInLocYr,"_",replicate),
           blockInRep=paste0(repInTrial,"_",blockNumber)) %>%
    nest(TrialData=-c(programName,locationName,studyYear,TrialType,studyName))

  # Define complete blocks
  nestedDBdata %>%
    mutate(Nobs=map_dbl(TrialData,~nrow(.)),
           MaxNOHAV=map_dbl(TrialData,~unique(.$MaxNOHAV)),
           Nrep=map_dbl(TrialData,~length(unique(.$replicate))),
           Nblock=map_dbl(TrialData,~length(unique(.$blockInRep))),
           Nclone=map_dbl(TrialData,~length(unique(.$germplasmName))),
           # median number of obs per clone
           medObsPerClone=map_dbl(TrialData,~count(.,germplasmName) %$% round(median(n),1)),
           # median number of obs per replicate
           medObsPerRep=map_dbl(TrialData,~count(.,replicate) %$% round(median(n),1)),
           # Define complete block effects based on the "replicate" variable
           CompleteBlocks=ifelse(Nrep>1 & medObsPerClone==Nrep & Nobs!=Nrep,TRUE,FALSE),
           # Additional trials with imperfect complete blocks
           CompleteBlocks=ifelse(Nrep>1 & medObsPerClone!=Nrep & medObsPerClone>1 & Nobs!=Nrep,TRUE,CompleteBlocks)) -> x
  x %>%
    # Some complete blocks may only be represented by the "blockNumber" column
    mutate(medBlocksPerClone=map_dbl(TrialData,~select(.,blockInRep,germplasmName) %>%
                                       # median number of blockInRep per clone
                                       distinct %>%
                                       count(germplasmName) %$%
                                       round(median(n))),
           # If CompleteBlocks==FALSE (complete blocks not detected based on replicate)
           # and if more than half the clones are represented in more than one block based on the blockInRep variable
           # Copy the blockInRep values into the repInTrial column
           # Recompute Nrep
           # and declare CompleteBlocks==TRUE
           TrialData=ifelse(medBlocksPerClone>1 & CompleteBlocks==FALSE,map(TrialData,~mutate(.,repInTrial=blockInRep)),TrialData),
           Nrep=map_dbl(TrialData,~length(unique(.$repInTrial))),
           CompleteBlocks=ifelse(medBlocksPerClone>1 & CompleteBlocks==FALSE,TRUE,CompleteBlocks)) -> y
  # Define incomplete blocks
  y %>%
    mutate(repsEqualBlocks=map_lgl(TrialData,~all(.$replicate==.$blockNumber)),
           NrepEqualNblock=ifelse(Nrep==Nblock,TRUE,FALSE),
           medObsPerBlockInRep=map_dbl(TrialData,~count(.,blockInRep) %$% round(median(n),1))) -> z
  # Define complete blocked trials with nested sub-blocks
  z %<>%
    mutate(IncompleteBlocks=ifelse(CompleteBlocks==TRUE & Nobs!=Nblock & Nblock>1 & medObsPerBlockInRep>1 & NrepEqualNblock==FALSE,TRUE,FALSE))
  # Define clearly unreplicated (CompleteBlocks==FALSE & Nrep==1) trials with nested sub-blocks
  z %<>%
    mutate(IncompleteBlocks=ifelse(CompleteBlocks==FALSE & Nobs!=Nblock & Nblock>1 & medObsPerBlockInRep>1 & Nrep==1,TRUE,IncompleteBlocks))
  # Define additional trials with incomplete blocks (blockInRep) where CompleteBlocks==FALSE but Nrep>1 and Nrep==Block
  z %<>%
    mutate(IncompleteBlocks=ifelse(CompleteBlocks==FALSE & IncompleteBlocks==FALSE &
                                     Nobs!=Nblock & Nblock>1 &  Nobs!=Nrep &
                                     medObsPerBlockInRep>1 & Nrep>1 & NrepEqualNblock==TRUE,TRUE,IncompleteBlocks))
  # Last few cases (2 trials actually) where Nrep>1 and Nblock>1 and Nrep!=Nblock but CompleteBlocks==FALSE
  z %<>%
    mutate(IncompleteBlocks=ifelse(CompleteBlocks==FALSE & IncompleteBlocks==FALSE &
                                     Nobs!=Nblock & Nobs!=Nrep &
                                     medObsPerBlockInRep>1 & Nrep>1,TRUE,IncompleteBlocks))
  z %<>%
    dplyr::select(-MaxNOHAV) %>%
    unnest(TrialData)
  return(z)
}

#' Set-up for a by-trait analysis downstream
#'
#' Need to restructure the data from per-trial by regrouping by trait.
#'
#' @param indata data.frame with experimental designs detected, output from \code{\link{detectExptDesigns}}
#' @param traits character vec, names of traits to include
#'
#' @details \strong{NOTICE:} This function is part of a family of functions (\code{"cassavabase_pheno_pipeline"}) developed as part of the NextGen Cassava Breeding Project genomic selection pipeline.
#' For some examples of their useage:
#' \itemize{
#'  \item \href{https://wolfemd.github.io/IITA_2021GS/}{IITA_2021GS (\strong{\emph{should be first use of this package}})}
#'  \item \href{https://wolfemd.github.io/NRCRI_2021GS/}{NRCRI_2021GS}
#'  \item \href{https://wolfemd.github.io/TARI_2020GS/}{TARI_2020GS}
#'  \item \href{https://wolfemd.github.io/IITA_2020GS/}{IITA_2020GS}
#' }
#' @family cassavabase_pheno_pipeline
#' @return
#' @export
nestDesignsDetectedByTraits<-function(indata,traits){
  indata %<>%
    select(programName,locationName,studyYear,TrialType,studyName,
           CompleteBlocks,IncompleteBlocks,
           yearInLoc,trialInLocYr,repInTrial,blockInRep,observationUnitDbId,
           germplasmName,FullSampleName,GID,all_of(traits),PropNOHAV) %>%
    mutate(IncompleteBlocks=ifelse(IncompleteBlocks==TRUE,"Yes","No"),
           CompleteBlocks=ifelse(CompleteBlocks==TRUE,"Yes","No")) %>%
    pivot_longer(cols = all_of(traits), names_to = "Trait", values_to = "Value") %>%
    filter(!is.na(Value),
           !is.na(GID)) %>%
    nest(MultiTrialTraitData=c(-Trait))
  return(indata)
}

#' Nest trials by trait
#'
#' @param indata
#' @param traits character vec, names of traits to include
#'
#' @details \strong{NOTICE:} This function is part of a family of functions (\code{"cassavabase_pheno_pipeline"}) developed as part of the NextGen Cassava Breeding Project genomic selection pipeline.
#' For some examples of their useage:
#' \itemize{
#'  \item \href{https://wolfemd.github.io/IITA_2021GS/}{IITA_2021GS (\strong{\emph{should be first use of this package}})}
#'  \item \href{https://wolfemd.github.io/NRCRI_2021GS/}{NRCRI_2021GS}
#'  \item \href{https://wolfemd.github.io/TARI_2020GS/}{TARI_2020GS}
#'  \item \href{https://wolfemd.github.io/IITA_2020GS/}{IITA_2020GS}
#' }
#' @family cassavabase_pheno_pipeline
#' @return
#' @export
nestTrialsByTrait<-function(indata,traits){
  nested_trialdata<-dbdata %>%
    select(-MaxNOHAV) %>%
    unnest(TrialData) %>%
    pivot_longer(cols = any_of(traits),
                 names_to = "Trait",
                 values_to = "TraitValue") %>%
    nest(TraitByTrialData=-c(Trait,studyYear,programName,locationName,studyName,TrialType))
  return(nested_trialdata)
}

#' Calcluate the proportion missing data
#'
#' @param TraitValues numeric vector
#'
#' @details \strong{NOTICE:} This function is part of a family of functions (\code{"cassavabase_pheno_pipeline"}) developed as part of the NextGen Cassava Breeding Project genomic selection pipeline.
#' For some examples of their useage:
#' \itemize{
#'  \item \href{https://wolfemd.github.io/IITA_2021GS/}{IITA_2021GS (\strong{\emph{should be first use of this package}})}
#'  \item \href{https://wolfemd.github.io/NRCRI_2021GS/}{NRCRI_2021GS}
#'  \item \href{https://wolfemd.github.io/TARI_2020GS/}{TARI_2020GS}
#'  \item \href{https://wolfemd.github.io/IITA_2020GS/}{IITA_2020GS}
#' }
#' @family cassavabase_pheno_pipeline
#' @return
#' @export
calcPropMissing<-function(TraitValues){ length(which(is.na(TraitValues))) / length(TraitValues) }

#' Fit model and remove outliers from one trial for one trait
#'
#' \strong{NOT CURRENTLY IN USE / PROBABLY DEPRECATED}
#'
#' Fits mixed-models to each trial and removes outliers. Returns de-regressed BLUPs and more.
#'
#' @param Trait
#' @param TraitByTrialData
#' @param GID
#'
#' @details \strong{NOTICE:} This function is part of a family of functions (\code{"cassavabase_pheno_pipeline"}) developed as part of the NextGen Cassava Breeding Project genomic selection pipeline.
#' For some examples of their useage:
#' \itemize{
#'  \item \href{https://wolfemd.github.io/IITA_2021GS/}{IITA_2021GS (\strong{\emph{should be first use of this package}})}
#'  \item \href{https://wolfemd.github.io/NRCRI_2021GS/}{NRCRI_2021GS}
#'  \item \href{https://wolfemd.github.io/TARI_2020GS/}{TARI_2020GS}
#'  \item \href{https://wolfemd.github.io/IITA_2020GS/}{IITA_2020GS}
#' }
#' @family cassavabase_pheno_pipeline
#' @return
#' @export
curateTrialOneTrait<-function(Trait,TraitByTrialData,GID="GID"){
  require(lme4)

  modelFormula<-paste0("TraitValue ~ (1|",GID,")")
  modelFormula<-ifelse(all(TraitByTrialData$CompleteBlocks),
                       paste0(modelFormula,"+(1|repInTrial)"),modelFormula)
  modelFormula<-ifelse(all(TraitByTrialData$IncompleteBlocks),
                       paste0(modelFormula,"+(1|blockInRep)"),modelFormula)
  modelFormula<-ifelse(grepl("logRTNO",Trait) | grepl("logFYLD",Trait) | grepl("logTOPYLD",Trait),
                       paste0(modelFormula,"+PropNOHAV"),modelFormula)

  propMiss<-calcPropMissing(TraitByTrialData$TraitValue)
  fit_model<-possibly(function(modelFormula,TraitByTrialData){
    model_out<-lmer(as.formula(modelFormula),data=TraitByTrialData)
    if(!is.na(model_out)){
      outliers<-which(abs(rstudent(model_out))>=3.3)
      if(length(outliers)>0){
        model_out<-lmer(as.formula(modelFormula),data=TraitByTrialData,
                        subset=abs(rstudent(model_out))<3.3)
      }
    }
    return(list(model_out=model_out,outliers=outliers)) },
    otherwise = NA)
  model_out<-fit_model(modelFormula,TraitByTrialData)
  if(is.na(model_out)){
    out <-tibble(H2=NA,VarComps=list(NULL),BLUPs=list(NULL),Model=modelFormula,Noutliers=NA,Outliers=NA,propMiss=propMiss)
  } else {
    varcomps<-as.data.frame(VarCorr(model_out[["model_out"]]))[,c("grp","vcov")] %>%
      spread(grp,vcov)
    Vg<-varcomps$GID
    H2<-Vg/(Vg+varcomps$Residual)
    BLUP<-ranef(model_out[["model_out"]], condVar=TRUE)[[GID]]
    PEV <- c(attr(BLUP, "postVar"))
    blups<-tibble(GID=rownames(BLUP),BLUP=BLUP$`(Intercept)`,PEV=PEV) %>%
      mutate(REL=1-(PEV/Vg),
             drgBLUP=BLUP/REL,
             WT=(1-H2)/((0.1 + (1-REL)/REL)*H2))
    out <- tibble(H2=H2,
                  VarComps=list(varcomps),
                  BLUPs=list(blups),
                  Model=modelFormula,
                  Noutliers=length(model_out[["outliers"]]),
                  Outliers=list(model_out[["outliers"]]),
                  propMiss=propMiss) }
  return(out)
}

#' Fits mixed-models to each trial and all traits, removing outliers
#'
#'
#' \strong{NOT CURRENTLY IN USE / PROBABLY DEPRECATED}
#'
#' Fits mixed-models across trials and removes outliers. Returns de-regressed BLUPs and more.
#'
#' @param nestedTrialData
#' @param traits
#'
#' @details \strong{NOTICE:} This function is part of a family of functions (\code{"cassavabase_pheno_pipeline"}) developed as part of the NextGen Cassava Breeding Project genomic selection pipeline.
#' For some examples of their useage:
#' \itemize{
#'  \item \href{https://wolfemd.github.io/IITA_2021GS/}{IITA_2021GS (\strong{\emph{should be first use of this package}})}
#'  \item \href{https://wolfemd.github.io/NRCRI_2021GS/}{NRCRI_2021GS}
#'  \item \href{https://wolfemd.github.io/TARI_2020GS/}{TARI_2020GS}
#'  \item \href{https://wolfemd.github.io/IITA_2020GS/}{IITA_2020GS}
#' }
#' @family cassavabase_pheno_pipeline
#' @return
#' @export
curateTrialsByTrait<-function(nestedTrialData,traits){
  outdata<-nestedTrialData %>%
    mutate(modelOutput=map2(Trait,TraitByTrialData,~curateTrialOneTrait(Trait = .x,TraitByTrialData = .y))) %>%
    dplyr::select(-TraitByTrialData) %>%
    unnest(modelOutput)
  return(outdata)
}

#' Restructure BLUPs for downstream analysis
#'
#' \strong{NOT CURRENTLY IN USE / PROBABLY DEPRECATED}
#'
#' Actually expects input to be BLUPs from single-trial mixed-models.
#' Fits mixed-models across trials and removes outliers. Returns de-regressed BLUPs and more.
#'
#' @param curatedTrialData
#'
#' @return
#' @details \strong{NOTICE:} This function is part of a family of functions (\code{"cassavabase_pheno_pipeline"}) developed as part of the NextGen Cassava Breeding Project genomic selection pipeline.
#' For some examples of their useage:
#' \itemize{
#'  \item \href{https://wolfemd.github.io/IITA_2021GS/}{IITA_2021GS (\strong{\emph{should be first use of this package}})}
#'  \item \href{https://wolfemd.github.io/NRCRI_2021GS/}{NRCRI_2021GS}
#'  \item \href{https://wolfemd.github.io/TARI_2020GS/}{TARI_2020GS}
#'  \item \href{https://wolfemd.github.io/IITA_2020GS/}{IITA_2020GS}
#' }
#' @family cassavabase_pheno_pipeline
#' @return
#' @export
nestForMultiTrialAnalysis<-function(curatedTrialData){
  nested_trialdata<-curatedTrialData %>%
    # remove trait-trial models that failed
    filter(!is.na(H2)) %>%
    # remove some per-trial summaries we don't want at this stage
    select(-H2,-VarComps,-Model,-Noutliers,-propMiss) %>%
    unnest(BLUPs) %>%
    nest(MultiTrialTraitData=c(-Trait))
  return(nested_trialdata)
}

#' Fit model and remove outliers from a multi-trial dataset
#'
#' \strong{NOT CURRENTLY IN USE / PROBABLY DEPRECATED}
#'
#' Actually expects input to be BLUPs from single-trial mixed-models.
#' Fits mixed-models across trials and removes outliers. Returns de-regressed BLUPs and more.
#'
#' @param curatedTrialData
#' @param GID
#'
#' @details \strong{NOTICE:} This function is part of a family of functions (\code{"cassavabase_pheno_pipeline"}) developed as part of the NextGen Cassava Breeding Project genomic selection pipeline.
#' For some examples of their useage:
#' \itemize{
#'  \item \href{https://wolfemd.github.io/IITA_2021GS/}{IITA_2021GS (\strong{\emph{should be first use of this package}})}
#'  \item \href{https://wolfemd.github.io/NRCRI_2021GS/}{NRCRI_2021GS}
#'  \item \href{https://wolfemd.github.io/TARI_2020GS/}{TARI_2020GS}
#'  \item \href{https://wolfemd.github.io/IITA_2020GS/}{IITA_2020GS}
#' }
#' @family cassavabase_pheno_pipeline
#' @return
#' @export
fitMultiTrialModel<-function(curatedTrialData,GID="GID"){
  require(lme4)
  modelFormula<-paste0("drgBLUP ~ (1|",GID,")")
  fit_model<-possibly(function(modelFormula,curatedTrialData){
    model_out<-lmer(as.formula(modelFormula),
                    data=curatedTrialData,
                    weights = WT)
    return(model_out) },
    otherwise = NA)
  model_out<-fit_model(modelFormula,curatedTrialData)
  summary(model_out)
  if(is.na(model_out)){
    out <-tibble(H2=NA,VarComps=list(NULL),BLUPs=list(NULL),Model=modelFormula)
  } else {
    varcomps<-as.data.frame(VarCorr(model_out))[,c("grp","vcov")] %>%
      spread(grp,vcov)
    Vg<-varcomps$GID
    H2<-Vg/(Vg+varcomps$Residual)
    BLUP<-ranef(model_out, condVar=TRUE)[[GID]]
    PEV <- c(attr(BLUP, "postVar"))
    blups<-tibble(GID=rownames(BLUP),BLUP=BLUP$`(Intercept)`,PEV=PEV) %>%
      mutate(REL=1-(PEV/Vg),
             drgBLUP=BLUP/REL,
             WT=(1-H2)/((0.1 + (1-REL)/REL)*H2))
    out <- tibble(H2=H2,
                  VarComps=list(varcomps),
                  BLUPs=list(blups),
                  Model=modelFormula) }
  return(out)
}
