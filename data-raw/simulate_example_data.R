## code to prepare `simulate_example_data` dataset goes here
# usethis::use_data(simulate_example_data, overwrite = TRUE)

library(tidyverse); library(magrittr)
library(AlphaSimR)
library(genomicMateSelectR)

rm(list=ls())
# SET-UP FOUNDER POPULATION
## simulate founder haplotypes
founderHap <- runMacs2(nInd=100,nChr=2,segSites=100)

# Specify some simulation parameters
SP <- SimParam$new(founderHap)
SP$restrSegSites(minQtlPerChr=10,minSnpPerChr=50,overlap=FALSE)
SP$setSexes("no") #all individuals are hermaphrodites
SP$addSnpChip(nSnpPerChr=50) # Observed SNPs per chromosome
SP$setTrackPed(TRUE) #keeps pedigree information in slot SP@pedigree
SP$setTrackRec(TRUE) #keeps recomb. records of all individuals in slot of "SP"

# Set-up 2 genetically correlated traits
## with both additive and dominance effects
G = 1.5*diag(2)-0.5 #Genetic correlation matrix
SP$addTraitAD(nQtlPerChr = 10,
              mean=c(0,0), meanDD=c(0,0),
              var=c(1,1), varDD = c(1,1),
              corA=G, corDD=G)

# Now create an initial population, a set of founders
founders <- newPop(founderHap,simParam=SP)

# Initial founder phenotypes
founders <- setPheno(pop=founders,h2=0.5,reps=2)

# Simulate 2 cylces of selection
## generates a small pedigree and phenotype data to analyze in the demo.
nCycles<-2
# very simple container for each cycles sim output
simOutput<-list(founders)

for(cycle in 1:nCycles){
  # choose the best from last cycle
  ## select based only on Trait1, trait2 is along for the ride
  chosenParents<- selectInd(pop=simOutput[[cycle]],nInd=5,use="pheno")
  # make crosses
  offspringPop<-randCross(pop=chosenParents,
                          nCrosses=10, nProgeny = 10)
  # phenotype  new offspring
  offspringPop<-setPheno(pop = offspringPop, h2=0.5, reps=2)
  # add new offspring to simOutput list
  simOutput[[cycle+1]]<-offspringPop
}

# Tidy up to output
tidySimOutput<-tibble(Cycle=0:nCycles,
                      Sims=simOutput) %>%
  dplyr::mutate(meanG=map(Sims,~colMeans(.@gv)),
                varG=map(Sims,~var(.@gv)))

save(tidySimOutput,SP,
     file = here::here("data-raw","tidySimOutput_and_SP.Rdata"))

# Create inputs for genomicMateSelectR

## Pedigree
ped<-SP$pedigree %>%
  as.data.frame %>%
  dplyr::mutate(GID=rownames(.)) %>%
  as_tibble %>%
  dplyr::select(GID,father,mother)

## Phenotypes
phenos<-tidySimOutput %>%
  dplyr::select(Cycle,Sims) %>%
  dplyr::mutate(phenos=map(Sims,function(Sims){
    AlphaSimR::pheno(Sims) %>%
      as_tibble(.) %>%
      dplyr::mutate(GID=Sims@id) %>%
      dplyr::relocate(GID,.before = V1) })) %>%
  dplyr::select(-Sims) %>%
  unnest(phenos) %>%
  rename(Trait1=V1,Trait2=V2)

## Tidy (ready) phenotypes
### will label them 'blups' to match
#### argument's name in e.g. runGenomicPredictions() function
#### even though they are not BLUPs
#### currently, genomicMateSelectR set-up to integrate into a two-stage
#### genomic evaluation with 2nd stage taking de-regressed BLUPs and
#### weighting error variances
### for this example dataset, simply set all weights to be ==1
blups<-phenos %>%
  tidyr::pivot_longer(cols = contains("Trait"),
                      names_to = "Trait",
                      values_to = "drgBLUP") %>%
  dplyr::mutate(BLUP=drgBLUP,
                WT=1) %>%
  dplyr::select(-Cycle) %>%
  nest(TrainingData=-Trait)


# centiMorgan scale genetic map
## get the genetic map
genmap<-getSnpMap(simParam = SP)
m<-genmap$pos*100; # convert it to centimorgans
names(m)<-genmap$id ## as a named vector
genmap<-m; rm(m)

# construct the recombination frequency matrix
## actually 1-2*recomb. freq. matrix
recombFreqMat<-1-(2*genomicMateSelectR::genmap2recombfreq(genmap,nChr = SP$nChr))
# image(genomicMateSelectR::genmap2recombfreq(m,nChr = SP$nChr))
# image(recombFreqMat)


# haplotype matrix
haploMat<-pullSnpHaplo(mergePops(tidySimOutput$Sims))
# change haplotype tags in rownames of the haploMat to please genomicMateSelectR
rownames(haploMat) %<>%
  gsub("_1","_HapA",.) %>%
  gsub("_2","_HapB",.)

# dosage matrix
doseMat<-pullSnpGeno(mergePops(tidySimOutput$Sims))

# kinships
A<-genomicMateSelectR::kinship(doseMat,"add")
Dgeno<-genomicMateSelectR::kinship(doseMat,"domGenotypic")
Dclassic<-genomicMateSelectR::kinship(doseMat,"domClassic")
grms<-list(A=A,Dgeno=Dgeno,Dclassic=Dclassic)
SIwts<- c(0.5,0.5) %>% `names<-`(.,blups$Trait)
# runGenomicPredictions and get marker effects for downstream use
gpredsA<-runGenomicPredictions(modelType = "A",
                               selInd = T,
                               SIwts = SIwts,
                               getMarkEffs = T,returnPEV = F,
                               blups = blups,dosages = doseMat,
                               grms = list(A=grms[["A"]]))

gpredsAD<-runGenomicPredictions(modelType = "AD",
                                selInd = T,
                                SIwts = SIwts,
                                getMarkEffs = T,returnPEV = F,
                                blups = blups,dosages = doseMat,
                                grms = list(A=grms[["A"]],
                                            D=grms[["Dclassic"]]))
gpredsDirDom<-runGenomicPredictions(modelType = "AD",
                                    selInd = T,
                                    SIwts = SIwts,
                                    getMarkEffs = T,returnPEV = F,
                                    blups = blups,dosages = doseMat,
                                    grms = list(A=grms[["A"]],
                                                D=grms[["Dgeno"]]))

snpeffsA<-gpredsA$genomicPredOut[[1]]
snpeffsAD<-gpredsAD$genomicPredOut[[1]]
snpeffsDirDom<-gpredsDirDom$genomicPredOut[[1]]

# that's pretty much all the file types one might need

# store the data we want
usethis::use_data(ped, overwrite = TRUE)
usethis::use_data(phenos, overwrite = TRUE)
usethis::use_data(blups, overwrite = TRUE)
usethis::use_data(genmap, overwrite = TRUE)
usethis::use_data(recombFreqMat, overwrite = TRUE)
usethis::use_data(haploMat, overwrite = TRUE)
usethis::use_data(doseMat, overwrite = TRUE)
usethis::use_data(grms, overwrite = TRUE)
usethis::use_data(snpeffsA, overwrite = TRUE)
usethis::use_data(snpeffsAD, overwrite = TRUE)
usethis::use_data(snpeffsDirDom, overwrite = TRUE)
