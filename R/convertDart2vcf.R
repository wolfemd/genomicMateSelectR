#' Convert DArTseqLD reports to VCF
#'
#' This function converts a dual-file report format that we (NextGen Cassava Breeding) have arranged with DArT.
#'
#' The output of this step is written to disk as a \code{*.vcf.gz} completely compatible with downstream bioinformatics softwares!
#'
#' @dartvcfInput input name and path of "vcf" file from DArT
#' @dartcountsInput input name and path of counts file from DArT
#' @outName output path and name
#' @nskipvcf number of "VCF" rows to skip on read-in
#' @nskipcounts number of "counts file" rows to skip on read in
#' @ncores number of cores to use, could be VERY memory intensive
#'
#' @details \strong{NOTICE:} This function is part of a family of functions (\code{"imputation_functions"}) developed as part of the NextGen Cassava Breeding Project genomic selection pipeline.
#' For some examples of their useage:
#' \itemize{
#'  \item \href{https://wolfemd.github.io/IITA_2021GS/}{IITA_2021GS (\strong{\emph{should be first use of this package}})}
#'  \item \href{https://wolfemd.github.io/NRCRI_2021GS/}{NRCRI_2021GS}
#'  \item \href{https://wolfemd.github.io/TARI_2020GS/}{TARI_2020GS}
#' }
#' @family imputation_functions
#' @return The output of this step is written to disk as a \code{*.vcf.gz} completely compatible with downstream bioinformatics softwares!
#' @export
convertDart2vcf<-function(dartvcfInput,dartcountsInput,outName,
                          nskipvcf=2,nskipcounts=3,ncores){
  require(tidyverse); require(magrittr)
  # Read and format the counts/vcf files  -------------------

  ## Read in "VCF" (called genotypes from DArT)
  ## and "counts" files (read counts)
  vcf<-read.table(dartvcfInput,
                  stringsAsFactors = F,skip = nskipvcf, header = T, sep = "\t", comment.char = "")
  readCounts<-read.csv(dartcountsInput, stringsAsFactors = F,header = T,skip=nskipcounts)

  ## Formatting chrom and position info
  vcf %<>%
    mutate(X.CHROM=gsub("Chromosome","",X.CHROM)) %>%
    rename(Pos=POS) %>%
    mutate(Chr=as.numeric(gsub("Chromosome","",X.CHROM)),
           Pos=as.numeric(Pos))

  readCounts %<>%
    mutate(RefAlt=ifelse(SNP=="","REF","ALT"),
           Chr=as.numeric(gsub("Chromosome","",Chrom_Cassava_v61)),
           Pos=as.numeric(SNP_ChromPos_Cassava_v61))
  # Add a SNPindex ------------------
  # Add a unique value "SNPindex" to each SNP in the vcf and readCounts df's
  # For readCounts, this is going to be the best way to keep track of pairs of rows. Since in many cases, multiple CloneID can point to same Chr-Pos-Alleles and it's otherwise unclear which pair of rows should go together when subsetting downstream.
  vcf %<>%
    mutate(SNPindex=1:nrow(.))
  readCounts %<>%
    mutate(SNPindex=sort(rep(1:(nrow(.)/2),2)))

  readCounts %<>%
    separate(AlleleID,c("tmp","Alleles"),sep=-3,remove = F) %>%
    select(-tmp) %>%
    separate(Alleles,c("REF","ALT"),">",remove = F)
  vcf %<>%
    separate(ID,c("tmp","Alleles"),sep=-3,remove = F) %>%
    select(-tmp)
  vcf %<>%
    separate(ID,c("CloneID","tmp"),"[|]",remove = F,extra = 'merge') %>%
    mutate(CloneID=as.numeric(CloneID)) %>%
    select(-tmp) %>%
    rename(AlleleID=ID)

  # Add required VCF fields --------------------------
  ## First have to do some data transformation and
  ## create some of the meta-data fields in a VCF, e.g. QUAL, FILTER INFO.
  readCounts %<>%
    arrange(Chr,Pos,SNPindex,CloneID,AlleleID,RefAlt) %>%
    mutate(QUAL=".",
           FILTER=".",
           INFO=".",
           FORMAT="GT:AD:DP:PL",
           SNP_ID=paste0("S",Chr,"_",Pos))
  vcf %<>%
    arrange(Chr,Pos,SNPindex,CloneID,AlleleID) %>%
    mutate(QUAL=".",
           FILTER=".",
           INFO=".",
           FORMAT="GT:AD:DP:PL",
           SNP_ID=paste0("S",Chr,"_",Pos))

  # Prune duplicate and multi-allelic sites ---------------
  sites2keep<-vcf %>%
    count(Chr,Pos) %>%
    ungroup() %>%
    filter(n==1) %>%
    select(-n) %>%
    semi_join(
      readCounts %>%
        count(Chr,Pos) %>%
        arrange(desc(n)) %>%
        filter(n==2) %>%
        select(-n))
  vcf %<>% semi_join(sites2keep)
  readCounts %<>% semi_join(sites2keep)

  sampleIDsFromDartVCF<-colnames(vcf) %>%
    .[!. %in% c("X.CHROM","Pos","AlleleID","Alleles","REF","ALT","QUAL","FILTER","INFO","FORMAT",
                "Chr","SNPindex","CloneID","SNP_ID","SNP")]

  tmp_counts <-readCounts %>%
    .[,c("Chr","Pos","SNPindex","SNP_ID",
         "CloneID","AlleleID","RefAlt","Alleles","QUAL","FILTER","INFO","FORMAT",sampleIDsFromDartVCF)] %>%
    semi_join(sites2keep)

  # readCounts and vcf long by sample ------------------------------
  tmp_counts %<>%
    pivot_longer(cols = all_of(sampleIDsFromDartVCF),
                 names_to = "FullSampleName",
                 values_to = "ReadCount")
  #gather(FullSampleName,ReadCount,sampleIDsFromDartVCF)
  tmp_counts %<>%
    select(-AlleleID) %>%
    spread(RefAlt,ReadCount)
  tmp_counts %<>%
    rename(AltCount=ALT,
           RefCount=REF)
  vcf_long<-vcf %>%
    .[,c("Chr","Pos","SNPindex","SNP_ID",
         "CloneID","Alleles","REF","ALT","QUAL","FILTER","INFO","FORMAT",sampleIDsFromDartVCF)] %>%
    pivot_longer(cols = all_of(sampleIDsFromDartVCF), names_to = "FullSampleName", values_to = "GT")
  #    gather(FullSampleName,GT,sampleIDsFromDartVCF)

  # Add GT from vcf to the counts -----------------------------------
  tmp_counts %<>%
    inner_join(vcf_long)

  # Calc PL field  --------------------------------------------------
  # AD+DP fields
  ## Now we can calc DP and formate the VCF field "AD" (e.g. "21,0" for 21 reference reads and 0 alt. allele reads)
  tmp_counts %<>%
    mutate(DP=AltCount+RefCount,
           AD=paste0(RefCount,",",AltCount))
  tmp1<-tmp_counts %>%
    filter(!is.na(AltCount),
           !is.na(RefCount))
  tmp<-tmp1; rm(tmp1)

  # Calc. genotype likelihoods ----------------------------------------

  ## Genotype likelihoods calculated according to:
  ### http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000862#s4
  ## Converted to Normalized Phred Scores according to:
  ### https://gatkforums.broadinstitute.org/gatk/discussion/5913/math-notes-how-pl-is-calculated-in-haplotypecaller
  ## Truncate low Phred probabilities (high Phred scores) to
  ### 255 max according to TASSEL's convention (Jeff Glaubitz, pers. communication).

  #ref<-171; alt<-171; error<-0.001
  calcPL<-function(ref,alt,error=0.001){
    # ref and alt arguments are read counts for ref and alt allele, repsectively
    dp<-ref+alt
    # for values >170, factorial() returns 'inf'
    # Since it means essentially 100% probability of a genotype...
    # set DP to 169 cieling, ref/alt to equiv. allele proportions
    if(dp>=170){ ref<-169*(ref/dp); alt<-169*(alt/dp); dp<-169 }
    gl_RefRef<-(factorial(dp)/(factorial(ref)*factorial(alt)))*(1-(0.75*error))^ref*(error/4)^(alt)
    gl_RefAlt<-(factorial(dp)/(factorial(ref)*factorial(alt)))*(0.5-(0.25*error))^(ref+alt)*(error/4)^(0)
    gl_AltAlt<-(factorial(dp)/(factorial(ref)*factorial(alt)))*(1-(0.75*error))^alt*(error/4)^(ref)
    phredScale<--10*log10(c(gl_RefRef,gl_RefAlt,gl_AltAlt))
    minPhred<-min(phredScale)
    normPhred<-round(phredScale-minPhred,0)
    normPhred[which(normPhred>=255)]<-255
    normPhred<-paste0(normPhred,collapse = ",")
    if(dp==0){ normPhred<-"." }
    return(normPhred)
  }
  require(furrr); plan(multisession, workers = ncores)
  options(future.globals.maxSize=+Inf); options(future.rng.onMisuse="ignore")
  tmp %<>%
    mutate(PL=future_map2_chr(RefCount,AltCount,~calcPL(ref=.x,alt=.y)))
  plan(sequential)

  # Final VCF format ----------------------------------------
  tmp %<>%
    mutate(FORMATfields=paste(GT,AD,DP,PL,sep=":")) %>%
    select(Chr,Pos,SNP_ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,FullSampleName,FORMATfields) %>%
    spread(FullSampleName,FORMATfields) %>%
    arrange(Chr,Pos) %>%
    rename(`#CHROM`=Chr,
           POS=Pos,ID=SNP_ID)

  # Header ----------------------------------------
  header<-c("##fileformat=VCFv4.0",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
            "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the reference and alternate alleles in the order listed\">",
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth (only filtered reads used for calling)\">",
            "##FORMAT=<ID=PL,Number=3,Type=Float,Description=\"Normalized, Phred-scaled likelihoods for AA,AB,BB genotypes where A=ref and B=alt; not applicable if site is not biallelic\">")

  # Write to disk ----------------------------------------

  options("scipen"=1000, "digits"=4)
  # for a few SNPs, position kept printing in sci notation e.g. 1e3, screws up Beagle etc., this avoids that (I hope)
  write_lines(header,
              path=paste0(outName,".vcf"))
  write.table(tmp,
              paste0(outName,".vcf"),
              append = T,sep = "\t",row.names=F, col.names=T, quote=F)
  # Save sitesWithAlleles
  tmp %>%
    rename(CHROM=`#CHROM`) %>%
    select(CHROM:ALT) %>%
    write.table(.,file=paste0(outName,".sitesWithAlleles"),
                row.names=F)
  # Save sample list
  write.table(sampleIDsFromDartVCF,file=paste0(outName,".samples"),
              row.names = F, col.names = F, quote = F)

  # BGzip
  system(paste0("cat ",outName,".vcf ",
                "| bgzip -c > ",outName,".vcf.gz"))
  system(paste0("rm ",outName,".vcf"))
}
