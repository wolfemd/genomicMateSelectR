#' Split a genome-wide VCF into separate chromosome-wise VCFs
#'
#' Uses \href{https://vcftools.github.io/index.html}{vcftools}
#'
#' @param Chr
#' @param vcfIn
#' @param filters
#' @param outPath
#' @param outSuffix
#'
#' @details \strong{NOTICE:} This function is part of a family of functions (\code{"imputation_functions"}) developed as part of the NextGen Cassava Breeding Project genomic selection pipeline.
#' For some examples of their useage:
#' \itemize{
#'  \item \href{https://wolfemd.github.io/IITA_2021GS/}{IITA_2021GS (\strong{\emph{should be first use of this package}})}
#'  \item \href{https://wolfemd.github.io/NRCRI_2021GS/}{NRCRI_2021GS}
#'  \item \href{https://wolfemd.github.io/TARI_2020GS/}{TARI_2020GS}
#' }
#' @family imputation_functions
#' @return
#' @export
splitVCFbyChr<-function(Chr,vcfIn,filters=NULL,outPath,outSuffix){
  system(paste0("vcftools --gzvcf ",vcfIn," ",
                "--chr ",Chr," ",filters," ",
                "--recode --stdout | bgzip -c -@ 24 > ",
                outPath,"chr",Chr,"_",outSuffix,".vcf.gz")) }

#' Run Beagle5 to impute a target VCF with a reference VCF
#'
#' Impute with \href{https://faculty.washington.edu/browning/beagle/b5_0.html}{Beagle V5.0}. Use an "imputation reference panel".
#' Refer to Beagle documentation for meaning of arguments passed.
#'
#' @param targetVCF passes to Beagle `\code{gt=targetVCF}`
#' @param refVCF passes to Beagle `\code{ref=refVCF}`
#' @param mapFile passes to Beagle `\code{map=mapFile}`
#' @param outName passes to Beagle `\code{out=outName} (\strong{don't put file suffix, Beagle adds \code{*.vcf.gz}}).
#' @param nthreads passes to Beagle `\code{nthreads=nthreads}`
#' @param maxmem passes to java `\code{-Xmx<maxmem>}`
#' @param impute passes to Beagle `\code{impute=TRUE}`
#' @param ne passes to Beagle `\code{ne=ne}`
#' @param samplesToExclude
#'
#' @details \strong{NOTICE:} This function is part of a family of functions (\code{"imputation_functions"}) developed as part of the NextGen Cassava Breeding Project genomic selection pipeline.
#'
#' \code{java -Xms2g -Xmx [maxmem] -jar /programs/beagle/beagle.jar gt= [targetVCF] map= [mapFile] ref= [refVCF] out= [outName] nthreads= [nthreads] impute= [impute]  ne=  [ne]}
#'
#' For some examples of their useage:
#' \itemize{
#'  \item \href{https://wolfemd.github.io/IITA_2021GS/}{IITA_2021GS (\strong{\emph{should be first use of this package}})}
#'  \item \href{https://wolfemd.github.io/NRCRI_2021GS/}{NRCRI_2021GS}
#'  \item \href{https://wolfemd.github.io/TARI_2020GS/}{TARI_2020GS}
#' }
#' @family imputation_functions
#' @return
#' @export
runBeagle5<-function(targetVCF,refVCF=NULL,mapFile,outName,window=40,overlap=2,
                     nthreads,maxmem="500g",impute=TRUE,ne=100000,samplesToExclude=NULL){
  system(paste0("java -Xms2g -Xmx",maxmem," -jar /programs/beagle/beagle.jar ",
                "gt=",targetVCF," ",
                "map=",mapFile," ",
                ifelse(!is.null(refVCF),paste0("ref=",refVCF," "), ""),
                "out=",outName," ",
                "nthreads=",nthreads," impute=",impute," ne=",ne, " window=",window, " overlap=",overlap,
                ifelse(!is.null(samplesToExclude),paste0(" excludesamples=",samplesToExclude),""))) }

#' Apply quality filters after imputation by Beagle5.0
#'
#' Remove markers from the imputed VCF based on adjustable thresholds.
#' For VCFs imputed by Beagle5
#'
#' Uses \href{https://vcftools.github.io/index.html}{vcftools} and R.
#' of minor allee frequency \code{MAFthresh}, p-value measuring likelihood of divergence from hardy-weinberg equilibrium ()
#'
#' @param inPath path to input VCFs
#' @param inName imputed input VCF name
#' @param outPath path for output to be written to
#' @param outName name desired, don't add file suffix, will be *.vcf.gz
#' @param DR2thresh remove if imputation quality score \code{DR2} supplied in the VCF INFO field by Beagle is less than or equal to the threshold
#' @param HWEthresh remove if p-value HWE from \code{vcftools --hardy} less than threshold, smaller p-value means more departure from HWE
#' @param MAFthresh remove if minor allele frequency less than threshold
#'
#' @details \strong{NOTICE:} This function is part of a family of functions (\code{"imputation_functions"}) developed as part of the NextGen Cassava Breeding Project genomic selection pipeline.
#' For some examples of their useage:
#' \itemize{
#'  \item \href{https://wolfemd.github.io/IITA_2021GS/}{IITA_2021GS (\strong{\emph{should be first use of this package}})}
#'  \item \href{https://wolfemd.github.io/NRCRI_2021GS/}{NRCRI_2021GS}
#'  \item \href{https://wolfemd.github.io/TARI_2020GS/}{TARI_2020GS}
#' }
#' @family imputation_functions
#' @return
#' @export
postImputeFilter<-function(inPath=NULL,inName,outPath=NULL,outName,DR2thresh=0.75,HWEthresh=1e-20,MAFthresh=0.005){
  require(magrittr); require(dplyr)
  # Extract imputation quality scores (DR2 and AF) from VCF
  system(paste0("vcftools --gzvcf ",inPath,inName,".vcf.gz --get-INFO DR2 --get-INFO AF --out ",outPath,inName))
  system(paste0("vcftools --gzvcf ",inPath,inName,".vcf.gz --hardy --out ",outPath,inName))

  # Read scores into R
  INFO<-read.table(paste0(outPath,inName,".INFO"),
                   stringsAsFactors = F, header = T)
  hwe<-read.table(paste0(outPath,inName,".hwe"),
                  stringsAsFactors = F, header = T)
  stats2filterOn<-left_join(INFO,hwe %>% rename(CHROM=CHR))
  # Compute MAF from AF and make sure numeric
  stats2filterOn %<>%
    dplyr::mutate(DR2=as.numeric(DR2),
                  AF=as.numeric(AF)) %>%
    dplyr::filter(!is.na(DR2),
                  !is.na(AF)) %>%
    dplyr::mutate(MAF=ifelse(AF>0.5,1-AF,AF))
  # Identify sites passing filter
  sitesPassingFilters<-stats2filterOn %>%
    dplyr::filter(DR2>=DR2thresh,
                  P_HWE>HWEthresh,
                  MAF>MAFthresh) %>%
    dplyr::select(CHROM,POS)
  print(paste0(nrow(sitesPassingFilters)," sites passing filter"))

  # Write a list of positions passing filter to disk
  write.table(sitesPassingFilters,
              file = paste0(outPath,inName,".sitesPassing"),
              row.names = F, col.names = F, quote = F)
  # Apply filter to vcf file with vcftools
  system(paste0("vcftools --gzvcf ",inPath,inName,".vcf.gz"," ",
                "--positions ",outPath,inName,".sitesPassing"," ",
                "--recode --stdout | bgzip -c -@ 24 > ",
                outPath,outName,".vcf.gz"))
  print(paste0("Filtering Complete: ",outName))
}

#' Merge two VCFs into one
#'
#' Uses \href{https://samtools.github.io/bcftools/bcftools.html}{\code{bcftools merge}}. Also requires \code{tabix}.
#'
#' @param inPath
#' @param inVCF1
#' @param inVCF2
#' @param outPath
#' @param outName
#'
#' @details \strong{NOTICE:} This function is part of a family of functions (\code{"imputation_functions"}) developed as part of the NextGen Cassava Breeding Project genomic selection pipeline.
#' For some examples of their useage:
#' \itemize{
#'  \item \href{https://wolfemd.github.io/IITA_2021GS/}{IITA_2021GS (\strong{\emph{should be first use of this package}})}
#'  \item \href{https://wolfemd.github.io/NRCRI_2021GS/}{NRCRI_2021GS}
#'  \item \href{https://wolfemd.github.io/TARI_2020GS/}{TARI_2020GS}
#' }
#' @family imputation_functions
#' @return
#' @export
mergeVCFs<-function(inPath=NULL,inVCF1,inVCF2,outPath=NULL,outName){
  system(paste0("tabix -f -p vcf ",inPath,inVCF1,".vcf.gz"))
  system(paste0("tabix -f -p vcf ",inPath,inVCF2,".vcf.gz"))
  system(paste0("bcftools merge ",
                "--output ",outPath,outName,".vcf.gz ",
                "--merge snps --output-type z --threads 6 ",
                inPath,inVCF1,".vcf.gz"," ",
                inPath,inVCF2,".vcf.gz"))
}

#' Subset a VCF file based on provided positions
#'
#' Uses \href{https://vcftools.github.io/index.html}{\code{vcftools --positions}}
#'
#' @param inPath
#' @param inVCF
#' @param positionFile see VCFtools
#' @param outPath
#' @param outName
#'
#' @details \strong{NOTICE:} This function is part of a family of functions (\code{"imputation_functions"}) developed as part of the NextGen Cassava Breeding Project genomic selection pipeline.
#' For some examples of their useage:
#' \itemize{
#'  \item \href{https://wolfemd.github.io/IITA_2021GS/}{IITA_2021GS (\strong{\emph{should be first use of this package}})}
#'  \item \href{https://wolfemd.github.io/NRCRI_2021GS/}{NRCRI_2021GS}
#'  \item \href{https://wolfemd.github.io/TARI_2020GS/}{TARI_2020GS}
#' }
#' @family imputation_functions
#' @return
#' @export
filter_positions<-function(inPath=NULL,inVCF,positionFile,outPath=NULL,outName){
  system(paste0("vcftools --gzvcf ",inPath,inVCF," ",
                "--positions ",inPath,positionFile," ",
                "--recode --stdout | bgzip -c -@ 24 > ",
                outPath,outName,".vcf.gz"))
}

#' Convert a VCF file to dosage matrix
#'
#' The function below will
#'
#' (1) convert the input VCF to plink1.9 binary format and
#' (2) convert the plink binary to a dosage (0,1,2) matrix with special attention to which allele gets counted in the file.
#'
#'
#' Uses \href{https://www.cog-genomics.org/plink/}{\code{plink1.9}}. With \code{plink1.9}
#' there is some risk the counted allele could switch between
#' e.g. the reference panel and the progeny files because of allele freq. (see plink documentation).
#' To avoid this, went to extra trouble: write a file suffixed \code{*.alleleToCount}
#' listing SNP ID (column 1) and the ALT allele from the VCF (column 2).
#' Pass the file to \code{plink1.9} using the \code{--recode-allele} flag
#' to ensure all output dosages count the ALT allele consistent with the VCFs.
#' The reason to use \code{plink1.9} is that \code{Beagle5} imputed files
#' don't have a \strong{DS} (dosage) field that can be directly extracted.
#' Instead, phased genotypes e.g. \code{0|1} need to be converted to dosages
#' (e.g. \code{0|1 --> 1}, code{1|1 --> 2}).
#' An alternative might be to extract the haplotypes using code{vcftools} and
#' manually compute the dosages.

#' @param pathIn
#' @param pathOut
#' @param vcfName
#'
#' @details \strong{NOTICE:} This function is part of a family of functions (\code{"imputation_functions"}) developed as part of the NextGen Cassava Breeding Project genomic selection pipeline.
#' For some examples of their useage:
#' \itemize{
#'  \item \href{https://wolfemd.github.io/IITA_2021GS/}{IITA_2021GS (\strong{\emph{should be first use of this package}})}
#'  \item \href{https://wolfemd.github.io/NRCRI_2021GS/}{NRCRI_2021GS}
#'  \item \href{https://wolfemd.github.io/TARI_2020GS/}{TARI_2020GS}
#' }
#' @family imputation_functions
#' @return
#' @export
convertVCFtoDosage<-function(pathIn,pathOut,vcfName){
  #' @vcfName expects file to be *.vcf.gz, expects name _not_ to include extension .vcf.gz
  # Make binary plink
  system(paste0("export PATH=/programs/plink-1.9-x86_64-beta3.30:$PATH;",
                "plink --vcf ",pathIn,vcfName,".vcf.gz ",
                "--make-bed --const-fid --keep-allele-order ",
                "--out ",pathOut,vcfName))
  # Write first 5 columns of VCF file to txt file *.sitesWithAlleles
  system(paste0("zcat ",pathIn,vcfName,".vcf.gz ",
                "| cut -f1-5 > ",pathOut,vcfName,".sitesWithAlleles"))
  # Read *.sitesWithAlleles, take columns required by plink and write back to disk
  read.table(paste0(pathOut,vcfName,".sitesWithAlleles"),
             stringsAsFactors = F, header = F, sep = c("\t")) %>%
    select(V3,V5) %>%
    write.table(.,file=paste0(pathOut,vcfName,".alleleToCount"),
                row.names = F, sep = c("\t"), quote = F, col.names = F)
  # Recode to dosage
  system(paste0("export PATH=/programs/plink-1.9-x86_64-beta3.30:$PATH;",
                "plink --bfile ",pathOut,vcfName," --keep-allele-order --recode A ",
                "--recode-allele ",pathOut,vcfName,".alleleToCount ",
                "--out ",pathOut,vcfName))
}

#' Combine multiple dosage matrix to a single genome-wide dosage
#'
#' Inputs are VCF files, one for each chromosome presumably.
#'
#' Use common naming for VCF with prefix e.g. "chr12_".
#'
#' Combines "<chroms>_<nameOfchromWiseDosageFiles>".
#'
#' Output is an .rds (R dataset) format.
#'
#' @param pathIn
#' @param chroms
#' @param nameOfchromWiseDosageFiles
#'
#' @details \strong{NOTICE:} This function is part of a family of functions (\code{"imputation_functions"}) developed as part of the NextGen Cassava Breeding Project genomic selection pipeline.
#' For some examples of their useage:
#' \itemize{
#'  \item \href{https://wolfemd.github.io/IITA_2021GS/}{IITA_2021GS (\strong{\emph{should be first use of this package}})}
#'  \item \href{https://wolfemd.github.io/NRCRI_2021GS/}{NRCRI_2021GS}
#'  \item \href{https://wolfemd.github.io/TARI_2020GS/}{TARI_2020GS}
#' }
#' @family imputation_functions
#' @return Output is an .rds (R dataset) format.
#' @export
createGenomewideDosage<-function(pathIn,chroms,nameOfchromWiseDosageFiles){
  snps<-future_map(chroms,~read.table(paste0(pathIn,"chr",.,nameOfchromWiseDosageFiles,".raw"), stringsAsFactor=F, header = T) %>%
                     dplyr::select(-FID,-PAT,-MAT,-SEX,-PHENOTYPE) %>%
                     column_to_rownames(var = "IID") %>%
                     as.matrix()) %>%
    reduce(.,cbind)
  saveRDS(snps,file = paste0(pathIn,"DosageMatrix",nameOfchromWiseDosageFiles,".rds"))
}

#' Run Beagle4.1 to impute a target VCF with a reference VCF
#'
#' \href{https://faculty.washington.edu/browning/beagle/b4_1.html}{Beagle4.1} is
#' slower than Beagle5 by far. However, it can use genotype-likelihoods (the GL
#' VCF field) for its algorithm, which is potentially more accurate. In cases
#' like imputing GBS data, or if you have time to wait around, this should be
#' better than Beagle5.0 for a first pass imputation of observed sites. Impute
#' with \href{https://faculty.washington.edu/browning/beagle/b5_0.html}{Beagle
#' V5.0}. Use an "imputation reference panel". Refer to Beagle documentation for
#' meaning of arguments passed.
#'
#' @param targetVCF passes to Beagle `\code{gl=targetVCF}`
#' @param refVCF passes to Beagle `\code{ref=refVCF}`
#' @param mapFile passes to Beagle `\code{map=mapFile}`
#' @param outName passes to Beagle `\code{out=outName} (\strong{don't put file
#'   suffix, Beagle adds \code{*.vcf.gz}}).
#' @param nthreads passes to Beagle `\code{nthreads=nthreads}`
#' @param maxmem passes to java `\code{-Xmx<maxmem>}`
#' @param impute passes to Beagle `\code{impute=TRUE}`
#' @param ne passes to Beagle `\code{ne=ne}`
#' @param samplesToExclude
#' @param niter passes to Beagle `\code{niterations=niter}`
#'
#' @details \strong{NOTICE:} This function is part of a family of functions (\code{"imputation_functions"}) developed as part of the NextGen Cassava Breeding Project genomic selection pipeline.
#'
#' \code{java -Xms2g -Xmx [maxmem] -jar /programs/beagle41/beagle41.jar gl= [targetVCF] map= [mapFile] ref= [refVCF] out= [outName] niterations= [niter] nthreads= [nthreads] impute= [impute]  ne=  [ne]}
#'
#' For some examples of their useage:
#' \itemize{
#'  \item \href{https://wolfemd.github.io/IITA_2021GS/}{IITA_2021GS (\strong{\emph{should be first use of this package}})}
#'  \item \href{https://wolfemd.github.io/NRCRI_2021GS/}{NRCRI_2021GS}
#'  \item \href{https://wolfemd.github.io/TARI_2020GS/}{TARI_2020GS}
#' }
#' @family imputation_functions
#' @return
#' @export
runBeagle4pt1GL<-function(targetVCF,refVCF,mapFile,outName,
                          nthreads,maxmem="500g",impute=TRUE,ne=100000,samplesToExclude=NULL,niter=10){
  system(paste0("java -Xms2g -Xmx",maxmem," -jar /programs/beagle41/beagle41.jar ",
                "gl=",targetVCF," ",
                "map=",mapFile," ",
                ifelse(!is.null(refVCF),paste0("ref=",refVCF," "), ""),
                "out=",outName," ",
                "niterations=",niter," ",
                "nthreads=",nthreads," impute=",impute," ne=",ne,
                ifelse(!is.null(samplesToExclude),paste0(" excludesamples=",samplesToExclude),""))) }

#' Apply quality filters after imputation by Beagle4.1
#'
#' Remove markers from the imputed VCF based on adjustable thresholds.
#' For VCFs imputed by Beagle4.1.
#'
#' Uses \href{https://vcftools.github.io/index.html}{vcftools} and R.
#' of minor allee frequency \code{MAFthresh}, p-value measuring likelihood of divergence from hardy-weinberg equilibrium ()
#'
#' @param inPath path to input VCFs
#' @param inName imputed input VCF name
#' @param outPath path for output to be written to
#' @param outName name desired, don't add file suffix, will be *.vcf.gz
#' @param DR2thresh remove if imputation quality score \code{DR2} supplied in the VCF INFO field by Beagle is less than or equal to the threshold
#' @param HWEthresh remove if p-value HWE from \code{vcftools --hardy} less than threshold, smaller p-value means more departure from HWE
#' @param MAFthresh remove if minor allele frequency less than threshold
#'
#' @details \strong{NOTICE:} This function is part of a family of functions (\code{"imputation_functions"}) developed as part of the NextGen Cassava Breeding Project genomic selection pipeline.
#' For some examples of their useage:
#' \itemize{
#'  \item \href{https://wolfemd.github.io/IITA_2021GS/}{IITA_2021GS (\strong{\emph{should be first use of this package}})}
#'  \item \href{https://wolfemd.github.io/NRCRI_2021GS/}{NRCRI_2021GS}
#'  \item \href{https://wolfemd.github.io/TARI_2020GS/}{TARI_2020GS}
#' }
#' @family imputation_functions
#' @return
#' @export
postImputeFilterBeagle4pt1<-function(inPath=NULL,inName,outPath=NULL,outName,AR2thresh=0.75,HWEthresh=1e-20,MAFthresh=0.005){
  require(magrittr); require(dplyr)
  # Extract imputation quality scores (DR2 and AF) from VCF
  system(paste0("vcftools --gzvcf ",inPath,inName,".vcf.gz --get-INFO AR2 --get-INFO AF --out ",outPath,inName))
  system(paste0("vcftools --gzvcf ",inPath,inName,".vcf.gz --hardy --out ",outPath,inName))

  # Read scores into R
  INFO<-read.table(paste0(outPath,inName,".INFO"),
                   stringsAsFactors = F, header = T)
  hwe<-read.table(paste0(outPath,inName,".hwe"),
                  stringsAsFactors = F, header = T)
  stats2filterOn<-left_join(INFO,hwe %>% rename(CHROM=CHR))
  # Compute MAF from AF and make sure numeric
  stats2filterOn %<>%
    dplyr::mutate(AR2=as.numeric(AR2),
                  AF=as.numeric(AF)) %>%
    dplyr::filter(!is.na(AR2),
                  !is.na(AF)) %>%
    dplyr::mutate(MAF=ifelse(AF>0.5,1-AF,AF))
  # Identify sites passing filter
  sitesPassingFilters<-stats2filterOn %>%
    dplyr::filter(AR2>=AR2thresh,
                  P_HWE>HWEthresh,
                  MAF>MAFthresh) %>%
    dplyr::select(CHROM,POS)
  print(paste0(nrow(sitesPassingFilters)," sites passing filter"))

  # Write a list of positions passing filter to disk
  write.table(sitesPassingFilters,
              file = paste0(outPath,inName,".sitesPassing"),
              row.names = F, col.names = F, quote = F)
  # Apply filter to vcf file with vcftools
  system(paste0("vcftools --gzvcf ",inPath,inName,".vcf.gz"," ",
                "--positions ",outPath,inName,".sitesPassing"," ",
                "--recode --stdout | bgzip -c -@ 24 > ",
                outPath,outName,".vcf.gz"))
  print(paste0("Filtering Complete: ",outName))
}
