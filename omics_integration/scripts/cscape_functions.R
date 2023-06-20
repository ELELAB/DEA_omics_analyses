#' This function lifts a MAF file to a different genomic build.
#' @param Infile A tibble of MAF. 
#' @param Current_Build A charcter string, either \code{GRCh38} or \code{GRCh37}.
#' @import dplyr  
#' @importFrom magrittr "%>%"
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom rtracklayer import.chain liftOver
#' @return MAF tibble with positions lifted to another build 
#' @export
#' @examples
#' 
#' LiftMAF(Infile, Current_Build = 'GRCh38')


LiftMAF <- function(Infile, Current_Build){
  #The input file is assumed to be maf_tibble, this file is then lifted to 
  # either 38 or 37 and return as a tibble
  flag <- FALSE
  if(Current_Build == 'GRCh38'){
    chainBuild <- "hg38ToHg19.over.chain"
    #Make chain from 38 to 19
    path = system.file(package="liftOver", "extdata", chainBuild)
    flag <- TRUE
  } else if(Current_Build == "GRCh37"){
    path <- "../data/hg19ToHg38.over.chain"
    flag <- TRUE
  } else {
    print("Error: Build must be either GRCh38 or GRCh37")
  }
  
  if(flag == TRUE){
    #Change to Grange format
    infile_GRange <- makeGRangesFromDataFrame(Infile, 
                                              start.field = "Start_Position", 
                                              end.field = "End_Position",
                                              seqnames.field = "Chromosome",
                                              keep.extra.columns = TRUE)
    #Import chain
    chain <- import.chain(path)
    
    #Do liftover
    infile_GRange_lifted <- liftOver(x = infile_GRange, chain = chain)
    
    #recreate maf column names and tibble format
    outfile_tibble_lifted <- as_tibble(infile_GRange_lifted) %>% 
      dplyr::rename(Start_Position = start,
                    End_Position = end, 
                    Chromosome = seqnames,
                    Strand = strand) %>% 
      dplyr::select(-c(group, group_name, width))
    
    return(outfile_tibble_lifted)
  }
}




#' This function extracts columns from a MAF tibble to fit CScape input format
#' @param MAF tibble of MAF 
#' @import dplyr  
#' @importFrom magrittr "%>%"
#' @importFrom tidyr separate
#' @return tibble of cscape-somatic input
#' @export
#' @examples
#' 
#' MAFtoCscape(MAF)


MAFtoCscape <- function(MAF){
  cscape <- MAF %>% 
    filter(Variant_Type == 'SNP') %>% 
    dplyr::select(Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele1, Tumor_Seq_Allele2) %>% 
    mutate(Mutant = case_when(Reference_Allele == Tumor_Seq_Allele1 ~ Tumor_Seq_Allele2,
                              Reference_Allele == Tumor_Seq_Allele2 ~ Tumor_Seq_Allele1)) %>% 
    separate(Chromosome, into = c(NA, "Chr"), sep = 3) %>% 
    dplyr::select(Chr, Start_Position, Reference_Allele, Mutant)
  return(cscape)
}





#' This function retrives the individial score for a SNP
#' @param Ranges The position
#' @param Reference_Allele The reference nucleotide
#' @param Mutant The mutant nucleotide
#' @param file_coding cscape_table with coding scores
#' @param file_noncoding cscape_table with noncoding scores
#' @import dplyr  
#' @importFrom magrittr "%>%"
#' @importFrom tidyr unite nest unnest
#' @importFrom stringr str_split str_replace_all
#' @importFrom seqminer tabix.read.table
#' @importFrom tibble tibble
#'
#' @return returns the score
#' @export
#' @examples
#' 
#' data <- tabix_func(Ranges, Reference_Allele, Mutant, file_coding, file_noncoding)



tabix_func <- function(Ranges, Reference_Allele, Mutant, file_coding, file_noncoding){
  
  flag <- FALSE
  bases <- c("A", "T","G", "C")
  chromosomes <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13",
                   "14", "15", "16", "17", "18", "19", "20", "21", "22") #X and Y not avail
  chromosome <- str_split(Ranges, pattern = ':',simplify = TRUE) %>% .[[1]]
  remark <- 'std'
  type <- ''
  
  reference <- as.character(Reference_Allele)
  mutant <- as.character(Mutant)
  # Check Reference, Mutant and Chromosome are correct format and possible 
  if( !(reference %in% bases) | !(mutant %in% bases)){
    remark <- "Error: Reference and mutant nucleotides must both be ACGT"
    score <- NA
  } else if( !(chromosome %in% chromosomes)){
    remark <- "Error: Unexpected chromosome (X and Y not available)"
    score <- NA
    # If all okay, continue:
  } else{
    #Look for score in coding region
    x <- as_tibble(tabix.read.table(tabixFile = file_coding,
                                    tabixRange = Ranges)) %>%
      mutate(across(where(is.logical),as.character)) %>%
      mutate(across(.cols = everything(),
                    .fns =~ str_replace_all(string =., pattern = "TRUE", replacement = "T")))
    
    # Has the data been found
    if(dim(x)[1] != 0){
      flag <- TRUE
      type <- 'Coding'
    } else {
      # If no annotation is found in this position try the noncoding file
      x <- as_tibble(tabix.read.table(tabixFile = file_noncoding,
                                      tabixRange = Ranges)) %>%
        mutate(across(where(is.logical),as.character)) %>%
        mutate(across(.cols = everything(),
                      .fns =~ str_replace_all(string =., pattern = "TRUE", replacement = "T")))
      
      
      if(dim(x)[1] != 0){
        flag <- TRUE
        type <- 'Noncoding'
      } else{
        # if there still is no data found at this position
        flag <- FALSE
        remark <- 'Error: Unexcepted position'
        score <- NA
        print('here 3')
      }
    }
    
    # If the position has been found
    if(flag == TRUE){
      if (reference == mutant){
        remark <- 'Error: Reference and mutant must be different'
        score <- NA
      }
      else if (reference == x$V3[1]) {
        score <- x %>% filter(V4 == mutant) %>% dplyr::select(V5) %>% pull()
        remark <- confidence(score, type)
      } else {
        remark <- 'Error: Unexpected reference'
        score <- NA
      }
    }}
  
  if (type == 'Coding'){
    data <- tibble(Coding_score = score, Remark = remark)
  } else if (type == 'Noncoding'){
    data <- tibble(Noncoding_score = score, Remark = remark)
  } else {
    data <- tibble(Remark = remark) #tibble(Coding_score = '', Noncoding_score = '', Remark = remark)
  }
  return(data)
}






#' This function retrive cscape-scores to SNPs
#' @param input Input matching cscape input
#' @param coding_file cscape_table with coding scores
#' @param noncoding_file cscape_table with noncoding scores
#' @import dplyr  
#' @importFrom magrittr "%>%"
#' @importFrom tidyr unite nest unnest
#' @importFrom purrr pmap
#' @importFrom readr parse_guess
#'
#' @return returns a tibble with a score and remark for each SNP
#' @export
#' @examples
#' 
#' cscape_out <- RunCscape_somatic(input, coding_file, noncoding_file)


RunCscape_somatic <- function(input, coding_file, noncoding_file){
  cscape_in <- input %>%
    unite(col = Ranges, c("Chr", "Start_Position"), sep = ":", remove = FALSE) %>%
    unite(col = Ranges, c("Ranges","Start_Position"), sep = "-", remove = FALSE)
  
  #When using Pmap important to call the list and varibales the same name
  cscape_annot <- cscape_in %>% 
    mutate(file_coding = coding_file,
           file_noncoding = noncoding_file) %>%
    dplyr::select(-Start_Position) %>% 
    mutate(Mydata = pmap(.l= list(Ranges, Reference_Allele, Mutant, file_coding, file_noncoding),
                         .f = tabix_func)) %>% 
    unnest(cols = Mydata)
  
  # Make standard Cscape output format
  cscape_out <- cscape_annot %>% 
    separate(col = Ranges, into = c('Chr', 'Position', NA)) %>% 
    dplyr::rename(Reference = Reference_Allele) %>% 
    dplyr::select(-c(file_coding, file_noncoding)) %>% 
    relocate(Remark, .after = last_col()) %>% 
    mutate(across(.fns = parse_guess)) 
  
  return(cscape_out)
}





#' This function annotated a confidence level to the score
#' @param s the score
#' @param type coding or noncoding
#'
#' @return returns a confidence level or remark/error message 
#' @export
#' @examples
#' 
#' remark <- confidence(s, type)
#' 
confidence <- function(s, type){
  #   Returns the confidence level for (neutral, low-confidence or high-confidence)
  #     for the given p-value string using the given threshold.  If the string is invalid,
  #     returns 'no prediction'
  
  # Initialization: 
  #Confidence thresholds
  DEFAULT_THRESH <- 0.5
  CODING_THRESH <- 0.89
  NONCODING_THRESH <- 0.7
  if(type == 'Coding'){thresh <- CODING_THRESH
  }else{thresh <- NONCODING_THRESH}
  
  # Function: 
  x <- as.numeric(s)
  if (x >= thresh){
    return('High-confidence')
  }
  else if (x >= DEFAULT_THRESH){
    return('Low-confidence')
  }
  else {
    return('Neutral')
  }
  #Something on exception return 'No prediction'
}