#' This function looks for replicates in the omics data frame and replaces replicates with sample mean 
#'
#' @param omics_df is a data frame containing only samples as rows
#'
#' @return
#' @export
#'
#' @examples
#' dataFilt_hugo <- find_replicates(dataFilt_hugo)
find_replicates <- function(omics_df){
  # Find total sample codes
  total_samples <- str_sub(colnames(omics_df), 1, 15)
  
  # Find replicated sample codes
  rep_samples <- unique(total_samples[duplicated(total_samples)])
  
  for (sample in rep_samples){
    temp <- omics_df %>% 
      select(starts_with(sample))
    
    # Calculate mean of the replicated samples
    omics_df <- cbind(omics_df, apply(temp, 1, mean, na.rm=TRUE))
    
    # Rename last column
    colnames(omics_df)[length(colnames(omics_df))] <- sample
    # Remain integers
    omics_df[,length(colnames(omics_df))] <- round(omics_df[,length(colnames(omics_df))])
    # Delete replicates
    omics_df <- omics_df[, !colnames(omics_df) %in% colnames(temp)]
  }
  return(omics_df)
}








#' This function wrangle the MAF file to the total number of mutations for each gene-sample pair
#'
#' @param mu_df should be a data frame of the mutation data
#' @param samples is a vector of samples for which you want to build the mutation table
#'
#' @return
#' @export
#'
#' @examples 
#' Mu_wrangle_number(mu_df = BRCA_mu,
#' samples = common_samples)
Mu_wrangle_number <- function(mu_df, samples){
  # Preprocessing
  mu_df <- mu_df %>% 
    # Add sample code information
    mutate(Sample = str_sub(Tumor_Sample_Barcode, 1, 15)) %>% 
    # Keep input samples
    filter(Sample %in% samples)
  
  BRCA_mu_number <- mu_df %>%
    # Find total mutation number for each gene-patient pair
    group_by(Hugo_Symbol, Sample) %>% 
    summarise(Total_mu = n()) %>% 
    pivot_wider(id_cols = Hugo_Symbol,
                names_from = Sample,
                values_from = Total_mu)
  
  # Replace NA with 0
  BRCA_mu_number[is.na(BRCA_mu_number)] <- 0
  
  return(BRCA_mu_number)
}






#' This function wrangles the MAF file to indicate the mutation type information of the genes for each sample
#'
#' @param mu_df should be a data frame of the mutation data
#' @param samples is a vector of samples for which you want to build the mutation table
#' @param size is a character with two options "large" or "small" that indicates the size of the mutations types to include
#'
#' @return a large list of data frames
#' @export
#'
#' @examples
#' mutation_list <- Mu_wrangle_type(mu_df = BRCA_mu,
#' samples = common_samples, cscape = FALSE)
Mu_wrangle_type <- function(mu_df, samples, cscape){
  
  # Preprocessing
  mu_df <- mu_df %>% 
    # Add sample code information
    mutate(Sample = str_sub(Tumor_Sample_Barcode, 1, 15)) %>% 
    # Keep input samples
    filter(Sample %in% samples)
  
  # Make empty list
  sample_mutation_list <- list()
  
  for (sample in samples){
    
    # Find unique mutated genes in sample
    sample_df <- mu_df %>% 
      filter(Sample == sample)
    
    # Extract unique features
    genes <- unique(sample_df$Hugo_Symbol)
    variant_class <- unique(sample_df$Variant_Classification)
    
    # Generate empty mutation matrix with only variant classification
    mutation <- matrix(data = 0,
                       nrow = length(genes), 
                       ncol = 9,
                       dimnames = list(genes,
                                       c("5'Flank", "5'UTR", "Nonsense_Mutation",
                                         "Splice_Site", "Splice_Region", "Intron",
                                         "3'UTR", "3'Flank", "RNA"))) # according to consequence table
    # Add matrix values
    for (gene in genes){
      
      # Count driving mutations based on CScape score for each gene
      if (cscape == TRUE){
        driving_mu <- sample_df %>% 
          filter(Hugo_Symbol == gene) %>% 
          filter(CScape_Mut_Class == "Driver")
        
        if (nrow(driving_mu) != 0){
          # Assign 1 to certain mutation type with driving mutations
          mutation_name <- driving_mu$Variant_Classification
          for (name in mutation_name){
            if (name %in% colnames(mutation)){
              mutation[gene, name] <- 1
            }
          }
        }
      }
      # Count all mutation numbers
      else if (cscape == FALSE){
        driving_mu <- sample_df %>% 
          filter(Hugo_Symbol == gene)
        
        mutation_name <- driving_mu$Variant_Classification
        for (name in mutation_name){
          if (name %in% colnames(mutation)){
            mutation[gene, name] <- mutation[gene, name] + 1
          }
        }
      }
      else {
        print("Please assign a logical value to cscape variable.")
      }
    }
    # Convert matrix and add to list
    sample_mutation_list[[sample]] <- tibble::rownames_to_column(as.data.frame(mutation), "Hugo_Symbol")
  }
  return(sample_mutation_list)
}




#' The function wrangles the CNV file from GISTIC 2.0 with certain samples
#'
#' @param CNV_gistic_data is the resulting data (all_data_by_gene.txt) from GISTIC 2.0
#' @param samples is a vector of samples for which you want to build the CNV table
#'
#' @return
#' @export
#'
#' @examples
#' BRCA_CNV_gistic <- CNV_wrangle(CNV_gistic_data = BRCA_CNV_gistic,
#' samples = common_samples, size = "large")
CNV_wrangle <- function(CNV_gistic_data, samples){
  # Rename sample columns
  colnames(CNV_gistic_data)[-c(1:3)] <- gsub("\\.", "-", colnames(CNV_gistic_data)[-c(1:3)])
  
  # Find replicates
  non_replicate_df <- find_replicates(CNV_gistic_data[, -c(1:3)])
  
  # Add other columns
  CNV_gistic_data <- cbind(CNV_gistic_data[, c(1:3)], non_replicate_df)
  
  # Extract only sample code
  colnames(CNV_gistic_data)[-c(1:3)] <- str_sub(colnames(CNV_gistic_data)[-c(1:3)], 1, 15)
  
  # Keep overlapping samples
  CNV_gistic_data <- CNV_gistic_data %>% 
    select(Gene.Symbol, colnames(.)[which(colnames(.) %in% samples)])
  
  return(CNV_gistic_data)
}





#' The function produces the gene-level methylation data
#'
#' @param methy_data is Methylation data with beta values for samples 
#' @param probe_data is the corresponding probe file for methylation data
#' @param samples is a vector of samples for which you want to build the methylation table
#'
#' @return
#' @export
#'
#' @examples
#' BRCA_met_wrangled <- Methylation_wrangle(methy_data = BRCA_met,
#' probe_data = probes,
#' samples = common_samples)
Methylation_wrangle <- function(methy_data, probe_data, samples){
  # Rename sample columns
  colnames(methy_data) <- gsub("\\.", "-", colnames(methy_data))
  
  # find replicates
  methy_data <- find_replicates(methy_data)
  
  # Extract only sample code
  colnames(methy_data) <- str_sub(colnames(methy_data), 1, 15)
  
  BRCA_met <- methy_data %>%
    # Keep only expected samples
    select(colnames(.)[which(colnames(.) %in% samples)]) %>% 
    # Convert row names to column
    rownames_to_column(var = "probeID")
  
  BRCA_met_wrangled <- probe_data %>% 
    select(probeID, genesUniq, transcriptTypes, CGIposition) %>% 
    # Merge probe data and met data
    inner_join(BRCA_met, by = "probeID") %>% 
    # Delete NA
    drop_na() %>% 
    # Keep only protein coding genes(why??)
    # filter(grepl('protein_coding', transcriptTypes)) %>%
    # Keep only genesUniq with one gene/transcript
    filter(!grepl(';', genesUniq)) %>%
    # Calculate average beta value for each CGI position on the gene
    group_by(genesUniq, CGIposition) %>% 
    summarise_at(vars(-c(transcriptTypes, probeID)), mean) %>% 
    pivot_wider(names_from = CGIposition,
                values_from = -c(genesUniq, CGIposition),
                values_fill = 0)
  
  return(BRCA_met_wrangled)
  
}