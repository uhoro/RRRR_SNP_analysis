
library(tidyverse)
library(dplyr)
library(stringr)
library(optparse)



get_corrected_gap_gene_IPS <- function(flexa_IPS_shaped_gap_gene, gene_id, gap_info_gap_gene_only) {
  
  if (nrow(flexa_IPS_shaped_gap_gene) != 0) {
    REF_gap_str <- gap_info_gap_gene_only$REF[gap_info_gap_gene_only$gene_id == gene_id]
    gap_blocks <- str_split(REF_gap_str, ",")[[1]]
    gap_starts <- c()
    gap_stops <- c()
    for (gap_i in 1:length(gap_blocks)) {
      gap_starts <- c(gap_starts, round(as.numeric(str_split(gap_blocks[gap_i], "-")[[1]][1])/3,0)) # divide by three to convert bp to codon/AA
      gap_stops<- c(gap_stops, round(as.numeric(str_split(gap_blocks[gap_i], "-")[[1]][2])/3,0))
    }
    # print(paste("gap_starts"))
    # print(gap_starts)
    # print(paste("gap_stops"))
    # print(gap_stops)
    for (ips_i in 1:nrow(flexa_IPS_shaped_gap_gene)) {
      start_location_before_ips_i <- as.numeric(flexa_IPS_shaped_gap_gene$start_location_before_correction[ips_i])
      start_location_after_ips_i <- start_location_before_ips_i
      stop_location_before_ips_i <- as.numeric(flexa_IPS_shaped_gap_gene$stop_location_before_correction[ips_i])
      stop_location_after_ips_i <- stop_location_before_ips_i
      # print(start_location_before_ips_i)
      # print(stop_location_before_ips_i)
      for (gap_i in 1:length(gap_blocks)) {
        if (stop_location_before_ips_i <= gap_starts[gap_i]) {
          break
        } else {
          if (start_location_before_ips_i <= gap_starts[gap_i]) {
            stop_location_after_ips_i <- stop_location_after_ips_i + (gap_stops[gap_i] - gap_starts[gap_i] + 1)
            
          } else {
            start_location_after_ips_i <- start_location_after_ips_i + (gap_stops[gap_i] - gap_starts[gap_i] + 1)
            stop_location_after_ips_i <- stop_location_after_ips_i + (gap_stops[gap_i] - gap_starts[gap_i] + 1)
          }
          
        }
      }
      flexa_IPS_shaped_gap_gene$start_location[ips_i] <- start_location_after_ips_i
      flexa_IPS_shaped_gap_gene$stop_location[ips_i] <- stop_location_after_ips_i
    }
  }
  
  return(flexa_IPS_shaped_gap_gene)
}


# Define options
option_list <- list(
  make_option(c("-u", "--corrected_ipr_fname"), type = "character", 
              default = "flexa_IPS_corrected.csv",
              help = "corrected ipr output file name. default: flexa_IPS_corrected.csv", metavar = "MESSAGE"),
  make_option(c("-d", "--out_dir"), type = "character",
              help = "output directory", metavar = "MESSAGE"),
  make_option(c("-i", "--original_ipr_path"), type = "character",
              help = "Path to the original ipr path", metavar = "MESSAGE"),
  make_option(c("-r", "--REF_gap_data"), type = "character", default = NA,
              help = "Path to the reference gap information csv file", metavar = "MESSAGE")
)

# Create the parser object
opt_parser <- OptionParser(option_list = option_list)
# Parse the arguments
opt <- parse_args(opt_parser)

# Access the arguments
cat("Corrected ipr file name:", opt$corrected_ipr_fname, "\n")
cat("out_dir: ", opt$out_dir, "\n")
cat("original interpro path: ", opt$original_ipr_path, "\n")
cat("Reference gap data: ", opt$REF_gap_data, "\n")


# Update accordingly the PATH to interpro scan output in csv format
flexa_IPS <- read.csv(opt$original_ipr_path)
flexa_IPS %>%
  separate(protein_id, c("geneid", "T1"), "-", remove = FALSE) -> flexa_IPS_shaped

# Update accordingly the PATH to gap information in ref sequence
gap_info_ref <- read.csv(opt$REF_gap_data)

# gap info in base pair whereas start_location and stop_location are in AA

flexa_IPS_shaped$start_location_before_correction <- flexa_IPS_shaped$start_location
flexa_IPS_shaped$stop_location_before_correction <- flexa_IPS_shaped$stop_location

gap_info_gap_gene_only <- gap_info_ref %>%
  filter(REF != "") 

flexa_IPS_shaped_corrected <- flexa_IPS_shaped %>%
  filter(!geneid %in% gap_info_gap_gene_only$gene_id)

n_processed_gene <- 0
print(paste("Total number of genes to process: ", length(gap_info_gap_gene_only$gene_id)))
for (gene_id in gap_info_gap_gene_only$gene_id) {
  n_processed_gene<- n_processed_gene +1
  if (n_processed_gene %% 100 == 0) {
    print(paste("processed ", n_processed_gene, " genes"))
  }
  flexa_IPS_shaped_gap_gene <- flexa_IPS_shaped %>%
    filter(geneid == gene_id)
  flexa_IPS_shaped_gap_gene_corrected <- get_corrected_gap_gene_IPS(flexa_IPS_shaped_gap_gene,
                                                                    gene_id,
                                                                    gap_info_gap_gene_only)
  flexa_IPS_shaped_corrected <- rbind(flexa_IPS_shaped_corrected, flexa_IPS_shaped_gap_gene_corrected)
}  

out_full_path <- paste(opt$out_dir, opt$corrected_ipr_fname, sep = "/")
write.csv(flexa_IPS_shaped_corrected,out_full_path, row.names = FALSE) 


