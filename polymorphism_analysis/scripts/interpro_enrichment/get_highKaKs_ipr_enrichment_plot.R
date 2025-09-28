library(optparse)

get_overlap2 <- function(Feature_set, set_KaKs_pos) {
  print(paste("total number of features: ", length(Feature_set$start_location)))
  res <- rep(FALSE, length(Feature_set$start_location))
  for (i in 1:length(Feature_set$start_location)) {
    pos_info_vec <- set_KaKs_pos$pos_info2[set_KaKs_pos$geneid == Feature_set$geneid[i]]
    if (i %% 10000 == 0) {
      print(paste("processing ", i, " th feature"))
    }
    if (length(pos_info_vec) == 0) {
      next
    }
    #print(pos_info_vec)
    feature_start <- as.numeric(Feature_set$start_location[i])
    feature_end <- as.numeric(Feature_set$stop_location[i])
    for (j in 1:(length(pos_info_vec)/2)) {
      start_pos_AA <- as.numeric(substring(pos_info_vec[2*j - 1], 2))/3 # convert to AA count
      end_pos_AA <- as.numeric(substring(pos_info_vec[2*j], 2))/3 # convert to AA count
      # print(paste("feature_end: ", feature_end))
      # print(paste("start_pos_AA: ", start_pos_AA))
      # print(paste("feature_start: ", feature_start))
      # print(paste("end_pos_AA: ", end_pos_AA))
      if (!(feature_end < start_pos_AA | end_pos_AA < feature_start)) {
        res[i] <- TRUE
        break
      }
      
    }
  }
  return(res)
}

count_no_overlap2 <- function(Feature_set, set_KaKs_pos) {
  print(paste("total number of KaKspos: ", length(set_KaKs_pos$geneid)/2))
  res <- rep(FALSE, length(set_KaKs_pos$geneid)/2)
  for (i in 1:(length(set_KaKs_pos$geneid)/2)) {
    start_pos_AA <- as.numeric(substring(set_KaKs_pos$pos_info2[2*i - 1], 2))/3 # convert to AA count
    end_pos_AA <- as.numeric(substring(set_KaKs_pos$pos_info2[2*i], 2))/3 # convert to AA count
    
    feature_start_vec <- as.numeric(Feature_set$start_location[Feature_set$geneid == set_KaKs_pos$geneid[2*i]])
    feature_stop_vec <- as.numeric(Feature_set$stop_location[Feature_set$geneid == set_KaKs_pos$geneid[2*i]])
    if (i %% 100 == 0) {
      print(paste("processing ", i, " th KaKs"))
    }
    if (length(feature_start_vec) == 0) {
      next
    }
    
    for (j in 1:(length(feature_start_vec))) {
      if (!(feature_stop_vec[j] < start_pos_AA | end_pos_AA < feature_start_vec[j])) {
        res[i] <- TRUE
        break
      }
      
    }
  }
  return(res)
}


get_highKaKs_ipr_enrichment <- function(sample1, 
                                        sample2, 
                                        KaKs_high_path,
                                        out_dir,
                                        genes_with_errors_path,
                                        interpro_path, suffix = "") {

  

  print(paste("KaKs high path: ", KaKs_high_path))
  
  library(tidyverse)
  library(dplyr)
  if (is.na(genes_with_errors_path)){
    print("Use all genes")
    genes_with_errors <- vector()
  } else {
    print("Use genes in-frame")
    genes_with_errors_raw <- readLines(genes_with_errors_path)
    genes_with_errors <- vector()
    for (i in 1:length(genes_with_errors_raw)) {
      vec <- str_split(genes_with_errors_raw[i], "_")[[1]]
      geneid <- paste("FUN", tail(vec, n=1), sep = "_")
      genes_with_errors <- c(genes_with_errors, geneid)
      
    }
    print(paste("genes_with_errors example: ", genes_with_errors[1]))
  }
  
  KaKs_high_df <- read.csv(KaKs_high_path)
  
  KaKs_high_df %>%
    dplyr::select(geneid, pos_info2) %>%
    separate_longer_delim(pos_info2, delim = "_") %>%
    filter(pos_info2 != "pos") -> Set_KaKs2_pos_original
  
  KaKs_high_df %>%
    dplyr::select(geneid, pos_info2) %>%
    separate_longer_delim(pos_info2, delim = "_") %>%
    filter(pos_info2 != "pos") %>%
    filter(! geneid %in% genes_with_errors) -> Set_KaKs2_pos
  
  print(paste("Nunber of KaKs2 positions removed: ",
              length(Set_KaKs2_pos_original$geneid) - length(Set_KaKs2_pos$geneid)))
  
  flexa_IPS <- read.csv(interpro_path)
  flexa_IPS %>%
    separate(protein_id, c("geneid", "T1"), "-", remove = FALSE) -> flexa_IPS_shaped
  
  flexa_IPS_shaped$In_KaKs2_region <- get_overlap2(flexa_IPS_shaped, Set_KaKs2_pos)
  res <- count_no_overlap2(flexa_IPS_shaped, Set_KaKs2_pos)
  Set_KaKs2_pos$with_feature <- rep(res, each = 2)
  ratio_withFeature <- sum(Set_KaKs2_pos$with_feature)/length(Set_KaKs2_pos$with_feature)
  print(paste("Ratio of high KaKs region with ipr feature: ", ratio_withFeature))
  # Filter only interpro annotation 
  flexa_IPS_shaped %>%
    filter(interpro_accession != "-") -> flexa_IPS_shaped_filt1
  
  # Filter genes with errors
  flexa_IPS_shaped_filt1 %>%
    filter(!geneid %in% genes_with_errors) ->flexa_IPS_shaped_filt
  
  # Filter redundant interpro annotation
  flexa_IPS_shaped_filt$rowid <- 1:length(flexa_IPS_shaped_filt$protein_accession)
  geneid_vec <- unique(flexa_IPS_shaped_filt$geneid)
  rowid_toremove <- c()
  i <- 0
  print(paste("Length of geneid_vec: ", length(geneid_vec)))
  for (gi in 1:length(geneid_vec)) {
    if (gi %% 1000 == 0) {
      print(paste("processing ", gi, " th gene"))
    }
    #gi <- 1
    gene <- geneid_vec[gi]
    flexa_IPS_shaped_filt %>%
      dplyr::filter(geneid == gene) -> flexa_IPS_shaped_filt_gene
    IPR_vec <- unique(flexa_IPS_shaped_filt_gene$interpro_accession)
    for (IPRi in 1:length(IPR_vec)) {
      #IPRi <- 1
      IPRid <- IPR_vec[IPRi]
      flexa_IPS_shaped_filt_gene %>%
        dplyr::filter(interpro_accession == IPRid) -> flexa_IPS_shaped_filt_gene_IPR
      if (length(flexa_IPS_shaped_filt_gene_IPR$start_location) == 1) {
        next
      } else {
        for (i in 1:(length(flexa_IPS_shaped_filt_gene_IPR$start_location)-1)) {
          #print(paste("i ", i))
          start <- flexa_IPS_shaped_filt_gene_IPR$start_location[i]
          stop <- flexa_IPS_shaped_filt_gene_IPR$stop_location[i]
          for (j in (i+1):length(flexa_IPS_shaped_filt_gene_IPR$start_location)) {
            #print(paste("j ", j))
            start2 <- flexa_IPS_shaped_filt_gene_IPR$start_location[j]
            stop2 <- flexa_IPS_shaped_filt_gene_IPR$stop_location[j]
            if (!(as.numeric(stop2) < as.numeric(start) | as.numeric(stop) < as.numeric(start2))) {
              rowid_toremove <- c(rowid_toremove, flexa_IPS_shaped_filt_gene_IPR$rowid[j])
            }
          }
        }
      }
    }
    
  }
  
  flexa_IPS_shaped_filt %>%
    filter(!(rowid %in% rowid_toremove)) -> flexa_IPS_shaped_filt_nr
  feature_n <- length(flexa_IPS_shaped_filt_nr$start_location)
  print(paste("number of non redundant features: ", feature_n))
  inKaKs2_n <- sum(flexa_IPS_shaped_filt_nr$In_KaKs2_region)
  print(paste("number of non redundant features in KaKs2: ", inKaKs2_n))
  
  flexa_IPS_shaped_filt_nr %>%
    mutate(
      geneid2 = geneid
    ) %>%
    group_by(interpro_accession) %>%
    reframe(
      total_freq = n(),
      inKaKs2 = sum(In_KaKs2_region),
      geneid_vec = paste0(geneid, collapse = ","),
      interpro_description = interpro_description[1]
    ) -> flexa_IPS_shaped_filt_nr_reframed
  
  flexa_IPS_shaped_filt_nr %>% 
    filter(In_KaKs2_region) %>%
    mutate(geneid2 = geneid) %>%
    group_by(interpro_accession) %>%
    reframe(
      geneid_vec_wKaKs2 = paste0(geneid, collapse = ",")
    ) -> flexa_IPS_shaped_filt_nr_overlap_reframed
  
  flexa_IPS_shaped_filt_nr_reframed %>%
    left_join(flexa_IPS_shaped_filt_nr_overlap_reframed,
              by = join_by(interpro_accession)) -> flexa_IPS_shaped_filt_nr_reframed
  
  flexa_IPS_shaped_filt_nr_reframed %>%
    mutate(
      expected_freq = total_freq * sum(flexa_IPS_shaped_filt_nr$In_KaKs2_region) / length(flexa_IPS_shaped_filt_nr$start_location),
      IDlong = paste(interpro_description, interpro_accession)
    ) -> flexa_IPS_KaKs2_freq
  
  fisher_res <- rep(1, length(flexa_IPS_KaKs2_freq$interpro_accession))
  outKaKs2_n <- feature_n - inKaKs2_n
  for (i in 1:length(flexa_IPS_KaKs2_freq$interpro_accession)) {
    dat <- data.frame(
      Background = c(flexa_IPS_KaKs2_freq$total_freq[i] - flexa_IPS_KaKs2_freq$inKaKs2[i], 
                     outKaKs2_n - (flexa_IPS_KaKs2_freq$total_freq[i] - flexa_IPS_KaKs2_freq$inKaKs2[i])),
      Set = c(flexa_IPS_KaKs2_freq$inKaKs2[i], inKaKs2_n - flexa_IPS_KaKs2_freq$inKaKs2[i]))
    fisher_res[i] <- fisher.test(dat, alternative = "less")$p.value
  }
  flexa_IPS_KaKs2_freq$fisher_res <- fisher_res
  
  flexa_IPS_KaKs2_freq %>%
    mutate(
      IDlong = fct_reorder(IDlong, -log10(fisher_res))
    ) -> flexa_IPS_KaKs2_freq
  
  flexa_IPS_KaKs2_freq <- flexa_IPS_KaKs2_freq[order(flexa_IPS_KaKs2_freq$fisher_res),] 
  flexa_IPS_KaKs2_freq %>%
    mutate(
      fisher_res_BF = p.adjust(fisher_res, method="bonferroni"),
      fisher_res_Holm = p.adjust(fisher_res, method="holm"),
      fisher_res_BH = p.adjust(fisher_res, method="BH"),
      fisher_res_BY = p.adjust(fisher_res, method="BY"),
    ) -> flexa_IPS_KaKs2_freq_stat
  
  KaKs2_freq_stat_fname <- paste(out_dir, "/flexa_IPS_KaKs2_freq_stat_", sample1, 
                                 "_", sample2, suffix, ".csv",
                                 sep = "")
  KaKs2_freq_stat_fname
  write.csv(flexa_IPS_KaKs2_freq_stat, KaKs2_freq_stat_fname, row.names = FALSE)
  
  # 
  flexa_IPS_KaKs2_freq_sample_pair <- read.csv(KaKs2_freq_stat_fname)
  
  high_KaKs_ipr_enrich_df <- flexa_IPS_KaKs2_freq_sample_pair %>%
    arrange(fisher_res)
  
  if (!file.exists(paste(out_dir, "/enrichment_plots", sep = ""))) {
    dir.create(paste(out_dir, "/enrichment_plots", sep = ""))
  }
  
  title_full <- paste("Enriched non redundant interpro feature in KaKs > 2 regions", 
                      sample1, sample2)
  p <- high_KaKs_ipr_enrich_df[1:20,] %>%
    # filter(IPR_ID %in% c("IPR000242",  
    #                      "IPR013783", "IPR000742", 
    #                      "IPR003410", "IPR015919", 
    #                      "IPR003961", "IPR000980")) %>%
    ggplot(aes(x = -log10(fisher_res), y = reorder(IDlong, -log10(fisher_res)), fill = inKaKs2)) +
    geom_bar(stat = "identity", aes(fill=inKaKs2), width = 0.5) + 
    #scale_fill_gradient(limits = c(0,14), low='#FFFFFF', high='#FF0000') + 
    scale_fill_gradient(limits = c(0,10), low='#FFFFFF', high='#40e0d0') + 
    #geom_point(shape = 1,size = 12,colour = "black")+ 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    ylab("InterProScan feature") +
    #xlim(0,7) +
    xlab("- log10 p value") +
    ggtitle(title_full) +
    theme(plot.title = element_text(hjust = 1))
  
  top20_plot_name <- paste(out_dir, "/enrichment_plots/KaKs_IPS_overlap_top20_", 
                           sample1, "_", sample2, "_",  suffix, ".svg",
                           sep = "")
  ggsave(top20_plot_name, width = 8, height = 8)
  
  high_KaKs_ipr_enrich_df_sort <- high_KaKs_ipr_enrich_df %>%
    arrange(desc(inKaKs2), fisher_res) %>%
    filter(fisher_res < 0.05)
  high_KaKs_ipr_enrich_df_sort$order <- 1:nrow(high_KaKs_ipr_enrich_df_sort)
  
  library(forcats)
  p <- high_KaKs_ipr_enrich_df_sort[1:20,] %>%
    # filter(IPR_ID %in% c("IPR000242",  
    #                      "IPR013783", "IPR000742", 
    #                      "IPR003410", "IPR015919", 
    #                      "IPR003961", "IPR000980")) %>%
    ggplot(aes(x = inKaKs2, y = reorder(IDlong, desc(order)), fill = -log10(fisher_res))) +
    geom_bar(stat = "identity", aes(fill=-log10(fisher_res)), width = 0.5) + 
    #scale_fill_gradient(limits = c(0,14), low='#FFFFFF', high='#FF0000') + 
    scale_fill_gradient(limits = c(0,6), low='#FFFFFF', high='#40e0d0') + 
    #geom_point(shape = 1,size = 12,colour = "black")+ 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    ylab("InterProScan feature") +
    #xlim(0,7) +
    xlab("Number of features with high KaKs value") +
    ggtitle("Enriched interpro feature in KaKs > 2 regions") +
    theme(plot.title = element_text(hjust = 1))
  

  top20_plot_name <- paste(out_dir, 
                           "/enrichment_plots/KaKs_IPS_overlap_top20_by_DomainCount_", 
                           sample1, "_", sample2, "_", suffix, ".svg",
                           sep = "")

  ggsave(top20_plot_name, width = 8, height = 8)
  
}

# Define options
option_list <- list(
  make_option(c("-x", "--sample1"), type = "character", default = "C1G11",
              help = "sample 1 id", metavar = "NUMBER"),
  make_option(c("-y", "--sample2"), type = "character", default = "C1C5v2",
              help = "sample 2 id", metavar = "MESSAGE"),
  make_option(c("-d", "--out_dir"), type = "character",
              help = "output directory", metavar = "MESSAGE"),
  make_option(c("-k", "--KaKs_high_path"), type = "character",
              help = "Path to high KaKs csv data", metavar = "MESSAGE"),
  make_option(c("-e", "--error_txt"), type = "character", default = NA,
              help = "Gene list with errors", metavar = "MESSAGE"),
  make_option(c("-p", "--ipr_path"), type = "character", 
              help = "Path to the updated interpro data", metavar = "MESSAGE"),
  make_option(c("-f", "--suffix"), type = "character",
              default = "",
              help = "additional suffix for output figure (optional)", metavar = "MESSAGE")
)

# Create the parser object
opt_parser <- OptionParser(option_list = option_list)
# Parse the arguments
opt <- parse_args(opt_parser)

# Access the arguments
cat("Sample 1:", opt$sample1, "\n")
cat("Sample 2:", opt$sample2, "\n")
cat("KaKs_high_path:", opt$KaKs_high_path, "\n")
cat("out_dir: ", opt$out_dir, "\n")
cat("gene error list: ", opt$error_txt, "\n")
cat("interpro path: ", opt$ipr_path, "\n")
cat("figure suffix: ", opt$suffix, "\n")

get_highKaKs_ipr_enrichment(opt$sample1, opt$sample2, 
                            opt$KaKs_high_path, opt$out_dir, 
                            opt$error_txt, opt$ipr_path, opt$suffix)



