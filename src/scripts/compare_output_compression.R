# This was run on my local with an sshfs connection

library(tidyverse)

base_dir <- "~/Documents/sshfs/Analysis/Lori_Sussel/Nicole_Moss/eclip_normalization"

my_bed <- file.path(base_dir, "python_results",
                    "001_01.basedon_001_01.peaks.l2inputnormnew.bed.compressed.bed")

their_bed <- file.path(base_dir, "data",
                       "test_compression_r2", 
                       "001_01.basedon_001_01.peaks.l2inputnormnew.bed.compressed.bed")

my_bed_file <- read.table(my_bed)

their_bed_file <- read.table(their_bed)

all_colnames <- c("chromosome", "start", "end", "log_p", "logfc", "strand")

colnames(my_bed_file) <- all_colnames

colnames(their_bed_file) <- all_colnames

my_bed_file <- my_bed_file %>%
  dplyr::mutate(position = paste0(chromosome, "_", start, "_", end, "_", strand)) %>%
  dplyr::mutate(log_p_round = round(log_p, 3))


their_bed_file <- their_bed_file %>%
  dplyr::mutate(position = paste0(chromosome, "_", start, "_", end, "_", strand))%>%
  dplyr::mutate(log_p_round = round(log_p, 3))

lapply(unique(my_bed_file$chromosome), function(x){
  print(x)
  new_my_bed <- my_bed_file %>% dplyr::filter(chromosome == all_of(x))
  
  new_their_bed <- their_bed_file %>% dplyr::filter(chromosome == all_of(x))
  
  if (!setequal(new_my_bed$log_p_round, new_their_bed$log_p_round)){
    print(new_my_bed)
    print(new_their_bed)
  }


})
