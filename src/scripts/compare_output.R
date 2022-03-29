# Run on my local with sshfs to the server

library(tidyverse)

base_dir <- "~/Documents/sshfs/Analysis/Lori_Sussel/Nicole_Moss/eclip_normalization"

my_bed <- file.path(base_dir, "python_results",
                    "001_01.basedon_001_01.peaks.l2inputnormnew.bed")

their_bed <- file.path(base_dir, "data",
                       "test_compression_r2", 
                       "001_01.basedon_001_01.peaks.l2inputnormnew.bed")

my_bed_file <- read.table(my_bed)

their_bed_file <- read.table(their_bed)

all_colnames <- c("chromosome", "start", "end", "log_p", "logfc", "strand")

colnames(my_bed_file) <- all_colnames

colnames(their_bed_file) <- all_colnames

my_bed_file <- my_bed_file %>%
  dplyr::mutate(position = paste0(chromosome, "_", start, "_", end, "_", strand))

their_bed_file <- their_bed_file %>%
  dplyr::mutate(position = paste0(chromosome, "_", start, "_", end, "_", strand))

their_bed_file <- their_bed_file[match(my_bed_file$position,
                                               their_bed_file$position),]

all.equal(my_bed_file$position, their_bed_file$position)

all.equal(my_bed_file$log_p, their_bed_file$log_p)

all.equal(my_bed_file$logfc, their_bed_file$logfc)

my_bed_file$log_p_round <- round(my_bed_file$log_p, 3)
their_bed_file$log_p_round <- round(their_bed_file$log_p, 3)

all.equal(my_bed_file$log_p_round, their_bed_file$log_p_round)