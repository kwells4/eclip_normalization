library(tidyverse)

base_dir <- "~/Documents/sshfs/Analysis/Lori_Sussel/Nicole_Moss/eclip_normalization"

my_bed <- file.path(base_dir, "python_results_all",
                    "001_01.basedon_001_01.peaks.l2inputnormnew.bed")

their_bed <- file.path(base_dir, "data",
                       "samtools_sort_r2", 
                       "001_01.basedon_001_01.peaks.l2inputnormnew.bed")

my_bed_file <- read.table(my_bed)

their_bed_file <- read.table(their_bed)

all_colnames <- c("chromosome", "start", "end", "log_p", "logfc", "strand")

my_colnames <- c(all_colnames, "clip_reads", "input_reads", "total_clip_reads",
                 "total_input_reads")

colnames(my_bed_file) <- my_colnames

colnames(their_bed_file) <- all_colnames

my_bed_file <- my_bed_file %>%
  dplyr::mutate(position = paste0(chromosome, "_", start, "_", end, "_", strand))

their_bed_file <- their_bed_file %>%
  dplyr::mutate(position = paste0(chromosome, "_", start, "_", end, "_", strand))

my_bed_file <- my_bed_file %>%
  dplyr::filter(log_p != "Inf")

their_bed_file <- their_bed_file[match(my_bed_file$position,
                                               their_bed_file$position),]

all.equal(my_bed_file$position, their_bed_file$position)

all.equal(my_bed_file$log_p, their_bed_file$log_p)

all.equal(my_bed_file$logfc, their_bed_file$logfc)

my_bed_file$log_p_round <- round(my_bed_file$log_p, 2)
their_bed_file$log_p_round <- round(their_bed_file$log_p, 2)

all.equal(my_bed_file$log_p_round, their_bed_file$log_p_round)

my_bed_file$their_log_p <- their_bed_file$log_p_round

my_bed_file$different <- my_bed_file$their_log_p != my_bed_file$log_p_round

unmatched <- my_bed_file %>% dplyr::filter(different == TRUE & their_log_p != 15.66 &
                              their_log_p != 400)


write.table(unmatched, file.path(base_dir, "data/unmatched.bed"), sep = "\t",
            quote = FALSE)