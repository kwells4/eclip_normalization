# This was run on my local using sshfs to connect to the server

library(tidyverse)

base_dir <- "~/Documents/sshfs/Analysis/Lori_Sussel/Nicole_Moss/eclip_normalization/data/samtools_sort_r2"

file_one <- file.path(base_dir, "001_01.basedon_001_01.peaks.l2inputnormnew.bed")
file_two <- file.path(base_dir, "001_01.basedon_001_01.peaks.l2inputnormnew.bed.compressed.bed")

bed_one <- read.table(file_one)
bed_two <- read.table(file_two)

colnames <- c("chromosome", "start", "end", "p_val", "logfc", "strand")

colnames(bed_one) <- colnames
colnames(bed_two) <- colnames

bed_one <- bed_one %>%
  dplyr::mutate(position = paste0(chromosome, "_", start, "_", end, "_", strand))

bed_two <- bed_two %>%
  dplyr::mutate(position = paste0(chromosome, "_", start, "_", end, "_", strand))

# Try with all:
positions <- setdiff(bed_one$position, bed_two$position)

all_overlaps <- lapply(positions, function(x){
  all_info = strsplit(x, split = "_")
  start_pos = as.integer(all_info[[1]][2])
  end_pos = as.integer(all_info[[1]][3])
  chromosome_keep = all_info[[1]][1]
  
  # Ways it works:
  # start_pos <= start <= end_pos (the start is between the start and end)
  # start_pos <= end <= end_pos (the end is between the start and end)
  # start <= start_pos and end >= end_pos(the new read is longer than the original)
  # start >= start_pos and end <= end_pos (the new read is shorter than the original)
  
  test_bed_one <- bed_one %>%
    dplyr::filter(chromosome == chromosome_keep) %>%
    dplyr::filter(start >= start_pos & start <= start_pos)
  
  test_bed_two <- bed_one %>%
    dplyr::filter(chromosome == chromosome_keep) %>%
    dplyr::filter(end >= start_pos & end <= start_pos)
  
  test_bed_three <- bed_one %>%
    dplyr::filter(chromosome == chromosome_keep) %>%
    dplyr::filter(start <= start_pos & end >= end_pos)
  
  test_bed_four <- bed_one %>%
    dplyr::filter(chromosome == chromosome_keep) %>%
    dplyr::filter(start >= start_pos & end <= end_pos)
  
  
  overlaps <- do.call(rbind, list(test_bed_one, test_bed_two, test_bed_three, test_bed_four)) %>%
    distinct()
  
  return(overlaps)
  
})

all_overlaps <- do.call(rbind, all_overlaps)

write.table(all_overlaps, file.path(base_dir, "all_overlaps.bed"))
