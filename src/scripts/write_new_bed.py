import os
from collections import defaultdict

output_directory = "/beevol/home/wellskri/Analysis/Lori_Sussel/Nicole_Moss/eclip_normalization/data/test_compression_r2"
input_directory = "/beevol/home/wellskri/Analysis/Lori_Sussel/Nicole_Moss/eclip_normalization/data/samtools_sort_r2"

overlapping_bed = os.path.join(input_directory, "all_overlaps.bed")
input_bed = os.path.join(input_directory, "Rbfox2_CLIP_S24_r2.peaks.bed")
output_bed = os.path.join(output_directory, "Rbfox2_CLIP_S24_r2.peaks.bed")

bed_dict = defaultdict(int)
with open(overlapping_bed, "r") as in_bed:
	next(in_bed)
	for line in in_bed:
		line = line.strip().split("\t")
		position = line[1] + "_" + line[2] + "_" + line[3]
		bed_dict[position] += 1

with open(input_bed, "r") as new_in_bed, open(output_bed, "w") as out_bed:
	for line in new_in_bed:
		test_line = line.strip().split("\t")
		position = test_line[0] + "_" + test_line[1] + "_" + test_line[2]
		if bed_dict[position] >= 1:
			out_bed.write(line)