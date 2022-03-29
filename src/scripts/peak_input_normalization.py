import numpy as np
from scipy.stats import chi2_contingency
import math
import pysam
import argparse
from collections import defaultdict
import re
import os
import subprocess
import sys


def main():
	# Initialize dictionaries and lists
    options = setup()

    make_output_dir(options.output_dir)
    sample_dict = {} # Outer level is sample id (save name), inner levels are input_bam, bed, clip_bam
    manifest_file = options.manifest_file
    read_count_file = options.manifest_file + ".mapped_read_num"

    make_input_dict(manifest_file, sample_dict)

    total_reads_dict = make_reads_dict(sample_dict, read_count_file)


    # For each sample in the input dict
    for sample in sample_dict:
        save_file = os.path.join(options.output_dir, sample +
            "_01.basedon_001_01.peaks.l2inputnormnew.bed")
        bed_file = sample_dict[sample]["bedfile"]
        input_bam = sample_dict[sample]["input_bam"]
        clip_bam = sample_dict[sample]["clip_bam"]

        with open(save_file, "w") as write_file, open(bed_file, "r") as in_bed:
            for line in in_bed:
                bed_dict = make_bed_dict(line)
                read_count_clip = count_bam_reads(bed_dict, clip_bam, options)
                read_count_input = count_bam_reads(bed_dict, input_bam, options)
                total_clip_reads = total_reads_dict[clip_bam]
                total_input_reads = total_reads_dict[input_bam]

                print(bed_dict)
                p_val, logfc = chi_square_or_fisher(read_count_clip, read_count_input, total_clip_reads, total_input_reads)
                write_file.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(bed_dict["chromosome"], bed_dict["start"], bed_dict["end"], p_val, logfc, bed_dict["strand"]))
                # save_clip = os.path.join("test_bams", clip_bam.split("/")[-1])
                # save_input  = os.path.join("test_bams", input_bam.split("/")[-1])
                # new_bam_file(bed_dict, clip_bam, options, save_clip)
                # new_bam_file(bed_dict, input_bam, options, save_input)

    # Open bed file
    # For each line in bed file
    # Pass to function to get reads from bam for just that region
    # Decide if we need chisquare or fisher exact
    # Run test
    # 2 options. 1 - write to output file, 2 - keep a dictionary. Not sure what is best
    # Make compressed output where overlapping peaks are discarded (just keep the top hit)
    # This maybe is better in R?





#########
# Setup #
#########

def setup():
    """
    Gets command line arguments and returns a Namespace object
    """

    # House keeping to read in arguments from the command line
    parser = argparse.ArgumentParser()

    parser.add_argument("-m", "--manifest_file", dest = "manifest_file",
        help = "Path to the manifest file, should have the same format as clipper.",
        default = "manifest_file.txt",
        action = "store",
        metavar = "\b")

    parser.add_argument("-o", "--output_dir", dest = "output_dir",
        help = "Path to the output directory.",
        default = "results",
        action = "store",
        metavar = "\b")

    # If we are using r1, this should be stranded I think
    parser.add_argument("-s", "--strandedness", dest = "strandedness",
        help = "The strandedness of the library. If you expect reads in the bam file to be on the same strand as the gene, this should be 'stranded', if it is the opposite strand, this should be 'reverse_stranded'",
        default = "stranded",
        action = "store",
        metavar = "\b")

    # If we are using r1, this should be stranded I think
    parser.add_argument("-c", "--cutoff", dest = "cutoff_length",
        help = "The cutoff for mapped length to be included, default is 10000",
        default = 1000,
        action = "store",
        metavar = "\b")



    args = parser.parse_args()

    return(args)

def _check_path(path):
    """
    Make sure a path exists and throw an error if it doesn't
    """

    if os.path.exists(path):
        return os.path.abspath(path)
    else:
        sys.exit("ERROR: " + path + " does not exist.")

def make_output_dir(path):
    """
    Make the output directory if it doesn't already exist
    """

    if not os.path.exists(path):
        os.makedirs(path)

###################
# Get input files #
###################

def make_input_dict(manifest_file, sample_dict):
    """
    Make a dictionary of input files
    """

    # Go through each entry in the manifest file
    with open(manifest_file, "r") as manifest_input:
        for line in manifest_input:
            save_name, bedfile, clip_bam, input_bam = get_input_files(line)

            sample_dict[save_name] = {"bedfile":bedfile, "clip_bam":clip_bam, "input_bam":input_bam}

            for i in sample_dict[save_name]:
                sample_dict[save_name][i] = _check_path(sample_dict[save_name][i])

def get_input_files(input_info):
    """
    Pull out bam files, get bed file name from bam file input
    """

    save_name, sample, celltype, clip_bam, input_bam = input_info.strip().split("\t")
    bedfile = re.sub("bam", "peaks.bed", clip_bam)

    return([save_name, bedfile, clip_bam, input_bam])

###################
# Count bam reads #
###################

def make_reads_dict(sample_dict, read_count_file):
    """
    Make a dictionary of the number of reads in each bam file. Save to a
    file if it doesn't already exist
    """

    reads_dict = {}
    if os.path.exists(read_count_file):
        bam_files = []
        for sample in sample_dict:
            bam_files.append(sample_dict[sample]["clip_bam"])
            bam_files.append(sample_dict[sample]["input_bam"])

        with open(read_count_file, "r") as count_file:
            for line in count_file:

                try:
                    # Get information from file
                    file_name, read_count = line.strip().split("\t")
                    int(read_count)
                except:
                    sys.exit("check manifest mapped read num file to make sure counts are correct")

                if file_name in bam_files:
                    bam_files.remove(file_name)
                    # Write to dictionary
                    reads_dict[file_name] = read_count

        if len(bam_files) >= 1:
            sys.exit("Not all files have been counted. Remove file named " +
                read_count_file + " and rerun.")

    else:
        with open(read_count_file, "w") as count_file:
            for sample in sample_dict:

                # Get file names
                clip_bam = sample_dict[sample]["clip_bam"]
                input_bam = sample_dict[sample]["input_bam"]

                # Count number of reads
                total_clip_reads = bam_read_count(clip_bam)
                total_input_reads = bam_read_count(input_bam)

                # Add to dictionary
                reads_dict[clip_bam] = total_clip_reads
                reads_dict[input_bam] = total_input_reads

                # Write to file
                count_file.write("{}\t{}\n".format(clip_bam, total_clip_reads))
                count_file.write("{}\t{}\n".format(input_bam, total_input_reads))

    return(reads_dict)


def bam_read_count(bamfile):
    """
    Return  the number of reads in a bam file
    """

    try:
        # Count reads using samtools
        read_count = subprocess.check_output(['samtools view -c -F 4 ' + bamfile], shell=True)
        read_count = int(read_count)
    except:
        sys.exit("There was a problem counting the number of reads in the bam file, " +
            bamfile + ". Check that samtools is installed and that the bamfile looks correct.")
    return (read_count)

#################
# Make bed dict #
################# 
def make_bed_dict(line):
    chromosome, start, end, gene, p_val, strand, middle_down, middle_up = line.strip().split("\t")
    bed_dict = {"chromosome": chromosome,
                "start": int(start),
                "end": int(end),
                "gene": gene,
                "p_val": p_val,
                "strand": strand,
                "middle_down": middle_down,
                "middle_up": middle_up}
    return(bed_dict)

####################
# Read bam section #
####################

def new_bam_file(bed_dict, bam_file, options, save_bam):
    chromosome = bed_dict["chromosome"]
    start = bed_dict["start"]
    end = bed_dict["end"]
    strand = bed_dict["strand"]

    total_reads = 0

    alignment_file = pysam.AlignmentFile(bam_file, "rb")

    write_file = pysam.AlignmentFile(save_bam, "wb", template = alignment_file)

    for read in alignment_file.fetch(chromosome, start, end):
        if options.strandedness == "stranded" and strand == "+":
            if not read.is_reverse:
                write_file.write(read)

        # If the library is stranded and the gene is - strand,
        # the read must map to the reverse strand
        elif options.strandedness == "stranded" and strand == "-":
            if read.is_reverse:
                write_file.write(read)

        # If the library is reverse stranded and the gene is + strand,
        # the read must map to the reverse strand
        elif options.strandedness == "reverse_stranded" and strand == "+":
            if read.is_reverse:
                write_file.write(read)

        # If the library is reverse stranded and the gene is - strand,
        # the read must map to the forward strand
        elif options.strandedness == "reverse_stranded" and strand == "-":
            if not read.is_reverse:
                write_file.write(read)

def count_bam_reads(bed_dict, bam_file, options):
    chromosome = bed_dict["chromosome"]
    start = bed_dict["start"]
    end = bed_dict["end"]
    strand = bed_dict["strand"]

    total_reads = 0

    alignment_file = pysam.AlignmentFile(bam_file, "rb")

    for read in alignment_file.fetch(chromosome, start, end):
        read_length = read.reference_length
        # print(read.reference_length)
        # print()

        if read_length > int(options.cutoff_length):
            continue

        # If the library is stranded and the gene is + strand,
        # the read must map to the forward strand
        if options.strandedness == "stranded" and strand == "+":
            if not read.is_reverse:
                total_reads += 1

        # If the library is stranded and the gene is - strand,
        # the read must map to the reverse strand
        elif options.strandedness == "stranded" and strand == "-":
            if read.is_reverse:
                total_reads += 1

        # If the library is reverse stranded and the gene is + strand,
        # the read must map to the reverse strand
        elif options.strandedness == "reverse_stranded" and strand == "+":
            if read.is_reverse:
                total_reads += 1

        # If the library is reverse stranded and the gene is - strand,
        # the read must map to the forward strand
        elif options.strandedness == "reverse_stranded" and strand == "-":
            if not read.is_reverse:
                total_reads += 1

    alignment_file.close()
    print(total_reads)
    return(total_reads)

"""
Based on comparing my numbers to the numbers from the perl script, I don't
need to change the index position.

some notes
From pysam fetch()
fetch one or more rows in a region using 0-based indexing. The region is specified
by reference, start and end. Alternatively, a samtools region string can be supplied.

From the perl script notes
# Keep in mind:
# STAR sam output (including output from 'samtools view .bam') is 1-based, closed ended
# Bed files are 0-based, open-ended

I'll need to think about this when counting... maybe I will need to shift the bed file
regions by 1?


From pysam faq Confusion might arise as different file formats might have different conventions.
For example, the SAM format is 1-based while the BAM format is 0-based. It is important to
remember that pysam will always conform to the python convention and translate to/from the file
format automatically.

The only exception is the region string in the fetch() and pileup() methods. This string follows
the convention of the samtools command line utilities. The same is true for any coordinates passed
to the samtools command utilities directly, such as pysam.mpileup().


pysam.AlignmentFile.fetch() returns all reads overlapping a region sorted by the first aligned
base in the reference sequence. Note that it will also return reads that are only partially
overlapping with the region. Thus the reads returned might span a region that is larger than
the one queried.
"""

###########################################
# Perform chi square or fisher exact test #
###########################################

def chi_square_or_fisher(peak_clip, peak_input, clip_total, input_total):
    a = int(peak_clip)
    b = int(clip_total) - a
    c = int(peak_input)
    d = int(input_total) - c

    tot = a + b + c + d
    expa = (a+c)*(a+b)/tot
    expb = (b+d)*(a+b)/tot
    expc = (a+c)*(c+d)/tot
    expd = (b+d)*(c+d)/tot

    # Check if fisher exact should be run
    if expa < 5 or expb < 5 or expc < 5 or expd < 5 or a < 5 or b < 5 or c < 5 or d < 5:
        fisher_test()
    elif expa >= 5 or expb >= 5 or expc >= 5 or expd >= 5:
        return_list = chi_square_test(a, b, c, d)

    return(return_list)

def chi_square_test(a, b, c, d):

    obs = np.array([[a, b], [c, d]])

    stat, p, dof, expected = chi2_contingency(obs)

    # P val
    p_value = abs(math.log10(p))

    # logfc
    logfc = math.log((a/b)/(c/d))/math.log(2)

    print(p_value)
    print(logfc)

    return([p_value, logfc])


if __name__ == "__main__":
    main()