import numpy as np
from scipy.stats import chi2_contingency
from scipy.stats import fisher_exact
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
        # Empty dictionary for compressing later
        results_dict = defaultdict(list)
        save_file_all = os.path.join(options.output_dir, sample +
            "_01.basedon_001_01.peaks.l2inputnormnew.bed")
        bed_file = sample_dict[sample]["bedfile"]
        input_bam = sample_dict[sample]["input_bam"]
        clip_bam = sample_dict[sample]["clip_bam"]

        with open(save_file_all, "w") as write_file, open(bed_file, "r") as in_bed:
            for line in in_bed:
                # print()
                bed_dict = make_bed_dict(line)

                # Find number of reads in peaks and overall
                read_count_clip = count_bam_reads(bed_dict, clip_bam, options)
                read_count_input = count_bam_reads(bed_dict, input_bam, options)

                if options.flavor == "perl_script":
                    # To be consistent with the perl script use:
                    read_count_input += 1
                elif options.flavor == "default":
                    # Otherwise use
                    if read_count_input == 0:
                        read_count_input = 1
                else:
                    sys.exit("unrecognized 'flavor' argument! Use 'perl_script' or 'default'")

                # print(read_count_input)
                total_clip_reads = total_reads_dict[clip_bam]
                total_input_reads = total_reads_dict[input_bam]

                # Run tests
                p_val_log, p_val, logfc = chi_square_or_fisher(read_count_clip, read_count_input, total_clip_reads, total_input_reads)
                
                # Write to a file
                write_file.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(bed_dict["chromosome"], bed_dict["start"], bed_dict["end"], p_val_log, logfc, bed_dict["strand"]))
                
                # Save all output to a dictionary for compressing
                bed_dict["p_value"] = p_val
                bed_dict["log"] = logfc
                bed_dict["p_val_log"] = p_val_log
                save_val = bed_dict["chromosome"] + "_" + bed_dict["strand"]

                results_dict[save_val].append(bed_dict)

                # save_clip = os.path.join("test_bams", clip_bam.split("/")[-1])
                # save_input  = os.path.join("test_bams", input_bam.split("/")[-1])
                # new_bam_file(bed_dict, clip_bam, options, save_clip)
                # new_bam_file(bed_dict, input_bam, options, save_input)

        # print(results_dict)
        new_results_dict = compress_peaks(results_dict)

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
        help = "Path to the output directory. Default is 'results'. This will be made if it doesn't exist",
        default = "results",
        action = "store",
        metavar = "\b")

    # If we are using r1, this should be stranded I think
    parser.add_argument("-s", "--strandedness", dest = "strandedness",
        help = ("The strandedness of the library. If you expect reads in the bam file to be on" + 
            "the same strand as the gene, this should be 'stranded', if it is the opposite strand," +
            " this should be 'reverse_stranded'. If the library is unstranded this should be " +
            "'unstranded'. The default is 'stranded'"),
        default = "stranded",
        action = "store",
        metavar = "\b")

    # If we are using r1, this should be stranded I think
    parser.add_argument("-c", "--cutoff", dest = "cutoff_length",
        help = "The cutoff for mapped length to be included, default is 1000",
        default = 1000,
        action = "store",
        metavar = "\b")


    # If we are using r1, this should be stranded I think
    parser.add_argument("-f", "--flavor", dest = "flavor",
        help = ("If input reads should be counted like the original perl script (always add 1) " +
                "set to 'perl_script', otherwise set to 'default'"),
        default = "default",
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

    Returns - No return, but the sample dictionary is updated with paths
    to all files needed for analysis.
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

    Returns - a list of file names based on the input manifest file.
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
    file if it doesn't already exist. If the file already exists, the
    read number is extracted from the file, if it doesn't, samtools is
    used to count reads.

    Returns - a dictionary of each file and the number of reads.
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
                read_count_file + " and rerun.\n")

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
    Uses samtools to count the number of reads in a bam file. Requires samtools to be
    loaded

    Returns - the number of reads in a bam file.
    """

    try:
        # Count reads using samtools
        read_count = subprocess.check_output(['samtools view -c -F 4 ' + bamfile], shell=True)
        read_count = int(read_count)
    except:
        sys.exit("There was a problem counting the number of reads in the bam file located: " +
            bamfile + ".\n Check that samtools is installed and that the bamfile looks correct.\n")
    return (read_count)

#################
# Make bed dict #
################# 
def make_bed_dict(line):
    """
    Makes a dictionary of information about each line in a bed file.

    Returns - a dictionary with information from the bed file that can be
    easily accessed.
    """
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
    """
    Generates a new bam file consisting only of the reads that mapped to a 
    specific region. This function was written so that I could look at how
    the perl script was processing specific reads and should not be needed
    with the final script. It is a helpful function though, so I will keep
    it.

    Returns - Nothing, writes a new bam file consisting only of the reads
    of interest.
    """

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
    """
    Uses pysam.fetch to identify all reads that map to a region of the genome.
    The reads will be counted according to the strandedness of the library and
    the read used.

    There is currently also a cutoff for the length of the mapped read because
    we were seeing many reads that had "introns" in the cigar string that were
    thousands of base pairs long and the start and end of the read were not in
    the same gene.

    This currently largely counts the same reads as the original perl script, 
    but I have not completely been able to replicate it. In the original perl
    script, they split a read at any intron and determined if either of the
    reads from the split read were in the region. All of these intron reads 
    ended up being thrown out in my experience.

    To attempt to correctly count reads like the perl script, I skip reads
    that don't have any positions mapped within the read. I do this using
    pysams get_reference_positions. To save time and computing power, I 
    only look though positions until I find one position that is between the
    start and end of the peak. In testing, this mostly had the same number 
    as the perl script. The one concern is that the perl script said there
    were 25 counts in the input and I see 24. Interestingly, when I manually
    count the output from the perl script I also see 24.

    For the index, in the perl script, they noted that bed files are 1-based
    and bam files are 0-based. Pysam automatically accounts for this and I
    compared all of my regions to the regions in the perl script and found
    that they lined up perfectly without me changing the index position.

    Returns - the total reads that mapped to the appropriate strand in a region
    of the bam file.
    """

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

        # This does a check to make sure that a position mapped to
        # the peak. This is to fix those reads where the intron
        # is thousands of bp long and skips the peak.
        keep_read = False
        for position in read.get_reference_positions():
            if position >= start and position <= end:
                keep_read = True
                continue

        # If no positions mapped to the peak, move on to the next read
        if not keep_read:
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

        # If the library is not stranded, count all reads
        elif options.strandedness == "unstranded":
            total_reads += 1

    alignment_file.close()
    # print(chromosome + "_" + str(start) + "_" + str(end))
    # print(total_reads)
    return(total_reads)

###########################################
# Perform chi square or fisher exact test #
###########################################

def chi_square_or_fisher(peak_clip, peak_input, clip_total, input_total):
    """
    Decides if a chi-square or fiser exact test. The fisher exact test will be
    run if any of the values or expected values will be less than 5. Once this
    decision is made, it passes the values to either the fishers exact test or 
    chi-sequare test.

    It also calculates a log2 fold change.

    It will calculate the p value as 1 if the input is higher than the clip.

    Returns: a list consisting of the p-value and the logfc calculated.
    """
    a = int(peak_clip)
    b = int(clip_total) - a
    c = int(peak_input)
    d = int(input_total) - c

    tot = a + b + c + d
    expa = (a + c) * (a + b) / tot
    expb = (b + d) * (a + b) / tot
    expc = (a + c) * (c + d) / tot
    expd = (b + d) * (c + d) / tot

    # Make a contengency table
    obs = np.array([[a, b], [c, d]])

    # logfc
    logfc = math.log2((a / int(clip_total)) / (c / int(input_total)))

    # Set p value to 1 and log p value to 0 if input is higher than clip
    if logfc < 0:
        return_list = [0, 1]
    # Check if fisher exact should be run
    elif expa < 5 or expb < 5 or expc < 5 or expd < 5 or a < 5 or b < 5 or c < 5 or d < 5:
        return_list = fisher_test(obs)
    elif expa >= 5 or expb >= 5 or expc >= 5 or expd >= 5:
        return_list = chi_square_test(obs)
    else:
        sys.exit("Unclear if chi squared or fishers test should be done.")

    return_list.append(logfc)

    return(return_list)

def chi_square_test(obs):
    """
    Runs a chi square test given the number of reads in the peak for both input
    and clip and the number of reads outside of the peak for both input and clip.
    These read numbers should be proveded as a contengency table.

    Returns: A list consisting of the p value and log p value calculated
    """

    stat, p, dof, expected = chi2_contingency(obs)

    # P val
    p_value = abs(math.log10(p))

    return([p_value, p])

def fisher_test(obs):
    """
    Runs a fisher exact test given the number of reads in the peak for both
    input and clip and the number of reads outside of the peak for both input
    and clip. These read numbers should be provided as a contengency table.

    Returns: A list consisting of the p value and log p value calculated
    """

    odds, p = fisher_exact(obs)

    # P val
    p_value = abs(math.log10(p))

    return([p_value, p])

##################
# Compress peaks #
##################

def compress_peaks(results_dict):
    """
    This compresses peaks so that only one peak overlapping a region is in the final
    output file. The winning peak is based on the log fold enrichment (highest
    enrichemnt wins)

    Return - a new dictionary consisting of only the winning sections.
    """
    new_results_dict = defaultdict(list)

    for chromosome_strand in results_dict:
        peak_list = results_dict[chromosome_strand]

        keep_index = []
        index = 0
        for peak in peak_list:
            # Default is to keep unless otherwise decided
            add = True
            for index_val in keep_index:
                remove_index = False
                keep_peak = peak_list[index_val]

                # If the peaks overlap, keep the best
                # Start of read is inside existing
                start_inside = peak["start"] >= keep_peak["start"] and peak["start"] <= keep_peak["end"]
                
                # End of read is inside existing
                end_inside = peak["end"] >= keep_peak["start"] and peak["end"] <= keep_peak["end"]
                
                # Full read is inside existing - I think this isn't necessary, test later
                shorter_inside = peak["start"] >= keep_peak["start"] and peak["end"] <= keep_peak["end"]
                
                # Existing read is inside new
                longer_inside = peak["start"] <= keep_peak["start"] and peak["end"] >= keep_peak["end"]
                
                if start_inside or end_inside or shorter_inside or longer_inside:
                    if peak["log"] > keep_peak["log"]:
                        remove_index = True
                        add = True

                    elif peak["log"] == keep_peak["log"]:
                        add = True
                        remove_index = False

                    else:
                        remove_index = False
                        add = False

                if remove_index:
                    keep_index.remove(index_val)

            if add:
                keep_index.append(index)

            index += 1

        chromosome_list = []

        # Add the indexes that were kept
        for index_num in keep_index:
            chromosome_list.append(peak_list[index_num])

        # Add to final dictionary
        new_results_dict[chromosome_strand] = chromosome_list




if __name__ == "__main__":
    main()