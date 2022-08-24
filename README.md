# eClip normalization
Scripts to replicate the normalization performed after clipper. The original scripts are [here]( https://github.com/YeoLab/gscripts/tree/1.1/perl_scripts)

## Requirements

`numpy`

`scipy`

`pysam`

`samtools` - must be loaded in your working environment

`python` v3.X

## Description

The overall goal of the `peak_input_normalization.py` script is to normalize read counts in peaks identified by `clipper` and then find p-values using `chi-squared` or `fisher-exact` tests. Reads are only kept if part of the mapped read overlaps a peak (so a spliced read will not be counted for a peak that is in the intron). The output is the same as the output for `peak_input_normalization.pl`. It includes columns with information about the peaks (`chromosome`, `start`, `end`, `strand`), the p-value (`pvalue`), the log of the p-value (`log10p`), the log fold change (`log2fc`), and the counts for the clip and input so you can check or alter the p-value calculation (`peak_counts_clip`, `peak_counts_input`, `total_clip_counts`, `total_input_counts`).

This script follows the general logic of the original `perl` script `peak_input_normalization.pl`, but has been rewritten in python to allow for some modifications.

1. The most important modification is the inclusion of a `--strand` argument. While analyzing our own `eCLIP` data, we realized that the protocol we used did not have the same strand specificity as the original published protocol (the original was reverse stranded, the new was stranded). This meant that when running `clipper` using Read 2 as the input (as described in their analysis pipeline), we saw almost no mapping to genes, while using Read 1 worked. Although `clipper` does not include a strand argument, you can get around this by using only reads from Read 1.
2. Instead of using `perl` and parsing the cigar string, I used `pysam` to parse the cigar string for me.
3. When calculating log fold change and p-values, I largely followed the logic in the original `perl` script, but they add a pseudo count of 1 only to the number of input reads before finding logfc and running a statistical test. I have added the option to do this using `flavor`, but the default is to add a pseudocount of 1 before finding the log fold change, but not altering any of the counts to determine the p-value.

## Usage

We assume that you have largely followed the protocol to use `clipper` found [here](https://www.encodeproject.org/documents/3b1b2762-269a-4978-902e-0e1f91615782/@@download/attachment/eCLIP_analysisSOP_v2.0.pdf)

This script is meant to replace the analysis found in the section labeled "Peak normalization vs SMInput and reproducible peak / IDR analysis"

**Please note, we have found that the strandedness of the protocol has changed and you may need to run `clipper` with read 1 instead of read 2 as described in the linked document**

This is a stand alone script, you only need to download `peak_input_normalization.py` found in `src/scripts`. Other files in the `scripts` folder were used for testing my pipeline against the original `perl` script.

To run:

```bash
python /path/to/peak_input_normalization.py -m /path/to/manifest
```

There are a few other optional arguments:

```
usage: peak_input_normalization.py [-h] [-m] [-o] [-s] [-f]

optional arguments:
  -h, --help            show this help message and exit
  -m, --manifest_file
                        Path to the manifest file, should have the same format
                        as clipper.
  -o, --output_dir  Path to the output directory. Default is 'results'.
                        This will be made if it doesn't exist
  -s, --strandedness
                        The strandedness of the library. If you expect reads
                        in the bam file to be onthe same strand as the gene,
                        this should be 'stranded', if it is the opposite
                        strand, this should be 'reverse_stranded'. If the
                        library is unstranded this should be 'unstranded'. The
                        default is 'stranded'
  -f, --flavor      If input reads should be counted like the original
                        perl script (always add 1) set to 'perl_script',
                        otherwise set to 'default'
```

### Further descriptions of arguments

* `-m` `--manifest_file`: This was written to run almost identically to the orignal perl script. The manifest file is the same as the manifest file originally used. This file should be a tab delimited file with 5 columns. Examples of these manifest files are in the directory `manifest_files`
  * Column 1 - sample number, should be a unique value for each line in the file
  * Column 2 - sample name, this could be the name of the protein
  * Column 3 - cell type
  * Column 4 - Bam file (for just read 1) containing the clip reads. The peak file with the same name (but `.bam` replaced with `.peaks.bed`) is expected to be in the same folder.
  * Column 5 - Bam file (for just read 1) containing the input reads. The peak file with the same name (but `.bam` replaced with `.peaks.bed`) is expected to be in the same folder.
* `-o` `--output_dir`: Path to the output directory. The directory will be made. Default is `results`
* `-s` `--strandeness`: If it is forwarded stranded (`stranded`, read maps to the same strand as the gene), reverse stranded (`reverse_stranded`, read maps to the opposite strand as the gene) or unstranded (`unstranded`).
* `-f` `--flavor`: If you want the exact results from the `perl` script, choose `perl_script`. This adds a psudocount of 1 to the input in all cases. Otherwise, the default is to only use a pseudocount to determine the log fold changes.

Example scripts to run the files are in `src/shell_scripts`

### Saved files
The first time through the script, the number of reads in each bam file will be counted. This is a slow process, so the number of reads is saved with the same name as your manifest file appended with `.mappped_read_num`. If this file exists, the counts from it will be used instead of recounting.

### Some potential issues that need testing
This was developed using read 1 to run clipper (because our protocol was stranded with read 1 mapping to the forward stand). The output of clipper is a bed file with the same name but `.peaks.bed` in place of the `.bam`. This should work if you use read 2 to call peaks with clipper, but then you either need to use read 1 (and rename the files) or use read 2 and the stranded argument will be opposite of what it would be if you use read 1. This has not been tested, but I plan to test it in the future. If your data was generated using a reverse stranded protocol, you should be able to use the original clipper script.