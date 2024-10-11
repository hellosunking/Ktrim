<img src="https://github.com/hellosunking/hellosunking.github.io/blob/master/logos/Ktrim.png" width="50%" height="50%">

# Ktrim: an extra-fast and accurate adapter- and quality-trimmer for sequencing data
Version 1.6.0, Oct 2024<br />
Author: Kun Sun \(sunkun@szbl.ac.cn\)<br />
<br />
Distributed under the
[GNU General Public License v3.0 \(GPLv3\)](https://www.gnu.org/licenses/gpl-3.0.en.html "GPLv3")
for personal and academic usage only.<br />
For detailed information please refer to the license files under `license` directory.

---
## Release of version 1.6
The author is pleased to release version 1.6 of Ktrim with usability improvement.
This version added a "-R" option to output reads with adapters only. This could be useful for libraries
that (e.g., CLIP-seq, whose adapters are also included for built-in support).

## Major features of Ktrim
1. Fast, sensitive, and accurate
2. Supports both paired- and single-end data
3. Supports both Gzipped and plain text
4. Supports output to stdout to pipe with downstream software, e.g., aligners
5. Supports multi-threading for speed-up
6. Built-in support for common adapters; customized adapters are also supported

## Installation
`Ktrim` is written in `C++` for GNU Linux/Unix platforms. After uncompressing the source package, you can
find an executable file at `bin/ktrim` compiled using `g++ v4.8.5` and linked with `libz v1.2.7` for Linux
x64 system. If you could not run it (which is usually caused by low version of `libc++` or `libz` libraries)
or you want to build a version optimized for your system, you can re-compile the programs:
```
user@linux$ make clean && make
```

## Run Ktrim
The main program is `ktrim` under `bin/` directory. You can add its path to your `~/.bashrc` file under
the `PATH` variable to call it from anywhere; or you can run the following command to add it to your
current session:
```
user@linux$ PATH=$PATH:$PWD/bin/
```

You can also add `ktrim` to your system to call it from anywhere and share with other users (requires
root privilege):
```
user@linux# make install
```

Call `ktrim` without any parameters to see the usage (or use '-h' option):
```
Usage: Ktrim [options] -f fq.list {-1/-U Read1.fq [-2 Read2.fq ]} -o out.prefix

Author : Kun Sun (sunkun@szbl.ac.cn)
Version: 1.6.0 (Oct 2024)

Ktrim is designed to perform adapter- and quality-trimming of FASTQ files.

Compulsory parameters:

  -f fq.list      Specify the path to a file containing path to read 1/2 fastq files

OR you can specify the fastq files directly:

  -1/-U Read1.fq  Specify the path to the files containing read 1
                  If your data is Paired-end, use '-1' and specify read 2 files using '-2' option
                  Note that if '-U' is used, specification of '-2' is invalid
                  If you have multiple files for your sample, use ',' to separate them
                  Note that Gzip-compressed files (with .gz suffix) are also supported

  -o out.prefix   Specify the prefix of the output files
                  Note that output files include trimmed reads in FASTQ format and statistics

Optional parameters:

  -2 Read2.fq     Specify the path to the file containing read 2
                  Use this parameter if your data is generated in paired-end mode
                  If you have multiple files for your sample, use ',' to separate them
                  and make sure that all the files are well paired in '-1' and '-2' options

  -c              Write the trimming results to stdout (default: not set)
                  Note that the interleaved fastq format will be used for paired-end data.
  -R                Only output reads with adapter (default: not set)
  -t threads      Specify how many threads should be used (default: 6)
                  You can set '-t' to 0 to use all threads (automatically detected)
                  2-8 threads are recommended, as more threads would not benefit the performance

  -p phred-base   Specify the baseline of the phred score (default: 33)
  -q score        The minimum quality score to keep the cycle (default: 20)
                  Note that 20 means 1% error rate, 30 means 0.1% error rate in Phred
  -w window       Set the window size for quality check (default: 5)
                  Ktrim will stop when all the bases in a consecutive window pass the quality threshold

                  Phred 33 ('!') and Phred 64 ('@') are the most widely used scoring system
                  while quality scores start from 35 ('#') in the FASTQ files is also common

  -s size         Minimum read size to be kept after trimming (default: 36)

  -k kit          Specify the sequencing kit to use built-in adapters
                  Currently supports 'Illumina' (default), 'Nextera', 'Transposase' and 'BGI'
  -a sequence     Specify the adapter sequence in read 1
  -b sequence     Specify the adapter sequence in read 2
                  If '-a' is set while '-b' is not, I will assume that read 1 and 2 use same adapter
                  Note that '-k' option has a higher priority (when set, '-a'/'-b' will be ignored)

  -m proportion   Set the proportion of mismatches allowed during index and sequence comparison
                  Default: 0.125 (i.e., 1/8 of compared base pairs)

  -h              Show this help information and quit
  -v              Show the software version and quit

Please refer to README.md file for more information (e.g., setting adapters).
Citation: Sun K. Bioinformatics 2020 Jun 1; 36(11):3561-3562. (PMID: 32159761)

Ktrim: extra-fast and accurate adapter- and quality-trimmer.
```

Please note that from version 1.2.0, Ktrim supports Gzip-compressed files as input. If you have multiple
lanes of FASTQ files, Ktrim even supports that some lanes are compressed while others are in plain text, as
long as READ1 and READ2 are the same (either both compressed or plain text) for each lane of paired-end data.

`Ktrim` contains built-in adapter sequences used by Illumina TruSeq kits, Nextera kits, Nextera transposase
adapters and BGI sequencing kits within the package. However, customized adapter sequences are also allowed
by setting '-a' (for read 1) and '-b' (for read 2; if it is the same as read 1, you can left it blank)
options. You may need to refer to the manual of your library preparation kit for the adapter sequences.
Note that in the current version of `Ktrim`, only 1 pair of adapters is allowed.

Here are the built-in adapter sequences (the copyright should belong to the corresponding companies):

```
Illumina TruSeq kits:
AGATCGGAAGAGC (for both read 1 and read 2)

Nextera kits (suitable for ATAC-seq, Cut & tag data):
CTGTCTCTTATACACATCT (for both read 1 and read 2)

BGI adapters:
Read 1: AAGTCGGAGGCCAAGCGGTC
Read 2: AAGTCGGATCGTAGCCATGT

CLIP-seq adapters:
Read 1: TGGAATTCTCGGGTGCCAAGG
Read 2: GATCGTCGGACTGTAGAACTCTGAAC
```

### Example 1
Your data is generated using Illumina TruSeq kit in Single-end mode, then you can run:
```
user@linux$ ktrim -U /path/to/read1.fq -o /path/to/output/dir
```

### Example 2
Your data is generated using a kit with customized adapters; your data is composed of 3 lanes in Paired-end
mode and you have prepared a `file.list` to record the paths as follows:
```
/path/to/lane1.read1.fq.gz	/path/to/lane1.read2.fq.gz
/path/to/lane2.read1.fq.gz	/path/to/lane2.read2.fq.gz
/path/to/lane3.read1.fq	/path/to/lane3.read2.fq
```
in addition, your Phred scoring system starts from 64; you want to keep the high quality (Phred score >=30)
bases and reads longer than 50 bp after trimming; and you want to use 8 threads to speed-up the analysis,
then you can run:
```
user@linux$ ktrim -f file.list -t 8 -p 64 -q 30 -s 50 -o /path/to/output/dir \
                  -a READ1_ADAPTER_SEQUENCE -b READ2_ADAPTER_SEQUENCE
```

## Outputs explanation
`Ktrim` outputs the trimmed reads in FASTQ format and key statistics (e.g., the numbers of reads that
contains adapters and the number of reads in the trimmed files).

## Testing dataset and benchmark evaluation
Under the `testing_dataset/` directory, a script named `simu.reads.pl` is provided to generate *in silico*
reads for testing purpose only. **Note that the results in the paper is based on the data generated by this
script.** Another script `check.accuracy.pl` is designed to evaluate the accuracies of the trimming tools.

Please refer to Supplementary Method in the paper for reproducing the results (using Ktrim v1.1.0).

## Citation
When referencing, please cite "Sun K: **Ktrim: an extra-fast and accurate adapter- and quality-trimmer
for sequencing data.** *Bioinformatics* 2020 Jun 1; 36(11):3561-3562."
[PubMed](https://www.ncbi.nlm.nih.gov/pubmed/32159761 "Ktrim@PubMed")
[Full Text](https://doi.org/10.1093/bioinformatics/btaa171 "Full text on Bioinformatics journal")

---
Please send bug reports to Kun Sun \(sunkun@szbl.ac.cn\).<br />
Ktrim is freely available at
[https://github.com/hellosunking/Ktrim/](https://github.com/hellosunking/Ktrim/ "Ktrim @ Github").

