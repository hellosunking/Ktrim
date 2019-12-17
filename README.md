
# Ktrim: an extra-fast and accurate adapter- and quality-trimmmer for sequencing data
Version 1.0.0, Dec 2019<br />
Author: Kun Sun \(sunkun@szbl.ac.cn\)<br />
<br />
Distributed under the
[GNU General Public License v3.0 \(GPLv3\)](https://www.gnu.org/licenses/gpl-3.0.en.html "GPLv3")
for personal and academic usage only.<br />
For detailed information please refer to the license files under `license` directory.

---

## Installation
`Ktrim` is written in `C++` for GNU Linux/Unix platformw. After uncompressing the source package, you
can find a pre-compiled executable under `bin/` directory for Linux x86_64 system compiled using
`G++` v7.4.0. If you could not run it (which usually caused by low version of `libc++` library) or you
want to build an optimized version for your machine, you can re-compile the programs:
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

You can also add `ktrim` to the system to call it from anywhere and share with other users (requires
root privilege):
```
user@linux# make install
```

Call `ktrim` without any parameters to see the usage (or use '-h' option):
```
Usage: Ktrim [options] -1/-U Read1.fq [ -2 Read2.fq ] -o out.prefix

Author : Kun Sun (sunkun@szbl.ac.cn)
Version: 1.0.0 (Dec 2019)

Ktrim is designed to perform adapter- and quality-trimming of FASTQ files.

Compulsory parameters:

  -1/-U Read1.fq   Specify the path to the files containing read 1
                   If your data is Paired-end, use '-1' and specify read 2 files using '-2' option
                   Note that if '-U' is used, specification of '-2' is invalid
                   If you have multiple files for your sample, use ',' to separate them

  -o out.prefix    Specify the prefix of the output files
                   Note that output files include trimmed reads in FASTQ format and statistics

Optional parameters:

  -2 Read2.fq      Specify the path to the file containing read 2
                   Use this parameter if your data is generated in paired-end mode
                   If you have multiple files for your sample, use ',' to separate them
                   and make sure that all the files are well paired in '-1' and '-2' options

  -t threads       Specify how many threads should be used (default: 1, single-thread)
                   You can set '-t' to 0 to use all threads (automatically detected)

  -p phred-base    Specify the baseline of the phred score (default: 33)
  -q score         The minimum quality score to keep the cycle (default: 20)
                   Note that 20 means 1% error rate, 30 means 0.1% error rate in Phred

                   Phred 33 ('!') and Phred 64 ('@') are the most widely used scoring system
                   while quality scores start from 35 ('#') in the FASTQ files is also common

  -s size          Minimum read size to be kept for alignment (default: 36)

  -k kit           Specify the sequencing kit to use built-in adapters
                   Currently supports 'Illumina' (default), 'Nextera', 'Transposase' and 'BGI'
  -a sequence      Specify the adapter sequence in read 1
  -b sequence      Specify the adapter sequence in read 2
                   If '-a' is set while '-b' is not, I will assume that read 1 and 2 use same adapter
                   Note that '-k' option has a higher priority (when set, '-a'/'-b' will be ignored)

  -h/--help        Show this help information and quit
  -v/--version     Show the software version and quit

Please refer to README.md file for more information (e.g., setting adapters).

Ktrim: extra-fast and accurate adapter- and quality-trimmer.
```

`Ktrim` contains built-in adapter sequences used by Illumina TruSeq kits, Nextera kits, Nextera transposase
adapters and BGI sequencing kits within the package. However, customized adapter sequences are also allowed
by setting '-a' (for read 1) and '-b' (for read 2; if it is the same as read 1, you can left it blank)
options. You may need to refer to the manual of your library preparation kit for the adapter sequences.
Note that in the current version of `Ktrim`, only 1 pair of adapters is allowed.

Here are the built-in adapter sequnces:

```
Illumina TruSeq kits:
AGATCGGAAGAGC (for both read 1 and read 2)

Nextera kits:
CTGTCTCTTATACACATCT (for both read 1 and read 2)

Nextera transposase adapters:
Read 1: TCGTCGGCAGCGTC
Read 2: GTCTCGTGGGCTCG

BGI adapters:
Read 1: AAGTCGGAGGCCAAGCGGTC
Read 2: AAGTCGGATCGTAGCCATGT
```

### Example 1
Your data is generated using Illumina TruSeq kits in Single-End mode, then you can run:
```
user@linux$ ktrim -U /path/to/read1.fq -o /path/to/output/dir
```

### Example 2
Your data is generated using a kit with customized adapters; your data is composed of 3 lanes in Paired-end
mode and uses Phred 35 scoring system; you want to keep the high quality (Phred score >=30) bases and reads
longer than 50 bp after trimming; and you want to use 4 threads to speed-up the analysis, then you can run:
```
user@linux$ ktrim -1 /path/to/lane1.read1.fq,/path/to/lane2.read1.fq,/path/to/lane3.read1.fq \
                  -2 /path/to/lane1.read2.fq,/path/to/lane2.read2.fq,/path/to/lane3.read2.fq \
                  -t 4 -p 35 -q 30 -s 50 -o /path/to/output/dir \
                  -a SEQUENCE_READ1_ADAPTER -b SEQUENCE_READ2_ADAPTER
```

## Testing dataset
Under the `testing_dataset/` directory, a script named `simu.reads.pl` is provided to generated *in silico*
reads for testing purpose only. Note that the comparison results used in the paper is based on the data
generated by this script. Another script `check.accuracy.pl` is designed to evaluate the accuracy of the
trimming tool.

The following are the information and codes that you can use to reproduce the results in the paper
(you may need to modify the codes by updating the paths of `Ktrim`, `Trim Galore` and `Trimmomatic`):

```
## Hardware:
## CPU: Intel(R) Core(TM) i5-6200U CPU @ 2.30GHz
## Memory: 8 GB DDR4 2133MHz
## Storage: 256 GB PCIe SSD
## 
## Operating System:
## OS version: Ubuntu 18.04.1 x86_64
## Linux kernal version: 5.0.0-37-generic
##
## Software versions:
## Trimmomatic: version 0.39; Java: openjdk version 11.0.4
## Trim Galore: version 0.6.5; cutadpat version 2.7; Perl version v5.26.1; Python version 3.6.9
## Ktrim; version 1.0.0; compiled using G++ version 7.4.0

## Note that 1 and 4 threads are both used for performance comparisons
THREAD=1

## generate in silico reads with adapters
## adapter sequence used: CTGTCTCTTATACACATCT, AGATGTGTATAAGAGACAG
## Output files: inSilicoReads.read1.fq, inSilicoReads.read2.fq
perl simu.reads.pl inSilicoReads

############################# Paired-end Mode #####################################
## run the 3 trimmers
## the 'adapter.fa' file is also stored under the "testing_dataset" directory
time java -jar trimmomatic-0.39.jar PE -threads $THREAD inSilicoReads.read1.fq inSilicoReads.read2.fq \
               Trimmomatic.read1.fq Trimmomatic.read1.unpaired.fq \
               Trimmomatic.read2.fq Trimmomatic.read2.unpaired.fq \
               ILLUMINACLIP:adapter.fa:2:30:10 MINLEN:36

time trim_galore inSilicoReads.read1.fq inSilicoReads.read2.fq \
                 --base TrimGalore --paired --length 36 --dont_gzip -j $THREAD \
                 -a CTGTCTCTTATACACATCT -a2 AGATGTGTATAAGAGACAG

time ktrim -1 inSilicoReads.read1.fq -2 inSilicoReads.read2.fq -o Ktrim \
           -t $THREAD -s 36 -a CTGTCTCTTATACACATCT -b AGATGTGTATAAGAGACAG

## compare the accuracies of the 3 trimmers
perl check.accuracy.pl inSilicoReads.read1.fq Trimmomatic.read1.fq
perl check.accuracy.pl inSilicoReads.read1.fq TrimGalore_val_1.fq
perl check.accuracy.pl inSilicoReads.read1.fq Ktrim.read1.fq


############################# Single-end Mode #####################################
time java -jar trimmomatic-0.39.jar SE -threads $THREAD inSilicoReads.read1.fq \
               Trimmomatic.SE.read1.fq ILLUMINACLIP:adapter.fa:2:30:10 MINLEN:36

time trim_galore inSilicoReads.read1.fq --base TrimGalore.SE \
                 --length 36 --dont_gzip -j $THREAD -a CTGTCTCTTATACACATCT

time ktrim -U inSilicoReads.read1.fq -o Ktrim.SE -t $THREAD -s 36 -a CTGTCTCTTATACACATCT

perl check.accuracy.pl inSilicoReads.read1.fq Trimmomatic.SE.read1.fq
perl check.accuracy.pl inSilicoReads.read1.fq TrimGalore.SE_trimmed.fq
perl check.accuracy.pl inSilicoReads.read1.fq Ktrim.SE.read1.fq
```

## Outputs explanation
`Ktrim` outputs the trimmed reads in FASTQ format and key statistics, including the number of reads with
adapters and the number of reads that are kept in the trimmed files.

---
Please send bug reports to Kun Sun \(sunkun@szbl.ac.cn\).<br />
`Ktrim` is freely available at
[https://github.com/hellosunking/Ktrim/](https://github.com/hellosunking/Ktrim/ "Ktrim @ Github").

