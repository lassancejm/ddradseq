# **ddRadSeq**

**Double Digest RADseq pipeline**

## * * *


## Table of contents

[[TOC]]

## * * *


## Overview

The ddRadSeq pipeline comprises a control script and three "worker" programs. Each worker program performs a separate step in the data processing. The three worker programs are

1. **parse_pool**: this program parses the raw fastQ files according to the standard Illumina multiplexing index sequences.

2. **trim_barcode**: the second program inputs paired fastQ files from each pool and sorts the forward sequences according to the custom 5' barcode, the barcode is then trimmed from the sequence and output.

3. **trim_3prime**: this final worker program checks for instances of the barcode on the 3' end of the reverse reads and, if they are present, trims them.

The behavior of the pipeline and worker programs is controlled by a python script called **ddradseq.py**. This document will provide an introduction to using the pipeline.

## The ddradseq.py control script

Running the script with the "--help" flag produces the following message:

* * *


./ddradseq.py --help

usage: use "python ./ddradseq.py --help" for more information

Control script for ddRadSeq pipeline

optional arguments:

  -h, --help            show this help message and exit

  -s N, --stage N       Start pipeline at a particular stage

                        0: Start pipeline at parse_pool and run to end [default]

                        1: Start pipeline at trim_barcode and run to end

                        2: Start pipeline at trim_3prime and run to end

                        3: Start pipeline at bwa and run to end

                        4: Run only parse_pool stage and exit

                        5: Run only trim_barcode stage and exit

                        6: Run only trim_3prime state and exit

                        7: Run only bwa stage and exit

  -d DIR, --dir DIR     Directory containing input fastQ files

  -o DIR, --out DIR     Parent directory for output files

  -r FILE, --ref FILE   Name of file with reference sequence for read mapping

  -p FILE, --pool FILE  CSV file with Illumina multiplex indices and pool names

  -b FILE, --barcode FILE

                        CSV file with custom barcodes and individual identifiers

  -c, --cluster         If this switch is present, script runs in cluster (slurm) mode

  -m N, --map N         Number of threads available for bwa read mapping

  -t N, --threads N     Number of threads available for concurrency

  -v, --version         show program's version number and exit

For more detailed information on how to run this script, see the overview document.

* * *


One of the most crucial input parameters for the control script is specified by the ‘-d’ or ‘--dir’ flag. This option takes as an argument the full path to the directory that contains the raw fastQ input files. For example, let us imagine that the directory "/opt/data/fastq/" contains all of our raw read data. The script will search both the “/opt/data/fastq/” directory, as well as all of the subdirectories contained within it, for example “/opt/data/fastq/run1/” and “/opt/data/fastq/run2/”. The script will consider files as input if the file ends with the extension “*.fastq.gz”. Thus, to insure that all of your data files are identified by the script, be sure that the files are gzipped and end with the proper extension. This is a mandatory option for the script to run.

Another import parameter is specified by the ‘-o’ or ‘--out’ flag. This option takes an argument that tells the script where to write all of the output files. For example, if one wanted to write all of the output files to the user’s home directory, one could type ‘--out /home/dgarriga’. There will be four subdirectories that contain both intermediate and final output files. The ‘pool/’ directory will contain the output from the first stage of the pipeline (the **parse_pool** step). In the case of our example, the ‘/home/dgarriga/pool’ directory will contain all of the fastQ files, sorted by the Illumina multiplexing index sequence. These output files represent the cumulative output across all input files (e.g., separate lanes or flow cells). Similarly, the second stage of the pipeline will produce output in a subdirectory called ‘trim/’, the third stage will produce output in a subdirectory called ‘final/’ and the last bwa mapping step will produce *.bam files in a subdirectory called ‘bam/’.

Certain stages of the pipeline will require additional input files as well. The control script is fairly smart about checking which stage of the pipeline you want to run and asking for these required auxiliary files. For example, the **parse_pool** stage requires the ‘-p’ or ‘--pool’ option to be specified. The expected argument is a file that contains a comma-separated value (CSV) list of the pool identifier associated with each Illumina index. For example, a file called ‘illumdb.csv’ may contain the following entries:

* * *


TTAGGC,CamembertA

TGACCA,CamembertB

* * *


Note that no header row is provided. Likewise, the **trim_barcode** step requires the ‘-b’ or ‘--barcode’ option to be specified, this should refer to another CSV file with the individual identifier for each custom barcode adapter. The format of this file should be the same as that specified with the ‘-p’ flag. And finally, if the **bwa** mapping step is being run, the ‘-r’ or ‘--ref’ flag should be used to input the reference fastA genome sequence.

The pipeline can optionally be run from any stage. This behavior is governed by the ‘-s’ or ‘--stage’ flag. It is also important to note that the default behavior of the script is to overwrite any existing results directory.

Lastly, the pipeline implements three levels of parallel processing. The first is that jobs can be spread across nodes in a **slurm** cluster. If this level of parallelization is desired, the user should invoke the ‘-c’ or ‘--cluster’ flag, with no arguments. Additionally, both the **trim_3prime** and **bwa** stages of the pipeline can be executed in parallel, the user may specify the number of threads to use by invoking the ‘-t’ or ‘--threads’ flag followed by an integer argument indicating the number of desired threads. Additionally, the ‘-m’ or ‘--map’ flag can be used to specify the number of threads to be used in each bwa instance for read mapping. Therefore maximum number of threads per node, will be the argument of ‘--threads’ multiplied by the argument of ‘--map’.

## The ddradseq.log file

The **ddradseq.py** control script will write a detailed log file ("ddradseq.log") to the specified output directory after any significant level of activity is invoked by the script. This log file contains informational messages, warnings, or errors that can be used for troubleshooting subsequent runs of the pipeline. Users are encouraged the read it, regardless of the outcome of any given run.

## The parse_pool program

The **parse_pool** program is written in the C language and successful compilation of the program depends on the **libfasta** library being installed (see **Installation** section below). The **parse_pool** program operates on a single "raw" fastQ input file at a time. The program will read the standard Illumina multiplexing index and create a separate output fastQ file for each index sequence. All files will be written to a new directory called “pool/”, which will be located in a parent directory that is specified by the user with the “-o” or “--out” flags,

* * *


./parse_pool -h

Usage : parse_pool [OPTIONS] [FASTQ] [CSV]

Parse fastQ file into separate files by Illumina index

Available options

  -o  DIR    Parent directory to write output files [default: same as input fastQ]

  -h         Display this help message

* * *


In addition to the input fastQ file, you may also notice that the **parse_pool** program also requires a "csv" file. This is a comma-separated file with two columns: the sequence of the Illumina index and the associated multiplexing pool identifier string. It is this pool identifier (“POOLID”) that will be used to name the resulting output fastQ files, e.g., “pool_POOLID.R1.fq.gz”. These output files will be appended, and not overwritten, when separate “raw” fastQ files are processed.

Pipeline end-users will not directly call any of the worker programs themselves, instead all worker programs are called from the python control script. In this way, the user must specify the highest level directory that contains all of the "raw" input fastQ files that are to be processed with **parse_pool**.

## The trim_barcode program

Like **parse_pool**, the **trim_barcode** worker program is also written in C.

* * *


./trim_barcode -h

Usage : trim_barcode [OPTIONS] [FASTQ.R1] [FASTQ.R2] [CSV]

Trim custom barcode sequences from fastQ file

Available options

  -d  INT    Edit distance [default: 1]

  -o  DIR    Parent directory to write output files [default: same as input fastQ]

  -h         Display this help message

* * *


## The trim_3prime program

* * *


./trim_3prime -h

Usage : trim_3prime [OPTIONS] [FASTQ.R1] [FASTQ.R2]

Aligns mate pairs in fastQ files and trims 3' end of reverse sequence

Available options

  -a  INT    Log odds score for identical bases [default: 1]

  -b  INT    Log odds score for different bases [default: 3]

  -o  DIR    Parent directory to write output files [default: same as input fastQ]

  -q  INT    Gap open cost [default: 5]

  -r  INT    Gap extend cost [default: 2]

  -t  INT    Minimum alignment score [default: 100]

  -h         Display this help message

* * *


## Installation

The ddRadSeq pipeline worker programs depend on the **libfasta** C library. Either .rpm packages for redhat-based systems or .deb packages for Debian-based systems can be downloaded from [https://lummei.net](https://lummei.net). Once the **libfasta** library is properly installed, the worker programs and python control script can be downloaded as a tar bundle from [https://lummei.net](https://lummei.net). 

