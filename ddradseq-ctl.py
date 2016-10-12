#!/usr/bin/env python2.7
#----------------------------------------------------------
# File: ddradseq.py
#
# Author: Lummei Analytics LLC
# Last updated: August 2016
#
# Description: A control script for running the ddRadSeq
# pipeline in the Hoekstra lab
#----------------------------------------------------------
import os
import sys
import time
import glob
import csv
import re
import logging
import argparse
import textwrap
import subprocess
import multiprocessing
# Globally scoped initialization of logging class
logger = logging.getLogger("ddradseq-ctl")
# Globally scoped arrays of stages
stageParsePool = [0, 4]
stageTrimBarcode = [0, 1, 5]
stageTrimThreePrime = [0, 1, 2, 6]
stageBWA = [0, 1, 2, 3, 7]
def main(args):
    # Parse command line arguments
    parameters = getCommandLine(args)
    # Get user name
    userName = os.environ['USER']
    # Start pipeline logger
    initializeLog(userName)
    # Check for write permissions on output directory
    checkPermissions(parameters.outputDir)
    # Check for executable worker programs in user's path
    #checkWorkers(parameters.stage)
    # Check CSV database files
    if parameters.stage in stageParsePool:
        checkDatabase(parameters.poolCSV)
    if parameters.stage in stageTrimBarcode:
        checkDatabase(parameters.barcodeCSV)
    # Check for files in existing output directories, if any
    checkExistingDirs(parameters.outputDir, parameters.stage)
    # Run the parse_pool stage of the pipeline
    if parameters.stage in stageParsePool:
        runParsePool(
            parameters.inputDir,
            parameters.outputDir,
            parameters.poolCSV,
			parameters.numThreads)
    # Run the trim_barcode stage of the pipeline
    if parameters.stage in stageTrimBarcode:
        runTrimBarcode(parameters.outputDir, parameters.barcodeCSV, parameters.numThreads)
    # Run the trim_3prime stage of the pipeline
    if parameters.stage in stageTrimThreePrime:
        runTrimThreePrime(parameters.outputDir, parameters.numThreads)
    # Run the bwa stage of the pipeline
    if parameters.stage in stageBWA:
        runBWA(
            parameters.outputDir,
            parameters.referenceFile,
            parameters.numThreads,
            parameters.threadsBWA)
    # Finish pipeline
    logger.info("ddradseq.py pipeline done")
"""
------------------------------------------------------------
execBWA()
------------------------------------------------------------
Executes the bwa application in parallel
Takes the directory of the output directory and the input
fastA reference genome file as arguments
Return value is trivial
Exits the script if bwa encounters run-time error
(e.g., non-zero return value to shell)
------------------------------------------------------------
"""
def execBWA(threadID, nthreadsBWA, outputDir, referenceFile,
            fileStart, fileEnd, filenamesForwardSort, filenamesReverseSort):
    for i in range(fileStart, fileEnd):
        fileParse = re.split('[_.]', filenamesForwardSort[i])
        outfileBWA = "bam/individ_" + fileParse[1] + ".bam"
        fulloutfileBWA = os.path.join(outputDir, outfileBWA)
        cmdBWA = "bwa mem -t {:d} {} {} {} | samtools view -bS -T {} -o {} -".format(nthreadsBWA, referenceFile,
                                                                                     filenamesForwardSort[i], filenamesReverseSort[i], referenceFile, fulloutfileBWA)
        logger.info(
            "Thread {:d}: Running command: {}".format(threadID, cmdBWA))
        p = subprocess.Popen(
            cmdBWA,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
        # Catch any run time errors
        exitCode = p.returncode
        if exitCode:
            logger.error(
                "thread {:d}: bwa command encountered run-time error".format(threadID))
            logger.error(stderr)
        else:
            logger.info(
                "thread {:d}: bwa command executed successfully".format(threadID))
    return 0
"""
------------------------------------------------------------
runBWA()
------------------------------------------------------------
Function to run bwa to map reads
Takes the directory of the output directory and the input
fastA reference genome file as arguments
Return value is trivial
Exits the script if bwa encounters run-time error
(e.g., non-zero return value to shell)
------------------------------------------------------------
"""
def runBWA(outputDir, referenceFile, numThreads, nthreadsBWA):
    logger.info("Running the bwa stage of the pipeline")
    logger.info(
        "Starting to map reads to reference genome in {}.".format(referenceFile))
    # Get the stage start time
    startTime = time.time()
    # Construct sorted arrays of mate pair input file names
    fullPathForward = os.path.join(outputDir, 'final/final_*.R1.fq.gz')
    filenamesForward = glob.glob(fullPathForward)
    if not filenamesForward:
        logger.error("bwa error: no input R1 fastQ files found")
        sys.exit('FATAL ERROR: no R1 fastQ files for bwa input found')
    fullPathReverse = os.path.join(outputDir, 'final/final_*.R2.fq.gz')
    filenamesReverse = glob.glob(fullPathReverse)
    if not filenamesReverse:
        logger.error("bwa error: no input R2 fastQ files found")
        sys.exit('FATAL ERROR: no R2 fastQ files for bwa input found')
    filenamesForwardSort = sorted(filenamesForward)
    filenamesReverseSort = sorted(filenamesReverse)
    # Make a bwa subdirectory in outdir
    outSubDir = os.path.join(outputDir, "bam")
    try:
        os.makedirs(outSubDir)
    except OSError:
        pass
    # Iterate through input fastQ files
    numFiles = len(filenamesForwardSort)
    fileStart = [0] * numThreads
    fileEnd = [0] * numThreads
    jobs = []
    filesPerThread = (numFiles + numThreads - 1) / numThreads
    for tid in range(numThreads):
        fileStart[tid] = tid * filesPerThread
        fileEnd[tid] = (tid + 1) * filesPerThread
    fileEnd[numThreads - 1] = numFiles
    for tid in range(numThreads):
        proc = multiprocessing.Process(
            target=execBWA,
            args=(tid,
                  nthreadsBWA,
                  outputDir,
                  referenceFile,
                  fileStart[tid],
                  fileEnd[tid],
                  filenamesForwardSort,
                  filenamesReverseSort))
        jobs.append(proc)
    # Run jobs
    for pr in jobs:
        pr.start()
    # Wait for join to finish
    for pr in jobs:
        pr.join()
    # Print information about stage run time
    elapsedTime = time.time() - startTime
    min, sec = divmod(int(elapsedTime), 60)
    hour, min = divmod(min, 60)
    logger.info(
        "Total elapsed bwa time: {0:02d}:{1:02d}:{2:02d}".format(hour, min, sec))
    return 0
"""
------------------------------------------------------------
execTrimThreePrime()
------------------------------------------------------------
Function to execute the trim_3prime program in parallel
Takes the directory of the output fastQ files from
trim_barcode work program as an argument
Return value is trivial
Exits the script if trim_3prime encounters run-time error
(e.g., non-zero return value to shell)
------------------------------------------------------------
"""
def execThreePrime(threadID, outputDir, fileStart,
                   fileEnd, filenamesForwardSort, filenamesReverseSort):
    for i in range(fileStart, fileEnd):
        cmdTrimThreePrime = "trim_3prime -o {} {} {}".format(
            outputDir,
            filenamesForwardSort[i],
            filenamesReverseSort[i])
        logger.info(
            "Thread {:}: Running command: {}".format(
                threadID,
                cmdTrimThreePrime))
        proc = subprocess.Popen(
            cmdTrimThreePrime,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        # Catch any run time errors
        exitCode = proc.returncode
        if exitCode:
            logger.error(
                "Thread {:d}: trim_3prime command encountered run-time error".format(threadID))
            logger.error(stderr)
        else:
            logger.info(
                "Thread {:d}: trim_3prime command executed successfully".format(threadID))
    return 0
"""
------------------------------------------------------------
runTrimThreePrime()
------------------------------------------------------------
Function to run the trim_3prime worker program
Takes the directory of the output fastQ files from
trim_barcode work program as an argument
Return value is trivial
Exits the script if trim_3prime encounters run-time error
(e.g., non-zero return value to shell)
------------------------------------------------------------
"""
def runTrimThreePrime(outputDir, numThreads):
    logger.info("Running the trim_3prime stage of the pipeline")
    # Get the stage start time
    startTime = time.time()
    # Construct sorted arrays of mate pair input file names
    fullPathForward = os.path.join(outputDir, 'trim/trim_*.R1.fq.gz')
    filenamesForward = glob.glob(fullPathForward)
    if not filenamesForward:
        logger.error("trim_3prime error: no input R1 fastQ files found")
        sys.exit('FATAL ERROR: no R1 fastQ files for trim_3prime input found')
    fullPathReverse = os.path.join(outputDir, 'trim/trim_*.R2.fq.gz')
    filenamesReverse = glob.glob(fullPathReverse)
    if not filenamesReverse:
        logger.error("trim_3prime error: no input R2 fastQ files found")
        sys.exit('FATAL ERROR: no R2 fastQ files for trim_3prime input found')
    filenamesForwardSort = sorted(filenamesForward)
    filenamesReverseSort = sorted(filenamesReverse)
    numFiles = len(filenamesForwardSort)
    fileStart = [0] * numThreads
    fileEnd = [0] * numThreads
    jobs = []
    filesPerThread = (numFiles + numThreads - 1) / numThreads
    for tid in range(numThreads):
        fileStart[tid] = tid * filesPerThread
        fileEnd[tid] = (tid + 1) * filesPerThread
    fileEnd[numThreads - 1] = numFiles
    for tid in range(numThreads):
        proc = multiprocessing.Process(
            target=execThreePrime,
            args=(tid,
                  outputDir,
                  fileStart[tid],
                  fileEnd[tid],
                  filenamesForwardSort,
                  filenamesReverseSort))
        jobs.append(proc)
    # Run the jobs
    for pr in jobs:
        pr.start()
    # Wait for all jobs to finish
    for pr in jobs:
        pr.join()
    # Print information about stage run time
    elapsedTime = time.time() - startTime
    min, sec = divmod(int(elapsedTime), 60)
    hour, min = divmod(min, 60)
    logger.info(
        "Total elapsed trim_3prime time: {0:02d}:{1:02d}:{2:02d}".format(hour, min, sec))
    return 0
"""
------------------------------------------------------------
runTrimBarcode()
------------------------------------------------------------
Function to run the trim_barcode worker program
Takes the name of the directory containing the output fastQ
files from the parse_pool worker program and the
comma-separated text file with custom barcodes and individual
sample names as arguments
Return value is trivial
Exits the script if trim_barcode encounters run-time error
(e.g., non-zero return value to shell)
------------------------------------------------------------
"""
def runTrimBarcode(outputDir, barcodeCSV, numThreads):
    logger.info("Running the trim_barcode stage of the pipeline")
    # Get the stage start time
    startTime = time.time()
    # Construct sorted arrays of mate pair input file names
    fullPathForward = os.path.join(outputDir, 'pool/pool_*.R1.fq.gz')
    filenamesForward = glob.glob(fullPathForward)
    if not filenamesForward:
        logger.error("trim_barcode error: no input R1 fastQ files found")
        sys.exit('FATAL ERROR: no R1 fastQ files for trim_barcode input found')
    fullPathReverse = os.path.join(outputDir, 'pool/pool_*.R2.fq.gz')
    filenamesReverse = glob.glob(fullPathReverse)
    if not filenamesReverse:
        logger.error("trim_barcode error: no input R2 fastQ files found")
        sys.exit('FATAL ERROR: no R2 fastQ files for trim_barcode input found')
    filenamesForwardSort = sorted(filenamesForward)
    filenamesReverseSort = sorted(filenamesReverse)
    # Iterate through input fastQ files
    for index, f in enumerate(filenamesForwardSort):
        cmdTrimBarcode = "trim_barcode -t {} -o {} {} {} {}".format(
            numThreads, outputDir, f, filenamesReverseSort[index], barcodeCSV)
        logger.info("Running command: {}".format(cmdTrimBarcode))
        proc = subprocess.Popen(
            cmdTrimBarcode,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        # Catch any run time errors
        exitCode = proc.returncode
        if exitCode:
            logger.error("trim_barcode command encountered run-time error")
            logger.error(stderr)
        else:
            logger.info("trim_barcode command executed successfully")
    # Print information about stage run time
    elapsedTime = time.time() - startTime
    min, sec = divmod(int(elapsedTime), 60)
    hour, min = divmod(min, 60)
    logger.info(
        "Elapsed trim_barcode time: {0:02d}:{1:02d}:{2:02d}".format(hour, min, sec))
    return 0
"""
------------------------------------------------------------
runParsePool()
------------------------------------------------------------
Function to run the parse_pool worker program
Takes the list of infile names and the comma-separated text
file with Illumina indices and pool names as arguments
Return value is trivial
Exits the script if parse_pool encounters run-time error
(e.g., non-zero return value to shell)
------------------------------------------------------------
"""
def runParsePool(inputDir, outputDir, poolCSV, numThreads):
    # Get a list of infile fastQ file names
    infileNames = getFastqFilenames(inputDir, "*.fastq.gz")
    # Count the number of input fastQ files
    numInfiles = len(infileNames)
    logger.info("Running the parse_pool stage of the pipeline")
    logger.info("Read {} input fastQ files".format(numInfiles))
    # Get the stage start time
    startTime = time.time()
    # Iterate through input fastQ files
    for f in infileNames:
        sizeInfile = os.path.getsize(f) / (1024 * 1024.0)
        logger.info("File: {}  size: {:.2f} Mb".format(f, sizeInfile))
        cmdParsePool = "parse_pool -t {:d} -o {} {} {}".format(numThreads, outputDir, f, poolCSV)
        logger.info("Running command: {}".format(cmdParsePool))
        proc = subprocess.Popen(
            cmdParsePool,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        # Catch any run time errors
        exitCode = proc.returncode
        if exitCode:
            logger.error("parse_pool command encountered run-time error")
            logger.error(stderr)
        else:
            logger.info("parse_pool command executed successfully")
    # Print information about stage run time
    elapsedTime = time.time() - startTime
    min, sec = divmod(int(elapsedTime), 60)
    hour, min = divmod(min, 60)
    logger.info(
        "Elapsed parse_pool time: {0:02d}:{1:02d}:{2:02d}".format(hour, min, sec))
    return 0
"""
------------------------------------------------------------
checkWorkers()
------------------------------------------------------------
Function to check for executable worker programs in the
user's path
Exits the script if any of the three worker programs are
not found
Return value is trivial
------------------------------------------------------------
"""
def checkWorkers(stage):
    # If needed, check for the parse_pool program
    if stage in stageParsePool:
        try:
            logger.info(
                "Checking whether parse_pool executable is in user PATH")
            pathParsePool = subprocess.check_output(
                'which parse_pool', shell=True)
        except:
            logger.error("parse_pool executable not found in user PATH")
            sys.exit('FATAL ERROR: parse_pool executable not found')
        logger.info("parse_pool executable found at: {}".format(pathParsePool))
    # If needed, check for the trim_barcode program
    if stage in stageTrimBarcode:
        try:
            logger.info(
                "Check whether the trim_barcode executable is in user PATH")
            pathTrimBarcode = subprocess.check_output(
                'which trim_barcode', shell=True)
        except:
            logger.error("trim_barcode executable not found in user PATH")
            sys.exit('FATAL ERROR: trim_barcode executable not found')
        logger.info(
            "trim_barcode executable found at: {}".format(pathTrimBarcode))
    # If needed, check for the trim_3prime program
    if stage in stageTrimThreePrime:
        try:
            logger.info(
                "Checking whether the trim_3prime executable is in user PATH")
            pathTrimThreePrime = subprocess.check_output(
                'which trim_3prime', shell=True)
        except:
            logger.error("trim_3prime executable not found in user PATH")
            sys.exit('FATAL ERROR: trim_3prime executable not found')
        logger.info(
            "trim_3prime executable found at: {}".format(pathTrimThreePrime))
    # If needed, check for both bwa and samtools programs
    if stage in stageBWA:
        try:
            logger.info("Checking whether the bwa executable is in user PATH")
            pathBWA = subprocess.check_output('which bwa', shell=True)
        except:
            logger.error("bwa executable not found in user PATH")
            sys.exit('FATAL ERROR: bwa executable not found')
        logger.info("bwa executable found at: {}".format(pathBWA))
        try:
            logger.info(
                "Checking whether the samtools executable is in user PATH")
            pathSAMtools = subprocess.check_output(
                'which samtools', shell=True)
        except:
            logger.error("samtools executable not found in user PATH")
            sys.exit('FATAL ERROR: samtools executable not found')
        logger.info("samtools executable found at: {}".format(pathSAMtools))
    return 0
"""
------------------------------------------------------------
checkDatabase()
------------------------------------------------------------
Function to check for a CSV database for duplicate keys
Exits the script if any duplicate keys are found in the DB
Return value is trivial
------------------------------------------------------------
"""
def checkDatabase(CSVfile):
    logger.info(
        "Checking CSV database file {} for duplicate keys".format(CSVfile))
    with open(CSVfile) as inputCSV:
        reader = csv.DictReader(inputCSV, fieldnames=('KEY', 'VALUE'))
        d = {}
        numDuplicates = 0
        for row in reader:
            if row['KEY'] in d:
                numDuplicates += 1
            else:
                d[row['KEY']] = row['VALUE']
        if numDuplicates > 0:
            logger.error(
                "CSV database file {} contains {:d} duplicate entries".format(CSVfile, numDuplicates))
            sys.exit('FATAL ERROR: CSV database contains duplicate entries')
        else:
            logger.info(
                "CSV database file {} has no duplicate keys".format(CSVfile))
    return 0
"""
------------------------------------------------------------
checkPermissions()
------------------------------------------------------------
Function to check for output directory write permissions
Return value is trivial
------------------------------------------------------------
"""
def checkPermissions(outputDir):
    isDir = os.path.isdir(outputDir)
    if isDir:
        logger.info(
            "Confirmed that specified path to output directory {} exists".format(outputDir))
        logger.info(
            "Now testing for user write permissions on {}".format(outputDir))
        try:
            fileName = os.path.join(outputDir, "test")
            f = open(fileName, "w")
            f.close()
            os.remove(fileName)
        except Exception as e:
            logger.error("{}".format(e))
            sys.exit('FATAL ERROR: cannot write to specified output directory')
        logger.info(
            "Confirmed that the user is able to write to {}".format(outputDir))
    return 0
"""
------------------------------------------------------------
checkExistingDirs()
------------------------------------------------------------
Function to check for output directory write permissions
Return value is trivial
------------------------------------------------------------
"""
def checkExistingDirs(outputDir, stage):
    logger.info(
        "Checking for existing output directories in {}".format(outputDir))
    # Log the results of the subdirectory scan
    if stage in stageParsePool:
        poolDir = os.path.join(outputDir, 'pool')
        pool_isdir = os.path.isdir(poolDir)
        if pool_isdir:
            logger.info(
                "parse_pool output directory {} already exists on disk".format(poolDir))
            logger.warning("Files will be overwritten")
            # Remove all existing files in the pool subdirectory
            for root, dirs, files in os.walk(poolDir):
                for name in files:
                    os.remove(os.path.join(root, name))
            # Remove the existing pool subdirectory
            os.rmdir(poolDir)
        else:
            logger.info(
                "parse_pool output directory {} does not exist-- it will be created".format(poolDir))
    if stage in stageTrimBarcode:
        trimDir = os.path.join(outputDir, 'trim')
        trim_isdir = os.path.isdir(trimDir)
        if trim_isdir:
            logger.info(
                "trim_barcode output directory {} already exists on disk".format(trimDir))
            logger.warning("Files will be overwritten")
            # Remove all existing files in the trim subdirectory
            for root, dirs, files in os.walk(trimDir):
                for name in files:
                    os.remove(os.path.join(root, name))
            # Remove existing trim subdirectory
            os.rmdir(trimDir)
        else:
            logger.info(
                "trim_barcode output directory {} does not exist-- it will be created".format(trimDir))
    if stage in stageTrimThreePrime:
        finalDir = os.path.join(outputDir, 'final')
        final_isdir = os.path.isdir(finalDir)
        if final_isdir:
            logger.info(
                "trim_3prime output directory {} already exists on disk".format(finalDir))
            logger.warning("Files will be overwritten")
            # Remove all existing files in the final subdirectory
            for root, dirs, files in os.walk(finalDir):
                for name in files:
                    os.remove(os.path.join(root, name))
            # Remove existing final subdirectory
            os.rmdir(finalDir)
        else:
            logger.info(
                "trim_3prime output directory {} does not exist-- it will be created".format(finalDir))
    if stage in stageBWA:
        bamDir = os.path.join(outputDir, 'bam')
        bam_isdir = os.path.isdir(bamDir)
        if bam_isdir:
            logger.info(
                "bwa output directory {} already exists on disk".format(bamDir))
            logger.warning("Files will be overwritten")
            # Remove all existing files in the bam subdirectory
            for root, dirs, files in os.walk(bamDir):
                for name in files:
                    os.remove(os.path.join(root, name))
            # Remove existing bam subdirectory
            os.rmdir(bamDir)
        else:
            logger.info(
                "bwa output directory {} does not exist-- it will be created".format(bamDir))
    return 0
"""
------------------------------------------------------------
getFastqFilenames()
------------------------------------------------------------
Function to read in names of all files in a specified
directory, or its child folders, that end in the suffix .fastq.gz
Takes the name of the parent directory as argument and the
file name extension pattern to search for
Returns a list of input fastQ filenames
------------------------------------------------------------
"""
def getFastqFilenames(inputDir, pattern):
    fileList = list()
    for root, directs, filenames in os.walk(inputDir):
        searchPattern = os.path.join(root, pattern)
        fileList.extend(glob.glob(searchPattern))
        for direct in directs:
            dirPath = os.path.join(root, direct)
            searchPattern = os.path.join(dirPath, pattern)
            fileList.extend(glob.glob(searchPattern))
    if not fileList:
        logger.error("No input fastQ files found in {}".format(inputDir))
        sys.exit('FATAL ERROR: no input fastQ files found')
    return fileList
"""
------------------------------------------------------------
initializeLog()
------------------------------------------------------------
Function to setup the log file
Takes the login name of the user running the script
Return value is trivial
------------------------------------------------------------
"""
def initializeLog(userName):
    logger.setLevel(logging.INFO)
    fh = logging.FileHandler("ddradseq.log")
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - ' +
        '%(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    logger.info("ddradseq.py program started by user {}".format(userName))
    return 0
"""
------------------------------------------------------------
getCommandLine()
------------------------------------------------------------
Function to read command line arguments using the argparse
module
Returns a Namespace object with command line parameters
------------------------------------------------------------
"""
def getCommandLine(args):
    parser = argparse.ArgumentParser(
        __file__,
            description="Control script for ddRadSeq pipeline",
            usage='use "python %(prog)s --help" for more information',
            epilog='For more detailed information on how to run this script, see the overview document.',
            formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-s', '--stage',
                        default=0,
                        type=int,
                        required=False,
                        dest='stage',
                        help='''\
Start pipeline at a particular stage
0: Start pipeline at parse_pool and run to end [default]
1: Start pipeline at trim_barcode and run to end
2: Start pipeline at trim_3prime and run to end
3: Start pipeline at bwa and run to end
4: Run only parse_pool stage and exit
5: Run only trim_barcode stage and exit
6: Run only trim_3prime state and exit
7: Run only bwa stage and exit''',
                        metavar='N'
                        )
    parser.add_argument('-d', '--dir',
                        required=True,
                        dest='inputDir',
                        help="Directory containing input fastQ files",
                        metavar='DIR'
                        )
    parser.add_argument('-o', '--out',
                        required=True,
                        dest='outputDir',
                        help="Parent directory for output files",
                        metavar='DIR'
                        )
    parser.add_argument('-r', '--ref',
                        required=False,
                        dest='referenceFile',
                        help="Name of file with reference sequence for read mapping",
                        metavar='FILE'
                        )
    parser.add_argument('-p', '--pool',
                        required=False,
                        dest='poolCSV',
                        help="CSV file with Illumina multiplex indices and pool names",
                        metavar='FILE'
                        )
    parser.add_argument('-b', '--barcode',
                        required=False,
                        dest='barcodeCSV',
                        help="CSV file with custom barcodes and individual identifiers",
                        metavar='FILE'
                        )
    parser.add_argument('-c', '--cluster',
                        default=False,
                        required=False,
                        dest='cluster',
                        action='store_true',
                        help="If this switch is present, script runs in cluster (slurm) mode"
                        )
    parser.add_argument('-m', '--map',
                        default=1,
                        type=int,
                        required=False,
                        dest='threadsBWA',
                        help="Number of threads available for bwa read mapping",
                        metavar='N'
                        )
    parser.add_argument('-t', '--threads',
                        default=1,
                        type=int,
                        required=False,
                        dest='numThreads',
                        help="Number of threads available for concurrency",
                        metavar='N'
                        )
    parser.add_argument('-v', '--version',
                        action='version',
                        version='%(prog)s 1.0'
                        )
    parameters = parser.parse_args(args)
    # Check if there are enough available CPU threads
    availableCores = multiprocessing.cpu_count()
    requestedCores = parameters.numThreads * parameters.threadsBWA
    if requestedCores > availableCores:
        parser.error(
            'FATAL ERROR: {:d} threads requested, only {:d} threads available on host'.format(
                requestedCores, availableCores))
    # Check if pool database CSV file is needed
    if parameters.stage in stageParsePool:
        if not parameters.poolCSV:
            parser.error('FATAL ERROR: pool CSV file needed as input')
        else:
            if not os.path.isfile(parameters.poolCSV):
                parser.error(
                    'FATAL ERROR: multiplexing/pool CSV file {} not found'.format(parameters.poolCSV))
    # Check if barcode database CSV file is needed
    if parameters.stage in stageTrimBarcode:
        if not parameters.barcodeCSV:
            parser.error('FATAL ERROR: barcode CSV file needed as input')
        else:
            if not os.path.isfile(parameters.barcodeCSV):
                parser.error(
                    'FATAL ERROR: custom barcode CSV file {} not found'.format(
                        parameters.barcodeCSV))
    # Check if reference genome file is needed
    if parameters.stage in stageBWA:
        if not parameters.referenceFile:
            parser.error(
                'FATAL ERROR: reference genome fastA file needed as input')
        else:
            if not os.path.isfile(parameters.referenceFile):
                parser.error(
                    'FATAL ERROR: reference fastA file {} not found'.format(
                        parameters.referenceFile))
    return parameters
if __name__ == '__main__':
    main(sys.argv[1:])
