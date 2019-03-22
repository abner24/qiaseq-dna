import glob
import gzip
import os
import os.path
import math
import subprocess
from multiprocessing.dummy import Pool as ThreadPool

# our modules
import primer_trim

#-------------------------------------------------------------------------------------
def runShellCommand(cmd):
    # run shell command, capture stdout and stderror (assumes log is not large and not redirected to disk)
    try:
        log = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
        error = False
    except subprocess.CalledProcessError, ex:
        log = ex.output
        error = True
        
    # print command log output
    for line in log.split("\n"):
        print("prep_trim: " + line.strip())
     
    # re-raise exception now that error detail has been printed
    if error:
        raise(ex)

#-------------------------------------------------------------------------------------
def worker(cmd):
    try:
        subprocess.check_call(cmd, shell=True)
        return True
    except:
        return False
  
#-------------------------------------------------------------------------------------
def splitReadFile(readFile,filePrefixOut,readSide,numBatchesMax,deleteLocalFiles):

    # validate existence of input read file
    if not os.path.isfile(readFile):
        raise UserWarning("input read file does not exist: " + readFile)
           
    # get number of reads in the fastq file
    isGzFile = readFile.endswith(".fastq.gz")
    if isGzFile:
        line = subprocess.check_output("zcat {} | wc -l".format(readFile), shell=True)
        numLines = int(line)
    else:
        line = subprocess.check_output("wc -l " + readFile, shell=True)
        vals = line.split(" ")
        numLines = int(vals[0])
    if numLines % 4 != 0 or numLines == 0:
        raise UserWarning("truncated or empty fastq read file {}".format(readFile))
    numReads = numLines / 4
    
    # set batch size
    batchSize = 4 * int(math.ceil(1.00 * numReads / numBatchesMax))
    
    # open input fastq read file
    if isGzFile:
        fileIn = gzip.open(readFile,"rb")
    else:
        fileIn = open(readFile,"r")
  
    # write fastq batches to disk - one batch to be done on each CPU core
    numLinesOut = 0
    batchNum = 0
    fileOut = None
    for line in fileIn:
        # open new file if needed
        if numLinesOut == 0:
            fileOut = open("{}.{:04d}.{}.fastq".format(filePrefixOut,batchNum,readSide), "w")
            
        # write fastq line to disk
        fileOut.write(line)
        numLinesOut += 1
           
        # close output file if done with batch
        if numLinesOut == batchSize:
            fileOut.close()
            numLinesOut = 0
            batchNum += 1
   
    # close out last batch
    if numLinesOut != 0:
        fileOut.close()
        batchNum += 1
    numBatches = batchNum
    
    # delete local input file if no longer needed
    if deleteLocalFiles and len(os.path.dirname(readFile)) == 0:
        os.remove(readFile)
     
    # done
    return numReads, numBatches  

#-------------------------------------------------------------------------------------
def runReadTrimmer(cfg):
    ''' Use paired-end read trimmer for Illumina sequencing
    '''
    cmd = "PYTHONPATH= python3 {trimmer} " \
          "--r1 {R1} --r2 {R2} --primer-file {primer} " \
          "--out-metrics {summary} --out-r1 {outR1} --out-r2 {outR2} " \
          "--check-primer-side --primer3-bases-R1 {primer3R1} " \
          "--primer3-bases-R2 {primer3R2} --ncpu {ncpu} --primer-col 3 " \
          "--seqtype dna --umi-len 12 --common-seq-len 11 " \
          "--min-primer-side-len 25 --min-umi-side-len 1 --tagname-primer {pr} " \
          "--tagname-primer-error {pe} --tagname-umi {mi} "
    if cfg.duplex:
        cmd = cmd + "--tagname-duplex {du} --is-duplex --is-phased-adapters ".format(du = cfg.tagNameDuplex)
    if cfg.multimodal:
        cmd = cmd + "--is-multimodal "
    cmd = cmd + "> {log} 2>&1"
    
    cmd = cmd.format(
        trimmer = os.path.join(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__))),'read-trimmer/trimmer/run.py'),
        R1 = cfg.readFile1, R2 = cfg.readFile2, primer = cfg.primerFile,
        summary = cfg.readSet + '.prep.detail.summary.txt',
        outR1 = cfg.readSet + '.prep.R1.fastq',
        outR2 = cfg.readSet + '.prep.R2.fastq',
        primer3R1 = cfg.primer3Bases, primer3R2 = cfg.primer3Bases,
        ncpu = int(cfg.numCores),
        pr = cfg.tagNamePrimer, pe = cfg.tagNamePrimerErr, mi = cfg.tagNameUmiSeq,
        log = cfg.readSet + '.prep.log')
    subprocess.check_call(cmd, shell=True)

    # create a less dense summary file
    minimal_metrics = [
        lambda x: x == "read fragments total", lambda x: x.startswith("read fragments dropped,"),
        lambda x: x.startswith("read fragments with duplex tag"), lambda x: x.startswith("read fragments after trimming")]
    
    with open(cfg.readSet + '.prep.detail.summary.txt','r') as IN, \
         open(cfg.readSet + '.prep.summary.txt','w') as OUT:
        for line in IN:
            val, metricname = line.strip('\n').split('\t')
            for comparison in minimal_metrics:
                if comparison(metricname):
                    OUT.write(line)
    
#-------------------------------------------------------------------------------------
def run(cfg):
    # report start
    print("prep: starting read prep - trimming ends and UMI extraction")
    
    # get params
    readSet          = cfg.readSet
    readFile1        = cfg.readFile1
    readFile2        = cfg.readFile2
    trimScript       = cfg.trimScript
    cutadaptDir      = cfg.cutadaptDir
    deleteLocalFiles = cfg.deleteLocalFiles
    numCores         = int(cfg.numCores)
    tagNameUmiSeq    = cfg.tagNameUmiSeq
    tagNamePrimer    = cfg.tagNamePrimer
    tagNamePrimerErr = cfg.tagNamePrimerErr
    primerFile       = cfg.primerFile
    primer3Bases     = int(cfg.primer3Bases)
    platform = cfg.platform

    # Use new read-trimmer for Illumina reads.
    if cfg.platform.lower() == 'illumina':
        runReadTrimmer(cfg)
        return
    
    # The code below is now executed for IonTorrent Reads alone
    # It is kept as-is from before, so there are references to illumina 
    
    # debug check
    if cfg.readFile1 == cfg.readFile2:
        raise UserWarning("R1 and R2 have the same filename. Please fix the file paths for the input files.")

    # set sequencing type for trimming
    if cfg.platform.lower() == "illumina": 
        if cfg.duplex.lower() == "true": ## Duplex sequencing run
            seqtype = "illumina_duplex"
        else:
            seqtype = "illumina"
    else:
        seqtype = "iontorrent"

    # set output file prefix
    filePrefixOut = readSet + ".prep"
    
    # check for adapter on primer side - customer forgot to use the custom sequencing primer, or are a small part of a large HiSeq run
    if cfg.platform.lower() == "illumina":
        if readFile1.endswith(".fastq"):  
            cmd = "head -n 4000 " + readFile1 + " > " 
        else:
            cmd = "zcat " + readFile1 + " | head -n 4000 > " 
        cmd += filePrefixOut + ".temp1000.R1.fastq"
        subprocess.check_call(cmd, shell=True)
        cmd = cutadaptDir + "cutadapt -e 0.18 -O 18" \
            + " -g ^AATGTACAGTATTGCGTTTTG -n 1" \
            + " -o /dev/null " \
            +         filePrefixOut + ".temp1000.R1.fastq" \
            + " > " + filePrefixOut + ".cutadapt.5.R1.log 2>&1"
        subprocess.check_call(cmd, shell=True)
        os.remove(filePrefixOut + ".temp1000.R1.fastq")
        pctWithAdapter = None
        for line in open(filePrefixOut + ".cutadapt.5.R1.log", "r"):
            if line.startswith("Reads with adapters:"):
                idx1 = line.find("(")
                idx2 = line.find("%)")
                pctWithAdapter = float(line[idx1+1:idx2])
                break
        print("prep: check for wrong sequencing primer using first 1,000 reads: % wrong is {}".format(pctWithAdapter))
        if pctWithAdapter > 95.0:
            print("WARNING: R1 reads start with PCR adapter - custom sequencing primer was not used on primer side!")
            cmd = cutadaptDir + "cutadapt -e 0.18 -O 18" \
                + " -g ^AATGTACAGTATTGCGTTTTG -n 1" \
                + " -o " + filePrefixOut + ".fixed.R1.fastq " \
                         + readFile1  \
                + " > "  + filePrefixOut + ".cutadapt.5.R1.log 2>&1"
            subprocess.check_call(cmd, shell=True)
            readFile1 = filePrefixOut + ".fixed.R1.fastq"
     
    # split both read files into chunks to be processed in parallel
    numReads1, numBatches = splitReadFile(readFile1,filePrefixOut,"R1",numCores,deleteLocalFiles)
    if cfg.platform.lower() == "illumina": # ion-torrent reads are single end
        numReads2, numBatches = splitReadFile(readFile2,filePrefixOut,"R2",numCores,deleteLocalFiles)    
        # debug check
        if numReads2 != numReads1:
            raise UserWarning("prep: R1 and R2 read count not equal!")        
    # debug check
    if numReads1 == 0:
        raise UserWarning("prep: input read files are empty!")  
   
    # run cd-hit to cluster close primer sequences; creates the file {primerFile}.clusters
    primer_trim.cluster_primer_seqs(primerFile)
    # cache primer search datastructure to avoid repeated compute over batches; creates the file {primerFile}.kmer.cache
    #primer_trim.create_primer_search_datastruct(primerFile,primerFile+'.clusters',cache=True,cache_file=primerFile+".kmer.cache")
    
    # set up trimming work to be run in parallel sub-processes, using another python script
    workIn = []
    for batchNum in range(numBatches):
        filePrefixBatch = "{}.{:04d}".format(filePrefixOut,batchNum)
        cmd = "python {0} {1} {2} {3} {4} {5} {6} {7} {8} > {6}.log 2>&1 ".format(trimScript,cutadaptDir,tagNameUmiSeq,tagNamePrimer,tagNamePrimerErr,primer3Bases,filePrefixBatch,primerFile,seqtype)
        workIn.append(cmd)
        
    # run cutadapt and UMI extraction in parallel sub-processes
    print("prep: starting parallel trimming batches")
    pool = ThreadPool(min(numCores,len(workIn)))
    workOut = pool.map(worker, workIn)
    pool.close()
    pool.join()
    print("prep: completed parallel trimming batches")
    
    # make sure all batches of work completed successfully
    for batchNum in range(len(workOut)):
        if not workOut[batchNum]:
            raise Exception("read trimming failed for batch: {:04d}".format(batchNum))
   
    # concatenate the read files back into one file
    for readEnd in ("R1","R2"):
    
        # delete output file if it exists
        readFileOut = "{}.{}.fastq".format(filePrefixOut,readEnd)
        if os.path.isfile(readFileOut):
            os.remove(readFileOut)
      
        # concatenate read file and delete (Linux cat probaby faster than Python line reads)
        for batchNum in range(numBatches):
            readFileIn = "{}.{:04d}.{}.fastq".format(filePrefixOut, batchNum, readEnd)
            cmd = "cat {} >> {} ".format(readFileIn,readFileOut)
            subprocess.check_call(cmd,shell=True)
            os.remove(readFileIn)
   
    # concatenate the log files - note dangerous wildcards here!
    logFileOut = open(filePrefixOut + ".log","w")
    for logFileIn in glob.glob(filePrefixOut + ".0*.log"):
        IN = open(logFileIn,"r")
        logFileOut.write(IN.read())
        IN.close()
        os.remove(logFileIn)
    logFileOut.close()

    # For ion-torrent reads concatenate the umi and primer tag files        
    if seqtype == 'iontorrent':
        OUT1 = open(readSet + ".umi.tag.txt","w")
        OUT2 = open(readSet + ".primer.tag.txt","w")
        OUT3 = open(readSet + ".cutadapt.5.R1.txt","w")
        OUT4 = open(readSet + ".cutadapt.3.R1.txt","w")
        
        for umiTagFileTemp in sorted(glob.glob(filePrefixOut + "*.umi.tag.txt")):
            IN1 =  open(umiTagFileTemp,"r")
            OUT1.write(IN1.read()) # reading into memory
            IN1.close()
            os.remove(umiTagFileTemp)
        for primerTagFileTemp in sorted(glob.glob(filePrefixOut + "*.primer.tag.txt")):
            IN2 = open(primerTagFileTemp,"r")
            OUT2.write(IN2.read()) # reading into memory
            IN2.close()
            os.remove(primerTagFileTemp)
        for cutadapt5FileTemp in sorted(glob.glob(filePrefixOut + "*.cutadapt.5.R1.txt")):
            IN3 = open(cutadapt5FileTemp,"r")
            OUT3.write(IN3.read()) # reading into memory
            IN3.close()
            os.remove(cutadapt5FileTemp)
        for cutadapt3FileTemp in sorted(glob.glob(filePrefixOut + "*.cutadapt.3.R1.txt")):
            IN4 = open(cutadapt3FileTemp,"r")
            OUT4.write(IN4.read()) # reading into memory
            IN4.close()
            os.remove(cutadapt3FileTemp)
            
        OUT1.close()
        OUT2.close()
        OUT3.close()
        OUT4.close()
        
    # aggregate summary read count files - for some trim scripts these files contain important read count metrics
    output = []
    firstFile = True
    for sumFileIn in glob.glob(filePrefixOut + ".0*.summary.txt"):
        iLine = 0
        for line in open(sumFileIn,"r"):
            metricVal,metricName = line.strip().split("\t")
            metricVal = int(metricVal)
            if firstFile:
                output.append([metricVal,metricName])
            else:
                output[iLine][0] += metricVal
            iLine += 1
        firstFile = False
        os.remove(sumFileIn)
    sumFileOut = open(filePrefixOut + ".summary.txt","w")
    for row in output:
        sumFileOut.write("\t".join((str(x) for x in row)))
        sumFileOut.write("\n")
    sumFileOut.close()
    
    # report completion
    print("prep: done")
    
#-------------------------------------------------------------------------------------
if __name__ == "__main__":
    cfg = lambda:0
    cfg.platform = "illumina"
    cfg.readSet = "NEB_S2"
    cfg.readFile1 = "/mnt/webserver/datadisk/resources/jdicarlo/NEB_S2_L001_R1_001.fastq.gz"
    cfg.readFile2 = "/mnt/webserver/datadisk/resources/jdicarlo/NEB_S2_L001_R2_001.fastq.gz"
    cfg.cutadaptDir = "/srv/qgen/bin/cutadapt-1.10/"
    cfg.trimScript  = "prep_trim.py"
    cfg.numCores  = "32"
    cfg.deleteLocalFiles = True
    cfg.primer3Bases = -1
    run(cfg)
