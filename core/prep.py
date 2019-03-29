import glob
import gzip
import os
import os.path
import math
import subprocess

# our modules
import misc.prep_ion

#-------------------------------------------------------------------------------------
def runReadTrimmer(cfg):
    ''' Use paired-end read trimmer for Illumina sequencing
    '''
    readTrimmerPath = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        'read-trimmer/trimmer/run.py')

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
        trimmer = readTrimmerPath,
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
    deleteLocalFiles = cfg.deleteLocalFiles
    numCores         = int(cfg.numCores)
    tagNameUmiSeq    = cfg.tagNameUmiSeq
    tagNamePrimer    = cfg.tagNamePrimer
    tagNamePrimerErr = cfg.tagNamePrimerErr
    primerFile       = cfg.primerFile
    primer3Bases     = int(cfg.primer3Bases)
    platform = cfg.platform
    
    # debug check
    if cfg.readFile1 == cfg.readFile2:
        raise UserWarning("R1 and R2 have the same filename. Please fix the file paths for the input files.")

    # Use new read-trimmer for Illumina reads.
    if cfg.platform.lower() == 'illumina':
        runReadTrimmer(cfg)   
    else: # legacy code for Ion Torrent reads
        misc.prep_ion.run(cfg)
    
    # report completion
    print("prep: done")
    
#-------------------------------------------------------------------------------------
if __name__ == "__main__":
    cfg = lambda:0
    cfg.platform = "illumina"
    cfg.readSet = "NEB_S2"
    cfg.readFile1 = "/mnt/webserver/datadisk/resources/jdicarlo/NEB_S2_L001_R1_001.fastq.gz"
    cfg.readFile2 = "/mnt/webserver/datadisk/resources/jdicarlo/NEB_S2_L001_R2_001.fastq.gz"
    cfg.numCores  = "32"
    cfg.deleteLocalFiles = True
    cfg.primer3Bases = -1
    run(cfg)
