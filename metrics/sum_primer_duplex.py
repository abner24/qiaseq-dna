from collections import defaultdict
import sys

import pysam


def run(cfg):
    '''
    '''
    print('sum_primer_duplex starting...')
    
    # input
    readset       = cfg.readSet
    inbam         = readset + '.umi_merge.bam'
    umi_mark_file = readset + '.umi_mark.for.sum.primer.txt'

    # output
    outfile       = readset + '.sum_primer_duplex.txt'
    # output metric file header
    out_header = "|".join(
        ["read set", "primer", "strand", "chrom", "loc5", "loc3", "Total unique UMIs", \
         "Total reads","Total CC-UMIs","Total reads for CC-UMIs","Total TT-UMIs", \
         "Total reads for TT-UMIs", "Total NN-UMIs", "Total reads for NN-UMIs", \
         "Only NN-UMIs", "Duplex-UMIs"])    

    # store duplex tag for each UMI
    IN = pysam.AlignmentFile(inbam,"rb")
    duplex_by_umi = defaultdict(lambda:defaultdict(int))    
    for read in IN:
        umi = read.get_tag(cfg.tagNameUmi)
        duplex_id = read.get_tag(cfg.tagNameDuplex) 
        duplex_by_umi[umi][duplex_id] += 1
        
    IN.close()

    primer_metrics = defaultdict(lambda:defaultdict(int))
    primers = set()
    primer_info = {}

    # store duplex UMIs for each primer
    with open(umi_mark_file, "r") as IN:
        for line in IN:
            chrom, strand, umi_loc, umi, num_reads, num_alignments, mt_read_idx, \
                is_resample,frag_len,primer,primer_loc5 = line.strip().split("|")            
            primers.add(primer)
            loc3 = int(primer_loc5) + len(primer) - 1 if strand == "0" \
                   else int(primer_loc5) - len(primer) + 1
            primer_info[primer] = "|".join(
                [readset, primer,strand,chrom,primer_loc5,str(loc3)])
            molecule = "-".join([chrom, strand, umi_loc, umi])

            if duplex_by_umi[molecule]["CC"] > 0:
                primer_metrics[primer]["CC_duplex"] += 1

            elif duplex_by_umi[molecule]["TT"] > 0:
                primer_metrics[primer]["TT_duplex"] += 1

            elif duplex_by_umi[molecule]["NN"] > 0:
                primer_metrics[primer]["NN_duplex"] += 1
                

            if duplex_by_umi[molecule]["NN"] > 0 and \
               duplex_by_umi[molecule]["CC"] == 0 and \
               duplex_by_umi[molecule]["TT"] == 0: # UMIs with only NN tags
                primer_metrics[primer]["only_NN"] += 1            

            primer_metrics[primer]["CC"] += duplex_by_umi[molecule]["CC"]        
            primer_metrics[primer]["TT"] += duplex_by_umi[molecule]["TT"]
            primer_metrics[primer]["NN"] += duplex_by_umi[molecule]["NN"]

            # for debug
            primer_metrics[primer]["UMI"] += 1
            primer_metrics[primer]["reads"] += int(num_reads)

    # write to output file
    with open(outfile, 'w') as OUT:
        OUT.write(out_header)
        OUT.write("\n")
        for primer in primers:    
            CC_umi     = primer_metrics[primer]["CC_duplex"]
            TT_umi     = primer_metrics[primer]["TT_duplex"]
            NN_umi     = primer_metrics[primer]["NN_duplex"]
            only_NN_umi =  primer_metrics[primer]["only_NN"]
            CC_reads   = primer_metrics[primer]["CC"]
            TT_reads   = primer_metrics[primer]["TT"]
            NN_reads   = primer_metrics[primer]["NN"]    

            num_reads  = primer_metrics[primer]["reads"]
            unique_umi = primer_metrics[primer]["UMI"]
            duplex_umi = (CC_umi + TT_umi) - (unique_umi - only_NN_umi)

            out = [primer_info[primer], str(unique_umi), str(num_reads), str(CC_umi), \
                   str(CC_reads), str(TT_umi), str(TT_reads), str(NN_umi), str(NN_reads), \
                   str(only_NN_umi), str(duplex_umi)]
            OUT.write("|".join(out))
            OUT.write("\n")


if __name__ == '__main__':
    readset   = sys.argv[1]
    paramfile = sys.argv[2]
    import core.run_config
    cfg = core.run_config.run(readset, paramfile)
    run(cfg)
