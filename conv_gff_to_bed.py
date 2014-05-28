#!/usr/bin/env python
"""
Convert genome annotation data in GFF/GTF to a 12 column BED format. 
BED format typically represents the transcript models. 

Usage: python gff_to_bed_conv.py in.gff > out.bed  
"""

import re
import sys
import GFFParser

def writeBED(tinfo):
    """
    writing result files in bed format 
    """
    for ent1 in tinfo:
        for idx, tid in enumerate(ent1['transcripts']):
            exon_cnt = len(ent1['exons'][idx])
            exon_len = ''
            exon_cod = '' 
            rel_start = None 
            rel_stop = None 
            for idz, ex_cod in enumerate(ent1['exons'][idx]):#check for exons of corresponding transcript  
                exon_len += str(int(ex_cod[1])-int(ex_cod[0])+1) + ','
                if idz == 0: #calculate the relative start position 
                    exon_cod += '0,'
                    rel_start = int(ex_cod[0])
                    rel_stop = int(ex_cod[1])
                else:
                    exon_cod += str(int(ex_cod[0])-rel_start) + ','
                    rel_stop = int(ex_cod[1])
            
            if exon_len:
                out_print = [ent1['chr'],
                            str(rel_start),
                            str(rel_stop),
                            tid[0],
                            ent1['score'][0], 
                            ent1['strand'], 
                            str(rel_start),
                            str(rel_stop),
                            '0',
                            str(exon_cnt),
                            exon_len,
                            exon_cod]
                print '\t'.join(out_print)  
    
def __main__():
    try:
        query_file = sys.argv[1]
    except:
        print __doc__
        sys.exit(-1)

    Transcriptdb = GFFParser.Parse(query_file)  
    writeBED(Transcriptdb)

if __name__ == "__main__": 
    __main__() 
