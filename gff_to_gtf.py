#!/usr/bin/env python 
"""
Program to convert data from GFF to GTF 

Usage: python gff_to_gtf.py in.gff > out.gtf 
"""

import re
import sys

import GFFParser

def printGTF(tinfo):
    """
    writing result file in GTF format

    @args tinfo: parsed object from gff file
    @type tinfo: numpy array 
    """

    for ent1 in tinfo:
        for idx, tid in enumerate(ent1['transcripts']):
            #print 
            #print tid 
            #print ent1['name']
            #print ent1['source']
            #print ent1['gene_info']
            #print ent1['chr']
            #print 
            
            exons = ent1['exons'][idx]
            cds_exons = ent1['cds_exons'][idx]

            stop_codon = start_codon = ()

            if ent1['strand'] == '+':
                start_codon = (cds_exons[0][0], cds_exons[0][0]+2) 
                stop_codon = (cds_exons[-1][1]-2, cds_exons[-1][1]) 
            elif ent1['strand'] == '-':
                start_codon = (cds_exons[-1][1]-2, cds_exons[-1][1])
                stop_codon = (cds_exons[0][0], cds_exons[0][0]+2)
            else:
                print 'STRAND information %s, skip the transcript' % ent1['strand']
                pass 
                
            print start_codon, stop_codon

            break

    
if __name__ == "__main__": 

    try:
        gff_fname = sys.argv[1]
    except:
        print __doc__
        sys.exit(-1)

    Transcriptdb = GFFParser.Parse(gff_fname)  

    printGTF(Transcriptdb) 
