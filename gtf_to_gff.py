#!/usr/bin/env python
"""
Convert Gene Transfer Format [GTF] to Generic Feature Format Version 3 [GFF3].

Usage: python gtf_to_gff.py in.gtf > out.gff3  
    
Requirement:
    GFFParser.py: https://github.com/vipints/GFFtools-GX/blob/master/GFFParser.py    
    helper.py : https://github.com/vipints/GFFtools-GX/blob/master/helper.py
    
Copyright (C) 
    2009-2012 Friedrich Miescher Laboratory of the Max Planck Society, Tubingen, Germany.
    2012-2014 Memorial Sloan Kettering Cancer Center New York City, USA.
"""

import re
import sys
import GFFParser
import helper

def GFFWriter(gtf_content):
    """
    write the feature information to GFF format

    @args gtf_content: Parsed object from gtf file 
    @type gtf_content: numpy array
    """

    print '##gff-version 3'

    for ent1 in gtf_content:

        chr_name = ent1['chr']
        strand = ent1['strand']
        start = ent1['start']
        stop = ent1['stop']
        source = ent1['source']
        ID = ent1['name']
        Name = ent1['gene_info']['Name']

        print '%s\t%s\tgene\t%d\t%d\t.\t%s\t.\tID=%s;Name=%s' % (chr_name, source, start, stop, strand, ID, Name) 

        for idx, tid in enumerate(ent1['transcripts']):

            t_start = ent1['exons'][idx][0][0]
            t_stop = ent1['exons'][idx][-1][-1]
            t_type = ent1['transcript_type'][idx]

            if ent1['exons'][idx].any() and ent1['cds_exons'][idx].any():
                utr5_exons, utr3_exons = helper.buildUTR(ent1['cds_exons'][idx], ent1['exons'][idx], strand)

            print '%s\t%s\t%s\t%d\t%d\t.\t%s\t.\tID=%s;Parent=%s' % (chr_name, source, t_type, t_start, t_stop, strand, tid[0], ID) 

            for ex_cod in utr5_exons:
                print '%s\t%s\tfive_prime_UTR\t%d\t%d\t.\t%s\t.\tParent=%s' % (chr_name, source, ex_cod[0], ex_cod[1], strand, tid[0]) 

            for ex_cod in ent1['cds_exons'][idx]:
                print '%s\t%s\tCDS\t%d\t%d\t.\t%s\t%d\tParent=%s' % (chr_name, source, ex_cod[0], ex_cod[1], strand, ex_cod[2], tid[0]) 

            for ex_cod in utr3_exons:
                print '%s\t%s\tthree_prime_UTR\t%d\t%d\t.\t%s\t.\tParent=%s' % (chr_name, source, ex_cod[0], ex_cod[1], strand, tid[0]) 

            for ex_cod in ent1['exons'][idx]:
                print '%s\t%s\texon\t%d\t%d\t.\t%s\t.\tParent=%s' % (chr_name, source, ex_cod[0], ex_cod[1], strand, tid[0]) 
            

def __main__():

    try:
        gtf_fname = sys.argv[1]
    except:
        print __doc__
        sys.exit(-1)

    gtf_file_content = GFFParser.Parse(gtf_fname)  

    GFFWriter(gtf_file_content)

if __name__ == "__main__": 
    __main__()
