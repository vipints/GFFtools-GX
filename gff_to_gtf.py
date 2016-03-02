#!/usr/bin/env python 
"""
Program to convert data from GFF to GTF 

Usage: python gff_to_gtf.py in.gff > out.gtf 

Requirement:
    GFFParser.py: https://github.com/vipints/GFFtools-GX/blob/master/GFFParser.py    

Copyright (C) 
    2009-2012 Friedrich Miescher Laboratory of the Max Planck Society, Tubingen, Germany.
    2012-2015 Memorial Sloan Kettering Cancer Center New York City, USA.
"""

import re
import sys
import numpy
import GFFParser

def printGTF(tinfo):
    """
    writing result file in GTF format

    @args tinfo: parsed object from gff file
    @type tinfo: numpy array 
    """

    for ent1 in tinfo:
        for idx, tid in enumerate(ent1['transcripts']):
            
            exons = ent1['exons'][idx]
            if numpy.isnan(numpy.min(exons)):
                continue 
            cds_exons = ent1['cds_exons'][idx]

            stop_codon = start_codon = ()

            if ent1['strand'] == '+':
                if cds_exons.any():
                    start_codon = (cds_exons[0][0], cds_exons[0][0]+2) 
                    stop_codon = (cds_exons[-1][1]-2, cds_exons[-1][1]) 
            elif ent1['strand'] == '-':
                if cds_exons.any():
                    start_codon = (cds_exons[-1][1]-2, cds_exons[-1][1])
                    stop_codon = (cds_exons[0][0], cds_exons[0][0]+2)
            else:
                sys.stdout.write('STRAND information missing - %s, skip the transcript - %s\n' % (ent1['strand'], tid[0]))
                pass 
                
            last_cds_cod = 0 
            for idz, ex_cod in enumerate(exons):

                sys.stdout.write('%s\t%s\texon\t%d\t%d\t.\t%s\t.\tgene_id "%s"; transcript_id "%s"; exon_number "%d"; gene_name "%s"; \n' % (ent1['chr'], ent1['source'], ex_cod[0], ex_cod[1], ent1['strand'], ent1['name'], tid[0], idz+1, ent1['gene_info']['Name']))

                if cds_exons.any():
                    try:
                        sys.stdout.write('%s\t%s\tCDS\t%d\t%d\t.\t%s\t%d\tgene_id "%s"; transcript_id "%s"; exon_number "%d"; gene_name "%s"; \n' % (ent1['chr'], ent1['source'], cds_exons[idz][0], cds_exons[idz][1], ent1['strand'], cds_exons[idz][2], ent1['name'], tid[0], idz+1, ent1['gene_info']['Name']))
                        last_cds_cod = idz 
                    except:
                        pass 

                    if idz == 0:
                        sys.stdout.write('%s\t%s\tstart_codon\t%d\t%d\t.\t%s\t%d\tgene_id "%s"; transcript_id "%s"; exon_number "%d"; gene_name "%s"; \n' % (ent1['chr'], ent1['source'], start_codon[0], start_codon[1], ent1['strand'], cds_exons[idz][2], ent1['name'], tid[0], idz+1, ent1['gene_info']['Name']))

            if stop_codon:
                sys.stdout.write('%s\t%s\tstop_codon\t%d\t%d\t.\t%s\t%d\tgene_id "%s"; transcript_id "%s"; exon_number "%d"; gene_name "%s"; \n' % (ent1['chr'], ent1['source'], stop_codon[0], stop_codon[1], ent1['strand'], cds_exons[last_cds_cod][2], ent1['name'], tid[0], idz+1, ent1['gene_info']['Name']))

    
if __name__ == "__main__": 

    try:
        gff_fname = sys.argv[1]
    except:
        print __doc__
        sys.exit(-1)

    Transcriptdb = GFFParser.Parse(gff_fname)  

    printGTF(Transcriptdb) 
