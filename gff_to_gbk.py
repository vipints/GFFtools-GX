#!/usr/bin/env python 
"""
Convert data from GFF and associated genome sequence in fasta file into GenBank.

Usage: 
python gff_to_gbk.py in.gff in.fasta out.gbk 

Requirements:
    BioPython:- http://biopython.org/
    helper.py : https://github.com/vipints/GFFtools-GX/blob/master/helper.py

Copyright (C) 
    2010-2012 Friedrich Miescher Laboratory of the Max Planck Society, Tubingen, Germany.
    2012-2014 Memorial Sloan Kettering Cancer Center New York City, USA.
"""

import sys
import helper

from Bio import SeqIO
from Bio.Alphabet import generic_dna

def __main__():
    """
    """

    try:
        gff_fname = sys.argv[1]
        fasta_fname = sys.argv[2]
        gb_fname = sys.argv[3]
    except: 
        print __doc__
        sys.exit(-1)

    fasta_fh = helper.open_file(fasta_fname) 

    fasta_rec = SeqIO.to_dict(SeqIO.parse(fasta_fh, "fasta", generic_dna))
    fasta_fh.close()

    #
    #from Core import GFF
    #gff_rec = GFF.parse(gff_fname, fasta_rec)
    #gb_fh=open(gb_fname, "w")
    #SeqIO.write(gff_rec, gb_fh, "genbank")
    #gb_fh.close()
    #

if __name__=="__main__":
    __main__()
