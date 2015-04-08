#!/usr/bin/env python
"""
Convert genome annotation data in a 12 column BED format to GFF3. 

Usage: 
    python bed_to_gff.py in.bed > out.gff

Requirement:
    helper.py : https://github.com/vipints/GFFtools-GX/blob/master/helper.py

Copyright (C) 
    2009-2012 Friedrich Miescher Laboratory of the Max Planck Society, Tubingen, Germany.
    2012-2015 Memorial Sloan Kettering Cancer Center New York City, USA.
"""

import re
import sys
import helper 

def __main__():
    """
    main function 
    """

    try:
        bed_fname = sys.argv[1]
    except:
        print __doc__
        sys.exit(-1)

    bed_fh = helper.open_file(bed_fname)

    for line in bed_fh: 
        line = line.strip( '\n\r' )

        if not line or line[0] in  ['#']:
            continue 

        parts = line.split('\t') 
        assert len(parts) >= 12, line

        rstarts = parts[-1].split(',')
        rstarts.pop() if rstarts[-1] == '' else rstarts

        exon_lens = parts[-2].split(',')
        exon_lens.pop() if exon_lens[-1] == '' else exon_lens
        
        if len(rstarts) != len(exon_lens):
            continue # checking the consistency col 11 and col 12 

        if len(rstarts) != int(parts[-3]): 
            continue # checking the number of exons and block count are same
        
        if not parts[5] in ['+', '-']:
            parts[5] = '.' # replace the unknown strand with '.' 

        # bed2gff result line 
        print '%s\tbed2gff\tgene\t%d\t%s\t%s\t%s\t.\tID=Gene:%s;Name=Gene:%s' % (parts[0], int(parts[1])+1, parts[2], parts[4], parts[5], parts[3], parts[3])
        print '%s\tbed2gff\ttranscript\t%d\t%s\t%s\t%s\t.\tID=%s;Name=%s;Parent=Gene:%s' % (parts[0], int(parts[1])+1, parts[2], parts[4], parts[5], parts[3], parts[3], parts[3])

        st = int(parts[1])
        for ex_cnt in range(int(parts[-3])):
            start = st + int(rstarts[ex_cnt]) + 1
            stop = start + int(exon_lens[ex_cnt]) - 1
            print '%s\tbed2gff\texon\t%d\t%d\t%s\t%s\t.\tParent=%s' % (parts[0], start, stop, parts[4], parts[5], parts[3])

    bed_fh.close()


if __name__ == "__main__": 
    __main__()
