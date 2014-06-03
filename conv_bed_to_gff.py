#!/usr/bin/env python

#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Written (W) 2010 Vipin T Sreedharan, Friedrich Miescher Laboratory of the Max Planck Society
# Copyright (C) 2010 Max Planck Society
#
# Description : Convert a BED format file to GFF3 format

import re, sys

def __main__():

    try:
        bed_fh = open(sys.argv[1], 'rU')
    except:
        sys.stderr.write('BED format file fail to open, Cannot continue...\n')
        sys.stderr.write('USAGE: bed_to_gff3_converter.py <bed file> > *.gff3\n')
        sys.exit(-1)
    print '##gff-version 3'
    for line in bed_fh: 
        line = line.strip( '\n\r' ).split( '\t' )
        if re.match('#', line[0]):continue
        if len(line) != 12: # considering BED lines with 12 fields
            line = '\t'.join(line)
            sys.stdout.write('Warning: Invalid BED line found- ' + line + '\n') 
            continue
        if len(line[-1].split(',')) != len(line[-2].split(',')):continue # checking the consistency b/w relative start of exon and its length
        rstart = line[-1].split(',')
        if rstart[-1] == '': rstart.pop()
        exon_len = line[-2].split(',')
        if exon_len[-1] == '': exon_len.pop()
        if len(rstart) != int(line[-3]): continue # checking the number of exons and block count are same
        if line[5] != '+' and line[5] != '-':line[5] = '.' # replace the unknown starnd with '.' 
        # write feature lines to the result file 
        print line[0] + '\tbed2gff\tgene\t' + str(int(line[1]) + 1) + '\t' + line[2] + '\t' + line[4] + '\t' + line[5] + '\t.\t' + 'ID=Gene:' + line[3] + ';Name=Gene:' + line[3] 
        print line[0] + '\tbed2gff\ttranscript\t' + str(int(line[1]) + 1) + '\t' + line[2] + '\t' + line[4] + '\t' + line[5] + '\t.\t' + 'ID=' + line[3] + ';Name=' + line[3] + ';Parent=Gene:' + line[3]
        st = int(line[1])
        for ex_cnt in range(int(line[-3])):
            start = st + int(rstart[ex_cnt]) + 1
            stop = start + int(exon_len[ex_cnt]) - 1
            print line[0] + '\tbed2gff\texon\t' + str(start) + '\t' + str(stop) + '\t' + line[4] + '\t' + line[5] + '\t.\t' + 'Parent=' + line[3]
    bed_fh.close()

if __name__ == "__main__": __main__()
