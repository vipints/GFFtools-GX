#!/usr/bin/env python
"""
Convert data from Genbank format to GFF. 

Usage: 
python gbk_to_gff_conv.py in.gbk > out.gff 

Requirements:
    BioPython:- http://biopython.org/
"""

import os, sys, re
import collections
from Bio import SeqIO
from helper import open_file

def feature_table(chr_id, source, orient, genes, transcripts, cds, exons, unk):
    """
    Write the feature information
    """
    for gname, ginfo in genes.items():
        line = [str(chr_id), 
                'gbk_to_gff',
                ginfo[3],
                str(ginfo[0]),
                str(ginfo[1]),
                '.',
                ginfo[2],
                '.',
                'ID=%s;Name=%s;Note=%s' % (str(gname), str(gname), ginfo[-1])]
        print '\t'.join(line) 
        ## construct the transcript line is not defined in the original file 
        t_line = [str(chr_id), 'gbk_to_gff', source, 0, 1, '.', ginfo[2], '.'] 

        if not transcripts:
            t_line.append('ID=Transcript:%s;Parent=%s' % (str(gname), str(gname)))

            if exons: ## get the entire transcript region  from the defined feature
                t_line[3] = str(exons[gname][0][0])
                t_line[4] = str(exons[gname][0][-1])
            elif cds:
                t_line[3] = str(cds[gname][0][0])
                t_line[4] = str(cds[gname][0][-1])
            print '\t'.join(t_line) 

            if exons:
                exon_line_print(t_line, exons[gname], 'Transcript:'+str(gname), 'exon')

            if cds:
                exon_line_print(t_line, cds[gname], 'Transcript:'+str(gname), 'CDS')
                if not exons:
                    exon_line_print(t_line, cds[gname], 'Transcript:'+str(gname), 'exon')

        else: ## transcript is defined 
            for idx in transcripts[gname]: 
                t_line[2] = idx[3]
                t_line[3] = str(idx[0])
                t_line[4] = str(idx[1])
                t_line.append('ID='+str(idx[2])+';Parent='+str(gname))
                print '\t'.join(t_line) 
                
                ## feature line print call 
                if exons:
                    exon_line_print(t_line, exons[gname], str(idx[2]), 'exon')
                if cds:
                    exon_line_print(t_line, cds[gname], str(idx[2]), 'CDS')
                    if not exons:
                        exon_line_print(t_line, cds[gname], str(idx[2]), 'exon')

    if len(genes) == 0: ## feature entry with fragment information 
        
        line = [str(chr_id), 'gbk_to_gff', source, 0, 1, '.', orient, '.'] 
        fStart = fStop = None 

        for eid, ex in cds.items(): 
            fStart = ex[0][0] 
            fStop = ex[0][-1]

        for eid, ex in exons.items(): 
            fStart = ex[0][0] 
            fStop = ex[0][-1]

        if fStart or fStart:

            line[2] = 'gene'
            line[3] = str(fStart)
            line[4] = str(fStop)
            line.append('ID=Unknown_Gene_' + str(unk) + ';Name=Unknown_Gene_' + str(unk))
            print "\t".join(line)

            if not cds:
                line[2] = 'transcript'
            else:
                line[2] = 'mRNA'
            line[8] = 'ID=Unknown_Transcript_' + str(unk) + ';Parent=Unknown_Gene_' + str(unk)
            print "\t".join(line)
           
            if exons:
                exon_line_print(line, cds[None], 'Unknown_Transcript_' + str(unk), 'exon')
                
            if cds:
                exon_line_print(line, cds[None], 'Unknown_Transcript_' + str(unk), 'CDS')
                if not exons:
                    exon_line_print(line, cds[None], 'Unknown_Transcript_' + str(unk), 'exon')
                
            unk +=1 

    return unk

def exon_line_print(temp_line, trx_exons, parent, ftype):
    """
    Print the EXON feature line 
    """
    for ex in trx_exons:
        temp_line[2] = ftype
        temp_line[3] = str(ex[0])
        temp_line[4] = str(ex[1])
        temp_line[8] = 'Parent=%s' % parent
        print '\t'.join(temp_line)

def gbk_parse(fname):
    """
    Extract genome annotation recods from genbank format 
    """
    fhand = open_file(gbkfname)
    unk = 1 

    for record in SeqIO.parse(fhand, "genbank"):

        gene_tags = dict()
        tx_tags = collections.defaultdict(list) 
        exon = collections.defaultdict(list) 
        cds = collections.defaultdict(list) 
        mol_type, chr_id = None, None 

        for rec in record.features:

            if rec.type == 'source':
                mol_type = rec.qualifiers['mol_type'][0]
                try:
                    chr_id = rec.qualifiers['chromosome'][0]
                except:
                    chr_id = record.name 
                continue 

            strand='-'
            strand='+' if rec.strand>0 else strand
            
            fid = None 
            try:
                fid = rec.qualifiers['gene'][0]
            except:
                pass

            transcript_id = None
            try:
                transcript_id = rec.qualifiers['transcript_id'][0]
            except:
                pass 

            if re.search(r'gene', rec.type):
                gene_tags[fid] = (rec.location._start.position+1, 
                                    rec.location._end.position, 
                                    strand,
                                    rec.type
                                    )
            elif rec.type == 'exon':
                exon[fid].append((rec.location._start.position+1, 
                                    rec.location._end.position))
            elif rec.type=='CDS':
                cds[fid].append((rec.location._start.position+1, 
                                    rec.location._end.position))
            else: 
                # get all transcripts 
                if transcript_id: 
                    tx_tags[fid].append((rec.location._start.position+1,
                                    rec.location._end.position, 
                                    transcript_id,
                                    rec.type))
        # record extracted, generate feature table
        unk = feature_table(chr_id, mol_type, strand, gene_tags, tx_tags, cds, exon, unk)
        
    fhand.close()

if __name__=='__main__': 

    try:
        gbkfname = sys.argv[1]
    except:
        print __doc__
        sys.exit(-1)

    ## extract gbk records  
    gbk_parse(gbkfname) 
