#!/usr/bin/env python
"""
Convert Gene Transfer Format [GTF] to Generic Feature Format Version 3 [GFF3].

Usage: python gtf_to_gff.py in.gtf > out.gff3  
    
Requirement:
    
Copyright (C) 
    2009-2012 Friedrich Miescher Laboratory of the Max Planck Society, Tubingen, Germany.
    2012-2014 Memorial Sloan Kettering Cancer Center New York City, USA.
"""

import re
import sys
import GFFParser

def GFFWriter(gtf_file_cont):
    """
    Write feature details to GFF3 
    """

    print '##gff-version 3'

    for contig, contig_info in sorted(gtf_file_cont.items()): # chromosome 
        for feature, details in contig_info.items(): # gene with source 
            
            gene_start= gene_stop = []
            gnames = None
            tnames= transcript_details = dict()

            for ftid, tinfo in details.items(): # transcripts 
                tinfo['exon'].sort() # coordinate system in ascending order. 
                tinfo['CDS'].sort()
                if tinfo['exon']:
                    gene_start.append(tinfo['exon'][0][0])
                    gene_stop.append(tinfo['exon'][-1][1])
                if not gene_start:
                    continue
                orient = tinfo['info'][0]
                tnames[ftid]=tinfo['info'][-1]
                gnames=tinfo['info'][-2]
                if len(tinfo['CDS']) == 0: # non coding transcript 
                    transcript_details[ftid] = dict(info = tinfo['info'], 
                                                exon = tinfo['exon'], 
                                                tpe = 'transcript')
                else:
                    if tinfo['sp_cod']: # stop codon are seperated from CDS, add the coordinates based on strand 
                        if orient == '+':
                            if tinfo['sp_cod'][0][0]-tinfo['CDS'][-1][1] == 1:
                                tinfo['CDS'][-1] = (tinfo['CDS'][-1][0], tinfo['sp_cod'][0][1])
                            else:
                                tinfo['CDS'].append(tinfo['sp_cod'][0])
                        if orient == '-':
                            if tinfo['CDS'][0][0]-tinfo['sp_cod'][0][1] == 1:
                                tinfo['CDS'][0] = (tinfo['sp_cod'][0][0], tinfo['CDS'][0][1])
                            else:
                                tinfo['CDS'].insert(0, tinfo['sp_cod'][0])
                    if tinfo['exon']:
                        utr5, utr3 = buildUTR(tinfo['CDS'], tinfo['exon'], orient) # getting UTR info from CDS and exon.
                        transcript_details[ftid] = dict(info = tinfo['info'], 
                                                    exon = tinfo['exon'], 
                                                    utr5 = utr5, 
                                                    utr3 = utr3, 
                                                    cds = tinfo['CDS'], 
                                                    tpe = 'mRNA')
            if gene_start and gene_stop: # displying Gene, transcript and subfeatures
                gene_start.sort()
                gene_stop.sort()
                if gnames == None:
                    gnames = feature[0] # assign gene name as gene id, if not defined 
                pline = [str(contig),
                        feature[1],
                        'gene',
                        str(gene_start[0]),
                        str(gene_stop[-1]),
                        '.',
                        orient,
                        '.',
                        'ID=' + feature[0] + ';Name=' + gnames]
                print '\t'.join(pline)

                for dtid, dinfo in transcript_details.items():
                    if dinfo['info'][3]:
                        pline = [str(contig),
                                feature[1],
                                dinfo['tpe'],
                                str(dinfo['exon'][0][0]),
                                str(dinfo['exon'][-1][1]),
                                dinfo['info'][1],
                                orient,
                                '.',
                                'ID=' + str(dtid) + ';Parent=' + feature[0] + ';Name=' + str(dinfo['info'][3]) ]
                    else:
                        pline = [str(contig),
                                feature[1],
                                dinfo['tpe'],
                                str(dinfo['exon'][0][0]),
                                str(dinfo['exon'][-1][1]),
                                dinfo['info'][1],
                                orient,
                                '.',
                                'ID=' + dtid + ';Parent=' + feature[0]]
                    print '\t'.join(pline) 

                    if 'utr5' in dinfo:
                        for ele in dinfo['utr5']:
                            pline = [str(contig),
                                    feature[1], 
                                    'five_prime_UTR',
                                    str(ele[0]), 
                                    str(ele[1]),
                                    '.',
                                    orient,
                                    '.',
                                    'Parent=' + dtid]
                            print '\t'.join(pline) 

                    if 'cds' in dinfo:
                        cds_w_phase = addCDSphase(orient, dinfo['cds'])
                        for ele in cds_w_phase:
                            pline = [str(contig),
                                    feature[1],
                                    'CDS',
                                    str(ele[0]),
                                    str(ele[1]),
                                    '.',
                                    orient,
                                    str(ele[-1]),
                                    'Parent=' + dtid]
                            print '\t'.join(pline) 

                    if 'utr3' in dinfo:
                        for ele in dinfo['utr3']:
                            pline = [str(contig),
                                    feature[1],
                                    'three_prime_UTR',
                                    str(ele[0]),
                                    str(ele[1]),
                                    '.',
                                    orient,
                                    '.',
                                    'Parent=' + dtid]
                            print '\t'.join(pline)

                    if 'exon' in dinfo:
                        intron_start = 0
                        for xq, ele in enumerate(dinfo['exon']):
                            
                            if xq > 0:
                                pline = [str(contig),
                                    feature[1],
                                    'intron',
                                    str(intron_start),
                                    str(ele[0]-1),
                                    '.',
                                    orient,
                                    '.',
                                    'Parent=' + dtid]
                                print '\t'.join(pline)

                            pline = [str(contig),
                                    feature[1],
                                    'exon',
                                    str(ele[0]),
                                    str(ele[1]),
                                    '.',
                                    orient,
                                    '.',
                                    'Parent=' + dtid]
                            print '\t'.join(pline) 

                            intron_start = ele[1]+1 


def getGTFcontent(gtf_file):
    """
    Extract GTF features 
    """

    GFH = _open_file(gtf_file)
    gtf_content, recall = dict(), None

    for rec in GFH:
        rec = rec.strip('\n\r')

        #skip empty line fasta identifier and commented line
        if not rec or rec[0] in  ['#', '>']:
            continue
        #skip the genome sequence 
        if not re.search('\t', rec):
            continue

        parts = rec.split('\t')
        assert len(parts) >= 8, rec 
    
        if re.search(r'^(start_codon|start-codon|startcodon)$', parts[2], re.IGNORECASE):
            continue

        gid= tid= gname= tname= ttype = None

        for attb in parts[-1].split(';'):
            if re.search(r'^\s?$', attb):
                continue

            attb = re.sub('"', '', attb).strip()
            attb = attb.split()

            if re.search(r'^(gene_id|geneid|name)$', attb[0], re.IGNORECASE): 
                gid = attb[1]
            elif re.search(r'^(transcript_id|transcriptId)$', attb[0], re.IGNORECASE):
                tid = attb[1]
            elif re.search(r'^(gene_name|genename)$', attb[0], re.IGNORECASE):
                gname = attb[1]
            elif re.search(r'^(transcript_name|transcriptname)$', attb[0], re.IGNORECASE):
                tname = attb[1]
            elif re.search(r'^(transcript_type)$', attb[0], re.IGNORECASE):
                ttype = attb[1]

        if gid == tid: #UCSC GTF files, gene & transcript have same identifier 
            gid = 'Gene:'+str(gid) 
            tid = 'Transcript:'+str(tid)

        if tid == None: #JGI GTF file dont have transcript ID for CDS line
            tid = recall 

        exon= cds= sp_cod= st_cod = []

        if re.search(r'^exon$', parts[2], re.IGNORECASE): 
            exon = [(int(parts[3]), int(parts[4]))]
        elif re.search(r'^CDS$', parts[2], re.IGNORECASE):
            cds = [(int(parts[3]), int(parts[4]))]
        elif re.search(r'^(stop_codon|stop-codon|stopcodon)$', parts[2], re.IGNORECASE):
            sp_cod = [(int(parts[3]), int(parts[4]))]
        else: #other lines are not required to GFF line 
            continue

        #creating feature connections 
        if parts[0] in gtf_content: # adding to existing chromosome
            if (gid, parts[1]) in gtf_content[parts[0]].keys(): # adding to existing gene 
                if tid in gtf_content[parts[0]][(gid, parts[1])].keys(): # adding to existing transcript
                    if exon:
                        gtf_content[parts[0]][(gid, parts[1])][tid]['exon'].append(exon[0])
                    elif cds:
                        gtf_content[parts[0]][(gid, parts[1])][tid]['CDS'].append(cds[0])
                    elif sp_cod:    
                        gtf_content[parts[0]][(gid, parts[1])][tid]['sp_cod'].append(sp_cod[0])
                else: # inserting new transcript
                    gtf_content[parts[0]][(gid, parts[1])][tid] = dict(exon = exon, 
                                                            CDS = cds, 
                                                            sp_cod = sp_cod, 
                                                            info = [parts[6], parts[5], gname, tname, ttype])
            else: # inserting new gene 
                gtf_content[parts[0]][(gid, parts[1])] = {tid : dict(exon = exon, 
                                                    CDS = cds,
                                                    sp_cod = sp_cod, 
                                                    info = [parts[6], parts[5], gname, tname, ttype])}
        else: # inserting new chromosome identifier 
            gtf_content[parts[0]] = {(gid, parts[1]) : {tid : dict(exon = exon, 
                                            CDS = cds,
                                            sp_cod = sp_cod, 
                                            info = [parts[6], parts[5], gname, tname, ttype])}}
        recall = tid #set previous id for CDS line 

    GFH.close()
    return gtf_content

def __main__():

    try:
        gtf_fname = sys.argv[1]
    except:
        print __doc__
        sys.exit(-1)

    gtf_file_content = getGTFcontent(gtf_fname)

    GFFWriter(gtf_file_content)

if __name__ == "__main__": 
    __main__()
