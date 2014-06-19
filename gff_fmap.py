#!/usr/bin/env python
"""
GFF feature mapping program, to find the relation between different features described in a given GFF file. 

Usage: 
python gff_fmap.py in.gff > out.txt 

Courtesy: Brad Chapman 
    Few functions are inherited from bcbio-GFFParser program. 
"""

import re
import sys 
import urllib
import collections
from helper import open_file

def _gff_line_map(line):
    """Parses a line of GFF into a dictionary.
    Given an input line from a GFF file, this:
        - breaks it into component elements
        - determines the type of attribute (flat, parent, child or annotation)
        - generates a dictionary of GFF info 
    """
    gff3_kw_pat = re.compile("\w+=")
    def _split_keyvals(keyval_str):
        """Split key-value pairs in a GFF2, GTF and GFF3 compatible way.
        GFF3 has key value pairs like:
          count=9;gene=amx-2;sequence=SAGE:aacggagccg
        GFF2 and GTF have:           
          Sequence "Y74C9A" ; Note "Clone Y74C9A; Genbank AC024206"
          name "fgenesh1_pg.C_chr_1000003"; transcriptId 869
        """
        quals = collections.defaultdict(list)
        if keyval_str is None:
            return quals
        # ensembl GTF has a stray semi-colon at the end
        if keyval_str[-1] == ';':
            keyval_str = keyval_str[:-1]
        # GFF2/GTF has a semi-colon with at least one space after it.
        # It can have spaces on both sides; wormbase does this.
        # GFF3 works with no spaces.
        # Split at the first one we can recognize as working
        parts = keyval_str.split(" ; ")
        if len(parts) == 1:
            parts = keyval_str.split("; ")
            if len(parts) == 1:
                parts = keyval_str.split(";")
        # check if we have GFF3 style key-vals (with =)
        is_gff2 = True
        if gff3_kw_pat.match(parts[0]):
            is_gff2 = False
            key_vals = [p.split('=') for p in parts]
        # otherwise, we are separated by a space with a key as the first item
        else:
            pieces = []
            for p in parts:
                # fix misplaced semi-colons in keys in some GFF2 files
                if p and p[0] == ';':
                    p = p[1:]
                pieces.append(p.strip().split(" "))
            key_vals = [(p[0], " ".join(p[1:])) for p in pieces]
        for key, val in key_vals:
            # remove quotes in GFF2 files
            if (len(val) > 0 and val[0] == '"' and val[-1] == '"'):
                val = val[1:-1] 
            if val:
                quals[key].extend(val.split(','))
            # if we don't have a value, make this a key=True/False style
            # attribute
            else:
                quals[key].append('true')
        for key, vals in quals.items():
            quals[key] = [urllib.unquote(v) for v in vals]
        return quals, is_gff2

    def _nest_gff2_features(gff_parts):
        """Provide nesting of GFF2 transcript parts with transcript IDs.

        exons and coding sequences are mapped to a parent with a transcript_id
        in GFF2. This is implemented differently at different genome centers
        and this function attempts to resolve that and map things to the GFF3
        way of doing them.
        """
        # map protein or transcript ids to a parent
        for transcript_id in ["transcript_id", "transcriptId", "proteinId"]:
            try:
                gff_parts["quals"]["Parent"] = \
                        gff_parts["quals"][transcript_id]
                break
            except KeyError:
                pass
        # case for WormBase GFF -- everything labelled as Transcript or CDS
        for flat_name in ["Transcript", "CDS"]:
            if gff_parts["quals"].has_key(flat_name):
                # parent types
                if gff_parts["type"] in [flat_name]:
                    if not gff_parts["id"]:
                        gff_parts["id"] = gff_parts["quals"][flat_name][0]
                        gff_parts["quals"]["ID"] = [gff_parts["id"]]
                # children types
                elif gff_parts["type"] in ["intron", "exon", "three_prime_UTR",
                        "coding_exon", "five_prime_UTR", "CDS", "stop_codon",
                        "start_codon"]:
                    gff_parts["quals"]["Parent"] = gff_parts["quals"][flat_name]
                break

        return gff_parts

    line = line.strip()
    if line == '':return [('directive', line)] # sometimes the blank lines will be there 
    if line[0] == '>':return [('directive', '')] # sometimes it will be a FATSA header
    if line[0] == "#":
        return [('directive', line[2:])]
    elif line:
        parts = line.split('\t')
        if len(parts) == 1 and re.search(r'\w+', parts[0]):return [('directive', '')] ## GFF files with FASTA sequence together 
        assert len(parts) == 9, line
        gff_parts = [(None if p == '.' else p) for p in parts]
        gff_info = dict()
            
        # collect all of the base qualifiers for this item
        quals, is_gff2 = _split_keyvals(gff_parts[8])

        gff_info["is_gff2"] = is_gff2

        if gff_parts[1]:quals["source"].append(gff_parts[1])
        gff_info['quals'] = dict(quals)

        # if we are describing a location, then we are a feature
        if gff_parts[3] and gff_parts[4]:
            gff_info['type'] = gff_parts[2]
            gff_info['id'] = quals.get('ID', [''])[0]
            
            if is_gff2:gff_info = _nest_gff2_features(gff_info)
            # features that have parents need to link so we can pick up
            # the relationship
            if gff_info['quals'].has_key('Parent'):
                final_key = 'child'
            elif gff_info['id']:
                final_key = 'parent'
            # Handle flat features
            else:
                final_key = 'feature'
        # otherwise, associate these annotations with the full record
        else:
            final_key = 'annotation'
        return [(final_key, gff_info)]
    
def parent_child_id_map(gff_handle):
    """
    Provide a mapping of parent to child relationships in the file.
    Gives a dictionary of parent child relationships:

    keys -- tuple of (source, type) for each parent
    values -- tuple of (source, type) as children of that parent
    """
    # collect all of the parent and child types mapped to IDs
    parent_sts = dict()
    child_sts = collections.defaultdict(list)
    for line in gff_handle:
        line_type, line_info = _gff_line_map(line)[0]
        if (line_type == 'parent' or (line_type == 'child' and line_info['id'])):
            parent_sts[line_info['id']] = (line_info['quals']['source'][0], line_info['type'])
        if line_type == 'child':
            for parent_id in line_info['quals']['Parent']:
                child_sts[parent_id].append((line_info['quals']['source'][0], line_info['type']))
    gff_handle.close()
    # generate a dictionary of the unique final type relationships
    pc_map = collections.defaultdict(list)
    for parent_id, parent_type in parent_sts.items():
        for child_type in child_sts[parent_id]:
            pc_map[parent_type].append(child_type)
    pc_final_map = dict()
    for ptype, ctypes in pc_map.items():
        unique_ctypes = list(set(ctypes))
        unique_ctypes.sort()
        pc_final_map[ptype] = unique_ctypes
    # some cases the GFF file represents a single feature type 
    if not pc_final_map:
        for fid, stypes in parent_sts.items():
            pc_final_map[stypes] = dict()
    # generate a report on feature id mapping in the file 
    print '+---------------------+---------------------------------+'
    print '| Parent feature type | Associated child feature type(s)|'
    print '+---------------------+---------------------------------+'
    for key, value in pc_final_map.items():
        print key[0], key[1]
        for child_to in value:
            print '\t\t\t|-',child_to[0], child_to[1]
        print '+---------------------+---------------------------------+'


if __name__=='__main__':

    try:
        gff_file = sys.argv[1]
    except:
        print __doc__
        sys.exit(-1)
    
    gff_handle = open_file(gff_file)
    parent_child_id_map(gff_handle)
