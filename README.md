GFFtools-GX 
===========

A collection of tools for converting genome annotation between [GTF](https://genome.ucsc.edu/FAQ/FAQformat.html#format4), [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) and [GFF](https://genome.ucsc.edu/FAQ/FAQformat.html#format3).

### INTRODUCTION

Several genome annotation centers provide their data in GTF, BED, GFF3 etc. I have few programs 
they mainly deals with converting between GTF, BED and GFF3 formats. They are extensively tested 
with files from different centers like ENSEMBL, UCSC, JGI and NCBI AceView. Please follow the 
instructions below to clone these tools into your galaxy instance.

CONTENTS

Tool configuration files in *.xml format. 

    gtf_to_gff.xml
    gff_to_gtf.xml
    bed_to_gff.xml
    gff_to_bed.xml
    gbk_to_gff.xml
    bed_to_gff.xml
    
Python based scripts. 

    gtf_to_gff.py: convert data from GTF to valid GFF3.
    gff_to_gtf.py: convert data from GFF3 to GTF.
    bed_to_gff.py: convert data from a 12 column UCSC wiggle BED format to GFF3.
    gff_to_bed.py: convert gene transcript annotation from GFF3 to UCSC wiggle 12 column BED format.
    gbk_to_gff.py: convert data from genbank format to GFF. 
    GFFParser.py: Parse GFF/GTF files.  
    helper.py: Utility functions.

test-data: Test data set. (move to your galaxy_root_folder/test-data/)
    
    You may need to move the test files into your test-data directory so galaxy can find them. 
    If you want to run the functional tests eg as: 

    exmaple: 
    sh run_functional_tests.sh -id fml_gtf2gff

REQUIREMENTS

    python 

COMMENTS/QUESTIONS 

I can be reached at vipin [at] cbio.mskcc.org 

LICENSE

Copyright (C) 2009-2012 Friedrich Miescher Laboratory of the Max Planck Society
              2013-2015 Memorial Sloan Kettering Cancer Center

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

COURTESY

To the Galaxy Team.
