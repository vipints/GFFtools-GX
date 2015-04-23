GFFtools-GX 
===========

A collection of tools for converting genome annotation between [GTF](https://genome.ucsc.edu/FAQ/FAQformat.html#format4), [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1), [GenBank](http://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html) and [GFF](https://genome.ucsc.edu/FAQ/FAQformat.html#format3).

##### INTRODUCTION

Several genome annotation centers provide their data in GTF, BED, GFF and GenBank format. I have few programs, they mainly deals with converting between GTF, BED GenBank and GFF formats. They are extensively tested with files from different centers like [ENSEMBL](http://www.ensembl.org), [UCSC](https://genome.ucsc.edu/), [JGI](http://genome.jgi.doe.gov/) and [NCBI AceView](http://www.ncbi.nlm.nih.gov/IEB/Research/Acembly/HelpJan.html). These programs can be easily integrated into your galaxy instance.

##### CONTENTS

Included utilities are: 

 BED-to-GFF: convert data from a 12 column UCSC wiggle BED format to GFF
 GBK-to-GFF: convert data from genbank format to GFF
 GFF-to-BED: convert data from GFF to 12 column BED format
 GFF-to-GTF: convert data from GFF to GTF 
 GTF-to-GFF: convert data from GTF to valid GFF

test-data: Test data set. (move to your galaxy-root-folder/test-data/)
    
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
