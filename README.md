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

##### REQUIREMENTS

    python2.6 or 2.7 and biopython  

    Galaxy should be able to automatically install biopython via Galaxy toolshed.

##### COMMENTS/QUESTIONS 

I can be reached at vipin [at] cbio.mskcc.org 

##### LICENSE

Copyright (c) 2009-2012, Friedrich Miescher Laboratory of the Max Planck Society

              2013-2015, Memorial Sloan Kettering Cancer Center

              Vipin T Sreedharan <vipin@cbio.mskcc.org>  
All rights reserved.

Licensed under the BSD 2-Clause License: <http://opensource.org/licenses/BSD-2-Clause>
    
    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:
    
        * Redistributions of source code must retain the above copyright notice,
          this list of conditions and the following disclaimer.
    
        * Redistributions in binary form must reproduce the above copyright notice,
          this list of conditions and the following disclaimer in the documentation
          and/or other materials provided with the distribution.
    
    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
    ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
    FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

##### COURTESY

To the Galaxy Team.
