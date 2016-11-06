# NCCU bioinformatrics HomeWork 2 
## Desciption
Homeword 2 is a project about implementation of __global alignment__ and __local alignment__ of __protein sequence__ with different scoring scheme.<br>
In this folder, __hw2_104753035_ywliu.R__ is the main program, and user have to provide the following argument:<br>
    input:  path of input file (e.g. ./test.fasta.txt)
    score:  path of scoring file (e.g. pam250.txt, pam100.txt)
    aln: global/local
    gap_open: any numerical value (e.g. -10, -5)
    gap_extend: any numerical value (e.g. -2, -4)
    output: name of output file (e.g. ./result.fasta)
    gap: open/extend
Sxample Script:
    R hw2_104753035_ywliu.R --input test.fasta.txt --score pam250.txt --aln global --gap_open -10 --gap_extend -2 --output result.fasta --gap open

## Require R Package
- Biostrings

