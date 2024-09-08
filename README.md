<code> 
usage: sashimi.py [-h] -c COORDINATES -b BAM -g GTF [-o OUTPUT] [-s STRAND] [-v VCF] [-sa SPLICEAI] [-gb {grch37,grch38}] [-t TEMP]
                  [-r REFERENCE]

Browsable html sashimi

options:
  -h, --help            show this help message and exit
  -c COORDINATES, --coordinates COORDINATES
  -b BAM, --bam BAM
  -g GTF, --gtf GTF
  -o OUTPUT, --output OUTPUT
  -s STRAND, --strand STRAND
  -v VCF, --vcf VCF
  -sa SPLICEAI, --spliceai SPLICEAI {False, True}
  -gb {grch37,grch38}, --genomebuild {grch37, grch38} // needed for SpliceAI
  -t TEMP, --temp TEMP // temp folder needed to store subset of vcf and SpliceAI annotated vcf.
  -r REFERENCE, --reference REFERENCE // genome fasta, see 'spliceai' python package for download instructions.
<\code>
