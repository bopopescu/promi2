_For more information, look in the supplement of PROmiRNA_

## [miRBase](http://www.mirbase.org)
- For a list of (e.g. position) information regarding the `miRNA_primary_transcript`
- The version of miRBase included (v20) is not the latest release -- v20 is the release used in the PROmiRNA publication

## PhastCons
- Files: `track_000.chr*.wib`
- Average PhastCons conservation score of the promoter region across vertebrates
 on the 46-way alignment (http://hgdownload.cse.ucsc.edu/goldenPath/hg19/phastCons46way/vertebrate/)
- You can find the files included in the ExternalData of http://promirna.molgen.mpg.de

## TATA TRAP affinity
- File: `TATA_box_jaspar.psem`
- position-specific mismatch energy matrix (psem)
- information derived from the Jaspar database (http://jaspar.genereg.net/, Jaspar ID: MA0108.2)
- You can find the file included in the ExternalData of http://promirna.molgen.mpg.de

## Files related to the human genome
- Files: `hg19.fa`, `hg19_repeats.bed`, `hg19.chrom.sizes`, `Homo_sapiens.GRCh37.75.gtf`
- You can find the files included in the ExternalData of http://promirna.molgen.mpg.de


* * *
## [miRIAD](http://www.bioinfo.mochsl.org.br/miriad): labelling miRNA as	intragenic/intergenic
- File: [miriad_human_labels_v2014.gff](miriad_human_labels_v2014.gff)
- The annotation file was [download](http://www.bioinfo.mochsl.org.br/miriad/downloads)
  and processed to gff format for ...
 1. use in [label.py](../code/label.py) and
 2. subsequent use for plotting in [plots.py](../code/plots.py)

- **Note**: it is acceptible to use another annotation file for [label.py](../code/label.py) given that the new file is in the same format as [miriad_human_labels_v2014.gff](miriad_human_labels_v2014.gff)