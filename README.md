- This is a component of my rotation project (Jan-Apr 2015)
- It is a reimplementation of [PROmiRNA](http://promirna.molgen.mpg.de)

## Dependencies
- Python >= 2.7
- R >= 2.12.1
- Perl >= 5.12
- [Bedtools](http://bedtools.readthedocs.org/en/latest/)
- [ANNOTATE 3.04 (TRAP: Transcription Factor Affinity Prediction)](http://www.mybiosoftware.com/trap-3-04-transcription-factor-affinity-prediction.html)
- BIO::Graphics perl module (for e.g. Bio/Graphics/Wiggle.pm)

### Additional input files
- miRBase (I have provided v20): `hsa.gff2`, `mirna.txt`, `mirna_context.txt`

These are large files; you can find them included in the [ExternalData](http://promirna.molgen.mpg.de) of PROmiRNA
- PhastCons: `track_000.chr*.wib`
- TATA box affinity: `TATA_box_jaspar.psem`
- Files related to the human genome:
  - `hg19.fa`
  - `hg19_repeats.bed`
  - `hg19.chrom.sizes`
  - `Homo_sapiens.GRCh37.75.gtf`

## Input files
### 1. File containing the putative TSSs you want to classify
This file is in gff format -- so the following 9 tab-delimited columns:

1. chromosome
2. source (field is not used)
3. feature (field is not used)
4. start of 1kb region around TSS (-500 of mid putative TSS)
5. stop of 1kb region around TSS (+500 of mid putative TSS)
6. normalized tag counts
7. strand
8. '.' (field is not used)
9. 'start:<start>;stop:<stop>' of putative TSS

Example shown [here](test/test.gff).

### 2. Config file
There is also a [config](code/config.ini) file you will need to update.

## Reference
1. Marsico A, Huska MR, Lasserre J, Hu H, Vucicevic D, Musahl A, Orom U, Vingron M. (2013). [PROmiRNA: a new miRNA promoter recognition method uncovers the complex regulation of intronic miRNAs](http://genomebiology.com/2013/14/8/R84). Genome Biology, 14(8):R84. doi: 10.1186/gb-2013-14-8-r84.
2. Kozomara A, Griffiths-Jones S. (2014). [miRBase: annotating high confidence microRNAs using deep sequencing data](http://nar.oxfordjournals.org/content/42/D1/D68). Nucleic Acids Research, 42(Database issue):D68-73. doi: 10.1093/nar/gkt1181.
3. Thomas-Chollier M, Hufton A, Heinig M, O'Keeffe S, Masri NE, Roider HG, Manke T, Vingron M. (2011). [Transcription factor binding predictions using TRAP for the analysis of ChIP-seq data and regulatory SNPs](http://www.nature.com/nprot/journal/v6/n12/full/nprot.2011.409.html). Nature Protocols, 3;6(12):1860-9. doi: 10.1038/nprot.2011.409.
