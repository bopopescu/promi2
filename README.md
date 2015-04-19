- This is a component of my rotation project (Jan-Apr 2015)
- It is a reimplementation of [PROmiRNA](http://promirna.molgen.mpg.de)

## Dependencies
- Python >= 2.7
- R
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

## Reference
1. Marsico A, Huska MR, Lasserre J, Hu H, Vucicevic D, Musahl A, Orom U, Vingron M. (2013). [PROmiRNA: a new miRNA promoter recognition method uncovers the complex regulation of intronic miRNAs](http://genomebiology.com/2013/14/8/R84). Genome Biology, 14(8):R84. doi: 10.1186/gb-2013-14-8-r84.
