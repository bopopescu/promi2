# Scripts in this repo

## Multipurpose
- `utils.py` - module of utility functions
- `label.py` - adds miRNA labelling: intragenic/intergenic/unknown/NA based on an annotation file
- `plots.py` - generate plots based on miRNA labelling

## Creating the training set
1. `tc_normalization.py` - performs rle normalization on tag counts (used for creating the negative training set)
2. `create-training-set.py` - using script 1, creates the positive and negative training set

```
## You need to first update "config.ini" and "tcnorm.ini"

python2.7 create-training-set.py -f '/path/to/*.gff' -o ../Testout-tset -c config.ini
```
**Note:** `-f` globs all sample.gff files from which to pick the positive training set; you will need quotations around `*`; and an example of the input format is found [here](../test/test.gff).

## Training the model
1. `training.py` - performs training to get model parameters

```
## You need to first update "config.ini"

python2.7 training.py -i ../Test-tset/TrainingSet.gff -c config.ini
```
**Note:** See below if you want to enable training with the correlation feature

## Using the model
1. `features.py` - performs sequence feature extraction of CpG content, conservation, and TATA box affinity
2. `mirna_proximity.py` - performs miRNA proximity feature extraction (requires the use of miRBase files)
3. `correlation.py` - **new feature**; performs correlation feature extraction
                      (requires 'srnaseqmatrix' and 'cageseqmatrix', in the config file, to be filled;
                       see instructions below for enabling this feature)
4. `gff_unify_features.py` - performs `bedtools intersect ... -s -f 1 -r -wao`
5. `promirna.py` - module of core mathematical formulas from PROmiRNA
6. `promi2.py` - runs scripts 1-5 of "using the model"

```
## You need to first update "config.ini"

python2.7 promi2.py -i ../test/test.gff -o ../Testout-promi2 -c config.ini
```
**Note:** `-p` enables plotting of results with `plot.py` (plotting dependencies need to first be installed)

- **About: the gff input file**
 - This is the file containing the putative TSSs you want to classify
 - The format consists of 9 tab-delimited columns:
  - 1. chromosome
  - 2. source (field is not used)
  - 3. feature (field is not used)
  - 4. start of 1kb region around TSS (-500 of mid putative TSS)
  - 5. stop of 1kb region around TSS (+500 of mid putative TSS)
  - 6. normalized tag counts
  - 7. strand
  - 8. '.' (field is not used)
  - 9. 'start:<start>;stop:<stop>' of putative TSS
 - An example of the format is shown [here](../test/test.gff)

- **About: the config.ini file**
 - You will need to update the paths in the [config](config.ini) file

**Tip:** When you already have the extracted-features.gff file ([example](../test/test-features.gff)), you can run PROmiRNA directly with `-f`.
```
python2.7 promi2.py -i ../test/test-features.gff -f -o ../Testout-promi2predict
```
_(Note: In the config file, you will still need to fill in `params =` and `features =`)_


## Summarizing results
1. `summarize-results.py` - summarize promi2 predictions from a series of samples

```
python2.7 summarize-results.py -f '../test/test-summary/Predictions.*' -o ../Testout-summary.gff -s
```
**Note:** `-s` is to print only results where mirna_prox score is > 0


* * *
## Using the model2: custom
- `custom.py` - given a gff input file ([example](../test/test-custom.gff)) containing a list of positions (and possibly miRNA partners in column 9), performs feature extraction, calculation of probabilities, labelling of positions, and plotting of results

```
## You need to first update "config.ini" as well as "tcnorm.ini"
python2.7 custom.py -i ../test/test-custom.gff -o ../Testout-custom -p
```
**Note:** `-p` enables plotting of results with `plot.py` (plotting dependencies need to first be installed)

**Note:** `-m` designate that the infile has miRNA partners in column 9
(and that miRNA partners do not been to be found)


* * *
One main difference between `promi2.py` and `custom.py` is that `promi2.py` assumes you start with normalized CAGE tag counts in a matrix, whereas `custom.py` pulls the tag counts from raw data and normalizes them.

* * *
# Enabling Correlation
## Using correlation feature with promi2/custom
1. At the "**Creating the training set**" step, you will need to create a dataset which contains the extracted correlation
 - Make sure `srnaseqmatrix =` and `cageseqmatrix =`, in the [config](config.ini) file, is filled in
2. At the "**Training the model**" step, you will neet to tell the program to include correlation as a feature
 - In `trainingfeatures =`, change `cpg,cons,tata,mirna_prox` to `cpg,cons,tata,mirna_prox,corr`
 - Also make sure to train with training set with the extracted correlations
   (e.g. "TrainingSet-corr.gff" -- not "TrainingSet.gff")
3. At the "**Using the model**" step, you will need to tell promi2 (or custom) to use the corrlelation feature as well as the new parameters
 - Update `params =` to the resulting `*finalparams` file generated in the previous step
 - In `features =`, change `cpg,cons,tata,mirna_prox` to `cpg,cons,tata,mirna_prox,corr`
 - **Note:** promi2 also expects `cageseqmatrix =` to be filled out


