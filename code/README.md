# Scripts in this repo

## Generic
- `utils.py` - module of utility functions

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

python2.7 training.py -i ../test/test-trainset.gff
```

## Using the model
1. `features.py` - performs sequence feature extraction of CpG content, conservation, and TATA box affinity
2. `mirna_proximity.py` - performs miRNA proximity feature extraction (requires the use of miRBase files)
3. `gff_unify_features.py` - performs `bedtools intersect ... -s -f 1 -r -wao`
4. `promirna.py` - module of core mathematical formulas from PROmiRNA
5. `promi2.py` - runs scripts 1-4 of "using the model"

```
## You need to first update "config.ini"

python2.7 promi2.py -i ../test/test.gff -o ../Testout-promi2 -c config.ini
```

## Summarizing results
1. summarize-results.py
