Written by: Robin Andersson

## Input
The `tc_quantify.sh` takes six arguments:

1. TC file (.txt file downloaded from FANTOM, I have tried it on tc.decompose_smoothing_merged.ctssMaxCounts3.clustername_update__description.txt (it only uses the **chr:start..stop,strand info** from the first column).
2. A file with library ids (not really important that it is the correct ids, you could use a unique integer per row, a remnant from my [Robin's] usage)
3. A file with the files (ctss.bed files) to be used including their paths
4. A file with the tag counts per library
5. A file with normalization factors for RLE normalization
6. out path

2-5 should be in the same order, one library per line.

## Output
It produces five files:

1. A bed file transformed from the TC .txt file
2. A *_expression.matrix file with the first six columns as the bed file (but resorted) and then the tag counts per library in the same order as in your library file
3. A *_expression_tpm_rle.matrix file with the first six columns as the bed file (but resorted) and then the tpm and rle normalized tag counts per library in the same order as in your library file
4. A *_expression_max_count.bed file that is the bed file plus a seventh column with the max counts across all libraries
5. A *_expression_max_count.bed file that is the bed file plus a seventh column with the max tpm (rle normalized) across all libraries

The fourth column of all out files holds the TC position id so you can map it back to the original data.

## Example
```
nohup nice ./tc_quantify.sh ${txt} ${libs} ${files} ${counts} ${norm_factors} ${outdir}/ > ${log_nohup.out} 2>&1 &
```
