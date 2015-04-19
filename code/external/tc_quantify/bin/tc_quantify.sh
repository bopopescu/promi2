#!/bin/bash

TCS=$1
LIBRARIES=$2
LIBRARYFILES=$3
LIBRARYCOUNTS=$4
LIBRARYNORMFACTORS=$5
OPATH=$6

STUB=$( echo $TCS | sed -e "s/\.txt//" )
STUB=`basename $STUB`
STUB=${OPATH}${STUB}
BINDIR="`dirname \"$0\"`"
REGIONS=${STUB}.bed

## Needed scripts in the same directory as this file:
## rle_tpm_transform.pl
## to_bed.pl

## transform tc file to bed file
${BINDIR}/to_bed.pl $TCS > $REGIONS

## Quantify expression
cp ${REGIONS} ${STUB}_expression.matrix
for FILE in $( cat ${LIBRARYFILES} ); do
	echo ${FILE}
	intersectBed -wao -s -a ${REGIONS} -b ${FILE} | sort -k 1,1 -k 2,2n | cut -f -6,11 | awk 'BEGIN{OFS="\t"} {v = $7; if (v <= 0) v = 0; print $1,$2,$3,$4,$5,$6,v}' | groupBy -g 1,2,3,4,5,6 -c 7 -o sum | sort -k 1,1 -k 2,2n | cut -f 7 > ${STUB}_expression.tmp
	paste ${STUB}_expression.matrix ${STUB}_expression.tmp > ${STUB}_expression.tmp2
	mv ${STUB}_expression.tmp2 ${STUB}_expression.matrix
done
rm ${STUB}_expression.tmp


## TPM normalize
${BINDIR}/rle_tpm_transform.pl ${STUB}_expression.matrix ${LIBRARIES} ${LIBRARYCOUNTS} ${LIBRARYNORMFACTORS} ${LIBRARIES} 6 | paste ${REGIONS} - > ${STUB}_expression_tpm_rle.matrix

## max values (extra column added)
awk 'function max(a,b) {if (a > b) {return a} else {return b}} BEGIN{OFS="\t"}{m=0; for (i = 7; i <= NF; i++) {if ($i != "NA") m = max(m,$i)}; print $1,$2,$3,$4,$5,$6,m}' ${STUB}_expression.matrix > ${STUB}_expression_max_count.bed
awk 'function max(a,b) {if (a > b) {return a} else {return b}} BEGIN{OFS="\t"}{m=0; for (i = 7; i <= NF; i++) {if ($i != "NA") m = max(m,$i)}; print $1,$2,$3,$4,$5,$6,m}' ${STUB}_expression_tpm_rle.matrix > ${STUB}_expression_max_tpm.bed

