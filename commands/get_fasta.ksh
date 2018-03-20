#!/usr/bin/env ksh93
set -o errexit
set -o nounset
INFASTA=$1
LISTDIR=$2
OUTDIR=$3
JOBMAX=$4
for path in $LISTDIR/*; do
  file=`basename $path`
  #echo "get_fasta_subset.pl -in $INFASTA -out $OUTDIR/$file -lis $path -c"
  get_fasta_subset.pl -in $INFASTA -out $OUTDIR/$file -lis $path -c &
done
wait

