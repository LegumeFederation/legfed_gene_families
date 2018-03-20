#!/usr/bin/env ksh93
set -o errexit
set -o nounset

HMMDIR=$1
FASDIR=$2
TBLOUTDIR=$3
JOBMAX=$4
for path in $FASDIR/*; do
  file=`basename $path`;
  hmmscan -o /dev/null --tblout $TBLOUTDIR/$file -E 1e-4 --noali --cpu 1 $HMMDIR/$file $FASDIR/$file &
done
wait

