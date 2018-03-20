#!/usr/bin/env ksh93
set -o errexit
set -o nounset

HMMDIR=$1
HMMEMITDIR=$2
JOBMAX=$3

for path in $HMMDIR/*; do
  file=`basename $path`;
  hmmemit -c -o $HMMEMITDIR/$file $path &
done
wait

