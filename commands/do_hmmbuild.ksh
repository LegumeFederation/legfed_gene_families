#!/usr/bin/env ksh93
set -o errexit
set -o nounset

INDIR=$1
OUTDIR=$2
JOBMAX=$3

for path in $INDIR/*; do
  file=`basename $path`;
  hmmbuild -n $file $OUTDIR/$file $path &
done
wait

