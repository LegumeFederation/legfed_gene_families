#!/usr/bin/env ksh93
set -o errexit
set -o nounset

FASTADIR=$1
HMMDIR=$2
HMMALIGNDIR=$3
JOBMAX=$4

for path in $FASTADIR/*; do
  file=`basename $path`;
  hmmalign --trim --outformat A2M --amino -o $HMMALIGNDIR/$file $HMMDIR/$file $FASTADIR/$file &
done
wait

