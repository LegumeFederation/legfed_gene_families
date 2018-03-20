#!/usr/bin/env ksh93

set -o errexit
set -o nounset 

HMMTRIMDIR1=$1
HMMTRIMDIR2=$2
LOGDIR=$3
MIN_DEPTH=$4
MIN_PCT_DEPTH=$5
MIN_PCT_ALIGNED=$6
JOBMAX=$7

echo "PARAMETERS: HMMTRIMDIR1 HMMTRIMDIR2 LOGDIR MIN_DEPTH MIN_PCT_DEPTH MIN_PCT_ALIGNED JOBMAX"
echo "PARAMETERS: $HMMTRIMDIR1 $HMMTRIMDIR2 $LOGDIR $MIN_DEPTH $MIN_PCT_DEPTH $MIN_PCT_ALIGNED $JOBMAX"

for path in $HMMTRIMDIR1/*; do
  file=`basename $path`
  filter_align.pl -in $path -out $HMMTRIMDIR2/$file -log $LOGDIR/$file \
                  -depth $MIN_DEPTH -pct_depth $MIN_PCT_DEPTH -min_pct_aligned $MIN_PCT_ALIGNED &
done
wait

