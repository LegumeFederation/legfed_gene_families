#!/usr/bin/env ksh93

set -o errexit
set -o nounset 

HMMALIGNTRIMMEDDIR=$1
TREEDIR=$2
JOBMAX=$3

# By default, FastTreeMP uses all machine cores, e.g. 64 on lathyrus. 
# It is more efficient to run more jobs on one core each by setting an environment variable.
OMP_NUM_THREADS=1
export OMP_NUM_THREADS
for path in $HMMALIGNTRIMMEDDIR/*; do
  file=`basename $path`
  FastTreeMP -quiet $path > $TREEDIR/$file &
done
wait


