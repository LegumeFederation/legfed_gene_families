#!/usr/bin/env ksh93
set -o errexit
set -o nounset

HMMALIGNDIR=$1
HMMTRIMDIR=$2
JOBMAX=$3

for path in $HMMALIGNDIR/*; do
  file=`basename $path`
  perl -ne 'if ($_ =~ />/) {print $_} else {$line = $_; $line =~ s/[a-z]//g; print $line}' $path |
    sed '/^$/d' > $HMMTRIMDIR/$file &
done
wait

