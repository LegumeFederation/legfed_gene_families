#!/usr/bin/env ksh93
set -o errexit
set -o nounset

BLASTDBDIR=blastdb
BLASTOUTDIR=blastout
FASTADIR1=$1
FASTADIR2=$2
JOBMAX=$3

for path1 in $FASTADIR1/*; do
  file1=`basename $path1`
  for path2 in $FASTADIR2/*; do
    file2=`basename $path2`
    #if [[ ! $file1 > $file2 ]]; then  # A x B, A x A, B x B -- but not B x A
    #if [[ $file1 > $file2 ]]; then     # B x A -- but not A x B, A x A, B x B
      echo "WORKING ON $file1 $file2";
     blastp -query $path1 -db $BLASTDBDIR/$file2 -out $BLASTOUTDIR/$file1.x.$file2.blp -evalue 1e-5 -outfmt 6 &
    #fi
  done
done

