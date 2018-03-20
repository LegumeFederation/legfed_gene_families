#!/usr/bin/env ksh93

set -o errexit
set -o nounset 

JOBMAX=$1   # e.g. 20
WORKDIR=$2  # e.g. /scratch/scannon/legume_genefams3
ALIGNDIR=$3 # e.g. 33_hmmalign_trim2
TREEDIR=$4  # e.g. 36_trees_RAxML
ROOTS=$5    # e.g. "prupe,cucsa,arath,vitvi,solly" or "orysa,sorbi,ambtr"
FILEPAT=$6  # e.g. "*" or "L.*" or "E.H*"

for path in $ALIGNDIR/$FILEPAT; do
  file=`basename $path`;
  echo "WORKING ON $file";
  #echo root_tree_by_species.pl -align $path -out $WORKDIR/$TREEDIR -roots $ROOTS -exe &
  root_tree_by_species.pl -align $path -out $WORKDIR/$TREEDIR -roots $ROOTS -exe &
done
wait

