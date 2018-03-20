#!/usr/bin/env ksh93
set -o errexit
set -o nounset

PEPDIR=$1 # 04_pair_peptide_fasta
CDSDIR=$2 # 04_pair_CDS_fasta
WORKDIR=$3 # 05_ks
MAINDIR=$4 # /scratch/scannon/lgf4
JOBMAX=$5 # 60

for path in $CDSDIR/*; do
  file=`basename $path`
  mkdir $MAINDIR/$WORKDIR/$file
  cd $MAINDIR/$WORKDIR/$file
  synonymous_calc.py $MAINDIR/$PEPDIR/$file $MAINDIR/$CDSDIR/$file \
    2>> $MAINDIR/$WORKDIR/$file.stderr 1>>$MAINDIR/$WORKDIR/$file.ks &
done
wait

