#!/usr/bin/env ksh93

# Wraps add_align_pct_to_blast.pl, which reads a BLAST file and tables of query and subject lengths.
# Adds four fields to BLAST output file: query len, subj len, pct of query aligned, pct of subj aligned. 
# NOTE: Filenames in the BLAST results are assumed to have the format medtr.x.lotja.blp
# where the query and subject species names each have five characters and are separated by .x.
set -o errexit
set -o nounset

LEN_DIR=$1
COV_DIR=$2
COV_PCT=$3
ID_PCT=$4
JOBMAX=$5

for path in blastout/*; do
  file=`basename $path .blp`;
  qry=${file:0:5};
  subj=${file:8:5};
  add_align_pct_to_blast.pl -bl $path -q 02_pep_len/$qry -s 02_pep_len/$subj -cov $COV_PCT -id $ID_PCT \
    > blastout_cov/$qry.x.$subj &
done
wait

