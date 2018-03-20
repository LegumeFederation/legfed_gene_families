#!/usr/bin/awk -f
#
# NAME
#   top_hmmscan_line.awk - Filter tabular hmmscan output to top hit per query.
# 
# SYNOPSIS
#   ./top_hmmscan_line.awk [INPUT_FILE(S)]
#     
# OPERANDS
#     INPUT_FILE
#         A file containing tabular blast output (or similar); the important
#         field is 3 (query ID in hmmscan output).
#
# EQUIVALENT TO
#   awk '$3==prev && ct<2 {print; ct++} $3!=prev {print; ct=1; prev=$3}' FILE(s)

BEGIN { MAX = 1 } 

$3 == prev && count < MAX { print; count++ }
$3 != prev { print; count = 1; prev = $3 }
