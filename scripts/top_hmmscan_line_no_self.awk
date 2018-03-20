#!/usr/bin/awk -f
#
# NAME
#   top_hmmscan_line_no_self.awk - Filter tabular hmmscan output to top hit per query
#                             ... but skipping self-matches.
# 
# SYNOPSIS
#   top_hmmscan_line_no_self.awk [INPUT_FILE(S)]
#     
# OPERANDS
#     INPUT_FILE
#         A file containing tabular HMM output (or similar); the important
#         fields are 1 and 3 (subject and query IDs in hmmscan output).

BEGIN { MAX = 1 } 

$3 == prev && count < MAX && $1 != $3 && $1 !~ /^#/ { print; count++ }
$3 != prev && $1 != $3 && $1 !~ /^#/ { print; count = 1; prev = $3 }
