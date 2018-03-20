#!/usr/bin/awk -f
#
# NAME
#   top_2_lines_no_self.awk - Filter tabular blast output to top 2 hits per query
#                          ... but skipping gene self-matches
# 
# SYNOPSIS
#   ./top_2_lines_no_self.awk [INPUT_FILE(S)]
#     
# OPERANDS
#     INPUT_FILE
#         A file containing tabular blast output (or similar); the important
#         fields are 1 and 2 (subject and query IDs in blast output).

BEGIN { MAX = 2 } 

$1 == prev && count < MAX && $1 != $2 { print; count++ }
$1 != prev && $1 != $2 { print; count = 1; prev = $1 }
