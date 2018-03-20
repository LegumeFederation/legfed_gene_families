#!/usr/bin/awk -f
#
# NAME
#   top_line_no_self.awk - Filter tabular blast output to top hit per query
#                          ... but skipping gene self-matches
# 
# SYNOPSIS
#   ./top_line_no_self.awk [INPUT_FILE(S)]
#     
# OPERANDS
#     INPUT_FILE
#         A file containing tabular blast output (or similar); the important
#         field is 0 (query ID in blast output).

BEGIN { MAX = 1 } 

$1 == prev && count < MAX && $1 != $2 { print; count++ }
$1 != prev && $1 != $2 { print; count = 1; prev = $1 }
