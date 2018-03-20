#!/usr/bin/awk -f
#
# NAME
#   top_4_lines.awk - Filter tabular blast output to top two hits per query.
# 
# SYNOPSIS
#   top_4_lines.awk [INPUT_FILE(S)]
#     
# OPERANDS
#     INPUT_FILE
#         A file containing tabular blast output (or similar); the important
#         field is 0 (query ID in blast output).

BEGIN { MAX = 4 } 

$1 == prev && count < MAX { print; count++ }
$1 != prev { print; count = 1; prev = $1 }
