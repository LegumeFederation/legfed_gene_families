#!/usr/bin/awk -f
#
# NAME
#   top_line_by_col2.awk - Filter tabular output to first instance of ID in 2nd column (e.g. blast target).
# 
# SYNOPSIS
#   ./top_line_by_col2.awk [INPUT_FILE(S)]
#     
# OPERANDS
#     INPUT_FILE
#         A file containing tabular blast output (or similar); the important
#         field is the second one (target ID in blast output).

BEGIN { MAX = 1 } 

$2 == prev && count < MAX { print; count++ }
$2 != prev { print; count = 1; prev = $2 }
