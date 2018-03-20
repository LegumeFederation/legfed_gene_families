#!/usr/bin/perl -w

# Program: split_table_to_files.pl

# S. Cannon Jan 2008

use strict;

my $usage_string = <<EOS;
Usage: perl $0  input_table   out_dir   field_number
  Input is a table to be split into sub-tables, with the sub-tables 
  named by values in the sorted input table.
  
  The field_number is the field to use for naming output files
  (one-indexed). Values in this field should be suitable for file-naming, 
  and the file should be sorted by this column.
  
EOS

my ($input_table, $out_dir, $field_number) = @ARGV;

die "\n$usage_string\n" unless (scalar(@ARGV) == 3);

my $params = <<EOS;
input_table:   $input_table
out_dir:       $out_dir
field_number:  $field_number
EOS

print "$params\n";

$out_dir =~ s{/$}{};
opendir (DIR, "$out_dir") or mkdir($out_dir); close DIR;

open (my $IN, '<', $input_table) or die "can't open in $input_table:$!";

my %seen;

while (my $line = <$IN>) {
  chomp $line;  
  next if $line =~ m/^#|^$/;
  
  my @elements = split /\t/, $line;
  my $file_out = $elements[$field_number-1];
  
  if ($seen{$file_out}) {
    print OUT "$line\n";
  }
  else {
    $seen{$file_out}++;
    my $out_file = "$out_dir/$file_out";
    open (OUT, "> $out_file") or die "can't open out $out_file: $!";
    # print "$out_file\n";
    print OUT "$line\n";
  }
}

# Version
# 0.01 Jan03'08 basic; works OK.
# 0.02 Feb25'08 Change @ARGV input test to allow $field_number = 0
# 0.03 Apr04'08 Change from zero-indexing in the field number to one-indexing
# 0.04 Nov06'08 Add chomp, in case file name is taken from last field
