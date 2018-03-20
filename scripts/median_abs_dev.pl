#!/usr/bin/env perl

# PROGRAM: median_abs_dev.pl
# S. Cannon 2017
# See description under Usage, or with median_abs_dev.pl -h

use strict;
use warnings;
use Getopt::Long;

my $scale_factor = 1.4826; # see https://en.wikipedia.org/wiki/Median_absolute_deviation
my $filename;
my $columns_in = 1;
my $help;
GetOptions (
  "filename:s"   => \$filename,
  "columns_in:i" => \$columns_in,
  "help"         => \$help,
  );

my $usage = <<EOS;
Given a stream of numbers OR a stream of IDs and numbers (sorted by ID), return the median 
and scaled median absolute deviation (MAD) - for all numbers in the case of the input 
consisting of a stream of numbers, or per ID in the case of input consisting of IDs and numbers.
  Example of the median and MAD output format:
    #filename-or-ID  med-MAD  median   med+MAD  100*(med-MAD)/med  count_of_nums
    L_473CLT         438.52   466.10   493.68   94.08              10
    L_475ZRT         955.44   989.10   1022.76  96.60              55

Usage: cat STRING_OF_NUMS or STRING_OF_IDS_AND_NUMS | median_abs_dev.pl -filename FILE

Options:
  -columns_in  1 or 2; default 1 (string of numbers); 2 requires input such as IDENTIFIER 123, sorted on col1.
  -filename    FILE -- NOTE: the name may be passed in for use in the output, for the one-column option
               (i.e. string-of-numbers), but data must come in via STDIN. The filename is not used for 
               opening the file; just for information in the output. This is because the data is likely
               generated via a unix pipeline rather than in a static file.
  -help        This message
EOS

my $wrong_input = "columns_in was specified as $columns_in but must be 1 or 2.";

die "$usage\n" if $help; 
die "$usage\n$wrong_input\n" if ( $columns_in!=1 and $columns_in!=2 );

my @nums;
my %nums_HoA;

while (<>) {
  chomp;
  next if $_ =~ /^$|^#/;
  if ($columns_in == 1) {
    push @nums, $_;
  }
  elsif ($columns_in == 2) {
    my ($ID,$num) = split(/\s+/,$_);
    # print "$ID, $num\n";
    push @{$nums_HoA{$ID}}, $num;
  }
}

if ($columns_in == 1) {
  my $median = median(@nums);

  my $count_of_nums = scalar(@nums);
  
  my @diffs_fr_median;
  for my $elt (@nums){
    push @diffs_fr_median, abs($median-$elt);
  }
  
  if (not defined($filename)) { $filename="NONE"}
  
  my $scaled_MAD = $scale_factor * median(@diffs_fr_median);
  my $median_minus = $median-$scaled_MAD;
  my $median_plus = $median+$scaled_MAD;
  my $pct_MAD_vs_median = 100*($median-$scaled_MAD)/$median;
  
  printf "%s\t%.2f\t%.2f\t%.2f\t%.2f\t%u\n", 
          $filename, $median_minus, $median, $median_plus, $pct_MAD_vs_median, $count_of_nums;
}
elsif ($columns_in == 2) {
  foreach my $ID ( keys %nums_HoA ) {
    my @nums_for_this_ID = @{ $nums_HoA{$ID} }; # give array a nicer name, since we use it a few times
    my $count_of_nums = scalar(@nums_for_this_ID);

    my $median = median(@nums_for_this_ID);
    
    my @diffs_fr_median;
    for my $elt (@nums_for_this_ID){
      push @diffs_fr_median, abs($median-$elt);
    }
    
    my $scaled_MAD = $scale_factor * median(@diffs_fr_median);
    my $median_minus = $median-$scaled_MAD;
    my $median_plus = $median+$scaled_MAD;
    my $pct_MAD_vs_median = 100*($median-$scaled_MAD)/$median;
    
    printf "%s\t%.2f\t%.2f\t%.2f\t%.2f\t%u\n", 
            $ID, $median_minus, $median, $median_plus, $pct_MAD_vs_median, $count_of_nums;
  }
}

##########
sub median
{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}

__END__
VERSIONS
v01 2017-12-31 initial; works OK on a stream of incoming numbers, e.g. for match scores for a gene family
v02 2018-02-22 add option for two-column input, e.g. for gene family ID & match score


