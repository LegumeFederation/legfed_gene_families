#!/usr/bin/env perl

# Program: filter_by_median_abs_dev.pl
# Program description: see usage message.
# S. Cannon 2018

use strict;
use warnings;
use strict;
use Getopt::Long;
use List::Util qw(max min);

my ($hmm_file, $median_file, $no_MAD);
my $ceil_pct_of_median = 95;
my $floor_pct_of_median = 80;
my $hmmsearch=0;
my $help;

GetOptions (
  "hmm_file=s" =>          \$hmm_file,
  "median_file=s" =>       \$median_file,
  "cell_pct_of_median:i" =>  \$ceil_pct_of_median,
  "floor_pct_of_median:i" => \$floor_pct_of_median,
  "no_MAD"                => \$no_MAD,
  "hmmsearch"             => \$hmmsearch,
  "help:s" =>              \$help
);

my $usage = <<EOS;
  Simplified concept: given a file of best hmmscan matches of a set of genes to a set of HMMs, 
  use a per-family score threshold to decide whether each gene matches its candidate HMM.

  Read a tabular hmmscan output file and a file of median and scaled median-absolute-deviation 
  (MAD) values calculated for sequences that were used to build the HMM for each family. 
  The median and MAD values are calculated from a family HMM file by median_abs_dev.pl
  Report sequences belonging to each family, with scores at or above the lesser of 
  (the median - scaled MAD) and (a provided percentage of the median).

  Example of the median and MAD format from median_abs_dev.pl:
    #filename-or-ID  med-MAD  median   med+MAD  100*(med-MAD)/med  count_of_nums
    L_473CLT         438.52   466.10   493.68   94.08              10
    L_475ZRT         955.44   989.10   1022.76  96.60              55

  Usage: perl $0 -hmm_file FILE -median_file FILE [options] > out_file
  
  -hmm_file:     (required) tabular HMM file; output of hmmscan
  -median_file:  (required) tabular file with five fields; ID 1st, median 3rd (above)
  -floor_pct_of_median (default $floor_pct_of_median) 
        Percent of median to use as lower cutoff, if the scaled minimum absolute deviation 
        (MAD) from second column of median_file is lower than this.
  -ceil_pct_of_median  (default $ceil_pct_of_median) 
        Percent of median to use as upper cutoff, to avoid cases in which 
        the scaled MAD is equal to (100% of) the median (i.e. too stringent).
  -no_MAD        If set, use only the floor_pct_of_median and median to set the cutoff
                 ... rather than median-(scaled MAD)
  -help:          display this help message
EOS

die "\n$usage\n" if ($help);
die "\n$usage\nPlease provide hmm_file\n" unless (defined($hmm_file));
die "\n$usage\nPlease provide median_file\n" unless (defined($median_file));

open(my $MEDFH, "<", $median_file) or die "can't open $median_file: $!";
open(my $HMMFH, "<", $hmm_file) or die "can't open $hmm_file: $!";

my (%medians, %med_less_MADs, %floors, %ceils, %obs_pct_of_med);

while (<$MEDFH>) {
  chomp;
  next if $_ =~ /^$/;
  if (scalar(split(/\s+/, $_)) < 5) {
    warn "fewer than five elements in median_file: [$_]\n";
  }
  my ($ID, $med_less_MAD, $median, $med_plus_MAD, $obs_pct_of_med) = split(/\s+/, $_);
  $medians{$ID} = sprintf "%.1f", $median;
  $med_less_MADs{$ID} = sprintf "%.1f", $med_less_MAD;
  $floors{$ID} = sprintf "%.1f", (($floor_pct_of_median/100) * $median);
  $ceils{$ID} = sprintf "%.1f", (($ceil_pct_of_median/100) * $median);
  $obs_pct_of_med{$ID} = sprintf "%.1f", $obs_pct_of_med;
  #print "$ID\t[$floors{$ID}\t$ceils{$ID}]\t$med_less_MADs{$ID}" . 
  #      "\t$medians{$ID}\t$obs_pct_of_med{$ID}\n";
}
close $MEDFH;

# read hmm file and filter based on medians file 
print "#S_ID____\tQ_ID_________\tscore\tcutoff\t[floors\tceils]" . 
      "\tadjMAD\tfam_pct\tverdict\tpctOfCutoff\tpctOfMedian\n";
while (<$HMMFH>) {
  chomp;
  my $line = $_;
  next if $line =~ /^$/;
  next if $line =~ /^#/;
  my ($S_ID, $S_acc, $Q_ID, $Q_acc, $E_val, $score, @rest);
  if ($hmmsearch) {
    ($Q_ID, $Q_acc, $S_ID, $S_acc, $E_val, $score, @rest) = split(/\s+/, $line);
  }
  else {
    ($S_ID, $S_acc, $Q_ID, $Q_acc, $E_val, $score, @rest) = split(/\s+/, $line);
  }
  # Set cutoff as the median-MAD unless median-MAD is too low or too high.
  # If $no_MAD is set, just use $floor_pct_of_median*$medians{$S_ID})/100
  my $cutoff;
  if (defined($medians{$S_ID})) {
    if ( $no_MAD || $med_less_MADs{$S_ID} < ($floor_pct_of_median*$medians{$S_ID})/100 ) {
      $cutoff = sprintf "%.1f", ($floor_pct_of_median*$medians{$S_ID})/100;
    }
    elsif ($med_less_MADs{$S_ID} > ($ceil_pct_of_median*$medians{$S_ID})/100 ) {
      $cutoff = sprintf "%.1f", ($ceil_pct_of_median*$medians{$S_ID})/100;
    }
    else {
      $cutoff = $med_less_MADs{$S_ID};
    }
  
    my $verdict;
    if ( $score >= $cutoff ) { $verdict = "MATCHES" }
    else { $verdict = "doesnt" }
  
    print "$S_ID\t$Q_ID\t$score\t$cutoff\t$floors{$S_ID}\t$ceils{$S_ID}" . 
          "\t$med_less_MADs{$S_ID}\t$obs_pct_of_med{$S_ID}\t$verdict\t" .
          sprintf("%.1f", 100*$score/$cutoff), "\t" . 
          sprintf("%.1f", 100*$score/$medians{$S_ID}), "\n";
  }
  else {
    warn "WARNING: median not defined for $S_ID\n";
  }
}
close $HMMFH;

__END__

Versions
v01 2018-01-02



