#!/usr/bin/env perl

# Program: filter_by_ks.pl
# Program description: see usage message.
# S. Cannon 2017

use strict;
use warnings;
use strict;
use Getopt::Long;

my ($ks_file, $cutoff_file, $out_file, $verbose);
my $scale_factor=1.5;
my $min_transform=10;
my $max_transform=100;
my $help;

GetOptions (
  "ks_file=s" =>       \$ks_file,
  "cutoff_file=s" =>   \$cutoff_file,
  "out_file:s" =>      \$out_file,
  "scale_factor:f" =>  \$scale_factor,
  "min_transform:i" => \$min_transform,
  "max_transform:i" => \$max_transform,
  "verbose" =>         \$verbose,
  "help:s" =>          \$help
);

my $usage = <<EOS;
  Read a tabular file of gene and Ks data and a file of Ks threshold values per species pair
  and filter the gene-pair file to a given multiple of the Ks threshold values.

  Output rows consist of a gene pair and a score between 10 and 100 (default values), 
  with higher numbers corresponding with lower Ks (and greater similarity):
  transformed_ks = 
      min_transform + (max_transform-min_transform)*((adjusted_cutoff-ks)/adjusted_cutoff);

  Example of the gene-pair format:
    aradu.vigun  aradu.Aradu.000JC  vigun.Vigun05g296900.1  0.8964
    aradu.vigun  aradu.Aradu.001N3  vigun.Vigun10g098600.1  0.5134
    aradu.vigun  aradu.Aradu.002J3  vigun.Vigun07g072400.1  0.5779

  Example of the Ks threshold format:
    aradu.vigun  0.65
    araip.cerca  0.65
    araip.chafa  0.7

  Usage: perl $0 -ks_file FILE -cutoff_file FILE [options] > out_file
  
  -ks_file:       (required) tabular gene-pair file, with Ks values (see format above)
  -cutoff_file:   (required) tabular file with species-pair in first column & Ks threshold in second

  -scale_factor:  factor for scaling the Ks threshold. Range typically 0.5-2.0. Default 1.
                  Useful if the Ks values in the threshold file are the Ks modes, and you
                  wish the cutoff to be a factor of the mode - say, 1.5x(mode)
  -min_transform: Transform Ks output to a new min-max range. This is the min value; default 10
  -max_transform: Transform Ks output to a new min-max range. This is the max value; default 100

  -v|verbose:     report debug-type info to stderr
  -h|help:        display this help message
EOS

die "\n$usage\n" if ($help);
die "\n$usage\nPlease provide ks_file\n" unless (defined($ks_file));
die "\n$usage\nPlease provide cutoff_file\n" unless (defined($cutoff_file));

open(my $CUTFH, "<", $cutoff_file) or die "can't open $cutoff_file: $!";
open(my $GENEPAIRSFH, "<", $ks_file) or die "can't open $ks_file: $!";

if ($verbose){
  print STDERR "Working on Ks file [$ks_file] with cutoff file [$cutoff_file]\n"; 
}

my %ks_cutoffs;
while (<$CUTFH>) {
  chomp;
  next if $_ =~ /^$/;
  my ($sp_pair, $ks_value) = split(/\s+/, $_);
  $ks_cutoffs{$sp_pair} = $ks_value;
  #print "$sp_pair, $ks_value\n";
}

# read Ks file and filter based on coords file 
my %seen_key;
while (<$GENEPAIRSFH>) {
  my $line = $_;
  chomp;
  next if $line =~ /^$/;
  my ($sp_pair_key, $gene1, $gene2, $ks) = split(/\t/, $_);

  # print "$sp_pair_key, $scale_factor, $adjusted_cutoff\n";
  if (defined($ks_cutoffs{$sp_pair_key})) {
    my $adjusted_cutoff = $scale_factor*($ks_cutoffs{$sp_pair_key});
    my $transformed_ks = 
      $min_transform + ($max_transform-$min_transform)*(($adjusted_cutoff-$ks)/$adjusted_cutoff);
    if ( $ks <= $adjusted_cutoff ) {
      print "$gene1\t$gene2\t";
      printf("%.4f\n",$transformed_ks);
      #print "$gene1\t$gene2\t$transformed_ks\n";
      # print "$gene1\t$gene2\t$ks\n";
      # print "YES:\t$ks\t<\t", $scale_factor * ($ks_cutoffs{$sp_pair_key}), "\n" ;
    }
    else {
      # print "NO:\t$ks\n";
    }
  }
  else {
    warn "WARNING: Ks cutoff not defined for species pair: [$sp_pair_key]\n";
  }
}

__END__

Versions
v01 2017-12-24 Start. Is functional.
v01 2018-02-19 Add some documentation and verbose option



