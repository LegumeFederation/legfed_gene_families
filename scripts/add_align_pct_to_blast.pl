#!/usr/bin/env perl

# Program: add_align_pct_to_blast.pl
# Program description: see usage message.
# S. Cannon 2017

use strict;
use warnings;
use strict;
use Getopt::Long;

my ($blast_file, $qlens, $slens);
my $cov_pct=0;
my $id_pct=0;
my $help;

GetOptions (
  "blast_file=s" => \$blast_file,
  "qlens=s" =>      \$qlens,
  "slens=s" =>      \$slens,
  "cov_pct:i" =>    \$cov_pct,
  "id_pct:i" =>     \$id_pct,
  "help:s" =>       \$help
);

my $usage = <<EOS;
  Reads a tabular BLAST file and tables of query and subject sequence lengths.
  Adds four fields to BLAST output file: query len, subj len, pct of query aligned, pct of subj aligned. 

  Usage: perl $0 -blast_file FILE -qlens FILE -slens FILE [options]  > out_file
  
  -blast_file (required) tabular blast output
  -qlens      (required) file with lengths of query sequences
  -slens      (required) file with lengths of subject sequences
  -cov_pct    coverage threshold in percent (align/qry & align/subj). Default 0; max 100
  -id_pct     identity threshold in percent (BLAST col 3). Default 0; max 100

  -h|help:    display this help message
EOS

die "\n$usage\n" if ($help);
die "\n$usage\nPlease provide blast_file\n" unless (defined($blast_file));
die "\n$usage\nPlease provide qlens\n" unless (defined($qlens));
die "\n$usage\nPlease provide slens\n" unless (defined($slens));

open(my $BLASTFH, "<", $blast_file) or die "can't open $blast_file: $!";
open(my $QLENFH, "<", $qlens) or die "can't open $qlens: $!";
open(my $SLENFH, "<", $slens) or die "can't open $slens: $!";

my %QlensH;
while (<$QLENFH>) {
  chomp;
  my @bits = split(/\s+/, $_);
  my $key = $bits[0];
  $QlensH{$key} = $bits[1];
  #print "$key, $QlensH{$key}\n";
}
close $QLENFH;

my %SlensH;
while (<$SLENFH>) {
  chomp;
  my @bits = split(/\s+/, $_);
  my $key = $bits[0];
  $SlensH{$key} = $bits[1];
  #print "$key, $SlensH{$key}\n";
}
close $SLENFH;

# read BLAST file and add length info
while (<$BLASTFH>) {
  chomp;
  my @bits = split(/\t/, $_);
  my ($qry,$subj,$blast_pct_id,$aln_len) = ($bits[0],$bits[1],$bits[2],$bits[3]);
  my ($Qlen, $Slen) = ($QlensH{$qry}, $SlensH{$subj});
  my ($Qcov_pct, $Scov_pct) = (100*$aln_len/$QlensH{$qry}, 100*$aln_len/$SlensH{$subj});
  if ($Qcov_pct>=$cov_pct && $Scov_pct>=$cov_pct && $blast_pct_id>=$id_pct) {
    print join("\t", @bits);
    printf ("\t%d\t%d\t%.2f\t%.2f\n", $Qlen, $Slen, $Qcov_pct, $Scov_pct);
  }
}
close $BLASTFH;

__END__

Versions
v01 2017-12-19



