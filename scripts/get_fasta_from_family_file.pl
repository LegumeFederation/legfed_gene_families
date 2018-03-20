#!/usr/bin/env perl

# PROGRAM: get_fasta_from_family_file.pl
# VERSION: see version notes at bottom of file.
# S. Cannon 2017
# see description under Usage

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;

####### File IO ########

my ($input_fas, $family_file, $help, $verbose, $out_dir);

GetOptions (
  "input_fas=s" =>   \$input_fas,   # required
  "family_file=s" => \$family_file,   # required
  "out_dir=s" =>     \$out_dir,   # required
  "verbose+" =>      \$verbose,
  "help" =>          \$help
);

my $usage = <<EOS;
  Usage: perl $0 [-options]
  
  Read a file of families consisting of rows each beginning with a family ID and a list of 
  fasta IDs; and a fasta file. Print those sequences to new fasta files, one per family.
  Format of family file is like this:
    #Fam_ID gene_IDs
    L.v7lGl cicar.Ca_26226  cicar.Ca_27089  lotja.Lj2g3v2291660.1 lotja.Lj0g3v0113769.1
    L.V7ll5 glyma.Glyma.01G021500 glyma.Glyma.09G200100 phavu.Phvul.002G145500  vigun.Viungv11005429m
    L.V7mqj glyma.Glyma.01G098900 glyma.Glyma.02G175600 medtr.Medtr8g079750.1 medtr.Medtr4g062590.1
    L.V7P9b medtr.Medtr3g449470.1 medtr.Medtr6g029430.1 medtr.Medtr6g012600.1 medtr.Medtr7g045640.1
   
   -input_fas:   input fasta file **
   -out_dir:     directory for output (to contain one file per line of family_file) **
   -family_file: file with family ID and list of sequence IDs, described above **
   -verbose      (boolean) for more stuff to terminal; may call multiple times
   -help         for more info
   
   ** = required
EOS

die "\n$usage\n" 
  if ($help or !defined($input_fas) or !defined($out_dir) or !defined($family_file) );

# Read in the sequence using the Bioperl SeqIO;
my $seqio_obj  = Bio::SeqIO->new(-file => $input_fas , '-format' => 'Fasta');
my %seq_hsh;
while ( my $seq_obj = $seqio_obj->next_seq ) {
  my $display_id = $seq_obj->display_id();
  $seq_hsh{$display_id} = $seq_obj->seq();
}

# Read in family file
open( my $FAM_FH, '<', $family_file ) or die "can't open family_file $family_file: $!";

# Traverse family file and write one fasta file per family
while (<$FAM_FH>) {
  chomp;
  next if ( $_ =~ /^$/ or $_ =~ /^#/);

  my ($fam_ID, @seq_IDs) = split(/\s+/, $_);
  if ($verbose) { print "$fam_ID\n" }

  # Open file for writing new fasta out
  open( my $OUT_FH, '>', "$out_dir/$fam_ID" ) or die "can't open out $out_dir/$fam_ID: $!"; 
  
  for my $seq_ID (@seq_IDs){
    if (defined($seq_hsh{$seq_ID})) {
      print $OUT_FH ">$seq_ID\n$seq_hsh{$seq_ID}\n";
    }
    else {
      warn "WARNING: can't find sequence for $seq_ID\n";
    }
  }
}

__END__
VERSIONS

v01 2017-07-02 Initial version
