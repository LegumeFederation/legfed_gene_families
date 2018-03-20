#!/usr/bin/env perl

# PROGRAM: root_tree_by_species.pl
# VERSION: see version notes at bottom of file.
# S. Cannon 2017
# see description under Usage

use strict;
use warnings;
use Getopt::Long;
use File::Basename;

####### File IO ########

my $raxml_executable = "/opt/raxml/8.0.26/raxmlHPC-PTHREADS-AVX";
my ($align_in, $output_path, $roots, $help, $verbose, $execute);
my $threads=2;
my $method="PROTGAMMAAUTO";

GetOptions (
  "align_in=s" =>     \$align_in,   # required
  "output_path=s" =>  \$output_path,   # required
  "roots:s" =>        \$roots,   
  "verbose+" =>       \$verbose,
  "execute" =>        \$execute,
  "threads:i" =>      \$threads,
  "help" =>           \$help
);

my $usage = <<EOS;
  Usage: perl $0 [-options]
  
  Calculate a RAxML tree and root it by specified sequence patterns, if those sequences are present.
  Rooting is done by a sequence name or substring in the name, if the name is provided and present; 
  otherwise, the rooting uses the tree-height minimizing method of RAxML.

  For example, if there are outgroup species
    arath.AT1G77750.1
    arath.AT5G14320.1
    orysa.LOC_Os03g49710.1
  and the outgroup strings are supplied as follows:
    -roots "orysa,arath"
  then rooting will be done on orysa and then arath; and if orysa and arath are not 
  monophyletic, the first name in the list will be used for identifying an outgroup
  i.e. orysa.LOC_Os03g49710.1 will be used. If none is found, then 
  the tree will be rooted using the RAxML method "-f I":
    "It roots the tree by rooting it at the branch that best balances the subtree lengths
     (sum over branches in the subtrees) of the left and right subtree."

  Flags and parameters:
   -align_in:    Input alignment (only used for identifying sequence IDs) **
   -output_path: Output directory (FULL, ABSOLUTE PATH) for rooted trees
   -roots:       List of strings found in sequence names to be used as an outgroup,
                    comma-separated, within double quotes, e.g. -roots="orysa,arath"
                    If missing, then the tree will be rooted using the RAxML method "-f I"
   -execute     Execute the command
   -threads     Threads used by RAxML. Default 2
   -verbose     For more stuff to terminal
   -help        For more info
   
   ** = required
EOS

die "\n$usage\n" 
  if ($help or !defined($align_in) or !defined($output_path) );

warn "Note: -execute wasn't set, so the command will be printed but not executed.\n" if (!defined($execute));

# read alignment in
open( my $ALIGN_IN, '<', $align_in ) or die "can't open align_in $align_in: $!";

# tree file name
my $name = fileparse($align_in);

my %seq_IDs;
while (my $line = <$ALIGN_IN>) {
  next unless ($line =~ /^>(\S+) *.*/);
  my $ID=$1;
  #print "$ID\n";
  $seq_IDs{$ID} = $ID;
}

my $raxml_command;
my $seed=int(100*rand());

if (defined($roots)) { # List of roots is provided
  # Put elements of LIST into array in case we need to retain list order
  my @outgrp_strs = split(/,/, $roots);

  my %seen_seq_ID;
  my $set_of_out_IDs;
  for my $root_str (@outgrp_strs) {
    #print "===== ", $root_str, " =====\n";
    for my $seq_ID (%seq_IDs) {
      #print "  $root_str $seq_ID\n";
      if ($seq_ID =~ /$root_str/) {
        if ($seen_seq_ID{$seq_ID}) { next }
        else {
          #print "$seq_ID\n";
          $set_of_out_IDs .= "$seq_ID,";
          $seen_seq_ID{$seq_ID}++;
        }
      }
    }
  }
  $set_of_out_IDs =~ s/,$//;
  
  if (length($set_of_out_IDs) > 0) { # Do outgroup rooting
    $raxml_command = "$raxml_executable -T $threads -m $method -o $set_of_out_IDs -s $align_in -w $output_path -n $name -p $seed";
  }
  else {
    if ($verbose) { print "can't find outgroup" }
    $raxml_command = "$raxml_executable -T $threads -m $method -f I -s $align_in -w $output_path -n $name -p $seed";
  }

} 
else { # No list of roots provided, so do midpoint rooting
  $raxml_command = "$raxml_executable -T $threads -m $method -f I -s $align_in -w $output_path -n $name -p $seed";
}

if ($verbose) {
  print $raxml_command, "\n";
}

if ($execute) {
  system($raxml_command);
}

__END__
VERSIONS
SC = Steven Cannon
v0.01 17-03-28 SC Initial version
v0.02 17-07-20 SC Add warning if -execute isn't set. 
                    Don't print outgroup $seq_ID again if already seen. 
                    Add -threads
                    Implement midpoint rooting, if outgroup list isn't provided or no outgroup is found in tree
v0.03 17-07-20 SC Remove option for taking in a user tree. 


root_tree_by_species.pl -tree 09_trees/59026833_a -align 08_hmmalign_trimmed/59026833_a -out $PWD/09_trees_rt_by_outgrp -roots "ambtr,zeama,orysa,vitvi,solly,arath,prupe" -v -e

