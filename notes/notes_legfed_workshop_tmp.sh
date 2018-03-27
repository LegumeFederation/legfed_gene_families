
# This file contains some notes for exploring the legfed gene families,
# for the Sata Fe LegFed workshop, 2018-03-27. 
# The file will be removed from this repository after the workshop.

##########
# Gene families on temporary dev URL
  http://fisher.ncgr.org:50005

##########
# Gene family exploration
# Do the work on your ~/Desktop (or elsewhere if you prefer)
  cd ~/Desktop

# Get a couple of files from the Data Store, from public/Gene_families
#   legume.genefam.fam1.M65K.info_fam_composition_list.tsv.gz
#   legume.genefam.fam1.M65K.info_fam_composition_sum.tsv.gz

  git clone https://github.com/LegumeFederation/legfed_gene_families

# Move into the gene families directory and explore a bit
  cd legfed_gene_families
  head stats/*

#####
# Explore the "sequences per family" file
# Check numbers of sequences per family in stats/species_in_55_hmmalign_trim2_legonly.tsv
  cat stats/species_in_55_hmmalign_trim2_legonly.tsv |
    awk '{print NF-1}'

  cat stats/species_in_55_hmmalign_trim2_legonly.tsv |
    awk '{print NF-1}' | histogram -n -s1

  cat stats/species_in_55_hmmalign_trim2_legonly.tsv |
    awk '{print NF-1}' | histogram -n -s1 | histplot -d 25 | head -40

  # more advanced: the same thing, but with a perl one-liner instead of histplot:
  cat stats/species_in_55_hmmalign_trim2_legonly.tsv |
    awk '{print NF-1}' | histogram -n -s1 | 
    perl -lane 'print "$F[0]\t", "." x ($F[1]/30)' | head -40


#####
# Explore the "composition file
  head stats/composition_55_hmmalign_trim2.tsv
  
  # Find average count for some species
  cat stats/composition_55_hmmalign_trim2.tsv |
    awk 'NR>1 {ct++; sum+=$7} END{print sum, ct, sum/ct}'
  
  # Find average count for some species
  cat stats/composition_55_hmmalign_trim2.tsv |
    awk 'NR>1 {ct++; sum+=$8} END{print sum, ct, sum/ct}'
  
  # Wrap this in a loop
    for i in {7..20}; do
      cat stats/composition_55_hmmalign_trim2.tsv |
        awk -v COL=$i 'NR>1 {ct++; sum+=$COL} END{print head, sum/ct}'
    done > tmp.ave_per_sp
  
  # Get species IDs ...
    head -1 stats/composition_55_hmmalign_trim2.tsv | 
      cut -f7-20 | perl -pe 's/\t/\n/g' > tmp.sp_IDs
  # ... then combine with average counts per species
    paste tmp.sp_IDs tmp.ave_per_sp 
  
  


##########
# Synteny exploration
# Online:
#   https://legumeinfo.org/genomes/gbrowse/Va3.0
# From Data Store, get GFF files for Vigna angularis
#   https://legumefederation.org/data/public/Vigna_angularis/
#                                      Gyeongwon.gnm3.syn1.wCFF
  tar -xzf Gyeongwon.gnm3.syn1.wCFF.tar
  cd Gyeongwon.gnm3.syn1.wCFF
  gunzip *gz

# Check Ks peaks for synteny blocks (median values for gene pairs per block)
  grep -v "#" vigan.Gyeongwon.v3.x.glyma.Wm82.a2.gff | 
    cut -f4 -d'=' | histogram -n -s0.05 | histplot -d 10

# Or for all files
  for file in *gff; do 
    echo $file; 
    grep -v "^#" $file | cut -f4 -d'=' | histogram -z -n -s0.05 | histplot -d 10; 
    echo
  done



##########
Command-line tips for working with data

I make a project directory, and within that directory, I make subdirectories only one level deep
(no further nesting of subdirectories). Most subdirectories are named with a number and a short name,
e.g. 01_fasta . The number indicates the place of the data in the analysis project -- larger numbers
being later in the process. I also have some standard directories, e.g. notes, blastdb, blastout. Advantages:
  - Most commands can be executed from the level of the project directory.
  - Numbers help to order directories and give clues about the chronology of the analysis.
  - Directories can be accessed quickly, by typing the number and "tabbing-out" the rest of the name.

I maintain detailed notes for every project, and keep them always in notes/ in the project directory*.
  - I give the notes file the suffix ".sh" and use BASH syntax:
  - # Comments are noted with leading hash symbol
  - Executable commands are entered without a hash symbol. If the file is constructed carefull, 
    then the whole analysis can be re-executed.

I maintain small, project-specific shell scripts in commands/ in the project directory
  - When a project involves numerous files at a given step (e.g. gene families, or multiple species,
    or analyses that need to be manually parallelized), I write short shell scripts - often in order
    to run multiple jobs in the background, on a cluster. I typically put these "job scripts" in a
    "commands" directory within the project directory. These scripts typically have the same form 
    and patterns:
      #!/usr/bin/env ksh93     # A shell that allows simple job control, using JOBMAX

      set -o errexit           # The script should fail on an error
      set -o nounset           # An unset variable is an error

      INDIR=$1                 # Take directory names and parameters in as arguments to the script
      OUTDIR=$2
      JOBMAX=$3

      for path in $INDIR/*; do     # directory+file = "path"
        file=`basename $path`;     # file is the path without the directory
        muscle -in $path -out $OUTDIR/$file -quiet &   # jobs are run "in the background"
      done
      wait                     # "wait" is useful if I have one batch job followed by another batch

  - Calling the control script in the following way helps to self-document the parameters:
    # Do alignments
      mkdir 13_align

      INDIR=12_family_fasta
      OUTDIR=13_align
      JOBMAX=20

      nohup commands/do_align.ksh $INDIR $OUTDIR $JOBMAX &

I maintain scripts that are more substantial or likely to be reused in my $HOME/bin/
... and on github (am trying to be better about that :-)

For data exploration, cleanup, transformation, etc., the following tools are essential:
  - perl one-liners. Useful perl one-liner flags:
      perl -ne    # read every line of input (STDIN or file)
      perl -pe    # read and print every line
      perl -pi -e # modify a file in-place (useful but potentially dangerous)
      perl -le    # execute, and add a line return - e.g. perl -le "print 2+3"
      perl -lane  # read every line, auto-splitting input into array @F

  - awk, for simple filtering and transformations - especially, of tabular data
  - cut | sort | uniq -c
  - join, paste
  

Some of the more common unix commands I used during the gene family constructions 
(with counts of the commands in my notes)

     ll    # aliased to  ls -l
     ls
     rmdir
     mkdir
     rm
     cp
     mv
     cd
     chmod
     export

   1 xargs
   4 paste
   4 sed
   4 sh
   5 comm
   7 head
   7 join
   7 wc
  11 uniq
  19 cut
  20 echo
  31 grep
  49 cat
  49 for ... do ... done
  49 sort

  37 awk
  41 perl


