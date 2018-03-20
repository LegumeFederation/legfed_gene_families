# Objective: identify ortholog sets between species pairs, and WGD-derived paralogs within species

# Work described in these notes was carried out on a 64-processor linux machine, generally using 
# bash commands except in a set of small process-control scripts running under KornShell (ksh93),
# to take advantage of the JOBMAX variable, which defines the maximum number running background 
# jobs that can run at a time. For running comparable computations in an environment with 
# job schedulers, those scripts will need to be modified accordingly.
# Notes below have been simplified somewhat. They should capture the main tasks and processes involved
# in the tree calculation for the LegumeFederation 2018 legume gene families, but may differ in some
# particulars: side-tests for making parameter selections, etc.

# Requirements/dependencies:
#  ksh93  or other job management system (see comments above)
#  cluster environment for parallelization of expensive computations (search, alignment, tree calculation)
#  bioperl    sequence and alignment manipulation
#  blast      seqence search
#  muscle     alignments
#  fasttree   tree calculation
#  RAxML      tree calculation
#  hmmer package  sequence search and alignment
#  synonymous_calc.py  calculating synonymous sequence changes - by Haibo Tang
#  mcl        Markov clustering


# Raw ingredients for the particular legume family calculations
# Species lists (composed of the three first letters of the genus and two of the species, 
# e.g. glyma = Glycine max
  # Selected dicot outgroups to the legumes
    arath cucsa prupe solly vitvi
  # Legumes
    aradu araip cajca cicar glyma lotja lupan medtr phavu tripr vigan vigra vigun arahy

####################
# The general outline of the procedure is ...


# Get copies of proteome data sets. (Get copies rather than linking, because we will focus on just a subset).
# NOTE: changes were made in many of these on 2017-05-09 for correspondence with CDS files
  mkdir 01_proteomes
  pushd /data/gene_families/01_proteomes
  cp aradu.pep.fa  arath.pep.fa  glyma.pep.fa  lupan.pep.fa phavu.pep.fa  solly.pep.fa  vigan.pep.fa  vitvi.pep.fa \
    araip.pep.fa  cicar.pep.fa  lotja.pep.fa  medtr.pep.fa  prupe.pep.fa  tripr.pep.fa  vigra.pep.fa \
    /scratch/scannon/legume_genefams/01_proteomes
  popd
  
# Do all-by-all BLAST searches, for all proteomes.
# This step takes a while, and needs to be spread across numerous processors, under job management.

  # Make blast databases
    for path in 01_proteomes/*fa; do 
      file=`basename $path .pep.fa`; 
      makeblastdb -in $path -dbtype prot -hash_index -parse_seqids -out blastdb/$file &
    done

  # Blast each species pair. See code for doing the "upper triangle" or lower triangle of comparisons.
    FASTADIR1=01_proteomes
    FASTADIR2=01_proteomes
    JOBMAX=30

    nohup ./commands/do_blast_all.ksh &

# Get top two hits per query - except for self-comparisons, where we only want one top (non-self) match
# Also (new step 2017-12), filter matches by percent query alignment and percent subject alignment
  
  mkdir 02_pep_len
  for path in 01_proteomes/*; do 
    file=`basename $path .pep.fa`
    seqlen.awk $path > 02_pep_len/$file
  done &

  mkdir blastout_cov

  # Call script with percent coverage = 50 and percent identity 60
  # NOTE: Filenames in the BLAST results are assumed to have the format medtr.x.lotja.blp
  # where the query and subject species names each have five characters and are separated by .x.

    LEN_DIR=02_pep_len
    COV_DIR=blastout_cov
    COV_PCT=50
    ID_PCT=60
    JOBMAX=20

    nohup ../commands/do_blast_coverage_filter.ksh $LEN_DIR $COV_DIR $COV_PCT $ID_PCT $JOBMAX &

# Pick top two hits per query - except top non-self single hit for self-comparisons

  mkdir blastout_top2

    for path in blastout_cov/*; do 
      file=`basename $path`
      echo $file
      top_2_lines_no_self.awk $path > blastout_top2/$file.top 
    done &

  # Prune the self-hit files to one hit per query 
     
    for SP in aradu araip arath cajca cicar cucsa glyma lotja lupan medtr \
              phavu prupe solly tripr vigan vigra vigun vitvi ; do
      file="$SP.x.$SP"
      echo $file
      top_line_no_self.awk blastout_cov/$file > blastout_top2/$file.top 
    done &


# For all comparisons with glyma as subject, pick top FOUR non-self hits 
# (non-self only being relevant for glyma x glyma)
    # three lines for glyma.x.glyma
      top_3_lines_no_self.awk MAX=3 blastout_cov/glyma.x.glyma > blastout_top2/glyma.x.glyma.top &
    
    # For all comparisons with glyma as subject, pick top FOUR non-self hits 
    # (non-self only being relevant for glyma x glyma)
    for SP in aradu araip arath cajca cicar cucsa lotja lupan medtr \
              phavu prupe solly tripr vigan vigra vigun vitvi ; do
      file="$SP.x.glyma"
      echo $file
      top_4_lines.awk blastout_cov/$file > blastout_top2/$file.top &
    done

# In preparation for Ks calculations, concatenate CDS into one file, and similarly for peptide sequence
# CDS and peptides come from /data/gene_families/02_CDS/ and /data/gene_families/01_proteomes/
  cat 01_proteomes/*fa > tmp.all_peptides.fa
  cat 01_CDS/*fa > tmp.all_CDS.fa

# In preparation for calculating Ks, process top-two blast hits, reporting unique query-target list, where query is always lexically before target
  cat blastout_top2/* | awk '$1<=$2 {print $1 "\t" $2} $2<$1 {print $2 "\t" $1}' | 
    sort -u > lis.top_two_blast_hits &


# Find for what pairs Ks has already been calculated. 
# Skip amborella, rice, sorghum for now, because the eudicot outgroups should suffice for outgroup pruning in the legumes.
  nohup cat 05_ks/*ks | perl -pe 's/[;,]/\t/g' | awk 'NR>1 {print $1 "\t" $2}' | sort -u > lis.ks_done &

  comm -13 lis.ks_done lis.top_two_blast_hits | 
    grep -v ambtr | grep -v sorbi | grep -v orysa > lis.ks_to_be_calculated.L &
  
# Get CDS and peptide files for those gene lists. 
  mkdir 03_pair_lists
  perl -pe 's/(\w\w\w\w\w)\.(\S+)\t(\w\w\w\w\w)\.(\S+)$/$1.$3\t$1.$2\t$3.$4/' lis.ks_to_be_calculated.L | 
    sort > lis.ks_to_be_calc_for_splitting &
    #NOTE: lis.ks_to_be_calc_for_splitting should be a three-column file, with the first being e.g. glyma.medtr

  split_table_to_files.pl lis.ks_to_be_calc_for_splitting 03_pair_lists 1 &

# For synonymous_calc.py, print query and subject in one column: Q S Q S Q S ...
  perl -pi -e 's/^\S+\t(\S+)\t(\S+)/$1\n$2/' 03_pair_lists/* &

# Get fasta sequence, given the lists above. Do this in a script for parallelization
  mkdir 04_pair_peptide_fasta 04_pair_CDS_fasta 

  nohup commands/get_fasta.ksh tmp.all_CDS.fa      03_pair_lists  04_pair_CDS_fasta &
  nohup commands/get_fasta.ksh tmp.all_peptides.fa 03_pair_lists  04_pair_peptide_fasta &

# Calculate ks
  mkdir 05_ks
  
  PEPDIR=04_pair_peptide_fasta
  CDSDIR=04_pair_CDS_fasta
  WORKDIR=05_ks
  MAINDIR=/scratch/scannon/lgf4/LEFTOVERS
  JOBMAX=60

  nohup commands/syn_calc.ksh $PEPDIR $CDSDIR $WORKDIR $MAINDIR $JOBMAX &

#####
# Planning for next steps: we want to get Ks results for all in lis.top_two_blast_hits
# The file lis.top_two_blast_hits has two columns. To find genes in both columns, double the file,
# to ensure that both genes occur in the first column; then sort and join with the list of long-branch genes.
  LANG=en_EN # so that sorting is consistent and suitable for the join command
  awk '{print $1 "\t" $2; print $2 "\t" $1}' lis.top_two_blast_hits | sort -k1,1 -k2,2 > lis.top_two_blast_hits_x2 &

  perl -pi -e 's/$/\tkill/' lis.long_branches_ge1
  LANG=en_EN sort -o lis.long_branches_ge1 lis.long_branches_ge1

  join -a 2 lis.long_branches_ge1 lis.top_two_blast_hits_x2 |
    grep -v kill | awk '$1 <= $2' | sort -u > lis.top_two_blast_hits_no_long &

  # Intermediate product: lis.top_two_blast_hits_no_long is the set of pairs for which we want Ks values


# Join the (doubled, parsed) Ks results with the blast list. The blast list adds the requirements of:
# percent coverage (50) and percent identity (60). If we had calculated Ks from that list in the first place,
# we wouldn't need to do the following join step; but since some of the Ks values were calculated earlier,
# we do this additional join/filter here.

  # Fields are:
    #   1     2     3      4      5      6
    #   name  name  dS-yn  dN-yn  dS-ng  dN-ng
    # Use 5: dS-Nei-Gojobori

    cat 05_ks_combined/* | 
      perl -pe 's/[;,]/\t/g; s/^(\w\w\w\w\w)\.(\S+)\t(\w\w\w\w\w)\.(\S+)\t/$1\t$3\t$1.$2\t$3.$4\t/' | 
      awk -v OFS="\t" '$1 <= $2 && $7>=0 && $7<2.5 {print $1 "." $2, $3, $4, $7 }
                       $2 <  $1 && $7>=0 && $7<2.5 {print $2 "." $1, $4, $3, $7 }' |
       LC_COLLATE=C sort -u -k1,1 -k2,2 > tab.all_Ks_before_splitting_by_sp_pair &

    split_table_to_files.pl tab.all_Ks_before_splitting_by_sp_pair 05_ks_combined_uniq 1 &

    # then remove the leading gensp.gensp\t  and join the two gene fields
    perl -pi -e 's/^\S+\t//; s/\t/___/' 05_ks_combined_uniq/* &

    # then sort on the combined gene___gene field, with LC_COLLATE=C
    for file in 05_ks_combined_uniq/*; do LC_COLLATE=C sort -k1,1 -o $file $file; done &

  
# Make the sorted BLAST file with two (joined) names, and Ks file, with two (joined) names and a Ks column
    mkdir 05_blast_uniq
    
    cat lis.top_two_blast_hits_no_long | 
      perl -pe 's/^(\w\w\w\w\w)\.(\S+)\s(\w\w\w\w\w)\.(\S+)/$1\t$3\t$1.$2\t$3.$4\t/' |
      awk -v OFS="\t" '$1 <= $2 {print $1 "." $2, $3, $4 }
                       $2 <  $1 {print $2 "." $1, $4, $3 }' |
       LC_COLLATE=C sort -k1,1 -k2,2 > lis.top_two_blast_hits_no_long_prefixed &

    split_table_to_files.pl lis.top_two_blast_hits_no_long_prefixed  05_blast_uniq 1 &

    # then remove the leading gensp.gensp\t  and join the two gene fields
    perl -pi -e 's/^\S+\t//; s/\t/___/' 05_blast_uniq/* &

    # then sort on the combined gene___gene field, with LC_COLLATE=C
    for file in 05_blast_uniq/*; do LC_COLLATE=C sort -k1,1 -o $file $file; done &
    
# Join the BLAST files and the Ks files - each of which was split into species-pair files above.
    mkdir 06_blast_ks_joined

    for path in 05_ks_combined_uniq/*; do 
      file=`basename $path`
      join 05_blast_uniq/$file 05_ks_combined_uniq/$file | 
        perl -pe 's/^(\w\w\w\w\w)\.(.+)___(\w\w\w\w\w)\.(\S+) (\d.+)/$1.$3\t$1.$2\t$3.$4\t$5/' \
        > 06_blast_ks_joined/$file
    done &


# Make histograms and find modal WGD Ks values
  mkdir histograms

  TAB=$'\t'
  DIR_OUT="histograms"
  for path in 05_ks_combined_uniq/*; do
    LINES=`wc -l $path | cut -f1 -d" "`;
    HIST_DIVISOR=$((10+$LINES/1000));
    FILE=`basename $path`
    cat /dev/null > $DIR_OUT/$FILE
    SPECIES=`head -1 $path | perl -pe 's/^(\w\w\w\w\w)\.\S+\t(\w\w\w\w\w)\.\S+\t.+/$1\t$2/'`
    echo "$FILE$TAB$SPECIES$TAB$LINES" >> $DIR_OUT/$FILE
    cut -f3 $path | histogram -z -n -s0.05 | histplot -d $HIST_DIVISOR >> $DIR_OUT/$FILE
  done  &

  # Get modal value, for legume WGD. This generally isn't the first (speciation) peak, so ignore 
  # values <= 0.4
  TAB=$'\t'
  for path in 05_ks_combined_uniq/*; do
    FILE=`basename $path`
    LINES=`wc -l $path | cut -f1 -d" "`;
    SPECIES=`head -1 $path | perl -pe 's/^(\w\w\w\w\w)\.\S+\t(\w\w\w\w\w)\.\S+\t.+/$1\t$2/'`
    PEAK=`awk '$3>=0.4 && $3<2 {print $3}' $path | histogram -z -n -s0.05 |
            sort -k2nr | head -1 | cut -f1` # end execute into $PEAK
    export SPECIES; export PEAK;
    perl -le 'print $ENV{"SPECIES"}, "\t", $ENV{"PEAK"}'
  done > WGD_peaks2.tab

  # Filter to legume only
    cat WGD_peaks2.tab | awk '!/arath/ && !/ambtr/ && !/solly/ && !/orysa/ && 
                              !/sorbi/ && !/prupe/ && !/vitvi/ && !/cucsa/' \
      > WGD_peaks2_leg_only.tab



# Filter the Ks output by Ks. Filter separately for each species pair, because the expected Ks ranges
# are different for different species pairs. The file with manually-assessed Ks thresholds is
# WGD_peaks2_leg_manual.tab, with format 
  aradu.lotja 0.55

  Example of the gene-pair format:
    aradu.vigun  aradu.Aradu.000JC  vigun.Vigun05g296900.1  0.8964
    aradu.vigun  aradu.Aradu.001N3  vigun.Vigun10g098600.1  0.5134
    aradu.vigun  aradu.Aradu.002J3  vigun.Vigun07g072400.1  0.5779

  # Example for one species-pair: 
    filter_by_ks.pl -ks 06_blast_ks_joined/araip.vigra -cut WGD_peaks2_leg_manual.tab -scale_factor 1.5

  # For all species-pairs:
    for path in 06_blast_ks_joined/* ; do 
      sp_pair=`basename $path`
      filter_by_ks.pl -ks $path -cut WGD_peaks2_leg_manual.tab -scale_factor 1.5
    done > tab.all_pairs_filtered_by_ks  &


  # Cluster the results using markov clustering, testing several inflation values
    nohup mcl tab.all_pairs_filtered_by_ks -I 1.1 -te 24 --abc &
    nohup mcl tab.all_pairs_filtered_by_ks -I 1.2 -te 12 --abc &
    nohup mcl tab.all_pairs_filtered_by_ks -I 1.4 -te 12 --abc &
    nohup mcl tab.all_pairs_filtered_by_ks -I 2.0 -te 12 --abc &
    nohup mcl tab.all_pairs_filtered_by_ks -I 3.0 -te 12 --abc &
    nohup mcl tab.all_pairs_filtered_by_ks -I 4.0 -te 12 --abc &


# Check composition summaries (e.g. counts per species, and numbers of species counts = 1 or 2)
# Also prune small and non-diverse families, using MIN_FAM_SIZE and MIN_SP_COUNT.
# Print to three files: KEEP, CUT, and ct_per_sp

    export MIN_FAM_SIZE=6
    export MIN_SP_COUNT=6
    for I in 11 12 14 20 30 40; do
      export file="out.tab.all_pairs_filtered_by_ks.I$I"
      cat $file | perl -ane '
          BEGIN{
            $MIN_FAM_SIZE=$ENV{"MIN_FAM_SIZE"};
            $MIN_SP_COUNT=$ENV{"MIN_SP_COUNT"};
            $FILE=$ENV{"file"};
            my $fam_num = 0;
            @species = qw(aradu araip cajca cicar glyma lotja lupan 
                          medtr phavu tripr vigan vigra vigun);
            open ($OUTCUT, ">", "$FILE.CUT") or die "cant open out lis.$FILE.CUT, $!\n";
            open ($OUTKEEP, ">", "$FILE.KEEP") or die "cant open out $FILE.KEEP, $!\n";
            open ($OUTCOUNTS, ">", "ct_per_sp.$FILE") or die "cant open out ct_per_sp.$FILE, $!\n";

            print $OUTCOUNTS "num.counts____\t", join("\t", @species),"\n";
          }
          my ($name, $ct_present, $ct_below, $ct_eq, $ct_above, $gene_ct) = (0,0,0,0,0,0);
          my $out_line;
          my %ct_sp;
          for my $sp (@species){
            $ct_sp{$sp} = 0;
            for my $seq (@F){ if ($seq =~ /$sp/) {$ct_sp{$sp}++ } }
            $out_line .= "$ct_sp{$sp}\t";
            if ($ct_sp{$sp}>0){$ct_present++}
            if ($ct_sp{$sp}<2){$ct_below++}
            if ($ct_sp{$sp}==2){$ct_eq++}
            if ($ct_sp{$sp}>2){$ct_above++}
          }
          if (scalar(@F) >= $MIN_FAM_SIZE && $ct_present >= $MIN_SP_COUNT) {
            $gene_ct = scalar(@F);
            $fam_num++;
            my $fam_id_and_stats = "$fam_num.$gene_ct.$ct_present.$ct_below.$ct_eq.$ct_above";
            print $OUTCOUNTS "$fam_id_and_stats\t$out_line\n";
            print $OUTKEEP $fam_id_and_stats, "\t", join("\t",@F),"\n";
          }
          else { # print IDs from small families
            print $OUTCUT join("\n",@F),"\n";
          } ' # end of perl script
      done &

# Get summary stats: histograms of family sizes - summarizing from stats/ct_per_sp.$file above
    base="out.tab.all_pairs_filtered_by_ks.I"
    for I in 11 12 14 20 30 40; do
      grep -v counts "stats/ct_per_sp.$base$I" | cut -f1 | cut -f2 -d'.' | 
        histogram -n -s1 | histplot -d 30 | head -50 > stats/hist_ct_per_sp.$base$I
    done

    for file in stats/hist_ct_per_sp.out.tab.all_pairs_filtered_by_ks*; do 
      echo; echo $file; cat $file; 
    done > stats/summary_hist_ct_per_sp.out.tab.all_pairs_filtered_by_ks_f$MIN_FAM_SIZE.s$MIN_SP_COUNT.txt 


  # Get coarse summary by some meaningful categories - on the pruned families (.KEEP)

    export MIN_FAM_SIZE=6
    export MIN_SP_COUNT=6
    base="out.tab.all_pairs_filtered_by_ks.I"
    for I in 11 12 14 20 30 40; do
      cat $base$I.KEEP |
      perl -lane '$ct=scalar(@F)-1; # subtract family ID, e.g. 16200.6.6.17.0.0
                  if($ct>$max){$max=$ct}
                  if($ct<10){$_lt_10++} 
                  elsif($ct>=10 && $ct<20){$_10_to_20++} 
                  elsif($ct>=20 && $ct<40){$_20_to_40++} 
                  elsif($ct>=40 && $ct<80){$_40_to_80++} 
                  elsif($ct>=80 && $ct<200){$_80_to_200++} 
                  else{$_ge_200++} 
                  $tot_genes += $ct;
                  $tot_fams ++;
                  END{print "_0 to 10\t$_lt_10"; 
                      print "10 to 20\t$_10_to_20"; 
                      print "20 to 40\t$_20_to_40"; 
                      print "40 to 80\t$_40_to_80"; 
                      print "80 to 200\t$_80_to_200"; 
                      print "_ >= 200\t$_ge_200";
                      print "_    max\t$max";
                      print "_  >= 10\t", $_lt_10+$_10_to_20+$_20_to_40+$_40_to_80+$_80_to_200+$_ge_200;
                      print "tot_genes\t$tot_genes";
                      print "tot_fams\t$tot_fams";
                      }
                  ' > stats/ct_coarse.$base$I.KEEP.f$MIN_FAM_SIZE.s$MIN_SP_COUNT.txt
    done
  
    for file in stats/ct_coarse.out.tab.all_pairs_filtered_by_ks.I*.f$MIN_FAM_SIZE.s$MIN_SP_COUNT.txt; do 
      echo; echo $file; cat $file; 
    done > stats/summary.ct_coarse.out.tab.all_pairs_filtered_by_ks.f$MIN_FAM_SIZE.s$MIN_SP_COUNT.txt


##### MAJOR PRODUCT: initial gene families
# Give this file a new name:
  cp stats/out.tab.all_pairs_filtered_by_ks.I11.KEEP legume_gene_families_v1_I11


# Get fasta sequence from the labeled mcl family file. Do this with get_fasta_from_family_file.pl
# in order not to have to re-load the source family file more than once.

  mkdir 12_family_fasta

  nohup get_fasta_from_family_file.pl -in tmp.all_peptides.fa  -out 12_family_fasta  -fam legume_gene_families_v1_I11 &

# Remove terminal X and internal XXs, to avoid pollution of alignments and HMM alignment models 
# (there are quite a few in these sequences). Leave one internal X 
  cd 12_family_fasta
    perl -pi -e 's/^([A-Z][^X]+)X+$/$1/; s/^([A-Z].+[^X]+)XX+([A-Z].+[^X]+)$/$1X$2/g' * &
    mkdir ../12_family_fasta_tmp
    for file in *; do sed '/^$/d' $file > ../12_family_fasta_tmp/$file; done &
    mv ../12_family_fasta_tmp/* .
    cd ..
  # NOTE: This didn't work fully as intended. Some XXs were removed but not all. Did some additional work

# Do alignments
  mkdir 13_align
  
  nohup commands/do_align.ksh 12_family_fasta 13_align &

# separately retrieved and aligned  12_family_fasta/11246.17.16.16.1.0

#####
# Build hmms
  mkdir 14_hmm
        
  nohup commands/do_hmmbuild.ksh 13_align 14_hmm 20 &


# Generate a consensus sequence for each HMM

  mkdir 15_hmmemit
  nohup commands/do_hmmemit.ksh 14_hmm 15_hmmemit 20 &


##########

  # Some stats
    ls 22_family_fasta | wc -l
      15256
    grep -ch '>' 22_family_fasta/* | histogram -n | head -50 | histplot -d 20
    # 6   .......
    # 7   ...................
    # 8   ...........
    # 9   ..........
    # 10  ...........
    # 11  ............
    # 12  .............
    # 13  ................
    # 14  .....................
    # 15  ..............................
    # 16  ..........................................
    # 17  .......................................................
    # 18  .............................................................
    # 19  ..................................................
    # 20  ..........................................
    # 21  ...............................
    # 22  .........................
    # 23  ....................
    # 24  ..................
    # 25  ..................
    # 26  ...............
    # 27  ................
    # 28  ...............
    # 29  ................
    # 30  .................
    # 31  ...............
    # 32  .............
    # 33  ...............
    # 34  ..............
    # 35  ............
    # 36  ...........
    # 37  .......
    # 38  .......
    # 39  .....
    # 40  .....
    # 41  ....
    # 42  ...
    # 43  ...
    # 44  ..
    # 45  ..
    # 46  ..
    # 47  ..
    # 48  ..
    # 49  ..
    # 50  ..
    # 51  
    # 52  .
    # 53  .
    # 54  .
    # 55  .
    # 

# Assign IDs

  COUNT=`ls 22_family_fasta | wc -l`
  echo $COUNT
    15256

  make_random_IDs.pl -count $COUNT -w 6 | sort | perl -pe 's/(\w+)/L_$1/' > lis.fam_names.L



  # Modify temporary names, prefixing them with gene count number, 
  # to help assign random IDs in order by family size
  mkdir 22_family_fasta.tmp
  cd 22_family_fasta
    for file in *; do  
      count=`grep -c '>' $file`; 
      echo $count.$file; 
      cp $file ../22_family_fasta.tmp/$count.$file
    done > ../lis.count_prefix.22_family_fasta &

    cd ..
  
  ls 22_family_fasta.tmp > lis.22_family_fasta.tmp

  mkdir 22_family_fasta.tmp2
  sort -k1nr,1nr -t'.' -o lis.count_prefix.22_family_fasta lis.count_prefix.22_family_fasta
  paste lis.count_prefix.22_family_fasta lis.fam_names.L | 
    perl -pe 's{(\S+)\t(\S+)}{mv 22_family_fasta.tmp/$1 22_family_fasta.tmp2/$2}' > cmd.rename

  sh cmd.rename
  rm 22_family_fasta/*
  mv 22_family_fasta.tmp2/* 22_family_fasta/
  rmdir 22_family_fasta.tmp*
  
#####
# Do alignments
  
  mkdir 23_align
  nohup commands/do_align.ksh 22_family_fasta 23_align &

#####
# Build hmms
        
  mkdir 24_hmm
  nohup commands/do_hmmbuild.ksh 23_align 24_hmm 20 &

# Generate a consensus sequence for each HMM

  mkdir 25_hmmemit
  nohup commands/do_hmmemit.ksh 24_hmm 25_hmmemit 20 &

  # Fix IDs in hmmemit files
    perl -pi -e 's/-consensus//' 25_hmmemit/*

# Add back sequences omitted through initial filtering steps

  # Get list of all legume sequences
    SP_DIR=01_proteomes
    OUT_DIR=tmp_all_species
    mkdir $OUT_DIR

    for SP in aradu araip cajca cicar glyma lotja lupan medtr phavu tripr vigan vigra vigun ; do
      grep '>' $SP_DIR/$SP.pep.fa | cut -f1 -d' ' | sed 's/>//' | sort > $OUT_DIR/$SP &
    done &
  
    cat $OUT_DIR/* > lis.genes_in_legume_proteomes
    
    rm -rf $OUT_DIR

  # Get list of genes in families
    cat 22_family_fasta/* | grep '>' | cut -f1 -d' ' | sed 's/>//' | sort > lis.genes_in_22_family_fasta &

  # Get list of genes not in families
    comm -23 lis.genes_in_legume_proteomes lis.genes_in_22_family_fasta > lis.genes_not_in_fams

    wc -l lis.g*
      361051 lis.genes_in_22_family_fasta
      608925 lis.genes_in_legume_proteomes
      247874 lis.genes_not_in_fams

  # To parallelize, split lis.genes_in_22_family_fasta by species, and get those sequences
    mkdir 22_lis_not_in_fams 22_fas_not_in_fams
    for SP in aradu araip cajca cicar glyma lotja lupan medtr phavu tripr vigan vigra vigun ; do
      grep $SP lis.genes_not_in_fams > 22_lis_not_in_fams/$SP 

      get_fasta_subset.pl -in tmp.all_peptides.fa -lis 22_lis_not_in_fams/$SP -out 22_fas_not_in_fams/$SP 
    done &
    
### Search the leftover genes using hmmscan

  # press hmm database
    cat 24_hmm/* > 24_hmm_ALL
    hmmpress 24_hmm_ALL

  # Search hmm database

    mkdir 27_leftovers_hmmscan
    for SP in aradu araip cajca cicar glyma lotja lupan medtr phavu tripr vigan vigra vigun ; do
        nohup hmmscan --noali -E 1e-4 --cpu 2 --tblout 27_leftovers_hmmscan/$SP.x.24_hmm_ALL.pftbl \
          24_hmm_ALL 22_fas_not_in_fams/$SP &
    done


  # Also compare sequences in families against their HMM, to get a measure of expected similarity
      
    mkdir 26_fam_vs_hmm
    nohup commands/do_hmmsearch_of_fam_vs_hmm.ksh \
      24_hmm 22_family_fasta 26_fam_vs_hmm 26_fam_vs_hmm_humanout &

# Check families against families to see what can be merged
  cat 25_hmmemit/* > 25_hmmemit.fa
  nohup hmmsearch -o /dev/null --tblout 25_hmmemit.x.24_hmm_ALL.pftbl --noali --cpu 12 24_hmm_ALL 25_hmmemit.fa &

# script to calculate median and scaled median-absolute-deviation
# called thus: 
  # awk '$1!~/^#/ {print $6}' L_BB0LK7 | ~/bin/median_abs_dev.pl 
  
  cat /dev/null > tab.26_fam_vs_hmm.MAD
  for path in 26_fam_vs_hmm/* ; do
    file=`basename $path`
    cat $path | awk '$1!~/^#/ {print $6}' | median_abs_dev.pl -file $file >> tab.26_fam_vs_hmm.MAD
  done &

# Get top hmmscan line for hmmemit sequences vs. HMMs
  cat 25_hmmemit.x.24_hmm_ALL.pftbl | top_hmmscan_line_no_self.awk \
    > 25_hmmemit.x.24_hmm_ALL.pftbl.top

# Generate table of candidate merges and associated statistics:
  filter_by_median_abs_dev.pl -hmm 25_hmmemit.x.24_hmm_ALL.pftbl.top \
                              -med tab.26_fam_vs_hmm.MAD -floor 80 | 
                awk '$9~/MATCHES|verdict/' > 25_hmmemit.x.24_hmm_ALL.candidate_merges
    # This generates 550 family pairs, which would reduce to 249 families if using 
    # mcl - -I 2 --abc   (237 two-member familes; 269 at two-or-three-member families)
    # Based on some spreadhsheet work, merge 53 small-ish familes (with <= 25 members) into 24.

# Redo some HMMs and consensus sequences
  ll 22_family_fasta | awk '$7==2 {print $9}' | xargs -I{} cp 22_family_fasta/{} 22_family_fasta_merge/
  mkdir 22_family_fasta_merge 23_align_merge 24_hmm_merge 25_hmmemit_merge
  nohup commands/do_align.ksh 22_family_fasta_merge 23_align_merge &
  nohup commands/do_hmmbuild.ksh 23_align_merge 24_hmm_merge 20 &
  nohup commands/do_hmmemit.ksh 24_hmm_merge 25_hmmemit_merge 20 &

  nohup commands/do_hmmsearch_of_fam_vs_hmm.ksh 24_hmm_merge 22_family_fasta_merge 26_fam_vs_hmm_merge &


  cat 24_hmm_merge/* > 24_hmm_merge_all
  hmmpress 24_hmm_merge_all
  mkdir 27_leftovers_hmmscan_merge

  HMMOUTPUT=27_leftovers_hmmscan_merge/$SP.x.24_hmm_ALL.pftbl

    for SP in aradu araip cajca cicar glyma lotja lupan medtr phavu tripr vigan vigra vigun ; do
      nohup hmmscan --noali -E 1e-4 --cpu 1 --tblout $HMMOUTPUT 24_hmm_merge_all 22_fas_not_in_fams/$SP &
    done


  filter_by_median_abs_dev.pl -hmm 27_leftovers_hmmscan_top.pftbl.top \
                              -med tab.26_fam_vs_hmm.MAD -floor 75 -no_MAD |
                awk '$9~/MATCHES|verdict/' > 27_leftovers_hmmscan_top_MATCHES_fl75noMAD
    # Stats:
        wc -l 27_leftovers_hmmscan_top_MATCHES_fl75noMAD
          21512  # <== use this


# Get top HMM matches for each gene to its gene family. 
# NOTE: I should have done hmmscan rather than hmmsearch. Switch the query and subject fields in the output:
  cat 26_fam_vs_hmm/* | perl -pi -e 's/ +/\t/g' | 
    awk -v OFS="\t" '$1 !~ /^#/ {print $3,$2,$1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19}' \
    > 26_fam_vs_hmm.pftbl.top
  

# Filter by median and ajdusted median absolute deviation (MAD)
  filter_by_median_abs_dev.pl -hmm 26_fam_vs_hmm.pftbl.top \
                              -med tab.26_fam_vs_hmm.MAD -floor 75 -no_MAD |
                awk '$9~/MATCHES|verdict/' > 26_fam_vs_hmm_MATCHES_fl75noMAD

  filter_by_median_abs_dev.pl -hmm 26_fam_vs_hmm.pftbl.top \
                              -med tab.26_fam_vs_hmm.MAD -floor 75 -no_MAD |
                awk '$9~/doesnt|verdict/' > 26_fam_vs_hmm_doesnt_fl75noMAD


  filter_by_median_abs_dev.pl -hmm 26_fam_vs_hmm.pftbl.top \
                              -med tab.26_fam_vs_hmm.MAD -floor 65 -no_MAD |
                awk '$9~/MATCHES|verdict/' > 26_fam_vs_hmm_MATCHES_fl65noMAD

  filter_by_median_abs_dev.pl -hmm 26_fam_vs_hmm.pftbl.top \
                              -med tab.26_fam_vs_hmm.MAD -floor 65 -no_MAD |
                awk '$9~/doesnt|verdict/' > 26_fam_vs_hmm_doesnt_fl65noMAD

  # Stats:
    # using MAD as cutoff; floor 75
      perl -le 'print 280603/(280603+80392)'
      # proportion matching: 0.7773
    # using floor as cutoff (-no_MAD); floor 75
      perl -le 'print 325930/(325930+35066)'
      # proportion matching: 0.9029 
    # using floor as cutoff (-no_MAD); floor 65
      perl -le 'print 336886/(336886+24111)'
      # proportion matching: 0.9332  <== use this

# Combine the "matching" genes from the fam_vs_hmm search, and the matching leftovers
  cat 26_fam_vs_hmm_MATCHES_fl65noMAD 27_leftovers_hmmscan_top_MATCHES_fl75noMAD | 
    awk '$1!~/^#/ {print $1 "\t" $2}' | sort -k1,1 -k2,2 > tab.28_all_MATCHES_via_hmmscan_median

cut -f1 tab.28_all_MATCHES_via_hmmscan_median | uniq -c | awk '{print $1}' | histogram -n -s1 | head -50| histplot -d 20
bin abcdefghijKLMNOPQRSTabcdefghijKLMNOPQRSTabcdefghijKLMNOPQRSTabcdefghijKLMNOPQRST
4   
5   .
6   ......
7   .............
8   ........
9   .......
10  ........
11  ..........
12  .............
13  .................
14  ...........................
15  .......................................
16  ......................................................
17  .................................................................
18  ..............................................................
19  ................................................
20  .................................
21  ..........................
22  .....................
23  ...................
24  .................
25  ................
26  ................
27  ................
28  ................
29  .................
30  ................
31  .................
32  ...............
33  ...............
34  .............
35  ...........
36  .........
37  ......
38  .....
39  ....
40  ....
41  ....
42  ..
43  ..
44  ...
45  ..
46  ..
47  .
48  ...
49  ..
50  ..
51  .
52  .
53  .


# Combine the core sequence set and the unplaced sequences and the top outgroup matches.
# To do this, first reshape output into a family file with format linke that generated by mcl:
# each line begins with a family ID, and the gene list follows, space-separated.
# then use get_fasta_from_family_file.pl to retrieve sequences for each family.

  cat ~/bin/hash_to_rows_by_1st_col.awk
    BEGIN{ORS=""; OFS="\t"}
    NR==1 { print $1 "\t" $2; count = 1; prev = $1 }
    NR>1 && $1 == prev { print "\t" $2; count++ }
    NR>1 && $1 != prev { print "\n" $1 "\t" $2; count = 1; prev = $1 }
    END{print "\n"}

  hash_to_rows_by_1st_col.awk tab.28_all_MATCHES_via_hmmscan_median \
      > 28_all_MATCHES_via_hmmscan_median.fam

# Get fasta sequences for combined (core + leftover) families
  
  mkdir 32_family_fasta
  nohup get_fasta_from_family_file.pl -in tmp.all_peptides.fa -fam 28_all_MATCHES_via_hmmscan_median.fam -out 32_family_fasta &

# do alignments
  mkdir 33_align
  nohup commands/do_align.ksh 32_family_fasta 33_align &

# build HMMs
  mkdir 34_hmm
  nohup commands/do_hmmbuild.ksh 33_align 34_hmm 20 &

  #make into pressed hmm database
  cat 34_hmm/* > 34_hmm_ALL
  hmmpress 34_hmm_ALL


# Add outgroup sequences: one sequence per eudicot outgroup, to help with rooting.
# Pick reciprocal top hits per species and family. 
      
  mkdir 35_outgrp_vs_hmm
 
  for sp in arath cucsa prupe solly vitvi; do
    nohup hmmscan -o /dev/null --cpu 4 -E 1e-4 \
      --tblout 35_outgrp_vs_hmm/$sp.x.34_hmm_ALL.tblout 34_hmm_ALL 01_proteomes/$sp.pep.fa &
  done &

# Find top outgroup match per species. 
# For each outgroup species, first find each gene's best match to a family; 
# then for each family, find the best-matching gene from that species.

    mkdir 36_top_outgroup_seqs
    INDIR=35_outgrp_vs_hmm
    OUTDIR=36_top_outgroup_seqs

    for species in prupe arath vitvi cucsa solly ; do 
      top_hmmscan_line.awk $INDIR/$species.x.34_hmm_ALL.tblout | 
        perl -pe 's/ +/\t/g' | awk -v OFS="\t" '$1!~/^#/ {print $1, $2, $3, $4, $5, $6, $7}' | 
        sort -k1,1 -k6nr,6nr | top_line.awk > $OUTDIR/$species
    done


    INDIR=36_top_outgroup_seqs
    OUTDIR=36_top_outgroup_seqs
    for species in prupe arath vitvi cucsa solly ; do 
      filter_by_median_abs_dev.pl -hmm $OUTDIR/$species \
          -med tab.26_fam_vs_hmm.MAD -floor 50 -no_MAD > $OUTDIR/$species.fl50noMAD
    done
          #-med tab.26_fam_vs_hmm.MAD -floor 65 -no_MAD > $OUTDIR/$species.fl65noMAD
    
    # check stats (proportion of top-matches retained)
      cd 36_top_outgroup_seqs
      for species in prupe arath vitvi cucsa solly; do 
        file=$species.fl65noMAD; lines=`awk 'END{print NR}' $file`; 
        file=$species.fl50noMAD; lines=`awk 'END{print NR}' $file`; 
        #matches=`awk '$9~/MATCHES/ {ct++} END{print ct}' $file`; 
        echo $species $matches $lines $((100*$matches/$lines)); 
      done
      # at 65% of median:
        # prupe  9461 10857 89
        # arath  8383 11496 74
        # vitvi  8756 11327 78
        # cucsa 10537 13041 80
        # solly  8119 10946 76
      # Proportion of families with at least one outgroup sequence:
        cut -f1 *50* | grep -v "#" | sort -u | wc -l
          14723

      # at 50% of median:      <== use this
        # prupe 10102 10857 95
        # arath  9967 11496 87
        # vitvi  9938 11327 89
        # cucsa 11815 13041 90
        # solly  9418 10946 88
      # Proportion of families with at least one outgroup sequence:
        cut -f1 *65* | grep -v "#" | sort -u | wc -l
          14723

     cd /scratch/scannon/lgf4/
    

  # Outgroup seqs to add to the legume fams. Report in two columns: family gene
    awk '$9~/MATCHES/ {print $1 "\t" $2}' 36_top_outgroup_seqs/*.fl50noMAD | 
      sort > hsh.outgroups_to_add_to_legfams
   
  # For current legume families, get list of all sequence IDs in each family. 
  # Report in two columns: family gene
    export FASDIR=32_family_fasta
    grep '>' $FASDIR/* | perl -pe '$DIR=$ENV{"FASDIR"}; s/$DIR\/([^:]+):>(.+)/$1\t$2/' \
      > hsh.core_legfam_IDs_in_legfams

  # Combine the core sequence set and the unplaced sequences and the top outgroup matches.
  # To do this, first reshape output into a family file with format linke that generated by mcl:
  # each line begins with a family ID, and the gene list follows, space-separated.
  # then use get_fasta_from_family_file.pl to retrieve sequences for each family.

    cat hsh.core_legfam_IDs_in_legfams hsh.outgroups_to_add_to_legfams |
      sort | hash_to_rows_by_1st_col.awk  > fam_lists.all_core_plus_out
       # Major product fam_lists.all_core_plus_out ^^ 


# Get fasta sequences for combined (core + leftover + outgroup) families
  
  mkdir 42_family_fasta
  get_fasta_from_family_file.pl -in tmp.all_peptides.fa \
    -fam fam_lists.all_core_plus_out -out 42_family_fasta &

# Align each family to HMM
    
  mkdir 43_hmmalign
  nohup commands/hmm_realign.ksh 42_family_fasta 34_hmm 43_hmmalign 20 &


# Remove non-match characters from A2M file

  mkdir 43_hmmalign_trim1
  nohup commands/hmm_trim1.ksh 43_hmmalign 43_hmmalign_trim1 20 &


# Remove columns, then cut sequences spanning less than 20% of the alignment

  mkdir 43_hmmalign_trim2 43_hmmalign_trim2_log
  nohup commands/hmm_trim2.ksh 43_hmmalign_trim1 43_hmmalign_trim2 43_hmmalign_trim2_log 3 20 20 &
                     # hmmtrimdir1   hmmtrimdir2  logdir  min_depth   min_pct_aligned   JOBMAX

  # An alignment to check for trimming effects: 43_hmmalign/L_3LTHV5
    # Manual work on
    # then run separately:
      filter_align.pl -in 43_hmmalign_trim1/L_3LTHV5 -out 43_hmmalign_trim2/L_3LTHV5 \
        -log 43_hmmalign_trim2_log/L_3LTHV5 -depth 3 -min 20


#####
# Build hmms again, on the trimmed files - mostly to generate better consensus sequences
          
    mkdir 44_hmm
    nohup commands/do_hmmbuild.ksh 43_hmmalign_trim2 44_hmm 20 &


# Generate a consensus sequence for each HMM

  mkdir 45_hmmemit
  nohup commands/do_hmmemit.ksh 44_hmm 45_hmmemit 20 &


#####
# Make trees. For the largest of them (630 of L.0*), use FastTree. For the rest, use RAxML. 
# Create some temp directories. Combine files later.

    mkdir 46_trees 46_trees_RAxML 46_trees_FT
    mkdir 43_hmmalign_trim2_big
    mv 43_hmmalign_trim2/L_0* 43_hmmalign_trim2_big/
    
    chmod u+x commands/do_fasttree.ksh
    nohup commands/do_fasttree.ksh 43_hmmalign_trim2_big 46_trees_FT 20 &


    nohup commands/do_fasttree.ksh 43_hmmalign_trim2 46_trees_FT 20 &


# For remaining families, calculate trees using RAxML with outgroup rooting
# Tree for evaluating: L.NHRL8
  raxml -T 2 -m PROTGAMMAAUTO \
    -o "solly.Solyc02g061940.2.1,vitvi.GSVIVT01035986001,arath.AT2G24960.1,cucsa.Cucsa.271410.1,prupe.7G091300.1" \
    -w /scratch/scannon/legume_genefams3/46_trees_RAxML -s 43_hmmalign_trim2/L.NHRL8 -n L.NHRL8 -p 123

  mkdir 46_trees_RAxML 

  JOBMAX=20
  WORKDIR=/scratch/scannon/lgf4
  ALIGNDIR=43_hmmalign_trim2
  TREEDIR=46_trees_RAxML
  ROOTS="prupe,cucsa,arath,vitvi,solly"
  FILEPAT="L_*"

  nohup commands/do_raxml_outgrp_root.ksh $JOBMAX $WORKDIR $ALIGNDIR $TREEDIR $ROOTS $FILEPAT &

  # Manual cleanup of RAxML files (incomplete notes)
    cd 46_trees_RAxML
      rm *_log* *_pars* *_best*
      mkdir ../46_trees_RAxML_info
      mv *info* ../46_trees_RAxML_info/
      rename 's/RAxML_result.(L.+)/$1.tree/' RAxML_result.*
      cd ..
  
  # Get list of trees not yet calculated
    ls 46_trees_RAxML | perl -pe 's/.tree//' > lis.46_trees_RAxML
    ls 46_trees_FT > lis.46_trees_FT
    comm -23 lis.46_trees_FT lis.46_trees_RAxML > lis.trees_not_calcd

  # Recalculate those 1562 trees
    mkdir 43_hmmalign_trim2_redo 46_trees_RAxML_redo
    perl -pe 's{(.+)}{cp 43_hmmalign_trim2/$1 43_hmmalign_trim2_redo/}' lis.trees_not_calcd \
      > cmd.cp_trees_not_calcd
    sh cmd.cp_trees_not_calcd 


  JOBMAX=10
  WORKDIR=/scratch/scannon/lgf4
  ALIGNDIR=43_hmmalign_trim2_redo
  TREEDIR=46_trees_RAxML_redo
  ROOTS="prupe,cucsa,arath,vitvi,solly"
  FILEPAT="L_*"

  nohup commands/do_raxml_outgrp_root.ksh $JOBMAX $WORKDIR $ALIGNDIR $TREEDIR $ROOTS $FILEPAT &

  # After another round and ~100 trees calculated, 1257 remain. 
  # These have no outgroups, so copy in the FastTree trees.

  mkdir 37_trees_combined
  cp 46_trees_FT/* 37_trees_combined/
  cd 37_trees_combined
    rename 's/(L.*)/$1.tree/' L.*
    cd ..
  cp 46_trees_RAxML/L.* 37_trees_combined/


########## START general procedure for adding or removing taxa from a set of trees
# Notes on pruning/removing taxa, using R ape package and/or phytools (Liam Revel):
# http://blog.phytools.org/2011/03/prune-tree-to-list-of-taxa.html
# http://blog.phytools.org/2014/08/remove-set-of-tips-matching-regular.html

# Starting materials: 
#   proteomes, e.g. 01_proteomes
#   HMMs, e.g. 44_hmm
  
# Add peanut genes. Copy from lis-stage:
# /usr/local/www/data/private/Arachis_hypogaea/arahy.Tifrunner.gnm1.ann1.CCJH
# Files: arahy.Tifrunner.gnm1.ann1.CCJH.protein_primaryTranscript.faa
#        arahy.Tifrunner.gnm1.ann1.CCJH.cds_primaryTranscript.fna
  mv 01_proteomes/arahy.Tifrunner.gnm1.ann1.CCJH.protein_primaryTranscript.faa 01_proteomes/arahy.pep.fa
  mv arahy.Tifrunner.gnm1.ann1.CCJH.cds_primaryTranscript.fna 01_CDS/arahy.CDS.fa
    
# Simplify the deflines:
  perl -pi -e 's/Tifrunner.gnm1.ann1.//' 01_proteomes/arahy.pep.fa
  perl -pi -e 's/Tifrunner.gnm1.ann1.//' 01_CDS/arahy.CDS.fa


##### Add one or more species
  # press hmm database
    cat 44_hmm/* > blastdb/44_hmm_ALL
    hmmpress blastdb/44_hmm_ALL &
    rm blastdb/44_hmm_ALL

  # Search hmm database
    mkdir 60_hmmscan
    SP=arahy
    nohup hmmscan --noali -E 1e-4 --cpu 12 --tblout 60_hmmscan/$SP.x.44_hmm_ALL.pftbl \
      blastdb/44_hmm_ALL 01_proteomes/$SP.pep.fa &

# For the query species, find each gene's best match to a family

  SP=arahy
  INDIR=60_hmmscan

  top_hmmscan_line.awk $INDIR/$SP.x.44_hmm_ALL.pftbl | 
    perl -pe 's/ +/\t/g' | awk -v OFS="\t" '$1!~/^#/ {print $1, $3}' \
    > tmp.lis_top_arahy_in_fams

  # Compare sequences in families calculated earlier, against their HMM, 
  # to get a measure of expected similarity
     
    mkdir 45_fam_vs_hmm
    chmod u+x commands/do_hmmsearch_of_fam_vs_hmm.ksh
    nohup commands/do_hmmsearch_of_fam_vs_hmm.ksh 44_hmm 42_family_fasta 45_fam_vs_hmm 20 &

# Then filter out those not matching with at least 50% of the median match score for the family
  
  cat /dev/null > tab.45_fam_vs_hmm.MAD
  for path in 26_fam_vs_hmm/* ; do
    file=`basename $path`
    cat $path | awk '$1!~/^#/ {print $6}' | median_abs_dev.pl -file $file >> tab.45_fam_vs_hmm.MAD
  done &

  filter_by_median_abs_dev.pl -hmm_file 60_hmmscan/arahy.x.44_hmm_ALL.pftbl \
      -median_file tab.45_fam_vs_hmm.MAD -floor 50 -no_MAD |
      grep -v "#" | sort -k1,1 -k3nr,3nr > 60_hmmscan/arahy.fl50noMAD &

# Take top match - of the hmmscan matches filtered for MAD
  top_line.awk 60_hmmscan/arahy.fl50noMAD | cut -f1,2 > tmp.lis_top_arahy_in_fams


# Combine all into a single list, per gene family
# To do this, first reshape output into a family file with format linke that generated by mcl:
# each line begins with a family ID, and the gene list follows, space-separated.
# then use get_fasta_from_family_file.pl to retrieve sequences for each family.
  cat tmp.lis_top_arahy_in_fams tmp.lis.some_genes_in_fams | sort -k1,1 -k2,2 |
    hash_to_rows_by_1st_col.awk > fam_lists.arahy_and_some_genes_in_fams


# Get fasta sequences for combined (core + leftover + outgroup) families
  
  mkdir 62_family_fasta
  get_fasta_from_family_file.pl -in tmp.all_peptides.fa \
    -fam fam_lists.arahy_and_some_genes_in_fams -out 62_family_fasta &

# Align each family to HMM
  mkdir 63_hmmalign
  nohup commands/hmm_realign.ksh 62_family_fasta 44_hmm 63_hmmalign 20 &

# Remove non-match characters from A2M file
  
  mkdir 63_hmmalign_trim1
  nohup commands/hmm_trim1.ksh 63_hmmalign 63_hmmalign_trim1 20 &

# Remove columns, then cut sequences spanning less than 20% of the alignment
  # NOTE change in alignment-cleaning parameters here: depth 4 and coverage 15, rather than 3 and 20
  mkdir 63_hmmalign_trim2 63_hmmalign_trim2_log
  nohup commands/hmm_trim2.ksh 63_hmmalign_trim1 63_hmmalign_trim2 63_hmmalign_trim2_log 4 15 20 &
                     # hmmtrimdir1   hmmtrimdir2  logdir  min_depth   min_pct_aligned   JOBMAX

#####
# Make trees. Use FastTree on all trees, but planning to replace all but the largest with RAxML 
# for outgroup rooting.
# Create some temp directories. Combine files later.

    mkdir 66_trees_FT
    nohup commands/do_fasttree.ksh 63_hmmalign_trim2 66_trees_FT 20 &

# For remaining families, calculate trees using RAxML with outgroup rooting

  mkdir 66_trees_RAxML 

  JOBMAX=20
  WORKDIR=/scratch/scannon/lgf4
  ALIGNDIR=63_hmmalign_trim2
  TREEDIR=66_trees_RAxML
  ROOTS="prupe,cucsa,arath,vitvi,solly"
  FILEPAT="L_*"

  nohup commands/do_raxml_outgrp_root.ksh $JOBMAX $WORKDIR $ALIGNDIR $TREEDIR $ROOTS $FILEPAT &

########## END general procedure for adding or removing taxa from a set of trees

# Check the proportion of genes for several species that are in gene families
  
  for species in arahy araip cicar glyma lotja medtr phavu vigun; do 
    echo
    echo $species
    grep -c '>' 01_proteomes/$species.pep.fa
  done &
    # arahy 67124
    # araip 41840
    # cicar 28269
    # glyma 56044
    # lotja 39734
    # medtr 50894
    # phavu 27197
    # vigun 29773

  for species in arahy araip cicar glyma lotja medtr phavu vigun; do 
    echo
    echo $species
    cat 62_family_fasta/L_* | grep -c $species
  done &
    # arahy 15025/67124 = 0.2238
    # araip 18747/41840 = 0.4481 
    # cicar 18104/28269 = 0.6404
    # glyma 40550/56044 = 0.7235
    # lotja 18287/39734 = 0.4602
    # medtr 24132/50894 = 0.4742 
    # phavu 22538/27197 = 0.8287 
    # vigun 24236/29773 = 0.8140 

# Observe: this isn't good - especially for arahy
# To try: match the leftover genes against the current gene families.

##### Add leftover genes to families #####
  # Get list of genes in the original proteomes
    mkdir 01_original_proteome_lists
    LC_COLLATE=C
    for species in aradu arahy araip cajca cicar glyma lotja lupan medtr phavu tripr vigan vigra vigun; do
      grep ">" 01_proteomes/$species.pep.fa | perl -pe 's/^>//; s/(^\S+)\s.+/$1/; s/ +//' |
        sort > 01_original_proteome_lists/$species
    done

  # Get list of genes in the current families
    mkdir 01_genes_in_62_lists
    LC_COLLATE=C
    for species in aradu arahy araip cajca cicar glyma lotja lupan medtr phavu tripr vigan vigra vigun; do
      cat 62_family_fasta/* | grep $species | perl -pe 's/^>//; s/ +//' | 
        sort > 01_genes_in_62_lists/$species
    done

  # Get list of genes MISSING FROM the current families
    mkdir 01_genes_NOT_in_62_lists
    for species in aradu arahy araip cajca cicar glyma lotja lupan medtr phavu tripr vigan vigra vigun; do
      comm -23 01_original_proteome_lists/$species 01_genes_in_62_lists/$species \
        > 01_genes_NOT_in_62_lists/$species
    done

  # Get leftover genes
    mkdir 60_leftover_genes
    for path in 01_genes_NOT_in_62_lists/*; do 
      file=`basename $path` 
      echo $file
      get_fasta_subset.pl -input_fas 01_proteomes/$file.pep.fa \
        -output_fas 60_leftover_genes/$file -lis 01_genes_NOT_in_62_lists/$file 
    done

  # Search hmm database
    mkdir 60_hmmscan

    HMMDB=blastdb/44_hmm_ALL
    FASDIR=60_leftover_genes
    TBLOUTDIR=60_hmmscan
    JOBMAX=14
    nohup commands/do_hmmscan_of_proteomes_vs_hmm.ksh $HMMDB $FASDIR $TBLOUTDIR $JOBMAX &

  # get top hits
    mkdir 60_hmmscan_top
    for path in 60_hmmscan/*; do 
      file=`basename $path`
      top_hmmscan_line.awk $path > 60_hmmscan_top/$file
    done

  # filter by median absolute deviation per family
    mkdir 60_hmmscan_top_MAD
    for path in 60_hmmscan_top/*; do 
      file=`basename $path`
      filter_by_median_abs_dev.pl -hmm_file $path -floor 50 -no_MAD \
      -median_file tab.45_fam_vs_hmm.MAD > 60_hmmscan_top_MAD/$file
    done

  # Lists of genes in each family in 62
    
    for path in 62_family_fasta/*; do
      FAM=`basename $path`
      awk -v family=$FAM '$1~/>/ {print family "\t" substr($1,2)}' $path
    done > hsh.70_genes_in_fams_62


  # Lists of genes to add to families, per 60_hmmscan_top_MAD
    cat 60_hmmscan_top_MAD/* | awk '$9~/MATCHES/ {print $1 "\t" $2}' > hsh.70_genes_to_add

  # Combine hash-lists to get larger set of genes for augmented families
  # The output is a family file with format linke that generated by mcl.
    cat hsh.70_genes_in_fams_62 hsh.70_genes_to_add | sort -u |
      hash_to_rows_by_1st_col.awk > fam_lists.70_genes_in_fams_augmented

  # then use get_fasta_from_family_file.pl to retrieve sequences for each family.
    mkdir 62_family_fasta
    get_fasta_from_family_file.pl -in tmp.all_peptides.fa \
      -fam fam_lists.70_genes_in_fams_augmented -out 62_family_fasta &
    # NOTE: The prior fasta files in 62_family_fasta were overwritten here with the augmented files.
    # From this point forward, redo alignment and tree work in the 6X_ series.

# Align each family to HMM
  mkdir 63_hmmalign
  nohup commands/hmm_realign.ksh 62_family_fasta 44_hmm 63_hmmalign 20 &

# Remove non-match characters from A2M file
  
  mkdir 63_hmmalign_trim1
  nohup commands/hmm_trim1.ksh 63_hmmalign 63_hmmalign_trim1 20 &

# Remove columns, then cut sequences spanning less than 20% of the alignment
  # NOTE change in alignment-cleaning parameters here: depth 4 and coverage 15, rather than 3 and 20
  mkdir 63_hmmalign_trim2 63_hmmalign_trim2_log
  nohup commands/hmm_trim2.ksh 63_hmmalign_trim1 63_hmmalign_trim2 63_hmmalign_trim2_log 4 15 20 &
                     # hmmtrimdir1   hmmtrimdir2  logdir  min_depth   min_pct_aligned   JOBMAX


##############################
# Get leftover genes

  # Get list of genes in the original proteomes
    mkdir 01_original_proteome_lists
    LC_COLLATE=C
    for species in aradu arahy araip cajca cicar glyma lotja lupan medtr phavu tripr vigan vigra vigun; do
      grep ">" 01_proteomes/$species.pep.fa | perl -pe 's/^>//; s/(^\S+)\s.+/$1/; s/ +//' |
        sort > 01_original_proteome_lists/$species
    done

  # Get list of genes in the current families
    mkdir 01_genes_in_62_lists
    LC_COLLATE=C
    for species in aradu arahy araip cajca cicar glyma lotja lupan medtr phavu tripr vigan vigra vigun; do
      cat 62_family_fasta/* | grep $species | perl -pe 's/^>//; s/ +//' | 
        sort > 01_genes_in_62_lists/$species
    done

  # Get list of genes MISSING FROM the current families
    mkdir 01_genes_NOT_in_62_lists
    for species in aradu arahy araip cajca cicar glyma lotja lupan medtr phavu tripr vigan vigra vigun; do
      comm -23 01_original_proteome_lists/$species 01_genes_in_62_lists/$species \
        > 01_genes_NOT_in_62_lists/$species
    done

  # Get leftover genes
    mkdir 60_leftover_genes
    for path in 01_genes_NOT_in_62_lists/*; do 
      file=`basename $path` 
      echo $file
      get_fasta_subset.pl -input_fas 01_proteomes/$file.pep.fa \
        -output_fas 60_leftover_genes/$file -lis 01_genes_NOT_in_62_lists/$file 
    done &

    mv 60_leftover_genes 01_proteomes_leftover

# Run blastp for all pairs
  # Make blast databases
  mkdir blastdb blastout
  for path in 01_proteomes_leftover/*; do 
    file=`basename $path`; 
    makeblastdb -in $path -dbtype prot -hash_index -parse_seqids -out blastdb/$file &
  done

  for path in 01_proteomes_out/*; do 
    file=`basename $path`; 
    makeblastdb -in $path -dbtype prot -hash_index -parse_seqids -out blastdb/$file &
  done


# Blast all by all
# assumes existing standard directories: blastdb blastout
# params: $FASTADIR $JOBMAX
  FASTADIR1=01_proteomes_leftover
  FASTADIR2=01_proteomes_leftover
  JOBMAX=30
  nohup ../commands/do_blast_all.ksh $FASTADIR1 $FASTADIR2 $JOBMAX &

# Blast legumes against outgroups
# assumes existing standard directories: blastdb blastout
# params: $FASTADIR $JOBMAX
  FASTADIR1=01_proteomes_leftover
  FASTADIR2=01_proteomes_out
  JOBMAX=20
  nohup ../commands/do_blast_all.ksh $FASTADIR1 $FASTADIR2 $JOBMAX &

# Get top two hits per query - except for self-comparisons, where we only want one top (non-self) match
# Also: filter matches by percent query alignment and percent subject alignment
  
  mkdir 02_pep_len blastout_cov

  mv 01_proteomes_out/* 01_proteomes_leftover/

  for path in 01_proteomes_leftover/*; do 
    file=`basename $path`
    seqlen.awk $path > 02_pep_len/$file
  done &


  # Call script with percent coverage = 50 and percent identity 60
  # NOTE: Filenames in the BLAST results are assumed to have the format medtr.x.lotja.blp
  # where the query and subject species names each have five characters and are separated by .x.

    LEN_DIR=02_pep_len
    COV_DIR=blastout_cov
    COV_PCT=50
    ID_PCT=60
    JOBMAX=20

    nohup ../commands/do_blast_coverage_filter.ksh $LEN_DIR $COV_DIR $COV_PCT $ID_PCT $JOBMAX &

# Pick top two hits per query - except top non-self single hit for self-comparisons

  mkdir blastout_top2

    for path in blastout_cov/*; do
      file=`basename $path`
      echo $file
      top_2_lines_no_self.awk $path > blastout_top2/$file.top &
    done &
     
  # Prune the self-hit files to one hit per query, for the legumes (skipping the outgroups)
    for SP in aradu araip arahy cajca cicar glyma lotja lupan medtr phavu tripr vigan vigra vigun ; do
      file="$SP.x.$SP"
      echo $file
      top_line_no_self.awk blastout_cov/$file > blastout_top2/$file.top &
    done &


# For all comparisons with glyma as subject, pick top FOUR non-self hits 
# (non-self only being relevant for glyma x glyma)
  nohup commands/get_top4_nonself_blast_hits.ksh &
    # three lines for glyma.x.glyma
      top_3_lines_no_self.awk MAX=3 blastout_cov/glyma.x.glyma > blastout_top2/glyma.x.glyma.top &
    
    # skips outgroups - and also glyma as query
    # For all comparisons with glyma as subject, pick top FOUR non-self hits 
    # (non-self only being relevant for glyma x glyma)
    for SP in aradu araip arahy cajca cicar lotja lupan medtr phavu tripr vigan vigra vigun ; do
      file="$SP.x.glyma"
      echo $file
      top_4_lines.awk blastout_cov/$file > blastout_top2/$file.top &
    done

# CDS and peptides come from /data/gene_families/02_CDS/ and /data/gene_families/01_proteomes/
  cat ../01_proteomes/*fa > ../tmp.all_peptides.fa
  cat ../01_CDS/*fa > ../tmp.all_CDS.fa

# In preparation for calculating Ks, process top-two blast hits, reporting unique query-target list, where query is always lexically before target
  cat blastout_top2/* | awk '$1<=$2 {print $1 "\t" $2} $2<$1 {print $2 "\t" $1}' | 
    sort -u > lis.top_two_blast_hits &

# Get CDS and peptide files for those gene lists. 
  mkdir 03_pair_lists
  perl -pe 's/(\w\w\w\w\w)\.(\S+)\t(\w\w\w\w\w)\.(\S+)$/$1.$3\t$1.$2\t$3.$4/' lis.top_two_blast_hits | 
    sort > lis.ks_to_be_calc_for_splitting &
    #NOTE: lis.ks_to_be_calc_for_splitting should be a three-column file, with the first being e.g. glyma.medtr

  split_table_to_files.pl lis.ks_to_be_calc_for_splitting 03_pair_lists 1 &

# For synonymous_calc.py, print query and subject in one column: Q S Q S Q S ...
  perl -pi -e 's/^\S+\t(\S+)\t(\S+)/$1\n$2/' 03_pair_lists/* &

# Get fasta sequence, given the lists above. Do this in a script for parallelization
  mkdir 04_pair_peptide_fasta 04_pair_CDS_fasta 

  INFASTA=../tmp.all_CDS.fa
  LISTDIR=03_pair_lists
  OUTDIR=04_pair_CDS_fasta
  JOBMAX=20
  nohup ../commands/get_fasta.ksh $INFASTA $LISTDIR $OUTDIR $JOBMAX &

  INFASTA=../tmp.all_peptides.fa
  LISTDIR=03_pair_lists
  OUTDIR=04_pair_peptide_fasta
  JOBMAX=20
  nohup ../commands/get_fasta.ksh $INFASTA $LISTDIR $OUTDIR $JOBMAX &

# Calculate ks
  mkdir 05_ks
  
  PEPDIR=04_pair_peptide_fasta
  CDSDIR=04_pair_CDS_fasta
  WORKDIR=05_ks
  MAINDIR=/scratch/scannon/lgf4/LEFTOVERS
  JOBMAX=60

  nohup ../commands/syn_calc.ksh $PEPDIR $CDSDIR $WORKDIR $MAINDIR $JOBMAX &


  # Fields in the Ks output are:
    #   1     2     3      4      5      6
    #   name  name  dS-yn  dN-yn  dS-ng  dN-ng
    # Use 5: dS-Nei-Gojobori

  # Parse into four-column gene-pair format:
    #species-pair   gene             gene                   dS-ng 
    aradu.vigun  aradu.Aradu.000JC  vigun.Vigun05g296900.1  0.8964
    aradu.vigun  aradu.Aradu.001N3  vigun.Vigun10g098600.1  0.5134
    aradu.vigun  aradu.Aradu.002J3  vigun.Vigun07g072400.1  0.5779

  # Also, skip outgroup comparisons: arath, cucsa, prupe, solly, vitvi

    cat 05_ks/*.ks | 
      grep -v arath | grep -v cucsa | grep -v prupe | grep -v solly | grep -v vitvi |
      perl -pe 's/[;,]/\t/g; s/^(\w\w\w\w\w)\.(\S+)\t(\w\w\w\w\w)\.(\S+)\t/$1\t$3\t$1.$2\t$3.$4\t/' | 
      awk -v OFS="\t" '$1 <= $2 && $7>=0 && $7<2.5 {print $1 "." $2, $3, $4, $7 }
                       $2 <  $1 && $7>=0 && $7<2.5 {print $2 "." $1, $4, $3, $7 }' |
       LC_COLLATE=C sort -u -k1,1 -k2,2 > tab.all_Ks_before_splitting_by_sp_pair &

  # Split into species pairs
    mkdir 05_ks_parsed
    split_table_to_files.pl tab.all_Ks_before_splitting_by_sp_pair 05_ks_parsed 1 &

# Filter the Ks output relative to Ks peaks. Filter separately for each species pair, because the expected 
# Ks ranges are different for different species pairs. The file with manually-assessed Ks thresholds is
# WGD_peaks2_leg_manual.tab, with format 
    aradu.lotja 0.55

  # Output rows consist of a gene pair and a score between 10 and 100 (default values), 
  # with higher numbers corresponding with lower Ks (and greater similarity):
  # transformed_ks = 
  #     min_transform + (max_transform-min_transform)*((adjusted_cutoff-ks)/adjusted_cutoff);

  # Example for one species-pair: 
    filter_by_ks.pl -ks 05_ks_parsed/araip.vigra -cut ../stats/WGD_peaks2_leg_manual.tab -scale_factor 1.5

  # Do the filtering for all species-pairs:
    for path in 05_ks_parsed/* ; do 
      sp_pair=`basename $path`
      filter_by_ks.pl -ks $path -cut ../stats/WGD_peaks2_leg_manual.tab -scale_factor 1.5
    done > tab.all_pairs_filtered_by_ks  &


# Cluster - testing with several inflation values
    nohup mcl tab.all_pairs_filtered_by_ks -I 1.1 -te 24 --abc &
    nohup mcl tab.all_pairs_filtered_by_ks -I 1.2 -te 12 --abc &
    nohup mcl tab.all_pairs_filtered_by_ks -I 1.4 -te 12 --abc &
    nohup mcl tab.all_pairs_filtered_by_ks -I 2.0 -te 12 --abc &
    nohup mcl tab.all_pairs_filtered_by_ks -I 3.0 -te 12 --abc &
    nohup mcl tab.all_pairs_filtered_by_ks -I 4.0 -te 12 --abc &

# Check composition summaries (e.g. counts per species, and numbers of species counts = 1 or 2)
# Also prune small and non-diverse families, using MIN_FAM_SIZE and MIN_SP_COUNT.
# Print to three files: KEEP, CUT, and ct_per_sp
    
    mkdir stats

    export MIN_FAM_SIZE=4
    export MIN_SP_COUNT=3
    for I in 11 12 14 20 30 40; do
      export file="out.tab.all_pairs_filtered_by_ks.I$I"
      cat $file | perl -ane '
          BEGIN{
            $MIN_FAM_SIZE=$ENV{"MIN_FAM_SIZE"};
            $MIN_SP_COUNT=$ENV{"MIN_SP_COUNT"};
            $FILE=$ENV{"file"};
            my $fam_num = 0;
            @species = qw(aradu araip cajca cicar glyma lotja lupan medtr phavu tripr vigan vigra vigun);
            open ($OUTCUT, ">", "$FILE.CUT") or die "cant open out lis.$FILE.CUT, $!\n";
            open ($OUTKEEP, ">", "$FILE.KEEP") or die "cant open out $FILE.KEEP, $!\n";
            open ($OUTCOUNTS, ">", "ct_per_sp.$FILE") or die "cant open out ct_per_sp.$FILE, $!\n";

            print $OUTCOUNTS "num.counts____\t", join("\t", @species),"\n";
          my $out_line;
          my %ct_sp;
          for my $sp (@species){
            $ct_sp{$sp} = 0;
            for my $seq (@F){ if ($seq =~ /$sp/) {$ct_sp{$sp}++ } }
            $out_line .= "$ct_sp{$sp}\t";
            if ($ct_sp{$sp}>0){$ct_present++}
            if ($ct_sp{$sp}<2){$ct_below++}
            if ($ct_sp{$sp}==2){$ct_eq++}
            if ($ct_sp{$sp}>2){$ct_above++}
          }
          if (scalar(@F) >= $MIN_FAM_SIZE && $ct_present >= $MIN_SP_COUNT) {
            $gene_ct = scalar(@F);
            $fam_num++;
            my $fam_id_and_stats = "$fam_num.$gene_ct.$ct_present.$ct_below.$ct_eq.$ct_above";
            print $OUTCOUNTS "$fam_id_and_stats\t$out_line\n";
            print $OUTKEEP $fam_id_and_stats, "\t", join("\t",@F),"\n";
          }
          else { # print IDs from small families
            print $OUTCUT join("\n",@F),"\n";
          } ' # end of perl script
      done &


# Examine family size histograms:
  for file in *KEEP; do echo $file; awk '{print NF}' $file | histogram -n -s2 | histplot -d 10 | head -40 ; done

  # Based on these, use -I 1.2 for these leftover genes (up from 1.1 for original families)
out.tab.all_pairs_filtered_by_ks.I12.KEEP
bin   abcdefghijKLMNOPQRSTabcdefghijKLMNOPQRSTabcdefghijKLMNOPQRSTabcdefghijKLMNOPQRST
4.00  ......................................................
6.00  .....................................................................................
8.00  .....................................................
10.00 .................................
12.00 ......................
14.00 ..................
16.00 .............
18.00 ............
20.00 .........
22.00 ......
24.00 .....
26.00 .....
28.00 ...
30.00 ...
32.00 ..
34.00 ...
36.00 .
38.00 ..
40.00 .
42.00 .
44.00 .
46.00 .
48.00 .
50.00 .


# Get summary stats: histograms of family sizes - summarizing from stats/ct_per_sp.$file above
    base="out.tab.all_pairs_filtered_by_ks.I"
    for I in 11 12 14 20 30 40; do
      grep -v counts "stats/$base$I.KEEP" | cut -f1 | cut -f2 -d'.' | 
        histogram -n -s1 | histplot -d 10 | head -50 > stats/hist_ct_per_sp.$base$I
    done

    for file in stats/hist_ct_per_sp.out.tab.all_pairs_filtered_by_ks*; do 
      echo; echo $file; cat $file; 
    done > stats/summary_hist_ct_per_sp.out.tab.all_pairs_filtered_by_ks_f$MIN_FAM_SIZE.s$MIN_SP_COUNT.txt 


  # Get coarse summary by some meaningful categories - on the pruned families (.KEEP)

    export MIN_FAM_SIZE=4
    export MIN_SP_COUNT=3
    base="out.tab.all_pairs_filtered_by_ks.I"
    for I in 11 12 14 20 30 40; do
      cat stats/$base$I.KEEP |
      perl -lane '$ct=scalar(@F)-1; # subtract family ID, e.g. 16200.6.6.17.0.0
                  if($ct>$max){$max=$ct}
                  if($ct<10){$_lt_10++} 
                  elsif($ct>=10 && $ct<20){$_10_to_20++} 
                  elsif($ct>=20 && $ct<40){$_20_to_40++} 
                  elsif($ct>=40 && $ct<80){$_40_to_80++} 
                  elsif($ct>=80 && $ct<200){$_80_to_200++} 
                  else{$_ge_200++} 
                  $tot_genes += $ct;
                  $tot_fams ++;
                  END{print "_0 to 10\t$_lt_10"; 
                      print "10 to 20\t$_10_to_20"; 
                      print "20 to 40\t$_20_to_40"; 
                      print "40 to 80\t$_40_to_80"; 
                      print "80 to 200\t$_80_to_200"; 
                      print "_ >= 200\t$_ge_200";
                      print "_    max\t$max";
                      print "_  >= 10\t", $_lt_10+$_10_to_20+$_20_to_40+$_40_to_80+$_80_to_200+$_ge_200;
                      print "tot_genes\t$tot_genes";
                      print "tot_fams\t$tot_fams";
                      }
                  ' > stats/ct_coarse.$base$I.KEEP.f$MIN_FAM_SIZE.s$MIN_SP_COUNT.txt
    done
  
    for file in stats/ct_coarse.out.tab.all_pairs_filtered_by_ks.I*.f$MIN_FAM_SIZE.s$MIN_SP_COUNT.txt; do 
      echo; echo $file; cat $file; 
    done > stats/summary.ct_coarse.out.tab.all_pairs_filtered_by_ks.f$MIN_FAM_SIZE.s$MIN_SP_COUNT.txt

##### MAJOR PRODUCT: initial gene families
# Give this file a new name:
  cp stats/out.tab.all_pairs_filtered_by_ks.I12.KEEP legume_gene_families_LEFTOVERS_v1_I12


# Get fasta sequence from the labeled mcl family file. Do this with get_fasta_from_family_file.pl
# in order not to have to re-load the source family file more than once.

  mkdir 12_family_fasta

  nohup get_fasta_from_family_file.pl -in ../tmp.all_peptides.fa  -out 12_family_fasta  -fam legume_gene_families_LEFTOVERS_v1_I12 &

# Remove terminal X and internal XXs (there are quite a few in these sequences). Leave one internal X 
  cd 12_family_fasta
    perl -pi -e 's/^([A-Z][^X]+)X+$/$1/; s/^([A-Z].+[^X]+)XX+([A-Z].+[^X]+)$/$1X$2/g' * &
    mkdir ../12_family_fasta_tmp
    for file in *; do sed '/^$/d' $file > ../12_family_fasta_tmp/$file; done &
    mv ../12_family_fasta_tmp/* .
    cd ..
  # NOTE: This didn't work fully as intended. Some XXs were removed but not all. Did some additional work

# Assign IDs

  # Find number of IDs to generate
    COUNT=`ls 12_family_fasta | wc -l`
    echo $COUNT
      3569
  # Get list of IDs already in use:
    ls ../62_family_fasta > lis.fam_IDs_round1

  make_random_IDs.pl -count $COUNT -width 6 -in_file lis.fam_IDs_round1 | 
    sort | perl -pe 's/(\w+)/L_$1/' > lis.fam_names_round2

  # Modify temporary names, prefixing them with gene count number, 
  # to help assign random IDs in order by family size
  mkdir 12_family_fasta.tmp
  cd 12_family_fasta
    for file in *; do  
      count=`grep -c '>' $file`; 
      echo $count.$file; 
      cp $file ../12_family_fasta.tmp/$count.$file
    done > ../lis.count_prefix.12_family_fasta 

    cd ..
  
  ls 12_family_fasta.tmp > lis.12_family_fasta.tmp

  mkdir 12_family_fasta.tmp2
  sort -k1nr,1nr -t'.' -o lis.count_prefix.12_family_fasta lis.count_prefix.12_family_fasta
  paste lis.count_prefix.12_family_fasta lis.fam_names_round2 | 
    perl -pe 's{(\S+)\t(\S+)}{mv 12_family_fasta.tmp/$1 12_family_fasta.tmp2/$2}' > cmd.rename

  sh cmd.rename
  rm 12_family_fasta/*
  mv 12_family_fasta.tmp2/* 12_family_fasta/
  rmdir 12_family_fasta.tmp*


# Do alignments
  mkdir 13_align
  
  INDIR=12_family_fasta
  OUTDIR=13_align
  JOBMAX=40

  nohup ../commands/do_align.ksh $INDIR $OUTDIR $JOBMAX &


#####
# Build HMMs
  mkdir 14_hmm

  INDIR=13_align
  OUTDIR=14_hmm
  JOBMAX=40
        
  nohup ../commands/do_hmmbuild.ksh $INDIR $OUTDIR $JOBMAX &

# Copy HMMs from first round into new combined hmm directory
  mkdir 44_hmm
  cp ../44_hmm/* 44_hmm/
  cp 14_hmm/* 44_hmm/

# Search all proteomes against combined HMMs, in order to allow free assortment among them

  # press hmm database
  mkdir hmmdb
  cat 44_hmm/* > hmmdb/44_hmm_ALL
  hmmpress hmmdb/44_hmm_ALL 
  rm hmmdb/44_hmm_ALL


  mkdir 50_hmmscan

  HMMDB=hmmdb/44_hmm_ALL
  FASDIR=../01_proteomes
  TBLOUTDIR=50_hmmscan
  JOBMAX=20

  nohup ../commands/do_hmmscan_of_proteomes_vs_hmm.ksh $HMMDB $FASDIR $TBLOUTDIR $JOBMAX &

  # tracking output for each output file, while the hmmscan is running:
    for file in *fa; do 
      FILE=$file; 
      COUNT=`top_hmmscan_line.awk $file | wc -l`; 
      echo "$FILE $COUNT"; 
    done

# For each family, calculate the hmmscan median and median-absolute-deviation (MAD) scores
  cat 50_hmmscan/* | top_hmmscan_line.awk | 
    awk '$1!~/^#/ {print $1 "\t" $6}' | sort -k1,1 -k2nr,2nr | 
    median_abs_dev.pl -col 2 | sort -k1,1 > tab.50_MAD

# Take the top hmmscan match
  cat 50_hmmscan/* | grep -v "#" | top_hmmscan_line.awk > tab.50_top_hmmscan

# Filter by percent of median for the family
  # DON'T use this one:
  filter_by_median_abs_dev.pl -hmm tab.50_top_hmmscan -med tab.50_MAD -floor 50 -no_MAD |
         sort -k1,1 -k10nr,10nr  > tab.50_top_hmmscan_filt_median_fl50 &

  # Go with this one:
  filter_by_median_abs_dev.pl -hmm tab.50_top_hmmscan -med tab.50_MAD -floor 40 -no_MAD |
         sort -k1,1 -k10nr,10nr  > tab.50_top_hmmscan_filt_median_fl40 &

  filter_by_median_abs_dev.pl -hmm tab.50_top_hmmscan_lgf5_legumeonly -med tab.50_MAD -floor 40 -no_MAD |
         sort -k1,1 -k10nr,10nr  > tab.50_top_hmmscan_filt_median_fl40_lgf5_legumeonly &

# Check proportions matching families at this cutoff, by species
  awk '$9!~/verdict/ {print $2 "\t" $9}' tab.50_top_hmmscan_filt_median_fl50 | 
    perl -pe 's/(^\w\w\w\w\w)\.\S+\t(\S+)/$1\t$2/' | sort | uniq -c | 
    awk -v OFS="\t" '$3~/doesnt/ {no=$1} $3~/MATCHES/ {yes=$1; sp=$2; print sp, no, yes, no/yes}'
      aradu 5496  27175 0.202245
      arahy 9455  53539 0.1766
      araip 7438  28992 0.256554
      arath 3880  21143 0.183512
      cajca 4248  31188 0.136206
      cicar 2349  23641 0.0993613
      cucsa 2403  17596 0.136565
      glyma 4931  47043 0.104819
      lotja 8493  24392 0.348188
      lupan 1982  29756 0.0666084
      medtr 7247  32994 0.219646
      phavu 1007  25463 0.0395476
      prupe 3121  21332 0.146306
      solly 6849  22198 0.308541
      tripr 8982  27711 0.324131
      vigan 2614  22132 0.11811
      vigra 2298  19616 0.117149
      vigun 881 27101 0.032508
      vitvi 5201  17741 0.293163

  awk '$9!~/verdict/ {print $2 "\t" $9}' tab.50_top_hmmscan_filt_median_fl40 | 
    perl -pe 's/(^\w\w\w\w\w)\.\S+\t(\S+)/$1\t$2/' | sort | uniq -c | 
    awk -v OFS="\t" '$3~/doesnt/ {no=$1} $3~/MATCHES/ {yes=$1; sp=$2; print sp, no, yes, no/yes}'
      aradu 4493  28178 0.159451
      arahy 7797  55197 0.141258
      araip 6165  30265 0.203701
      arath 2814  22209 0.126705
      cajca 3080  32356 0.095191
      cicar 1758  24232 0.0725487
      cucsa 1699  18300 0.0928415
      glyma 3930  48044 0.0818
      lotja 6899  25986 0.265489
      lupan 1502  30236 0.0496759
      medtr 6253  33988 0.183977
      phavu 783 25687 0.0304823
      prupe 2413  22040 0.109483
      solly 5506  23541 0.23389
      tripr 7578  29115 0.260278
      vigan 2015  22731 0.0886455
      vigra 1784  20130 0.0886239
      vigun 689 27293 0.0252446
      vitvi 4230  18712 0.226058

# Suppress sequences not matching -- i.e. $9~/doesnt/ -- and cases where there are fewer than 3 matching legume sequences.
  
  # get count of legumes per family 
  cat tab.50_top_hmmscan_filt_median_fl40 | 
    awk '$9~/MATCHES/ && $2~/aradu|arahy|araip|cajca|cicar|glyma|lotja|lupan|medtr|phavu|tripr|vigan|vigra|vigun/ {print $1}' | 
    uniq -c > tmp.ct_by_species.top_hmmscan_filt_median_fl40

  # Join the family count to the HMM match lines, then exclude outgroups in the output (add the top outgroup later)
  # and filter by number of legume sequences in the family (require at least 3)
  sort tab.50_top_hmmscan_filt_median_fl40 | join -t $'\t' - tmp.ct_by_species.top_hmmscan_filt_median_fl40 | 
    awk '$9~/MATCHES/ && $2!~/arath|cucsa|prupe|solly|vitvi/ && $12>=3' > tmp.ct_by_species.top_hmmscan_filt_median_fl40.filt.leg

  # Extract the top outgroup sequences per family. The perl regex joins the species to the gene family name;
  # then sort by the first column (gene_fam--species) and score, and take the top line per first column.
  sort tab.50_top_hmmscan_filt_median_fl40 | join -t $'\t' - tmp.ct_by_species.top_hmmscan_filt_median_fl40 | 
    awk '$9~/MATCHES/ && $2~/arath|cucsa|prupe|solly|vitvi/ && $12>=3' | # require match to outgroup
    perl -pe 's/(^\S+)\t(\w\w\w\w\w)(\.\S+)/$1--$2\t$2$3/' | # add species to gene-family name, e.g. L_001QTQ--arath
    sort -k1,1 -k2,2 -k3nr,3nr | top_line.awk | # pick top-scoring outgroup per gene family
    perl -pe 's/(^\w+)--(\w{5})\t/$1\t/' | # remove species from gene-family name, e.g. L_001QTQ--arath --> L_001QTQ
    sed '' > tmp.ct_by_species.top_hmmscan_filt_median_fl40.filt.out

  # Combine legume-only and top-outgroup files, then
  # Reshape output into a family file with format like that generated by mcl:
  # each line begins with a family ID, and the gene list follows, space-separated.
  # then use get_fasta_from_family_file.pl to retrieve sequences for each family.

  cat tmp.ct_by_species.top_hmmscan_filt_median_fl40.filt.leg \
    tmp.ct_by_species.top_hmmscan_filt_median_fl40.filt.out |
    sort -k1,1 -k3nr,3nr | cut -f1,2 | 
    hash_to_rows_by_1st_col.awk > tab.50_top_hmmscan_filt_median_fl40.fam

    # MAJOR PRODUCT: gene family file, ready for getting fasta sequences
    # tab.50_top_hmmscan_filt_median_fl40.fam

# Get fasta sequences for combined (core + leftover) families
  mkdir 52_family_fasta
  nohup get_fasta_from_family_file.pl -in ../tmp.all_peptides.fa -fam tab.50_top_hmmscan_filt_median_fl40.fam -out 52_family_fasta &

# Check gene counts per species
  cat 52_family_fasta/* | grep '>' | cut -b2-6 | sort | uniq -c | awk '{print $2 "\t" $1}' > tab.in_families
  cat ../01_proteomes/* | grep '>' | cut -b2-6 | sort | uniq -c | awk '{print $2 "\t" $1}' > tab.in_proteomes

  paste tab.in_families tab.in_proteomes | 
    awk '$1!~/arath|cucsa|prupe|solly|vitvi/ {ctFam+=$2; ctProt+=$4; printf "%s\t%i\t%i\t%2.1f\n", $1, $2, $4, 100*$2/$4} 
         END{printf "%s\t%i\t%i\t%2.1f\n", "total", ctFam, ctProt, 100*ctFam/ctProt }'
      # aradu  28154   36734  76.6
      # arahy  55147   67124  82.2
      # araip  30225   41840  72.2
      # cajca  32341   40071  80.7
      # cicar  24210   28269  85.6
      # glyma  48008   56044  85.7
      # lotja  25953   39734  65.3
      # lupan  30222   32865  92.0
      # medtr  33935   50894  66.7
      # phavu  25677   27197  94.4
      # tripr  29078   39947  72.8
      # vigan  22707   26857  84.5
      # vigra  20112   22368  89.9
      # vigun  27286   29773  91.6
      # total  433055 539717  80.2

# Irritatingly, family sizes don't correspond with list order of filenames, because
# the late-added families have a different size distribution than the initial set, 
# and because the re-assortment of all sequences during the hmmsearch changed family sizes.
# Fix this by renaming
  mkdir 52_family_fasta_tmp
  ls 52_family_fasta > lis.52_family_fasta
  grep -c '>' 52_family_fasta/* | perl -pe 's{^\w+/(\w+):(\d+)}{$1\t$2}' | sort -k2nr > lis.50_fam_fasta_reordered
  paste lis.50_fam_fasta_reordered lis.52_family_fasta |
    awk '{print "cp 52_family_fasta/" $1 " 52_family_fasta_tmp/" $3}' > cmd.rename_52_family_fasta
  sh cmd.rename_52_family_fasta
  rm -rf 52_family_fasta
  mv 52_family_fasta_tmp 52_family_fasta

# do alignments
  mkdir 53_align

  INDIR=52_family_fasta
  OUTDIR=53_align
  JOBMAX=30

  nohup ../commands/do_align.ksh $INDIR $OUTDIR $JOBMAX &

# build HMMs. Note: these HMMs should be essentially the same as 44_hmm, but with the new names.
  mkdir 54_hmm

  INDIR=53_align
  OUTDIR=54_hmm
  JOBMAX=30

  nohup ../commands/do_hmmbuild.ksh $INDIR $OUTDIR $JOBMAX &

# Generate a consensus sequence for each HMM

  mkdir 54_hmmemit

  HMMDIR=54_hmm
  HMMEMITDIR=54_hmmemit
  JOBMAX=10

  nohup ../commands/do_hmmemit.ksh $HMMDIR $HMMEMITDIR $JOBMAX &

  perl -pi -e 's/-consensus//' 54_hmmemit/*


# Align each family to HMM
  mkdir 55_hmmalign

  FASTADIR=52_family_fasta
  HMMDIR=54_hmm
  HMMALIGNDIR=55_hmmalign
  JOBMAX=30

  nohup ../commands/hmm_realign.ksh $FASTADIR $HMMDIR $HMMALIGNDIR $JOBMAX &


# Remove non-match characters from A2M file
  
  mkdir 55_hmmalign_trim1

  HMMALIGNDIR=55_hmmalign
  HMMTRIMDIR=55_hmmalign_trim1
  JOBMAX=30

  nohup ../commands/hmm_trim1.ksh $HMMALIGNDIR $HMMTRIMDIR $JOBMAX &


# Remove columns, then cut sequences spanning less than 20% of the alignment
  # NOTE change in alignment-cleaning parameters here: depth 3 and coverage 20 AND pct_depth 20
  mkdir 55_hmmalign_trim2 55_hmmalign_trim2_log

  HMMTRIMDIR1=55_hmmalign_trim1
  HMMTRIMDIR2=55_hmmalign_trim2
  LOGDIR=55_hmmalign_trim2_log
  MIN_DEPTH=3
  MIN_PCT_DEPTH=20
  MIN_PCT_ALIGNED=20
  JOBMAX=20

  nohup ../commands/hmm_trim2.ksh $HMMTRIMDIR1 $HMMTRIMDIR2 $LOGDIR $MIN_DEPTH $MIN_PCT_DEPTH $MIN_PCT_ALIGNED $JOBMAX &

#####
# Make trees. Use FastTree on all trees, but planning to replace all but the largest with RAxML 
# for outgroup rooting.
# Create some temp directories. Combine files later.

  mkdir 56_trees_FT
  
  HMMALIGNTRIMMEDDIR=55_hmmalign_trim2
  TREEDIR=56_trees_FT
  JOBMAX=20

  nohup ../commands/do_fasttree.ksh $HMMALIGNTRIMMEDDIR $TREEDIR $JOBMAX &

# For remaining families, calculate trees using RAxML with outgroup rooting

  mkdir 56_trees_RAxML 

  JOBMAX=20
  WORKDIR=/scratch/scannon/lgf4/LEFTOVERS
  ALIGNDIR=55_hmmalign_trim2
  TREEDIR=56_trees_RAxML
  ROOTS="prupe,cucsa,arath,vitvi,solly"
  FILEPAT="L_*"

  nohup ../commands/do_raxml_outgrp_root.ksh $JOBMAX $WORKDIR $ALIGNDIR $TREEDIR $ROOTS $FILEPAT &


  mkdir 56_trees_RAxML_big

  JOBMAX=20
  WORKDIR=/scratch/scannon/lgf4/LEFTOVERS
  ALIGNDIR=55_hmmalign_trim2_big
  TREEDIR=56_trees_RAxML_big
  ROOTS="prupe,cucsa,arath,vitvi,solly"
  FILEPAT="L_*"

  nohup ../commands/do_raxml_outgrp_root.ksh $JOBMAX $WORKDIR $ALIGNDIR $TREEDIR $ROOTS $FILEPAT &


# Clean up the RAxML files
  cd 56_trees_RAxML
    rm *_bestTree* *_info* *_log* *_parsimony*
    rename 's/RAxML_result.//' RAxML_result.*
    cd ..

# Combine RAxML and FastTree trees, replacing the FastTree instances RAxMA where both were calculated. 
# FastTree handles families that couldn't be rooted and that were skipped because of size.

  mkdir 56_trees_combined
  cp 56_trees_FT/L_* 56_trees_combined/
  cp 56_trees_RAxML/L_* 56_trees_combined/


# Check composition summaries (e.g. counts per species, and numbers of species counts = 1 or 2)
# Also prune small and non-diverse families, using MIN_FAM_SIZE and MIN_SP_COUNT.
# Print to three files: KEEP, CUT, and ct_per_sp
    
    OUTFILE=summary.55_hmmalign_trim2.tab

    perl -e '
        @species = qw(aradu arahy araip cajca cicar glyma lotja lupan medtr phavu tripr vigan vigra vigun); 
        print "num.counts____\ttotal\tspecies\tbelow\tequal\tabove", join("\t", @species),"\n";
      ' > $OUTFILE

    for path in 55_hmmalign_trim2/*; do 
      export file=`basename $path`; 
      cat $path | awk -v ORS=" " '$1~/>/ {print substr($1,2,5)}' | 
        perl -ane '
          BEGIN{
            @species = qw(aradu arahy araip cajca cicar glyma lotja lupan medtr phavu tripr vigan vigra vigun); 
            $FILE=$ENV{"file"}
          }
          my $out_line;
          my %ct_sp;
          my ($ct_present, $ct_below, $ct_eq, $ct_above) = (0,0,0,0);
          for my $sp (@species){
            $ct_sp{$sp} = 0;
            for my $seq (@F){ if ($seq =~ /$sp/) {$ct_sp{$sp}++ } }
            $out_line .= "$ct_sp{$sp}\t";
            if ($ct_sp{$sp}>0){$ct_present++}
            if ($ct_sp{$sp}<2){$ct_below++}
            if ($ct_sp{$sp}==2){$ct_eq++}
            if ($ct_sp{$sp}>2){$ct_above++}
          }
          $gene_ct = 0;
          foreach $sp (keys %ct_sp){ $gene_ct+=$ct_sp{$sp}}
          $fam_num++;
          my $fam_id_and_stats = "$FILE\t$gene_ct\t$ct_present\t$ct_below\t$ct_eq\t$ct_above";
          print "$fam_id_and_stats\t$out_line\n";
        ' # end of perl script
    done >> $OUTFILE

  # get histogram
    cat $OUTFILE | cut -f2 | awk 'NR>1' | histogram -n -s1 | 
      histplot -d 20 | sed 's/\.00//' | head -75 > $OUTFILE.hist

  # Check the proportion of genes in families, of all included legume species
    cat 55_hmmalign_trim2/* | grep -c '>'
      493967
    cat ../01_proteomes/* | grep -c '>'
      539717
    perl -le 'print 100*493967/539717'
      91.52  # Combined attrition proportion is 8.5% (i.e. 8.5% of genes not finding a family) - which seems OK.


# Generate AHRD annotations for consensus sequences
  # Get files of consensus sequences. Run this in chunks, because AHRD is a memory hog.
    mkdir 57_hmmemit_chunks
    for part in 0 1 2 3 4 5 6 7 8 9 B C D F G H J K L M N P Q R S T V W X Y Z ; do
      cat 54_hmmemit/L_$part* > 57_hmmemit_chunks/L_$part.fa 
    done

  # Run Will Nelson's standalone annotation pipeline (as modified by me). 
  # The executable and databases are at /home/scannon/standalone_annot/


# Files for loading at legumeinfo
  
#####
# For comparison, place proteomes into Phytozome12 families
# Get pz12 families: cluster_HMMs_Fabidae_5205.hmm

  cd /scratch/scannon/phytozome12

  # press hmm database
  mkdir hmmdb
  mv cluster_HMMs_Fabidae_5205.hmm hmmdb/
  hmmpress hmmdb/cluster_HMMs_Fabidae_5205.hmm

  
  mkdir 60_hmmscan_pz12

  HMMDB=hmmdb/cluster_HMMs_Fabidae_5205.hmm
  FASDIR=../lgf5/01_proteomes
  TBLOUTDIR=60_hmmscan_pz12
  JOBMAX=20

  nohup ./commands/do_hmmscan_of_proteomes_vs_hmm.ksh $HMMDB $FASDIR $TBLOUTDIR $JOBMAX &


  # tracking output for each output file, while the hmmscan is running:
    for file in *fa; do 
      FILE=$file; 
      COUNT=`top_hmmscan_line.awk $file | wc -l`; 
      echo "$FILE $COUNT"; 
    done


# For each family, calculate the hmmscan median and median-absolute-deviation (MAD) scores
  cat 60_hmmscan_pz12/* | top_hmmscan_line.awk | 
    awk '$1!~/^#/ {print $1 "\t" $6}' | sort -k1,1 -k2nr,2nr | 
    median_abs_dev.pl -col 2 | sort -k1,1 > tab.60_MAD

# Take the top hmmscan match
  cat 60_hmmscan_pz12/* | grep -v "#" | top_hmmscan_line.awk > tab.60_top_hmmscan

# Filter by percent of median for the family

  # Go with this one:
  filter_by_median_abs_dev.pl -hmm tab.60_top_hmmscan -med tab.60_MAD -floor 40 -no_MAD |
         sort -k1,1 -k10nr,10nr  > tab.60_top_hmmscan_filt_median_fl40 &

# Check proportions matching families at this cutoff, by species
  awk '$9!~/verdict/ {print $2 "\t" $9}' tab.60_top_hmmscan_filt_median_fl40 | 
    perl -pe 's/(^\w\w\w\w\w)\.\S+\t(\S+)/$1\t$2/' | sort | uniq -c | 
    awk -v OFS="\t" '$3~/doesnt/ {no=$1} $3~/MATCHES/ {yes=$1; sp=$2; print sp, no, yes, no/yes}'
      aradu 5480  26656 0.205582
      arahy 9642  52359 0.184152
      araip 7377  28593 0.258
      cajca 3503  35237 0.0994125
      cicar 1822  24189 0.0753235
      glyma 3551  49011 0.0724531
      lotja 7282  26164 0.278321
      lupan 1704  30003 0.0567943
      medtr 5481  39106 0.140158
      phavu  743  25815 0.0287817
      tripr 6296  31377 0.200657
      vigan 2280  22642 0.100698
      vigra 2015  19896 0.101277
      vigun  807  26940 0.0299555

# Compare hmmsearch matches for lgf5 and the phytozome fabid models
  wc -l tab.60_top_hmmscan  # Phytozome12 Fabid
    495971
  wc -l tab.50_top_hmmscan  # lgf5
    609628
  perl -le 'print 100*495971/609628'
    81.35  # The phytozome families capture only 81% as many sequences as the lgf5 families

# Suppress sequences not matching -- i.e. $9~/doesnt/ -- and cases where there are fewer than 3 matching legume sequences.
  
  # get count of legumes per family 
  cat tab.60_top_hmmscan_filt_median_fl40 | 
    awk '$9~/MATCHES/ && $2~/aradu|arahy|araip|cajca|cicar|glyma|lotja|lupan|medtr|phavu|tripr|vigan|vigra|vigun/ {print $1}' | 
    uniq -c | awk '{print $2 "\t" $1}' | sort > tmp.ct_by_species.top_hmmscan_filt_median_fl40

  # Join the family count to the HMM match lines, then exclude outgroups in the output (add the top outgroup later)
  # and filter by number of legume sequences in the family. Unlike for lgf5, DON'T suppress the smallest families (<3) 
  sort tab.60_top_hmmscan_filt_median_fl40 | join -t $'\t' - tmp.ct_by_species.top_hmmscan_filt_median_fl40 | 
    awk '$9~/MATCHES/' > tmp.ct_by_species.top_hmmscan_filt_median_fl40.filt.leg
  
  # Reshape output into a family file with format like that generated by mcl:
  # each line begins with a family ID, and the gene list follows, space-separated.
  # then use get_fasta_from_family_file.pl to retrieve sequences for each family.

  cat tmp.ct_by_species.top_hmmscan_filt_median_fl40.filt.leg |
    sort -k1,1 -k3nr,3nr | cut -f1,2 | 
    hash_to_rows_by_1st_col.awk > tab.50_top_hmmscan_filt_median_fl40.fam

    # MAJOR PRODUCT: gene family file, ready for getting fasta sequences
    # tab.50_top_hmmscan_filt_median_fl40.fam
  # This phytozome file has 26131 families - in contrast with 18549 in lgf5


  TO DO: find a way to compare gene memberships between the phytozome and lgf5 families
# Suppress sequences not matching -- i.e. $9~/doesnt/ -- and cases where there are fewer than 3 matching legume sequences.
  
  # get count of legumes per family 
  cat tab.60_top_hmmscan_filt_median_fl40 | 
    awk '$9~/MATCHES/ && $2~/aradu|arahy|araip|cajca|cicar|glyma|lotja|lupan|medtr|phavu|tripr|vigan|vigra|vigun/ {print $1}' | 
    uniq -c | awk '{print $2 "\t" $1}' | sort > tmp.ct_by_species.top_hmmscan_filt_median_fl40

#####
  # Compare these from lgf5:     tab.50_top_hmmscan_filt_median_fl40_lgf5_legumeonly
  # with these from phytozome12: tab.60_top_hmmscan_filt_median_fl40
  # Strategy: for each family set, print output in the format
  #   gene       family
  # then sort and join on families

    cat tab.50_top_hmmscan_filt_median_fl40_lgf5_legumeonly |
      awk '$9~/MATCHES/ {print $2 "\t" $1}' | LC_COLLATE=C sort > tab.lgf5_gene_and_family &
    cat tab.60_top_hmmscan_filt_median_fl40 |
      awk '$9~/MATCHES/ {print $2 "\t" $1}' | LC_COLLATE=C sort > tab.pz12_gene_and_family &

    join tab.lgf5_gene_and_family tab.pz12_gene_and_family | awk '{print $2 "\t" $3}' | 
      sort | uniq -c | awk -v OFS="\t" '{print $1, $2, $3}' > tab.lgf5_pz12_corr

    join tab.pz12_gene_and_family tab.lgf5_gene_and_family | awk '{print $2 "\t" $3}' | 
      sort | uniq -c | awk -v OFS="\t" '{print $1, $2, $3}' > tab.pz12_lgf5_corr
  
# Calculate some basic stats
  cut -f2 tab.lgf5_gene_and_family | sort | uniq -c | awk '{print $1}' | histogram | sort -n
    # largest five: 811, 959, 974, 1112, 1260
    # mode at 16 (peak with counts>1000 spanning 15,16,17)
    # bins with counts > 100 span 2 - 38
    # there are 950 families of size 1-4 ("small")
    # there are 20 families of size >= 500 ("large")
  cut -f2 tab.pz12_gene_and_family | sort | uniq -c | awk '{print $1}' | histogram | sort -n
    # largest six: 1017, 1059, 1076, 1177, 1214, 1447
    # mode at 16 (peak with counts>1000 spanning 14,15,16,17,18)
    # bins with counts > 100 span 1-36
    # there are 5401 families of size 1-4 ("small")
    # there are 13131313131313131313131313 families of size >= 500 ("large")

