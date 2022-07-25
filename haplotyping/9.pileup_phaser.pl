use strict;
use warnings;
use POSIX;

# REQUIREMENTS
my $basedir  = '/t1-data/project/pregcare/sbush';
my $samtools = '/t1-data/project/GorielyLab2021/sbush/programs/samtools-1.12/bin/samtools';
my $tabix    = '/t1-data/project/GorielyLab2021/sbush/programs/htslib-1.12/tabix';
my $dbsnp    = '/t1-data/project/pregcare/sbush/dbSNP/dbSNP_153.hg38.vcf.gz';
my $coords   = "$basedir/ont_coords.tsv";
my $fatal    = 0;
if (!(-e($dbsnp)))    { $fatal++; print "ERROR: cannot find $dbsnp\n";    }
if (!(-e($tabix)))    { $fatal++; print "ERROR: cannot find $tabix\n";    }
if (!(-e($coords)))   { $fatal++; print "ERROR: cannot find $coords\n";   }
if (!(-e($samtools))) { $fatal++; print "ERROR: cannot find $samtools\n"; }
exit 1 if ($fatal > 0);

# PARAMETERS
my $min_depth 	     = 30;
my @families_to_run  = (qw/FAM44/);
my %families_to_run  = map {$_ => 1} @families_to_run;
my %informative_snps = ('FAM44' => 'rs1635560'); # this script will iterate through every position in a given sequenced region to find a potential SNP for phasing the DNM. This is of course time consuming. If we already know of a potential SNP, though, we'll use that instead. This hash is used to quickly skip to the specified SNP and ignore all others. We do this should we wish to re-run the script and get results quickly.

# OUTPUT
my $out_file1 = "$basedir/2x2_tables_for_pileup_inheritance_prediction.txt";
my $out_file2 = "$basedir/inheritance.pileup.tsv";
open(OUT1,'>',$out_file1) or die $!; open(OUT2,'>',$out_file2) or die $!;
print OUT2 "Family\tRegion sequenced\tLength of region sequenced (bp)\tDe novo mutation (chr:position)\tReference allele (for DNM)\tVariant allele (for DNM)\tRead counts for DNM SNP or number of reads supporting DNM indel (for SNPs, count are given as A|C|G|T), obtained from pileup\tTotal depth of coverage (for DNM)\tdbSNP ID (for DNM)\tHow many SNPs are in the same phase set as the DNM, not including the DNM itself?\tSNPs in the proband's phase set, not including the DNM (listed in order of distance from DNM and reported as dbSNP IDs if available, else 'chr:pos')\tHow many INFORMATIVE SNPs (= different between parents) are in the same phase set, not including the DNM?\tInformative SNPs (= different between parents) in the proband's phase set, not including the DNM (listed in order of distance from DNM and reported as dbSNP IDs if available, else 'chr:pos')\tTotal depth of coverage for informative SNPs in the proband's phase set, listed in the order previously stated\tBases that phase with the DNM (in the proband, listed in the order previously stated)\tGenotypes, listed in order previously stated (mother)\tGenotypes, listed in order previously stated (father)\tGenotypes, listed in order previously stated (proband)\tPattern of inheritance\tMin. distance of an informative SNP from DNM (bp)\tMax. distance of an informative SNP from DNM (bp)\n";

# STORE DNM VARIANTS FOR EACH FAMILY AND THE BARCODES AND MinION RUNS OF EACH MOTHER/FATHER/CHILD
my %snps_per_fam = (); my %bcs_per_fam = (); my %regions_sequenced = ();
open(IN,$coords) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line); $line =~ s/\r$//g; $line =~ s/\s$//g;
	  my @line = split(/\t/,$line);
	  my $onts = $line[0]; my $bc_line = $line[1]; my $family = $line[2]; my $region_sequenced = $line[3]; my $dbsnp_id = $line[4]; my $dnm_loc = $line[5]; my $ref = $line[6]; my $alt = $line[7]; my $dnm_type = $line[8]; my $indel_in_pileup_notation = $line[9];
	  my $dnm_chr = ''; my $dnm_pos = '';
	  if ($dnm_loc =~ /^(.*?)\:(\d+)$/) { $dnm_chr = $1; $dnm_pos = $2; } else { print "ERROR: unable to parse $dnm_loc\n"; exit 1; }
	  $regions_sequenced{$family} = $region_sequenced;
	  $snps_per_fam{$family}{dnm_as_indel} = $indel_in_pileup_notation;
	  $snps_per_fam{$family}{ont} = $onts;
	  $snps_per_fam{$family}{chr} = $dnm_chr;
	  $snps_per_fam{$family}{pos} = $dnm_pos;
	  $snps_per_fam{$family}{ref} = $ref;
	  $snps_per_fam{$family}{alt} = $alt;
	  $snps_per_fam{$family}{id}  = $dbsnp_id;
	   $bcs_per_fam{$family} 	  = $bc_line;
	}
close(IN) or die $!;

my @family = ();
while((my $family,my $irrel)=each(%snps_per_fam))
	{ push(@family,$family); }
my @sorted_family = sort {$a cmp $b} @family;
foreach my $family (@sorted_family)
	{ next if (!(exists($families_to_run{$family}))); # CHECKPOINT: restrict analysis only to specific families of interest (i.e. those for which inheritance has not been resolved by Medaka already)
	  next if (!(exists($bcs_per_fam{$family})));
	  my $region_sequenced = $regions_sequenced{$family};
	  my $length_of_region_sequenced = 0;
	  if ($region_sequenced =~ /^.*?\:(\d+)\-(\d+)$/)
		{ $length_of_region_sequenced = ($2-$1)+1; }

	  # check to see if the de novo mutation is an indel. We will make use of this variable later if so.
	  my $dnm_indel = $snps_per_fam{$family}{dnm_as_indel};
	  
	  # determine which MinION run is associated with which member of the trio (and their barcode)
	  my $onts = $snps_per_fam{$family}{ont};
	  next if ($onts !~ /^(.*?)\|(.*?)\|(.*?)$/);
	  my $mother_ont = ''; my $father_ont = ''; my $child_ont = '';
	  if ($onts =~ /^(.*?)\|(.*?)\|(.*?)$/)
		{ $mother_ont = $1; $father_ont = $2; $child_ont = $3; }
	  
	  # 1. in the child BAM, which reads cover the DNM, i.e. position $snps_per_fam{$family}{pos}?
	  # (notes: see an alternative pileup parser at https://raw.githubusercontent.com/Niknafs/NGSTools/master/baseParser.py and associated discussion: https://www.biostars.org/p/76829/)
	  my $bc_line   = $bcs_per_fam{$family};
	  my @bc_line   = split(/\|/,$bc_line);
	  my $child_bc  = $bc_line[2];
	  my $child_bam = "$basedir/$child_ont/6.bam/$child_bc.bam";
	  next if (!(-e($child_bam)));
	  my $pileup = ''; my $read_ids = '';
	  open(IN,"$samtools mpileup -r $snps_per_fam{$family}{chr}:$snps_per_fam{$family}{pos}-$snps_per_fam{$family}{pos} -aa $child_bam --output-QNAME 2> /dev/null |") or die $!;
	  while(<IN>)
		{ my $line = $_; chomp($line);
		  my @line = split(/\t/,$line);
		  $pileup = $line[4]; $read_ids = $line[6];
		}
	  close(IN) or die $!;
	  my @pileup   = split(//,$pileup);			  
	  my @read_ids = split(",",$read_ids);
	
	  # iterate through @pileup and make a new pileup array, one where each entry has 1-to-1 correspondence with @read_ids. To do this, we must set to one side all sequences beginning \+[0-9]+[ACGTNacgtn]+. These refer to indels that *START FROM THE NEXT POSITION* - which means their QNAME will NOT be output HERE.
	  my @new_pileup = ();
	  my $skip = 0; my $entry_num_in_read_id_array = 0;
	  for(my $x=0;$x<@pileup;$x++)
		{ my $n = $pileup[$x];
		  my $n1; if (defined($pileup[$x+1])) { $n1 = $pileup[$x+1]; }
		  my $skip_this = 0;
		  if ($skip > 0)
			{ $skip--; $skip_this++; }
		  elsif ($n =~ /^(\d+)$/)
			{ $skip = $1; $skip_this++;
			  if ( (defined($n1)) and ($n1 =~ /^(\d+)$/) )
				{ $skip .= $1; $skip++; }
			}
		  next if ($skip_this > 0);
		  next if (($n eq '+') or ($n eq '-')); # a sequence matching the regular expression \+[0-9]+[ACGTNacgtn]+ denotes an insertion of one or more bases *STARTING FROM THE NEXT POSITION* (if starting "\-", a deletion). For example, +2AG means insertion of AG in the forward strand and -2ct means deletion of CT in the reverse strand
		  next if (($n eq '^') or ($n eq '$') or ($n eq '[') or ($n eq ']')); # ^ (caret) marks the start of a read segment. $ (dollar) marks the end of a read segment
		  push(@new_pileup,$n); # most (but not all) $n at this point will be these: AGTCN (upper case), which denotes a base that did not match the reference on the forward strand and agtcn (lower case), which denotes a base that did not match the reference on the reverse strand. That these values of $n 'do not match the reference' is because in the call to mpileup, above, we have saved a little runtime by NOT specifying the reference genome (this means that for each base, the program won't check to see if it matches the reference or not). If we had specified then reference, then we would instead see values of $n that are either . (dot), which means a base that matched the reference on the forward strand, or , (comma), which means a base that matched the reference on the reverse strand.
		  my $read_id = $read_ids[$entry_num_in_read_id_array];
		  $entry_num_in_read_id_array++;
		}
	  if ($#new_pileup != $#read_ids)
		{ print "FATAL ERROR with parsing the pileup of the child DNM (family $family): the size of the pileup line ($#pileup) and the number of read IDs ($#read_ids) does not match. pileup line: $pileup\n"; exit 1; }
	  next if ($#new_pileup != $#read_ids);
	  print "$family: @pileup\n";
	  
	  # count how many reads call the DNM as REF and the DNM as ALT
	  my %nuc_counts = (A => 0, G => 0, C => 0, T => 0);
	  my %dnm_read_ids_REF = (); my %dnm_read_ids_ALT = ();
	  my $num_reads_containing_the_dnm_del = 0; my $num_reads_calling_a_base_not_a_del = 0;
	  my $num_reads_containing_the_dnm_ins = 0; my $num_reads_calling_a_base_not_an_ins = 0;
	  if ($dnm_indel eq '.') # THE DNM IS A SNP
		{ for(my $x=0;$x<@new_pileup;$x++)
			{ my $nt 	  = uc($new_pileup[$x]);
			  my $read_id = $read_ids[$x];
			  $nuc_counts{$nt}++;
			  if    ($nt eq ($snps_per_fam{$family}{ref}))
				{ $dnm_read_ids_REF{$read_id}++; }	
			  elsif ($nt eq ($snps_per_fam{$family}{alt}))
				{ $dnm_read_ids_ALT{$read_id}++; }			  
			}
		}
	  elsif ($dnm_indel =~ /^\-(\d+)/) # THE DNM IS A DELETION
		{ my $length_of_deletion = $1;
		  my $start_of_deletion = $snps_per_fam{$family}{pos}+1;
		  my $end_of_deletion = $snps_per_fam{$family}{pos}+$length_of_deletion;
		  my %reads_with_a_del = (); my %reads_without_a_del = ();
		  for(my $pos=($start_of_deletion-1);$pos<=($end_of_deletion+1);$pos++)
			{ my $pileup = ''; my $read_ids = '';
			  open(IN,"$samtools mpileup -r $snps_per_fam{$family}{chr}:$pos-$pos -aa $child_bam --output-QNAME 2> /dev/null |") or die $!;
			  while(<IN>)
				{ my $line = $_; chomp($line);
				  my @line = split(/\t/,$line);
				  $pileup = $line[4]; $read_ids = $line[6];
				}
			  close(IN) or die $!;
			  my @pileup   = split(//,$pileup);			  
			  my @read_ids = split(",",$read_ids);
			  my $skip = 0; my $entry_num_in_read_id_array = 0;
			  for(my $x=0;$x<@pileup;$x++)
				{ my $n = $pileup[$x];
				  my $n1; if (defined($pileup[$x+1])) { $n1 = $pileup[$x+1]; }
				  my $skip_this = 0;
				  if ($skip > 0)
					{ $skip--; $skip_this++; }
				  elsif ($n =~ /^(\d+)$/)
					{ $skip = $1; $skip_this++;
					  if ( (defined($n1)) and ($n1 =~ /^(\d+)$/) )
						{ $skip .= $1; $skip++; }
					}
				  next if ($skip_this > 0);
				  next if (($n eq '+') or ($n eq '-')); # a sequence matching the regular expression \+[0-9]+[ACGTNacgtn]+ denotes an insertion of one or more bases *STARTING FROM THE NEXT POSITION* (if starting "\-", a deletion). For example, +2AG means insertion of AG in the forward strand and -2ct means deletion of CT in the reverse strand
				  next if (($n eq '^') or ($n eq '$') or ($n eq '[') or ($n eq ']')); # ^ (caret) marks the start of a read segment. $ (dollar) marks the end of a read segment
				  my $read_id = $read_ids[$entry_num_in_read_id_array];
				  if ($n eq '*') # * (asterisk) is a placeholder for a deleted base in a multiple basepair deletion that was mentioned in a previous line by the -[0-9]+[ACGTNacgtn]+ notation. IMPORTANT: the length of the deletion cannot be inferred from the presence of the *
					{ $reads_with_a_del{$read_id}{$pos}++; }
				  else
					{ $reads_without_a_del{$read_id}{$pos}++; }
				  $entry_num_in_read_id_array++;
				}
			}
		  # reads that contain the DNM deletion must, by definition, have a * in the pileup for each base of the deletion AND no * in the first base of the deletion-1 AND no * in the last base of the deletion +1
		  my %reads_containing_the_dnm_del = ();
		  while((my $read_id,my $irrel)=each(%reads_with_a_del))
			{ my $failure = 0;
			  for(my $pos=$start_of_deletion;$pos<=$end_of_deletion;$pos++)
				{ if (!(exists($reads_with_a_del{$read_id}{$pos})))
					{ $failure++; }
				}
			  my $start_minus1 = $start_of_deletion-1;
			  my $end_plus1    = $end_of_deletion+1;
			  if (exists($reads_with_a_del{$read_id}{$start_minus1})) { $failure++; }
			  if (exists($reads_with_a_del{$read_id}{$end_plus1}))    { $failure++; }
			  if ($failure == 0)
				{ $reads_containing_the_dnm_del{$read_id}++; }
			}
		  $num_reads_containing_the_dnm_del = scalar keys %reads_containing_the_dnm_del;
		  $num_reads_calling_a_base_not_a_del = scalar keys %reads_without_a_del;
		  %dnm_read_ids_ALT = %reads_containing_the_dnm_del;
		  %dnm_read_ids_REF = %reads_without_a_del;
		}
	  elsif ($dnm_indel =~ /^\+(\d+)(\w+)/) # THE DNM IS AN INSERTION
		{ my $length_of_insertion = $1;
		  my $inserted_bases 	  = $2;
		  my @inserted_bases 	  = split(//,$inserted_bases);
		  my $base_post_insertion = '';
		  if ($snps_per_fam{$family}{ref} =~ /^.+?(\w{1})$/)
			{ $base_post_insertion = $1; }
		  my $start_of_insertion  = $snps_per_fam{$family}{pos};
		  my $end_of_insertion    = $snps_per_fam{$family}{pos}+$length_of_insertion;
		  my %bases_per_read_id   = ();
		  for(my $pos=$start_of_insertion;$pos<=$end_of_insertion+1;$pos++)
			{ my $pileup = ''; my $read_ids = '';
			  open(IN,"$samtools mpileup -r $snps_per_fam{$family}{chr}:$pos-$pos -aa $child_bam --output-QNAME 2> /dev/null |") or die $!;
			  while(<IN>)
				{ my $line = $_; chomp($line);
				  my @line = split(/\t/,$line);
				  $pileup = $line[4]; $read_ids = $line[6];
				}
			  close(IN) or die $!;
			  my @pileup   = split(//,$pileup);			  
			  my @read_ids = split(",",$read_ids);
			  my $skip = 0; my $entry_num_in_read_id_array = 0;
			  for(my $x=0;$x<@pileup;$x++)
				{ my $n = $pileup[$x];
				  my $n1; if (defined($pileup[$x+1])) { $n1 = $pileup[$x+1]; }
				  my $skip_this = 0;
				  if ($skip > 0)
					{ $skip--; $skip_this++; }
				  elsif ($n =~ /^(\d+)$/)
					{ $skip = $1; $skip_this++;
					  if ( (defined($n1)) and ($n1 =~ /^(\d+)$/) )
						{ $skip .= $1; $skip++; }
					}
				  next if ($skip_this > 0);
				  next if (($n eq '+') or ($n eq '-'));
				  next if (($n eq '^') or ($n eq '$') or ($n eq '[') or ($n eq ']'));
				  my $read_id = $read_ids[$entry_num_in_read_id_array];
				  if ($n =~ /^[ATCG]$/) { $bases_per_read_id{$read_id}{$pos} = $n; }
				  $entry_num_in_read_id_array++;
				}
			}
		  # reads that contain the DNM insertion must, by definition, call all the bases of the insertion in the appropriate order. We can manually check this by using Data::Dumper to dump %bases_per_read_id, cross-referencing its content with both %reads_without_an_ins and %reads_containing_the_dnm_ins
		  my %reads_containing_the_dnm_ins = (); my %reads_without_an_ins = ();
		  while((my $read_id,my $irrel)=each(%bases_per_read_id))
			{ my $success = 0; my $num_tests = 0; my $base_num = 0;
			  for(my $pos=($start_of_insertion+1);$pos<=$end_of_insertion+1;$pos++) # in pileup notation the 'start' is the base BEFORE the inserted base
				{ if (!(defined($inserted_bases[$base_num])))
					{ print "FATAL ERROR with parsing the pileup of the child DNM (family $family): unable to synchronize the length of the insertion as given by its coords ($start_of_insertion-$end_of_insertion) with the number of bases it is said to contain ($inserted_bases)\n"; exit 1; }
			      my $inserted_base = $inserted_bases[$base_num];
				  if (exists($bases_per_read_id{$read_id}{$pos}))
					{ my $base_in_this_read_at_this_pos = $bases_per_read_id{$read_id}{$pos};
					  if ($pos <= $end_of_insertion)
						{ if ($base_in_this_read_at_this_pos eq $inserted_base) # we expect the read to contain every inserted base, in that order
							{ $success++; }
						}
					  elsif ($pos == ($end_of_insertion+1))
						{ if ($base_in_this_read_at_this_pos eq $base_post_insertion) # we also expect to contain the appropriate call for the first base after the insertion
							{ $success++; }
						}
					}
				  $num_tests++;
				}
			  $base_num++;
			  if ($success == $num_tests)
				{ $reads_containing_the_dnm_ins{$read_id}++; }
			  else
				{ $reads_without_an_ins{$read_id}++; }
			}
		  $num_reads_containing_the_dnm_ins    = scalar keys %reads_containing_the_dnm_ins;
		  $num_reads_calling_a_base_not_an_ins = scalar keys %reads_without_an_ins;
		  %dnm_read_ids_ALT = %reads_containing_the_dnm_ins;
		  %dnm_read_ids_REF = %reads_without_an_ins;
		}

	  # create summary counts of the number and % of reads supporting the DNM REF/ALT in the child
	  my $dnm_A = '.'; my $dnm_T = '.'; my $dnm_C = '.'; my $dnm_G = '.'; my $dnm_total = '.';
	  if ($dnm_indel eq '.')
		{ $dnm_A  = $nuc_counts{'A'};
		  $dnm_T  = $nuc_counts{'T'};
		  $dnm_C  = $nuc_counts{'C'};
		  $dnm_G  = $nuc_counts{'G'};
		  $dnm_total = $dnm_A+$dnm_T+$dnm_C+$dnm_G;
		}
	  elsif ($dnm_indel =~ /^\-\d+/)
		{ $dnm_total = $num_reads_containing_the_dnm_del+$num_reads_calling_a_base_not_a_del; }
	  elsif ($dnm_indel =~ /^\+\d+\w+$/)
		{ $dnm_total = $num_reads_containing_the_dnm_ins+$num_reads_calling_a_base_not_an_ins; }
	  my $num_of_reads_calling_DNM_REF_in_child = scalar keys %dnm_read_ids_REF;
	  my $num_of_reads_calling_DNM_ALT_in_child = scalar keys %dnm_read_ids_ALT;
	  my $pc_of_reads_calling_DNM_REF_in_child  = 0; if ($dnm_total > 0) { $pc_of_reads_calling_DNM_REF_in_child = sprintf("%.2f",(($num_of_reads_calling_DNM_REF_in_child/$dnm_total)*100)); }
	  my $pc_of_reads_calling_DNM_ALT_in_child  = 0; if ($dnm_total > 0) { $pc_of_reads_calling_DNM_ALT_in_child = sprintf("%.2f",(($num_of_reads_calling_DNM_ALT_in_child/$dnm_total)*100)); }
	  
	  # 2. extract each REF- and ALT-containing read from the child BAM. We do this because want to know, for the ALT-containing reads *specifically*, which genomic region this entire set of alignments covers. This is because we will later be looking for informative (=we can phase using them) SNPs within it.
	  open (TMP,'>','TEMP_read_ids_REF.txt') or die $!;
	  while((my $read_id,my $irrel)=each(%dnm_read_ids_REF))
		{ print TMP "$read_id\n"; } # REMINDER: these are the REF-containing reads IN THE CHILD. They all contain base $snps_per_fam{$family}{ref} at position $snps_per_fam{$family}{pos}
	  close(TMP) or die $!;
	  my %dnm_read_ids_REF_with_coords = ();
	  open(IN,"$samtools view -N TEMP_read_ids_REF.txt $child_bam 2> /dev/null |") or die $!; # you can't just open this "samtools view -N <(echo $read_id) $child_bam |" because you're not allowed to open a command that both pipes in AND out
	  while(<IN>)
		{ my $line = $_; chomp($line);
		  my @line = split(/\t/,$line);
		  my $id = $line[0]; my $chr = $line[2]; my $start = $line[3]; my $cigar = $line[5]; my $seq = $line[9];
		  next if ($chr !~ /^$snps_per_fam{$family}{chr}$/); # CHECKPOINT: discard data not associated with the chromosome on which the DNM is found. This should never be triggered, but it's here just to make sure.
		  # $start is the start coordinate for the alignment, but how do we get the end coordinate? We parse the CIGAR string, like so: https://www.biostars.org/p/1680/
		  # we increment $end (i.e. the position of the alignment on the reference) when any of the 'consume reference' codes are raised. See page 8 of https://samtools.github.io/hts-specs/SAMv1.pdf
		  my $end = $start;
		  while ($cigar =~ /(\d+)(\w{1})/g)
			{ my $num = $1; my $type = $2;
			  if (($type eq 'M') or ($type eq 'D') or ($type eq 'N') or ($type eq '=') or ($type eq 'X'))
				{ $end += $num; }
			}
		  $dnm_read_ids_REF_with_coords{$id}{start} = $start;
		  $dnm_read_ids_REF_with_coords{$id}{end}   = $end;
		}
	  close(IN) or die $!;
	  unlink 'TEMP_read_ids_REF.txt' or die $!;
	  
	  open (TMP,'>','TEMP_read_ids_ALT.txt') or die $!;
	  while((my $read_id,my $irrel)=each(%dnm_read_ids_ALT))
		{ print TMP "$read_id\n"; } # REMINDER: these are the ALT-containing reads IN THE CHILD. They all contain base $snps_per_fam{$family}{alt} at position $snps_per_fam{$family}{pos}
	  close(TMP) or die $!;
	  my @map_coords = (); my %dnm_read_ids_ALT_with_coords = ();
	  open(IN,"$samtools view -N TEMP_read_ids_ALT.txt $child_bam 2> /dev/null |") or die $!; # you can't just open this "samtools view -N <(echo $read_id) $child_bam |" because you're not allowed to open a command that both pipes in AND out
	  while(<IN>)
		{ my $line = $_; chomp($line);
		  my @line = split(/\t/,$line);
		  my $id = $line[0]; my $chr = $line[2]; my $start = $line[3]; my $cigar = $line[5]; my $seq = $line[9];
		  next if ($chr !~ /^$snps_per_fam{$family}{chr}$/); # CHECKPOINT: discard data not associated with the chromosome on which the DNM is found. This should never be triggered, but it's here just to make sure.
		  # $start is the start coordinate for the alignment, but how do we get the end coordinate? We parse the CIGAR string, like so: https://www.biostars.org/p/1680/
		  # we increment $end (i.e. the position of the alignment on the reference) when any of the 'consume reference' codes are raised. See page 8 of https://samtools.github.io/hts-specs/SAMv1.pdf
		  my $end = $start;
		  while ($cigar =~ /(\d+)(\w{1})/g)
			{ my $num = $1; my $type = $2;
			  if (($type eq 'M') or ($type eq 'D') or ($type eq 'N') or ($type eq '=') or ($type eq 'X'))
				{ $end += $num; }
			}
		  #print "$id\t$chr\t$start\t$end\t$seq\n"; # print this line if you want to see what the sequence of each ALT-containing read in the child is
		  push(@map_coords,$start,$end);
		  $dnm_read_ids_ALT_with_coords{$id}{start} = $start;
		  $dnm_read_ids_ALT_with_coords{$id}{end}   = $end;
		}
	  close(IN) or die $!;
	  unlink 'TEMP_read_ids_ALT.txt' or die $!;
	  
	  # 3. now we need to work out which SNPs - in either the mother or father - have been called within the region covered by the proband's DNM ALT-containing reads
	  # we do this by parsing dbSNP to identify all biallelic SNPs within the region
	  my @sorted_coords = sort {$a <=> $b} @map_coords;
	  my $min_coord = '.'; my $max_coord = '.'; my $length = '.';
	  my %possible_phasing_snps = ();
	  if ($#sorted_coords != -1)
		{ $min_coord = $sorted_coords[0]; $max_coord = $sorted_coords[$#sorted_coords];
		  $length = ($max_coord-$min_coord)+1;
		  print OUT1 "*****\n";
		  print OUT1 "$family:\nthe DNM occurs at $snps_per_fam{$family}{chr}:$snps_per_fam{$family}{pos} and is $snps_per_fam{$family}{ref}/$snps_per_fam{$family}{alt}\nthere are $num_of_reads_calling_DNM_ALT_in_child reads which call the DNM ALT, spanning $snps_per_fam{$family}{chr}:$min_coord-$max_coord ($length bp)\n";
		  open(IN,"$tabix $dbsnp $snps_per_fam{$family}{chr}:$min_coord-$max_coord |") or die $!;
		  while(<IN>)
			{ my $line = $_; chomp($line);
			  my @line = split(/\t/,$line);
			  my $pos = $line[1]; my $dbsnp_id = $line[2]; my $ref = $line[3]; my $alt = $line[4];
			  my $len_ref = length($ref); my $len_alt = length($alt);
			  next if  ($len_ref != 1);
			  next if (($len_alt > 1) and ($alt !~ /\,/)); # this is an indel, not a multiallelic SNP
			  next if  ($ref eq $alt);
			  $possible_phasing_snps{$dbsnp_id}{pos} = $pos;
			  $possible_phasing_snps{$dbsnp_id}{ref} = $ref;
			  $possible_phasing_snps{$dbsnp_id}{alt} = $alt;
			}
		  close(IN) or die $!;
		}
	  
	  # 4. we then need to sort this list of possible candidate SNPs by order of distance from the DNM
	  # reminder: the DNM ALT is $snps_per_fam{$family}{alt} and there are $num_of_reads_calling_DNM_ALT_in_child ($pc_of_reads_calling_DNM_ALT_in_child%) reads supporting this
	  # we are trying to find an informative (= heterozygous) SNP within the parent, one allele of which is disproportionately present on the set of $num_of_reads_calling_DNM_ALT_in_child reads. This is a clear indication that derive from the same haplotype, i.e. are in phase.
	  my @dist_to_dnm = ();
	  while((my $dbsnp_id,my $irrel)=each(%possible_phasing_snps))
		{ my $pos = $possible_phasing_snps{$dbsnp_id}{pos};
		  my $dist = 'NA';
		  if ($pos >= $snps_per_fam{$family}{pos})
			{ $dist = $pos-$snps_per_fam{$family}{pos}; }
		  else
			{ $dist = $snps_per_fam{$family}{pos}-$pos; }
		  if ($dist > 0)
			{ push(@dist_to_dnm,[$dist,$dbsnp_id]); }
		}
	  my $printed_something = 0;
	  my @sorted_dist_to_dnm = map { $_->[0] } sort { $a->[1] <=> $b->[1] } map { [$_, $_->[0]] } @dist_to_dnm;
	  
	  # we're going to iterate through the entire array to find a potential candidate SNP. If we already know of one, though, let's just use that instead.
	  my $preselected_dbsnp_id = '';		  
	  if (exists($informative_snps{$family}))
		{ $preselected_dbsnp_id = $informative_snps{$family}; }
	  
	  for(my $x=0;$x<@sorted_dist_to_dnm;$x++)
		{ my $dist = $sorted_dist_to_dnm[$x][0]; my $dbsnp_id = $sorted_dist_to_dnm[$x][1];
		  last if  ($preselected_dbsnp_id eq 'none'); # we have already run this script on this family before and we know there is no informative SNP to be found - so if we run this script again, let's not waste trying to look for one and just skip out of this loop at the first opportunity
		  next if (($preselected_dbsnp_id ne '') and ($preselected_dbsnp_id =~ /^rs/) and ($dbsnp_id ne $preselected_dbsnp_id));
		  my $pos  = $possible_phasing_snps{$dbsnp_id}{pos};
		  my $ref  = $possible_phasing_snps{$dbsnp_id}{ref};
		  my $alt  = $possible_phasing_snps{$dbsnp_id}{alt};
		  my @alt  = split(",",$alt); # some SNPs in dbSNP, e.g., rs3733875 are multiallelic: G > A,C,T. If this is the case, we first need to confirm which allele we are looking at.
		  my %alt  = map {$_ => 1} @alt;
		  print OUT1 "there is a candidate phasing SNP $dist bases from the DNM. it is $dbsnp_id, a $ref/$alt\n";
		  
		  # 5. for each candidate SNP, we look up what call has been made in both the mother and father. We require our of candidate SNP that it is a heterozygote in at least one of the parents.
		  my $genotype_mother = ''; my $genotype_father = ''; my $genotype_child = '';
		  for(my $x=0;$x<=1;$x++)
			{ my $parent_bc = $bc_line[$x]; my $person = ''; my $ont = '';
			  if 	($x == 0) { $person = 'mother'; $ont = $mother_ont; }
			  elsif ($x == 1) { $person = 'father'; $ont = $father_ont; }
			  my $parent_bam = "$basedir/$ont/6.bam/$parent_bc.bam";
			  if (!(-e($parent_bam))) { print "ERROR: unable to find the BAM for the $person of $family at $parent_bam\n"; exit 1; }
			  my $pileup = ''; my $read_ids = '';
			  open(IN,"$samtools mpileup -r $snps_per_fam{$family}{chr}:$pos-$pos -aa $parent_bam --output-QNAME 2> /dev/null |") or die $!;
			  while(<IN>)
				{ my $line = $_; chomp($line);
				  my @line = split(/\t/,$line);
				  $pileup = $line[4]; $read_ids = $line[6];
				}
			  close(IN) or die $!;
			  my @pileup   = split(//,$pileup);
			  my @read_ids = split(",",$read_ids);
			  
			  # iterate through the parent's @pileup and make a new pileup array, one where each entry has 1-to-1 correspondence with @read_ids. To do this, we must set to one side all sequences beginning \+[0-9]+[ACGTNacgtn]+. These refer to indels that *START FROM THE NEXT POSITION* - which means their QNAME will NOT be output HERE.
			  my @new_pileup = ();
			  my $skip = 0; my $entry_num_in_read_id_array = 0;
			  for(my $x=0;$x<@pileup;$x++)
				{ my $n = $pileup[$x];
				  my $n1; if (defined($pileup[$x+1])) { $n1 = $pileup[$x+1]; }
				  my $skip_this = 0;
				  if ($skip > 0)
					{ $skip--; $skip_this++; }
				  elsif ($n =~ /^(\d+)$/)
					{ $skip = $1; $skip_this++;
					  if ( (defined($n1)) and ($n1 =~ /^(\d+)$/) )
						{ $skip .= $1; $skip++; }
					}
				  next if ($skip_this > 0);
				  next if (($n eq '+') or ($n eq '-')); # a sequence matching the regular expression \+[0-9]+[ACGTNacgtn]+ denotes an insertion of one or more bases *STARTING FROM THE NEXT POSITION* (if starting "\-", a deletion). For example, +2AG means insertion of AG in the forward strand and -2ct means deletion of CT in the reverse strand
				  next if (($n eq '^') or ($n eq '$') or ($n eq '[') or ($n eq ']')); # ^ (caret) marks the start of a read segment. $ (dollar) marks the end of a read segment
				  push(@new_pileup,$n); # most (but not all) $n at this point will be these: AGTCN (upper case), which denotes a base that did not match the reference on the forward strand and agtcn (lower case), which denotes a base that did not match the reference on the reverse strand. That these values of $n 'do not match the reference' is because in the call to mpileup, above, we have saved a little runtime by NOT specifying the reference genome (this means that for each base, the program won't check to see if it matches the reference or not). If we had specified then reference, then we would instead see values of $n that are either . (dot), which means a base that matched the reference on the forward strand, or , (comma), which means a base that matched the reference on the reverse strand.
				  my $read_id = $read_ids[$entry_num_in_read_id_array];
				  $entry_num_in_read_id_array++;
				}
			  if ($#new_pileup != $#read_ids) { print "ERROR with parsing the pileup of the $person at the phasingSNP locus: the size of the pileup line ($#pileup) and the number of read IDs ($#read_ids) does not match; discarding this SNP...\n"; }
			  next if ($#new_pileup != $#read_ids);
		 
			  # count how many of the parent's reads call the phasingSNP as REF and the phasingSNP as ALT
			  my %nuc_counts = (A => 0, G => 0, C => 0, T => 0);
			  my %parent_read_ids_REF = (); my %parent_read_ids_ALT_per_nt = ();
			  for(my $x=0;$x<@new_pileup;$x++)
				{ my $nt 	  = uc($new_pileup[$x]);
				  my $read_id = $read_ids[$x];
				  $nuc_counts{$nt}++;
				  if    ($nt eq $ref)
					{ $parent_read_ids_REF{$read_id}++; }
				  elsif (exists($alt{$nt}))
					{ $parent_read_ids_ALT_per_nt{$nt}{$read_id}++; }
				}
			  
			  # if dbSNP shows this is a multiallelic SNP, we will oblige it to be 'be biallelic' by considering only the nt with the greatest number of reads assigned to it
			  my @reads_per_nt = ();
			  while((my $nt,my $irrel)=each(%parent_read_ids_ALT_per_nt))
				{ my $no_reads = scalar keys %{$parent_read_ids_ALT_per_nt{$nt}};
				  push(@reads_per_nt,[$no_reads,$nt]);
				}
			  my %parent_read_ids_ALT = ();
			  if ($#reads_per_nt == -1)
				{ $alt = $alt[0]; }
			  else
				{ my @sorted_reads_per_nt = map { $_->[0] } sort { $b->[1] <=> $a->[1] } map { [$_, $_->[0]] } @reads_per_nt;
				  $alt = $sorted_reads_per_nt[0][1];
			      %parent_read_ids_ALT = %{$parent_read_ids_ALT_per_nt{$alt}};
				}
			  
			  my $parent_A  = $nuc_counts{'A'};
			  my $parent_T  = $nuc_counts{'T'};
			  my $parent_C  = $nuc_counts{'C'};
			  my $parent_G  = $nuc_counts{'G'};
			  my $parent_total = $parent_A+$parent_T+$parent_C+$parent_G;
			  my $num_of_reads_calling_INFOSNP_REF_in_parent = scalar keys %parent_read_ids_REF;
			  my $num_of_reads_calling_INFOSNP_ALT_in_parent = scalar keys %parent_read_ids_ALT;
			  my $pc_of_reads_calling_INFOSNP_REF_in_parent  = 0; if ($parent_total > 0) { $pc_of_reads_calling_INFOSNP_REF_in_parent = sprintf("%.2f",(($num_of_reads_calling_INFOSNP_REF_in_parent/$parent_total)*100)); }
			  my $pc_of_reads_calling_INFOSNP_ALT_in_parent  = 0; if ($parent_total > 0) { $pc_of_reads_calling_INFOSNP_ALT_in_parent = sprintf("%.2f",(($num_of_reads_calling_INFOSNP_ALT_in_parent/$parent_total)*100)); }
			  print OUT1 "in the $person, this position has ACGT counts of $parent_A|$parent_C|$parent_G|$parent_T. The ALT base, $alt, has a MAF of $pc_of_reads_calling_INFOSNP_ALT_in_parent%\n";
			  my $genotype = "$ref"."$ref"; # default
			  my $genotype_if_on_x = $ref;
			  if ($parent_total > $min_depth)
				{ if (($pc_of_reads_calling_INFOSNP_ALT_in_parent > 30) and ($pc_of_reads_calling_INFOSNP_ALT_in_parent < 70))
					{ $genotype = "$ref"."$alt";
					  $genotype_if_on_x = $alt;
					}
				  elsif ($pc_of_reads_calling_INFOSNP_ALT_in_parent >= 70)
					{ $genotype = "$alt"."$alt";
					  $genotype_if_on_x = $alt;
					}
				}
			  else
				{ print OUT1 "the min. coverage depth of $min_depth was not met so we revert to a default position of assuming this genotype is homozygous for the reference allele\n"; }
			  if    ($person eq 'mother') { $genotype_mother = $genotype; }
			  elsif ($person eq 'father')
				{ if ($snps_per_fam{$family}{chr} =~ /^X$/)
					{ $genotype_father = $genotype_if_on_x; }
				  else
					{ $genotype_father = $genotype; }
				}
			}
		  if ($genotype_mother eq $genotype_father) { print OUT1 "genotype identical in both parents: discarding this SNP as a candidate\n"; }
		  next if  ($genotype_mother eq $genotype_father); # if the parents have the same genotype, the SNP is of no value for phasing
		  next if (($genotype_mother eq '') or ($genotype_father eq '')); # if no genotype can be called, the SNP is of no value for phasing
		  $alt = $possible_phasing_snps{$dbsnp_id}{alt};
		  
		  # 6. now that we have a candidate, what is the child's genotype for this position?
		  # we'll look up which reads cover this locus, i.e. position $pos, in the child. We want to know whether they contain the same call as either $mother_call or $father_call
		  print OUT1 "genotype of mother: $genotype_mother\n";
		  print OUT1 "genotype of father: $genotype_father\n";
		  my $pileup = ''; my $read_ids = '';
		  open(IN,"$samtools mpileup -r $snps_per_fam{$family}{chr}:$pos-$pos -aa $child_bam --output-QNAME 2> /dev/null |") or die $!;
		  while(<IN>)
			{ my $line = $_; chomp($line);
			  my @line = split(/\t/,$line);
			  $pileup = $line[4]; $read_ids = $line[6];
			}
		  close(IN) or die $!;
		  my @pileup   = split(//,$pileup);			  
		  my @read_ids = split(",",$read_ids);
		  
		  # iterate through the child's @pileup - for the phasingSNP locus - and make a new pileup array, one where each entry has 1-to-1 correspondence with @read_ids. To do this, we must set to one side all sequences beginning \+[0-9]+[ACGTNacgtn]+. These refer to indels that *START FROM THE NEXT POSITION* - which means their QNAME will NOT be output HERE.
		  my @new_pileup = ();
		  my $skip = 0; my $entry_num_in_read_id_array = 0;
		  for(my $x=0;$x<@pileup;$x++)
			{ my $n = $pileup[$x];
			  my $n1; if (defined($pileup[$x+1])) { $n1 = $pileup[$x+1]; }
			  my $skip_this = 0;
			  if ($skip > 0)
				{ $skip--; $skip_this++; }
			  elsif ($n =~ /^(\d+)$/)
				{ $skip = $1; $skip_this++;
				  if ( (defined($n1)) and ($n1 =~ /^(\d+)$/) )
					{ $skip .= $1; $skip++; }
				}
			  next if ($skip_this > 0);
			  next if (($n eq '+') or ($n eq '-')); # a sequence matching the regular expression \+[0-9]+[ACGTNacgtn]+ denotes an insertion of one or more bases *STARTING FROM THE NEXT POSITION* (if starting "\-", a deletion). For example, +2AG means insertion of AG in the forward strand and -2ct means deletion of CT in the reverse strand
			  next if (($n eq '^') or ($n eq '$') or ($n eq '[') or ($n eq ']')); # ^ (caret) marks the start of a read segment. $ (dollar) marks the end of a read segment
			  push(@new_pileup,$n); # most (but not all) $n at this point will be these: AGTCN (upper case), which denotes a base that did not match the reference on the forward strand and agtcn (lower case), which denotes a base that did not match the reference on the reverse strand. That these values of $n 'do not match the reference' is because in the call to mpileup, above, we have saved a little runtime by NOT specifying the reference genome (this means that for each base, the program won't check to see if it matches the reference or not). If we had specified then reference, then we would instead see values of $n that are either . (dot), which means a base that matched the reference on the forward strand, or , (comma), which means a base that matched the reference on the reverse strand.
			  my $read_id = $read_ids[$entry_num_in_read_id_array];
			  $entry_num_in_read_id_array++;
			}
		  if ($#new_pileup != $#read_ids) { print "ERROR with parsing the pileup of the child at the phasingSNP locus: the size of the pileup line ($#pileup) and the number of read IDs ($#read_ids) does not match; discarding this SNP...\n"; }
		  next if ($#new_pileup != $#read_ids);
		  
		  # count how many of the child's reads call the phasingSNP as REF and the phasingSNP as ALT
		  my %nuc_counts = (A => 0, G => 0, C => 0, T => 0);
		  my %informative_pos_read_ids_REF = (); my %informative_pos_read_ids_ALT_per_nt = ();
		  for(my $x=0;$x<@new_pileup;$x++)
			{ my $nt 	  = uc($new_pileup[$x]);
			  my $read_id = $read_ids[$x];
			  $nuc_counts{$nt}++;
			  if    ($nt eq $ref)
				{ $informative_pos_read_ids_REF{$read_id}++; }
			  elsif (exists($alt{$nt}))
				{ $informative_pos_read_ids_ALT_per_nt{$nt}{$read_id}++; }
			}
			
		  # if dbSNP shows this is a multiallelic SNP, we will oblige it to be 'be biallelic' by considering only the nt with the greatest number of reads assigned to it
		  my @reads_per_nt = ();
		  while((my $nt,my $irrel)=each(%informative_pos_read_ids_ALT_per_nt))
			{ my $no_reads = scalar keys %{$informative_pos_read_ids_ALT_per_nt{$nt}};
			  push(@reads_per_nt,[$no_reads,$nt]);
			}
		  my %informative_pos_read_ids_ALT = ();
		  if ($#reads_per_nt == -1)
			{ $alt = $alt[0]; }
		  else
			{ my @sorted_reads_per_nt = map { $_->[0] } sort { $b->[1] <=> $a->[1] } map { [$_, $_->[0]] } @reads_per_nt;
			  $alt = $sorted_reads_per_nt[0][1];
		      %informative_pos_read_ids_ALT = %{$informative_pos_read_ids_ALT_per_nt{$alt}};
			}
		  
		  my $num_A    = $nuc_counts{'A'};
		  my $num_T    = $nuc_counts{'T'};
		  my $num_C    = $nuc_counts{'C'};
		  my $num_G    = $nuc_counts{'G'};
		  my $num_total = $num_A+$num_T+$num_C+$num_G;
		  my $num_of_reads_calling_INFOSNP_REF_in_child = scalar keys %informative_pos_read_ids_REF;
		  my $num_of_reads_calling_INFOSNP_ALT_in_child = scalar keys %informative_pos_read_ids_ALT;
		  my $pc_of_reads_calling_INFOSNP_REF_in_child  = sprintf("%.2f",(($num_of_reads_calling_INFOSNP_REF_in_child/$num_total)*100));
		  my $pc_of_reads_calling_INFOSNP_ALT_in_child  = sprintf("%.2f",(($num_of_reads_calling_INFOSNP_ALT_in_child/$num_total)*100));
		  my $child_is_het = 0;
		  my $genotype = "$ref"."$ref"; # default
		  if ($num_total > $min_depth)
			{ if (($pc_of_reads_calling_INFOSNP_ALT_in_child >= 30) and ($pc_of_reads_calling_INFOSNP_ALT_in_child < 70))
				{ $genotype = "$ref"."$alt"; $child_is_het++; }
			  elsif ($pc_of_reads_calling_INFOSNP_ALT_in_child >= 70)
				{ $genotype = "$alt"."$alt"; }
			}
		  else
			{ print OUT1 "the min. coverage depth of $min_depth was not met so we revert to a default position of assuming this genotype is homozygous for the reference allele\n"; }
		  $genotype_child = $genotype;
		  print OUT1 "in the child, this position has ACGT counts of $num_A|$num_C|$num_G|$num_T. The ALT base, $alt, has a MAF of $pc_of_reads_calling_INFOSNP_ALT_in_child%\n";
		  print OUT1 "genotype of child : $genotype_child\n";
		  if ($child_is_het == 0)
			{ print OUT1 "unable to proceed as at this position the child is not considered a heterozygote (we require > $min_depth reads [here we see $num_total] and that the % of SNPs calling the SNP is > 30 and < 70; here we see $pc_of_reads_calling_INFOSNP_ALT_in_child%)\n"; }
		  next if ($child_is_het == 0); # skip if the child genotype is not heterozygous, as this will not help us phase things
		  
		  # 7. create a 2x2 count table
		  # it will show for the reads that cover both positions, i.e. the DNM and the 'informative' (= can phase with it) position, the following. Variable names (defined later in this script) are shown here for reference.
		  # 									 	no. of reads calling REF at informative SNP pos		no. of reads calling ALT at informative SNP pos
		  # no. of reads calling REF at DNM pos  	$num_reads_DNM_REF_infosnp_REF						$num_reads_DNM_REF_infosnp_ALT
		  # no. of reads calling ALT at DNM pos		$num_reads_DNM_ALT_infosnp_REF						$num_reads_DNM_ALT_infosnp_ALT
		  # to generate this table we will cross-reference read IDs from the child at the DNM locus and the child at the informative locus
		  
		  # of those reads that contain REF at the DNM pos *AND* cover the informative location, how many call the informative SNP as REF and how many as ALT?
		  my $num_reads_DNM_REF_infosnp_REF = 0; my $num_reads_DNM_REF_infosnp_ALT = 0;
		  while((my $dnm_read_id,my $irrel)=each(%dnm_read_ids_REF_with_coords))
			{ my $start = $dnm_read_ids_REF_with_coords{$dnm_read_id}{start};
			  my $end   = $dnm_read_ids_REF_with_coords{$dnm_read_id}{end};
			  next if (($pos <= $start) or ($pos >= $end)); # CHECKPOINT: we don't consider this DNM-covering read because it doesn't cover the informative base (i.e. the one which will allow us to phase it)
			  while((my $variant_read_id,my $irrel)=each(%informative_pos_read_ids_REF))
				{ if ($dnm_read_id eq $variant_read_id)
					{ $num_reads_DNM_REF_infosnp_REF++;
					}
				}
			  while((my $variant_read_id,my $irrel)=each(%informative_pos_read_ids_ALT))
				{ if ($dnm_read_id eq $variant_read_id)
					{ $num_reads_DNM_REF_infosnp_ALT++;
					}
				}
			}
		
		  # of those reads that contain ALT at the DNM pos *AND* cover the informative location, how many call the informative SNP as REF and how many as ALT?
		  my $num_reads_DNM_ALT_infosnp_REF = 0; my $num_reads_DNM_ALT_infosnp_ALT = 0;
		  while((my $dnm_read_id,my $irrel)=each(%dnm_read_ids_ALT_with_coords))
			{ my $start = $dnm_read_ids_ALT_with_coords{$dnm_read_id}{start};
			  my $end   = $dnm_read_ids_ALT_with_coords{$dnm_read_id}{end};
			  next if (($pos <= $start) or ($pos >= $end)); # CHECKPOINT: we don't consider this DNM-covering read because it doesn't cover the informative base (i.e. the one which will allow us to phase it)
			  while((my $variant_read_id,my $irrel)=each(%informative_pos_read_ids_REF))
				{ if ($dnm_read_id eq $variant_read_id)
					{ $num_reads_DNM_ALT_infosnp_REF++;
					}
				}
			  while((my $variant_read_id,my $irrel)=each(%informative_pos_read_ids_ALT))
				{ if ($dnm_read_id eq $variant_read_id)
					{ $num_reads_DNM_ALT_infosnp_ALT++;
					}
				}
			}
		  
		  my $num_of_reads_calling_DNM_ALT_and_infosnp    = $num_reads_DNM_ALT_infosnp_REF+$num_reads_DNM_ALT_infosnp_ALT;
		  my $pc_of_reads_calling_DNM_ALT_and_infosnp 	  = sprintf("%.2f",(($num_of_reads_calling_DNM_ALT_and_infosnp/$num_of_reads_calling_DNM_ALT_in_child)*100));
		  my $pc_of_reads_calling_DNM_ALT_and_infosnp_REF = sprintf("%.2f",(($num_reads_DNM_ALT_infosnp_REF/$num_of_reads_calling_DNM_ALT_and_infosnp)*100));
		  my $pc_of_reads_calling_DNM_ALT_and_infosnp_ALT = sprintf("%.2f",(($num_reads_DNM_ALT_infosnp_ALT/$num_of_reads_calling_DNM_ALT_and_infosnp)*100));
		  
		  my $num_of_reads_calling_DNM_REF_and_infosnp    = $num_reads_DNM_REF_infosnp_REF+$num_reads_DNM_REF_infosnp_ALT;
		  my $pc_of_reads_calling_DNM_REF_and_infosnp 	  = sprintf("%.2f",(($num_of_reads_calling_DNM_REF_and_infosnp/$num_of_reads_calling_DNM_REF_in_child)*100));
		  my $pc_of_reads_calling_DNM_REF_and_infosnp_REF = sprintf("%.2f",(($num_reads_DNM_REF_infosnp_REF/$num_of_reads_calling_DNM_REF_and_infosnp)*100));
		  my $pc_of_reads_calling_DNM_REF_and_infosnp_ALT = sprintf("%.2f",(($num_reads_DNM_REF_infosnp_ALT/$num_of_reads_calling_DNM_REF_and_infosnp)*100));
		  
		  print OUT1 "of the $num_of_reads_calling_DNM_ALT_in_child reads which call the DNM ALT, $num_of_reads_calling_DNM_ALT_and_infosnp ($pc_of_reads_calling_DNM_ALT_and_infosnp%) also make a call at the informative (= can phase with it) locus\n";
		  print OUT1 "the number of reads that make biallelic calls, i.e. REF or ALT, at both loci ('DNM' and 'phasingSNP') are as follows:\n\n";
		  print OUT1 "\tphasingSNP REF ($ref)\tphasingSNP ALT ($alt)\n";
		  print OUT1 "DNM REF ($snps_per_fam{$family}{ref})\t$num_reads_DNM_REF_infosnp_REF\t$num_reads_DNM_REF_infosnp_ALT\n";
		  print OUT1 "DNM ALT ($snps_per_fam{$family}{alt})\t$num_reads_DNM_ALT_infosnp_REF\t$num_reads_DNM_ALT_infosnp_ALT\n";
		  print OUT1 "\nof the $num_of_reads_calling_DNM_ALT_and_infosnp reads that call the DNM ALT (i.e. $snps_per_fam{$family}{alt}), $pc_of_reads_calling_DNM_ALT_and_infosnp_REF% also call the phasingSNP REF and $pc_of_reads_calling_DNM_ALT_and_infosnp_ALT% also call the phasingSNP ALT\n";
		  print OUT1 "of the $num_of_reads_calling_DNM_REF_and_infosnp reads that call the DNM REF (i.e. $snps_per_fam{$family}{ref}), $pc_of_reads_calling_DNM_REF_and_infosnp_REF% also call the phasingSNP REF and $pc_of_reads_calling_DNM_REF_and_infosnp_ALT% also call the phasingSNP ALT\n";
		  
		  # does DNM REF occur more often with phasingSNP REF or phasingSNP ALT, and from which parent is this base sourced?
		  my @pc_of_reads_mapping_with_DNM_REF        = ([$pc_of_reads_calling_DNM_REF_and_infosnp_REF,$ref,'REF'],[$pc_of_reads_calling_DNM_REF_and_infosnp_ALT,$alt,'ALT']);
		  my @sorted_pc_of_reads_mapping_with_DNM_REF = map { $_->[0] } sort { $b->[1] <=> $a->[1] } map { [$_, $_->[0]] } @pc_of_reads_mapping_with_DNM_REF;
		  my $snp_which_phases_with_DNM_REF 		  = $sorted_pc_of_reads_mapping_with_DNM_REF[0][1];
		  my $type_of_snp_which_phases_with_DNM_REF   = $sorted_pc_of_reads_mapping_with_DNM_REF[0][2];
		  
		  # does DNM ALT occur more often with phasingSNP REF or phasingSNP ALT, and from which parent is this base sourced?
		  my @pc_of_reads_mapping_with_DNM_ALT        = ([$pc_of_reads_calling_DNM_ALT_and_infosnp_REF,$ref,'REF'],[$pc_of_reads_calling_DNM_ALT_and_infosnp_ALT,$alt,'ALT']);
		  my @sorted_pc_of_reads_mapping_with_DNM_ALT = map { $_->[0] } sort { $b->[1] <=> $a->[1] } map { [$_, $_->[0]] } @pc_of_reads_mapping_with_DNM_ALT;
		  my $snp_which_phases_with_DNM_ALT 		  = $sorted_pc_of_reads_mapping_with_DNM_ALT[0][1];
		  my $type_of_snp_which_phases_with_DNM_ALT   = $sorted_pc_of_reads_mapping_with_DNM_ALT[0][2];
		
		  my $mother_REF = ''; my $father_REF = ''; my $mother_ALT = ''; my $father_ALT = '';
		  if ($genotype_mother =~ /^(\w{1})(\w{1})$/) { $mother_REF = $1; $mother_ALT = $2; }
		  if ($genotype_father =~ /^(\w{1})(\w{1})$/) { $father_REF = $1; $father_ALT = $2; } elsif ($genotype_father =~ /^(\w{1})$/) { $father_REF = $1; $father_ALT = $1; }
		  
		  # to resolve the pattern of inheritance we will first look to see which base the DNM ALT phases with
		  my $inheritance = 'unresolved'; my $in_phase_base = '';
		  if ($type_of_snp_which_phases_with_DNM_ALT eq 'REF') # i.e. the DNM ALT phases with phasingSNP REF
			{ $in_phase_base = $snp_which_phases_with_DNM_ALT;
			  if 	(($mother_REF eq $snp_which_phases_with_DNM_ALT) and ($father_REF ne $snp_which_phases_with_DNM_ALT)) { $inheritance = 'maternal'; }
			  elsif (($mother_REF ne $snp_which_phases_with_DNM_ALT) and ($father_REF eq $snp_which_phases_with_DNM_ALT)) { $inheritance = 'paternal'; }
			}
		  elsif ($type_of_snp_which_phases_with_DNM_ALT eq 'ALT') # i.e. the DNM ALT phases with phasingSNP ALT
		    { $in_phase_base = $snp_which_phases_with_DNM_ALT;
			  if 	(($mother_ALT eq $snp_which_phases_with_DNM_ALT) and ($father_ALT ne $snp_which_phases_with_DNM_ALT)) { $inheritance = 'maternal'; }
			  elsif (($mother_ALT ne $snp_which_phases_with_DNM_ALT) and ($father_ALT eq $snp_which_phases_with_DNM_ALT)) { $inheritance = 'paternal'; }
			}
		  
		  # we won't necessarily have resolved the matter using the above if/elsif, though.
		  # consider FAM44 which has genotypes of mother (GA), father (GG) and child (GA) and a 2x2 table as follows:
		  #				phasingSNP REF (G)	phasingSNP ALT (A)
		  # DNM REF (C)	10					321
		  # DNM ALT (G)	297					2
		  # of the 299 reads that call the DNM ALT (i.e. G), 99.33% also call the phasingSNP REF and 0.67% also call the phasingSNP ALT
		  # i.e. $type_of_snp_which_phases_with_DNM_ALT eq 'REF' and $snp_which_phases_with_DNM_ALT eq 'G'
		  # at this point $inheritance eq 'unresolved' because $mother_REF eq 'G' *AND* $father_REF eq 'G'
		  # so, let's take at a look at what the DNM REF phases with. It is VERY IMPORTANT to note that compared to the above if/elsif, the order in which $inheritance options are listed is reversed. Don't forget we only care what DNM REF phases with because from this we can infer what DNM ALT is in phase with - and this is all that we wish to report
		  # we will conclude for FAM44 that because the DNM REF is on the same haplotype as the phasingSNP A, and that A can only come from the mother, then by extension the G that is the DNM ALT can only have come from the father
		  if ($inheritance eq 'unresolved')
			{ if ($type_of_snp_which_phases_with_DNM_REF eq 'REF') # i.e. the DNM REF phases with phasingSNP REF, which means by implication the DNM ALT phases with phasingSNP ALT
				{ $in_phase_base = $snp_which_phases_with_DNM_ALT;
				  if 	(($mother_REF eq $snp_which_phases_with_DNM_REF) and ($father_REF ne $snp_which_phases_with_DNM_REF)) { $inheritance = 'paternal'; }
				  elsif (($mother_REF ne $snp_which_phases_with_DNM_REF) and ($father_REF eq $snp_which_phases_with_DNM_REF)) { $inheritance = 'maternal'; }
				}
			  elsif ($type_of_snp_which_phases_with_DNM_REF eq 'ALT') # i.e. the DNM REF phases with phasingSNP ALT, which means by implication the DNM ALT phases with phasingSNP REF
				{ $in_phase_base = $snp_which_phases_with_DNM_ALT;
				  if 	(($mother_ALT eq $snp_which_phases_with_DNM_REF) and ($father_ALT ne $snp_which_phases_with_DNM_REF)) { $inheritance = 'paternal'; }
				  elsif (($mother_ALT ne $snp_which_phases_with_DNM_REF) and ($father_ALT eq $snp_which_phases_with_DNM_REF)) { $inheritance = 'maternal'; }
				}
			}
		  
		  if ($pc_of_reads_calling_DNM_ALT_and_infosnp_REF == $pc_of_reads_calling_DNM_ALT_and_infosnp_ALT)
			{ print OUT1 "conclusion: reads containing the DNM ALT are as equally likely to derive from the strand containing the phasingSNP REF as the phasingSNP ALT\n";
			  $inheritance = 'unresolved';
			}
		  else
			{ print OUT1 "conclusion 1: reads containing the DNM REF are more likely to derive from the strand containing the phasingSNP $type_of_snp_which_phases_with_DNM_REF\n";
			  print OUT1 "conclusion 2: reads containing the DNM ALT are more likely to derive from the strand containing the phasingSNP $type_of_snp_which_phases_with_DNM_ALT\n";
			}
		  print OUT1 "predicted inheritance: $inheritance\n";
		  my $total_depth_of_cov_phasing_snp = $num_reads_DNM_REF_infosnp_REF+$num_reads_DNM_REF_infosnp_ALT+$num_reads_DNM_ALT_infosnp_REF+$num_reads_DNM_ALT_infosnp_ALT;

		  $printed_something++;
		  my $read_counts = '.';
		  if    ($dnm_indel eq '.') 	 	{ $read_counts = "$dnm_A|$dnm_C|$dnm_G|$dnm_T"; 	}
		  elsif ($dnm_indel =~ /^\-\d+/) 	{ $read_counts = $num_reads_containing_the_dnm_del; }
		  elsif ($dnm_indel =~ /^\+\d+\w+/) { $read_counts = $num_reads_containing_the_dnm_ins; }
		  print OUT2 "$family\t$region_sequenced\t$length_of_region_sequenced\t$snps_per_fam{$family}{chr}:$snps_per_fam{$family}{pos}\t$snps_per_fam{$family}{ref}\t$snps_per_fam{$family}{alt}\t$read_counts\t$dnm_total\t$snps_per_fam{$family}{id}\t1\t$dbsnp_id\t1\t$dbsnp_id\t$total_depth_of_cov_phasing_snp\t$in_phase_base\t$genotype_mother\t$genotype_father\t$genotype_child\t$inheritance\t$dist\t$dist\n";
		  last unless ($inheritance eq 'unresolved');
		}
	 if ($printed_something == 0)
		{ my $read_counts = '.';
		  if    ($dnm_indel eq '.') 	 	{ $read_counts = "$dnm_A|$dnm_C|$dnm_G|$dnm_T"; 	}
		  elsif ($dnm_indel =~ /^\-\d+/) 	{ $read_counts = $num_reads_containing_the_dnm_del; }
		  elsif ($dnm_indel =~ /^\+\d+\w+/) { $read_counts = $num_reads_containing_the_dnm_ins; }
		  print OUT2 "$family\t$region_sequenced\t$length_of_region_sequenced\t$snps_per_fam{$family}{chr}:$snps_per_fam{$family}{pos}\t$snps_per_fam{$family}{ref}\t$snps_per_fam{$family}{alt}\t$read_counts\t$dnm_total\t$snps_per_fam{$family}{id}\t0\t.\t0\t.\t.\t.\t.\t.\t.\tunresolved\t.\t.\n";
		}
	}

close(OUT1) or die $!; close(OUT2) or die $!;
exit 1;
