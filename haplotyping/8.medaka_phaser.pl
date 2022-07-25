use strict;
use warnings;

# REQUIREMENTS
my $basedir  = '/t1-data/project/pregcare/sbush';
my $samtools = '/t1-data/project/GorielyLab2021/sbush/programs/samtools-1.12/bin/samtools';
my $coords   = "$basedir/ont_coords.tsv";
my $fatal    = 0;
if (!(-e($coords)))   { $fatal++; print "ERROR: cannot find $coords\n";   }
if (!(-e($samtools))) { $fatal++; print "ERROR: cannot find $samtools\n"; }
exit 1 if ($fatal > 0);

# PARAMETERS
my $known_only   = 'no'; # 'no' or 'yes' | This parameter is used to restrict the phase set only to SNPs considered already known, i.e. in dbSNP
my $pass_only    = 'no'; # 'no' or 'yes' | This parameter is used to restrict the phase set only to SNPs flagged by Medaka as 'PASS', i.e. meeting Medaka's internal QC filters
my $min_depth    = 10;   # numeric       | This parameter is used to restrict the phase set only to SNPs with a total depth of coverage at that locus greater than this threshold

# OUTPUT
my $out_file = "$basedir/inheritance.medaka.known_only_$known_only.pass_only_$pass_only.min_dp_$min_depth.tsv";
open(OUT,'>',$out_file) or die $!;
print OUT "Family\tRegion sequenced\tLength of region sequenced (bp)\tDe novo mutation (chr:position)\tReference allele (for DNM)\tVariant allele (for DNM)\tRead counts (for DNM, given as A|C|G|T only if called by Medaka in the proband)\tTotal depth of coverage (for DNM)\t";
print OUT "dbSNP ID (for DNM)\tHas Medaka detected the DNM in the proband? (expectation: yes)\tHas Medaka detected the DNM in the mother? (expectation: no)\tHas Medaka detected the DNM in the father? (expectation: no)\tHas Medaka phased the DNM in the proband?\tIf yes, how many SNPs are in the same phase set, not including the DNM?\tSNPs in the proband's phase set, not including the DNM (listed in order of distance from DNM and reported as dbSNP IDs if available, else 'chr:pos')\tHow many INFORMATIVE SNPs (= different between parents) are in the same phase set, not including the DNM?\tInformative SNPs (= different between parents) in the proband's phase set, not including the DNM (listed in order of distance from DNM and reported as dbSNP IDs if available, else 'chr:pos')\tTotal depth of coverage for informative SNPs in the proband's phase set, listed in the order previously stated\tBases that phase with the DNM (in the proband, listed in the order previously stated)\tGenotypes, listed in order previously stated (mother)\tGenotypes, listed in order previously stated (father)\tGenotypes, listed in order previously stated (proband)\tPattern of inheritance\tMin. distance of an informative SNP from DNM (bp)\tMax. distance of an informative SNP from DNM (bp)\tReason the pattern of inheritance was unresolved (if applicable)\n";

# STORE DNM SNPs FOR EACH FAMILY AND THE BARCODES OF EACH MOTHER/FATHER/CHILD
my %snps_per_fam = (); my %bcs_per_fam = (); my %regions_sequenced = ();
open(IN,$coords) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line); $line =~ s/\r$//g; $line =~ s/\s$//g;
	  my @line = split(/\t/,$line);
	  my $onts = $line[0]; my $bc_line = $line[1]; my $family = $line[2]; my $region_sequenced = $line[3]; my $dbsnp_id = $line[4]; my $dnm_loc = $line[5]; my $ref = $line[6]; my $alt = $line[7]; my $dnm_type = $line[8];
	  next if ($dnm_type ne 'SNP');
	  my $dnm_chr = ''; my $dnm_pos = '';
	  if ($dnm_loc =~ /^(.*?)\:(\d+)$/) { $dnm_chr = $1; $dnm_pos = $2; } else { print "ERROR: unable to parse $dnm_loc\n"; exit 1; }
	  $regions_sequenced{$family} = $region_sequenced;
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
	{ next if (!(exists($bcs_per_fam{$family})));
	
	  # determine which MinION run is associated with which member of the trio (and their barcode)
	  my $onts = $snps_per_fam{$family}{ont};
	  next if ($onts !~ /^(.*?)\|(.*?)\|(.*?)$/);
	  my $mother_ont = ''; my $father_ont = ''; my $child_ont = '';
	  if ($onts =~ /^(.*?)\|(.*?)\|(.*?)$/)
		{ $mother_ont = $1; $father_ont = $2; $child_ont = $3; }
	  
	  # PARSE VCFs FOR MOTHER, FATHER, AND CHILD
	  my %coord_lookup = ();
	  my $bc_line = $bcs_per_fam{$family};
	  my @bc_line = split(/\|/,$bc_line);
	  my $id_of_dnm = '.';
	  my $is_dnm_in_child = 'no'; my $is_dnm_in_mother = 'no'; my $is_dnm_in_father = 'no';
	  my $position_number = 0;
	  my $dnm_in_child_is_phased = 'no';
	  my $dnm_A = 0; my $dnm_T = 0; my $dnm_C = 0; my $dnm_G = 0; my $dnm_total = 0;
	  my $child_bam = '';
	  my @phase_set_snps = ();
	  for(my $x=0;$x<@bc_line;$x++)
		{ my $bc = $bc_line[$x]; my $person = ''; my $ont = '';
		  if 	($x == 0) { $person = 'mother'; $ont = $mother_ont; }
		  elsif ($x == 1) { $person = 'father'; $ont = $father_ont; }
		  elsif ($x == 2) { $person = 'child';  $ont = $child_ont;  }
		  my $dnm = "$snps_per_fam{$family}{chr}.$snps_per_fam{$family}{pos}.$snps_per_fam{$family}{ref}.$snps_per_fam{$family}{alt}";
		  my $bam = "$basedir/$ont/6.bam/$bc.bam";
		  my $vcf = "$basedir/$ont/7.vcf/$bc.vcf";
		  next if ( (!(-e($bam))) or (!(-e($vcf))) );
		  if ($person eq 'child')
			{ $child_bam = $bam; }
		  
		  # (1) do we detect the DNM in any member of the trio?
		  my %dnm_phase_set = (); my $gt_of_dnm = '';			  
		  open(VCF,$vcf) or die $!;
		  while(<VCF>)
			{ next if ($_ =~ /^\#/);
			  my $line = $_; chomp($line);
			  my @line = split(/\t/,$line);
			  my $chr = $line[0]; my $pos = $line[1]; my $dbsnp_id = $line[2]; my $ref = $line[3]; my $alt = $line[4]; my $filter = $line[6]; my $info_line = $line[7]; my $phase_line = $line[9];
			  my $len_ref = length($ref); my $len_alt = length($alt);
			  next if (($len_ref != 1) or ($len_alt != 1));
			  if ($dnm eq "$chr.$pos.$ref.$alt") # we are NOT having the line "next if ($filter ne 'PASS')". This is because, for the purpose of the DNM *ONLY*, we will permit lower-quality/non-passing calls because we expect to find this - we know it is there.
				{ $id_of_dnm = $dbsnp_id;
				  if 	($person eq 'mother') { $is_dnm_in_mother = 'yes'; }
				  elsif ($person eq 'father') { $is_dnm_in_father = 'yes'; }
				  elsif ($person eq 'child')
					{ $is_dnm_in_child  = 'yes';
					  if ($phase_line =~ /^(.*?)\:.+?\:(\d+)$/) # a $phase_line looks like this: "0|1:22:21813596"
						{ my $gt = $1; my $phase_set = $2;
						  if (($gt ne '0|1') and ($gt ne '1|0'))
							{ print "ERROR ($family): the DNM ($chr:$pos) is not considered by Medaka to be a heterozygote - unable to proceed\n"; exit 1; }
						  $dnm_in_child_is_phased = 'yes';
						  $dnm_phase_set{$phase_set} = $gt;
						  $gt_of_dnm = $gt;
						}
					  # calculate, from the pileup, the nucleotide counts for the DNM
					  my $pileup = '';
					  open(IN,"$samtools mpileup -r $chr:$pos-$pos -aa $bam |") or die $!;
					  while(<IN>)
						{ my $line = $_; chomp($line);
						  my @line = split(/\t/,$line);
						  $pileup = $line[4];
						}
					  close(IN) or die $!;
					  my %nuc_counts = (A => 0, G => 0, C => 0, T => 0);
					  my $total = 0; my $skip = 0;
					  for my $n (split(//, $pileup))
						{ my $skip_this = 0;
						  if ($n =~ /^(\d+)$/)
							{ $skip = $1; $skip_this++; }
						  elsif ($skip > 0)
							{ $skip--; $skip_this++; }
						  next if ($skip_this > 0);
						  if (exists $nuc_counts{uc $n})
							{ $nuc_counts{uc $n}++;
							  $total++;
							}
						}
					  $dnm_A     = $nuc_counts{'A'};
					  $dnm_T     = $nuc_counts{'T'};
					  $dnm_C     = $nuc_counts{'C'};
					  $dnm_G     = $nuc_counts{'G'};
					  $dnm_total = $total;
					}
				}
			}
		  close(VCF) or die $!;
		  
		  # we need to identify which variants phase with the DNM. This will be every variant in the phase set on the same side of the pipe.
		  if 	($gt_of_dnm eq '1|0') { $position_number = 1; } # i.e. the "1" is in position 1
		  elsif ($gt_of_dnm eq '0|1') { $position_number = 2; } # i.e. the "1" is in position 2
		  
		  # (2) for the child sample, identify which SNPs are in the phase set
		  if ($person eq 'child')
			{ open(VCF,$vcf) or die $!;
			  while(<VCF>)
				{ next if ($_ =~ /^\#/);
				  my $line = $_; chomp($line);
				  my @line = split(/\t/,$line);
				  my $chr = $line[0]; my $pos = $line[1]; my $dbsnp_id = $line[2]; my $ref = $line[3]; my $alt = $line[4]; my $filter = $line[6]; my $info_line = $line[7]; my $phase_line = $line[9];
				  my $len_ref = length($ref); my $len_alt = length($alt);
				  next if (($len_ref != 1) or ($len_alt != 1));
				  next if (($dnm ne "$chr.$pos.$ref.$alt") and ($pass_only eq 'yes') and ($filter ne 'PASS')); # CHECKPOINT: restrict analysis only to passing SNPs. An exception is made for the DNM.
				  next if (($dnm ne "$chr.$pos.$ref.$alt") and ($known_only eq 'yes') and ($dbsnp_id !~ /^rs/)); # CHECKPOINT: restrict analysis only to known (=dbSNP) SNPs. An exception is made for the DNM.
				  my $depth = 0;
				  if ($info_line =~ /^.*?DP\=(\d+)\;.*?$/) { $depth = $1; }
				  next if (($dnm ne "$chr.$pos.$ref.$alt") and ($depth < $min_depth)); # CHECKPOINT: restrict analysis only to SNPs with total depth of coverage at that locus > 10. An exception is made for the DNM.
				  if ($dbsnp_id eq '.') { $dbsnp_id = "$chr:$pos"; }
				  if ($phase_line =~ /^(.*?)\:.+?\:(\d+)$/) # a $phase_line looks like this: "0|1:22:21813596"
					{ my $gt = $1; my $phase_set = $2;
					  if (exists($dnm_phase_set{$phase_set}))
						{ my $dist = 'NA';
						  if ($pos >= $snps_per_fam{$family}{pos})
							{ $dist = $pos-$snps_per_fam{$family}{pos}; }
						  else
							{ $dist = $snps_per_fam{$family}{pos}-$pos; }
						  if ($dist =~ /NA/) { print "ERROR: unable to calculate distance between $pos and $snps_per_fam{$family}{pos}\n"; exit 1; }
						  if ($dbsnp_id eq '.')	{ $dbsnp_id = "$chr:$pos"; }
						  my $genotype = '';
						  if (($gt eq '0|1') or ($gt eq '1|0'))
							{ $genotype = "$ref"."$alt"; }
						  elsif ($gt eq '0|0')
							{ $genotype = "$ref"."$ref"; }
						  elsif ($gt eq '1|1')
							{ $genotype = "$alt"."$alt"; }
						  push(@phase_set_snps,[$dist,$dbsnp_id,$gt,$ref,$alt,$genotype,"$chr:$pos",$depth]);
						}
					}
				}
			  close(VCF) or die $!;
			}
		}
	  next if ($child_bam eq '');
	  next if (!(-e($child_bam)));
	  my @sorted_phase_set_snps = map { $_->[0] } sort { $a->[1] <=> $b->[1] } map { [$_, $_->[0]] } @phase_set_snps;
	  my $num_of_phase_set_snps = @sorted_phase_set_snps; $num_of_phase_set_snps-- unless ($num_of_phase_set_snps == 0); # we make a single subtraction so that when reporting the results, we do not count the DNM in the total
	  my $phase_set_snps = ''; my $genotypes_child = ''; my $bases_that_phase_with_the_dnm = ''; my $depths_of_phasing_snps = '';
	  my $min_distance_to_dnm = 'NA'; my $max_distance_to_dnm = 'NA';
	  if (($num_of_phase_set_snps > 0) and ($position_number > 0))
		{ # we want to know which of the phase set SNPs phase together with this DNM. To do this, we need to look at where we find the 1. Is it for the $ref or $alt position - positions 1 and 2, respectively?
		  for(my $x=1;$x<@sorted_phase_set_snps;$x++)
			{ my $dist		= $sorted_phase_set_snps[$x][0];
			  if ($x == 1) { $min_distance_to_dnm = $dist; }
			  $max_distance_to_dnm = $dist; # will auto-update and always be the last value
			  my $dbsnp_id  = $sorted_phase_set_snps[$x][1];
			  my $gt_of_snp = $sorted_phase_set_snps[$x][2];
			  my $ref       = $sorted_phase_set_snps[$x][3];
			  my $alt       = $sorted_phase_set_snps[$x][4];
			  my $genotype  = $sorted_phase_set_snps[$x][5];
			  my $coord		= $sorted_phase_set_snps[$x][6];
			  my $depth     = $sorted_phase_set_snps[$x][7];
			  my $phasing_base = '';
			  if ($gt_of_snp eq '1|0')
				{ if    ($position_number == 1) { $phasing_base = $alt; } # if the DNM's "1" is in position 1, then the DNM has phased with whichever base is on the left hand side of the pipe, in this case base 1 out of 1|0, i.e. ALT
				  elsif ($position_number == 2) { $phasing_base = $ref; } # if the DNM's "1" is in position 2, then the DNM has phased with whichever base is on the right hand side of the pipe, in this case base 0 out of 1|0, i.e. REF
				}
			  elsif ($gt_of_snp eq '0|1')
				{ if    ($position_number == 1) { $phasing_base = $ref; } # if the DNM's "1" is in position 1, then the DNM has phased with whichever base is on the left hand side of the pipe, in this case base 0 out of 0|1, i.e. REF
				  elsif ($position_number == 2) { $phasing_base = $alt; } # if the DNM's "1" is in position 2, then the DNM has phased with whichever base is on the right hand side of the pipe, in this case base 1 out of 0|1, i.e. ALT
				}
			  $phase_set_snps  				 .= "$dbsnp_id, ";
			  $genotypes_child 				 .= "$genotype, ";
			  $depths_of_phasing_snps 		 .= "$depth, ";
			  $bases_that_phase_with_the_dnm .= "$phasing_base, ";
			  $coord_lookup{$coord}++;
			}
		}
	  $phase_set_snps =~ s/\, $//; $depths_of_phasing_snps =~ s/\, $//; $genotypes_child =~ s/\, $//; $bases_that_phase_with_the_dnm =~ s/\, $//;
	  my $original_num_of_phase_set_snps = $num_of_phase_set_snps;
	  my $original_phase_set_snps = $phase_set_snps;
	  
	  # (3) now that we know which bases phase with the DNM in the child, we will look up what calls have been made at these positions in the mother and father
	  my %mother_snps = (); my %father_snps = ();
	  for(my $x=0;$x<=1;$x++)
		{ my $bc = $bc_line[$x]; my $person = ''; my $ont = '';
		  if 	($x == 0) { $person = 'mother'; $ont = $mother_ont; }
		  elsif ($x == 1) { $person = 'father'; $ont = $father_ont; }
		  my $vcf = "$basedir/$ont/7.vcf/$bc.vcf";
		  next if (!(-e($vcf)));
		  open(VCF,$vcf) or die $!;
		  while(<VCF>)
			{ next if ($_ =~ /^\#/);
			  my $line = $_; chomp($line);
			  my @line = split(/\t/,$line);
			  my $chr = $line[0]; my $pos = $line[1]; my $dbsnp_id = $line[2]; my $ref = $line[3]; my $alt = $line[4]; my $filter = $line[6]; my $gt_line = $line[9];
			  my $len_ref = length($ref); my $len_alt = length($alt);
			  next if (($len_ref != 1) or ($len_alt != 1)); # CHECKPOINT: exclude non-SNPs
			  next if (($pass_only eq 'yes') and ($filter ne 'PASS')); # CHECKPOINT: exclude non-passing SNPs. We can relax this requirement in this case. Why? Because we are looking for phase-set SNPs we have already identified in the child. So, there is a prior they must be in at least one parent too as de novo mutations (an alternative explanation for why a SNP is found in the child but not a parent too) are very rare.
			  next if (($known_only eq 'yes') and ($dbsnp_id !~ /^rs/)); # CHECKPOINT: restrict analysis only to known (=dbSNP) SNPs
			  my $gt = '';
			  if ($gt_line =~ /^(.*?)\:.+$/) # a $gt_line looks like this: "0/1:22"
				{ $gt = $1; }
			  my $genotype = "$ref"."$ref"; # default is WT...
			  if (($person eq 'father') and ($chr =~ /^X$/)) # which, if we're looking at the X chromosome of the father, will only contain a single base
				{ $genotype = $ref; }
			  if (($gt eq '0/1') or ($gt eq '0|1') or ($gt eq '1|0'))
				{ $genotype = "$ref"."$alt"; }
			  if (($gt eq '1/1') or ($gt eq '1|1'))
				{ $genotype = "$alt"."$alt";
				  if (($person eq 'father') and ($chr =~ /^X$/))
					{ $genotype = $alt; }
				}
			  if (exists($coord_lookup{"$chr:$pos"}))
				{ if ($person eq 'mother')
					{ $mother_snps{"$chr:$pos"} = $genotype; }
				  elsif ($person eq 'father')
					{ $father_snps{"$chr:$pos"} = $genotype; }
				}
			}
		  close(VCF) or die $!;
		}
	  my %parent_genotypes = ();
	  my $genotypes_mother = ''; my $genotypes_father = '';
	  for(my $x=1;$x<@sorted_phase_set_snps;$x++)
		{ my $dbsnp_id = $sorted_phase_set_snps[$x][1];
		  my $coord    = $sorted_phase_set_snps[$x][6];
		  my $chr      = '';
		  if ($coord   =~ /^(.*?)\:\d+$/) { $chr = $1; }
		  my $ref_call = $sorted_phase_set_snps[$x][3];
		  my $wt = "$ref_call"."$ref_call";
		  my $gt_m = ''; my $gt_f = '';
		  if (exists($mother_snps{$coord}))
			{ $genotypes_mother .= "$mother_snps{$coord}, "; $gt_m = $mother_snps{$coord};   }
		  else
			{ $genotypes_mother .= "$wt, "; 				 $gt_m = $wt;				     }
		  if (exists($father_snps{$coord}))
			{ $genotypes_father .= "$father_snps{$coord}, "; $gt_f = $father_snps{$coord};   }
		  else
			{ if    ($chr !~ /^X$/) { $genotypes_father .= "$wt, "; 	  $gt_f = $wt; 		 }
			  elsif ($chr =~ /^X$/) { $genotypes_father .= "$ref_call, "; $gt_f = $ref_call; }
			}
		  $parent_genotypes{$dbsnp_id}{mother} = $gt_m;
		  $parent_genotypes{$dbsnp_id}{father} = $gt_f;
		}
	  $genotypes_mother =~ s/\, $//; $genotypes_father =~ s/\, $//;
	  
	  # (4) we now have (a) a list of bases which phase with the DNM, (b) the mother's genotype, and (c) the father's genotype
	  # we can iterate through each of A to work out which base can ONLY derive from B or C. In this way, we can predict the pattern of inheritance.
	  # note that a haplotype will not (in human) last the length of the entire chromosome! At some point, it will break down due to recombination
	  # at that point, we would expect to see a shift in inheritance pattern with respect to what has phased with the DNM, e.g. SNPs which phase like so: paternal-paternal-paternal-maternal-maternal
	  # if such a shift occurs, what we do is give greater weight to those SNPs which are closest to the DNM, discarding the others, e.g. in this case, we retain only paternal-paternal-paternal
	  my @bases_that_phase_with_the_dnm = split(", ",$bases_that_phase_with_the_dnm);
	  my @phasing_base_ids = split(", ",$phase_set_snps);
	  my @dp_of_phase_snps = split(", ",$depths_of_phasing_snps);
	  my @genotypes_mother = split(", ",$genotypes_mother);
	  my @genotypes_father = split(", ",$genotypes_father);
	  my @genotypes_child  = split(", ",$genotypes_child);
	  my $new_num_of_phase_set_snps = 0; my $new_bases_that_phase_with_the_dnm = ''; my $new_phase_set_snps = ''; my $new_depths_of_phasing_snps = ''; my $new_genotypes_mother = ''; my $new_genotypes_father = ''; my $new_genotypes_child = '';
	  my %inheritances = (); my $last_inheritance_seen = '';
	  for(my $x=0;$x<@bases_that_phase_with_the_dnm;$x++)
		{ my $phase_snp_id = $phasing_base_ids[$x];
		  my $phasing_base = $bases_that_phase_with_the_dnm[$x];
		  my $gt_mother	   = $genotypes_mother[$x];
		  my $gt_father    = $genotypes_father[$x];
		  my $gt_child	   = $genotypes_child[$x];
		  my $depth		   = $dp_of_phase_snps[$x];
		  my @mother_bases = split(//,$gt_mother); my %mother_bases = map {$_ => 1} @mother_bases;
		  my @father_bases = split(//,$gt_father); my %father_bases = map {$_ => 1} @father_bases;
		  my @child_bases  = split(//,$gt_child);
		  my %potential_origin_of_bases = ();
		  foreach my $child_nt (@child_bases) # because the child is a heterozygote, there will always be two different bases. One must come from the mother, the other the feather.
			{ if (exists($mother_bases{$child_nt})) { $potential_origin_of_bases{$child_nt}{'mother'}++; } # print "$child_nt is potentially ma\n";
			  if (exists($father_bases{$child_nt})) { $potential_origin_of_bases{$child_nt}{'father'}++; } # print "$child_nt is potentially pa\n";
			}
		  my $maternal_base = ''; my $paternal_base = '';
		  while((my $child_nt,my $irrel)=each(%potential_origin_of_bases))
			{ my $no_of_potential_origins = scalar keys %{$potential_origin_of_bases{$child_nt}};
			  if ($no_of_potential_origins == 1) # one of the two alleles can ONLY come from one parent...
				{ while((my $origin,my $irrel)=each(%{$potential_origin_of_bases{$child_nt}}))
					{ if    ($origin eq 'mother') { $maternal_base = $child_nt; }
					  elsif ($origin eq 'father') { $paternal_base = $child_nt; }
					  # print "$phase_snp_id = $child_nt is from $origin\n";
					}
				}
			}
		  my %nts_remaining = ();
		  foreach my $child_nt (@child_bases) # ... which means the other allele can only come from the other one.
			{ next if (($maternal_base eq $child_nt) or ($paternal_base eq $child_nt));
			  $nts_remaining{$child_nt}++;
			}
		  my $no_of_nts_remaining = scalar keys %nts_remaining;
		  if ($no_of_nts_remaining == 1)
			{ while((my $nt,my $irrel)=each(%nts_remaining))
				{ if    ($maternal_base eq '') { $maternal_base = $nt; }
				  elsif ($paternal_base eq '') { $paternal_base = $nt; }
				}
			}
		  my $inheritance = 'unresolved';
		  if ($maternal_base ne $paternal_base)
			{ if    ($maternal_base eq $phasing_base) { $inheritance = 'maternal'; }
			  elsif ($paternal_base eq $phasing_base) { $inheritance = 'paternal'; }
			}
		  #print "$phase_snp_id --> child genotype: $gt_child. base $phasing_base phases with the DNM --> ma: $maternal_base ($gt_mother), pa: $paternal_base ($gt_father) --> $inheritance\n";
		  next if ($inheritance eq 'unresolved');
		  last if (($last_inheritance_seen ne '') and ($inheritance ne $last_inheritance_seen)); # we break out of this loop if we cease to see SNPs whose pattern of inheritance is consistent with the previous
		  $inheritances{$inheritance}++;
		  $last_inheritance_seen = $inheritance;
		  $new_bases_that_phase_with_the_dnm .= "$phasing_base, ";
		  $new_depths_of_phasing_snps .= "$depth, ";
		  $new_genotypes_father .= "$gt_father, ";
		  $new_genotypes_mother .= "$gt_mother, ";
		  $new_genotypes_child  .= "$gt_child, ";
		  $new_phase_set_snps   .= "$phase_snp_id, ";
		  $new_num_of_phase_set_snps++;
		}
	  $num_of_phase_set_snps = $new_num_of_phase_set_snps;
	  $new_bases_that_phase_with_the_dnm =~ s/\, $//; $new_depths_of_phasing_snps =~ s/\, $//; $new_genotypes_mother =~ s/\, $//; $new_genotypes_father =~ s/\, $//; $new_genotypes_child =~ s/\, $//; $new_phase_set_snps =~ s/\, $//;
	  $bases_that_phase_with_the_dnm = $new_bases_that_phase_with_the_dnm;
	  $depths_of_phasing_snps = $new_depths_of_phasing_snps;
	  $genotypes_father = $new_genotypes_father;
	  $genotypes_mother = $new_genotypes_mother;
	  $genotypes_child  = $new_genotypes_child;
	  $phase_set_snps   = $new_phase_set_snps;
	  my $inheritance   = 'unresolved';
	  my $no_possible_inheritances = scalar keys %inheritances;
	  if ($no_possible_inheritances == 1)
		{ while((my $source,my $irrel)=each(%inheritances))
			{ $inheritance = $source; }
		}
	  
	  print "SNPs in the phase set: $phase_set_snps\n";
	  print "phasing bases: $bases_that_phase_with_the_dnm\n";
	  print "ma: $genotypes_mother\npa: $genotypes_father\nch: $genotypes_child\n";
	  print "inheritance for $family: $inheritance\n";
	  if ($phase_set_snps 	eq '') { $phase_set_snps   = '.'; }
	  if ($genotypes_mother eq '') { $genotypes_mother = '.'; }
	  if ($genotypes_father	eq '') { $genotypes_father = '.'; }
	  if ($genotypes_child	eq '') { $genotypes_child  = '.'; }
	  if ($bases_that_phase_with_the_dnm eq '') { $bases_that_phase_with_the_dnm = '.'; }
	  if ($depths_of_phasing_snps eq '') 		{ $depths_of_phasing_snps 		 = '.'; }
	  if ($is_dnm_in_child eq 'no')
		{ $dnm_A = '.'; $dnm_C = '.'; $dnm_G = '.'; $dnm_T = '.'; $dnm_total = '.'; }
	  my $reason_for_lack_of_resolution = '.';
	  if ($inheritance eq 'unresolved')
		{ if ($is_dnm_in_child eq 'no')
			{ $reason_for_lack_of_resolution = 'DNM not called in the child'; }
		  elsif (($is_dnm_in_child eq 'yes') and ($dnm_in_child_is_phased eq 'no'))
			{ $reason_for_lack_of_resolution = 'DNM called in the child but no in-phase SNPs identified'; }
		  elsif (($is_dnm_in_child eq 'yes') and ($dnm_in_child_is_phased eq 'yes') and ($original_num_of_phase_set_snps == 0))
			{ $reason_for_lack_of_resolution = "DNM called in the child and in-phase SNPs identified but the workflow discards them all because none meet user-specified QC requirements"; }
		  elsif (($is_dnm_in_child eq 'yes') and ($dnm_in_child_is_phased eq 'yes') and ($original_num_of_phase_set_snps > 0) and ($num_of_phase_set_snps == 0)) # we get to this point if there are heterozygous in-phase SNPs in the child but these SNPs are (a) homozygotes in both parents (conclusion: v. likely to be an erroneous call in at least one parent as otherwise this indicates a [very rare] DNM in the child), or (b) heterozygotes in both parents (conclusion: you would need to manually resolve this)
			{ # what is the reason a phase set SNP was NOT informative, i.e. in $original_phase_set_snps but not $phase_set_snps?
			  my @original_phase_set_snps    = split(", ",$original_phase_set_snps);
			  my @informative_phase_set_snps = split(", ",$phase_set_snps);
			  my %informative_phase_set_snps = map {$_ => 1} @informative_phase_set_snps;
			  my $reason_line = '';
			  foreach my $dbsnp_id (@original_phase_set_snps)
				{ if (!(exists($informative_phase_set_snps{$dbsnp_id}))) # i.e. this SNP, a heterozygote in the child and in-phase with the DNM, is not informative for phasing. Why is this?
					{ my $gt_mother = $parent_genotypes{$dbsnp_id}{'mother'};
					  my $gt_father = $parent_genotypes{$dbsnp_id}{'father'};
					  if ( ((length($gt_mother))==(length($gt_father))) and ($gt_mother ne $gt_father) ) { print "FATAL ERROR: the SNP $dbsnp_id is thought non-informative with respect to phasing the DNM of $family but it is $gt_mother in mother and $gt_father in father\n"; exit 1; }
					  my $mother_status = ''; my $father_status = '';
					  if ($gt_mother =~ /^(\w{1})(\w{1})$/)
						{ if ($1 eq $2) { $mother_status = 'hom'; } else { $mother_status = 'het'; } }
					  if ($gt_father =~ /^(\w{1})(\w{1})$/)
						{ if ($1 eq $2) { $father_status = 'hom'; } else { $father_status = 'het'; } }
					  elsif ($gt_father =~ /^(\w{1})$/) # this will be the case if we are looking at the X chromosome
						{ $father_status = 'hom'; }
					  if ($mother_status ne $father_status) { print "FATAL ERROR: the SNP $dbsnp_id is thought non-informative with respect to phasing the DNM of $family but it is $mother_status in mother and $father_status in father\n"; exit 1; }
					  if    ($mother_status eq 'hom')
						{ $reason_line .= "$dbsnp_id is homozygous in both parents ($gt_mother), "; }
					  elsif ($mother_status eq 'het')
						{ $reason_line .= "$dbsnp_id is heterozygous in both parents ($gt_mother), "; }
					}
				}
			  $reason_line =~ s/\, $//;
			  $reason_for_lack_of_resolution = "DNM called in the child and in-phase SNPs identified, but none are informative - $reason_line";
			}
		}
	  my $region_sequenced = $regions_sequenced{$family};
	  my $length_of_region_sequenced = 0;
	  if ($region_sequenced =~ /^.*?\:(\d+)\-(\d+)$/)
		{ $length_of_region_sequenced = ($2-$1)+1; }
	  print OUT "$family\t$region_sequenced\t$length_of_region_sequenced\t$snps_per_fam{$family}{chr}:$snps_per_fam{$family}{pos}\t$snps_per_fam{$family}{ref}\t$snps_per_fam{$family}{alt}\t$dnm_A|$dnm_C|$dnm_G|$dnm_T\t$dnm_total\t$id_of_dnm\t$is_dnm_in_child\t$is_dnm_in_mother\t$is_dnm_in_father\t$dnm_in_child_is_phased\t$original_num_of_phase_set_snps\t$original_phase_set_snps\t$num_of_phase_set_snps\t$phase_set_snps\t$depths_of_phasing_snps\t$bases_that_phase_with_the_dnm\t$genotypes_mother\t$genotypes_father\t$genotypes_child\t$inheritance\t$min_distance_to_dnm\t$max_distance_to_dnm\t$reason_for_lack_of_resolution\n";
	}
close(OUT) or die $!;
exit 1;
