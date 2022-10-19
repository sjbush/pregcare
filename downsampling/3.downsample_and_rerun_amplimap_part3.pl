use strict;
use warnings;

# REQUIREMENTS
my $home 		 = '/t1-data/project/GorielyLab2021/sbush/pregcare';
my $lookup_table = "$home/miseq_fam_to_code_lookup.txt";
my $fatal		 = 0;
if (!(-e($lookup_table))) { $fatal++; print "ERROR: cannot find $lookup_table\n"; }
exit 1 if ($fatal > 0);

# PARAMETERS
my @families  = (qw/FAM27 FAM34/);
my @fold 	  = (qw/25 100 500 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 15000 20000/);
my @seeds 	  = (qw/123 456 789 234 567 891 321 654 987 432/);

# OUTPUT
my $out_dir = "$home/pileups";
if (!(-d($out_dir))) { mkdir $out_dir or die $!; }

# STORE LOOKUP TABLE
my %lookup = (); my %print_order = ();
open(IN,$lookup_table) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $family = $line[0]; my $bc = $line[1]; my $code = $line[2];
	  $lookup{$family}{$bc} = $code;
	  push(@{$print_order{$family}},$bc);
	}
close(IN) or die $!;

# FOR EACH FAMILY, PRINT A VERSION OF AMPLIMAP'S PILEUPS_LONG_DETAILED.CSV FILE, RESTRICTED TO THE DNM
foreach my $family (@families)
	{ my $dnm_chr; my $dnm_pos;
	  if    ($family eq 'FAM27') { $dnm_chr = 'chrX';  $dnm_pos = 154031022; }
	  elsif ($family eq 'FAM34') { $dnm_chr = 'chr20'; $dnm_pos = 63407094;  }
	  my $in_dir  = "$home/$family";
	  next if (!(-d($in_dir)));
	  
	  # print output for raw (not downsampled) BAM
	  my $in_file = "$in_dir/raw/analysis/pileups/pileups_long_detailed.csv";
	  next if (!(-e($in_file)));
	  my $out_file = "$out_dir/$family-no_downsampling.csv";
	  open(OUT,'>',$out_file) or die $!;
	  print OUT "label,sample,chr,pos,target_id,target_type,ref,alts,probes,count_hq_A,count_hq_C,count_hq_G,count_hq_T,count_hq_INS,count_hq_DEL,count_hq_SKIP,number_called_hq,ref_hq_count_fraction,nonref_hq_count_fraction,filter\n";
	  my %lines = ();
	  open(IN,$in_file) or die $!;
	  while(<IN>)
		{ next if ($. == 1);
		  my $line = $_; chomp($line);
		  my @line = split(/\,/,$line);
		  my $bc = $line[0]; my $chr = $line[1]; my $pos = $line[2];
		  if (($chr eq $dnm_chr) and ($pos == $dnm_pos))
			{ $lines{$bc} = $line; }
		}
	  close(IN) or die $!;
	  foreach my $bc (@{$print_order{$family}})
		{ my $line = $lines{$bc};
		  my $code = $lookup{$family}{$bc};
		  print OUT "$code,$line\n";
		}
	  close(OUT) or die $!;
	  
	  # print output for downsampled BAM
	  $out_file = "$out_dir/$family-downsampled.csv";
	  open(OUT,'>',$out_file) or die $!;
	  print OUT "label,group,avg_depth_to_which_reads_downsampled,seed,sample,chr,pos,target_id,target_type,ref,alts,probes,count_hq_A,count_hq_C,count_hq_G,count_hq_T,count_hq_INS,count_hq_DEL,count_hq_SKIP,number_called_hq,ref_hq_count_fraction,nonref_hq_count_fraction,DEL/ALL,count_hq_A.ALLreplicates,count_hq_C.ALLreplicates,count_hq_G.ALLreplicates,count_hq_T.ALLreplicates,count_hq_INS.ALLreplicates,count_hq_DEL.ALLreplicates,count_hq_SKIP.ALLreplicates\n";
	  
  	  my %counts_across_all_replicates = (); # i.e. summed across all members of "$group, depth $num, seed $seed"
	  foreach my $num (@fold)
		{ for(my $x=0;$x<@seeds;$x++)
			{ my $seed_num = $x+1; my $seed = $seeds[$x];
			  my $in_file = "$in_dir/downsampled/$num/$seed/analysis/pileups/pileups_long_detailed.csv";
	          next if (!(-e($in_file)));
			  print "1. $family, depth $num, seed $seed\n";
			  my %lines = ();
			  open(IN,$in_file) or die $!;
			  while(<IN>)
				{ next if ($. == 1);
				  my $line = $_; chomp($line);
				  my @line = split(/\,/,$line);
				  my $bc = $line[0]; my $chr = $line[1]; my $pos = $line[2]; my $count_hq_A = $line[8]; my $count_hq_C = $line[9]; my $count_hq_G = $line[10]; my $count_hq_T = $line[11]; my $count_hq_INS = $line[12]; my $count_hq_DEL = $line[13]; my $count_hq_SKIP = $line[14];
				  my $code  = $lookup{$family}{$bc};
				  my $group = 'NA';
				  if ($code =~ /^(.*?)Co\d+\-\d+\-I{1,3}$/) { $group = "$1"."Co"; } elsif ($code =~ /^(.*?)\-I{1,3}$/) { $group = $1; }
				  $counts_across_all_replicates{$group}{$num}{$chr}{$pos}{$seed}{'A'}    += $count_hq_A;
				  $counts_across_all_replicates{$group}{$num}{$chr}{$pos}{$seed}{'C'}    += $count_hq_C;
				  $counts_across_all_replicates{$group}{$num}{$chr}{$pos}{$seed}{'G'}    += $count_hq_G;
				  $counts_across_all_replicates{$group}{$num}{$chr}{$pos}{$seed}{'T'}    += $count_hq_T;
				  $counts_across_all_replicates{$group}{$num}{$chr}{$pos}{$seed}{'INS'}  += $count_hq_INS;
				  $counts_across_all_replicates{$group}{$num}{$chr}{$pos}{$seed}{'DEL'}  += $count_hq_DEL;
				  $counts_across_all_replicates{$group}{$num}{$chr}{$pos}{$seed}{'SKIP'} += $count_hq_SKIP;
				}
			  close(IN) or die $!;
			}
		}
	  
	  foreach my $num (@fold)
		{ for(my $x=0;$x<@seeds;$x++)
			{ my $seed_num = $x+1; my $seed = $seeds[$x];
			  my $in_file = "$in_dir/downsampled/$num/$seed/analysis/pileups/pileups_long_detailed.csv";
	          next if (!(-e($in_file)));
			  print "2. $family, depth $num, seed $seed\n";
			  my %lines = ();
			  open(IN,$in_file) or die $!;
			  while(<IN>)
				{ next if ($. == 1);
				  my $line = $_; chomp($line);
				  my @line = split(/\,/,$line);
				  pop(@line); # remove because we don't need 'filter'
				  my $bc = $line[0]; my $chr = $line[1]; my $pos = $line[2]; my $count_hq_DEL = $line[13]; my $number_called_hq = $line[15];
				  my $DEL_ALL = 'NA';
				  if ($number_called_hq > 0) { $DEL_ALL = sprintf("%.2f",(($count_hq_DEL/$number_called_hq)*100)); }
				  push(@line,$DEL_ALL);
				  $line = join(",",@line);
				  if (($chr eq $dnm_chr) and ($pos == $dnm_pos))
					{ $lines{$bc} = $line; }
				}
			  close(IN) or die $!;
			  foreach my $bc (@{$print_order{$family}})
				{ my $line  = 'NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA';
				  if (exists($lines{$bc})) { $line = $lines{$bc}; }
				  my $code  = $lookup{$family}{$bc};
				  my $group = 'NA';
				  if ($code =~ /^(.*?)Co\d+\-\d+\-I{1,3}$/) { $group = "$1"."Co"; } elsif ($code =~ /^(.*?)\-I{1,3}$/) { $group = $1; }
				  my $count_A_all    = 'NA'; if (exists($counts_across_all_replicates{$group}{$num}{$dnm_chr}{$dnm_pos}{$seed}{'A'}))    { $count_A_all    = $counts_across_all_replicates{$group}{$num}{$dnm_chr}{$dnm_pos}{$seed}{'A'};    }
				  my $count_C_all    = 'NA'; if (exists($counts_across_all_replicates{$group}{$num}{$dnm_chr}{$dnm_pos}{$seed}{'C'}))    { $count_C_all    = $counts_across_all_replicates{$group}{$num}{$dnm_chr}{$dnm_pos}{$seed}{'C'};    }
				  my $count_G_all    = 'NA'; if (exists($counts_across_all_replicates{$group}{$num}{$dnm_chr}{$dnm_pos}{$seed}{'G'}))    { $count_G_all    = $counts_across_all_replicates{$group}{$num}{$dnm_chr}{$dnm_pos}{$seed}{'G'};    }
				  my $count_T_all    = 'NA'; if (exists($counts_across_all_replicates{$group}{$num}{$dnm_chr}{$dnm_pos}{$seed}{'T'}))    { $count_T_all    = $counts_across_all_replicates{$group}{$num}{$dnm_chr}{$dnm_pos}{$seed}{'T'};    }
				  my $count_INS_all  = 'NA'; if (exists($counts_across_all_replicates{$group}{$num}{$dnm_chr}{$dnm_pos}{$seed}{'INS'}))  { $count_INS_all  = $counts_across_all_replicates{$group}{$num}{$dnm_chr}{$dnm_pos}{$seed}{'INS'};  }
				  my $count_DEL_all  = 'NA'; if (exists($counts_across_all_replicates{$group}{$num}{$dnm_chr}{$dnm_pos}{$seed}{'DEL'}))  { $count_DEL_all  = $counts_across_all_replicates{$group}{$num}{$dnm_chr}{$dnm_pos}{$seed}{'DEL'};  }
				  my $count_SKIP_all = 'NA'; if (exists($counts_across_all_replicates{$group}{$num}{$dnm_chr}{$dnm_pos}{$seed}{'SKIP'})) { $count_SKIP_all = $counts_across_all_replicates{$group}{$num}{$dnm_chr}{$dnm_pos}{$seed}{'SKIP'}; }
				  print OUT "$code,$group,$num,$seed_num,$line,$count_A_all,$count_C_all,$count_G_all,$count_T_all,$count_INS_all,$count_DEL_all,$count_SKIP_all\n";
				}
			}
		}
      close(OUT) or die $!;
	}

exit 1;
