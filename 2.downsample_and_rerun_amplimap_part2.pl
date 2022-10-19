=head

BEFORE USAGE:

module add python-base/3.6.10
module load python-cbrg
module add amplimap/0.4.19
module add samtools/1.5
module add bcftools/1.5

=cut

use strict;
use warnings;
use Acme::Tools qw(median);

# REQUIREMENTS
my $home 		 = '/t1-data/project/GorielyLab2021/sbush/pregcare';
my $in_dir_fam27 = '/t1-data/project/GorielyLab2021/shared/mgiebler/200303_Analysis/FAM_027/reads_in';
my $in_dir_fam34 = '/t1-data/project/GorielyLab2021/shared/mgiebler/200303_Analysis/FAM_034/reads_in';
my $picard		 = '/t1-data/project/GorielyLab2021/sbush/programs/picard.jar';
my $config		 = "$home/amplimap_prerequisites/config_for_downsampling.yaml";
my $probes_fam27 = "$home/amplimap_prerequisites/FAM27/probes.csv";
my $target_fam27 = "$home/amplimap_prerequisites/FAM27/targets.csv";
my $probes_fam34 = "$home/amplimap_prerequisites/FAM34/probes.csv";
my $target_fam34 = "$home/amplimap_prerequisites/FAM34/targets.csv";
my $lookup_table = "$home/miseq_fam_to_code_lookup.txt";
my $fatal = 0;
if (!(-d($in_dir_fam27))) { $fatal++; print "ERROR: cannot find $in_dir_fam27\n"; }
if (!(-d($in_dir_fam34))) { $fatal++; print "ERROR: cannot find $in_dir_fam34\n"; }
if (!(-e($probes_fam27))) { $fatal++; print "ERROR: cannot find $probes_fam27\n"; }
if (!(-e($target_fam27))) { $fatal++; print "ERROR: cannot find $target_fam27\n"; }
if (!(-e($probes_fam34))) { $fatal++; print "ERROR: cannot find $probes_fam34\n"; }
if (!(-e($target_fam34))) { $fatal++; print "ERROR: cannot find $target_fam34\n"; }
if (!(-e($lookup_table))) { $fatal++; print "ERROR: cannot find $lookup_table\n"; }
if (!(-e($picard))) 	  { $fatal++; print "ERROR: cannot find $picard\n"; 	  }
if (!(-e($config))) 	  { $fatal++; print "ERROR: cannot find $config\n"; 	  }
exit 1 if ($fatal > 0);

# PARAMETERS
my $num_procs = 10;
my @families  = (qw/FAM27 FAM34/);
my @fold 	  = (qw/25 100 500 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 15000 20000/);
my @seeds 	  = (qw/123 456 789 234 567 891 321 654 987 432/);
my $min_fold  = 25;

# OUTPUT
my $sh_file = "$home/run_downsampler_and_amplimap_part2.sh";
open(SH,'>',$sh_file) or die $!;
print SH "#!/bin/bash\n";

# STORE LOOKUP TABLE
my %lookup = ();
open(IN,$lookup_table) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $family = $line[0]; my $bc = $line[1]; my $code = $line[2];
	  $lookup{$family}{$bc} = $code;
	}
close(IN) or die $!;

# FOR EACH FAMILY, LOOK UP WHICH BARCODES WE REQUIRE, DOWNSAMPLE, AND THEN RE-RUN AMPLIMAP
foreach my $family (@families)
	{ print "$family\n";
	  my $in_dir; my $probes; my $targets;
	  if 	($family eq 'FAM27') { $in_dir = $in_dir_fam27; $probes = $probes_fam27; $targets = $target_fam27; }
	  elsif ($family eq 'FAM34') { $in_dir = $in_dir_fam34; $probes = $probes_fam34; $targets = $target_fam34; }
	  my $out_dir = "$home/$family";
	  
	  # obtain coords of target region
	  my $region;
	  open(IN,$targets) or die $!;
	  while(<IN>)
		{ next if ($. == 1);
		  my $line = $_; chomp($line);
		  my @line = split(/\,/,$line);
		  my $chr = $line[0]; my $start = $line[1]; my $end = $line[2];
		  $region = "$chr:$start-$end";
		}
	  close(IN) or die $!;
	  if (!(defined($region))) { print "ERROR: unable to define target region for $family\n"; exit 1; }
	  
	  # for each Amplimap BAM, run samtools depth so that we can derive an average depth across all positions
	  opendir(DIR,"$out_dir/raw/analysis/bams") or die $!;
	  my @files = readdir(DIR);
	  closedir(DIR) or die $!;
	  my @sorted_files = sort {$a cmp $b} @files;
	  foreach my $file (@sorted_files)
		{ if ($file =~ /^(.*?)\.bam$/)
			{ my $bc = $1;
			  if (!(-e("$out_dir/raw/analysis/bams/$bc.depth")))
				{ print SH "samtools depth -d 0 -r $region $out_dir/raw/analysis/bams/$file > $out_dir/raw/analysis/bams/$bc.depth\n";
				}
			}
		}
	  
	  # for each Amplimap BAM, derive an average (median) depth across all positions
	  my %depth_per_bc = ();
	  foreach my $file (@sorted_files)
		{ if ($file =~ /^(.*?)\.depth$/)
			{ my $bc = $1;
			  my @depth = ();
			  open(IN,"$out_dir/raw/analysis/bams/$bc.depth") or die $!;
			  while(<IN>)
				{ my $line = $_; chomp($line);
				  my @line = split(/\t/,$line);
				  my $depth = $line[2];
				  push(@depth,$depth);
				}
			  close(IN) or die $!;
			  my $median_depth = median(@depth);
			  $depth_per_bc{$bc} = $median_depth;
			}
		}
	  
	  # for each BAM, downsample reads 10 times and then re-run Amplimap
	  # we'll downsample to x-fold coverage per base, across the target region
	  # IMPORTANT: note that not every Amplimap BAM is bona fide, because there is one BAM per sequencing barcode and not all barcodes are necessary for each experiment
	  foreach my $num (@fold)
		{ if (!(-d("$out_dir/downsampled/$num"))) 					 	   { print SH "mkdir $out_dir/downsampled/$num\n"; 						}
		  foreach my $seed (@seeds)
			{ if (!(-d("$out_dir/downsampled/$num/$seed"))) 		 	   { print SH "mkdir $out_dir/downsampled/$num/$seed\n"; 		        }
			  if (!(-d("$out_dir/downsampled/$num/$seed/mapped_bams_in"))) { print SH "mkdir $out_dir/downsampled/$num/$seed/mapped_bams_in\n"; }
			  if (!(-e("$out_dir/downsampled/$num/$seed/analysis/pileups/pileups_long_detailed.csv")))
				{ my @bcs = ();
				  while((my $bc,my $irrel)=each(%{$lookup{$family}}))
					{ push(@bcs,$bc); }
				  my @sorted_bcs = sort {$a cmp $b} @bcs;
				  foreach my $bc (@sorted_bcs)
					{ next if (!(exists($lookup{$family}{$bc}))); # CHECKPOINT: there is no point running this procedure on barcodes that haven't been used for the experiment
					  my $bam     = "$out_dir/raw/analysis/bams/$bc.bam";
					  my $out_bam = "$out_dir/downsampled/$num/$seed/mapped_bams_in/$bc.bam";
					  
					  # we will downsample to $num/$median_depth reads, e.g. if average read depth is 50,000 and we want to downsample to 1000x then the probability of keeping any given read will be 1000/50000 = 0.02. In this case, then, we are downsampling to 2% of the original total
					  my $median_depth  = $depth_per_bc{$bc};
					  next if ($median_depth < $min_fold);
					  my $downsample_to = ($num/$median_depth);
					  next if (($downsample_to < 0) or ($downsample_to > 1));
					  print "$family. downsample to: $num-fold. seed: $seed. barcode: $bc. median depth of original BAM: $median_depth. probability of keeping a read from this BAM: $downsample_to\n";
					  print SH "java -jar $picard DownsampleSam -I $bam -O $out_bam -P $downsample_to --RANDOM_SEED $seed --STRATEGY HighAccuracy --TMP_DIR $out_dir/downsampled/$num/$seed/mapped_bams_in\n";
					  print SH "samtools index $out_bam\n";
					}
				  print SH "cd $out_dir/downsampled/$num/$seed\n";
				  print SH "cp $config config.yaml\n";
				  print SH "cp $probes .\n";
				  print SH "cp $targets .\n";
				  print SH "amplimap pileups --ncores $num_procs --run\n";
				  print SH "cd $out_dir\n";
				}
			}
		}
	}
close(SH) or die $!;
exit 1;