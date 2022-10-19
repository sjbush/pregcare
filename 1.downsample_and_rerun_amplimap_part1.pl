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

# REQUIREMENTS
my $home 		 = '/t1-data/project/GorielyLab2021/sbush/pregcare';
my $in_dir_fam27 = '/t1-data/project/GorielyLab2021/shared/mgiebler/200303_Analysis/FAM_027/reads_in';
my $in_dir_fam34 = '/t1-data/project/GorielyLab2021/shared/mgiebler/200303_Analysis/FAM_034/reads_in';
my $config		 = "$home/amplimap_prerequisites/config.yaml";
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
if (!(-e($config))) 	  { $fatal++; print "ERROR: cannot find $config\n"; 	  }
exit 1 if ($fatal > 0);

# PARAMETERS
my $num_procs = 10;
my @families  = (qw/FAM27 FAM34/);

# OUTPUT
my $sh_file = "$home/run_downsampler_and_amplimap_part1.sh";
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

# FOR EACH FAMILY, LOOK UP WHICH BARCODES WE REQUIRE AND THEN RUN AMPLIMAP
foreach my $family (@families)
	{ print "$family\n";
	  my $in_dir; my $probes; my $targets;
	  if 	($family eq 'FAM27') { $in_dir = $in_dir_fam27; $probes = $probes_fam27; $targets = $target_fam27; }
	  elsif ($family eq 'FAM34') { $in_dir = $in_dir_fam34; $probes = $probes_fam34; $targets = $target_fam34; }
	  
	  # make an output dir
	  my $out_dir = "$home/$family";
	  if (!(-d($out_dir))) 		 		  { print SH "mkdir $out_dir\n"; 			  }
	  if (!(-d("$out_dir/raw"))) 		  { print SH "mkdir $out_dir/raw\n"; 		  }
	  if (!(-d("$out_dir/raw/reads_in"))) { print SH "mkdir $out_dir/raw/reads_in\n"; }
	  if (!(-d("$out_dir/downsampled")))  { print SH "mkdir $out_dir/downsampled\n";  }
	  print SH "cd $out_dir\n";
	  
	  # iterate through input dir and act only on those fastqs corresponding to the appropriate barcodes
	  opendir(DIR,$in_dir) or die $!;
	  my @files = readdir(DIR);
	  closedir(DIR) or die $!;
	  foreach my $file (@files)
		{ next if (($file eq '.') or ($file eq '..'));
		  if ($file =~ /^(.*?)\_L.+?\.fastq\.gz$/)
			{ my $bc = $1;
			  if (exists($lookup{$family}{$bc}))
				{ # copy this file into a storage directory, if we have not already done so
				  print SH "cp $in_dir/$file $out_dir/raw/reads_in/$file\n" unless (-e("$out_dir/raw/reads_in/$file"));
				}
			}
		}
	  
	  # run Amplimap on the original data
	  if (!(-e("$out_dir/raw/analysis/pileups/pileups_long_detailed.csv")))
		{ print SH "cd $out_dir/raw\n";
		  print SH "cp $config .\n";
		  print SH "cp $probes .\n";
		  print SH "cp $targets .\n";
		  print SH "amplimap pileups variants --ncores $num_procs --run\n";
		  print SH "cd $out_dir\n";
		}
	}
close(SH) or die $!;
exit 1;