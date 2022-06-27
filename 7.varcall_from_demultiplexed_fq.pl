=head
BEFORE USAGE:

# make reference index
cd /t1-data/project/GorielyLab2021/sbush/genomes
wget http://ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
/home/s/sbush/programs/minimap2-2.18_x64-linux/minimap2 -x map-ont -I 100G -t 10 -d Homo_sapiens.GRCh38.dna.primary_assembly.idx Homo_sapiens.GRCh38.dna.primary_assembly.fa

# download dbSNP SNPs
curl https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.38.gz --output dbSNP_153.hg38.vcf.gz
curl https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.38.gz.tbi --output dbSNP_153.hg38.vcf.gz.tbi

# we need to rename dbSNP chromosomes as those in dbSNP are given RefSeq IDs whereas those in the Medaka VCFs have the more straightforward "1", "2", etc. To do this we need the RefSeq assembly report. See https://www.ncbi.nlm.nih.gov/assembly/GCA_000001405.28; more specifically, https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_report.txt
# make a two-column old\tnew file from this and use as input to bcftools annotate
bcftools annotate --threads 15 --rename-chrs RefSeqAcc_to_Chr_lookup.txt dbSNP_153.hg38.vcf.gz | bgzip > dbSNP_153.hg38.renamed.vcf.gz
mv dbSNP_153.hg38.renamed.vcf.gz dbSNP_153.hg38.vcf.gz
tabix dbSNP_153.hg38.vcf.gz

=cut

use strict;
use warnings;

# REQUIREMENTS
my $homedir = '/home/s/sbush';
my $basedir = '/t1-data/project/pregcare/sbush';
my @ont_dir = (qw/MinION1 MinION2 MinION3 MinION4 MinION5 MinION6 MinION7 MinION8/);
my $minimap = '/t1-data/project/GorielyLab2021/sbush/programs/minimap2-2.18_x64-linux/minimap2';
my $hs_idx  = '/t1-data/project/GorielyLab2021/sbush/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.idx'; # for indexing commands, see above
my $hs_fa   = '/t1-data/project/GorielyLab2021/sbush/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa';
my $dbsnp   = '/t1-data/project/pregcare/sbush/dbSNP/dbSNP_153.hg38.vcf.gz'; # for source of this file, see above
my $primers = "$basedir/ont_coords.tsv";
my $fatal   = 0;
if (!(-e($minimap))) { $fatal++; print "ERROR: cannot find $minimap\n"; }
if (!(-e($hs_idx)))  { $fatal++; print "ERROR: cannot find $hs_idx\n";  }
if (!(-e($hs_fa)))   { $fatal++; print "ERROR: cannot find $hs_fa\n";   }
if (!(-e($dbsnp)))   { $fatal++; print "ERROR: cannot find $dbsnp\n";   }
if (!(-e($primers))) { $fatal++; print "ERROR: cannot find $primers\n"; }
exit 1 if ($fatal > 0);

# PARAMETERS
my $num_procs = 10;

# OUTPUT
my $sh_file = "run_ont_varcaller.sh";
open(SH,'>',$sh_file) or die $!;
print SH "#!/bin/bash\n";

# STORE COORDS OF EACH TARGETED SEQUENCING REGION
my %coords = ();
open(IN,$primers) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line); $line =~ s/\r$//g; $line =~ s/\s$//g;
	  my @line = split(/\t/,$line);
	  my $ont_line = $line[0]; my $bc_line = $line[1]; my $fam = $line[2]; my $coords = $line[3]; $coords =~ s/\s$//g;
	  my @ont_line = split(/\|/,$ont_line);
	  my @bc_line  = split(/\|/,$bc_line);
	  for(my $x=0;$x<@bc_line;$x++)
		{ my $bc = $bc_line[$x];
		  my $ont = $ont_line[$x];
		  next if ($bc !~ /^BC\d+$/);
		  my $chr; my $start; my $end;
		  if ($coords =~ /^(.*?)\:(\d+)\-(\d+)$/)
			{ $chr = $1; $start = $2; $end = $3; }
		  if ( (!(defined($start))) or (!(defined($end))) )
			{ print "ERROR in $primers, ONT $ont $bc_line, line $. --> $coords\n"; exit 1; }
		  if ($start > $end)
			{ print "ERROR in $coords, ONT $ont $bc_line, line $. --> $coords --> $start > $end\n"; exit 1; }
		  $coords{$ont}{$bc}{$coords}++;
		}
	}
close(IN) or die $!;

foreach my $ont (@ont_dir)
	{ next if (!(-d("$basedir/$ont")));
	  if (!(exists($coords{$ont})))
		{ print "ERROR: unable to find primer coords for $ont\n"; }
	  next if (!(exists($coords{$ont})));
	  my $fq_dir  = "$basedir/$ont/5.final_fq"; # from 6.nanostat_on_demultiplexed_fqs.sh
	  my $bam_dir = "$basedir/$ont/6.bam";
	  my $vcf_dir = "$basedir/$ont/7.vcf";	
	  next if (!(-d($fq_dir)));
	  if (!(-d($bam_dir))) { mkdir $bam_dir or die $!; }
	  if (!(-d($vcf_dir))) { mkdir $vcf_dir or die $!; }
	  opendir(FQ,$fq_dir) or die $!;
	  my @fqs = readdir(FQ);
	  closedir(FQ) or die $!;
	  my @sorted_fqs = sort {$a cmp $b} @fqs;
	  foreach my $fq (@sorted_fqs) # one fq per barcode
		{ next if (($fq eq '.') or ($fq eq '..'));
		  next if ($fq !~ /^(.*?)\.fq\.gz$/);
		  my $bc = '';
		  if ($fq =~ /^(.*?)\.fq\.gz$/) { $bc = $1; }
		  next if ($bc eq '');
		  if (!(exists($coords{$ont}{$bc})))
			{ print "ERROR: unable to find primer coords for $ont, barcode $bc\n"; }
		  next if (!(exists($coords{$ont}{$bc})));
		  
		  # align reads to reference, outputting a sorted, indexed, BAM
		  print SH "cd $bam_dir\n";
		  my $bam = "$bam_dir/$bc.bam";
		  if (!(-e($bam)))
			{ print SH "$minimap -ax map-ont --sam-hit-only -t $num_procs $hs_idx $fq_dir/$fq | samtools view -h -b -q 20 -F 256 -F 2048 - | samtools sort -T tmp -o $bam -\n"; }
		  if (!(-e("$bam.bai")))
			{ print SH "samtools index $bam\n"; }
		  
		  # identify the coords of the target regions
		  my @coords = ();
		  while((my $coords,my $irrel)=each(%{$coords{$ont}{$bc}}))
			{ push(@coords,$coords); }
		  my @sorted_coords = sort {$a cmp $b} @coords;
		  
		  # move to the output directory and prepare to make $final_vcf
		  # we do this by calling and phasing using Medaka for each targeted region, following recommended practice: https://labs.epi2me.io/notebooks/Human_Variant_Calling_with_Medaka.html
		  my $final_vcf = "$vcf_dir/$bc.vcf";
		  print SH "cd $vcf_dir\n";
		  if (!(-d("$vcf_dir/$bc"))) { print SH "mkdir $vcf_dir/$bc\n"; }
		  print SH "cd $vcf_dir/$bc\n";
		  my $vcf_line = ''; my $called_something = 0;
		  foreach my $coords (@sorted_coords)
			{ if ( (!(-e("$vcf_dir/$bc/$coords.vcf"))) or (!(-d("$vcf_dir/$bc/$coords"))) )
				{ print SH "medaka_variant -f $hs_fa -i $bam -o $coords -t $num_procs -r $coords -p -d\n";
				  print SH "if [ -e $vcf_dir/$bc/$coords/round_1_phased.vcf ]; then cp $vcf_dir/$bc/$coords/round_1_phased.vcf $vcf_dir/$bc/$coords.vcf; else cp $vcf_dir/$bc/$coords/round_1.vcf $vcf_dir/$bc/$coords.vcf; fi\n";
				  $called_something++;
				}
			  $vcf_line .= "$vcf_dir/$bc/$coords.vcf ";
			}
		  if ( (!(-e($final_vcf))) or ($called_something > 0) )
			{ # concatenate the set of Medaka VCFs
		      my $concat_vcf 	  		  = "$vcf_dir/$bc.concat.vcf";
			  my $concat_sort_vcf 		  = "$vcf_dir/$bc.concat_sorted.vcf";
			  my $concat_sort_vcf_with_dp = "$vcf_dir/$bc.concat_sorted_with_DP.vcf";
			  print SH "bcftools concat $vcf_line > $concat_vcf\n";
			  print SH "bcftools sort $concat_vcf > $concat_sort_vcf\n";
			  print SH "rm $concat_vcf\n";
			  
			  # add read depth to concatenated VCF, then sort and index
			  print SH "medaka tools annotate --dpsp $concat_sort_vcf $hs_fa $bam $concat_sort_vcf_with_dp\n";
			  print SH "rm $concat_sort_vcf\n";
			  print SH "bgzip $concat_sort_vcf_with_dp\n";
			  print SH "tabix $concat_sort_vcf_with_dp.gz\n";
			  
			  # annotate concatenated+depth VCF with known SNP IDs, then sort and index. This is the final VCF.
			  print SH "bcftools annotate --threads $num_procs -a $dbsnp -c ID -o $final_vcf $concat_sort_vcf_with_dp.gz\n";
			  print SH "rm $concat_sort_vcf_with_dp.gz $concat_sort_vcf_with_dp.gz.tbi\n";
			}
		}
	}

close(SH) or die $!;
exit 1;