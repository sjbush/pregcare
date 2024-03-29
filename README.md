# pregcare

This repository contains scripts used in [“Personalised recurrence risk assessment following the birth of a child with a pathogenic _de novo_ mutation”](https://www.nature.com/articles/s41467-023-36606-w) by Bernkopf, Abdullah, _et al._

# Haplotyping

The scripts used for haplotyping the long-read data, which should be run in numbered order, perform the following steps and have the following software prerequisites, which should be accessible on the command line. More specific detail on prerequisites is available as comments within the scripts themselves.

**Script 1: split_fast5_among_subdirs.sh**

The raw data for this analysis takes the form of eight ONT (Oxford Nanopore) MinION runs, in which targeted sequencing was performed for each member of a family trio (mother/father/proband). Each ONT run contains multiple families which are distinguished by sequencing barcodes. The ONT run numbers, barcodes, sequencing coordinates, and associated family IDs are available in the metadata file ont_coords.tsv. We first need to convert all ONT fast5 files to fastqs, and then demultiplex them according to the barcodes present in this file. Scripts 1 to 5 are used to do this.
To begin, this script runs within a directory containing raw fast5 files. It creates a series of subdirectories (of the format “subdirX”, where X is numeric) each of which contain up to three fast5s. The purpose of this script is to speed up downstream processing by allowing a series of smaller subdirectories to be submitted to a SLURM cluster rather than submitting just the original large directory. (To submit individual fast5s to the cluster, you could edit this script to change ‘3’ to ‘1’).

**Script 2: submit_each_subdir_to_slurm.sh**

This script submits a series of batch jobs to a SLURM cluster, each one running script 3 (‘run_guppy.sh’) with a different fast5-containing subdirectory as input.

**Script 3: run_guppy.sh**

This script runs [Guppy basecaller](https://nanoporetech.com/) with the parameter --config dna_r9.4.1_450bps_hac.cfg, and produces one fastq file per fast5 input.

**Script 4: concat_all_guppy_fqs.sh**

This script concatenates all individual Guppy fastqs into one fastq per ONT run: minion1.fq.gz, minion2.fq.gz, and so on.

**Script 5: demultiplex_nanopore_fq.sh**

This script runs Guppy barcoder on each ONT-run-specific fastq, demultiplexing it to generate one fastq per barcode, i.e. one set of reads for each family member of a trio.

**Script 6: nanostat_on_demultiplexed_fqs.sh**

This script runs the QC tool [NanoStat](https://github.com/wdecoster/nanostat) each demultiplexed fastq, producing a summary stats report. This contains, among others, total read and base counts, average read length and quality, read N50, and the number, percentage and Mb of reads above a range of quality cutoffs.

**Script 7: varcall_from_demultiplexed_fq.pl**

This script takes as input the demultiplexed fastqs created by the previous script. It produces a shell script which contains the commands to (a) align the reads and filter the ensuing BAM for quality (using [minimap2](https://github.com/lh3/minimap2) and [samtools](http://www.htslib.org/)), (b) variant call per barcode, producing a series of VCFs (using [Medaka](https://github.com/nanoporetech/medaka)), (c) concatenate the Medaka VCFs, one per barcode, to produce one final VCF which (d) is sorted, indexed and annotated with [dbSNP](https://ftp.ncbi.nih.gov/snp/latest_release/VCF/) IDs (using [bcftools](https://samtools.github.io/bcftools/bcftools.html)).

**Script 8: medaka_phaser.pl**

This script takes as input the metadata file ont_coords.tsv plus the BAMs and VCFs generated by the previous script. It calls [samtools](http://www.htslib.org/).
It parses the VCFs to first determine whether Medaka had called and phased the de novo mutation (DNM) in the proband (but not in either parent, confirming their true de novo nature). For each DNM, it then obtains the associated ‘phase set’ SNPs, those on the same haplotype, retaining only those which had a total depth of coverage > 10x. It then cross-references the phase set SNPs with the VCFs from the mother and father to identify which calls (if any) had been made at those positions. This produces a set of three haplotypes from which the script classifies the inheritance of the DNM as either maternal or paternal (the SNPs in phase with the DNM could only be derived from the chromosome inherited from the mother or father, respectively), else unresolved (Medaka either did not call the DNM in the child, called it but did not construct a phase set, or, if it did construct a phase set, either did not call its constituent SNPs in the parents or made identical calls for both of them). The output is a table detailing, per family ID, the following:

* Family ID
* _De novo_ mutation (chr:position)
* Reference allele (for DNM)
* Variant allele (for DNM)
* Read counts (for DNM, given as A|C|G|T only if called by Medaka in the proband)
* Total depth of coverage (for DNM)
* dbSNP ID (for DNM)
* Has Medaka detected the DNM in the proband? (expectation: yes)
* Has Medaka detected the DNM in the mother? (expectation: no)
* Has Medaka detected the DNM in the father? (expectation: no)
* Has Medaka phased the DNM in the proband?
* If yes, how many SNPs are in the same phase set, not including the DNM?
* SNPs in the proband's phase set, not including the DNM (listed in order of distance from DNM and reported as dbSNP IDs if available, else 'chr:pos')
* How many informative SNPs (= different between parents) are in the same phase set, not including the DNM?
* Informative SNPs (= different between parents) in the proband’s phase set, not including the DNM (listed in order of distance from DNM and reported as dbSNP IDs if available, else 'chr:pos')
* Total depth of coverage for informative SNPs in the proband's phase set, listed in the order previously stated
* Bases that phase with the DNM (in the proband, listed in the order previously stated)
* Genotypes, listed in order previously stated (mother)
* Genotypes, listed in order previously stated (father)
* Genotypes, listed in order previously stated (proband)
* Pattern of inheritance (maternal, paternal, or unresolved)
* Min. distance of an informative SNP from DNM (bp)
* Max. distance of an informative SNP from DNM (bp)
* Reason the pattern of inheritance was unresolved (if applicable). The reasons are either:
  - DNM not called in the child
  - DNM called in the child but no in-phase SNPs identified
  - DNM called in the child and in-phase SNPs identified but the workflow discards them all because none meet user-specified QC requirements
  - DNM called in the child and in-phase SNPs identified, but none are informative - [SNP ID] is either homozygous or heterozygous in both parents

**Script 9: pileup_phaser.pl**

As with script 8, this script takes as input the metadata file ont_coords.tsv plus the BAMs generated by script 7. It calls [samtools](http://www.htslib.org/). This script implements an alternative pileup-based approach to phasing. It is used when a DNM could not be successfully phased using Medaka and implements the following approach:

* Using samtools mpileup with parameters -aa --output-QNAME, produce a ‘pileup’ string (format described [here](http://www.htslib.org/doc/samtools-mpileup.html)) and an associated list of read IDs covering the DNM locus in the child BAM.
* Parse the pileup and read ID strings to generate an associative array: a one-to-one correspondence of calls and read IDs.
* Summarise the number and proportion of reads supporting the DNM REF and ALT calls in the child. We require that the DNM locus be covered by > 50 reads, and consider homozygous reference calls to be those which are < 30% ALT, heterozygous calls to be 70% > ALT ≥ 30%, and homozygous variant calls to be ≥ 70% ALT. We considered a DNM to not be called if it was not heterozygous according to these criteria.
* Extract each REF- and ALT-containing read from the child BAM and identify, for the latter, the coordinates of the genomic region this set of alignments covers.
* Parse dbSNP to obtain a candidate list of SNPs which have previously been found within this region.
* After sorting this list in descending order of distance from the DNM, iterate through each candidate and determine its genotype in each parent and in child. We discard candidate SNPs which were not heterozygous in the child and at least one parent, according to the depth and percentage requirements described above. We also required that the minimum depth threshold be met by all three individuals.
* At this point, we have a set of reads which call the ALT-allele for the DNM and a candidate phasing SNP (the closest one to it and for which there is a prior, given its inclusion in dbSNP). The script then constructs a 2x2 count table, the rows of which are ‘number of reads calling REF at DNM position’ and ‘number of reads calling ALT at DNM position’ and columns ‘number of reads calling REF at phasing SNP position’ and ‘number of reads calling ALT at phasing SNP position.’
* Inheritance of the DNM is then resolved by considering which of the two alleles for the phasing SNP, REF or ALT, were disproportionately found on the same read as the DNM ALT.

The output of this script is the 2x2 count table (from which we perform Fisher’s exact test) and a table detailing, per family ID, the contents of each DNM’s phase set and predicted inheritance. The format of this table is as described for script 8.

# Downsampling

The scripts in this subdirectory are used for downsampling the short-read (MiSeq) data for two mosaic families with low-frequency de novo mutations (FAM27 and FAM34). Scripts should be run in numbered order and collectively perform the following analysis.

As described in [the paper](https://www.biorxiv.org/content/10.1101/2022.07.26.501520v1), the MiSeq data were originally analyzed using [Amplimap](https://github.com/koelling/amplimap) to obtain both the VAF of each family-specific mutation and the total count of >Q30 bases at the corresponding genomic position in each PCR replicate and sample.

For each replicate and sample, Amplimap produces a BAM (**script 1**). Using samtools ‘depth’ with parameters -d 0 -r, we determined the median sequencing depth across all positions of the target region (**script 2**).

We then downsampled each BAM to x-fold coverage per base, where x was 25, 100, 500, and 1000, every multiple of 1000 to 10,000, plus 15,000 and 20,000 (i.e. downsampling to 15 different depths) (**script 2**).

Downsampling was performed 10 times per depth using [Picard](https://broadinstitute.github.io/picard/) DownsampleSam with parameters --STRATEGY HighAccuracy and --P, where P (the probability of retaining a given read when iterating through the BAM) was equal to desired fold depth / median sequencing depth. Samples where P was either not a number or a number > 1 (indicative of a failed, low- or zero-depth, sample) were discarded. To facilitate reproducibility, seeds were not randomly generated but manually assigned: 123, 456, 789, 234, 567, 891, 321, 654, 987, and 432.

Finally, for each PCR replicate, sample, fold depth and seed, we re-ran Amplimap ‘pileups’ but in ‘mapped_bams_in’ mode (**script 2**).

This produced, for each PCR replicate, sample, fold depth and seed, an Amplimap output file, pileups_long_detailed.csv. These files were combined into a single output file per family (script 3), from which we corrected VAFs using **scripts 4 and 5**.

To visualise this data, boxplots were created using **script 6**. This reproduces Supplementary Figure S3 of the paper.
