#!/bin/bash
mkdir /t1-data/project/pregcare/sbush/MinION1/3.all_ont_fq
find /t1-data/project/pregcare/sbush/MinION1/2.fast5_to_fastq/subdir*/pass -name "*.fastq.gz" -exec cat {} > /t1-data/project/pregcare/sbush/MinION1/3.all_ont_fq/minion1.fq.gz \;
mkdir /t1-data/project/pregcare/sbush/MinION2/3.all_ont_fq
find /t1-data/project/pregcare/sbush/MinION2/2.fast5_to_fastq/subdir*/pass -name "*.fastq.gz" -exec cat {} > /t1-data/project/pregcare/sbush/MinION2/3.all_ont_fq/minion2.fq.gz \;
mkdir /t1-data/project/pregcare/sbush/MinION3/3.all_ont_fq
find /t1-data/project/pregcare/sbush/MinION3/2.fast5_to_fastq/subdir*/pass -name "*.fastq.gz" -exec cat {} > /t1-data/project/pregcare/sbush/MinION3/3.all_ont_fq/minion3.fq.gz \;
mkdir /t1-data/project/pregcare/sbush/MinION4/3.all_ont_fq
find /t1-data/project/pregcare/sbush/MinION4/2.fast5_to_fastq/subdir*/pass -name "*.fastq.gz" -exec cat {} > /t1-data/project/pregcare/sbush/MinION4/3.all_ont_fq/minion4.fq.gz \;
mkdir /t1-data/project/pregcare/sbush/MinION5/3.all_ont_fq
find /t1-data/project/pregcare/sbush/MinION5/2.fast5_to_fastq/subdir*/pass -name "*.fastq.gz" -exec cat {} > /t1-data/project/pregcare/sbush/MinION5/3.all_ont_fq/minion5.fq.gz \;
mkdir /t1-data/project/pregcare/sbush/MinION6/3.all_ont_fq
find /t1-data/project/pregcare/sbush/MinION6/2.fast5_to_fastq/subdir*/pass -name "*.fastq.gz" -exec cat {} > /t1-data/project/pregcare/sbush/MinION6/3.all_ont_fq/minion6.fq.gz \;
mkdir /t1-data/project/pregcare/sbush/MinION7/3.all_ont_fq
find /t1-data/project/pregcare/sbush/MinION7/2.fast5_to_fastq/subdir*/pass -name "*.fastq.gz" -exec cat {} > /t1-data/project/pregcare/sbush/MinION7/3.all_ont_fq/minion7.fq.gz \;
mkdir /t1-data/project/pregcare/sbush/MinION8/3.all_ont_fq
find /t1-data/project/pregcare/sbush/MinION8/2.fast5_to_fastq/subdir*/pass -name "*.fastq.gz" -exec cat {} > /t1-data/project/pregcare/sbush/MinION8/3.all_ont_fq/minion8.fq.gz \;
