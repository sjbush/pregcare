#!/bin/bash
mkdir /t1-data/project/pregcare/sbush/MinION2/5.final_fq
cp /t1-data/project/pregcare/sbush/MinION2/4.demultiplexed_fq/barcode01/*.fastq.gz /t1-data/project/pregcare/sbush/MinION2/5.final_fq/BC1.fq.gz
cp /t1-data/project/pregcare/sbush/MinION2/4.demultiplexed_fq/barcode02/*.fastq.gz /t1-data/project/pregcare/sbush/MinION2/5.final_fq/BC2.fq.gz
cp /t1-data/project/pregcare/sbush/MinION2/4.demultiplexed_fq/barcode03/*.fastq.gz /t1-data/project/pregcare/sbush/MinION2/5.final_fq/BC3.fq.gz
cp /t1-data/project/pregcare/sbush/MinION2/4.demultiplexed_fq/barcode04/*.fastq.gz /t1-data/project/pregcare/sbush/MinION2/5.final_fq/BC4.fq.gz
cp /t1-data/project/pregcare/sbush/MinION2/4.demultiplexed_fq/barcode05/*.fastq.gz /t1-data/project/pregcare/sbush/MinION2/5.final_fq/BC5.fq.gz
cp /t1-data/project/pregcare/sbush/MinION2/4.demultiplexed_fq/barcode06/*.fastq.gz /t1-data/project/pregcare/sbush/MinION2/5.final_fq/BC6.fq.gz
cp /t1-data/project/pregcare/sbush/MinION2/4.demultiplexed_fq/barcode07/*.fastq.gz /t1-data/project/pregcare/sbush/MinION2/5.final_fq/BC7.fq.gz
cp /t1-data/project/pregcare/sbush/MinION2/4.demultiplexed_fq/barcode08/*.fastq.gz /t1-data/project/pregcare/sbush/MinION2/5.final_fq/BC8.fq.gz
cp /t1-data/project/pregcare/sbush/MinION2/4.demultiplexed_fq/barcode09/*.fastq.gz /t1-data/project/pregcare/sbush/MinION2/5.final_fq/BC9.fq.gz
cp /t1-data/project/pregcare/sbush/MinION2/4.demultiplexed_fq/barcode10/*.fastq.gz /t1-data/project/pregcare/sbush/MinION2/5.final_fq/BC10.fq.gz
cp /t1-data/project/pregcare/sbush/MinION2/4.demultiplexed_fq/barcode11/*.fastq.gz /t1-data/project/pregcare/sbush/MinION2/5.final_fq/BC11.fq.gz
cp /t1-data/project/pregcare/sbush/MinION2/4.demultiplexed_fq/barcode12/*.fastq.gz /t1-data/project/pregcare/sbush/MinION2/5.final_fq/BC12.fq.gz
NanoStat --fastq BC1.fq.gz > BC1.stats
NanoStat --fastq BC2.fq.gz > BC2.stats
NanoStat --fastq BC3.fq.gz > BC3.stats
NanoStat --fastq BC4.fq.gz > BC4.stats
NanoStat --fastq BC5.fq.gz > BC5.stats
NanoStat --fastq BC6.fq.gz > BC6.stats
NanoStat --fastq BC7.fq.gz > BC7.stats
NanoStat --fastq BC8.fq.gz > BC8.stats
NanoStat --fastq BC9.fq.gz > BC9.stats
NanoStat --fastq BC10.fq.gz > BC10.stats
NanoStat --fastq BC11.fq.gz > BC11.stats
NanoStat --fastq BC12.fq.gz > BC12.stats

mkdir /t1-data/project/pregcare/sbush/MinION3/5.final_fq
cp /t1-data/project/pregcare/sbush/MinION3/4.demultiplexed_fq/barcode01/*.fastq.gz /t1-data/project/pregcare/sbush/MinION3/5.final_fq/BC1.fq.gz
cp /t1-data/project/pregcare/sbush/MinION3/4.demultiplexed_fq/barcode02/*.fastq.gz /t1-data/project/pregcare/sbush/MinION3/5.final_fq/BC2.fq.gz
cp /t1-data/project/pregcare/sbush/MinION3/4.demultiplexed_fq/barcode03/*.fastq.gz /t1-data/project/pregcare/sbush/MinION3/5.final_fq/BC3.fq.gz
cp /t1-data/project/pregcare/sbush/MinION3/4.demultiplexed_fq/barcode05/*.fastq.gz /t1-data/project/pregcare/sbush/MinION3/5.final_fq/BC5.fq.gz
cp /t1-data/project/pregcare/sbush/MinION3/4.demultiplexed_fq/barcode06/*.fastq.gz /t1-data/project/pregcare/sbush/MinION3/5.final_fq/BC6.fq.gz
cp /t1-data/project/pregcare/sbush/MinION3/4.demultiplexed_fq/barcode07/*.fastq.gz /t1-data/project/pregcare/sbush/MinION3/5.final_fq/BC7.fq.gz
cp /t1-data/project/pregcare/sbush/MinION3/4.demultiplexed_fq/barcode09/*.fastq.gz /t1-data/project/pregcare/sbush/MinION3/5.final_fq/BC9.fq.gz
cp /t1-data/project/pregcare/sbush/MinION3/4.demultiplexed_fq/barcode10/*.fastq.gz /t1-data/project/pregcare/sbush/MinION3/5.final_fq/BC10.fq.gz
cp /t1-data/project/pregcare/sbush/MinION3/4.demultiplexed_fq/barcode11/*.fastq.gz /t1-data/project/pregcare/sbush/MinION3/5.final_fq/BC11.fq.gz
NanoStat --fastq BC1.fq.gz > BC1.stats
NanoStat --fastq BC2.fq.gz > BC2.stats
NanoStat --fastq BC3.fq.gz > BC3.stats
NanoStat --fastq BC5.fq.gz > BC5.stats
NanoStat --fastq BC6.fq.gz > BC6.stats
NanoStat --fastq BC7.fq.gz > BC7.stats
NanoStat --fastq BC9.fq.gz > BC9.stats
NanoStat --fastq BC10.fq.gz > BC10.stats
NanoStat --fastq BC11.fq.gz > BC11.stats

mkdir /t1-data/project/pregcare/sbush/MinION4/5.final_fq
cp /t1-data/project/pregcare/sbush/MinION4/4.demultiplexed_fq/barcode01/*.fastq.gz /t1-data/project/pregcare/sbush/MinION4/5.final_fq/BC1.fq.gz
cp /t1-data/project/pregcare/sbush/MinION4/4.demultiplexed_fq/barcode02/*.fastq.gz /t1-data/project/pregcare/sbush/MinION4/5.final_fq/BC2.fq.gz
cp /t1-data/project/pregcare/sbush/MinION4/4.demultiplexed_fq/barcode03/*.fastq.gz /t1-data/project/pregcare/sbush/MinION4/5.final_fq/BC3.fq.gz
cp /t1-data/project/pregcare/sbush/MinION4/4.demultiplexed_fq/barcode05/*.fastq.gz /t1-data/project/pregcare/sbush/MinION4/5.final_fq/BC5.fq.gz
cp /t1-data/project/pregcare/sbush/MinION4/4.demultiplexed_fq/barcode06/*.fastq.gz /t1-data/project/pregcare/sbush/MinION4/5.final_fq/BC6.fq.gz
cp /t1-data/project/pregcare/sbush/MinION4/4.demultiplexed_fq/barcode07/*.fastq.gz /t1-data/project/pregcare/sbush/MinION4/5.final_fq/BC7.fq.gz
NanoStat --fastq BC1.fq.gz > BC1.stats
NanoStat --fastq BC2.fq.gz > BC2.stats
NanoStat --fastq BC3.fq.gz > BC3.stats
NanoStat --fastq BC5.fq.gz > BC5.stats
NanoStat --fastq BC6.fq.gz > BC6.stats
NanoStat --fastq BC7.fq.gz > BC7.stats

mkdir /t1-data/project/pregcare/sbush/MinION5/5.final_fq
cp /t1-data/project/pregcare/sbush/MinION5/4.demultiplexed_fq/barcode01/*.fastq.gz /t1-data/project/pregcare/sbush/MinION5/5.final_fq/BC1.fq.gz
cp /t1-data/project/pregcare/sbush/MinION5/4.demultiplexed_fq/barcode02/*.fastq.gz /t1-data/project/pregcare/sbush/MinION5/5.final_fq/BC2.fq.gz
cp /t1-data/project/pregcare/sbush/MinION5/4.demultiplexed_fq/barcode03/*.fastq.gz /t1-data/project/pregcare/sbush/MinION5/5.final_fq/BC3.fq.gz
cp /t1-data/project/pregcare/sbush/MinION5/4.demultiplexed_fq/barcode04/*.fastq.gz /t1-data/project/pregcare/sbush/MinION5/5.final_fq/BC4.fq.gz
cp /t1-data/project/pregcare/sbush/MinION5/4.demultiplexed_fq/barcode05/*.fastq.gz /t1-data/project/pregcare/sbush/MinION5/5.final_fq/BC5.fq.gz
cp /t1-data/project/pregcare/sbush/MinION5/4.demultiplexed_fq/barcode06/*.fastq.gz /t1-data/project/pregcare/sbush/MinION5/5.final_fq/BC6.fq.gz
cp /t1-data/project/pregcare/sbush/MinION5/4.demultiplexed_fq/barcode07/*.fastq.gz /t1-data/project/pregcare/sbush/MinION5/5.final_fq/BC7.fq.gz
cp /t1-data/project/pregcare/sbush/MinION5/4.demultiplexed_fq/barcode08/*.fastq.gz /t1-data/project/pregcare/sbush/MinION5/5.final_fq/BC8.fq.gz
cp /t1-data/project/pregcare/sbush/MinION5/4.demultiplexed_fq/barcode09/*.fastq.gz /t1-data/project/pregcare/sbush/MinION5/5.final_fq/BC9.fq.gz
NanoStat --fastq BC1.fq.gz > BC1.stats
NanoStat --fastq BC2.fq.gz > BC2.stats
NanoStat --fastq BC3.fq.gz > BC3.stats
NanoStat --fastq BC4.fq.gz > BC4.stats
NanoStat --fastq BC5.fq.gz > BC5.stats
NanoStat --fastq BC6.fq.gz > BC6.stats
NanoStat --fastq BC7.fq.gz > BC7.stats
NanoStat --fastq BC8.fq.gz > BC8.stats
NanoStat --fastq BC9.fq.gz > BC9.stats

mkdir /t1-data/project/pregcare/sbush/MinION6/5.final_fq
cp /t1-data/project/pregcare/sbush/MinION6/4.demultiplexed_fq/barcode01/*.fastq.gz /t1-data/project/pregcare/sbush/MinION6/5.final_fq/BC1.fq.gz
cp /t1-data/project/pregcare/sbush/MinION6/4.demultiplexed_fq/barcode02/*.fastq.gz /t1-data/project/pregcare/sbush/MinION6/5.final_fq/BC2.fq.gz
cp /t1-data/project/pregcare/sbush/MinION6/4.demultiplexed_fq/barcode03/*.fastq.gz /t1-data/project/pregcare/sbush/MinION6/5.final_fq/BC3.fq.gz
cp /t1-data/project/pregcare/sbush/MinION6/4.demultiplexed_fq/barcode05/*.fastq.gz /t1-data/project/pregcare/sbush/MinION6/5.final_fq/BC5.fq.gz
cp /t1-data/project/pregcare/sbush/MinION6/4.demultiplexed_fq/barcode06/*.fastq.gz /t1-data/project/pregcare/sbush/MinION6/5.final_fq/BC6.fq.gz
cp /t1-data/project/pregcare/sbush/MinION6/4.demultiplexed_fq/barcode07/*.fastq.gz /t1-data/project/pregcare/sbush/MinION6/5.final_fq/BC7.fq.gz
NanoStat --fastq BC1.fq.gz > BC1.stats
NanoStat --fastq BC2.fq.gz > BC2.stats
NanoStat --fastq BC3.fq.gz > BC3.stats
NanoStat --fastq BC5.fq.gz > BC5.stats
NanoStat --fastq BC6.fq.gz > BC6.stats
NanoStat --fastq BC7.fq.gz > BC7.stats

mkdir /t1-data/project/pregcare/sbush/MinION7/5.final_fq
cp /t1-data/project/pregcare/sbush/MinION7/4.demultiplexed_fq/barcode01/*.fastq.gz /t1-data/project/pregcare/sbush/MinION7/5.final_fq/BC1.fq.gz
cp /t1-data/project/pregcare/sbush/MinION7/4.demultiplexed_fq/barcode02/*.fastq.gz /t1-data/project/pregcare/sbush/MinION7/5.final_fq/BC2.fq.gz
cp /t1-data/project/pregcare/sbush/MinION7/4.demultiplexed_fq/barcode03/*.fastq.gz /t1-data/project/pregcare/sbush/MinION7/5.final_fq/BC3.fq.gz
NanoStat --fastq BC1.fq.gz > BC1.stats
NanoStat --fastq BC2.fq.gz > BC2.stats
NanoStat --fastq BC3.fq.gz > BC3.stats

mkdir /t1-data/project/pregcare/sbush/MinION8/5.final_fq
cp /t1-data/project/pregcare/sbush/MinION8/4.demultiplexed_fq/barcode10/*.fastq.gz /t1-data/project/pregcare/sbush/MinION8/5.final_fq/BC10.fq.gz
cp /t1-data/project/pregcare/sbush/MinION8/4.demultiplexed_fq/barcode11/*.fastq.gz /t1-data/project/pregcare/sbush/MinION8/5.final_fq/BC11.fq.gz
cp /t1-data/project/pregcare/sbush/MinION8/4.demultiplexed_fq/barcode12/*.fastq.gz /t1-data/project/pregcare/sbush/MinION8/5.final_fq/BC12.fq.gz
NanoStat --fastq BC10.fq.gz > BC10.stats
NanoStat --fastq BC11.fq.gz > BC11.stats
NanoStat --fastq BC12.fq.gz > BC12.stats
