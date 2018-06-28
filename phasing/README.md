# Phasing of given iSNVs
## The application in YFV amplicon-seq dataset

The phasing is implemented for selected serial iSNVs. As the distribution of serial iSNVs is limited, the phasing is performed by analyse the corresponding sequencing bases within reads.

The workflow of phasing analysis is as follows. 

* First, the sample IDs and serial iSNV positions of interested are manually selected and recorded in a file as `patientId_SerialGroupId_pos.txt`.
* Align the clean reads of YFV amplicon-seq fastq to reference YFV genomes, and obtain the SAM files.
* The script `phasing_check.pl` read the SAM files and serial iSNV information, and extract the corresponding bases at the iSNV positions from the reads, and output to a file with 'phasing' suffix. The output for YFV dataset are given in `sam.phasing.zip`.
* The script `related_phasing_pos.pl` load the files with a phasing suffix, and summarizes the status of phasing of each serial iSNV group. The output of YFV dataset is the `related_phasing_pos.out`. The ratio of phased reads, and corresponding  counts of phased and unphased reads are given.
