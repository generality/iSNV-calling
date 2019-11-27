iSNV-calling
===
For ebolavirus genome amplicon-sequencing protocol
---

Bbolavirus (EBOV) specific primers were designed to cover EBOV genome with overlaps, and the amplication products from one patient (or one blood sample) were pooled to prepare library and be sequenced by MiSeq. The scripts here are bioinformatics pipeline for the EBOV intrahost single nucleotide variation (iSNV) calling from the next-generation sequencing data.<br>  

The results of EBOV intrahost diversity generated from this pipeline were published in:
- Intra-host dynamics of Ebola virus during 2014, Nat Microbiol. 2016 Sep 5;1(11):16151. doi: 10.1038/nmicrobiol.2016.151.
The genome [1]

# Quality control of amplicon-sequencing data
We implemented a quality control and error correction strategy suggested by Schirmer et al based on their investigations on error patterns of amplicon sequencing data generated by Illumina’s MiSeq and Nextera XT Sample Preparation Kit, which was very similar with our experimental protocol [2] .Nucleotide-specific substitution errors could introduce false positive iSNVs, and these errors are more enriched at the start and end of reads (both read 1 and read 2). Therefore, the 10 bp at the start of all reads were trimmed, and we employed Sickle v1.3.3 for adaptively trimming of low quality bases at the end of reads with a threshold of Q20 and a minimum read length requirement of 100 bp. Next, Bayeshammer (included in SPAdes v3.5.0) was used for error correlation. Reads without their corresponding paired reads were disregarded. Then the properly paired reads were used as clean data for following bioinformatics analysis.<br>  

# Determination of EBOV iSNVs
The clean reads were aligned (pair-end alignment) to the earliest strain of EBOV 2014, H.sapiens-wt/GIN/2014/Makona-Kissidougou-C15 (Accession No. KJ660346.2), by using Bowtie2 v2.2.5  with default parameters. SAMtools v1.2 was employed to generate ‘mpileup’ files from SAM files produced by Bowtie2, with no limit of the maximum site depth. Home-made PERL scripts presented here were developed for iSNV calling using mpileup files as input. Briefly, first, at each EBOV reference genome site, the aligned low quality sequencing bases (< Q20) and insertions/deletions were excluded to reduce possible false positive, and the site depth as well as strand bias were re-calculated. Next, samples with no smaller than 3,000 sites which had a sequencing depth ≥ 100X or minor allele depth (disregarding the dominated allele) ≥ 5X were selected, and iSNVs were called from these valid sites of the selected samples. A series of criteria were used to detect a site with iSNVs:<br>  
(1)	The minor allele frequency ≥ 5%;<br>  
(2)	The depth of minor allele/alleles ≥ 5X;<br>  
(3)	The strand bias of minor allele/alleles was less than 10-fold;<br>  
(4)	The site did not locate in the EBOV-specific primer for the corresponding sample;<br>  
(5)	The site did not locate at the 1-100 bp downstream region of the first primer, and the 1-100 bp upstream of the last primer, for reducing potential bias induced by the transposon-based library preparation method.<br>  

# Reference
[1] Ni M, Chen C, Qian J, Xiao HX, Shi WF, Luo Y, Wang HY, Li Z, Wu J, Xu PS, Chen SH, Wong G, Bi YH, Xia ZP, Li W, Lu HJ, Ma JC, Tong YG, Zeng H, Wang SQ, Gao FG, Bo XC, Liu D. (2016). Intra-host dynamics of Ebola virus during 2014. Nature Microbiology 1(11), 16151.

[2]Schirmer, M., Ijaz, U.Z., D'Amore, R., Hall, N., Sloan, W.T., and Quince, C. (2015). Insight into biases and sequencing errors for amplicon sequencing with the Illumina MiSeq platform. Nucleic Acids Res 43, e37.
