![petanc](https://user-images.githubusercontent.com/47217665/192289530-60a0d40b-604c-44f7-9815-ceb7ddbe56c4.png)

**BLABLA : This pipeline clean reads and ...   for *Escherichia coli* strain** PE 

# Dependencies Installation
Petanc is a pipeline used on the cdc of IAME. Our cdc use Docker. Each software have been installed first in container.
In the file `parameter.py` , there is a list of the names, versions and articles to cite for a software in each images.

# Command line options
```
% petanc.py
Script usage :
	-h          : print this message and exit
	--fasta     : start the analyse from fasta sequences
	--directory : directory with fastq or fasta files
```

# Differents steps
- [Create a list of samples and clean names of fastq files](#Create-a-list-of-samples-and-clean-names-of-fastq-files)
- [Quality Control of the raw reads](#Quality-Control-of-the-raw-reads)
- [Trimming](#Trimming)
- [Quality Control of the clean reads](#Quality-Control-of-the-clean-reads)
- [Assembly](#Assembly)
- [Quality control of the assembly](#Quality-control-of-the-assembly)
- [Serotype](#Serotype)


## Create a list of samples and clean names of fastq files
A list of the samples is created from the directory (options). If the samples are fastq files, we clean their names. Depending on the sequencing technology, the reads files can be shortened. For exemple, reads files from MiSeq have a name Name_S[0-9]\*_L001_R1_001.fastq.gz and could be shortened to Name_R1.fastq.gz . This pipeline clean names for HiSeq, MiSeq and NextSeq.

## Quality Control of the raw reads
We used fastQC to generated a quality control of each file (with default parameters) of reads and MultiQC for visualize all fastQC output in one glance.

[Philip Ewels, Måns Magnusson, Sverker Lundin, Max Käller, MultiQC: summarize analysis results for multiple tools and samples in a single report, Bioinformatics, Volume 32, Issue 19, 1 October 2016, Pages 3047–3048, https://doi.org/10.1093/bioinformatics/btw354 ](https://doi.org/10.1093/bioinformatics/btw354)

:heavy_check_mark: FastQC 
:heavy_check_mark: QUAST: 5.0.2

## Trimming
Trimgalore is used to trim and filer the paired-end reads (--paired). We trim low-quality ends (-q 30) from reads in addition to adapter removal. We discard reads that became shorter (-t 50). We keep reads unpaired (--retain_umpaired).

## Quality control of the clean reads
The quality control after trimming is exactly the same before trimmming.

## Assembly
The assembly is done with Spades

## Quality control of the assembly
QUAST is used to check the quality of the assembly. 

[Alexey Gurevich, Vladislav Saveliev, Nikolay Vyahhi, Glenn Tesler, QUAST: quality assessment tool for genome assemblies, Bioinformatics, Volume 29, Issue 8, 15 April 2013, Pages 1072–1075, https://doi.org/10.1093/bioinformatics/btt086](https://doi.org/10.1093/bioinformatics/btt086)


:heavy_check_mark: QUAST: 5.0.2

:hourglass: Reflexion to use MultiQC again

## Serotype

We use abricate with the database ecoh (--db ecoh) and threshold for the percentage of identity (--minid 80) and coverage (--minicov 90)

[https://github.com/tseemann/abricate](https://github.com/tseemann/abricate)

:heavy_check_mark: abricate : 0.8.11

## MLST


:heavy_check_mark: abricate : 0.8.11

## fimH
## Clermontyping
## Virulence
## Resistance
## Plasmides
## Plascopes
## Annotation
### Clean gbk
## Capsule
## Shigella or EHEC
## Layout excel
