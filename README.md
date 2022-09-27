![petanc](https://user-images.githubusercontent.com/47217665/192289530-60a0d40b-604c-44f7-9815-ceb7ddbe56c4.png)

**The pipeline PETANC for “Plasmid-Exploration Typing Assembly N’Contig-ordering” ... This pipeline clean reads and ...   for *Escherichia coli* strain** PE 

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
- [MLST](#MLST)
- [*fimH*](#fimH)
- [Phylogroup](#Phylogroup)
- [Genes of virulences](#Genes-of-virulences)
- [Genes of resistance](#Genes-of-resistance)
- [Plasmides](#Plasmides)
- [Classification of contigs to plasmid or chomosome](#Classification-of-contigs-to-plasmid-or-chomosome)
- [Annotation](#Annotation)
- [Capsule systems](#Capsule-systems)
- [*Shigella or EHEC*](#Shigella-or-EHEC)
- [Layout Excel](#Layout-Excel)

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

:hourglass: to use MultiQC again

## Serotype

We use abricate with the database ecoh (--db ecoh) and threshold for the percentage of identity (--minid 80) and coverage (--minicov 90)

[https://github.com/tseemann/abricate](https://github.com/tseemann/abricate)

:heavy_check_mark: abricate : 0.8.11

## MLST

We use mlst to looking for Escherichia coli Warwick MLST (--scheme ecoli) and Escherichia coli Pasteur MLST (--scheme ecoli_2) 

*"This publication made use of the PubMLST website (https://pubmlst.org/) developed by Keith Jolley [(Jolley & Maiden 2010, BMC Bioinformatics, 11:595)](https://doi.org/10.1186/1471-2105-11-595) and sited at the University of Oxford.  The development of that website was funded by the Wellcome Trust".*

[https://github.com/tseemann/mlst](https://github.com/tseemann/mlst)

:heavy_check_mark: mlst : 2.16.2

## *fimH*

We use FimTyper to know the allele of *fimH* with a threshold for %identity (-k 95.00) and a minimum length for the overlap (-l 0.60)

[Roer L, Tchesnokova V, Allesoe R, Muradova M, Chattopadhyay S, Ahrenfeldt J, Thomsen MCF, Lund O, Hansen F, Hammerum AM, Sokurenko E, and Hasman H. J  Development of a web tool for Escherichia coli subtyping based on fimh alleles. Clin Microbiol. 2017. 55(8): 2538-2543. https://doi.org/10.1128/JCM.00737-17 ](https://doi.org/10.1128/JCM.00737-17)

[https://bitbucket.org/genomicepidemiology/fimtyper/src/master/](https://bitbucket.org/genomicepidemiology/fimtyper/src/master/)

:heavy_check_mark: FimTyper : 1.1

## Phylogroup

We use the ClermonTyping to know the Phylogroup of *Escherichia coli* with all the contigs (--threshold 0)

[Beghain J, Bridier-Nahmias A, Le Nagard H, E.Denamur and O.Clermont. Clermontyping: an easy-to-use and accurate in silico method for Escherichia genus strain phylotyping. Microb Genom. 2018 Jun 19. PMID: 29916797 PMCID: PMC6113867 https://doi.org/10.1099/mgen.0.000192 ](https://doi.org/10.1099/mgen.0.000192)

[Clermont O, Dixit OVA, Vangchhia B, Condamine B, Dion S, Bridier-Nahmias A, Denamur E, Gordon D. Characterization and rapid identification of phylogroup G in Escherichia coli, a lineage with high virulence and antibiotic resistance potential. Environ Microbiol. 2019 Jun 12. PMID: 31188527 https://doi.org/10.1111/1462-2920.14713](https://doi.org/10.1111/1462-2920.14713)

:heavy_check_mark: ClermonTyping : 21.03

## Genes of virulences

With Abricate, we are looking for genes of virulences with thresholds for the percentage of identity (--minid 80) and coverage (--minicov 90). The database was done by Guilhem Royer in Kieffer *et al*, 2019. 

[Kieffer N, Royer G, Decousser JW, Bourrel AS, Palmieri M, Ortiz De La Rosa JM, Jacquier H, Denamur E, Nordmann P, Poirel L. mcr-9, an Inducible Gene Encoding an Acquired Phosphoethanolamine Transferase in Escherichia coli, and Its Origin. Antimicrob Agents Chemother. 2019 Aug 23;63(9):e00965-19. doi: 10.1128/AAC.00965-19. Erratum in: Antimicrob Agents Chemother. 2019 Oct 22;63(11): PMID: 31209009; PMCID: PMC6709461. https://doi.org/10.1128/AAC.00965-19](https://doi.org/10.1128/AAC.00965-19)

[https://github.com/tseemann/abricate](https://github.com/tseemann/abricate)

:heavy_check_mark: abricate : 0.8.11

## Genes of resistance


[https://github.com/tseemann/abricate](https://github.com/tseemann/abricate)

:heavy_check_mark: abricate : 0.8.11

:hourglass: to change for resfinder

## Plasmid
## Classification of contigs to plasmid or chomosome

## Annotation
### Clean gbk
## Capsule systems
https://research.pasteur.fr/en/tool/capsulefinder/
## *Shigella* or EHEC
## Layout excel
