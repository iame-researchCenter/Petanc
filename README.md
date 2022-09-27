![petanc](https://user-images.githubusercontent.com/47217665/192289530-60a0d40b-604c-44f7-9815-ceb7ddbe56c4.png)

**The pipeline PETANC for “Plasmid-Exploration Typing Assembly N’Contig-ordering” is a pipeline for Illumina paired-end reads of *Escherichia coli* strain.**

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

[Philip Ewels, Måns Magnusson, Sverker Lundin, Max Käller, **MultiQC: summarize analysis results for multiple tools and samples in a single report**, Bioinformatics, Volume 32, Issue 19, 1 October 2016, Pages 3047–3048, https://doi.org/10.1093/bioinformatics/btw354 ](https://doi.org/10.1093/bioinformatics/btw354)

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

[Alexey Gurevich, Vladislav Saveliev, Nikolay Vyahhi, Glenn Tesler, **QUAST: quality assessment tool for genome assemblies**, Bioinformatics, Volume 29, Issue 8, 15 April 2013, Pages 1072–1075, https://doi.org/10.1093/bioinformatics/btt086](https://doi.org/10.1093/bioinformatics/btt086)


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

[Roer L, Tchesnokova V, Allesoe R, Muradova M, Chattopadhyay S, Ahrenfeldt J, Thomsen MCF, Lund O, Hansen F, Hammerum AM, Sokurenko E, and Hasman H. J  **Development of a web tool for Escherichia coli subtyping based on fimh alleles**, Clin Microbiol. 2017. 55(8): 2538-2543. https://doi.org/10.1128/JCM.00737-17 ](https://doi.org/10.1128/JCM.00737-17)

[https://bitbucket.org/genomicepidemiology/fimtyper/src/master/](https://bitbucket.org/genomicepidemiology/fimtyper/src/master/)

:heavy_check_mark: FimTyper : 1.1

## Phylogroup

We use the ClermonTyping to know the Phylogroup of *Escherichia coli* with all the contigs (--threshold 0)

[Beghain J, Bridier-Nahmias A, Le Nagard H, E.Denamur and O.Clermont. **Clermontyping: an easy-to-use and accurate in silico method for Escherichia genus strain phylotyping**. Microb Genom. 2018 Jun 19. PMID: 29916797 PMCID: PMC6113867 https://doi.org/10.1099/mgen.0.000192 ](https://doi.org/10.1099/mgen.0.000192)

[Clermont O, Dixit OVA, Vangchhia B, Condamine B, Dion S, Bridier-Nahmias A, Denamur E, Gordon D. **Characterization and rapid identification of phylogroup G in Escherichia coli, a lineage with high virulence and antibiotic resistance potential**. Environ Microbiol. 2019 Jun 12. PMID: 31188527 https://doi.org/10.1111/1462-2920.14713](https://doi.org/10.1111/1462-2920.14713)

:heavy_check_mark: ClermonTyping : 21.03

## Genes of virulences

With Abricate, we are looking for genes of virulences with thresholds for the percentage of identity (--minid 80) and coverage (--minicov 90). The database was done by Guilhem Royer in Kieffer *et al*, 2019. 

[Kieffer N, Royer G, Decousser JW, Bourrel AS, Palmieri M, Ortiz De La Rosa JM, Jacquier H, Denamur E, Nordmann P, Poirel L. **mcr-9, an Inducible Gene Encoding an Acquired Phosphoethanolamine Transferase in Escherichia coli, and Its Origin**. Antimicrob Agents Chemother. 2019 Aug 23;63(9):e00965-19. doi: 10.1128/AAC.00965-19. Erratum in: Antimicrob Agents Chemother. 2019 Oct 22;63(11): PMID: 31209009; PMCID: PMC6709461. https://doi.org/10.1128/AAC.00965-19](https://doi.org/10.1128/AAC.00965-19)

[https://github.com/tseemann/abricate](https://github.com/tseemann/abricate)

:heavy_check_mark: abricate : 0.8.11

## Genes of resistance

Currently, we use abricate with thresholds for the percentage of identity (--minid 80) and coverage (--minicov 90) and the database of ResFinder.

[Florensa AF, Kaas RS, Clausen PTLC, Aytan-Aktug D, Aarestrup FM. **ResFinder - an open online resource for identification of antimicrobial resistance genes in next-generation sequencing data and prediction of phenotypes from genotypes**. Microb Genom. 2022 Jan;8(1):000748. doi: 10.1099/mgen.0.000748. PMID: 35072601; PMCID: PMC8914360.https://doi.org/10.1099/mgen.0.000748](https://doi.org/10.1099/mgen.0.000748)

[https://github.com/tseemann/abricate](https://github.com/tseemann/abricate)

:heavy_check_mark: abricate : 0.8.11

:hourglass: To change for resfinder

## Plasmid

With Abricate, we are looking for genes of plasmids with thresholds for the percentage of identity (--minid 80) and coverage (--minicov 90). The database is PlasmidFinder_DB. 

[Carattoli A, Zankari E, García-Fernández A, Voldby Larsen M, Lund O, Villa L, Møller Aarestrup F, Hasman H. **In silico detection and typing of plasmids using PlasmidFinder and plasmid multilocus sequence typing**. Antimicrob Agents Chemother. 2014 Jul;58(7):3895-903. doi: 10.1128/AAC.02412-14. Epub 2014 Apr 28. PMID: 24777092; PMCID: PMC4068535. https://doi.org/10.1128/AAC.02412-14](https://doi.org/10.1128/AAC.02412-14)

[https://github.com/tseemann/abricate](https://github.com/tseemann/abricate)

:heavy_check_mark: abricate : 0.8.11

## Classification of contigs to plasmid or chomosome

We classify contigs from the assembly according to their location (i.e. plasmid or chromosome) with PlaScope. 

[Royer G, Decousser JW, Branger C, Dubois M, Médigue C, Denamur E, Vallenet D. **PlaScope: a targeted approach to assess the plasmidome from genome assemblies at the species level**. Microb Genom. 2018 Sep;4(9):e000211. doi: 10.1099/mgen.0.000211. PMID: 30265232; PMCID: PMC6202455.https://doi.org/10.1099/mgen.0.000211](https://doi.org/10.1099/mgen.0.000211)

[https://github.com/labgem/PlaScope](https://github.com/labgem/PlaScope)

:heavy_check_mark: PlaScope : 2018

## Annotation

We annote the bacteria assembly (--gcode 11) with prokka.

[Seemann T. **Prokka: rapid prokaryotic genome annotation**. Bioinformatics, 2014. PMID: 24642063 https://doi.org/10.1093/bioinformatics/btu153](https://doi.org/10.1093/bioinformatics/btu153)

[https://github.com/tseemann/prokka](https://github.com/tseemann/prokka)

:heavy_check_mark: Prokka : 1.14

Sometime, we have a problem with the first line of each contig of the genebank file from prokka. The locus name is crushed by length of the sequence, so the script `cleanGBK_locusProkka.py` rename the locus without the coverage of spades output.

## Capsule systems

We detects capsule systems with CapsuleFinder for the model "ABC", "GroupIV_e_stricte", "GroupIV_f", "GroupIV_s_stricte", "PGA", "Syn_cps3", "Syn_has" and "Wzy_stricte".

[Rendueles O, Garcia-Garcerà M, Néron B, Touchon M, Rocha EPC. **Abundance and co-occurrence of extracellular capsules increase environmental breadth: Implications for the emergence of pathogens**. PLoS Pathog. 2017 Jul 24;13(7):e1006525. doi: 10.1371/journal.ppat.1006525.https://doi.org/10.1371/journal.ppat.1006525](https://doi.org/10.1371/journal.ppat.1006525)

https://research.pasteur.fr/en/tool/capsulefinder/

:heavy_check_mark: CapsuleFinder : 02/02/2018

## *Shigella* or EHEC

With Abricate, we are looking for the gene *ipaH3* of *Shigella flexneri* Y strain PE57 (CP042980.1 from 1400701 to 1402416) with thresholds for the percentage of identity (--minid 95) and coverage (--minicov 95).

[https://github.com/tseemann/abricate](https://github.com/tseemann/abricate)

:heavy_check_mark: abricate : 0.8.11

## Layout excel

All the data are synthesized by the script `petanc_layout.py` to make an Excel file. 
