# Genome assembly and evaluation using long reads (at the example of the tree Quercus robur)

Quercus robur as a theme ?

Issues:

 - "large" genome. not possible to assemble
 - not enough coverage
 - tools such as BUSCO wont be meaningful on partial genome (too slow on entire genome)

Possible solutions:

 - assemble only 1 chromosome
 - ask Christoph for PacBio reads for 1 chromosome - and reference


New problem: As C. Plomion did not reply and does not participate at the course we might not have enough PacBio reads to actually assemble anything.

--> possible solution might be to take e.g. the longest scaffold/contig and the mapped PacBio read + simulate PacBio reads of the contig with existing programs to obtain enough coverage.




Advantage:

 - both C. Plomion and S. Salzberg assembled the genome and can spike in some experience
 - C. Plomion has chrom. level assembly and could present the genome, their strategy their way of evaluation
 - S. Salzberg has multiple trees and some crazy genome sizes in trees to present (MaSurCa as well)
 - R. Waterhouse has new plant categories for BUSCO now - which have not yet been tested on the oak

## Date

Week from 18-22 of June, ideally 2 days.

8.30  - 8.45 Welcoming
8.45  - 9.45 session 1
9.45  - 10.00 Coffee break
10.00 - 11.00 session 1 / second part
11.00 - 12.00 session 2
12.00 - 13.00 lunch break
13.00 - 14.00 session 2 / second part
14.00 - 15.00 session 3
15.00 - 15.15 coffee break
15.15 - 16.15 session 3 /second part
16.15 - 17.00 to determine, maybe first introduction to UNIX again

## 1st day: assembly of long reads

08.30 - 08.45 Welcoming
08.45 - 09.15 Philippe Reymond: Assembly of the Oak genome and assessment of somatic mutations
09.15 - 09.45 Emanuel Schmid: overview on single molecule sequencing and assembly of long reads
09.45 - 010.00 Coffee break
10.00 - 10.30 Sergey Koren: presentation of the assembler Canu and latest developments
10.30 - 11.30 demo of Sergey Koren, assembly of E. coli
11.30 - 12.30 lunch break
12.30 - 15.00 hands-on: assembling in small groups a contig from Quercus robur
15.00 - 15.15 coffee break
15.15 - 17.00 hands-on: continuing genome assembly

## 2nd day: evaluation of genome assemblies

08.30 - 08.45 Welcoming
08.45 - 09.15 Robert Waterhouse: BUSCO for evaluation of assemblies and more
09.15 - 09.45 Damien Liebherr: SwissProt, curation of (plant) proteins
09.45 - 010.00 Coffee break
10.00 - 10.30 Emanuel Schmid: comparison and evaluation of genome assemblies
10.30 - 12.00 hands-on: comparison of assembled genomes (& to reference)
12.00 - 13.00 lunch break
13.00 - 15.00 hands-on: comparison of assembled genomes
15.00 - 15.15 coffee break
15.15 - 17.00 hands-on: gene prediction



## 2nd day: workshop, hands-on

 - Assembly of a small organism (E. coli?) - or 1 Oak chromosome -  using minimap2 + miniasm.
 - maybe provide each group with a different coverage ?
 - Polishing will be skipped and polished version provided
 - Comparison of assembled/polished versions with reference (mummer)
 - evaluation of assemblies
   - 1 group using BUSCO
   - 1 group remapping reads
   - .....
 - prediction of genes --> proteins --> effect of polishing

We would need somebody to assist with the workshop (maybe 2, depending on the size and interest).


## Adresses

[1]Steven L. Salzberg, PhD
	Bloomberg Distinguished Professor
	Professor of Biomedical Engineering, Computer Science, and Biostatistics
	Director, Center for Computational Biology
	McKusick-Nathans Institute of Genetic Medicine, Johns Hopkins University School of Medicine
	Office: Welch Medical Library 107 | Malone 325
	Lab: Salzberg Lab
	410-614-6112
	salzberg@jhu.edu

[2] Christophe Plomion
	Plant Molecular Geneticist
	« Biodiversity Genes & Communities, BIOGECO» research unit, INRA
	69, route d'Arcachon, 33610 CESTAS– France
	Phone : +33 (0)5 57 12 27 65
	Email : plomion.christophe@inra.fr

[3] Robert Waterhouse
	Route de Chavannes 13
	1007 Lausanne
	Switzerland

[4] Adam M. Phillippy, Ph.D.
	 Investigator
	Computational and Statistical Genomics Branch
	Head
	Genome Informatics Section
	B.S. Loyola University Maryland, 2002
	M.S. University of Maryland, College Park, 2009
	Ph.D. University of Maryland, College Park, 2010
	T: (301) 451-8748
	F: (301) 402-2170
	E: adam.phillippy@nih.gov
	Building 49, Room 4A22
	49 CONVENT DR
	BETHESDA, MD 20892
