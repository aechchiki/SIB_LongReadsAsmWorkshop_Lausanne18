# Assembling a model plant genome

The oak genome is difficult to assemble in this session because of issues of biological (high level of heterozygosity >1%) and technical (PacBio coverage 20x) nature.

Instead, we will assemble another diploid plant: the model species *Arabidopsis thaliana*.
The raw data and the assembly are available from a very recent publication:

**High contiguity Arabidopsis thaliana genome assembly with a single nanopore flow cell** [(Michael et al., 2018, Nat. Comm)](https://www.nature.com/articles/s41467-018-03016-2#Sec9)

# Raw data

From this publication, we have the following data available from ENA:

 - [1 raw PacBio Sequel cell (`.bam`)](ftp://ftp.sra.ebi.ac.uk/vol1/ERA111/ERA1116568/bam/pb.bam)
 - [1 Nanopore run in (`.fq`)](ftp://ftp.sra.ebi.ac.uk/vol1/ERA111/ERA1116595/fastq/ont.fq.gz)


# Working data

We will try to assemble the chromosome 4 from *Arabidopsis* ( ~20 Mbp: the smallest chromosome).

We reported [here](documentation/preAssembly.md) the steps to subset reads from chromosome 4 from the whole raw dataset. If time allows, we will repeat the sub-setting of chromosomal reads together.

In short, the following steps have been executed:

 - raw data download and conversion
 - reads mapping onto chromosome 4
 - reads splitting for groups
 - reads filtering according to mapping
 - raw data subset to the reads of interest

 # Navigation

 Go to Next section: [Setting up working groups](documentation/work.md)

 Return to [Main page](README.md)
