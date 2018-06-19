# Getting ready for assembly

Remember from [Previous section](documentation/ass_main.md) that, for each group:

 - all of you can go through the Miniasm assembly
 - only one of you may proceed to the Canu assembly


# Initial setup

Before we start with the assembly process we have to provide some kind of structure to our data.

1. Navigate to your working directory:  `cd scratch/cluster/monthly/`

2. Generate a folder with the name of your userID: `mkdir <username>`

3. Navigate into that new directory: `cd <username>`

4. According your [group ID](documentation/work.md), copy the long reads for your group into that folder from:  `/scratch/beegfs/monthly/eschmid/workshop_assembly/ATH/PACBIO` or `/scratch/beegfs/monthly/eschmid/workshop_assembly/ATH/ONT`

5. You'll need also to copy some working scripts: `cp -R /scratch/beegfs/monthly/eschmid/workshop_assembly/ATH/software .`

## Note

In the next step, we will describe the assembly itself for both ONT and PacBio.
Take note that the variable `$i` represents the fragment on which you are working. Thus, you need to:
 - replace `${i}` by your group number
 - define globally i to your group number (e.g. `i=3`)


## Working example

We will describe here two working examples to fetch data for PacBio or Nanopore to your current working directory (point 4 here above):

- Group 3 (PacBio, fragment 3): `cp /scratch/beegfs/monthly/eschmid/workshop_assembly/ATH/PACBIO/PB_chrom4_S3.reads.fasta .`

- Group 5 (ONT, fragment 1): `cp /scratch/beegfs/monthly/eschmid/workshop_assembly/ATH/ONT/ONT_chrom4_S1.reads.fasta .`


# Tutorials

According your [group ID](documentation/work.md), follow the tutorial:
 - [Miniasm: PacBio data](documentation/miniasm_pb.md)
 - [Miniasm: Nanopore data](documentation/miniasm_nano.md)
 - [Canu: PacBio data](documentation/canu_pb.md)
 - [Canu: Nanopore data](documentation/canu_nano.md)


 # Navigation

 Go to Next section: Choose any of the above (Tutorials).

 Return to [Main page](README.md)
