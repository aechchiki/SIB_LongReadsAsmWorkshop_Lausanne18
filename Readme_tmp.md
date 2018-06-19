# Genome Assembly: *Arabidosis thaliana*, chr4



# Assembly of chrom4 - almost finished

Members of each group can either work together or do each the excersise themself for the miniasm assembly.
For the canu assembly we will do only 1 assembly/group as it is more demanding CPU wise for the cluster.

**[Please find here the instructions](documentation/Assembly.md)**





# Steps post-assembly - still work to be done

The assembly of the genome is often only half of the rent and one should invest
a reasonable amount of time as well for the post-assembly steps. Probably even the same equivalent of time.

In particular as we often need these post-assembly steps in order to validate the
suggested assembly by the assembly algorithm.

So, one might want to re-do the assembly steps then again if major flaws/problems appear.
Even better (but sometimes limited by resources) is the strategy to do multiple attempts in parallel and evaluate them accordingly.

The steps which we will execute in order to analyze our assemblies  **[are documented in here](documentation/postAssembly.md)**


## Prediction of coding genes

The prediction/annotation of genes is an important part in order to understand better an assembly but can be extremely time intensive - depending whether a good model is already available.
Together with annotation of other elements it helps one to understand better the architecture of the assembled contigs, e.g.

 - repeat elements (Lines, Sines, viral elements)
 - low complexity regions
 - centromeric & telomeric repeats
 - potential contamination with other organism


 One of the best pipelines is the MAKER pipeline which is based on the Augustus gene predictor, in particular for eukaryotic genomes.
 Here we will simply predict genes on our contigs and use this information for later analysis.


## Evaluation of genome completeness

 - remapping of raw reads & coverage Analysis
 - BUSCO in order to get a feeling for completeness
 - mapping known resources (proteins, transcripts,....)
 - comparison to a reference

### Remapping of raw reads

As already done for the polishing step, one can map back the original reads to the assembled pieces and then visualize the
coverage. Ideally this should be rather uniform over the contig often with a decreasing slope towards the edges.

Repeats are often an exception to this rule. And depending on the way a tool is mapping the reads one might see an increase in coverage or actively suppress this.
E.g. blasr allows to distribute randomly reads mapping in multiple places leading to an even coverage profile. If one changes this behavior we will see reads mapping into
multiple places and subsequently an increased coverage. Both ways have their advantages and disadvantages.


### BUSCO

**B**enchmarking **U**niversal **S**ingle- **C**opy **O**rthologs [(webpage)](http://busco.ezlab.org/) can be run on predicted genes (proteins) or the entire genome.
In general it is a great exercise to run it on both and compare the results.
If the gene prediction is accurate, the results should resemble a lot.

Running it onto the entire genome is much more time consuming though.

## Comparison with a reference (if existing)

If a reference is available, one can study similarities & differences at many levels - exceeding the goal of this course.
We will only use within this scope a comparison using the `Mummer`suite
