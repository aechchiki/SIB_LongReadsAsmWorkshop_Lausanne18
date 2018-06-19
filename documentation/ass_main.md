# Assembly

## Selected software

We will use two approaches for genome assembly from long reads:

- Minimap2 + Miniasm
- Canu

### Minimap2 + Miniasm

The [Minimap2](https://github.com/lh3/minimap2)+ [Miniasm](https://github.com/lh3/miniasm) approach is extremely rapid and allows us to obtain within short time the assembly of chromosome 4. It manages such an advantage in speed by assembling directly the raw reads without an hierarchical approach & correcting the reads.

There are limitations though:

- large repeats and repeats with very high similarity are much more difficult to resolve
- the assembly has the original long read error rate
- often many non-assembled fragments remain in the assembly

It is separated into 2 steps:

 1. overlapping of all reads with each other
 2. overlap-graph and steps to simplify it further (bubble popping, removing dead ends)

Suggested read:

- [Minimap + Miniasm](https://arxiv.org/pdf/1512.01801.pdf)
- [Minimap2 (major update)](https://arxiv.org/pdf/1708.01492v5.pdf)

### Canu

You had a pretty good introduction to [Canu](https://canu.readthedocs.io/en/latest/) by Sergey Koren himself this morning :)

Suggested read:

- [Canu](https://genome.cshlp.org/content/27/5/722)

# Navigation

Go to Next section: [Set up your assembly directory](assembling.md)

Return to [Main page](README.md)
