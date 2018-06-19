# Assemble your PacBio data: Miniasm

This section is relevant to [group ID](work.md) 1-4.

## Note

In the next steps, we will describe the assembly itself for both ONT and PacBio.
Take note that the variable `$i` represents the fragment on which you are working. Thus, you need to:
 - replace `${i}` by your group number
 - define globally i to your group number (e.g. `i=3`)

# Overlapping

`bsub -q lr_course -n 5 -R "span[hosts=1]" "module add UHTS/Analysis/minimap2/2.8; minimap2 -x ava-pb -t 5 PB_chrom4_S${i}.reads.fasta PB_chrom4_S${i}.reads.fasta > PB_chrom4_S${i}.overlap.paf"`

# Assembly

`bsub -q lr_course "module add UHTS/Analysis/miniasm/0.2.r159.dirty; miniasm -f PB_chrom4_S${i}.reads.fasta PB_chrom4_S${i}.overlap.paf -R  > PB_chrom4_S${i}.miniasm.gfa"`

# Format conversion

The `gfa` format is very useful as it keeps the assembly path and how elements are connected. Unfortunately most down-stream tools dont accept it yet and we need to convert the `gfa` -format to `fasta`-format.

### Note

The next script are very lightweight and can quickly be executed without bsub!

`awk '/^S/{print ">"$2"\n"$3}' PB_chrom4_S${i}.miniasm.gfa | fold > PB_chrom4_S${i}.miniasm.fasta`

`software/assesV4.2.pl -t 1 -s 1 -i PB_chrom4_S${i}.miniasm.fasta -l -p > PB_chrom4_S${i}.miniasm.assess.txt`


# Done!

# Navigation

Go to Next section: [].

Return to [Main page](README.md)
