# Organizing set of Raw Data

# our Qr3.5

After splitting into contigs, these are the 3 longest ones:

90997|>Contig889|4|90997
93565|>Contig86|7|93565
99588|>Contig2538|3|99588

Really, shabby :(

Try re-assembly with PacBio only:

    time minimap2 -x ava-pb -t 20 allCells.fasta allCells.fasta > allCells.overlap.paf
module add UHTS/Analysis/miniasm/0.2.r159.dirty;
time miniasm -f allCells.fasta allCells.fasta  allCells.overlap.paf -s 1000 -R > allCells.miniasm.gfa

Real time: 1138.095 sec; CPU: 13050.582 sec

real    19m1.481s
user    204m20.172s
sys     13m13.038s


Real time: 298.348 sec; CPU: 171.876 sec

real    4m58.358s
user    2m21.652s
sys     0m30.225s

Incredibly fast....



# Plomion
I downloaded the assembled genome from the French consortium of Quercus robur [http://www.oakgenome.fr/?page_id=587]
and extracted the first chromosome.
This chromosome has the size of 55'068'941 bp = 55Mbp - maybe too large ?


It has as well some gaps (#/bps/%):                               1488    1'597'513    2.90

Lets have a look at their asembly:

/home/eschmid/tools/perl_script_contigs/assesV4.2.pl -t 1 -s 720 -i quercus_robur.Plomion.fasta -l
quercus_robur.Plomion

    INFORMATION:
    A threshold of 1 bp was applied for ALL measures and genome size estimated 720000000 Bp
    Gaps are considered from a size of  10 bp on, below ignored
    lower case bases are often either low confidence bases (PB) or filled gaps (PBHoney,SSPACE)



    overall stats

    GC ratio (in %):                            35.65
    Shortest sequence length:                    2095
    Longest sequence length:                115639695
    bp-sum of all sequences :               814368469
    Number of sequences >1:                       550
    N50 sequence length :                    55068941
    NG50 sequence length :                   57352617
    lower Case (#/bps/%):                           0          0    0.00
   gaps (#/bps/%):                              22088   24054386    2.95

----------------------------------------------------------------------------------------------------


    overall sequence size distribution (# of contigs & scaffolds)

                              absolute    |      % total sequ.
    100-500:                       0                 0.00
    500-1000:                      0                 0.00
    1000-10000:                   96                17.45
    10000-100000:                248                45.09
    100000-1000000:              175                31.82
    >1000000:                     31                 5.64



    overall sequence size distribution (bp SUM of contigs & scaffolds)

                              absolute Bps |        % total Bps |   % genome size
    100-500:                       0                   0.00              0.00
    500-1000:                      0                   0.00              0.00
    1000-10000:               595845                   0.07              0.08
    10000-100000:            9491705                   1.17              1.32
    100000-1000000:         53001592                   6.51              7.36
    >1000000:              751279327                  92.25            104.34

Crazy what you can do with enough money...

## genome comparison

A first comparison with mashmap, as super fast:

![](QV35_vsPlomion.mashmap.png)

I think it looks better than expected but there are quite a few blue lines on the red diagonal.
This is less great meaning that maybe a few thing were wrongly combined (scaffolded).


With my blast-based approach

| comparison | cutoff                                               | prov..genome.size | genome.covered | genome.missing | cov.1.   | cov.2.  | cov.above. | median.pident | median.indels | median.mismatch | calc..genome.size. |            |
|------------|------------------------------------------------------|-------------------|----------------|----------------|----------|---------|------------|---------------|---------------|-----------------|--------------------|------------|
| blastComp/quercus_robur.V3.5.soft_masked.blastGenome | 0                 | 814368469      | 42.9905        | 57.0095  | 43.0103 | 2.2103     | 0.7214        | 98.519        | 0.00412         | 0.00943            | 813.640812 |
| blastComp/quercus_robur.V3.5.soft_masked.blastGenome | 500               | 814368469      | 41.1819        | 58.8181  | 41.7151 | 1.8618     | 0.5566        | 98.148        | 0.00945         | 0.0101             | 813.640812 |
| blastComp/quercus_robur.V3.5.soft_masked.blastGenome | 1000              | 814368469      | 35.5157        | 64.4843  | 36.8651 | 1.1657     | 0.4365        | 98.156        | 0.01162         | 0.00986            | 813.640812 |
| blastComp/quercus_robur.V3.5.soft_masked.blastGenome | 2500              | 814368469      | 19.654         | 80.346   | 22.3476 | 0.2482     | 0.0098        | 98.51         | 0.00943         | 0.00787            | 813.640812 |
| blastComp/quercus_robur.V3.5.soft_masked.blastGenome | 5000              | 814368469      | 6.831          | 93.169   | 9.7394  | 0.0398     | 0.0034        | 98.733        | 0.00818         | 0.0067             | 813.640812 |
| blastComp/quercus_robur.V3.5.soft_masked.blastGenome | 7500              | 814368469      | 1.4965         | 98.5035  | 4.4362  | 0.0097     | 0.0022        | 98.8315       | 0.00751         | 0.00592            | 813.640812 |
| blastComp/quercus_robur.V3.5.soft_masked.blastGenome | 10000             | 814368469      | -0.6958        | 100.6958 | 2.2534  | 0.0024     | 0             | 98.9035       | 0.00721         | 0.00563            | 813.640812 |

Now this is very strange. Essentially we cover only less than half their genome size with matches from our genome. As theirs is 813Mbp (-some gaps) it would mean that only 400 Mbp of our genome map. We have 720Mbp assembled with still quite some gaps.
But there is still a confusing gap. Maybe blast matches are too quickly abonned due to large differences and subsequently often <500bp long ?

The similarity is though remarkably high (here at 2.5kbp)

![](similarity2kbp.png)

## mummer

```bash
[REF]                [QRY]
[Sequences]
TotalSeqs                        550                85557
AlignedSeqs              498(90.55%)        84092(98.29%)
UnalignedSeqs              52(9.45%)          1465(1.71%)

[Bases]
TotalBases                 814368469            719779348
AlignedBases       579416100(71.15%)    600979405(83.49%)
UnalignedBases     234952369(28.85%)    118799943(16.51%)

[Alignments]
1-to-1                        380646               380646
TotalLength                461170810            460171622
AvgLength                    1211.55              1208.92
AvgIdentity                    96.25                96.25

M-to-M                       1040321              1040321
TotalLength                828428492            826994031
AvgLength                     796.32               794.94
AvgIdentity                    94.95                94.95

[Feature Estimates]
Breakpoints                  2080230              1961349
Relocations                    83188                79749
Translocations                193632                83520
Inversions                      4028                 9774

Insertions                    588237               799238
InsertionSum               416876384            305700264
InsertionAvg                  708.69               382.49

TandemIns                       1087                 1232
TandemInsSum                  176401               294451
TandemInsAvg                  162.28               239.00
```

This is now speaking a different language. Essentially 90% of their scaffolds aligne against ours and >98% of ours against theirs. If we break it down to nucleotides, they have >71% aligning against our and we have >83% aligning against theirs.
Considering that we have still 70Mbp gaps in our assembly and they have 24Mbp as well it becomes better understandable.

 Which is roughly 600Mbp common ground - not too bad, could be better though. Indeed, the length of aligned pieces is pretty short with 1,2kbp in average for single hits and 800bp for multiple ones. This will obviously be much worse for blast and explains the bad result in blast.

In general, the assemblies are very different - which is less great but we had unfortunately not the same ressources and techniques....

# Selecting a canidate for re-assembly

Scaffold  Qrob_H2.3_Sc0000306 looks good with 914kbp - that should go rather quickly.
2.63% being gaps. I did the entire excersise but it was not a good example - fragmenting too much.
If we look at the mapping on the scaffold it becomes obvious why:

![](Qrob_H2.3_Sc0000306_coverage.png)

Looking at coverage graph (log) we can see that actually many points exist with 0 coverage.
Likely they are then gaps where nothing is mapping - and indeed nothing from our RAW reads is mapping.

    awk '{if($4=="0") print }' readsMappedOnto.refQrob_H2.3_Sc0000306.coverage.bedGraph | wc -l
    87

Indeed 87 spots have no coverage at all (--> will be broken likely). If we look at the size:

    awk '{if($4=="0") print $3-$2 }' readsMappedOnto.refQrob_H2.3_Sc0000306.coverage.bedGraph  | sort | uniq -c
     56 1
      1 1064
      1 1196
      1 15
      1 153
      1 1786
      9 2
      1 2546
      1 2579
      1 280
      1 294
      4 3
      1 3111
      2 4
      2 7
      1 721
      1 753
      1 8
      1 821

There are some really big ones....

## systematic approach

The contigs are often heavily fragmented. To get the largest possible continious contig, I wrote a script
to split scaffolds into contigs:

    perl  /home/eschmid/tools/perl_script_contigs/fasta_scaffold2contig.pl -i quercus_robur.Plomion.fasta -n 50 > quercus_robur.Plomion.Contigs50.fasta

Taking the 40 largest
    grep "^>" quercus_robur.Plomion.Contigs50.fasta | awk '{FS="|"}{OFS="|"}{print $3 $1,$2,$3}' | sort -n | tail -40 | perl -npe 's/[0-9]+\>//' > top40.reads

    /home/eschmid/tools/perl_script_contigs/extract_scaffold_version4.pl -i top40.reads -f quercus_robur.Plomion.Contigs50.fasta -s > top40.fasta
  query: MyContigs.txt
	        FASTA: quercus_robur.Plomion.Contigs.fasta
		-d: 0
		-s: 1
		-c: 0
		-r: 0
		-v: 0


Then I explode it and rename the headers:

    /home/eschmid/tools/perl_script_contigs/correct_fasta_format.pl -i top40.fasta -r -p Contig -x
    renumbering chosen
    each FASTA will give rise to a new file
    manually prefix Contig chosen
    This should be first FASTA >Qrob_Chr08|577|621829


Then I took all our PacBio reads (fasta format) and aligned them with ``minimap2`` against it.

    module add UHTS/Analysis/minimap2/2.8;
    for contig in `ls Contig*.fasta`; do myContig=$(basename $contig .fasta); minimap2 $contig allCells.fasta -x map-pb > $myContig.paf ; done

The obvious candidates of interest are the Contig27 and 31 with 621829 and 535545 respectively.


Lets extract now for each of these 2 contigs the reads which map in order to see whether we have a "decent "   coverage.

    awk '{print $1}' Contig27.paf | sort | uniq > Contig27.reads &
    awk '{print $1}' Contig31.paf | sort | uniq > Contig31.reads &

    /home/eschmid/tools/perl_script_contigs/extract_scaffold_version4.pl -i Contig27.reads -f allCells.fasta -s > Contig27.reads.fasta &
    /home/eschmid/tools/perl_script_contigs/extract_scaffold_version4.pl -i Contig31.reads -f allCells.fasta -s > Contig31.reads.fasta &


Contig27.reads.fasta


    time minimap2 -x ava-pb -t 20 Contig27.reads.fasta Contig27.reads.fasta > Contig27.reads.paf &
    time minimap2 -x ava-pb -t 20 Contig31.reads.fasta Contig27.reads.fasta > Contig31.reads.paf &

    time miniasm -f simulatedAndRealReads.Contig31.fasta simulatedAndRealReads.Contig31.paf -s 1000 -R > simulatedAndRealReads.Contig31.miniasm.gfa &


This is actually pretty different:

 70290 Contig27.reads
 37333 Contig31.reads


So, on the contig27 we have actually almost the double of reads mapping which is surprising as the length is not that different. The average length of my reads was ~3000 and 200 reads would theoretically =~1x coverage. 10000x reads then would be already sufficient. This indicates already that it will be pretty tricky as we have very likely strange repeats here.


Lets map again onto reference all reads:

    module add UHTS/Analysis/samtools/1.4
    minimap2 Contig27.fasta Contig27.reads.fasta -a -t 30 | samtools view -bhS | samtools sort -o validateContig27.bam &

    minimap2 Contig31.fasta Contig31.reads.fasta -a -t 30 | samtools view -bhS | samtools sort -o validateContig31.bam &

    samtools index validateContig27.bam &
    samtools index validateContig31.bam &



     module add UHTS/Analysis/BEDTools/2.26.0;
     bedtools genomecov -ibam validateContig27.bam -bga -g Contig27.fasta > validateContig27.bedGraph &

     bedtools genomecov -ibam validateContig31.bam -bga -g Contig31.fasta > validateContig31.bedGraph &


Lets produce now 10'000 synthetic ones for each on top (~50x):

          module add UHTS/Analysis/BBMap/37.82;
          randomreads.sh ref=Contig27.fasta out=simulatedReads.Contig27.fq reads=10000 minlength=1000 maxlength=40000 pacbio=t pbmin=0.12 pbmax=0.17

          randomreads.sh ref=Contig31.fasta out=simulatedReads.Contig31.fq reads=10000 minlength=1000 maxlength=40000 pacbio=t pbmin=0.12 pbmax=0.17

Fastq is not useful for the further steps and larger - lets convert to FASTA:

         cat simulatedReads.Contig27.fq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > simulatedReads.Contig27.fasta
         cat simulatedReads.Contig31.fq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > simulatedReads.Contig31.fasta

Similar, re-mapping these and then mapping plots of coverage. Which should be pretty even in that case...

         minimap2 Contig31.fasta simulatedReads.Contig31.fasta -a -t 30 | samtools view -bhS | samtools sort -o validateContig31_sim.bam &
         minimap2 Contig27.fasta simulatedReads.Contig27.fasta -a -t 30 | samtools view -bhS | samtools sort -o validateContig27_sim.bam &
         samtools index validateContig31_sim.bam &
         samtools index validateContig27_sim.bam &
         bedtools genomecov -ibam validateContig27_sim.bam -bga -g Contig27.fasta > validateContig27_sim.bedGraph &
         bedtools genomecov -ibam validateContig31_sim.bam -bga -g Contig31.fasta > validateContig31_sim.bedGraph &

Here are the results:

![](coverage27.png)

In general it looks not too bad - maybe we will get that assembled as some of my reads seem to span the regions which have low coverage from the sampling. Likely towards the positions which have still small "N" - remember, I allowed <50bp gaps.

![](coverage31.png)

That one looks actually much better. Lets give it a first try to assemble now with miniasm the 2 contigs:

    cat Contig27.reads.fasta simulatedReads.Contig27.fasta > simulatedAndRealReads.Contig27.fasta &
    cat Contig31.reads.fasta simulatedReads.Contig31.fasta > simulatedAndRealReads.Contig31.fasta &




    module add UHTS/Analysis/minimap2/2.5;
    module add UHTS/Analysis/miniasm/0.2.r159.dirty;

    time minimap2 -x ava-pb -t 20 simulatedAndRealReads.Contig27.fasta simulatedAndRealReads.Contig27.fasta > simulatedAndRealReads.Contig27.paf &

    time minimap2 -x ava-pb -t 20 simulatedAndRealReads.Contig31.fasta simulatedAndRealReads.Contig31.fasta > simulatedAndRealReads.Contig31.paf &


    time miniasm -f simulatedAndRealReads.Contig31.fasta simulatedAndRealReads.Contig31.paf -s 1000 -R > simulatedAndRealReads.Contig31.miniasm.gfa &

    miniasm -f simulatedAndRealReads.Contig27.fasta simulatedAndRealReads.Contig27.paf -s 1000 -R > simulatedAndRealReads.Contig27.miniasm.gfa &


canu -d Contig27 -p test1 -assemble -pacbio-raw simulatedAndRealReads.Contig27.fasta genomeSize=600k useGrid=false minReadLength=500 minOverlapLength=200

canu -d Contig31 -p test1 -assemble -pacbio-raw simulatedAndRealReads.Contig31.fasta genomeSize=600k useGrid=false minReadLength=500 minOverlapLength=200

-----------------

    awk '/^S/{print ">"$2"\n"$3}' simulatedAndRealReads.Contig27.miniasm.gfa | fold > simulatedAndRealReads.Contig27.miniasm.fasta

    awk '/^S/{print ">"$2"\n"$3}' simulatedAndRealReads.Contig31.miniasm.gfa | fold > simulatedAndRealReads.Contig31.miniasm.fasta
