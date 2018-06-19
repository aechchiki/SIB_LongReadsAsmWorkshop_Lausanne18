## Working on the cluster

You can find [here](https://edu.sib.swiss/course/view.php?name=Vital-IT_Infrastructure_Usage) some documentation describing in detail the usage of the cluster.

A few **important** points:

 - do not execute programs directly on the cluster ! (exceptions are e.g. `wget; mkdir; mv; cp`)
 - running jobs can be checked via `bjobs` and `bjobs -a` to show finished jobs as well
 - unintended jobs can be killed via `bkill jobID`
 - if you want to kill all your jobs `bkill -u userID 0`
 - add `-o output.txt -e error.txt` to capture errors and additional information, but give meaningful names!
 - sending many jobs will decrease your priority, thus verify your commands before you send them

# Setting up your working directory

On Vital-IT, you must read and write files in the "scratch" directory:

```
cd /scratch/beegfs/weekly/
```

Be careful, files in this folder are erased after a month! At the end of the practicals, if you want to keep your results, you need to back-up the data, for example in a compressed tarball that you move to your home or archive folder, or to another computer.

Create your own directory:

```
mkdir <username>; cd <username>
```

You will always be working from this directory. Before launching commands please be sure that you are located in the right directory (`pwd`). The expected output should always be: `/scratch/beegfs/monthly/<username>` or any nested folders.

# Submitting commands to the cluster

It is **NOT** allowed to launch any calculation on the frontal machine where you are connected (`prd.vital-it.ch`), except for very light jobs.

You need to submit each job for batched execution through a job scheduler that will dispatch it on the cluster nodes (LSF).

Use for example:

```
bsub "<command line>"
```

Or, better, write your commands in a script and submit it with:

```
bsub < script.sh
```

## Useful links

Have a look at [this short tutorial](http://bioinfo.unil.ch/tp/SIB_LongReadsWorkshop_Bern16/vital-it-usage.html) to help you write such a script yourself.

To save you some typing, we also provide you a [skeleton](documentation/skeleton.sh) for a submission script.

You can specify the BSUB parameters according [LSF BSUB parameters](https://www.ibm.com/support/knowledgecenter/en/SSETD4_9.1.2/lsf_command_ref/bsub.1.html).

You can monitor jobs according [LSF rules](https://www.ibm.com/support/knowledgecenter/en/SSETD4_9.1.2/lsf_kc_cmd_ref.html).

# Navigation

Go to Next section: [Pre-assembly steps](documentation/bio.md)

Get a [Unix refresh](documentation/unix.md)

Return to [Main page](README.md)
