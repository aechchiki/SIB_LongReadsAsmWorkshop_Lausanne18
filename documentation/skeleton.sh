#!/bin/bash

#BSUB -L /bin/bash
#BSUB -J my_job_name

# Add the commands to be executed below:
# ...

# Global variables can be defined, and reused later in the script, for example:
work_dir=/scratch/beegfs/weekly/$USER
mkdir -p $work_dir/test/
