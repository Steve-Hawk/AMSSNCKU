#$Id: submitncku.sh,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
# usage: qsub submitncku.sh
#!/bin/sh
# -----------------------------------------------------------------------------
# Define the running shell to N1GE
#$ -S /bin/bash
# -----------------------------------------------------------------------------
# Define the running directory for your job.
# Cwd stands for current working directory
#$ -cwd
# -----------------------------------------------------------------------------
# Define the name of job (by default, the job name is set as the script name.
#$ -N zjcao
# -----------------------------------------------------------------------------
# Define the standard and error output files for your job.
# -e specifies the stderr file, -o for the stdout file, and
# -j merges both stderr and stdout to the stdout file.
#$ -j y
#$ -o out.log
# -----------------------------------------------------------------------------
# Specify memory to be required for this job. M stands for megabytes, or G for gigabytes.
# By default, N1GE will reserve 2G for this job.
#$ -l h_vmem=32G
# -----------------------------------------------------------------------------
# Request the Parallel Environment with 16 slots (computing cores)
#$ -pe pgi_mpich 4

#  -m be -M zjcao@amt.ac.cn
#$ -v MPI_HOME=/ap/pgi/linux86-64/6.2/mpi/mpich,PGI=/ap/pgi
#  -v MPI_HOME=/opt/vltmpi/OPENIB/mpi.icc.rsh , PGI=/ap/pgi
 
# needs in $NSLOTS         : the number of tasks to be used
#          $TMPDIR/machines: a valid machine file to be passed to mpirun
/bin/hostname
cd /home/lun1/hpc0001/zjcao/nn

/ap/pgi/linux86-64/6.2/mpi/mpich/bin/mpirun -np $NSLOTS -machinefile $TMPDIR/machines ABE xyzinput

#/bin/date
