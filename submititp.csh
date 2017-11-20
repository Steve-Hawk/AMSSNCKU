#$Id: submititp.csh,v 1.1.1.1 2016/05/30 06:22:13 zjcao Exp $
# usage: qsub submititp.csh
# for itp's cluster
#PBS -N zjcao
#PBS -l nodes=4:ppn=32
#PBS -l walltime=240:00:00
#PBS -q N_large
echo "This jobs is "$PBS_JOBID@$PBS_QUEUE
cd $PBS_O_WORKDIR
mpirun -np 128 -machinefile $PBS_NODEFILE ./ABE > out.log
