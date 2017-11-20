#!/bin/bash
#PBS -q n_cpu
#PBS -N test_bbh
#PBS -l nodes=4:ppn=32
#PBS -l walltime=142000 
#PBS -o test_bbh.out
#PBS -e test_bbh.err

cd $PBS_O_WORKDIR

OMP_NUM_THREADS=4
NC=`cat $PBS_NODEFILE | wc -l`
NN=`cat $PBS_NODEFILE | uniq | wc -l`
NP=$((NC/OMP_NUM_THREADS))

hosts=`cat $PBS_NODEFILE | uniq`
for i in $hosts
do
   for((j=1;j<=NP/NN;j++))
   do
       echo $i
   done
done > mymachinefile

#$mpirun -np $NP -machinefile mymachinefile -x OMP_NUM_THREADS=$OMP_NUM_THREADS ./cosmomc hede.ini$
mpirun -np 128 -machinefile mymachinefile ./ABE > out.log 




