#$Id: submitbnu.csh,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
# usage: qsub submitbnu.csh
# for altix330
##PBS -N gj
##PBS -l select=1:ncpus=1:cluster=false
#cd $PBS_O_WORKDIR
#mpirun -np 1 ./ABE
# for chess
#PBS -N e10
#PBS -l select=1:ncpus=4:cluster=true
cd $PBS_O_WORKDIR
mpirun.ch_gm -np 4 ./ABE > out.log
