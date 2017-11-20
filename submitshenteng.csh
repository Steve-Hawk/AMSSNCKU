#$Id: submitshenteng.csh,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
# usage: csh submitshenteng.csh
#bsub -W 6:00 -a intelmpi -n 32 -q x64_blades -o out.log -e err.log mpirun.lsf ./ABE > PI.log
bsub -W 144:00 -a intelmpi -n 8 -q x64_3950small -o out.log -e err.log mpirun.lsf ./ABE > PI.log
