#$Id: submitibm1350.csh,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
#!/bin/sh
# LoadLeveler JCF file for runing an MPICH job
# Command File name = submitibm1350.csh
# Usage : llsubmit submitibm1350.csh
# check job: llq -f %jn %id %o %st %dd -u u60hjy00
# delete job: llcancel 450383
#
#@ job_name = schwarzschild
#@ initialdir = /work3/u60hjy00/zjcao
#@ output = out.log
#@ error = err.log
#@ job_type = parallel
#@ node = 8
#@ tasks_per_node = 4
#@ class = 32cpu
#@ wall_clock_limit = 96:00:00
#@ notification = error
#@ notify_user = zjcao@amt.ac.cn
#@ queue
#
echo "------------------------------------------------------------"
echo LOADL_STEP_ID=$LOADL_STEP_ID
echo "------------------------------------------------------------"
# Make sure that the ll_get_machine_list binary is accessible on all machines
# in the LoadLeveler cluster.
/opt/LoadL/bin/ll_get_machine_list > /tmp/machinelist.$LOADL_STEP_ID
#
machine_count=`cat /tmp/machinelist.$LOADL_STEP_ID | wc -l`
echo $machine_count
echo MachineList:
cat /tmp/machinelist.$LOADL_STEP_ID
echo "------------------------------------------------------------"
mpirun -np $machine_count -machinefile \
/tmp/machinelist.$LOADL_STEP_ID ABE
#
rm /tmp/machinelist.$LOADL_STEP_ID
