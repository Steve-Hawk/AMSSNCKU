#$Id: submitnchc.csh,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
#!/bin/csh
# Command File name = submitnchc.csh
# Usage : llsubmit submitnchc.csh
# check job: llq -f %jn %id %o %st %dd -u u60hjy00
# delete job: llcancel 450383
#
#@ job_name = schwarzschild
#@ executable = /usr/bin/poe
#@ arguments = /work2/u60ccc01/run6/ABE
#@ output = out.log
#@ error = err.log
#@ job_type = parallel
#@ network.MPI = csss,shared,US
#@ node = 1
#@ tasks_per_node = 16
#@ class = 16cpu
#@ notify_user = zjcao@amt.ac.cn
#  notification = always
#  notification = never
#@ notification = error
#@ wall_clock_limit = 72:00:00
#@ queue
