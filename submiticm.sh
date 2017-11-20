#$Id: submiticm.sh,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
# usage: sh submiticm.sh
#!/bin/sh
NAME=`pwd`
bsub -np 256 -o out.log -N ${NAME##/*/} mpijob ABE xyzinput
