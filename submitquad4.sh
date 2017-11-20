#$Id: submitquad4.sh,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
# usage: sh submiquad4.sh
#!/bin/sh
mpirun -np 4 -machinefile ../machines ABE > out.log &
