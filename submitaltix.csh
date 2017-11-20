#$Id: submitaltix.csh,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
# usage: csh submitaltix.csh
bsub -W 6:00 -n 16 -q altix_s -o out.log -e err.log mpijob.sgi ./ABE > PI.log
