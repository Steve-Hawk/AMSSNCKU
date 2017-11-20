#$Id: pickdata.sh,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
#!/bin/bash

if [ "$#" -lt 3 ]
  then echo "usage: sh pickdata.sh everyline inputfile outputfile"
       exit
fi

j=$1
inputfile=$2
outputfile=$3

if [ ! -f "${inputfile}" ]
then echo "${inputfile} does not exits"
     exit
fi

k=0
while [ $k -lt $j ]; do
k=`expr $k + 1`
if [ -f "${outputfile}$k.dat" ]
then echo "${outputfile}$k.dat already exits"
     exit
fi
done
k=0;

while   read   data;   do
     if [ "${data}" = "" ]
       then k=`expr $k + 1`
	    while [ $k -le $j ]; do
	       echo $data >> ${outputfile}$k.dat
               k=`expr $k + 1`
            done
	    k=0;
     elif [ ${data:0:1} = "#" ] 
       then k=`expr $k + 1`
	    while [ $k -le $j ]; do
	       echo $data >> ${outputfile}$k.dat
               k=`expr $k + 1`
            done
	    k=0;
     else k=`expr $k + 1`
          echo $data >> ${outputfile}$k.dat
          if [ "$k" -eq "$j" ]
             then k=0
          fi
     fi
done   < ${inputfile}
