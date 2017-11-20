#$Id: splotalldat.sh,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
#!/bin/bash

file="temp.plt"

echo "set term gif" > ${file}
for f in *.dat; do
   echo "set output \"${f%.dat}.gif\"" >> ${file}
   echo "splot \"$f\" using 1:2:3 w l title \"\"" >> ${file}
done

gnuplot ${file}

rm ${file}
