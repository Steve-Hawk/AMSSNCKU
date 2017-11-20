#$Id: splotallbin.sh,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
#!/bin/bash

file="temp.plt"

echo "set term gif" > ${file}
echo "set xrange[-10:10]" >> ${file}
echo "set yrange[-10:10]" >> ${file}
echo "set zrange[-0.8:0.8]" >> ${file}
for f in *.bin; do
   ~/ReadData XY $f
done

for f in Lev06-00*.dat; do
   echo "set output \"${f%.dat}.gif\"" >> ${file}
   fn=${f#*Lev06-00}
   str="splot"
   for f1 in *${fn}; do
       str="${str} \"${f1}\" using 1:2:3 w l title \"\","
   done
   str=${str%,}
   echo ${str} >> ${file}
done


gnuplot ${file}

# thumb in rule: keep original files, delete all files created by this routine
rm ${file}

for f in *.bin; do
   rm -f ${f%.bin}.dat
done
