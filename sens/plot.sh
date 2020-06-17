#!/bin/bash
if [ $# -lt 1 ]; then
    echo "Usage :: $0 yyyymmddhh"
    exit
fi

yyyymmddhh=${1}
outfile=mTE-jma-${yyyymmddhh}.ps
intxt=mTE-1mode-jma-${yyyymmddhh}.txt
gmt psxy $intxt -R0.0/170.0/0.01/3000.0 -Jx0.1/1l -Bg6.0a12.0:"hour":/g2a1:"Total Energy":WeSn -W3,green -K >$outfile
for m in $(seq 1 1 26);do
    nm=$(printf %0.2d $m)
    intxt=mTE-${nm}-jma-${yyyymmddhh}.txt
    gmt psxy -R -J -O $intxt -W1 -K >> $outfile
done
intxt=mTE-mean-jma-${yyyymmddhh}.txt
gmt psxy -R -J -O $intxt -W3,red >> $outfile

pstopdf ${outfile}
pdfcrop ${outfile%.ps}.pdf
mv ${outfile%.ps}-crop.pdf ${outfile%.ps}.pdf
