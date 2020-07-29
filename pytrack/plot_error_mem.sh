#!/bin/sh
init=2019100912
orig=rjtd
outfile=$orig/error${init:4:9}.ps
vtime=$(date -jf "%Y%m%d%H%M%S" 20191012120000 +%s)

function plot_error() {
  MEM=1
  while test $MEM -le 26;do
    if [ $MEM -lt 10 ]; then
      MEM=0$MEM
    fi
    errortxt=jma/error${1}_${MEM}.txt
    if [ -f ${errortxt} ]; then
    time=$(date -jf "%Y%m%d%H%M%S" ${1}0000 +%s)
    maxdt=$(((${vtime}-${time})/3600))
    echo ${maxdt}
    awk -v max="${maxdt}" '{
      if ( $1 <= max ) print $1, $2; 
      }' ${errortxt} > tmp.txt
    echo $time
    pencol=lightgreen
    echo $pencol
    gmt psxy -R -J -O tmp.txt -i0,1 -W1p,${pencol} -K >> ${outfile} 
    gmt psxy -R -J -O tmp.txt -i0,1 -G${pencol} -Sc0.1 -Wfaint -K >> ${outfile} 
    fi
    MEM=`expr $MEM + 1`
  done
  errortxt=$orig/error${1}.txt
  if [ -f ${errortxt} ]; then
    time=$(date -jf "%Y%m%d%H%M%S" ${1}0000 +%s)
    maxdt=$(((${vtime}-${time})/3600))
    echo ${maxdt}
    awk -v max="${maxdt}" '{
      if ( $1 <= max ) print $1, $2; 
      }' ${errortxt} > tmp.txt
    echo $time
    pencol=red
    echo $pencol
    gmt psxy -R -J -O tmp.txt -i0,1 -W1p,${pencol} -K >> ${outfile} 
    gmt psxy -R -J -O tmp.txt -i0,1 -G${pencol} -Sc0.2 -Wfaint -K >> ${outfile} 
  fi
}


gmt psbasemap -R0/90/0/500 -JX15c/10c -Y8c -Bxa12df6hg6d -Bya100f2g50 -K > ${outfile}
plot_error ${init}

pstopdf ${outfile}
pdfcrop ${outfile%.ps}.pdf
mv -f ${outfile%.ps}-crop.pdf ${outfile%.ps}.pdf
