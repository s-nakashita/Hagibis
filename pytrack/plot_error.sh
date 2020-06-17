#!/bin/sh
init0=2019100900
init1=2019100912

dh=$((12 * 3600))
./makecpt.sh ${init0} ${init1} ${dh}
./init.sh ${init0} ${init1} ${dh} > init.txt
./legend.sh
bstfile=bst_hagibis.txt
orig=rjtd
outfile=$orig/error${init0:4:9}-${init1:4:9}.ps
vtime=$(date -jf "%Y%m%d%H%M%S" 20191012120000 +%s)

function plot_error() {
  errortxt=$orig/error${1}.txt
  if [ -f ${errortxt} ]; then
    time=$(date -jf "%Y%m%d%H%M%S" ${1}0000 +%s)
    maxdt=$(((${vtime}-${time})/3600))
    echo ${maxdt}
    awk -v max="${maxdt}" -v dt="${2}" '{
      if ( $1 <= max ) print $1+dt, $2; 
      }' ${errortxt} > tmp.txt
    #ft=$(awk '{print $1}' ${errortxt})
    #echo ${ft}
    #sec=`expr ${ft} \* 3600 + ${time}`
    #fdate=$(date -jf "%s" ${sec} +%Y-%m-%dT%H)
    #awk -v fd="${fdate}" '{print fd,$2}' ${errortxt} > tmp.txt
    #awk '{print $1 "-" $2 "-" $3 "T" $4, $7*0.01, '${time}'}' ${errortxt} > tmp.txt
    echo $time
    pencol=$(awk '$1~'${time}'{print $2}' track2.cpt)
    echo $pencol
    gmt psxy -R -J -O tmp.txt -i0,1 -W1p,${pencol} -K >> ${outfile} 
    gmt psxy -R -J -O tmp.txt -i0,1 -G${pencol} -Sc0.2 -Wfaint -K >> ${outfile} 
  fi
}


gmt psbasemap -R0/90/0/300 -JX15c/10c -Y8c -Bxa12df6hg6d -Bya100f2g50 -K > ${outfile}
dt=0
for init in $(cat init.txt); do
  plot_error ${init} ${dt}
  dt=`expr $dt + 12`
done
# best track
#awk '$4 % 6 == 0{print $1 "-" $2 "-" $3 "T" $4, $7*0.01}' ${bstfile} > tmp.txt
#gmt psxy -R -J -O tmp.txt -i0,1 -W3p -K >> ${outfile}
#gmt psxy -R -J -O tmp.txt -i0,1 -Sc0.1i -Gwhite -Wthin -K >> ${outfile}

#gmt pslegend -R -J -Dn0.05/0.6+w0.2/0.5 -K -O < legend.txt >> ${outfile}
#head -6 legend.txt  | gmt pslegend -Dx0c/-11c+w5c/10c -K -O >> ${outfile}
tail -n +2 legend.txt | head -6 | gmt pslegend -Dx5c/-11c+w5c/10c -K -O >> ${outfile}
#tail -n +13 legend.txt | gmt pslegend -Dx10c/-11c+w5c/10c -O >> ${outfile}

pstopdf ${outfile}
pdfcrop ${outfile%.ps}.pdf
mv -f ${outfile%.ps}-crop.pdf ${outfile%.ps}.pdf
