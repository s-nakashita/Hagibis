#!/bin/sh
init0=2019100900
init1=2019100912

dh=$((12 * 3600))
./makecpt.sh ${init0} ${init1} ${dh}
./init.sh ${init0} ${init1} ${dh} > init.txt
./legend.sh track.cpt init.txt
bstfile=bst_hagibis.txt
gsmfile=track2019100900_gsm.txt
orig=rjtd
outfile=$orig/slpc1009+gsm.ps

function plot_slpc() {
  tracktxt=$orig/track${1}.txt
  if [ -f ${tracktxt} ]; then
    time=$(date -jf "%Y%m%d%H%M%S" ${1}0000 +%s)
    awk '{print $1 "-" $2 "-" $3 "T" $4, $7*0.01, '${time}'}' ${tracktxt} > tmp.txt
    pencol=$(awk '$1~'${time}'{print $2}' track.cpt)
    echo $pencol
    gmt psxy -R -J -O tmp.txt -i0,1 -W1p,${pencol} -K >> ${outfile} 
    gmt psxy -R -J -O tmp.txt -i0,1,2 -Ctrack.cpt -Sc0.1i -Wfaint -K >> ${outfile} 
  fi
}


#gmt psbasemap -R2019-10-04T18/2019-10-15T00/912/1012 -JX15c/10c -Y8c -Bxa1df6hg1d -Bya4f2g4 -K > ${outfile}
gmt psbasemap -R2019-10-08T18/2019-10-13T00/912/1012 -JX15c/10c -Y8c -Bxa1df6hg1d -Bya4f2g4 -K > ${outfile}
# best track
awk '$4 % 6 == 0{print $1 "-" $2 "-" $3 "T" $4, $7*0.01}' ${bstfile} > tmp.txt
gmt psxy -R -J -O tmp.txt -i0,1 -W3p -K >> ${outfile}
gmt psxy -R -J -O tmp.txt -i0,1 -Sc0.1i -Gwhite -Wthin -K >> ${outfile}
for init in $(cat init.txt); do
  plot_slpc ${init}
done
# GSM analysis
awk '$4 % 6 == 0{print $1 "-" $2 "-" $3 "T" $4, $7*0.01}' ${gsmfile} > tmp.txt
gmt psxy -R -J -O tmp.txt -i0,1 -W1p,blue -K >> ${outfile}
gmt psxy -R -J -O tmp.txt -i0,1 -Sc0.1i -Gblue -Wfaint -K >> ${outfile}
echo "S 0c c 0.1i blue 1p,blue 0.5c GSM analysis" >> legend.txt

#gmt pslegend -R -J -Dn0.05/0.6+w0.2/0.5 -K -O < legend.txt >> ${outfile}
head -6 legend.txt  | gmt pslegend -Dx0c/-11c+w5c/10c -K -O >> ${outfile}
tail -n +7 legend.txt | head -6 | gmt pslegend -Dx5c/-11c+w5c/10c -K -O >> ${outfile}
tail -n +13 legend.txt | gmt pslegend -Dx10c/-11c+w5c/10c -O >> ${outfile}

pstopdf ${outfile}
pdfcrop ${outfile%.ps}.pdf
mv -f ${outfile%.ps}-crop.pdf ${outfile%.ps}.pdf
