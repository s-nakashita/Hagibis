#!/bin/sh
outfile=error-tigge.ps

gmt psbasemap -R2019-10-05T00/2019-10-13T00/0/600 -JX15c/10c -Y8c -BWeSn -Bxa1df6hg1d+l"initial date" -Bya100f100g50+l"error(km)" -K > ${outfile}
# mean error
#awk '$4 % 6 == 0{print $1 "-" $2 "-" $3 "T" $4, $7*0.01}' ${bstfile} > tmp.txt
ln -fs ./error-mean.txt tmp.txt
gmt psxy -R -J -O tmp.txt -i0,1 -W3p -K >> ${outfile}
gmt psxy -R -J -O tmp.txt -i0,1 -Sc0.1i -Gwhite -Wthin -K >> ${outfile}

line=1
for orig in ecmf rjtd kwbc egrr; do
  ln -fs ./error-${orig}.txt tmp.txt
  pencol=$(cat track_error_tigge.cpt | awk -v awk_var=$line '{if(NR == awk_var) {print $2}}' )
  gmt psxy -R -J -O tmp.txt -i0,1 -W1p,${pencol} -K >> ${outfile} 
  gmt psxy -R -J -O tmp.txt -i0,1 -G${pencol} -Sc0.2 -Wfaint -K >> ${outfile} 
  line=$((${line} + 1))
done

gmt pslegend -R -J -Dn0.7/0.7+w4c -F+gwhite+p -K -O < legend_error_tigge.txt >> ${outfile}
#head -6 legend.txt  | gmt pslegend -Dx0c/-11c+w5c/10c -K -O >> ${outfile}
#tail -n +2 legend.txt | head -6 | gmt pslegend -Dx5c/-11c+w5c/10c -K -O >> ${outfile}
#tail -n +13 legend.txt | gmt pslegend -Dx10c/-11c+w5c/10c -O >> ${outfile}

pstopdf ${outfile}
pdfcrop ${outfile%.ps}.pdf
mv -f ${outfile%.ps}-crop.pdf ${outfile%.ps}.pdf
