#!/bin/sh
head=error-slp
outfile=${head}-tigge.ps

if [ "$head" = "error" ]; then
gmt psbasemap -R2019-10-07T00/2019-10-13T00/0/600 -JX15c/10c -Y8c -BWeSn -Bxa1df6hg1d+l"initial date" -Bya100f100g50+l"error(km)" -K > ${outfile}
else
gmt psbasemap -R2019-10-07T00/2019-10-13T00/0/40 -JX15c/10c -Y8c -BWeSn -Bxa1df6hg1d+l"initial date" -Bya10f10g5+l"error(hPa)" -K > ${outfile}
fi
# mean error
#awk '$4 % 6 == 0{print $1 "-" $2 "-" $3 "T" $4, $7*0.01}' ${bstfile} > tmp.txt
ln -fs ./${head}-mean.txt tmp.txt
gmt psxy -R -J -O tmp.txt -i0,1 -W3p -K >> ${outfile}
gmt psxy -R -J -O tmp.txt -i0,1 -Sc0.1i -Gwhite -Wthin -K >> ${outfile}

line=1
for orig in ecmf rjtd kwbc egrr; do
  #if [ $orig = ecmf ]; then
  #ln -fs ./${head}-${orig}_mod.txt tmp.txt
  #else
  ln -fs ./${head}-${orig}.txt tmp.txt
  #fi
  pencol=$(cat track_error_tigge.cpt | awk -v awk_var=$line '{if(NR == awk_var) {print $2}}' )
  gmt psxy -R -J -O tmp.txt -i0,1 -W1p,${pencol} -K >> ${outfile} 
  gmt psxy -R -J -O tmp.txt -i0,1 -G${pencol} -Sc0.2 -Wfaint -K >> ${outfile} 
  line=$((${line} + 1))
  rm tmp.txt
done

gmt pslegend -R -J -Dn0.7/0.7+w4c -F+gwhite+p -K -O < legend_error_tigge.txt >> ${outfile}
#head -6 legend.txt  | gmt pslegend -Dx0c/-11c+w5c/10c -K -O >> ${outfile}
#tail -n +2 legend.txt | head -6 | gmt pslegend -Dx5c/-11c+w5c/10c -K -O >> ${outfile}
#tail -n +13 legend.txt | gmt pslegend -Dx10c/-11c+w5c/10c -O >> ${outfile}

convert -trim -density 300 -rotate 90 ${outfile} ${outfile%.ps}.png
pstopdf ${outfile}
pdfcrop ${outfile%.ps}.pdf
mv -f ${outfile%.ps}-crop.pdf ${outfile%.ps}.pdf
