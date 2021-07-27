#!/bin/zsh
alias gmt
CDIR=`pwd`
#datadir=rjtd
datadir=`pwd` #jma
bstfile=$CDIR/bst_hagibis.txt
outfile=track_gsm_0900-0912.ps

cd $datadir
pwd
function plot_track() {
    tracktxt=track${1}_gsm.txt
    #tracktxt=track_merra2.txt
    yyyy=`echo ${1:0:4}`
    mm=`echo ${1:4:2}`
    dd=`echo ${1:6:2}`
    hh=`echo ${1:8:2}`
    if [ -f ${tracktxt} ]; then
	echo "exist"
	time=$(date -jf "%Y/%m/%d %H:%M:%S" "$yyyy/$mm/$dd $hh:00:00" +%s)
	awk '{s=(100000-$7)/20000+0.1;print($5, $6,'${time}', s > 0.1 ? s : 0.1)}' ${tracktxt} > tmp.txt
	pencol=$(awk '$1~'${time}'{print $2}' $CDIR/track.cpt)
	gmt psxy -R -J -O tmp.txt -i0,1 -W1p,${pencol} -K >> ${outfile} 
	gmt psxy -R -J -O tmp.txt -i0,1,2,3 -C$CDIR/track.cpt -Scc -Wfaint -K >> ${outfile} 
    fi
#  rm -f tmp.txt
}

gmt pscoast -R125/165/12/45 -JQ12c -Bag -Dh -Gburlywood -Sazure -Wthinnest -A100 -P -K > ${outfile}
#for init in $(cat $CDIR/init.txt); do
#  plot_track ${init}
#done
plot_track 2019100900
plot_track 2019100912
#gmt psxy -R -J -O ${bstfile} -i4,5 -W3p -K >> ${outfile} 
#gmt psxy -R -J -O ${bstfile} -i4,5 -Sc0.1i -Gblack >> ${outfile} 
# best track
awk '$4 % 6 == 0{s=(100000-$7)/20000+0.1;print($5, $6, s > 0.1 ? s : 0.1)}' ${bstfile} > tmp.txt
gmt psxy -R -J -O tmp.txt -i0,1 -W3p -K >> ${outfile} 
gmt psxy -R -J -O tmp.txt -i0,1,2 -Scc -W1p >> ${outfile} 
#gmt psxy -R -J -O tmp.txt -i0,1,2 -Scc -W1p -K >> ${outfile} 
#gmt pslegend -R -J -Dg135/20+w10 -O < legend.txt >> ${outfile}

ps2pdfwr ${outfile}
pdfcrop ${outfile%.ps}.pdf
mv -f ${outfile%.ps}-crop.pdf ${outfile%.ps}.pdf
