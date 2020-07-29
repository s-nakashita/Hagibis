#!/bin/bash
init0=2019100900
CDIR=`pwd`
bstfile=$CDIR/bst_hagibis.txt
outfile=track_jma_${init0}_mem.ps

function plot_track() {
	orig=2
    datadir=jma
	MEM=1
  	while test $MEM -le 26;do
      if [ $MEM -lt 10 ]; then
      	MEM=0$MEM
	  fi
	  tracktxt=$datadir/track${1}_${MEM}.txt
	  yyyy=`echo ${1:0:4}`
	  mm=`echo ${1:4:2}`
	  dd=`echo ${1:6:2}`
	  hh=`echo ${1:8:2}`
	  if [ -f ${tracktxt} ]; then
	    echo "exist"
	    ptime=$(date -jf "%Y/%m/%d %H:%M:%S" "2019/10/12 12:00:00" +"%s")
	    awk '{s=(100000-$7)/20000+0.1;print($5, $6,'${orig}', s > 0.1 ? s : 0.1)}' ${tracktxt} > tmp.txt
	    awk '!($3~12 && $4~12){s=(100000-$7)/20000+0.1;print($5, $6,'${orig}', s > 0.1 ? s : 0.1)}' ${tracktxt} > tmp_n.txt
	    awk '$3~12 && $4~12{s=(100000-$7)/20000+0.1;print($5, $6,'${orig}', s > 0.1 ? s : 0.1)}' ${tracktxt} > tmp_p.txt
	    pencol=lightgreen
	    gmt psxy -R -J -O tmp.txt -i0,1 -W1p,${pencol} -K >> ${outfile} 
	    gmt psxy -R -J -O tmp_n.txt -i0,1,2,3 -C$CDIR/track_mem.cpt -Scc -Wfaint -K >> ${outfile}
	    gmt psxy -R -J -O tmp_p.txt -i0,1,2,3 -C$CDIR/track_mem.cpt -SA -Wfaint -K >> ${outfile} 
	  fi
	  MEM=`expr $MEM + 1`
    done
	datadir=rjtd
	orig=1
	tracktxt=$datadir/track${1}.txt
	yyyy=`echo ${1:0:4}`
	mm=`echo ${1:4:2}`
	dd=`echo ${1:6:2}`
	hh=`echo ${1:8:2}`
	if [ -f ${tracktxt} ]; then
	  echo "exist"
	  ptime=$(date -jf "%Y/%m/%d %H:%M:%S" "2019/10/12 12:00:00" +"%s")
	  awk '{s=(100000-$7)/20000+0.1;print($5, $6,'${orig}', s > 0.1 ? s : 0.1)}' ${tracktxt} > tmp.txt
	  awk '!($3~12 && $4~12){s=(100000-$7)/20000+0.1;print($5, $6,'${orig}', s > 0.1 ? s : 0.1)}' ${tracktxt} > tmp_n.txt
	  awk '$3~12 && $4~12{s=(100000-$7)/20000+0.1;print($5, $6,'${orig}', s > 0.1 ? s : 0.1)}' ${tracktxt} > tmp_p.txt
	  pencol=red
	  gmt psxy -R -J -O tmp.txt -i0,1 -W1p,${pencol} -K >> ${outfile} 
	  gmt psxy -R -J -O tmp_n.txt -i0,1,2,3 -C$CDIR/track_mem.cpt -Scc -Wfaint -K >> ${outfile}
	  gmt psxy -R -J -O tmp_p.txt -i0,1,2,3 -C$CDIR/track_mem.cpt -SA -Wfaint -K >> ${outfile} 
	fi
}

gmt pscoast -R125/150/12/45 -JQ12c -Bag -Dh -Gburlywood -Sazure -Wthinnest -A100 -P -K > ${outfile}
#for init in $(cat $CDIR/init.txt); do
plot_track ${init0}
#done
#gmt psxy -R -J -O ${bstfile} -i4,5 -W3p -K >> ${outfile} 
#gmt psxy -R -J -O ${bstfile} -i4,5 -Sc0.1i -Gblack >> ${outfile} 
# best track
awk '$4 % 6 == 0 && !($3~12 && $4~12){s=(100000-$7)/20000+0.1;print($5, $6, s > 0.1 ? s : 0.1)}' ${bstfile} > tmp.txt
awk '$3~12 && $4~12{s=(100000-$7)/20000+0.1;print($5, $6, s > 0.1 ? s : 0.1)}' ${bstfile} > tmp_p.txt
#gmt psxy -R -J -O tmp.txt -i0,1 -W3p -K >> ${outfile} 
gmt psxy -R -J -O tmp.txt -i0,1,2 -Scc -W1p -K >> ${outfile} 
gmt psxy -R -J -O tmp_p.txt -i0,1,2 -SA -W1p >> ${outfile} 
#gmt pslegend -R -J -Dg135/20+w10 -O < legend.txt >> ${outfile}

ps2pdfwr ${outfile}
pdfcrop ${outfile%.ps}.pdf
mv -f ${outfile%.ps}-crop.pdf ${outfile%.ps}.pdf
