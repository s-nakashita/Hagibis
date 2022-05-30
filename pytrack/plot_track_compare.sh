#!/bin/bash
#init0=2019100900
init0=${1:-2019100900}
CDIR=`pwd`
bstfile=$CDIR/bst_hagibis.txt
outfile=track_all_${init0}_mod.ps

function plot_track() {
    for orig in 0 1 2 3;do
	case $orig in
	    0 ) datadir=ecmf ;;
	    1 ) datadir=rjtd ;;
	    2 ) datadir=kwbc ;;
	    3 ) datadir=egrr ;;
	esac
	#if [ $orig -eq 0 ];then
	tracktxt=$datadir/track${1}_mod.txt
	#else
	#tracktxt=$datadir/track${1}.txt
	#fi
	yyyy=`echo ${1:0:4}`
	mm=`echo ${1:4:2}`
	dd=`echo ${1:6:2}`
	hh=`echo ${1:8:2}`
	if [ -f ${tracktxt} ]; then
	    echo "exist"
	    ptime=$(date -jf "%Y/%m/%d %H:%M:%S" "2019/10/12 12:00:00" +"%s")
	    awk '{s=(100000-$7)/20000+0.1;print($5, $6,'${orig}', s > 0.1 ? s : 0.1)}' ${tracktxt} > tmp.txt
	    awk '!($3~12 && $4~12){s=(100000-$7)/20000+0.1;print($5, $6, '${orig}', s > 0.1 ? s : 0.1)}' ${tracktxt} > tmp_n.txt
	    awk '$3~12 && $4~12{print($5, $6, '${orig}', 0.4)}' ${tracktxt} > tmp_p.txt
	    pencol=$(awk '$1~'${orig}'{print $2}' $CDIR/track_tigge.cpt)
	    gmt5 psxy -R -J -O tmp.txt -i0,1 -W1p,${pencol} -K >> ${outfile} 
	    gmt5 psxy -R -J -O tmp_n.txt -i0,1,2,3 -C$CDIR/track_tigge.cpt -Scc -Wfaint -K >> ${outfile}
	    gmt5 psxy -R -J -O tmp_p.txt -i0,1,2,3 -C$CDIR/track_tigge.cpt -SA -Wfaint -K >> ${outfile} 
	fi
#	datadir=rjtd
#	# 00 UTC
#	hh=00
#	orig=1
#	tracktxt=$datadir/track${1}${hh}.txt
#	yyyy=`echo ${1:0:4}`
#	mm=`echo ${1:4:2}`
#	dd=`echo ${1:6:2}`
#	#hh=`echo ${1:8:2}`
#	if [ -f ${tracktxt} ]; then
#	    echo "exist"
#	    ptime=$(date -jf "%Y/%m/%d %H:%M:%S" "2019/10/12 12:00:00" +"%s")
#	    awk '{s=(100000-$7)/20000+0.1;print($5, $6,'${orig}', s > 0.1 ? s : 0.1)}' ${tracktxt} > tmp.txt
#	    awk '!($3~12 && $4~12){s=(100000-$7)/20000+0.1;print($5, $6, '${orig}', s > 0.1 ? s : 0.1)}' ${tracktxt} > tmp_n.txt
#	    awk '$3~12 && $4~12{print($5, $6, '${orig}', 0.4)}' ${tracktxt} > tmp_p.txt
#	    pencol=$(awk '$1~'${orig}'{print $2}' $CDIR/track_tigge.cpt)
#	    gmt5 psxy -R -J -O tmp.txt -i0,1 -W1p,${pencol},- -K >> ${outfile} 
#	    gmt5 psxy -R -J -O tmp_n.txt -i0,1,3 -Scc -W1p,${pencol} -K >> ${outfile}
#	    gmt5 psxy -R -J -O tmp_p.txt -i0,1,3 -SA -W1p,${pencol} -K >> ${outfile} 
#	fi
#	# 12 UTC
#	hh=12
#	orig=1
#	tracktxt=$datadir/track${1}${hh}.txt
#	if [ -f ${tracktxt} ]; then
#	    echo "exist"
#	    ptime=$(date -jf "%Y/%m/%d %H:%M:%S" "2019/10/12 12:00:00" +"%s")
#	    awk '{s=(100000-$7)/20000+0.1;print($5, $6,'${orig}', s > 0.1 ? s : 0.1)}' ${tracktxt} > tmp.txt
#	    awk '!($3~12 && $4~12){s=(100000-$7)/20000+0.1;print($5, $6,'${orig}', s > 0.1 ? s : 0.1)}' ${tracktxt} > tmp_n.txt
#	    awk '$3~12 && $4~12{print($5, $6,'${orig}', 0.4)}' ${tracktxt} > tmp_p.txt
#	    pencol=$(awk '$1~'${orig}'{print $2}' $CDIR/track_tigge.cpt)
#	    gmt5 psxy -R -J -O tmp.txt -i0,1 -W1p,${pencol} -K >> ${outfile} 
#	    gmt5 psxy -R -J -O tmp_n.txt -i0,1,2,3 -C$CDIR/track_tigge.cpt -Scc -Wfaint -K >> ${outfile}
#	    gmt5 psxy -R -J -O tmp_p.txt -i0,1,2,3 -C$CDIR/track_tigge.cpt -SA -Wfaint -K >> ${outfile} 
#	fi
	#orig=0
	#tracktxt=jma/track${1}_mean.txt
	#yyyy=`echo ${1:0:4}`
	#mm=`echo ${1:4:2}`
	#dd=`echo ${1:6:2}`
	#hh=`echo ${1:8:2}`
	#if [ -f ${tracktxt} ]; then
	#    echo "exist"
	#    ptime=$(date -jf "%Y/%m/%d %H:%M:%S" "2019/10/12 12:00:00" +"%s")
	#    awk '{s=(100000-$7)/20000+0.1;print($5, $6,'${orig}', s > 0.1 ? s : 0.1)}' ${tracktxt} > tmp.txt
	#    awk '!($3~12 && $4~12){s=(100000-$7)/20000+0.1;print($5, $6,'${orig}', s > 0.1 ? s : 0.1)}' ${tracktxt} > tmp_n.txt
	#    awk '$3~12 && $4~12{s=(100000-$7)/20000+0.1;print($5, $6,'${orig}', s > 0.1 ? s : 0.1)}' ${tracktxt} > tmp_p.txt
	#    pencol=$(awk '$1~'${orig}'{print $2}' $CDIR/track_tigge.cpt)
	#    gmt5 psxy -R -J -O tmp.txt -i0,1 -W1p,${pencol} -K >> ${outfile} 
	#    gmt5 psxy -R -J -O tmp_n.txt -i0,1,2,3 -C$CDIR/track_tigge.cpt -Scc -Wfaint -K >> ${outfile}
	#    gmt5 psxy -R -J -O tmp_p.txt -i0,1,2,3 -C$CDIR/track_tigge.cpt -SA -Wfaint -K >> ${outfile} 
	#fi
	#orig=2
	#tracktxt=jma/gtrack${1}_mean.txt
	#yyyy=`echo ${1:0:4}`
	#mm=`echo ${1:4:2}`
	#dd=`echo ${1:6:2}`
	#hh=`echo ${1:8:2}`
	#if [ -f ${tracktxt} ]; then
	#    echo "exist"
	#    ptime=$(date -jf "%Y/%m/%d %H:%M:%S" "2019/10/12 12:00:00" +"%s")
	#    awk '{s=(100000-$7)/20000+0.1;print($5, $6,'${orig}', s > 0.1 ? s : 0.1)}' ${tracktxt} > tmp.txt
	#    awk '!($3~12 && $4~12){s=(100000-$7)/20000+0.1;print($5, $6,'${orig}', s > 0.1 ? s : 0.1)}' ${tracktxt} > tmp_n.txt
	#    awk '$3~12 && $4~12{s=(100000-$7)/20000+0.1;print($5, $6,'${orig}', s > 0.1 ? s : 0.1)}' ${tracktxt} > tmp_p.txt
	#    pencol=$(awk '$1~'${orig}'{print $2}' $CDIR/track_tigge.cpt)
	#    gmt5 psxy -R -J -O tmp.txt -i0,1 -W1p,${pencol} -K >> ${outfile} 
	#    gmt5 psxy -R -J -O tmp_n.txt -i0,1,2,3 -C$CDIR/track_tigge.cpt -Scc -Wfaint -K >> ${outfile}
	#    gmt5 psxy -R -J -O tmp_p.txt -i0,1,2,3 -C$CDIR/track_tigge.cpt -SA -Wfaint -K >> ${outfile} 
	#fi
	rm -f tmp.txt
	done
}

gmt5 pscoast -R130/150/12/45 -JQ12c -Bag -Dh -Gburlywood -Sazure -Wthinnest -A100 -P -K > ${outfile}
#gmt psxy -R -J -O ${bstfile} -i4,5 -W3p -K >> ${outfile} 
#gmt psxy -R -J -O ${bstfile} -i4,5 -Sc0.1i -Gblack >> ${outfile} 
# best track
awk '$4 % 6 == 0 && !($3~12 && $4~12){s=(100000-$7)/20000+0.1;print($5, $6, s > 0.1 ? s : 0.1)}' ${bstfile} > tmp.txt
awk '$3~12 && $4~12{print($5, $6, 0.4)}' ${bstfile} > tmp_p.txt
#gmt psxy -R -J -O tmp.txt -i0,1 -W3p -K >> ${outfile} 
gmt5 psxy -R -J -O tmp.txt -i0,1,2 -Scc -W1p -K >> ${outfile} 
gmt5 psxy -R -J -O tmp_p.txt -i0,1,2 -SA -W1p -K >> ${outfile} 
#gmt pslegend -R -J -Dg135/20+w10 -O < legend.txt >> ${outfile}
#for init in $(cat $CDIR/init.txt); do
plot_track ${init0}
#done

convert ${outfile} -trim ${outfile%.ps}.png
ps2pdfwr ${outfile}
pdfcrop ${outfile%.ps}.pdf
mv -f ${outfile%.ps}-crop.pdf ${outfile%.ps}.pdf
