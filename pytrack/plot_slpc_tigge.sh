#!/bin/bash
init0=2019101200
init1=2019101212
CDIR=`pwd`

rm init_*.txt
dh=$((12 * 3600))
./makecpt.sh ${init0} ${init1} ${dh}
./makecpt2.sh
./init.sh ${init0} ${init1} ${dh} > init_tmp.txt
./init2.sh ${init0} ${init1} ${dh} > init_c.txt
./legend2.sh
bstfile=bst_hagibis.txt
outfile=slpc_tigge_${init0:4:4}.ps

function plot_slpc() {
    line=$2
    for center in ecmf rjtd kwbc egrr; do
	ln -s $CDIR/$center/track${1}.txt track.txt
	tracktxt=track.txt
	yyyy=`echo ${1:0:4}`
	mm=`echo ${1:4:2}`
	dd=`echo ${1:6:2}`
	hh=`echo ${1:8:2}`
	echo $line
	if [ -f ${tracktxt} ]; then
	    echo "exist $center"
	    time=$(date -j -f "%Y%m%d%H%M%S" "${1}0000" "+%s")
	    awk '{print $1 "-" $2 "-" $3 "T" $4, $7*0.01, '${time}'}' ${tracktxt} > tmp.txt
	    pencol=$(cat track_tigge.cpt | awk -v awk_var=$line '{if(NR == awk_var) {print $2}}' )
	    echo $pencol
	    gmt psxy -R -J -O tmp.txt -i0,1 -W1p,${pencol} -K >> ${outfile} 
	    if [ $line -lt 5 ]; then
		gmt psxy -R -J -O tmp.txt -i0,1,2 -Sc0.1i -Wfaint -G${pencol} -K >> ${outfile}
	    else
		gmt psxy -R -J -O tmp.txt -i0,1,2 -St0.1i -Wfaint -G${pencol} -K >> ${outfile}
	    fi
	fi
	line=`expr $line + 1`
	rm track.txt
    done
}


gmt psbasemap -R2019-10-04T18/2019-10-15T00/912/1012 -JX15c/10c -Y8c -Bxa1df6hg1d -Bya4f2g4 -K > ${outfile}
count=1
for init in $(cat init_tmp.txt); do
    echo ${init} ${count}
    plot_slpc ${init} ${count}
    count=`expr $count + 4`
done
echo $?
# best track
awk '$4 % 6 == 0{print $1 "-" $2 "-" $3 "T" $4, $7*0.01}' ${bstfile} > tmp.txt
gmt psxy -R -J -O tmp.txt -i0,1 -W3p -K >> ${outfile}
echo $?
gmt psxy -R -J -O tmp.txt -i0,1 -Sc0.1i -Gwhite -Wthin -K >> ${outfile}
echo $?
gmt pslegend -Dx0c/-11c+w5c/10c -K -O < legend2.txt >> ${outfile}
#head -5 legend2.txt  | gmt pslegend -Dx0c/-11c+w5c/10c -K -O >> ${outfile}
#echo $?
#tail -n +4 legend2.txt | head -6 | gmt pslegend -Dx5c/-11c+w5c/10c -K -O >> ${outfile}
#echo $?
#tail -n +12 legend2.txt | gmt pslegend -Dx10c/-11c+w5c/10c -O >> ${outfile}
echo $?
pstopdf ${outfile}
pdfcrop ${outfile%.ps}.pdf
mv -f ${outfile%.ps}-crop.pdf ${outfile%.ps}.pdf
