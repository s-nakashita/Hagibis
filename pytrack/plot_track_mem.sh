#!/bin/bash
GMT=gmt6
init0=${1} # yyyymmddhh
CDIR=`pwd`
bstfile=$CDIR/bst_hagibis.txt
center=${2} # ecmwf, jma, ncep, ukmo
outfile=track_${center}_${init0}_mem.ps
rm -rf tmp.txt
function plot_track() {
	orig=2
	center=${2}
    datadir=${center}
	MEM=1
	M=${3}
	b=${4}
	w=${5}
  	while test $MEM -le $M;do
	  if [ $MEM -eq $b ]; then
	    orig=3
	  elif [ $MEM -eq $w ]; then
	    orig=4
	  else
		orig=2
	  fi
#		cptfile=$CDIR/track_mem_${center}.cpt
      if [ $MEM -lt 10 ]; then
      	MEM=0$MEM
	  fi
	  cptfile=$CDIR/track_mem.cpt
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
	    awk '$3~12 && $4~12{s=(100000-$7)/20000+0.1;print($5, $6,'${orig}', s > 0.2 ? s : 0.2)}' ${tracktxt} > tmp_p.txt
	    #awk '$3~12 && $4~12{print($5, $6,'${orig}', 0.5)}' ${tracktxt} > tmp_p.txt
	    pencol=$(awk '$1~'${orig}'{print $2}' ${cptfile})
		echo $pencol
	    #pencol=lightgreen
	    $GMT psxy -R -J -O tmp.txt -i0,1 -W1p,${pencol} -K >> ${outfile} 
	    #$GMT psxy -R -J -O tmp_n.txt -i0,1,2,3 -C$cptfile -Scc -Wfaint -K >> ${outfile}
	    $GMT psxy -R -J -O tmp_p.txt -i0,1,2,3 -C$cptfile -SA -Wfaint -K >> ${outfile} 
	  fi
	  MEM=`expr $MEM + 1`
    done
#	for MEM in $b $w; do
#	  if [ $MEM -eq $b ]; then 
#	    pencol=green4
#		cptfile=$CDIR/track_mem_b.cpt
#	  elif [ $MEM -eq $w ]; then
#	    pencol=purple
#		cptfile=$CDIR/track_mem_w.cpt
#	  fi
#      if [ $MEM -lt 10 ]; then
#      	MEM=0$MEM
#	  fi
#	  tracktxt=$datadir/track${1}_${MEM}.txt
#	  yyyy=`echo ${1:0:4}`
#	  mm=`echo ${1:4:2}`
#	  dd=`echo ${1:6:2}`
#	  hh=`echo ${1:8:2}`
#	  if [ -f ${tracktxt} ]; then
#	    echo "exist"
#	    ptime=$(date -jf "%Y/%m/%d %H:%M:%S" "2019/10/12 12:00:00" +"%s")
#	    awk '{s=(100000-$7)/20000+0.1;print($5, $6,'${orig}', s > 0.1 ? s : 0.1)}' ${tracktxt} > tmp.txt
#	    awk '!($3~12 && $4~12){s=(100000-$7)/20000+0.1;print($5, $6,'${orig}', s > 0.1 ? s : 0.1)}' ${tracktxt} > tmp_n.txt
#	    awk '$3~12 && $4~12{print($5, $6,'${orig}', 0.5)}' ${tracktxt} > tmp_p.txt
#	    #pencol=lightgreen
#	    $GMT psxy -R -J -O tmp.txt -i0,1 -W2p,${pencol} -K >> ${outfile} 
#	    #$GMT psxy -R -J -O tmp_n.txt -i0,1,2,3 -C$cptfile -Scc -Wfaint -K >> ${outfile}
#	    $GMT psxy -R -J -O tmp_p.txt -i0,1,2,3 -C$cptfile -SA -Wfaint -K >> ${outfile} 
#	  fi
#	done
	case $center in
		jma ) datadir=rjtd; orig=1 ;;
		ecmwf ) datadir=ecmf; orig=0 ;;
		ncep ) datadir=kwbc; orig=2 ;;
		ukmo ) datadir=egrr; orig=3 ;;
	esac
	#orig=1
	tracktxt=$datadir/track${1}.txt
	cptfile=${CDIR}/track_tigge.cpt
	#cptfile=${CDIR}/track_mem_${center}.cpt
	yyyy=`echo ${1:0:4}`
	mm=`echo ${1:4:2}`
	dd=`echo ${1:6:2}`
	hh=`echo ${1:8:2}`
	if [ -f ${tracktxt} ]; then
	  echo "exist"
	  ptime=$(date -jf "%Y/%m/%d %H:%M:%S" "2019/10/12 12:00:00" +"%s")
	  awk '{s=(100000-$7)/20000+0.1;print($5, $6,'${orig}', s > 0.1 ? s : 0.1)}' ${tracktxt} > tmp.txt
	  awk '!($3~12 && $4~12){s=(100000-$7)/20000+0.1;print($5, $6,'${orig}', s > 0.1 ? s : 0.1)}' ${tracktxt} > tmp_n.txt
	  awk '$3~12 && $4~12{s=(100000-$7)/20000+0.1;print($5, $6,'${orig}', s > 0.2 ? s : 0.2)}' ${tracktxt} > tmp_p.txt
	  #awk '$3~12 && $4~12{print($5, $6,'${orig}', 0.5)}' ${tracktxt} > tmp_p.txt
	  pencol=$(awk '$1~'${orig}'{print $2}' $cptfile)
	  echo $pencol
	  $GMT psxy -R -J -O tmp.txt -i0,1 -W1p,${pencol} -K >> ${outfile} 
	  $GMT psxy -R -J -O tmp_n.txt -i0,1,2,3 -C$cptfile -Scc -Wfaint -K >> ${outfile}
	  $GMT psxy -R -J -O tmp_p.txt -i0,1,2,3 -C$cptfile -SA -Wfaint -K >> ${outfile} 
	fi
}

$GMT pscoast -R130/150/10/40 -JQ12c -Bag -Dh -Gburlywood -Sazure -Wthinnest -A100 -P -K > ${outfile}
#for init in $(cat $CDIR/init.txt); do
b=-1
w=-1
case $center in
    ecmwf ) M=50 ;; # b=26 ; w=15 ;;
    jma   ) M=26 ;; # b=5  ; w=26 ;; #12UTC
    #jma   ) M=26 ; b=15 ; w=19 ;; #00UTC
    ncep  ) M=20 ;; # b=16 ; w=1  ;;
    ukmo  ) M=17 ;; # b=2  ; w=3  ;;
esac
plot_track ${init0} ${center} ${M} ${b} ${w}
#done
#$GMT psxy -R -J -O ${bstfile} -i4,5 -W3p -K >> ${outfile} 
#$GMT psxy -R -J -O ${bstfile} -i4,5 -Sc0.1i -Gblack >> ${outfile} 
# best track
awk '$4 % 6 == 0 && !($3~12 && $4~12){s=(100000-$7)/20000+0.1;print($5, $6, s > 0.1 ? s : 0.1)}' ${bstfile} > tmp.txt
awk '$3~12 && $4~12{s=(100000-$7)/20000+0.1;print($5, $6, s > 0.2 ? s : 0.2)}' ${bstfile} > tmp_p.txt
#$GMT psxy -R -J -O tmp.txt -i0,1 -W3p -K >> ${outfile} 
$GMT psxy -R -J -O tmp.txt -i0,1,2 -Scc -W1p -K >> ${outfile} 
$GMT psxy -R -J -O tmp_p.txt -i0,1,2 -SA -W1p >> ${outfile} 
#$GMT pslegend -R -J -Dg135/20+w10 -O < legend.txt >> ${outfile}

convert ${outfile} -trim -density 300 ${outfile%.*}.png
ps2pdfwr ${outfile}
pdfcrop ${outfile%.ps}.pdf
mv -f ${outfile%.ps}-crop.pdf ${outfile%.ps}.pdf
