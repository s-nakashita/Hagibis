#!/bin/bash
init0=2019100512
init1=2019101212
CDIR=`pwd`

dh=$((12 * 3600))
./init.sh ${init0} ${init1} ${dh} > init_tmp.txt

function select_error() {
  case ${2} in
    ecmwf ) orig=ecmf ;;
    jma ) orig=rjtd ;;
    ncep ) orig=kwbc ;;
    ukmo ) orig=egrr ;;
  esac
  ln -fs ./${orig}/error${1}.txt error.txt
  yyyy=`echo ${1:0:4}`
	mm=`echo ${1:4:2}`
	dd=`echo ${1:6:2}`
	hh=`echo ${1:8:2}`
  time=$(date -jf "%Y%m%d%H%M%S" "${1}0000" "+%s")
  vtime=$(date -jf "%Y%m%d%H%M%S" "20191012120000" "+%s")
  fsec=$(($vtime - $time))
  ft=$(($fsec / 3600))
  echo $fsec $ft
  slp=$(cat error.txt | awk -v ft=$ft '{if($1 == ft) {print $2}}')
  echo 00 ${slp} >>${4}
  M=${3}
  MEM=1
  while [ $MEM -le $M ];do
  if [ $MEM -lt 10 ]; then
    MEM=0$MEM
  fi
  ln -fs ./${2}/error${1}_${MEM}.txt error.txt
  yyyy=`echo ${1:0:4}`
	mm=`echo ${1:4:2}`
	dd=`echo ${1:6:2}`
	hh=`echo ${1:8:2}`
  time=$(date -jf "%Y%m%d%H%M%S" "${1}0000" "+%s")
  vtime=$(date -jf "%Y%m%d%H%M%S" "20191012120000" "+%s")
  fsec=$(($vtime - $time))
  ft=$(($fsec / 3600))
  echo $fsec $ft
  slp=$(cat error.txt | awk -v ft=$ft '{if($1 == ft) {print $2}}')
  echo ${MEM} ${slp} >>${4}
  MEM=`expr $MEM + 1`
  done
}

for center in ecmwf ncep ukmo; do
  for init in 2019100900 2019100912; do
    case $center in
      ecmwf ) M=50 ;;
      jma   ) M=26 ;;
      ncep  ) M=20 ;;
      ukmo  ) M=17 ;;
    esac
    outfile=error${init}-${center}.txt
    rm ${outfile}
    touch ${outfile}
    select_error ${init} ${center} ${M} ${outfile}
  done
done

