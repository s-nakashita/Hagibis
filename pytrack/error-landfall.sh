#!/bin/bash
init0=2019100700
init1=2019101212
CDIR=`pwd`

dh=$((12 * 3600))
./init.sh ${init0} ${init1} ${dh} > init_tmp.txt

function select_error() {
  #if [ -f ./${2}/error${1}_mod.txt ]; then
  #ln -fs ./${2}/error${1}_mod.txt error.txt
  #else
  ln -fs ./${2}/error${1}.txt error.txt
  #fi
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
  echo ${yyyy}-${mm}-${dd}T${hh} ${slp} >>${3} 
}

for center in ecmf rjtd kwbc egrr; do
  outfile=error-${center}.txt
  rm ${outfile}
  touch ${outfile}
  for init in $(cat init_tmp.txt); do
    select_error ${init} ${center} ${outfile}
  done
done

