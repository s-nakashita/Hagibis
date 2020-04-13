#!/bin/bash
if [ $# -lt 3 ]; then
  echo "$0 yyyymmddhh0 yyyymmddhh1 dtsec"
  exit
fi
yyyy0=`echo ${1:0:4}`
mm0=`echo ${1:4:2}`
dd0=`echo ${1:6:2}`
hh0=`echo ${1:8:2}`
yyyy1=`echo ${2:0:4}`
mm1=`echo ${2:4:2}`
dd1=`echo ${2:6:2}`
hh1=`echo ${2:8:2}`
fmt=%Y%m%d%H%M%S
outfmt=+%s
init0=${1}
init1=${2}
#echo $init0
#echo $init1
dt=${3}
#echo $dt
t0=$(date -j -f "${fmt}" "${init0}0000"  "$outfmt" )
t1=$(date -j -f "${fmt}" "${init1}0000"  "$outfmt" )
t=${t0}
while [ ${t} -le ${t1} ]; do
    for center in ecmwf jma ncep ukmo; do
    	vt=$(date -r ${t} "+%Y%m%d%H" )
    	echo ${vt}_$center
    done
    t=$((t + dt))
done

