#!/bin/bash
if [ $# -lt 3 ]; then
  echo "$0 yyyymmddhh0 yyyymmddhh1 dtsec"
  exit
fi
fmt=%Y%m%d%H%M%S
yyyy0=`echo ${1:0:4}`
mm0=`echo ${1:4:2}`
dd0=`echo ${1:6:2}`
yyyy1=`echo ${2:0:4}`
mm1=`echo ${2:4:2}`
dd1=`echo ${2:6:2}`
outfmt=+%s
init0=${1}
init1=${2}
dt=${3}
echo $init0 $init1 $dt
t0=$(date -j -f "$fmt" "${init0}0000" "${outfmt}")
t1=$(($(date -j -f "$fmt" "${init1}0000" "${outfmt}") + ${dt}))
echo $t0
echo $t1
palette=polar
cptfile=track2.cpt
gmt makecpt -C${palette} -T${t0}/${t1}/${dt} > ${cptfile}
