#!/bin/sh
model=${model:-40r1v2}
init=${init:-2018070112}
ntrunc=${ntrunc:-1023}
nlev=${nlev:-60}

ctldir=/Volumes/Pegasus/tctrack/201303Yagi/forecast/iECMWF-OpenIFS/40r1v2/TL1023L60/2013060912/rotate
datadir=/Volumes/Pegasus/OpenIFS/runs/${model}/TL${ntrunc}L${nlev}/${init}
dd=${init:6:2}
hh=${init:8:2}
grdate=$(LANG=C date -jf %Y%m%d%H%M%S ${init}0000 +%H:%MZ%d%b%Y) 

trackfile=${datadir}/track.txt
nl=$(wc -l ${trackfile} | awk '{print $1}')

for ctl in ${ctldir}/*.ctl; do
  sed -e "s/tdef 21/tdef ${nl}/" -e "s/12:00Z09JUN2013/${grdate}/" ${ctl} > ${datadir}/rotate/$(basename ${ctl})
done
