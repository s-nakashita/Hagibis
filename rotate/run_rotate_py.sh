#!/bin/sh
set -e
if [ $# -lt 2 ]; then
    echo "Usage :: $0 yyyymmddhh center"
    exit
fi

yyyymmddhh=${1}
yyyy=${yyyymmddhh:0:4}
nlon=${nlon:-360}
nlat=${nlat:-17}
latmax=${latmax:-8}

center=${2}
case $center in
    ecmwf ) orig=ecmf ; member=50;;
    jma ) orig=rjtd ; member=26;;
    ncep ) orig=kwbc ; member=20;;
    ukmo ) orig=egrr ; member=17;;
esac
#orig=${center}
trackf=../pytrack/$orig/track${yyyymmddhh}.txt
#trackf=../pytrack/$center/track${yyyymmddhh}_mean.txt
#trackf=../pytrack/$orig/track_anl.txt
datadir=/Volumes/dandelion/netcdf/tigge/${yyyy}/$center
outdir=/Volumes/dandelion/netcdf/tigge/${yyyy}/$center/rotate
if [ ! -d $outdir ]; then
    mkdir $outdir
fi
for mem in $(seq 1 $member); do
if [ $mem -lt 10 ];then
mem=0$mem
fi
echo $yyyymmddhh $datadir $trackf $nlon $nlat $latmax $center $mem | python rotate_scalar.py || exit 10

echo $yyyymmddhh $datadir $trackf $nlon $nlat $latmax $center $mem | python rotate_uv_v2.py || exit 11

echo $yyyymmddhh $center $mem | python create_netcdf.py || exit 12

mv np_${yyyymmddhh}_$mem.nc ${outdir}/
rm np_*_${yyyymmddhh}_$mem.nc
done
# mean
echo $yyyymmddhh $datadir $trackf $nlon $nlat $latmax $center | python rotate_scalar.py || exit 13

echo $yyyymmddhh $datadir $trackf $nlon $nlat $latmax $center | python rotate_uv_v2.py || exit 14

echo $yyyymmddhh $center | python create_netcdf.py || exit 15

mv np_${yyyymmddhh}_mean.nc ${outdir}/
rm np_*_${yyyymmddhh}_mean.nc