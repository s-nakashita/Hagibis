#!/bin/sh
orig=${1:-jma}
ntype=${2:-TE}
case $orig in
    ecmwf ) mem=50 ;;
    jma ) mem=26 ;;
    ncep ) mem=20 ;;
    ukmo ) mem=17 ;;
esac
mkdir -p $orig
cd $orig
idate=2019100912
edate=2019101212
smode=1
emode=1
nlev=6
CPPFLAGS="-cpp -E -P"
if [ $nlev -eq 6 ]; then
CPPFLAGS=${CPPFLAGS}" -Dlev6"
fi
echo $CPPFLAGS
cat > sens.nml << EOF
&sens_nml
orig = "${orig}",
mem = ${mem},
idate = ${idate},
edate = ${edate},
smode = ${smode},
emode = ${emode}
/
EOF
cat sens.nml
yyyy=${idate:0:4}
mm=${idate:4:2}
dd=${idate:6:2}
hh=${idate:8:2}
echo $yyyy $mm $dd $hh
em=${edate:4:2}
ed=${edate:6:2}
eh=${edate:8:2}
echo $em $ed $eh
datadir=/Volumes/dandelion/netcdf/tigge/${yyyy}/${orig}
if [ $smode -eq $emode ]; then
    header1=ensvsa-${ntype}-m${emode}
else
    header1=ensvsa-${ntype}-m${smode}-${emode}
fi
#header2=${orig}-${idate}-${edate}
#wfile=weight-${ntype}-${orig}-${idate}.grd
#logfile=${ntype}-${orig}-${idate}.log
header2=${orig}-${idate}-${edate}_nlev${nlev}
wfile=weight-${ntype}-${orig}-${idate}_nlev${nlev}.grd
logfile=${ntype}-${orig}-${idate}_nlev${nlev}.log
ofile=${header1}-${header2}-gr
cfile=${header1}-${header2}.ctl
ncfile=${header1}-${header2}.nc
mfile=glb_${idate}_mean.nc
rm -f glb_${idate}_*.nc
ln -s ${datadir}/glb_${idate}_mean.nc .

gfile=${header1}_t+wind_nlev${nlev}_${orig}_${mm}${dd}${hh}-${em}${ed}${eh}_isign
tefile=${header1}_EN_${orig}_${mm}${dd}${hh}-${em}${ed}${eh}
enfile=`echo ${header1} | sed -e "s/ensvsa-${ntype}//g"`-${header2}.txt
if [ ! -s $ofile ]; then
#    rm -f ${ofile} ${ncfile} ${wfile}
#fi
cd ..
if [ read_netcdf.F90 -nt tmp/read_netcdf.mod ]; then
    make clean
fi
make ${ntype} CPPFLAGS="${CPPFLAGS}" || exit 10
cd $orig
#rsync -auv ${datadir}/glb_${idate}_*.nc .
for m in $(seq 1 $mem);do
m2=`printf '%0.2d' $m`
ln -s ${datadir}/glb_${idate}_${m2}.nc
done

ln -s ../bin/ensvsa .
ln -s ../bin/grads-ensvsa .
ln -s ../bin/tevol-ensvsa .
./ensvsa                       || exit 11
./grads-ensvsa                 || exit 12
./tevol-ensvsa                 || exit 13
rm ensvsa grads-ensvsa tevol-ensvsa
# weights' ASCII text
awk '{if($1 == "vector="){$1="";print}}' $logfile > ${logfile%.*}.txt
fi
isec=$(date -jf "%Y%m%d%H" "${idate}" +"%s")
esec=$(date -jf "%Y%m%d%H" "${edate}" +"%s")
diff=$((${esec} - ${isec}))
#nd=$((${diff} / 43200 + 1))
nd=$((${diff} / 21600 + 1))
echo $nd
case $mm in
    01 ) mon=Jan ;;
    02 ) mon=Feb ;;
    03 ) mon=Mar ;;
    04 ) mon=Apr ;;
    05 ) mon=May ;;
    06 ) mon=Jun ;;
    07 ) mon=Jul ;;
    08 ) mon=Aug ;;
    09 ) mon=Sep ;;
    10 ) mon=Oct ;;
    11 ) mon=Nov ;;
    12 ) mon=Dec ;;
esac
if [ $nlev -eq 3 ]; then
nlev=5
plevs="200 250 300 500 850"
elif [ $nlev -eq 6 ]; then
nlev=8
plevs="200 250 300 500 700 850 925 1000"
fi
if [ ${ntype} = dTE ]; then
cat > $cfile << EOF
DSET ^${ofile}
UNDEF 9.999E+20
OPTIONS big_endian template
ydef 361 linear -90.0 0.5
xdef 720 linear 0.0 0.5
zdef ${nlev} levels ${plevs}
tdef ${nd} linear ${hh}Z${dd}${mon}${yyyy} 6hr
vars 9
ugrd ${nlev} 33,100,0 U component of wind [m/s]
vgrd ${nlev} 33,100,0 V component of wind [m/s]
t ${nlev} 33,100,0 Temperature [K]
gh ${nlev} 33,100,0 Geopotential Height [gpm]
te ${nlev} 33,100,0 Total moist energy [J/kg/m^2]
ke ${nlev} 33,100,0 Kinetic energy [J/kg/m^2]
pe ${nlev} 33,100,0 Potential energy (T) [J/kg/m^2]
pres_meansealev 0 33,100,0  Pressure at mean sea level [hPa]
peps 0 33,100,0 Potential energy (Ps) [J/kg/m^2]
endvars
EOF
else
cat > $cfile << EOF
DSET ^${ofile}
UNDEF 9.999E+20
OPTIONS big_endian template
ydef 361 linear -90.0 0.5
xdef 720 linear 0.0 0.5
zdef ${nlev} levels ${plevs}
tdef ${nd} linear ${hh}Z${dd}${mon}${yyyy} 6hr
vars 11
ugrd ${nlev} 33,100,0 U component of wind [m/s]
vgrd ${nlev} 33,100,0 V component of wind [m/s]
t ${nlev} 33,100,0 Temperature [K]
gh ${nlev} 33,100,0 Geopotential Height [gpm]
q ${nlev} 33,100,0 Specific Humidity [kg/kg]
te ${nlev} 33,100,0 Total moist energy [J/kg/m^2]
ke ${nlev} 33,100,0 Kinetic energy [J/kg/m^2]
pe ${nlev} 33,100,0 Potential energy (T) [J/kg/m^2]
le ${nlev} 33,100,0 Latent heat [J/kg/m^2]
pres_meansealev 0 33,100,0  Pressure at mean sea level [hPa]
peps 0 33,100,0 Potential energy (Ps) [J/kg/m^2]
endvars
EOF
fi
cat $cfile
cdo -f nc import_binary ${cfile} ${ncfile}
pwd
## plot
#if [ ${orig} = jma ]; then
    lfilter=False
#else
#    lfilter=True
#fi
dev=png
cat > config.ncl << EOF
ntype="${ntype}"
ofile="${ofile}"
ncfile="${ncfile}"
mfile="${mfile}"
gfile="${gfile}"
nlev=${nlev}
orig="${orig}"
yyyymmddhh="${idate}"
dev="${dev}"
lfilter=${lfilter}
EOF
cat config.ncl

### plotting perturbations
##for d in $(seq 0 $((${nd} - 1)));do
#ncl -nQ d=0 ../ensvsa.ncl
#ncl -nQ d=$((${nd} - 1)) ../ensvsa.ncl
##done
### plotting perturbation vorticity
#ncl -nQ d=0 lfilter=True ../ensvsa_vor.ncl
### plotting height-longitude sector
#for d in $(seq 0 3); do # 20.0 30.0 35.0; do
#latc=15.0
#ncl -nQ d=0 latc=${latc} ../ensvsa_h-lon.ncl
#latc=22.0
#ncl -nQ d=0 latc=${latc} ../ensvsa_h-lon.ncl
##ncl -nQ d=$((${nd} - 1)) latc=${latc} ../ensvsa_h-lon.ncl
#done
### plotting SLP & vertical-interpolated winds
ncl -nQ d=0 ../ensvsa_vint.ncl
ncl -nQ d=4 ../ensvsa_vint.ncl
ncl -nQ d=$((${nd} - 1)) ../ensvsa_vint.ncl
### plotting Energy distribution
for EN in te ke pe; do
out=`echo ${tefile} | sed -e "s/EN/${EN}/g"`
ncl -nQ nd=${nd} EN=\"${EN}\" out=\"${out}\" ../ensvsa-ENonly.ncl
done
if [ ${ntype} != dTE ]; then
EN=le
out=`echo ${tefile} | sed -e "s/EN/${EN}/g"`
ncl -nQ nd=${nd} EN=\"${EN}\" out=\"${out}\" ../ensvsa-ENonly.ncl
fi
### compute Energy vertical profile
#for EN in ke pe; do
#for latc in 15.0 25.0 35.0; do
#ncl -nQ latc=${latc} EN=\"${EN}\" ../ensvsa-ENavg.ncl
#done
#done
#if [ ${ntype} != dTE ]; then
#EN=le
#for latc in 15.0 25.0 35.0; do
#ncl -nQ latc=${latc} EN=\"${EN}\" ../ensvsa-ENavg.ncl
#done
#fi
#rm -f ${ofile} ${ncfile}
rm -f glb_${idate}_*.nc
ls -ltr | tail -9