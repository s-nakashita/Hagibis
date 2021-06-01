#!/bin/sh
orig=jma
orig=${1}
case $orig in
    ecmwf ) mem=50 ;;
    jma ) mem=26 ;;
    ncep ) mem=20 ;;
    ukmo ) mem=17 ;;
esac
idate=2019100912
edate=2019101212
smode=1
emode=1
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
ntype=TE
ntype=${2}
yyyy=${idate:0:4}
mm=${idate:4:2}
dd=${idate:6:2}
hh=${idate:8:2}
echo $yyyy $mm $dd $hh
em=${edate:4:2}
ed=${edate:6:2}
eh=${edate:8:2}
echo $em $ed $eh
if [ $smode -eq $emode ]; then
    ofile=ensvsa-${ntype}-m${emode}-${orig}-${idate}-${edate}_n-gr
    cfile=ensvsa-${ntype}-m${emode}-${orig}-${idate}-${edate}_n.ctl
    ncfile=ensvsa-${ntype}-m${emode}-${orig}-${idate}-${edate}_n.nc
    gfile=ensvsa-${ntype}_m${emode}_t+wind_n_${orig}_${mm}${dd}${hh}-${em}${ed}${eh}_mono
    tefile=ensvsa-${ntype}_m${emode}_te_${orig}_${mm}${dd}${hh}-${em}${ed}${eh}
else
    ofile=ensvsa-${ntype}-m${smode}-${emode}-${orig}-${idate}-${edate}_n-gr
    cfile=ensvsa-${ntype}-m${smode}-${emode}-${orig}-${idate}-${edate}_n.ctl
    ncfile=ensvsa-${ntype}-m${smode}-${emode}-${orig}-${idate}-${edate}_n.nc
    gfile=ensvsa-${ntype}_m${smode}-${emode}_t+wind_n_${orig}_${mm}${dd}${hh}-${em}${ed}${eh}
    tefile=ensvsa-${ntype}_m${smode}-${emode}_te_${orig}_${mm}${dd}${hh}-${em}${ed}${eh}
fi
if [ -e $ofile ]; then
    rm ${ofile}
    rm weight-${ntype}-${orig}-${idate}_n.grd
fi
if [ read_netcdf.f90 -nt read_netcdf.mod ]; then
    make clean
fi
make ${ntype}
./ensvsa-${ntype}
./grads-ensvsa-${ntype}

isec=$(date -jf "%Y%m%d%H" "${idate}" +"%s")
esec=$(date -jf "%Y%m%d%H" "${edate}" +"%s")
diff=$((${esec} - ${isec}))
nd=$((${diff} / 43200 + 1))
#nd=$((${diff} / 21600 + 1))
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
if [ ${ntype} = TE ]; then
cat > $cfile << EOF
DSET ^${ofile}
UNDEF 9.999E+20
OPTIONS big_endian template
ydef 361 linear -90.0 0.5
xdef 720 linear 0.0 0.5
zdef 3 levels 300 500 850
tdef ${nd} linear ${hh}Z${dd}${mon}${yyyy} 12hr
vars 6
ugrd 3 33,100,0 U component of wind [m/s]
vgrd 3 33,100,0 V component of wind [m/s]
t 3 33,100,0 Temperature [K]
q 3 33,100,0 Specific Humidity [kg/kg]
pres_meansealev 0 33,100,0  Pressure at mean sea level [hPa]
te 0 33,100,0 Total energy [J/kg]
endvars
EOF
elif [ ${ntype} = dTE ]; then
cat > $cfile << EOF
DSET ^${ofile}
UNDEF 9.999E+20
OPTIONS big_endian template
ydef 361 linear -90.0 0.5
xdef 720 linear 0.0 0.5
zdef 3 levels 300 500 850
tdef ${nd} linear ${hh}Z${dd}${mon}${yyyy} 12hr
vars 5
ugrd 3 33,100,0 U component of wind [m/s]
vgrd 3 33,100,0 V component of wind [m/s]
t 3 33,100,0 Temperature [K]
pres_meansealev 0 33,100,0  Pressure at mean sea level [hPa]
te 0 33,100,0 Total energy [J/kg]
endvars
EOF
fi
cat $cfile
cdo -f nc import_binary ${cfile} ${ncfile}
if [ ${orig} = jma ]; then
    lfilter=False
else
    lfilter=True
fi
dev=pdf
ncl -nQ d=0 dev=\"${dev}\" lfilter=${lfilter} ncfile=\"${ncfile}\" gfile=\"${gfile}\" orig=\"${orig}\" yyyymmddhh=\"${idate}\" ../ncl/ensvsa-${ntype}.ncl
#ncl -nQ d=$((${nd} - 1)) dev=\"${dev}\" lfilter=${lfilter} ncfile=\"${ncfile}\" gfile=\"${gfile}\" orig=\"${orig}\" yyyymmddhh=\"${idate}\" ../ncl/ensvsa-${ntype}.ncl
#ncl -nQ d=2 dev=\"${dev}\" lfilter=${lfilter} ncfile=\"${ncfile}\" gfile=\"${gfile}\" orig=\"${orig}\" yyyymmddhh=\"${idate}\" ../ncl/ensvsa-${ntype}.ncl
ncl -nQ d=0 dev=\"${dev}\" lfilter=${lfilter} ncfile=\"${ncfile}\" gfile=\"${gfile}\" orig=\"${orig}\" yyyymmddhh=\"${idate}\" ../ncl/ensvsa-${ntype}_vint.ncl
ncl -nQ d=2 dev=\"${dev}\" lfilter=${lfilter} ncfile=\"${ncfile}\" gfile=\"${gfile}\" orig=\"${orig}\" yyyymmddhh=\"${idate}\" ../ncl/ensvsa-${ntype}_vint.ncl
ncl -nQ d=$((${nd} - 1)) dev=\"${dev}\" lfilter=${lfilter} ncfile=\"${ncfile}\" gfile=\"${gfile}\" orig=\"${orig}\" yyyymmddhh=\"${idate}\" ../ncl/ensvsa-${ntype}_vint.ncl
ncl -nQ nd=${nd} dev=\"${dev}\" ncfile=\"${ncfile}\" out=\"${tefile}\" ../ncl/ensvsa-${ntype}only.ncl
ls -ltr | tail -9