#!/bin/sh
init="2019100912"
#exp=est
for exp in clim est mgdsst; do
    echo $exp $init
    python rotate_phy3mpp_gsm.py $exp $init
    python create_netcdf_phy3mpp.py $exp $init
    rm np_fcst_phy3mpp_${init}.nc np_fcst_phy2m_avr_${init}.nc
done