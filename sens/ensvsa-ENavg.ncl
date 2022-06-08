;----------------------------------------------------------------------
; vector_1.ncl
;
; Concepts illustrated:
;   - Overlaying vectors and filled contours on a map
;   - Manually attaching lat/lon coordinate arrays to a variable
;   - Changing the length of the smallest vector as a fraction of the reference vector
;   - Moving the vector reference annotation to the top right of the plot
;   - Drawing curly vectors
;   - Setting the color for vectors
;   - Making the labelbar be vertical
;   - Increasing the thickness of vectors
;----------------------------------------------------------------------
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
  load "config.ncl"
begin

;---Open netCDF file
  a = addfile(ncfile,"r")
;
; Read in Sea Surface Temperature Anomalies
; Read in U and V at 1000 mb [subscript 0]
; Determine the subscripts corresponding to Jan 88
;
  lon0 = 120.0
  lon1 = 160.0
  lat0 = latc-5.0
  lat1 = latc+5.0
  out = EN + "_avr_latc" + sprinti("%d",toint(latc)) + ".mtx"
  if (latc.eq.999) then ; targeted
  lon0 = 137.0
  lon1 = 142.0
  lat0 = 33.0
  lat1 = 37.0
  out = EN + "_avr_target.mtx"
  end if
  te3d     = a->$EN$(:,:,{lat0:lat1},{lon0:lon1})
  if (EN .eq. "pe") then
    ps = a->peps(:,{lat0:lat1},{lon0:lon1})
    te3d(:,{1000},:,:) = te3d(:,{1000},:,:) + ps 
  end if
  print(te3d&lev)
  lat      = a->lat({lat0:lat1})
  lon      = a->lon({lon0:lon1})
  plev     = a->lev
  time     = a->time
  ntime = dimsizes(time)
  nlev  = dimsizes(plev)
  nlon  = dimsizes(lon)
  wgt = cos(lat*3.1415/180.0)
  area = sum(wgt)*nlon
  print(area)
  te_ave = dim_sum_n(te3d*conform(te3d,wgt,2),(/2,3/))
  te_ave = te_ave / area
  dp   = dpres_plevel_Wrap(plev, plev(nlev-1), plev(0), 0)
  dpall= 1000.0
  sig  = dble2flt( dp / dpall )
  te_ave = te_ave*conform(te_ave,sig,1)
  te_int = dim_sum_n(te_ave,1)
  printMinMax(te_ave, 0)
  print(te_int)
  fmt_str = "%s" + str_concat(conform_dims(nlev,"%10.3f",-1)) + "%s"
  dlist = NewList("lifo")
  ListPush(dlist,"int")
  do j=nlev-1,0,1
    ListPush(dlist,(/plev(j)/))
  end do
  ListPush(dlist,"plev")
  write_table(out,"w",dlist,fmt_str)
  delete(dlist)
  fmt_str = "%s" + str_concat(conform_dims(nlev+1,"%10.3f",-1))
  dlist = NewList("lifo")
  ListPush(dlist,(/te_int(:)/))
  do j=nlev-1,0,1
    ListPush(dlist,(/te_ave(:,j)/))
  end do
  row_labels = sprinti("%4i",ispan(0,ntime-1,1)*6)
  ListPush(dlist,row_labels)
  write_table(out,"a",dlist,fmt_str)
  ;te = te3d(:,0,:,:)
  ;te3d = te3d * conform_dims(dimsizes(te3d),sig,1)
  ;te = dim_sum_n(te3d,1)
  ;printMinMax(te(0,{lat0:lat1},{lon0:lon1}),0)
  ;if ((EN.eq."pe").or.(EN.eq."te")) then 
  ;  peps   = a->peps(:,:,:)
  ;  printMinMax(peps(0,{lat0:lat1},{lon0:lon1}),0)
  ;  te     = te + peps
  ;end if
end
