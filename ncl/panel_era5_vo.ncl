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

begin

;---Open netCDF file
  ;hh = 12
  ;cmap = read_colormap_file("precip3_16lev") ; all level
  cmap = read_colormap_file("GMT_polar") ; low level
  ;cmap = cmap(::-1,:) ; for dyn
  wks = gsn_open_wks("png","VOR_ERA5-GSM_0900")
  plot = new(6,graphic)

  idy = 9
  ihr = 0
  d = 0
  id = 0
  ;do id = 0, 5
  idate = 2019100000 + idy*100 + ihr
  yyyy = str_get_cols(sprinti("%0.10d",idate), 0, 3)
  mm   = str_get_cols(sprinti("%0.10d",idate), 4, 5)
  dd   = str_get_cols(sprinti("%0.10d",idate), 6, 7)
  yyyymmddhh = sprinti("%0.10d",idate)
  print(yyyymmddhh)
  a1 = addfile("/Volumes/dandelion/netcdf/era5/"+yyyy+"/"+mm+"/"+dd+"/plev_"+dd+".nc", "r")
  lat     = a1->latitude
  lon     = a1->longitude
  time3d  = a1->time
  date3d    = cd_string(time3d, "%N%D %H:%M")
  print(date3d)
;
  ;pv = a1->pv(d,{500},:,:)
  u  = a1->u(d,:,:,:)
  v  = a1->v(d,:,:,:)
  ;t  = a1-> t(d,{500},:,:)
  vo = a1->vo(d,:,:,:)
  ;lev= dble2flt(a1->level({500}))
  ;lev=lev*1e2 
  ;lev@units = "Pa"
  ;pt = pot_temp(lev, t, -1, False)
  ;PV = pot_vort_isobaric(lev, u, v, t, t&latitude, 1, 0)
  ;pv = PV({500},:,:)
;
  a2 = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est/fcst_mdlpp_"+yyyymmddhh+".nc", "r")
  a3 = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est/fcst_mdl_pvpp_"+yyyymmddhh+".nc", "r")
  lat2     = a2->lat
  lon2     = a2->lon

  ;pv2 = a2->PV(d,{500},:,:)
  u2  = a2->U(d,:,:,:)
  v2  = a2->V(d,:,:,:)
  vo2 = a3->VOR(d,:,:,:)
;---Create plot
  res                  = True               ; plot mods desired
  res@gsnDraw	      = False
  res@gsnFrame	      = False
;     res@cnInfoLabelOn	      = False
     
; Set resources for tickmarks
  res@gsnStringFont              = "helvetica-bold"  ; will be used for gsnCenterString later
  res@gsnCenterStringFontHeightF = 0.04              ; ditto
  res@gsnCenterStringOrthogonalPosF = 0.04
  res@gsnCenterStringParallelPosF = 0.4
  res@gsnCenterStringFontHeightF = 0.03
  ;res@gsnCenterString	   = date3d(d)
  res@gsnLeftString         = ""
  res@gsnAddCyclic	      = False
  res@cnFillOn             = True               ; turn on color for contours
  res@cnLinesOn            = True              ; turn off contour lines
  res@cnLineLabelsOn       = False              ; turn off contour line labels
  ;res@gsnScalarContour     = True               ; (for gsn_csm_vector_scalar_map)draw contours
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnFillPalette        = cmap

  res@mpFillOn             = False 
;     res@mpLandFillColor      = "gray"            ; set land to be gray
  res@mpMinLonF            =  110.               ; select a subregion
  res@mpMaxLonF            =  150.
  res@mpMinLatF            =  10. 
  res@mpMaxLatF            =  40.

  ;res@lbLabelBarOn         = False
  ;res@lbOrientation            = "Vertical"     ; vertical label bar
  ;res@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer

; note, when doing a subregion, NCL determines the range of the data from
; the full domain. If you wish to just consider the domain you are plotting,
; you must manually set those levels.

  res@cnMinLevelValF       = -30.0               ; set min contour level
  res@cnMaxLevelValF       = 30.0                ; set max contour level
  res@cnLevelSpacingF      = 3.0                ; set contour spacing
  res@gsnCenterString = "ERA5"
  ;plot(id) = gsn_csm_contour_map_ce(wks,u,res)
  plot(id) = gsn_csm_pres_hgt(wks,u({1000:100},{20},{130:150}),res)
  res@gsnCenterString = "GSM Tl479"
  ;plot(id+1) = gsn_csm_contour_map_ce(wks,u2,res)
  plot(id+1) = gsn_csm_pres_hgt(wks,u2({1000:100},{20},{130:150}),res)
  res@cnMinLevelValF       = -50.0               ; set min contour level
  res@cnMaxLevelValF       = 50.0                ; set max contour level
  res@cnLevelSpacingF      = 5.0                ; set contour spacing
  res@gsnCenterString = "ERA5"
  ;plot(id+2) = gsn_csm_contour_map_ce(wks,v,res)
  plot(id+2) = gsn_csm_pres_hgt(wks,v({1000:100},{20},{130:150}),res)
  res@gsnCenterString = "GSM Tl479"
  ;plot(id+3) = gsn_csm_contour_map_ce(wks,v2,res)
  plot(id+3) = gsn_csm_pres_hgt(wks,v2({1000:100},{20},{130:150}),res)
  res@cnMinLevelValF       = -1.0e-3               ; set min contour level
  res@cnMaxLevelValF       = 1.0e-3                ; set max contour level
  res@cnLevelSpacingF      = 1.0e-4                ; set contour spacing
  res@gsnCenterString = "ERA5"
  ;plot(id+4) = gsn_csm_contour_map_ce(wks,vo,res)
  plot(id+4) = gsn_csm_pres_hgt(wks,vo({1000:100},{20},{130:150}),res)
  res@gsnCenterString = "GSM Tl479"
  ;plot(id+5) = gsn_csm_contour_map_ce(wks,vo2,res)
  plot(id+5) = gsn_csm_pres_hgt(wks,vo2({1000:100},{20},{130:150}),res)
  ;ihr = ihr + 3
  ;if (ihr.ge.24) then
  ;  ihr = ihr - 24
  ;  idy = idy + 1
  ;end if
  ;d = d + 1
  ;if (d.gt.7) then
  ;d = d - 8
  ;end if
  ;delete(pv)
  ;end do
  resP	                   = True
  ;resP@gsnPanelMainString	   = "MERRA-2 slp + divergence at 200hPa"
  resP@gsnPanelMainString	   = date3d(d)
  ;resP@gsnPanelLabelBar     = True
  gsn_panel(wks,plot,(/3,2/),resP)
  delete(plot)
end
