load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
begin
  dirname = "/Volumes/dandelion/grib/dias/data/jra55/Hist/Daily/anl_isentrop125/201910"
  vname = "pvor"
  date  = "2019100900"
  fname = dirname + "/anl_isentrop125_"+vname+"."+date
  a = addfile(fname,"r")
  ;dNames = getfiledimnames(a)
  ;print(dNames)
  ;vNames = getfilevarnames(a)
  ;print(vNames)
  level = a->lv_THEL0
  lat   = a->g0_lat_1
  lon   = a->g0_lon_2
  pv    = a->pVOR_GDS0_THEL
  nlev = dimsizes(level)
  print(level)
  pv = pv * 1.0e6
  pv@units = "PVU"
  printVarSummary(pv)
  printMinMax(pv, 0)

  cmap = read_colormap_file("GMT_polar")
  iplot = 2
  k_offset = (iplot-1)*12
  wks = gsn_open_wks("png","IPV_"+date+"_"+sprinti("%d",iplot))
  plots = new(12,graphic)
;---Create plot
  res                  = True               ; plot mods desired
  res@gsnDraw	      = False
  res@gsnFrame	      = False
  res@gsnAddCyclic	      = False
  res@cnFillOn             = True               ; turn on color for contours
  res@cnLinesOn            = True              ; turn off contour lines
  res@cnLineLabelsOn       = False              ; turn off contour line labels
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnFillPalette        = cmap

  res@gsnCenterStringParallelPosF = 0.6

  res@mpFillOn             = False 
  res@mpMinLonF            =  110.               ; select a subregion
  res@mpMaxLonF            =  160.
  res@mpMinLatF            =  10. 
  res@mpMaxLatF            =  60.

; note, when doing a subregion, NCL determines the range of the data from
; the full domain. If you wish to just consider the domain you are plotting,
; you must manually set those levels.

  do k=k_offset, min((/k_offset+11,nlev-1/))
  if (k.lt.8) then
  res@cnMinLevelValF       = -5.0               ; set min contour level
  res@cnMaxLevelValF       = 5.0                ; set max contour level
  res@cnLevelSpacingF      = .5                ; set contour spacing
  else if (k.lt.12) then
  res@cnMinLevelValF       = -10.0               ; set min contour level
  res@cnMaxLevelValF       = 10.0                ; set max contour level
  res@cnLevelSpacingF      = 1.                ; set contour spacing
  else if (k.lt.15) then
  res@cnMinLevelValF       = -20.0               ; set min contour level
  res@cnMaxLevelValF       = 20.0                ; set max contour level
  res@cnLevelSpacingF      = 2.                ; set contour spacing
  else if (k.lt.18) then
  res@cnMinLevelValF       = -50.0               ; set min contour level
  res@cnMaxLevelValF       = 50.0                ; set max contour level
  res@cnLevelSpacingF      = 5.                ; set contour spacing
  else
  res@cnMinLevelValF       = -200.0               ; set min contour level
  res@cnMaxLevelValF       = 200.0                ; set max contour level
  res@cnLevelSpacingF      = 20.                ; set contour spacing
  end if
  end if
  end if
  end if
  res@gsnCenterString = sprinti("%0.3d",level(k))+"K"
  plots(k-k_offset) = gsn_csm_contour_map_ce(wks,pv(k,:,:),res)
  end do
  resP	                   = True
  resP@gsnPanelMainString	   = pv@initial_time
  gsn_panel(wks,plots,(/3,4/),resP)
end