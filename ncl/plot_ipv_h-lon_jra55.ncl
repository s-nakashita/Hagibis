load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
begin
  dirname = "/Volumes/dandelion/grib/dias/data/jra55/Hist/Daily/anl_isentrop125/201910"
  vname = "pvor"
  date  = "2019100900"
  fname = dirname + "/anl_isentrop125_"+vname+"."+date
  a = addfile(fname,"r")
  vname = "pres"
  fname = dirname + "/anl_isentrop125_"+vname+"."+date
  b = addfile(fname,"r")
  ;dNames = getfiledimnames(a)
  ;print(dNames)
  vNames = getfilevarnames(b)
  print(vNames)
  level = a->lv_THEL0
  lat   = a->g0_lat_1
  lon   = a->g0_lon_2
  pv    = a->pVOR_GDS0_THEL
  pres  = b->PRES_GDS0_THEL
  nlev = dimsizes(level)
  print(level)
  pv = pv * 1.0e6
  pv@units = "PVU"
  printVarSummary(pv)
  printMinMax(pv, 0)
  pres = pres * 1.0e-2
  pres@units = "hPa"
  printVarSummary(pres)
  printMinMax(pres, 0)

  cmap = read_colormap_file("GMT_polar")
  latc = 20
  wks = gsn_open_wks("png","IPV_h-lon_"+date+"_lat"+sprinti("%d",latc))
;---Create plot
  res                  = True               ; plot mods desired
  res@gsnDraw	      = False
  res@gsnFrame	      = False
  res@gsnAddCyclic	      = False
  resp = res
  res@cnFillOn             = True               ; turn on color for contours
  res@cnLinesOn            = False              ; turn off contour lines
  res@cnLineLabelsOn       = False              ; turn off contour line labels
  res@cnFillPalette        = cmap

  res@gsnCenterStringParallelPosF = 0.6
  res@gsnCenterString = pv@initial_time
; note, when doing a subregion, NCL determines the range of the data from
; the full domain. If you wish to just consider the domain you are plotting,
; you must manually set those levels.

  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF       = -10.0               ; set min contour level
  res@cnMaxLevelValF       = 10.0                ; set max contour level
  res@cnLevelSpacingF      = 1.                ; set contour spacing
  plot = gsn_csm_contour(wks, pv({300:400},{latc},{110:160}), res)

  resp@gsnCenterString = ""
  resp@gsnRightString = ""
  resp@gsnLeftString = ""
  resp@cnFillOn             = False             ; turn on color for contours
  resp@cnLinesOn            = True              ; turn off contour lines
  resp@cnLineLabelsOn       = True              ; turn off contour line labels
  resp@cnLevelSelectionMode = "ManualLevels"
  resp@cnMinLevelValF       = 100.0               ; set min contour level
  resp@cnMaxLevelValF       = 1000.0              ; set max contour level
  resp@cnLevelSpacingF      = 100.0               ; set contour spacing
  plotp = gsn_csm_contour(wks, pres({300:400},{latc},{110:160}), resp)
  overlay(plot,plotp)
  draw(plot)
  frame(wks)
end