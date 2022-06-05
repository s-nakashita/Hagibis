load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
begin
  dirname = "/Volumes/dandelion/grib/dias/data/jra55/Hist/Daily/anl_isentrop125/201910"
  vname = "pvor"
  tunits = "hours since 2019-10-05 00:00:00"
  time = ispan(0, 210, 6)
  time@units = tunits
  ntime= dimsizes(time)
  level = 320 ;K
  cmap = read_colormap_file("GMT_polar")
  ;iplot = 3
  t_offset = 18;(iplot-1)*12
  wks = gsn_open_wks("png","IPV_"+sprinti("%0.3d",level)+"K");_"+sprinti("%d",iplot))
  plots = new(12,graphic)
;---Create plot
  res                  = True               ; plot mods desired
  res@gsnDraw	      = False
  res@gsnFrame	      = False
  res@gsnAddCyclic	      = False
  res@cnFillOn             = True               ; turn on color for contours
  res@cnLinesOn            = True              ; turn off contour lines
  res@cnLineLabelsOn       = False              ; turn off contour line labels
  res@cnFillPalette        = cmap
  res@lbLabelBarOn         = False

  res@mpFillOn             = False 
  res@mpMinLonF            =  110.               ; select a subregion
  res@mpMaxLonF            =  160.
  res@mpMinLatF            =  10. 
  res@mpMaxLatF            =  60.

  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF       = -5.0               ; set min contour level
  res@cnMaxLevelValF       = 5.0                ; set max contour level
  res@cnLevelSpacingF      = .5                ; set contour spacing
  
  do it=t_offset,min((/t_offset+11,ntime-1/))
  date = cd_string(time(it),"%Y%N%D%H")
  print(date)
  fname = dirname + "/anl_isentrop125_"+vname+"."+date
  a = addfile(fname,"r")
  ;dNames = getfiledimnames(a)
  ;print(dNames)
  ;vNames = getfilevarnames(a)
  ;print(vNames)
  lat   = a->g0_lat_1
  lon   = a->g0_lon_2
  pv    = a->pVOR_GDS0_THEL({level},:,:)
  pv = pv * 1.0e6
  pv@units = "PVU"
  printVarSummary(pv)
  printMinMax(pv, 0)

  res@gsnLeftString = pv@initial_time
  plots(it-t_offset) = gsn_csm_contour_map_ce(wks,pv,res)
  end do
  resP	                   = True
  resP@gsnPanelMainString	 = sprinti("%0.3d",level)+"K"
  resP@gsnPanelLabelBar    = True
  gsn_panel(wks,plots,(/3,4/),resP)
end