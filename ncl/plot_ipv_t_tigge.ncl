load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_inv_string.ncl"
begin
  dirname = "/Volumes/dandelion/grib/tigge"
  orig = "egrr" ; ecmf or kwbc or egrr
  fname = dirname + "/tigge_glb_"+orig+"_pt_2019-10-09_12_cntl.grb"
  a = addfile(fname,"r")
  dNames = getfiledimnames(a)
  print(dNames)
  vNames = getfilevarnames(a)
  print(vNames)
  lat   = a->lat_0
  lon   = a->lon_0
  time  = a->forecast_time0
  pv    = a->pv_P1_L107_GLL0(:,:,:)
  pv = pv * 1.0e6
  pv@units = "PVU"
  printVarSummary(pv)
  printMinMax(pv, 0)
  format = "%N/%D/%Y (%H:%M)"
  init = cd_inv_string(pv@initial_time, format)
  print(init)
  ntime= dimsizes(time)
  level = 320 ;K
  cmap = read_colormap_file("GMT_polar")
  iplot = 1
  t_offset = (iplot-1)*12
  wks = gsn_open_wks("png","IPV_"+orig+"_"+sprinti("%0.3d",level)+"K_"+sprinti("%d",iplot))
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
  res@gsnLeftString = "FT="+sprinti("%d",time(it))+" H"
  plots(it-t_offset) = gsn_csm_contour_map_ce(wks,pv(it,:,:),res)
  end do
  resP	                   = True
  resP@gsnPanelMainString	 = sprinti("%0.3d",level)+"K "+pv@initial_time
  resP@gsnPanelLabelBar    = True
  gsn_panel(wks,plots,(/3,4/),resP)
end