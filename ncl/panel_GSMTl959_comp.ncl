load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
begin
 yyyymmddhh="2019101100"
 yyyymmdd = str_get_cols(yyyymmddhh, 0, 7)
 yyyy = str_get_cols(yyyymmddhh, 0, 3)
 mm = str_get_cols(yyyymmddhh, 4, 5)
 dd = str_get_cols(yyyymmddhh, 6, 7)
 rdrdir = "/Volumes/dandelion/netcdf/JMAGPV/hourly/"+yyyy+"/"+mm
 init = "2019100900"
 iy = toint(str_get_cols(yyyymmddhh, 0, 3))
 im = toint(str_get_cols(yyyymmddhh, 4, 5))
 id = toint(str_get_cols(yyyymmddhh, 6, 7))
 ih = toint(str_get_cols(yyyymmddhh, 8, 9))
 fclim = "/Volumes/dandelion/GSMJob/Jobwk_Tl959L100/fcst_surf_asia_"+init+".nc"
 fest  = "/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_est/fcst_surf_asia_"+init+".nc"
 fmgd  = "/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_mgdsst/fcst_surf_asia_"+init+".nc"
 files = (/fclim, fest, fmgd/)
 sst   = (/"CLIM","EST","MGD"/)

 fgsm = addfile(fclim, "r")
 time = fgsm->time
 tunits = time@units
 delete(fgsm)

 int_h = 36
 wks = gsn_open_wks("pdf", "comp_prcp_"+sprinti("%d", int_h)+"h_"+init+"_"+yyyymmddhh)
 plot_region = True

 res = True
 res@gsnDraw = False
 res@gsnFrame = False

 res@cnFillOn = True
 res@cnLinesOn = False
 res@cnLineLabelsOn = False
 
 res@lbTitleFontHeightF = .03
 res@lbTitleString = "mm"
 res@lbTitlePosition            = "Bottom"            ; title position
 res@lbTitleDirection           = "Across"            ; title direction
 res@lbLabelStride              = 1
 res@lbLabelFontHeightF         = 0.015
 ;res@pmLabelBarHeightF          = 0.1
 ;res@lbBoxEndCapStyle           = "TriangleBothEnds"
 res@lbOrientation              = "Horizontal"        ; orientation

 res@mpProjection = "CylindricalEquidistant"
 res@mpLimitMode  = "Corners"
 res@mpLeftCornerLatF  = 23.5
 res@mpRightCornerLatF = 38.5
 res@mpLeftCornerLonF  = 127.5
 res@mpRightCornerLonF = 141.5
 ;if (yyyymmdd .eq. "20191012") then
  res@mpLeftCornerLatF = 33.0
  res@mpRightCornerLatF = 40.0
  res@mpLeftCornerLonF = 135.0
  res@mpRightCornerLonF = 143.0
 ;end if
 res@mpFillDrawOrder = "PostDraw"
 res@mpOutlineOn = True
 res@mpFillOn = False
 res@mpGeophysicalLineColor = "white"
 res@mpGeophysicalLineThicknessF = 3.0
 res@mpDataBaseVersion = "HighRes"

 res@gsnAddCyclic               = False            ; regional data, don't add cyclic point
 
 ; Set resources for tickmarks
 res@tmYLLabelsOn               = False            ; turn off lat labels
 res@tmXBLabelsOn               = False            ; turn off lon labels
 
 res@tmYROn                     = False
 res@tmYLOn                     = False
 res@tmXTOn                     = False
 res@tmXBOn                     = False

 res@gsnStringFont              = "helvetica-bold"  ; will be used for gsnCenterString later
 res@gsnCenterStringFontHeightF = 0.04              ; ditto
 res@gsnLeftString              = ""                ; remove the left and right substrings
 res@gsnRightString             = ""

 cmap_default = read_colormap_file("BkBlAqGrYeOrReViWh200")
 cmap = new((/10,4/), float)
 cmap(0,:) = cmap_default(16,:)
 cmap(1,:) = cmap_default(32,:)
 cmap(2,:) = cmap_default(48,:)
 cmap(3,:) = cmap_default(64,:)
 cmap(4,:) = cmap_default(96,:)
 cmap(5,:) = cmap_default(100,:)
 cmap(6,:) = cmap_default(116,:)
 cmap(7,:) = cmap_default(132,:)
 cmap(8,:) = cmap_default(148,:)
 cmap(9,:) = cmap_default(164,:)
 res@cnFillPalette             = cmap;_default(16:160,:)
 res@cnLevelSelectionMode      = "ExplicitLevels"  ; manually specify contour levels
 res@cnLevels = (/.5, 1., 2., 3., 5., 10., 15., 20., 25./)
 res@cnLevels = res@cnLevels * 18.0

 plot = new(4, "graphic")
 deg2rad = acos(-1.) / 180.
 ;Radar
 fname = rdrdir+"/"+dd + "/jrap"+yyyymmdd+sprinti("%0.2d", ih)+".nc"
 f = addfile(fname, "r")
 lon_in = f->longitude
 lon_out = lon_in(5::10)
 lat_in = f->latitude
 lat_out = lat_in(5::10)
 prcp_in = f->var0_1_201_surface(0,:,:)
 prcp = area_hi2lores_Wrap(lon_in, lat_in, prcp_in, False, cos(lat_in * deg2rad), lon_out, lat_out, False)
 delete(f)
 if (int_h .gt. 1) then
 do t = 1, int_h-1
  if (ih+t.lt.24) then
    fname = rdrdir+"/"+dd + "/jrap"+yyyymmdd+sprinti("%0.2d", ih+t)+".nc"
  else
    dd2 = sprinti("%0.2d", id+1)
    hh2 = sprinti("%0.2d", ih+t-24)
    fname = rdrdir+"/"+dd2 + "/jrap"+yyyy+mm+dd2+hh2+".nc"
  end if
  f = addfile(fname, "r")
  prcp_in = f->var0_1_201_surface(0,:,:)
  prcp1 = area_hi2lores_Wrap(lon_in, lat_in, prcp_in, False, cos(lat_in * deg2rad), lon_out, lat_out, False)
  prcp = prcp + prcp1
  delete(prcp1)
 end do 
 delete(f)
 end if
 printMinMax(prcp, 0)
 res@gsnCenterString = "Radar"
 res@lbLabelBarOn               = False
 res@tmYRLabelsOn               = True            ; turn on lat labels
 res@tmXBLabelsOn               = True
 res@tmYROn                     = True
 res@tmXBOn                     = True
 res@tmYRLabelFontHeightF       = 0.02
 res@tmXBLabelFontHeightF       = 0.02
 res@cnFillMode = "RasterFill"
 res@cnRasterSmoothingOn = True
 plot(0) = gsn_csm_contour_map(wks, prcp, res)
 ; polyline
 if (plot_region) then
 ; Tokai
 x = (/136., 138., 138., 136., 136./)
 y = (/ 36.,  36.,  34.,  34.,  36./)
 pgres = True
 pgres@gsLineColor = "Black"
 pgres@gsLineThicknessF = 2.0
 dum = new(4*3, graphic)
 i = 0
 do j = 0, 3
   dum(i) = gsn_add_polyline(wks, plot(0), x(j:j+1), y(j:j+1), pgres)
   i = i + 1
 end do
 ; Kanto
 x = (/ 138.,  141.,  141.,  138.,  138./)
 y = (/ 36.5,  36.5,  34.5,  34.5,  36.5/)
 do j = 0, 3
   dum(i) = gsn_add_polyline(wks, plot(0), x(j:j+1), y(j:j+1), pgres)
   i = i + 1
 end do
 ; Southern Tohoku
 x = (/ 139., 141.5, 141.5,  139.,  139./)
 y = (/ 38.5,  38.5,  36.5,  36.5,  38.5/)
 do j = 0, 3
   dum(i) = gsn_add_polyline(wks, plot(0), x(j:j+1), y(j:j+1), pgres)
   i = i + 1
 end do
 end if
 delete(prcp)
 delete(res@cnFillMode)
 res@cnRasterSmoothingOn = False
 ; GSM GPV
 t0h = cd_inv_calendar( iy, im, id, ih, 0, 0, tunits, 0)
 ih_pos = ih + int_h
 if (ih_pos .lt. 24) then
  tpos = cd_inv_calendar( iy, im, id, ih_pos, 0, 0, tunits, 0)
 else
  tpos = cd_inv_calendar( iy, im, id+1, ih_pos-24, 0, 0, tunits, 0)
 end if
 do j = 0, 2
 f = addfile(files(j), "r")
 prcp1 = f->RAIN({t0h}, :, :)
 prcp2 = f->RAIN({tpos}, :, :)
 prcp = prcp2
 prcp = prcp - prcp1
 printMinMax(prcp, 0)

 res@gsnCenterString = sst(j)
 res@lbLabelBarOn = False
 plot(j+1) = gsn_csm_contour_map(wks, prcp, res)
 delete(prcp1)
 delete(prcp2)
 delete(prcp)
 end do
 pmres = True
 pmres@gsMarkerIndex = 16
 pmres@gsMarkerSizeF = 10.0
 pmres@gsMarkerThicknessF = 3.0
 ;dum = gsn_add_polymarker(wks, plot(3), 139.0, 35.0, pmres)

 resP                = True
 resP@gsnMaximize    = True
 resP@gsnFrame       = False
 resP@gsnPanelLeft   = 0.005
 resP@gsnPanelRight  = 0.965
 resP@gsnPanelBottom = 0.1 
 resP@gsnPanelLabelBar   = True
 resP@lbTitleFontHeightF = .015
 resP@lbTitleString      = "mm"
 resP@lbTitlePosition    = "Bottom"            ; title position
 resP@lbTitleDirection   = "Across"            ; title direction
 resP@lbLabelStride      = 1
 resP@lbLabelFontHeightF = 0.01
 resP@gsnPanelSave   = True      ; Save the state of the paneled plots so we can
                                    ; query for text strings locations.
 ;resP@gsnPanelMainString = yyyy+"-"+mm+"-"+dd + " " + sprinti("%0.2d", ih)+"-"+sprinti("%0.2d", ih_pos)+" UTC"
 if (ih_pos .lt. 24) then
 resP@gsnPanelMainString = sprinti("%0.2d", ih)+"00UTC"+ " - "+ sprinti("%0.2d", ih_pos)+"00UTC"+dd
 else
 resP@gsnPanelMainString = sprinti("%0.2d", ih)+"00UTC"+dd + " " + " - "+ sprinti("%0.2d", ih_pos-24)+"00UTC"+dd2
 end if
 gsn_panel(wks, plot, (/1,4/), resP)  ; plots
 frame(wks)
end