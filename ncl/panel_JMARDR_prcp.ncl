begin
 yyyymmdd="20191012"
 yyyy = str_get_cols(yyyymmdd, 0, 3)
 mm = str_get_cols(yyyymmdd, 4, 5)
 dd = str_get_cols(yyyymmdd, 6, 7)
 dir = "/Volumes/dandelion/netcdf/JMAGPV/hourly/"+yyyy+"/"+mm+"/"+dd

 wks = gsn_open_wks("png", "JMARDR_prcp_3h_"+yyyymmdd)

 res = True
 res@gsnDraw = False
 res@gsnFrame = False

 res@cnFillOn = True
 res@cnLinesOn = False
 res@cnLineLabelsOn = False
 res@cnFillMode = "RasterFill"
 res@cnRasterSmoothingOn = True

 res@lbTitleFontHeightF = .03
 res@lbTitleString = "mm/3h"
 res@lbTitlePosition            = "Bottom"            ; title position
 res@lbTitleDirection           = "Across"            ; title direction
 res@lbLabelStride              = 2
 res@lbLabelFontHeightF         = 0.025
 res@pmLabelBarHeightF          = 0.1
 res@lbBoxEndCapStyle           = "TriangleBothEnds"
 res@lbOrientation              = "Horizontal"        ; orientation

 res@mpProjection = "CylindricalEquidistant"
 res@mpLimitMode  = "Corners"
 res@mpLeftCornerLatF = 33.0
 res@mpRightCornerLatF = 40.0
 res@mpLeftCornerLonF = 135.0
 res@mpRightCornerLonF = 143.0
 res@mpFillDrawOrder = "PostDraw"
 res@mpOutlineOn = True
 res@mpFillOn = False

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

 res@cnFillPalette             = "MPL_YlGnBu"
 res@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
 res@cnMinLevelValF            = 10             ; min level
 res@cnMaxLevelValF            = 90             ; max level
 res@cnLevelSpacingF           = 10.             ; contour interval

 plot = new(8, "graphic")
 deg2rad = acos(-1.) / 180.
 t_offset = 0
 do t = 0, 7
 fname = dir + "/jrap"+yyyymmdd+sprinti("%0.2d", 3*t+t_offset)+".nc"
 f = addfile(fname, "r")
 lon_in = f->longitude
 lon_out = lon_in(5::10)
 lat_in = f->latitude
 lat_out = lat_in(5::10)
 prcp_in = f->var0_1_201_surface(0,:,:)
 prcp1 = area_hi2lores_Wrap(lon_in, lat_in, prcp_in, False, cos(lat_in * deg2rad), lon_out, lat_out, False)
 prcp = prcp1
 fname = dir + "/jrap"+yyyymmdd+sprinti("%0.2d", 3*t+1+t_offset)+".nc"
 f = addfile(fname, "r")
 lon_in = f->longitude
 lon_out = lon_in(5::10)
 lat_in = f->latitude
 lat_out = lat_in(5::10)
 prcp_in = f->var0_1_201_surface(0,:,:)
 prcp2 = area_hi2lores_Wrap(lon_in, lat_in, prcp_in, False, cos(lat_in * deg2rad), lon_out, lat_out, False)
 prcp = prcp + prcp2
 fname = dir + "/jrap"+yyyymmdd+sprinti("%0.2d", 3*t+2+t_offset)+".nc"
 f = addfile(fname, "r")
 lon_in = f->longitude
 lon_out = lon_in(5::10)
 lat_in = f->latitude
 lat_out = lat_in(5::10)
 prcp_in = f->var0_1_201_surface(0,:,:)
 prcp3 = area_hi2lores_Wrap(lon_in, lat_in, prcp_in, False, cos(lat_in * deg2rad), lon_out, lat_out, False)
 prcp = prcp + prcp3
 printMinMax(prcp, 0)

 res@gsnCenterString = sprinti("%0.2d", 3*t+t_offset)+"-"+sprinti("%0.2d", 3*(t+1)+t_offset)+" UTC"
 if (t .eq. 0).or.(t .eq. 4) then
   res@lbLabelBarOn = True
 else
   res@lbLabelBarOn = False
 end if
 ;if (t .eq. 7) then
 if (t .eq. 0) then
    res@tmYRLabelsOn               = True            ; turn on lat labels
    res@tmXBLabelsOn               = True
    res@tmYROn                     = True
    res@tmXBOn                     = True
    res@tmYRLabelFontHeightF       = 0.02
    res@tmXBLabelFontHeightF       = 0.02
 end if
 plot(t) = gsn_csm_contour_map(wks, prcp, res)
 end do

 ; polygon
 ; Tokai
 x = (/136., 138., 138., 136., 136./)
 y = (/ 36.,  36.,  34.,  34.,  36./)
 pgres = True
 pgres@gsFillColor = "lightcoral"
 pgres@gsFillOpacityF = 0.3
 gon1 = gsn_add_polygon(wks, plot(0), x, y, pgres)
 ; Kanto
 x = (/ 138.,  141.,  141.,  138.,  138./)
 y = (/ 36.5,  36.5,  34.5,  34.5,  36.5/)
 pgres = True
 pgres@gsFillColor = "lightcoral"
 pgres@gsFillOpacityF = 0.3
 gon2 = gsn_add_polygon(wks, plot(1), x, y, pgres)
 ; Southern Tohoku
 x = (/ 139., 141.5, 141.5,  139.,  139./)
 y = (/ 38.5,  38.5,  36.5,  36.5,  38.5/)
 pgres = True
 pgres@gsFillColor = "lightcoral"
 pgres@gsFillOpacityF = 0.3
 gon3 = gsn_add_polygon(wks, plot(2), x, y, pgres)
 ; Northern Tohoku
 x = (/140.5, 142.5, 142.5, 140.5, 140.5/)
 y = (/ 40.5,  40.5,  38.5,  38.5,  40.5/)
 pgres = True
 pgres@gsFillColor = "lightcoral"
 pgres@gsFillOpacityF = 0.3
 gon4 = gsn_add_polygon(wks, plot(5), x, y, pgres)

 resP                = True
 resP@gsnMaximize    = True
 resP@gsnFrame       = False
 resP@gsnPanelLeft   = 0.005
 resP@gsnPanelRight  = 0.965
 resP@gsnPanelBottom = 0.08 
 resP@gsnPanelSave   = True      ; Save the state of the paneled plots so we can
                                    ; query for text strings locations.
 resP@gsnPanelMainString = yyyymmdd

 gsn_panel(wks, plot, (/2,4/), resP)  ; plots
 frame(wks)
end