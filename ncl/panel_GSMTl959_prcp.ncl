load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
begin
 yyyymmdd="20191012"
 init = "2019100900"
 Exp = "clim"
 iy = toint(str_get_cols(yyyymmdd, 0, 3))
 im = toint(str_get_cols(yyyymmdd, 4, 5))
 id = toint(str_get_cols(yyyymmdd, 6, 7))
 if (Exp .eq. "clim") then
   fname = "/Volumes/dandelion/GSMJob/Jobwk_Tl959L100/fcst_surf_"+init+".nc"
 else if (Exp .eq. "est") then
   fname = "/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_"+Exp+"/fcst_surf_"+init+".nc"
 else if (Exp .eq. "mgdsst") then
   fname = "/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_"+Exp+"/fcst_surf_"+init+".nc"
 end if
 end if
 end if
 f = addfile(fname, "r")
 time = f->time
 tunits = time@units
 
 wks = gsn_open_wks("png", "GSM_"+Exp+"_prcp_3h_"+init+"_"+yyyymmdd)

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
 res@mpLeftCornerLatF  = 30.0
 res@mpRightCornerLatF = 40.0
 res@mpLeftCornerLonF  = 134.0
 res@mpRightCornerLonF = 142.0
 if (yyyymmdd .eq. "20191012") then
  res@mpLeftCornerLatF = 33.0
  res@mpRightCornerLatF = 40.0
  res@mpLeftCornerLonF = 135.0
  res@mpRightCornerLonF = 143.0
 end if
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
 do t = 0, 3
 t0h = cd_inv_calendar( iy, im, id, 3*t, 0, 0, tunits, 0)
 if (t .lt. 7) then
  t3h = cd_inv_calendar( iy, im, id, 3*(t+1), 0, 0, tunits, 0)
 else
  t3h = cd_inv_calendar( iy, im, id+1, 0, 0, 0, tunits, 0)
 end if
 prcp1 = f->RAIN({t0h}, :, :)
 prcp2 = f->RAIN({t3h}, :, :)
 prcp = prcp2
 prcp = prcp - prcp1
 printMinMax(prcp, 0)

 res@gsnCenterString = sprinti("%0.2d", 3*t+t_offset)+"-"+sprinti("%0.2d", 3*(t+1)+t_offset)+" UTC"
 if (t .eq. 4) then
   res@lbLabelBarOn = True
 else
   res@lbLabelBarOn = False
 end if
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
 resP@gsnPanelBottom = 0.08 
 resP@gsnPanelSave   = True      ; Save the state of the paneled plots so we can
                                    ; query for text strings locations.
 resP@gsnPanelMainString = yyyymmdd+" "+Exp

 gsn_panel(wks, plot, (/2,4/), resP)  ; plots
 frame(wks)
end