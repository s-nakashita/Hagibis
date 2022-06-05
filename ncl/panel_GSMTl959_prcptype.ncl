load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
begin
 ;yyyymmdd="20191012"
 init = "2019100912"
 Exp = "clim"
 ;iy = toint(str_get_cols(yyyymmdd, 0, 3))
 ;im = toint(str_get_cols(yyyymmdd, 4, 5))
 ;id = toint(str_get_cols(yyyymmdd, 6, 7))
 if (Exp .eq. "clim") then
   fname = "/Volumes/dandelion/GSMJob/Jobwk_Tl959L100/fcst_phy2m_avr_"+init+".nc"
 else if (Exp .eq. "est") then
   fname = "/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_"+Exp+"/fcst_phy2m_avr_"+init+".nc"
 else if (Exp .eq. "mgdsst") then
   fname = "/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_"+Exp+"/fcst_phy2m_avr_"+init+".nc"
 end if
 end if
 end if
 f = addfile(fname, "r")
 time = f->time
 ntime = dimsizes(time)
 tunits = time@units
 date = cd_string(time, "%N%D %H:%M")
 
 ds = ind(date.eq."1012 06:00")
 de = min((/ds+4-1, ntime-1/))
 wks = gsn_open_wks("png", \
 "GSM_"+Exp+"_prcptype_6h_"+init+"_"+cd_string(time(ds), "%D%H")+"-"+cd_string(time(de), "%D%H"))

 res = True
 res@gsnDraw = False
 res@gsnFrame = False

 res@cnFillOn = True
 res@cnLinesOn = False
 res@cnLineLabelsOn = False
 ;res@cnFillMode = "RasterFill"
 ;res@cnRasterSmoothingOn = True

 res@lbTitleFontHeightF = .02
 res@lbTitleString = "mm/hr"
 res@lbTitlePosition            = "Right"            ; title position
 res@lbTitleDirection           = "Across"            ; title direction
 res@lbLabelStride              = 2
 res@lbLabelFontHeightF         = 0.02
 ;res@lbBoxEndCapStyle           = "TriangleBothEnds"
 res@lbOrientation              = "Horizontal"        ; orientation

 res@mpProjection = "CylindricalEquidistant"
 res@mpLimitMode  = "Corners"
 ;res@mpLeftCornerLatF  = 23.5
 ;res@mpRightCornerLatF = 38.5
 ;res@mpLeftCornerLonF  = 127.5
 ;res@mpRightCornerLonF = 141.5
 ;if (ds .eq. 12) then
  res@mpLeftCornerLatF = 33.0
  res@mpRightCornerLatF = 40.0
  res@mpLeftCornerLonF = 137.0
  res@mpRightCornerLonF = 143.0
 ;end if
 res@mpFillDrawOrder = "PostDraw"
 res@mpDataBaseVersion = "HighRes"
 res@mpOutlineOn = True
 res@mpFillOn = False

 res@gsnAddCyclic               = True            ; regional data, don't add cyclic point
 
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
 ;res@cnMinLevelValF            = 20.             ; min level
 ;res@cnMaxLevelValF            = 180.            ; max level
 ;res@cnLevelSpacingF           = 20.             ; contour interval

 plot = new(4*3, "graphic")
 deg2rad = acos(-1.) / 180.
 t_offset = 0
 do d = ds, de
 prcp2 = f->RAIN(d, :, :)
 prcp = prcp2
 if (d.gt.1) then
    prcp1 = f->RAIN(d-1, :, :)
    prcp = prcp - prcp1
 end if
 prcp = prcp / 6. ; per 6 hours -> per 1 hour
 printMinMax(prcp, 0)
 prc = f->aPRCd(d,:,:)
 prc = prc / 24. ; per day -> per 1 hour
 prl = f->aPRLd(d,:,:)
 prl = prl / 24. ; perday -> per 1 hour

 ;res@gsnLeftString = "mm/hr"
 if (d .eq. de) then
   res@lbLabelBarOn = True
 else
   res@lbLabelBarOn = False
 end if
 if (d .eq. ds) then
    res@tmYRLabelsOn               = True            ; turn on lat labels
    res@tmXBLabelsOn               = True
    res@tmYROn                     = True
    res@tmXBOn                     = True
    res@tmYRLabelFontHeightF       = 0.02
    res@tmXBLabelFontHeightF       = 0.02
 end if
 ;if (d .eq. ds) then
 ;   res@gsnCenterString = "total"
 ;end if
 res@gsnCenterString = date(d-1) + " - " + date(d)
 plot(3*(d-ds))   = gsn_csm_contour_map(wks, prcp, res)
 ;if (d .eq. ds) then
    res@gsnCenterString = "convective"
 ;else
 ;   res@gsnCenterString = ""
 ;end if
 res@gsnRightString = ""
 plot(3*(d-ds)+1) = gsn_csm_contour_map(wks, prc, res)
 ;if (d .eq. ds) then
    res@gsnCenterString = "stratiform"
 ;else
 ;   res@gsnCenterString = ""
 ;end if
 res@gsnRightString = ""
 ;res@lbLabelBarOn = True
 plot(3*(d-ds)+2) = gsn_csm_contour_map(wks, prl, res)
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
 ;resP@gsnPanelMainString = yyyymmdd+" "+Exp

 gsn_panel(wks, plot, (/4,3/), resP)  ; plots
 frame(wks)
end