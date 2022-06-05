;**********************************************************************
; This script creates a panel of 8 plots (4 rows and 4 columns):
;
;   - the 1st column of plots uses one color map and is represented by a single labelbar
;   - the 2nd and 3rd column of plots uses one color map and is represented by a single labelbar
;   - the 4th column of plots uses one color map and is represented by a single labelbar
;
; This script was written by Appolinaire Derbetini of the Laboratory 
; for Environmental Modelling and Atmospheric Physics, University of 
; Yaounde 1, Yaounde, Cameroon.
; 
; The idea was inspired by the paper of Pokam et al., 2018, Consequences
; of 1.5 C and 2 C global warming levels for temperature and
; precipitation changes over Central Africa, Environ. Res. Lett. 13
; (2018) 055011, https://doi.org/10.1088/1748-9326/aab048
;**********************************************************************

;**********************************************************************
; Simple function that takes the difference of two variables and
; copies the coordinate arrays from one of the variables to the
; difference variable for plotting later.
;**********************************************************************
function diff_vars(x1,x2)
local xdiff
begin
  xdiff = x1-x2
  copy_VarCoords(x1,xdiff)
  return(xdiff)
end

;**********************************************************************
; Function to retrieve the location of a plot so we can calculate
; a location for a title on the left side of the plot.
;**********************************************************************
function get_title_position(plot)
local vpy, vpx, vph
begin
  getvalues plot
      "vpYF"       : vpy
      "vpXF"       : vpx
      "vpHeightF"  : vph
  end getvalues
  pos = (/vpx*0.75,vpy-(vph/2.)/)
  return(pos)
end

;**********************************************************************
; Main code
;**********************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
begin 
   ;Reading netcdf data
   yyyymmddhh = 2019100900
   fhp1    = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_est/fcst_surf_asia_"+yyyymmddhh+".nc", "r")
   fh1     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_est/fcst_phy2m_snp_"+yyyymmddhh+".nc", "r")
   fhp2    = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100/fcst_surf_asia_"+yyyymmddhh+".nc", "r")
   fh2     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100/fcst_phy2m_snp_"+yyyymmddhh+".nc","r")
   fhp3    = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_mgdsst/fcst_surf_asia_"+yyyymmddhh+".nc", "r")
   fh3     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_mgdsst/fcst_phy2m_snp_"+yyyymmddhh+".nc", "r")
   
   pr_1 = fhp1->RAIN(:, :, :)
   lh_1 = fh1->TPW(:, :, :)
   u    = fh1->U(:, :, :)
   v    = fh1->V(:, :, :)
   pr_2 = fhp2->RAIN(:, :, :)
   lh_2 = fh2->TPW(:, :, :)
   pr_3 = fhp3->RAIN(:, :, :)
   lh_3 = fh3->TPW(:, :, :)
   t1h  = fhp1->time
   t6h  = fh1->time

   date = cd_string(t6h, "%N%D %H:%M")
   ntim = dimsizes(date)
   
   ;Reading track data
   data = readAsciiTable("../pytrack/track"+yyyymmddhh+"_gsm_tl959_est.txt", 7, "float", 0)
   npoints = dimsizes(data(:,0))
   dd = toint(data(:,2))
   hh = toint(data(:,3))
   lon1 = data(:,4)
   lat1 = data(:,5)
   delete(data)
   data = readAsciiTable("../pytrack/track"+yyyymmddhh+"_gsm_tl959.txt", 7, "float", 0)
   lon2 = data(:,4)
   lat2 = data(:,5)
   delete(data)
   data = readAsciiTable("../pytrack/track"+yyyymmddhh+"_gsm_tl959_mgdsst.txt", 7, "float", 0)
   lon3 = data(:,4)
   lat3 = data(:,5)
   
   pltDir = "./"                ; where output will be written
   pltName= "GSMtl959diff_tpw+rain_"+yyyymmddhh                                   
   pltType= "png"               ; ps, eps, png, pdf, x11 .....

;****************************************************************************
; Graphic resources
;****************************************************************************
    ;ver = 4
  nrow = 3
  ncol = 3
;  do ds = 1, 13, 4
  plot = new(nrow*ncol , "graphic")
  dum = new(nrow*ncol, graphic)
  ds = 10
  de = min((/ds+nrow-1,ntim-1/))
  i = 0
  pltPath = pltDir+pltName+"_"+cd_string(t6h(ds), "%D%H")+"-"+cd_string(t6h(de), "%D%H")
  wks     = gsn_open_wks(pltType, pltPath) 
  do d=ds, de
    res1                            = True     ; plot mods desired
    res1@gsnDraw                    = False    ; will panel plots later
    res1@gsnFrame                   = False

    res1@cnFillOn                   = True     ; turn on color fill
    res1@cnLinesOn                  = False    ; turn off contour lines
    res1@cnLineLabelsOn             = False    ; turn off contour line labels

;**********************************************************************
; We are turning the labelbar off for all but the bottom row of plots
;**********************************************************************
    res1@lbLabelBarOn               = False
    res1@lbTitleFontHeightF         = .03                 ; make title smaller
    res1@lbTitleString              = "kg/m^2"            ; title string
    res1@lbTitlePosition            = "Bottom"            ; title position
    res1@lbTitleDirection           = "Across"            ; title direction
    ;res1@lbLabelStride              = 2
    res1@lbLabelFontHeightF         = 0.025
    res1@pmLabelBarHeightF          = 0.1
    res1@lbBoxEndCapStyle           = "TriangleBothEnds"
    ; Set resources for map
    if(ds.lt.9) then
      res1@mpMinLatF                  = 15.
      res1@mpMaxLatF                  = 35.
      res1@mpMinLonF                  = 125.
    else
      res1@mpMinLatF                  = 25.
      res1@mpMaxLatF                  = 40.
      res1@mpMinLonF                  = 130.
    end if
    res1@mpMaxLonF                  = 150.

    res1@mpDataBaseVersion          = "MediumRes"
    res1@mpOutlineOn                = True
    res1@mpFillOn                   = False                ; turn off gray fill
    ;res@mpOutlineBoundarySets      = "National"

    res1@gsnAddCyclic               = False            ; regional data, don't add cyclic point
 
    ; Set resources for tickmarks
    res1@gsnStringFont              = "helvetica-bold"  ; will be used for gsnCenterString later
    res1@gsnCenterStringFontHeightF = 0.04              ; ditto
    res1@gsnLeftString              = ""                ; remove the left and right substrings
    res1@gsnRightString             = date(d)
    
   ; Set resources specifically for plots in column 1
    res1@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
    
    res2 = True
    res2@gsnDraw = False
    res2@gsnFrame = False
    res2@gsnLeftString              = ""
    res2@gsnRightString             = ""
    res2@cnLinesOn                  = True    ; turn off contour lines
    res2@cnLineLabelsOn             = True
    res2@cnLevelSelectionMode       = "ExplicitLevels"  ; manually specify contour levels
    res2@gsnContourNegLineDashPattern = 2
    res2@cnLineThicknessF           = 3.0

   ; For track
    pmres = True
    pmres@gsMarkerIndex = 16
    pmres@gsMarkerSizeF = 5.0
    pmres@gsMarkerColor = "red"
;****************************************************************************
; Create the plots for the panel
;****************************************************************************
    ; first column : raw values
    ind_p = min((/ind(t1h.eq.t6h(d))+3, dimsizes(t1h)-1/)) ; 3H after
    ind_t = 0
    day = toint(cd_string(t6h(d), "%D"))
    hour = toint(cd_string(t6h(d), "%H"))
    do ip = 0, npoints-1
    if(dd(ip).eq.day) .and. (hh(ip).eq.hour) then
    break
    end if
    ind_t = ind_t + 1
    end do
    print(ind_t)
    ;if (d .gt. 1) then
    ind_p2 = ind(t1h.eq.t6h(d-1))+3 ;3H before
    ;end if
    res1@cnMinLevelValF            = 20.           ; min level
    res1@cnMaxLevelValF            = 80.           ; max level
    res1@cnLevelSpacingF           = 10.            ; contour interval
    res1@cnFillPalette             = "CBR_wet"
    res1@gsnScalarContour          = True
  
    res1@vcRefMagnitudeF           = 20.0            ; define vector ref mag
    res1@vcRefLengthF              = 0.045           ; define length of vec ref
    res1@vcRefAnnoOrthogonalPosF   = -1.0            ; move ref vector
    res1@vcRefAnnoArrowLineColor   = "black"         ; change ref vector color
    res1@vcRefAnnoArrowUseVecColor = False           ; don't use vec color for ref

    res1@vcGlyphStyle = "LineArrow"
    res1@vcLineArrowColor        = "black"           ; change vector color
    res1@vcLineArrowThicknessF   = 2.5               ; change vector thickness
    res1@vcMinDistanceF	         = 0.05
    res1@vcMinMagnitudeF         = 5.0
    res1@vcVectorDrawOrder       = "Draw"
    if (d.eq.ds) then
      res1@gsnCenterString            = "EST"
    else
      res1@gsnCenterString            = ""
    end if
    if (d.eq.de) then
      res1@lbLabelBarOn = True
    end if
    ;plot(i) = gsn_csm_contour_map_ce(wks, lh_1(d,:,:), res1)
    plot(i) = gsn_csm_vector_scalar_map_ce(wks, u(d,:,:), v(d,:,:), lh_1(d,:,:), res1)
    res2@cnLevels = fspan(20, 180, 5)
    res2@cnLineColor = "green"
    pr = pr_1(ind_p,:,:)
    pr = pr - pr_1(ind_p2,:,:)
    plot2 = gsn_csm_contour(wks, pr, res2)
    overlay(plot(i),plot2)
    dum(i) = gsn_add_polymarker(wks, plot(i), lon1(ind_t), lat1(ind_t), pmres)
    delete(plot2)
    delete(res2@cnLevels)
    i = i + 1
    lh = lh_2(d,:,:)
    lh = lh - lh_1(d,:,:)
    if (ds .lt. 9) then
    res1@cnMinLevelValF            = -5.            ; min level
    res1@cnMaxLevelValF            = 5.         ; max level
    res1@cnLevelSpacingF           = 1.           ; contour interval
    else
    res1@cnMinLevelValF            = -10.            ; min level
    res1@cnMaxLevelValF            = 10.         ; max level
    res1@cnLevelSpacingF           = 2.           ; contour interval
    end if
    res1@cnFillPalette             = "CBR_drywet"
    res1@gsnScalarContour          = False
    res1@lbLabelBarOn = False
    if (d.eq.ds) then
      res1@gsnCenterString            = "CLIM-EST"
    else
      res1@gsnCenterString            = ""
    end if
    plot(i) = gsn_csm_contour_map_ce(wks, lh, res1)
    res2@cnLevels = (/-30., -20., -10., 10., 20., 30./)
    res2@cnLineColor = "black"
    pr = pr_2(ind_p,:,:)
    pr = pr - pr_1(ind_p,:,:)
    ;if (d.gt.1) then
    pr = pr - pr_2(ind_p2,:,:) + pr_1(ind_p2,:,:)
    ;end if
    plot2 = gsn_csm_contour(wks, pr, res2)
    overlay(plot(i),plot2)
    dum(i) = gsn_add_polymarker(wks, plot(i), lon2(ind_t), lat2(ind_t), pmres)
    delete(plot2)
    i = i + 1
    lh = lh_3(d,:,:)
    lh = lh - lh_1(d,:,:)
    if (d.eq.ds) then
      res1@gsnCenterString            = "MGD-EST"
    else
      res1@gsnCenterString            = ""
    end if
    if (d.eq.de) then
      res1@lbLabelBarOn = True
      res1@pmLabelBarWidthF          = 1.05
      res1@pmLabelBarParallelPosF     = -0.02
    end if
    res1@gsnScalarContour          = False
    plot(i) = gsn_csm_contour_map_ce(wks, lh, res1)
    pr = pr_3(ind_p,:,:)
    pr = pr - pr_1(ind_p,:,:)
    ;if (d.gt.1) then
    pr = pr - pr_3(ind_p2,:,:) + pr_1(ind_p2,:,:)
    ;end if
    plot2 = gsn_csm_contour(wks, pr, res2)
    overlay(plot(i),plot2)
    dum(i) = gsn_add_polymarker(wks, plot(i), lon3(ind_t), lat3(ind_t), pmres)
    delete(plot2)
    i = i + 1
    lh = lh_2(d,:,:)
    lh = lh - lh_3(d,:,:)
    if (d.eq.ds) then
      res1@gsnCenterString            = "CLIM-MGD"
    else
      res1@gsnCenterString            = ""
    end if
    res1@lbLabelBarOn = False
    res1@gsnScalarContour          = False
    ;plot(i) = gsn_csm_contour_map_ce(wks, lh, res1)
    pr = pr_2(ind_p,:,:)
    pr = pr - pr_3(ind_p,:,:)
    ;if (d.gt.1) then
    pr = pr - pr_2(ind_p2,:,:) + pr_3(ind_p2,:,:)
    ;end if
    ;plot2 = gsn_csm_contour(wks, pr, res2)
    ;overlay(plot(i),plot2)
    ;dum(i) = gsn_add_polymarker(wks, plot(i), lon2(ind_t), lat2(ind_t), pmres)
    ;delete(plot2)
    delete(res2@cnLevels)
    ;i = i + 1
;    delete(ud)
;    delete(vd)
  end do
;****************************************************************************
; Draw the panel and label bars
;****************************************************************************
    resP                = True
    resP@gsnMaximize    = True
    resP@gsnFrame       = False
    resP@gsnPanelLeft   = 0.005
    resP@gsnPanelBottom = 0.08 
    resP@gsnPanelSave   = True      ; Save the state of the paneled plots so we can
                                    ; query for text strings locations.
   ;resP@gsnPanelMainString = "Land Only: Tmax: 2 Colorbars"
    resP@gsnPanelLabelBar = False

    gsn_panel(wks, plot, (/nrow,ncol/), resP)  ; plots
    frame(wks)
    delete(wks)
    delete(res1)
    delete(res2)
;  end do
end
