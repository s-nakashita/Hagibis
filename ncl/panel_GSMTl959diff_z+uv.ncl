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
   ;Reading Historical data

   sst = "est"
   fh1     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100/fcst_p_asia_2019100900.nc", "r")
   fh2     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_"+sst+"/fcst_p_asia_2019100900.nc","r")
   ;fh3     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_mgdsst/fcst_surf_2019100900.nc", "r")
   
   ;level = 850
   z   = fh1->Z(:, {level}, :, :)
   zp  = fh2->Z(:, {level}, :, :)
   t3h = fh1->time
   printMinMax(z, 0)
   
   date = cd_string(t3h, "%N%D %H:%M")
   ntim = dimsizes(date)
   
   pltDir = "./"                ; where output will be written
   pltName= "GSMtl959_diff_"+sst+"_z"+sprinti("%3d",level)+"_2019100900"
   pltType= "png"               ; ps, eps, png, pdf, x11 .....

;****************************************************************************
; Graphic resources
;****************************************************************************
    ;ver = 4
  ;ds = 16
  de = min((/ds + 15, 28/)) ;max 28
  i = 0
  pltPath = pltDir+pltName+"_"+cd_string(t3h(ds), "%D%H")+"-"+cd_string(t3h(de), "%D%H")
  wks     = gsn_open_wks(pltType, pltPath) 
  gsn_define_colormap(wks, "GMT_polar")
  nrow = 4
  ncol = 4
  plot = new(nrow*ncol , "graphic")

  do d=ds, de
;;****************************************************************************
;; Calculate divergent winds
;;****************************************************************************
;   dv     = u(d,:,:)
;   dv     = uv2dvF(u(d,:,:), v(d,:,:))
;   
;   ud     = u(d,:,:)
;   vd     = v(d,:,:)
;   uvd    = dv2uvF(dv)
;   ud     = (/uvd(0,:,:)/)
;   vd     = (/uvd(1,:,:)/)
;   copy_VarCoords(u(d,:,:), ud)
;   copy_VarCoords(v(d,:,:), vd)
;   ud@long_name = "Zonal Divergent Wind"
;   ud@units     = "m/s" ;u@units
;   vd@long_name = "Meridional Divergent Wind"
;   vd@units     = "m/s" ;v@units

    res1                            = True     ; plot mods desired
    res1@gsnDraw                    = False    ; will panel plots later
    res1@gsnFrame                   = False

    res1@cnFillOn                   = True     ; turn on color fill
    res1@cnLinesOn                  = False    ; turn off contour lines
    res1@cnLineLabelsOn             = False    ; turn off contour line labels
    ;res1@cnFillPalette              = "GMT_polar"
    res1@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
    res1@cnMinLevelValF            = -50.            ; min level
    res1@cnMaxLevelValF            = 50.            ; max level
    res1@cnLevelSpacingF           = 5.              ; contour interval

;**********************************************************************
; We are turning the labelbar off for all but the bottom row of plots
;**********************************************************************
    res1@lbLabelBarOn               = False

    ; Set resources for map
    if (d.eq.0) then
    res1@mpMinLatF                  = 0.
    res1@mpMaxLatF                  = 50.
    else
    res1@mpMinLatF                  = 10.
    res1@mpMaxLatF                  = 60.
    end if
    res1@mpMinLonF                  = 110.
    res1@mpMaxLonF                  = 180.
    res1@mpDataBaseVersion          = "MediumRes"
    res1@mpOutlineOn                = True
    res1@mpFillOn                   = False                ; turn off gray fill
    ;res@mpOutlineBoundarySets      = "National"

    res1@gsnAddCyclic               = False            ; regional data, don't add cyclic point
    res1@cnFillDrawOrder            = "Predraw"
 
    ; Set resources for tickmarks
    res1@gsnStringFont              = "helvetica-bold"  ; will be used for gsnCenterString later
    res1@gsnCenterStringFontHeightF = 0.04              ; ditto
    res1@gsnLeftString              = ""                ; remove the left and right substrings
    res1@gsnRightString             = ""
    res1@gsnCenterString            = date(d)

   ; Set resources specifically for plots in column 1
    res2 = True
    res2@gsnAddCyclic               = False            ; regional data, don't add cyclic point
    res2@cnFillOn                   = False    ; turn on color fill
    res2@cnLinesOn                  = True     ; turn off contour lines
    res2@cnLineLabelsOn             = True     ; turn off contour line labels
    res2@cnLevelSelectionMode      = "ExplicitLevels"  ; manually specify contour levels
    if (level .eq. 250) then
    ;; Z250
      res2@cnLevels          = (/9900.,10000.,10100.,10200.,10300.,10400.,10500.,10600.,10700.,10800.,10900.,11000./)
    else if (level .eq. 300) then
    ;; Z300
      res2@cnLevels          = (/8600.,8700.,8800.,8900.,9000.,9100.,9200.,9300.,9400.,9500.,9600.,9700./)
    else if (level .eq. 500) then
    ;; Z500
      res2@cnLevels          = (/5500.,5550.,5600.,5650.,5700.,5750.,5800.,5850.,5900.,5950./)
    else if (level .eq. 700) then
    ; Z700
      res2@cnLevels          = (/3060.,3080.,3100.,3120.,3140.,3160.,3180.,3200.,3220./)
    else if (level .eq. 850) then
    ;; Z850
     res2@cnLevels          = (/1300.,1350.,1400.,1450.,1500.,1550.,1600.,1650.,1700./)
    end if
    end if
    end if
    end if
    end if
    res2@cnMonoLineThickness = False
    ;if (level .ne. 300) then
    ;  res2@cnLineThicknesses = (/   4.,   1.,   2.,   1.,   2.,   1.,   2.,   1.,   2.,   1./)
    ;else
      res2@cnLineThicknesses = (/   2.,   1.,   2.,   1.,   2.,   1.,   2.,   1.,   2.,   1.,   2.,   1.,   2.,   1./)
    ;end if

    ;res@gsnScalarContour = True
    ;res@vcRefMagnitudeF           = 10.0            ; define vector ref mag
    ;res@vcRefLengthF              = 0.045           ; define length of vec ref
    ;res@vcRefAnnoOrthogonalPosF   = -1.0            ; move ref vector
    ;res@vcRefAnnoArrowLineColor   = "black"         ; change ref vector color
    ;res@vcRefAnnoArrowUseVecColor = False           ; don't use vec color for ref

    ;res@vcGlyphStyle = "CurlyVector"
    ;res@vcLineArrowColor        = "black"           ; change vector color
    ;res@vcLineArrowThicknessF   = 2.0               ; change vector thickness
    ;res@vcMinDistanceF	      = 0.03
    ;res@vcVectorDrawOrder       = "PostDraw"        ; draw vectors last
;****************************************************************************
; Create the plots for the panel
;****************************************************************************
    ; Center strings will be set for top plots only
    ;ind_p = ind(t1h.eq.t3h(d))
    zd = zp(d,:,:)
    zd = zp(d,:,:) - z(d,:,:)
    plot(i) = gsn_csm_contour_map_overlay(wks, zd, z(d,:,:), res1, res2)
    i = i + 1
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
    resP@gsnPanelLabelBar = True

    gsn_panel(wks, plot, (/nrow,ncol/), resP)  ; plots
    frame(wks)
    delete(wks)
end
