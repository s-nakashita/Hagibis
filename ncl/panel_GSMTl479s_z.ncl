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

   fh1     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est/fcst_p_asia_2019100912.nc", "r")
   fh2     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est_0912+p/fcst_p_asia_2019100912.nc","r")
   fh3     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est_0912+p2/fcst_p_asia_2019100912.nc","r")
   fh4     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est_0912+pn/fcst_p_asia_2019100912.nc", "r")
   ;fh3     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est_0912+pf/fcst_p_asia_2019100912.nc", "r")
   level = 500
   z1   = fh1->Z(:, {level}, :, :)
   z2   = fh2->Z(:, {level}, :, :)
   z3   = fh3->Z(:, {level}, :, :)
   z4   = fh4->Z(:, {level}, :, :)
   t    = fh1->time
   
   date = cd_string(t, "%D-%HZ")
   ntim = dimsizes(date)
   
   pltDir = "./"                ; where output will be written
   pltName= "GSMtl479s_z"+sprinti("%3d", level)                                   
   pltType= "png"               ; ps, eps, png, pdf, x11 .....

;****************************************************************************
; Graphic resources
;****************************************************************************
    ;ver = 4
  ssec = cd_inv_calendar(2019,10,10,0,0,0,t@units,0)
  esec = cd_inv_calendar(2019,10,10,9,0,0,t@units,0)
  ds = ind(t.eq.ssec)
  de = ind(t.eq.esec)
  i = 0
  pltPath = pltDir+pltName+"_"+cd_string(t(ds), "%D%H")+"-"+cd_string(t(de), "%D%H")
  wks     = gsn_open_wks(pltType, pltPath) 
  nrow = 2
  ncol = 2
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

    res                            = True     ; plot mods desired
    res@gsnDraw                    = False    ; will panel plots later
    res@gsnFrame                   = False

    res@cnFillOn                   = False    ; turn on color fill
    res@cnLinesOn                  = True    ; turn off contour lines
    res@cnLineLabelsOn             = True    ; turn off contour line labels
    res@cnInfoLabelOn              = False
    ;res@cnFillPalette              = "GMT_polar"

;**********************************************************************
; We are turning the labelbar off for all but the bottom row of plots
;**********************************************************************
    res@lbLabelBarOn               = False

    res@gsnAddCyclic               = False            ; regional data, don't add cyclic point

    ; Set resources for tickmarks
    ;res@tmYLLabelsOn               = False            ; turn off lat labels
    ;res@tmXBLabelsOn               = False            ; turn off lon labels
  
    res@tmYROn                     = False
    ;res@tmYLOn                     = False
    res@tmXTOn                     = False
    ;res@tmXBOn                     = False

    res@gsnStringFont              = "helvetica-bold"  ; will be used for gsnCenterString later
    res@gsnCenterStringFontHeightF = 0.04              ; ditto
    res@gsnLeftString              = ""                ; remove the left and right substrings
    res@gsnRightString             = ""
    res@gsnCenterString            = ""

   ; Set resources specifically for plots in column 1
    ;res@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
    ;res@cnMinLevelValF            = 1000            ; min level
    ;res@cnMaxLevelValF            = 1550            ; max level
    ;res@cnLevelSpacingF           = 50              ; contour interval
    res@cnLevelSelectionMode     = "ExplicitLevels"
    ;res@cnLevels = (/1000., 1200., 1400., 1450., 1500.,1550./)
    ;res@cnLevels = (/2800., 2900., 2950., 3000., 3050., 3100., 3150./)
    res@cnLevels = (/5500., 5600., 5700., 5750., 5800., 5850., 5900./)
    ;res@cnLevels = (/8600.,8700.,8800.,8900.,9000.,9100.,9200.,9300.,9400.,9500.,9600.,9700./)
    res@cnLineLabelInterval       = 2               ; all contour has the label
    res@cnLineLabelFontHeightF    = 0.01
    res@cnLineLabelPerimOn        = False
    res@cnLineThicknessF          = 3.0

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
    res1 = res
    ; Set resources for map
    res1@mpProjection = "CylindricalEquidistant"
    res1@mpMinLatF                  = 5.
    res1@mpMaxLatF                  = 35.
    res1@mpMinLonF                  = 120.
    res1@mpMaxLonF                  = 160.

    res1@mpDataBaseVersion          = "MediumRes"
    res1@mpOutlineOn                = True
    res1@mpFillOn                   = False                ; turn off gray fill
    ;res@mpOutlineBoundarySets      = "National"
;****************************************************************************
; Create the plots for the panel
;****************************************************************************
    ; Center strings will be set for top plots only
    ;dir  = cd_string(t(d), "%Y/%N/%D")
    ;dd   = cd_string(t(d), "%D")
    ;print(dir)
    ;fh0  = addfile("/Volumes/dandelion/netcdf/era5/"+dir+"/plev_"+dd+".nc","r")
    ;t0   = fh0->time
    ;print(cd_string(time(t+t_offset),"%Y%N%D%H"))
    ;print(cd_string(t0,"%Y%N%D%H"))
    ;ind0 = ind(cd_string(t0,"%Y%N%D%H").eq.cd_string(t(d),"%Y%N%D%H"))
    ;print(ind0)
    ;z0   = fh0->z(ind0,{level},:,:)
    ;z0   = z0 / 9.81
    ;printMinMax(z0,0)
    res1@gsnCenterString = date(d)
    res1@cnLineColor = "black"
    res1@cnLineLabelFontColor = "black"
    ;plot(i) = gsn_csm_contour_map(wks, z0, res1)
    plot(i) = gsn_csm_contour_map(wks, z1(d,:,:), res1)
    ;delete(fh0)
    ;delete(z0)
    res@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
    res@cnMinLevelValF            = -50             ; min level
    res@cnMaxLevelValF            = 50              ; max level
    res@cnLevelSpacingF           = 5               ; contour interval
    res@gsnContourNegLineDashPattern = 2
    res@gsnContourZeroLineThicknessF = 0.
    ;res@cnLineColor = "red"
    ;res@cnLineLabelFontColor = "red"
    res@cnLineColor = "blue"
    res@cnLineLabelFontColor = "blue"
    zd = z2(d,:,:)
    zd = z2(d,:,:) - z1(d,:,:)
    plot1 = gsn_csm_contour(wks, zd, res)
    overlay(plot(i), plot1)
    res@cnLineColor = "green"
    res@cnLineLabelFontColor = "green"
    zd = z3(d,:,:) - z1(d,:,:)
    plot2 = gsn_csm_contour(wks, zd, res)
    overlay(plot(i), plot2)
    res@cnLineColor = "cyan"
    res@cnLineLabelFontColor = "cyan"
    ;res@cnLineColor = "magenta"
    ;res@cnLineLabelFontColor = "magenta"
    zd = z4(d,:,:) - z1(d,:,:)
    plot3 = gsn_csm_contour(wks, zd, res)
    overlay(plot(i), plot3)
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
    ;resP@gsnPanelLabelBar = True

    gsn_panel(wks, plot, (/nrow,ncol/), resP)  ; plots
    frame(wks)
    delete(wks)
end
