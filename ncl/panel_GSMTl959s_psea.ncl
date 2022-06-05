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

   fh1     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_mgdsst_0912/fcst_surf_asia_2019100912.nc", "r")
   fh2     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_mgdsst_0912/fcst_p_asia_2019100912.nc","r")
   ;fh3     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_mgdsst/fcst_surf_2019100900.nc", "r")
   
   psea = fh1->PSEA(:, :, :)
   u    = fh2->U(:, {200}, :, :)
   v    = fh2->V(:, {200}, :, :)
   t1h  = fh1->time
   t3h  = fh2->time

   date = cd_string(t3h, "%N%D %H:%M")
   ntim = dimsizes(date)
   
   pltDir = "./"                ; where output will be written
   pltName= "GSMtl959_mgdsst_psea_2019100912"                                   
   pltType= "png"               ; ps, eps, png, pdf, x11 .....

;****************************************************************************
; Graphic resources
;****************************************************************************
    ;ver = 4
  ds = 0
  de = 10
  i = 2
  pltPath = pltDir+pltName+"_"+cd_string(t3h(ds), "%D%H")+"-"+cd_string(t3h(de), "%D%H")
  wks     = gsn_open_wks(pltType, pltPath) 
  nrow = 2
  ncol = 4
  plot = new(nrow*ncol , "graphic")

  do d=ds, de, 2
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

    res@cnFillOn                   = True     ; turn on color fill
    res@cnLinesOn                  = True    ; turn off contour lines
    res@cnLineLabelsOn             = True    ; turn off contour line labels
    res@cnFillPalette              = "GMT_polar"

;**********************************************************************
; We are turning the labelbar off for all but the bottom row of plots
;**********************************************************************
    res@lbLabelBarOn               = False

    ; Set resources for map
    res@mpMinLatF                  = -10.
    res@mpMaxLatF                  = 30.
    res@mpMinLonF                  = 110.
    res@mpMaxLonF                  = 180.

    res@mpDataBaseVersion          = "MediumRes"
    res@mpOutlineOn                = True
    res@mpFillOn                   = False                ; turn off gray fill
    ;res@mpOutlineBoundarySets      = "National"

    res@gsnAddCyclic               = False            ; regional data, don't add cyclic point
 
    ; Set resources for tickmarks
    res@gsnStringFont              = "helvetica-bold"  ; will be used for gsnCenterString later
    res@gsnCenterStringFontHeightF = 0.04              ; ditto
    res@gsnLeftString              = ""                ; remove the left and right substrings
    res@gsnRightString             = ""
    res@gsnCenterString            = date(d)

   ; Set resources specifically for plots in column 1
    res@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
    res@cnMinLevelValF            = 1000            ; min level
    res@cnMaxLevelValF            = 1020            ; max level
    res@cnLevelSpacingF           = 1.              ; contour interval

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
    ind_p = ind(t1h.eq.t3h(d))
    ;plot(i) = gsn_csm_vector_scalar_map_ce(wks, ud, vd, psea(ind_p,:,:), res)
    plot(i) = gsn_csm_contour_map_ce(wks, psea(ind_p,:,:), res)
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
