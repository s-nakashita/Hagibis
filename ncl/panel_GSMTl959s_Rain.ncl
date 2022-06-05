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
   init = "2019100912"
   fh1     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100/fcst_surf_"+init+".nc", "r")
   fh2     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_est/fcst_surf_"+init+".nc","r")
   fh3     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_mgdsst/fcst_surf_"+init+".nc", "r")
   
   pr1 = fh1->RAIN(:, :, :)
   pr2 = fh2->RAIN(:, :, :)
   pr3 = fh3->RAIN(:, :, :)

   date = cd_string(fh1->time, "%D-%HZ")
   ntim = dimsizes(date)
   
   pltDir = "./"                ; where output will be written
   pltName= "GSMtl959_rain_"+init                                   
   pltType= "png"               ; ps, eps, png, pdf, x11 .....

;****************************************************************************
; Graphic resources
;****************************************************************************
    ;ver = 4
  do ver=1, 5
    pltPath = pltDir+pltName+sprinti("%d", ver)
    wks     = gsn_open_wks(pltType, pltPath) 

    res                            = True     ; plot mods desired
    res@gsnDraw                    = False    ; will panel plots later
    res@gsnFrame                   = False

    res@cnFillOn                   = True     ; turn on color fill
    res@cnLinesOn                  = False    ; turn off contour lines
    res@cnLineLabelsOn             = False    ; turn off contour line labels

;**********************************************************************
; We are turning the labelbar off for all but the bottom row of plots
;**********************************************************************
    res@lbLabelBarOn               = False

    ; Set resources for when we turn on labelbar later
    res@lbTitleFontHeightF         = .03                 ; make title smaller
    res@lbTitleString              = "mm/3h"            ; title string
    res@lbTitlePosition            = "Bottom"            ; title position
    res@lbTitleDirection           = "Across"            ; title direction
    res@lbLabelStride              = 2
    res@lbLabelFontHeightF         = 0.025
    res@pmLabelBarHeightF          = 0.1
    res@lbBoxEndCapStyle           = "TriangleBothEnds"
    res@lbOrientation              = "Horizontal"        ; orientation

    ; Set resources for map
    res@mpProjection               = "CylindricalEquidistant"       ; projection
    res@mpLimitMode                = "Corners"                      ; method to zoom
    if (ver .lt. 3) then
    res@mpLeftCornerLatF           = 10.0
    res@mpRightCornerLatF          = 30.0
    else if (ver .lt. 4) then
    res@mpLeftCornerLatF           = 15.0
    res@mpRightCornerLatF          = 35.0
    else
    res@mpLeftCornerLatF           = 20.0
    res@mpRightCornerLatF          = 40.0
    end if
    end if
    res@mpLeftCornerLonF           = 135.0
    res@mpRightCornerLonF          = 145.0
    res@mpFillDrawOrder            = "PostDraw"

    res@mpDataBaseVersion          = "MediumRes"
    res@mpFillOn                   = True    
    res@mpOutlineOn                = True
    res@mpFillOn                   = False                ; turn off gray fill
    res@mpOutlineBoundarySets      = "National"

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

   ; Set resources specifically for plots in column 1
    res1                           = res
    res1@cnFillPalette             = "MPL_YlGnBu"
    res1@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
    res1@cnMinLevelValF            = 10             ; min level
    res1@cnMaxLevelValF            = 90             ; max level
    res1@cnLevelSpacingF           = 10.             ; contour interval
   
;****************************************************************************
; Create the plots for the panel
;****************************************************************************
    nrow = 3
    ncol = 6                           
    t_offset = ncol*(ver-1)
    plot = new(nrow*ncol , "graphic")

    ; Center strings will be set for top plots only
    do t = 1, ncol
    if (t+t_offset .lt. ntim) then
    if (t .eq. ncol) then
    res1@tmYRLabelsOn               = True            ; turn on lat labels
    res1@tmXBLabelsOn               = True
    res1@tmYROn                     = True
    res1@tmXBOn                     = True
    res1@tmYRLabelFontHeightF       = 0.02
    res1@tmXBLabelFontHeightF       = 0.02
    end if
    res1@gsnCenterString = date(t+t_offset)
    pr = pr1(t+t_offset,:,:)
    if (t+t_offset .gt. 1) then
        pr = pr - pr1(t+t_offset-1,:,:)
    end if
    printMinMax(pr, 0)
    plot(t-1) = gsn_csm_contour_map(wks, pr, res1)
    end if
    end do
   
    ; Turn off center strings for rest of plots
    res1@gsnCenterString = ""
    res1@tmYRLabelsOn               = False            ; turn off lat labels
    res1@tmXBLabelsOn               = False
    res1@tmYROn                     = False
    res1@tmXBOn                     = False
    do t = 1, ncol
    if (t+t_offset .lt. ntim) then
    pr = pr2(t+t_offset,:,:)
    if (t+t_offset .gt. 1) then
        pr = pr - pr2(t+t_offset-1,:,:)
    end if
    plot(ncol+t-1) = gsn_csm_contour_map(wks, pr, res1)
    end if
    end do
    
    ; Turn on labelbar for bottom four plots
    res1@lbLabelBarOn          = True
    do t = 1, ncol
    if (t+t_offset .lt. ntim) then
    pr = pr3(t+t_offset,:,:)
    if (t+t_offset .gt. 1) then
        pr = pr - pr3(t+t_offset-1,:,:)
    end if
    plot(2*ncol+t-1) = gsn_csm_contour_map(wks, pr, res1)
    end if
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

    gsn_panel(wks, plot, (/nrow,ncol/), resP)  ; plots
  
;**********************************************************************
; Add rotated text strings on left side of paneled plots.
; Use the resized plots to determine the position of the titles
;**********************************************************************
    pos1 = get_title_position(plot(0))         ; 1st row leftmost plot
    pos2 = get_title_position(plot(ncol))         ; 2nd row leftmost plot
    pos3 = get_title_position(plot(2*ncol))         ; 3rd row leftmost plot
    ;pos4 = get_title_position(plot(12))        ; 4th row leftmost plot
    txres               = True
    txres@txFontHeightF = 0.012
    txres@txAngleF      = 90.
    gsn_text_ndc(wks, "climatology",   pos1(0), pos1(1), txres)
    gsn_text_ndc(wks, "estimate",   pos2(0), pos2(1), txres)
    gsn_text_ndc(wks, "MGDSST", pos3(0), pos3(1), txres)
    ;gsn_text_ndc(wks, "CCLM (MPI-ESM)",    pos4(0), pos4(1), txres)
     
    frame(wks)
    delete(wks)
  end do
end
