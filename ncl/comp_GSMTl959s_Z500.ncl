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

   fh1     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100/fcst_p_asia_2019100900.nc", "r")
   fh2     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_est/fcst_p_asia_2019100900.nc","r")
   fh3     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_mgdsst/fcst_p_asia_2019100900.nc", "r")
   
   pr1 = fh1->Z(:, {500}, :, :)
   pr2 = fh2->Z(:, {500}, :, :)
   pr3 = fh3->Z(:, {500}, :, :)

   date = cd_string(fh1->time, "%D-%HZ")
   ntim = dimsizes(date)
   
   pltDir = "./"                ; where output will be written
   pltName= "GSMtl959_est-mgd_z500"                                   
   pltType= "png"               ; ps, eps, png, pdf, x11 .....

;****************************************************************************
; Graphic resources
;****************************************************************************
    ;ver = 4
  do ver=1, 7
    pltPath = pltDir+pltName+sprinti("%d", ver)
    wks     = gsn_open_wks(pltType, pltPath) 
    gsn_define_colormap(wks, "GMT_Polar")

    res                            = True     ; plot mods desired
    res@gsnDraw                    = False    ; will panel plots later
    res@gsnFrame                   = False

    res@cnFillOn                   = True     ; turn on color fill
    res@cnLinesOn                  = True    ; turn off contour lines
    res@cnLineLabelsOn             = False    ; turn off contour line labels
    res@cnInfoLabelOn              = False

;**********************************************************************
; We are turning the labelbar off for all but the bottom row of plots
;**********************************************************************
    res@lbLabelBarOn               = False

    ; Set resources for when we turn on labelbar later
    res@lbTitleFontHeightF         = .03                 ; make title smaller
    res@lbTitleString              = "m"            ; title string
    res@lbTitlePosition            = "Bottom"            ; title position
    res@lbTitleDirection           = "Across"            ; title direction
    res@lbLabelStride              = 2
    res@lbLabelFontHeightF         = 0.025
    res@pmLabelBarHeightF          = 0.1
    res@lbBoxEndCapStyle           = "TriangleBothEnds"
    res@lbOrientation              = "Horizontal"        ; orientation

    res@gsnAddCyclic               = False            ; regional data, don't add cyclic point
 
    ; Set resources for tickmarks
    ;res@tmYLLabelsOn               = False            ; turn off lat labels
    ;res@tmXBLabelsOn               = False            ; turn off lon labels
  
    ;res@tmYROn                     = False
    ;res@tmYLOn                     = False
    ;res@tmXTOn                     = False
    ;res@tmXBOn                     = False

    res@gsnStringFont              = "helvetica-bold"  ; will be used for gsnCenterString later
    res@gsnCenterStringFontHeightF = 0.04              ; ditto
    res@gsnLeftString              = ""                ; remove the left and right substrings
    res@gsnRightString             = ""
    res@gsnCenterString             = ""

   ; Set resources specifically for plots in column 1
    res@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
    res@cnMinLevelValF            = -50             ; min level
    res@cnMaxLevelValF            = 50             ; max level
    res@cnLevelSpacingF           = 5.              ; contour interval
    res@cnLineLabelInterval       = 1               ; all contour has the label
    res@cnLineLabelFontHeightF    = 0.02
    res@gsnContourNegLineDashPattern = 1            ; negative contur is dashed
    res@cnLineLabelPerimOn = False
   ; Set resources for map
    res1 = res
    res1@mpProjection               = "CylindricalEquidistant"       ; projection
    res1@mpMinLatF = 30.0
    res1@mpMaxLatF = 50.0
    res1@mpMinLonF = 120.0
    res1@mpMaxLonF = 160.0
    res1@mpFillDrawOrder            = "PostDraw"

    res1@mpDataBaseVersion          = "MediumRes"
    res1@mpFillOn                   = False    
    res1@mpOutlineOn                = True
    res1@mpOutlineBoundarySets      = "National"

;****************************************************************************
; Create the plots for the panel
;****************************************************************************
    nrow = 2
    ncol = 2                           
    t_int = 1 ; every 3 hours
    t_offset = nrow*ncol*t_int*(ver-1)
    plot = new(nrow*ncol , "graphic")

    ; Center strings will be set for top plots only
    i = 0
    do t = 1, nrow*ncol*t_int, t_int
    if (t+t_offset .lt. ntim) then
    res1@gsnCenterString = date(t+t_offset)
    pr = diff_vars(pr2(t+t_offset,:,:), pr3(t+t_offset,:,:))
    printMinMax(pr, 0)
    if (i .eq. 3)
      res1@lbLabelBarOn               = True
    end if
    plot(i) = gsn_csm_contour_map(wks, pr, res1)
    i = i + 1
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
     
    frame(wks)
    delete(wks)
  end do
end
