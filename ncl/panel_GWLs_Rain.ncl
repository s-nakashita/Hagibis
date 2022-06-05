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
begin 
   ;Reading Historical data

   fh1     = addfile("mean_INT_pr_CLMcom-CCLM4_CNRM-CERFACS-CNRM-CM5_1971_2000.nc", "r")
   fh2     = addfile("mean_INT_pr_CLMcom-CCLM4_ICHEC-EC-EARTH_1971_2000.nc", "r")
   fh3     = addfile("mean_INT_pr_CLMcom-CCLM4_MOHC-HadGEM2-ES_1971_2000.nc", "r")
   fh4     = addfile("mean_INT_pr_CLMcom-CCLM4_MPI-M-MPI-ESM_1971_2000.nc", "r")
   
   hpr1 = fh1->pr(0, :, :)
   hpr2 = fh2->pr(0, :, :)
   hpr3 = fh3->pr(0, :, :)
   hpr4 = fh4->pr(0, :, :)
   
   ;Reading GWL  1.5°C
   
    gg1  = addfile("mean15_INT_pr_AFR-44_CNRM-CERFACS-CNRM-CM5_rcp85_r1i1p1_CLMcom-CCLM4-8-17_v1_day_20150101-20441231.nc", "r")
    gg2  = addfile("mean15_INT_pr_AFR-44_ICHEC-EC-EARTH_rcp85_r12i1p1_CLMcom-CCLM4-8-17_v1_day_20060101-20341231.nc", "r")
    gg3  = addfile("mean15_INT_pr_AFR-44_MOHC-HadGEM2-ES_rcp85_r1i1p1_CLMcom-CCLM4-8-17_v1_day_20100101-20391230.nc", "r")
    gg4  = addfile("mean15_INT_pr_AFR-44_MPI-M-MPI-ESM-LR_rcp85_r1i1p1_CLMcom-CCLM4-8-17_v1_day_20060101-20331231.nc", "r")
    
    gpr1 = gg1->pr(0, :, :)
    gpr2 = gg2->pr(0, :, :)
    gpr3 = gg3->pr(0, :, :)
    gpr4 = gg4->pr(0, :, :)
    
    ;Reading GWL  2.0°C 
    
    k1  = addfile("mean20_INT_pr_AFR-44_CNRM-CERFACS-CNRM-CM5_rcp85_r1i1p1_CLMcom-CCLM4-8-17_v1_day_20290101-20581231.nc", "r")
    k2  = addfile("mean20_INT_pr_AFR-44_ICHEC-EC-EARTH_rcp85_r12i1p1_CLMcom-CCLM4-8-17_v1_day_20210101-20501231.nc", "r")
    k3  = addfile("mean20_INT_pr_AFR-44_MOHC-HadGEM2-ES_rcp85_r1i1p1_CLMcom-CCLM4-8-17_v1_day_20230101-20521230.nc", "r")
    k4  = addfile("mean20_INT_pr_AFR-44_MPI-M-MPI-ESM-LR_rcp85_r1i1p1_CLMcom-CCLM4-8-17_v1_day_20210101-20501231.nc", "r")
       
    kpr1 = k1->pr(0, :, :)
    kpr2 = k2->pr(0, :, :)
    kpr3 = k3->pr(0, :, :)
    kpr4 = k4->pr(0, :, :)
    
    pltDir = "./"                ; where output will be written
    pltName= "GWL_rain"                                   
    pltType= "png"               ; ps, eps, png, pdf, x11 .....

;****************************************************************************
; Graphic resources
;****************************************************************************
    pltPath = pltDir+pltName
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
    res@lbTitleString              = "mm/day"            ; title string
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
    res@mpLeftCornerLatF           = min(kpr1&lat)
    res@mpLeftCornerLonF           = min(kpr1&lon)
    res@mpRightCornerLatF          = max(kpr1&lat)       
    res@mpRightCornerLonF          = max(kpr1&lon)
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
    res1@cnMinLevelValF            = 1             ; min level
    res1@cnMaxLevelValF            = 12             ; max level
    res1@cnLevelSpacingF           =  1.             ; contour interval
   
   ; Set resources specifically for plots in columns 2 & 3
    res2                           = res
    res2@cnFillPalette             = "MPL_PuOr"
    res2@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
    res2@cnMinLevelValF            = -1.5             ; min level
    res2@cnMaxLevelValF            =  1.5             ; max level
    res2@cnLevelSpacingF           =  0.2            ; contour interval
    res2@pmLabelBarWidthF          = 1.05
    res2@pmLabelBarParallelPosF     = -0.02

   ; Set resources specifically for plots in column 4
    cmap = read_colormap_file("MPL_bwr") 
    cmap = cmap(::-1,:) ; reverse the color map 

    res3                           = res
    res3@cnFillPalette             = cmap
    res3@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
    res3@cnMinLevelValF            = -1.             ; min level
    res3@cnMaxLevelValF            =  1.             ; max level
    res3@cnLevelSpacingF           =  0.2            ; contour interval

;****************************************************************************
; Create the plots for the panel
;****************************************************************************
    nrow = 4
    ncol = 4                           
    plot = new(nrow*ncol , "graphic")
  
;---Calculate the difference variables
    gpr1_hpr1_diff = diff_vars(gpr1,hpr1)
    kpr1_hpr1_diff = diff_vars(kpr1,hpr1)
    kpr1_gpr1_diff = diff_vars(kpr1,gpr1)
    gpr2_hpr2_diff = diff_vars(gpr2,hpr2)
    kpr2_hpr2_diff = diff_vars(kpr2,hpr2)
    kpr2_gpr2_diff = diff_vars(kpr2,gpr2)
    gpr3_hpr3_diff = diff_vars(gpr3,hpr3)
    kpr3_hpr3_diff = diff_vars(kpr3,hpr3)
    kpr3_gpr3_diff = diff_vars(kpr3,gpr3)
    kpr4_gpr4_diff = diff_vars(kpr4,gpr4)
    gpr4_hpr4_diff = diff_vars(gpr4,hpr4)
    kpr4_hpr4_diff = diff_vars(kpr4,hpr4)

    ; Center strings will be set for top plots only
    res1@gsnCenterString = "CTL:1971-2000"
    plot(0) = gsn_csm_contour_map(wks, hpr1, res1)
   
    res2@gsnCenterString = "1.5~S~o~N~ - CTL"
    plot(1) = gsn_csm_contour_map(wks, gpr1_hpr1_diff, res2)

    res2@gsnCenterString = "2.0~S~o~N~ - CTL"
    plot(2) = gsn_csm_contour_map(wks, kpr1_hpr1_diff, res2)

    res3@gsnCenterString = "2.0~S~o~N~ - 1.5~S~o~N~"
    plot(3) = gsn_csm_contour_map(wks, kpr1_gpr1_diff, res3)
   
    ; Turn off center strings for rest of plots
    res1@gsnCenterString = ""
    res2@gsnCenterString = ""
    res3@gsnCenterString = ""
     
    plot(4) = gsn_csm_contour_map(wks, hpr2, res1)
    plot(5) = gsn_csm_contour_map(wks, gpr2_hpr2_diff,  res2)
    plot(6) = gsn_csm_contour_map(wks, kpr2_hpr2_diff, res2)
    plot(7) = gsn_csm_contour_map(wks, kpr2_gpr2_diff,  res3)
   
    plot(8) = gsn_csm_contour_map(wks, hpr2, res1)
    plot(9) = gsn_csm_contour_map(wks, gpr3_hpr3_diff,  res2)
    plot(10) = gsn_csm_contour_map(wks, kpr3_hpr3_diff, res2)
    plot(11) = gsn_csm_contour_map(wks, kpr3_gpr3_diff,  res3)
   
    ; Turn on labelbar for bottom four plots
    res1@lbLabelBarOn          = True
   
    plot(12) = gsn_csm_contour_map(wks, hpr4, res1)
    plot(13) = gsn_csm_contour_map(wks, gpr4_hpr4_diff,  res2)
    
    res2@lbLabelBarOn       = True
    plot(14) = gsn_csm_contour_map(wks, kpr4_hpr4_diff, res2)
    
    res3@lbLabelBarOn       = True
    plot(15) = gsn_csm_contour_map(wks, kpr4_gpr4_diff, res3)
  
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
    pos2 = get_title_position(plot(4))         ; 2nd row leftmost plot
    pos3 = get_title_position(plot(8))         ; 3rd row leftmost plot
    pos4 = get_title_position(plot(12))        ; 4th row leftmost plot
    txres               = True
    txres@txFontHeightF = 0.012
    txres@txAngleF      = 90.
    gsn_text_ndc(wks, "CCLM (CNRM-CM5)",   pos1(0), pos1(1), txres)
    gsn_text_ndc(wks, "CCLM (EC-EARTH)",   pos2(0), pos2(1), txres)
    gsn_text_ndc(wks, "CCLM (HadGEM2-ES)", pos3(0), pos3(1), txres)
    gsn_text_ndc(wks, "CCLM (MPI-ESM)",    pos4(0), pos4(1), txres)
     
    frame(wks)
end
