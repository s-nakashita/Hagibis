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
   yyyymmddhh = 2019100912
   sst = ""
   fh1     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100"+sst+"/fcst_surf_asia_"+yyyymmddhh+".nc", "r")
   fh2     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100"+sst+"/fcst_phy2m_avr_"+yyyymmddhh+".nc","r")
   ;fh3     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_mgdsst/fcst_surf_2019100900.nc", "r")
   
   psea = fh1->PSEA(:, :, :)
   flux = fh2->aRLUT(:, :, :)
   t1h  = fh1->time
   t6h  = fh2->time

   date = cd_string(t6h, "%N%D %H:%M")
   ntim = dimsizes(date)
   
   pltDir = "./"                ; where output will be written
   pltName= "GSMtl959"+sst+"_psea+rlut_"+yyyymmddhh                                   
   pltType= "png"               ; ps, eps, png, pdf, x11 .....

;****************************************************************************
; Graphic resources
;****************************************************************************
    ;ver = 4
  ds = 1
  de = ntim-2
  i = 2
  pltPath = pltDir+pltName+"_"+cd_string(t6h(ds), "%D%H")+"-"+cd_string(t6h(de), "%D%H")
  wks     = gsn_open_wks(pltType, pltPath) 
  nrow = 4
  ncol = 4
  plot = new(nrow*ncol , "graphic")
  cmap = read_colormap_file("spread_15lev")

  do d=ds, de
    res1                            = True     ; plot mods desired
    res1@gsnDraw                    = False    ; will panel plots later
    res1@gsnFrame                   = False

    res1@cnFillOn                   = True     ; turn on color fill
    res1@cnLinesOn                  = False    ; turn off contour lines
    res1@cnLineLabelsOn             = False    ; turn off contour line labels
    res1@cnFillPalette              = cmap(::-1,:)

;**********************************************************************
; We are turning the labelbar off for all but the bottom row of plots
;**********************************************************************
    res1@lbLabelBarOn               = False

    ; Set resources for map
    res1@mpMinLatF                  = -10.
    res1@mpMaxLatF                  = 55.
    res1@mpMinLonF                  = 120.
    res1@mpMaxLonF                  = 180.

    res1@mpDataBaseVersion          = "MediumRes"
    res1@mpOutlineOn                = True
    res1@mpFillOn                   = True                ; turn off gray fill
    ;res@mpOutlineBoundarySets      = "National"

    res1@gsnAddCyclic               = False            ; regional data, don't add cyclic point
 
    ; Set resources for tickmarks
    ;res1@gsnStringFont              = "helvetica-bold"  ; will be used for gsnCenterString later
    ;res1@gsnCenterStringFontHeightF = 0.04              ; ditto
    res1@gsnLeftString              = ""                ; remove the left and right substrings
    res1@gsnRightString             = ""
    res1@gsnCenterString            = date(d)

   ; Set resources specifically for plots in column 1
    res1@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
    res1@cnMinLevelValF            = 100.            ; min level
    res1@cnMaxLevelValF            = 280.           ; max level
    res1@cnLevelSpacingF           = 20.            ; contour interval

    res2 = True
    res2@gsnDraw = False
    res2@gsnFrame = False
    res2@gsnLeftString              = ""
    res2@gsnRightString             = ""
    res2@cnLinesOn                  = True    ; turn off contour lines
    res2@cnLineLabelsOn             = True
    res2@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
    res2@cnMinLevelValF            = 1000            ; min level
    res2@cnMaxLevelValF            = 1020           ; max level
    res2@cnLevelSpacingF           = 4            ; contour interval
;****************************************************************************
; Create the plots for the panel
;****************************************************************************
    ; Center strings will be set for top plots only
    ind_p = ind(t1h.eq.t6h(d))
    plot(i) = gsn_csm_contour_map_ce(wks, flux(d,:,:), res1)
    plot2 = gsn_csm_contour(wks, psea(ind_p,:,:), res2)
    overlay(plot(i),plot2)
    delete(plot2)
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
