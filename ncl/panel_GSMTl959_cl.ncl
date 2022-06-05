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
   yyyymmddhh = 2019100900
   sst = "_est"
   fh1     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100"+sst+"/fcst_phy2m_snp_"+yyyymmddhh+".nc","r")
   fh2     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100"+sst+"/fcst_p_asia_"+yyyymmddhh+".nc","r")
   ;fh3     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_mgdsst/fcst_surf_2019100900.nc", "r")
   
   cla  = fh1->CLA(:, ::-1, :) ; cloud fraction [m^2/m^2]
   cll  = fh1->CLL(:, ::-1, :) ; cloud fraction (Low) [m^2/m^2]
   clm  = fh1->CLM(:, ::-1, :) ; cloud fraction (Middle) [m^2/m^2]
   clh  = fh1->CLH(:, ::-1, :) ; cloud fraction (High) [m^2/m^2]
   t6h  = fh1->time
   
   ;Grid information
   longs = fh1->lon
   latgs = fh1->lat(::-1) ;monotonically increasing
   lonll = fh2->lon
   latll = fh2->lat(::-1) ;monotonically increasing
   printVarSummary(longs)
   printVarSummary(latgs)
   printVarSummary(lonll)
   printVarSummary(latll)

   date = cd_string(t6h, "%Y%N%D %H:%M")
   ntim = dimsizes(date)
   
   pltDir = "./"                ; where output will be written
   pltName= "GSMtl959"+sst+"_cl_"+yyyymmddhh                                   
   pltType= "png"               ; ps, eps, png, pdf, x11 .....
   d = ind(date.eq."20191009 00:00")
   print(d)
   pltPath = pltDir+pltName+"_"+cd_string(t6h(d), "%D%H")
;****************************************************************************
; Horizontal interpolation
;****************************************************************************
   cla_ll = linint2_Wrap(longs, latgs, cla(d,:,:), True, lonll, latll, 0)
   cla_ll!0    = "lat"
   cla_ll!1    = "lon"
   cla_ll&lat  = latll
   cla_ll&lon  = lonll
   cla_ll@long_name = "cloud fraction"
   cla_ll@units     = "m^2/m^2"
   printMinMax(cla_ll, 0)
   cll_ll = linint2_Wrap(longs, latgs, cll(d,:,:), True, lonll, latll, 0)
   cll_ll!0    = "lat"
   cll_ll!1    = "lon"
   cll_ll&lat  = latll
   cll_ll&lon  = lonll
   cll_ll@long_name = "cloud fraction [low]"
   cll_ll@units     = "m^2/m^2"
   printMinMax(cll_ll, 0)
   clm_ll = linint2_Wrap(longs, latgs, clm(d,:,:), True, lonll, latll, 0)
   clm_ll!0    = "lat"
   clm_ll!1    = "lon"
   clm_ll&lat  = latll
   clm_ll&lon  = lonll
   clm_ll@long_name = "cloud fraction [middle]"
   clm_ll@units     = "m^2/m^2"
   printMinMax(clm_ll, 0)
   clh_ll = linint2_Wrap(longs, latgs, clh(d,:,:), True, lonll, latll, 0)
   clh_ll!0    = "lat"
   clh_ll!1    = "lon"
   clh_ll&lat  = latll
   clh_ll&lon  = lonll
   clh_ll@long_name = "cloud fraction [high]"
   clh_ll@units     = "m^2/m^2"
   printMinMax(clh_ll, 0)
;****************************************************************************
; Graphic resources
;****************************************************************************
  wks     = gsn_open_wks(pltType, pltPath) 
  plot = new(4,graphic)
  cmap = read_colormap_file("WhBlGrYeRe")
  res1                            = True     ; plot mods desired
  res1@gsnDraw                    = False    ; will panel plots later
  res1@gsnFrame                   = False

  res1@cnFillOn                   = True     ; turn on color fill
  res1@cnLinesOn                  = False    ; turn off contour lines
  res1@cnLineLabelsOn             = False    ; turn off contour line labels

  ; Set resources for map
  res1@mpMinLatF                  = -10.
  res1@mpMaxLatF                  = 55.
  res1@mpMinLonF                  = 120.
  res1@mpMaxLonF                  = 180.


  res1@mpDataBaseVersion          = "MediumRes"
  res1@mpOutlineOn                = True
  res1@mpFillOn                   = False
  ;res@mpOutlineBoundarySets      = "National"

  res1@gsnAddCyclic               = False            ; regional data, don't add cyclic point
 
  ; Set resources for tickmarks
  res1@gsnStringFont              = "helvetica-bold"  ; will be used for gsnCenterString later
  res1@gsnCenterStringFontHeightF = 0.02              ; ditto
  res1@gsnLeftString              = ""                ; remove the left and right substrings
  res1@gsnRightString             = ""

;****************************************************************************
; Create the plots for the panel
;****************************************************************************
  ; Set resources specifically for plots in column 1
  res1@cnFillPalette             = cmap
  res1@gsnCenterString           = cla_ll@long_name
  res1@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
  res1@cnMinLevelValF            = .1            ; min level
  res1@cnMaxLevelValF            = .9          ; max level
  res1@cnLevelSpacingF           = .1            ; contour interval
  plot(0) = gsn_csm_contour_map_ce(wks, cla_ll, res1)
  ; Set resources specifically for plots in column 1
  res1@gsnCenterString           = cll_ll@long_name
  plot(1) = gsn_csm_contour_map_ce(wks, cll_ll, res1)
  ; Set resources specifically for plots in column 1
  res1@gsnCenterString           = clm_ll@long_name
  plot(2) = gsn_csm_contour_map_ce(wks, clm_ll, res1)
  ; Set resources specifically for plots in column 1
  res1@gsnCenterString           = clh_ll@long_name
  plot(3) = gsn_csm_contour_map_ce(wks, clh_ll, res1)
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
  resP@gsnPanelMainString = date(d)
;  resP@gsnPanelLabelBar = True
;
  gsn_panel(wks, plot, (/2,2/), resP)  ; plots
  frame(wks)
;    delete(wks)
end
