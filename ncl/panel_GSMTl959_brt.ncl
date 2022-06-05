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
  sigma = 5.67e-8 ; stefan-boltzman constant [W/m^2/K^4]
  ;Reading Historical data
  yyyymmddhh = 2019100912
  sst = ""
  fh1     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100"+sst+"/fcst_phy2m_avr_"+yyyymmddhh+".nc","r")
   
  t6h  = fh1->time
  date = cd_string(t6h, "%Y%N%D %H:%M")
  ntim = dimsizes(date)
   
  pltDir = "./"                ; where output will be written
  pltName= "GSMtl959"+sst+"_tbb_"+yyyymmddhh                                   
  pltType= "png"               ; ps, eps, png, pdf, x11 .....
  d = ind(date.eq."20191012 00:00")
  print(d)
  pltPath = pltDir+pltName+"_"+cd_string(t6h(d), "%D%H")
  
;****************************************************************************
; Convert to brightness temperature
;****************************************************************************
  rlut = fh1->aRLUT(d, :, :) ; outgoing long wave radiation [W/m^2]
  tbb = new(dimsizes(rlut), typeof(rlut))
  tbb = sqrt(rlut / sigma)
  tbb = sqrt(tbb)
  tbb!0 = "lat"
  tbb!1 = "lon"
  tbb&lat = rlut&lat
  tbb&lon = rlut&lon
  tbb@long_name = "Brightness temperature"
  tbb@units     = "K"
;****************************************************************************
; Graphic resources
;****************************************************************************
  wks     = gsn_open_wks(pltType, pltPath) 
  ;plot = new(4,graphic)
  cmap = read_colormap_file("MPL_gist_yarg")
  res1                            = True     ; plot mods desired
  ;res1@gsnDraw                    = False    ; will panel plots later
  ;res1@gsnFrame                   = False

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
  ;res1@gsnStringFont              = "helvetica-bold"  ; will be used for gsnCenterString later
  ;res1@gsnCenterStringFontHeightF = 0.02              ; ditto
  res1@gsnLeftString              = ""                ; remove the left and right substrings
  res1@gsnRightString             = ""

;****************************************************************************
; Create the plots for the panel
;****************************************************************************
  ; Set resources specifically for plots in column 1
  res1@cnFillPalette             = cmap;(::-1,:)
  res1@gsnCenterString           = date(d)
  res1@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
  res1@cnMinLevelValF            = 200.            ; min level
  res1@cnMaxLevelValF            = 320.           ; max level
  res1@cnLevelSpacingF           = 10.            ; contour interval
  plot = gsn_csm_contour_map_ce(wks, tbb, res1)
;****************************************************************************
; Draw the panel and label bars
;****************************************************************************
;  resP                = True
;  resP@gsnMaximize    = True
;  resP@gsnFrame       = False
;  resP@gsnPanelLeft   = 0.005
;  resP@gsnPanelBottom = 0.08 
;  resP@gsnPanelSave   = True      ; Save the state of the paneled plots so we can
;                                    ; query for text strings locations.
;  resP@gsnPanelMainString = date(d)
;;  resP@gsnPanelLabelBar = True
;
;  gsn_panel(wks, plot, (/2,2/), resP)  ; plots
;  frame(wks)
;    delete(wks)
end
