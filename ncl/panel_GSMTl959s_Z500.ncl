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
   
   lev = 300
   var = "Z"
   var_l = str_lower(var)
   pr1 = fh1->$var$(:, {lev}, :, :)
   pr2 = fh2->$var$(:, {lev}, :, :)
   pr3 = fh3->$var$(:, {lev}, :, :)
   time = fh1->time

   date = cd_string(time, "%D-%HZ")
   ntim = dimsizes(date)
   
   pltDir = "./"                ; where output will be written
   pltName= "GSMtl959s_"+var_l+sprinti("%d",lev)+"_"                                   
   pltType= "png"               ; ps, eps, png, pdf, x11 .....

;****************************************************************************
; Graphic resources
;****************************************************************************
    ;ver = 4
  ;do ver=1, 7
    ;pltPath = pltDir+pltName+sprinti("%d", ver)
    t_offset = 20
    nrow = 3
    ncol = 3
    pltPath = pltDir+pltName+cd_string(fh1->time(t_offset), "%D%H")+"-"+cd_string(fh1->time(t_offset+nrow*ncol-1), "%D%H")
    wks     = gsn_open_wks(pltType, pltPath) 

    res                            = True     ; plot mods desired
    res@gsnDraw                    = False    ; will panel plots later
    res@gsnFrame                   = False

    res@cnFillOn                   = False     ; turn on color fill
    res@cnLinesOn                  = True    ; turn off contour lines
    res@cnLineLabelsOn             = True    ; turn off contour line labels
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
  
    res@tmYROn                     = False
    ;res@tmYLOn                     = False
    res@tmXTOn                     = False
    ;res@tmXBOn                     = False

    res@gsnStringFont              = "helvetica-bold"  ; will be used for gsnCenterString later
    res@gsnCenterStringFontHeightF = 0.04              ; ditto
    res@gsnLeftString              = ""                ; remove the left and right substrings
    res@gsnRightString             = ""
    res@gsnCenterString             = ""

   ; Set resources specifically for plots in column 1
    res@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
    ;res@cnMinLevelValF            = 5300.             ; min level
    ;res@cnMaxLevelValF            = 5900.             ; max level
    ;res@cnLevelSpacingF           = 50.              ; contour interval
    res@cnMinLevelValF            = 9000.             ; min level
    res@cnMaxLevelValF            = 9800.             ; max level
    res@cnLevelSpacingF           = 100.              ; contour interval
    ;res@cnLevelSelectionMode     = "ExplicitLevels"
    ;res@cnLevels = (/5500., 5600., 5700., 5740., 5780., 5820., 5860., 5900./)
    ;res@cnLevels = (/8600.,8700.,8800.,8900.,9000.,9100.,9200.,9300.,9400.,9500.,9600.,9700./)
    res@cnLineLabelInterval       = 2               ; all contour has the label
    res@cnLineLabelFontHeightF    = 0.01
    res@cnLineLabelPerimOn        = False
    res@cnLineThicknessF          = 3.0
   ; Set resources for map
    res1 = res
    res1@mpProjection               = "CylindricalEquidistant"       ; projection
    res1@mpMinLatF = 20.0
    res1@mpMaxLatF = 50.0
    res1@mpMinLonF = 120.0
    res1@mpMaxLonF = 150.0
    res1@mpFillDrawOrder            = "PostDraw"

    res1@mpDataBaseVersion          = "MediumRes"
    res1@mpFillOn                   = False    
    res1@mpOutlineOn                = True
    res1@mpOutlineBoundarySets      = "National"

;****************************************************************************
; Create the plots for the panel
;****************************************************************************
    t_int = 1 ; every 3 hours
    ;t_offset = nrow*ncol*t_int*(ver-1)
    plot = new(nrow*ncol , "graphic")

    ; Center strings will be set for top plots only
    i = 0
    do t = 0, nrow*ncol*t_int-1, t_int
    if (t+t_offset .lt. ntim) then
    ;dir  = cd_string(time(t+t_offset), "%Y/%N/%D")
    ;dd   = cd_string(time(t+t_offset), "%D")
    ;print(dir)
    ;fh0  = addfile("/Volumes/dandelion/netcdf/era5/"+dir+"/plev_"+dd+".nc","r")
    ;t0   = fh0->time
    ;;print(cd_string(time(t+t_offset),"%Y%N%D%H"))
    ;;print(cd_string(t0,"%Y%N%D%H"))
    ;ind0 = ind(cd_string(t0,"%Y%N%D%H").eq.cd_string(time(t+t_offset),"%Y%N%D%H"))
    ;print(ind0)
    ;pr   = fh0->$var_l$(ind0,{lev},:,:)
    ;pr   = pr / 9.81
    ;printMinMax(pr,0)
    res1@gsnCenterString = date(t+t_offset)
    ;res1@cnLineColor = "black"
    ;res1@cnLineLabelFontColor = "black"
    ;plot(i) = gsn_csm_contour_map(wks, pr, res1)
    ;delete(fh0)
    ;delete(pr)

    pr = pr1(t+t_offset,:,:)
    printMinMax(pr, 0)
    res1@cnLineColor = "red"
    res1@cnLineLabelFontColor = "red"
    plot(i) = gsn_csm_contour_map(wks, pr, res1)
    ;plot1 = gsn_csm_contour(wks, pr, res)
    ;overlay(plot(i), plot1)
    
    pr = pr2(t+t_offset,:,:)
    printMinMax(pr, 0)
    res@cnLineColor = "blue"
    res@cnLineLabelFontColor = "blue"
    plot2 = gsn_csm_contour(wks, pr, res)
    overlay(plot(i), plot2)
    
    pr = pr3(t+t_offset,:,:)
    res@cnLineColor = "green"
    res@cnLineLabelFontColor = "green"
    plot3 = gsn_csm_contour(wks, pr, res)
    overlay(plot(i), plot3)
    i = i + 1
    delete(pr)
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
  ;end do
end
