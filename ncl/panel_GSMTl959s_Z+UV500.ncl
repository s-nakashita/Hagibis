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
   
   t = fh1->time
   vt = cd_inv_calendar(2019, 10, 12, 0, 0, 0, t@units, 0)
   
   date = cd_string(vt, "%Y%N%D%H")
   
   pltDir = "./"                ; where output will be written
   pltName= "GSMtl959_z+uv500_"+date                                   
   pltType= "png"               ; ps, eps, png, pdf, x11 .....
   
;****************************************************************************
; Graphic resources
;****************************************************************************
    pltPath = pltDir+pltName
    wks     = gsn_open_wks(pltType, pltPath) 

    res                            = True     ; plot mods desired
    res@gsnDraw                    = False    ; will panel plots later
    res@gsnFrame                   = False

    res@cnFillOn                   = False     ; turn on color fill
    res@cnLinesOn                  = True    ; turn off contour lines
    res@cnLineLabelsOn             = True    ; turn off contour line labels
    res@cnInfoLabelOn              = False

    res@gsnScalarContour      = True              ; vectors over scalars
    res@vcRefMagnitudeF       = 20.0              ; vec magnitude
    res@vcRefLengthF          = 0.04              ; size of ref vector
    res@vcLineArrowColor      = "black"           ; white vectors
    res@vcMinDistanceF        = 0.02              ; thins arrows near pole  
    res@vcGlyphStyle          = "CurlyVector"     ; turns on curly vectors

    res@vcRefAnnoArrowLineColor   = "black"       ; black ref vector
    res@vcRefAnnoArrowUseVecColor = False
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
    res@tmYLLabelsOn               = False            ; turn off lat labels
    res@tmXBLabelsOn               = False            ; turn off lon labels
  
    res@tmYROn                     = False
    res@tmYLOn                     = False
    res@tmXTOn                     = False
    res@tmXBOn                     = False

    res@gsnStringFont              = "helvetica-bold"  ; will be used for gsnCenterString later
    res@gsnCenterStringFontHeightF = 0.02              ; ditto
    res@gsnLeftString              = ""                ; remove the left and right substrings
    res@gsnRightString             = ""
    res@gsnCenterString             = ""

   ; Set resources specifically for plots in column 1
    res@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
    res@cnMinLevelValF            = 5100             ; min level
    res@cnMaxLevelValF            = 5950             ; max level
    res@cnLevelSpacingF           = 50.              ; contour interval
    res@cnLineLabelInterval       = 1               ; all contour has the label
    res@cnLineLabelFontHeightF    = 0.02
    res@cnLineLabelPerimOn = False
   ; Set resources for map
    res1 = res
    res1@mpProjection               = "CylindricalEquidistant"       ; projection
    res1@mpMinLatF = 0.0
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
    nrow = 1
    ncol = 3                           
    plot = new(nrow*ncol , "graphic")

    ; Center strings will be set for top plots only
    i = 0
    res1@gsnCenterString = "clim"
    z = fh1->Z({vt},{500},:,:)
    u = fh1->U({vt},{500},:,:)
    v = fh1->V({vt},{500},:,:)
    plot(i) = gsn_csm_vector_scalar_map(wks, u, v, z, res1)
    
    i = i + 1
    res1@gsnCenterString = "est"
    z = fh2->Z({vt},{500},:,:)
    u = fh2->U({vt},{500},:,:)
    v = fh2->V({vt},{500},:,:)
    plot(i) = gsn_csm_vector_scalar_map(wks, u, v, z, res1)
    
    i = i + 1
    res1@gsnCenterString = "mgd"
    z = fh3->Z({vt},{500},:,:)
    u = fh3->U({vt},{500},:,:)
    v = fh3->V({vt},{500},:,:)
    plot(i) = gsn_csm_vector_scalar_map(wks, u, v, z, res1)
    
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
    resP@gsnPanelMainString = date

    gsn_panel(wks, plot, (/nrow,ncol/), resP)  ; plots
     
    frame(wks)
    delete(wks)
end
