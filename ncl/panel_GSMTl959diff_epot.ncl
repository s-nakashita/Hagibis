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
   fhp1    = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_est/fcst_surf_asia_"+yyyymmddhh+".nc", "r")
   fh1     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_est/fcst_phy2m_snp_"+yyyymmddhh+".nc", "r")
   fhp2    = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100/fcst_surf_asia_"+yyyymmddhh+".nc", "r")
   fh2     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100/fcst_phy2m_snp_"+yyyymmddhh+".nc","r")
   fhp3    = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_mgdsst/fcst_surf_asia_"+yyyymmddhh+".nc", "r")
   fh3     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_mgdsst/fcst_phy2m_snp_"+yyyymmddhh+".nc", "r")
   
   ps_1 = fhp1->PSEA(:, :, :)
   ps_2 = fhp2->PSEA(:, :, :)
   ps_3 = fhp3->PSEA(:, :, :)
   t1h  = fhp1->time
   t6h  = fh1->time

   date = cd_string(t6h, "%N%D %H:%M")
   ntim = dimsizes(date)
   
   pltDir = "./"                ; where output will be written
   pltName= "GSMtl959diff_epot_"+yyyymmddhh                                   
   pltType= "png"               ; ps, eps, png, pdf, x11 .....
;****************************************************************************
; Calculate equivalent potential temperature
;****************************************************************************
  q    = fh1->Q(:,:,:)
  p    = fh1->P(:,:,:)
  temp = fh1->T(:,:,:)
  w    = mixhum_convert(q, "w", (/0,0/))
  tlcl = tlcl_mixr_bolton(temp, w, p, (/0,0,1,1/))
  ept1 = pot_temp_equiv_tlcl(p, temp, tlcl, w, (/0,1,0,1/))
  printMinMax(ept1,0)
  delete(q)
  delete(p)
  delete(temp)
  delete(w)
  delete(tlcl)
  q    = fh2->Q(:,:,:)
  p    = fh2->P(:,:,:)
  temp = fh2->T(:,:,:)
  w    = mixhum_convert(q, "w", (/0,0/))
  tlcl = tlcl_mixr_bolton(temp, w, p, (/0,0,1,1/))
  ept2 = pot_temp_equiv_tlcl(p, temp, tlcl, w, (/0,1,0,1/))
  printMinMax(ept2,0)
  delete(q)
  delete(p)
  delete(temp)
  delete(w)
  delete(tlcl)
  q    = fh3->Q(:,:,:)
  p    = fh3->P(:,:,:)
  temp = fh3->T(:,:,:)
  w    = mixhum_convert(q, "w", (/0,0/))
  tlcl = tlcl_mixr_bolton(temp, w, p, (/0,0,1,1/))
  ept3 = pot_temp_equiv_tlcl(p, temp, tlcl, w, (/0,1,0,1/))
  printMinMax(ept3,0)
  delete(q)
  delete(p)
  delete(temp)
  delete(w)
  delete(tlcl)
;****************************************************************************
; Graphic resources
;****************************************************************************
    ;ver = 4
  nrow = 4
  ncol = 3
  do ds = 0, 12, 4
  plot = new(nrow*ncol , "graphic")
  ;ds = 12
  de = min((/ds+nrow-1,ntim-1/))
  i = 0
  pltPath = pltDir+pltName+"_"+cd_string(t6h(ds), "%D%H")+"-"+cd_string(t6h(de), "%D%H")
  wks     = gsn_open_wks(pltType, pltPath) 
  do d=ds, de
    res1                            = True     ; plot mods desired
    res1@gsnDraw                    = False    ; will panel plots later
    res1@gsnFrame                   = False

    res1@cnFillOn                   = True     ; turn on color fill
    res1@cnLinesOn                  = False    ; turn off contour lines
    res1@cnLineLabelsOn             = False    ; turn off contour line labels

;**********************************************************************
; We are turning the labelbar off for all but the bottom row of plots
;**********************************************************************
    res1@lbLabelBarOn               = False
    res1@lbTitleFontHeightF         = .03                 ; make title smaller
    res1@lbTitleString              = "K"            ; title string
    res1@lbTitlePosition            = "Bottom"            ; title position
    res1@lbTitleDirection           = "Across"            ; title direction
    ;res1@lbLabelStride              = 2
    res1@lbLabelFontHeightF         = 0.025
    res1@pmLabelBarHeightF          = 0.1
    res1@lbBoxEndCapStyle           = "TriangleBothEnds"
    ; Set resources for map
    res1@mpMinLatF                  = 5.
    res1@mpMaxLatF                  = 40.
    res1@mpMinLonF                  = 120.
    res1@mpMaxLonF                  = 170.

    res1@mpDataBaseVersion          = "MediumRes"
    res1@mpOutlineOn                = True
    res1@mpFillOn                   = True                ; turn off gray fill
    ;res@mpOutlineBoundarySets      = "National"

    res1@gsnAddCyclic               = False            ; regional data, don't add cyclic point
 
    ; Set resources for tickmarks
    res1@gsnStringFont              = "helvetica-bold"  ; will be used for gsnCenterString later
    res1@gsnCenterStringFontHeightF = 0.04              ; ditto
    res1@gsnLeftString              = ""                ; remove the left and right substrings
    res1@gsnRightString             = date(d)
    
   ; Set resources specifically for plots in column 1
    res1@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
    
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
    ; first column : raw values
    ind_p = ind(t1h.eq.t6h(d))
    res1@cnMinLevelValF            = 280.           ; min level
    res1@cnMaxLevelValF            = 310.           ; max level
    res1@cnLevelSpacingF           = 2.            ; contour interval
    res1@cnFillPalette             = "temp_19lev"
    if (d.eq.ds) then
      res1@gsnCenterString            = "EST"
    else
      res1@gsnCenterString            = ""
    end if
    if (d.eq.de) then
      res1@lbLabelBarOn = True
    end if
    plot(i) = gsn_csm_contour_map_ce(wks, ept1(d,:,:), res1)
    plot2 = gsn_csm_contour(wks, ps_1(ind_p,:,:), res2)
    overlay(plot(i),plot2)
    delete(plot2)
    i = i + 1
    eptd = ept2(d,:,:)
    eptd = eptd - ept1(d,:,:)
    res1@cnMinLevelValF            = -2.            ; min level
    res1@cnMaxLevelValF            = 2.           ; max level
    res1@cnLevelSpacingF           = .4            ; contour interval
    res1@cnFillPalette              = "CBR_drywet"
    res1@lbLabelBarOn = False
    if (d.eq.ds) then
      res1@gsnCenterString            = "CLIM-EST"
    else
      res1@gsnCenterString            = ""
    end if
    plot(i) = gsn_csm_contour_map_ce(wks, eptd, res1)
    plot2 = gsn_csm_contour(wks, ps_2(ind_p,:,:), res2)
    overlay(plot(i),plot2)
    delete(plot2)
    i = i + 1
    eptd = ept3(d,:,:)
    eptd = eptd - ept1(d,:,:)
    if (d.eq.ds) then
      res1@gsnCenterString            = "MGD-EST"
    else
      res1@gsnCenterString            = ""
    end if
    if (d.eq.de) then
      res1@lbLabelBarOn = True
      res1@pmLabelBarWidthF          = 1.05
      res1@pmLabelBarParallelPosF     = -0.02
    end if
    plot(i) = gsn_csm_contour_map_ce(wks, eptd, res1)
    plot2 = gsn_csm_contour(wks, ps_3(ind_p,:,:), res2)
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
    resP@gsnPanelLabelBar = False

    gsn_panel(wks, plot, (/nrow,ncol/), resP)  ; plots
    frame(wks)
    delete(wks)
  end do
end
