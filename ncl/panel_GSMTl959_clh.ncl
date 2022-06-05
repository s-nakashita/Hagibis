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
   sst = ""
   fh1     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100"+sst+"/fcst_phy2m_snp_"+yyyymmddhh+".nc","r")
   fh2     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100"+sst+"/fcst_p_asia_"+yyyymmddhh+".nc","r")
   ;fh3     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_mgdsst/fcst_surf_2019100900.nc", "r")
   
   ctop = fh1->CTOP(:, ::-1, :) ; cloud top pressure [hPa]
   cbase= fh1->CBASE(:, ::-1, :) ; cloud bottom pressure [hPa]
   plev = fh2->level         ; pressure levels [hPa]
   z    = fh2->Z(:, :, ::-1, :) ; geopotential heights [m]
   phi  = fh2->PHI(0, :, :) ; surface geopotential [m^2/s^2]
   zsfc = phi / 9.8065
   t6h  = fh1->time
   t3h  = fh2->time
   print(plev)

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
   pltName= "GSMtl959"+sst+"_clth_"+yyyymmddhh                                   
   pltType= "png"               ; ps, eps, png, pdf, x11 .....
   d = ind(date.eq."20191009 00:00")
   print(d)
   pltPath = pltDir+pltName+"_"+cd_string(t6h(d), "%D%H")
;****************************************************************************
; Horizontal interpolation
;****************************************************************************
   ctop_ll = linint2_Wrap(longs, latgs, ctop(d,:,:), True, lonll, latll, 0)
   ctop_ll!0    = "lat"
   ctop_ll!1    = "lon"
   ctop_ll&lat  = latll
   ctop_ll&lon  = lonll
   ctop_ll@long_name = "cloud top pressure"
   ctop_ll@units     = "hPa"
   printMinMax(ctop_ll, 0)
   cbase_ll = linint2_Wrap(longs, latgs, cbase(d,:,:), True, lonll, latll, 0)
   cbase_ll!0    = "lat"
   cbase_ll!1    = "lon"
   cbase_ll&lat  = latll
   cbase_ll&lon  = lonll
   cbase_ll@long_name = "cloud bottom pressure"
   cbase_ll@units     = "hPa"
   printMinMax(ctop_ll, 0)
  delete(ctop)
  delete(cbase)
;****************************************************************************
; Convert pressure to height
;****************************************************************************
  clth = new(dimsizes(ctop_ll), typeof(z))
  clbh = new(dimsizes(cbase_ll), typeof(z))
  tmp_z    = z(ind(t3h.eq.t6h(d)),:,:,:)
  printMinMax(tmp_z,0)
  do j = 0, dimsizes(latll)-1
    do i = 0, dimsizes(lonll)-1
      if ismissing(ctop_ll(j,i)) then
        clth(j,i) = 0.0
        continue
      end if
      pj = ctop_ll(j,i)
      do k = 0, dimsizes(plev)-1
        if (pj.gt.plev(k)) then
          break
        end if
      end do
      if (k.gt.0) then
        wgt = (log(pj) - log(plev(k-1)))/(log(plev(k)) - log(plev(k-1)))
        clth(j,i) = ((1.0 - wgt)*tmp_z(k-1,j,i) + wgt*tmp_z(k,j,i) - zsfc(j,i)) * 0.001 ; m->km
      else
        clth(j,i) = 0.0
      end if
    end do
  end do
  do j = 0, dimsizes(latll)-1
    do i = 0, dimsizes(lonll)-1
      if ismissing(cbase_ll(j,i)) then
        clbh(j,i) = 0.0
        continue
      end if
      pj = cbase_ll(j,i)
      do k = 0, dimsizes(plev)-1
        if (pj.gt.plev(k)) then
          break
        end if
      end do
      if (k.gt.0) then
        wgt = (log(pj) - log(plev(k-1)))/(log(plev(k)) - log(plev(k-1)))
        clbh(j,i) = ((1.0 - wgt)*tmp_z(k-1,j,i) + wgt*tmp_z(k,j,i) - zsfc(j,i)) * 0.001 ; m->km
      else
        clbh(j,i) = 0.0
      end if
    end do
  end do
  delete(tmp_z)
  clth!0 = "lat"
  clth!1 = "lon"
  clth&lat  = ctop_ll&lat
  clth&lon  = ctop_ll&lon
  clth@long_name = "cloud top height"
  clth@units     = "km"
  printVarSummary(clth)
  printMinMax(clth, 0)
  clbh!0 = "lat"
  clbh!1 = "lon"
  clbh&lat  = cbase_ll&lat
  clbh&lon  = cbase_ll&lon
  clbh@long_name = "cloud bottom height"
  clbh@units     = "km"
  printVarSummary(clbh)
  printMinMax(clbh, 0)
  delete(z)
;  exit()
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
  res1@cnFillPalette              = cmap
  res1@gsnCenterString           = clth@long_name
  res1@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
  res1@cnMinLevelValF            = 0.            ; min level
  res1@cnMaxLevelValF            = 18.          ; max level
  res1@cnLevelSpacingF           = 1.            ; contour interval
  plot(0) = gsn_csm_contour_map_ce(wks, clth, res1)
  ; Set resources specifically for plots in column 1
  res1@cnFillPalette              = cmap(::-1,:)
  res1@gsnCenterString           = ctop_ll@long_name
  res1@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
  res1@cnMinLevelValF            = 100.            ; min level
  res1@cnMaxLevelValF            = 1000.           ; max level
  res1@cnLevelSpacingF           = 50.            ; contour interval
  plot(1) = gsn_csm_contour_map_ce(wks, ctop_ll, res1)
  ; Set resources specifically for plots in column 1
  res1@cnFillPalette              = cmap
  res1@gsnCenterString           = clbh@long_name
  res1@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
  res1@cnMinLevelValF            = 0.            ; min level
  res1@cnMaxLevelValF            = 18.           ; max level
  res1@cnLevelSpacingF           = 1.            ; contour interval
  plot(2) = gsn_csm_contour_map_ce(wks, clbh, res1)
  ; Set resources specifically for plots in column 1
  res1@cnFillPalette              = cmap(::-1,:)
  res1@gsnCenterString           = cbase_ll@long_name
  res1@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
  res1@cnMinLevelValF            = 100.            ; min level
  res1@cnMaxLevelValF            = 1000.           ; max level
  res1@cnLevelSpacingF           = 50.            ; contour interval
  plot(3) = gsn_csm_contour_map_ce(wks, cbase_ll, res1)
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
