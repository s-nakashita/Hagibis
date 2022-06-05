;----------------------------------------------------------------------
; vector_1.ncl
;
; Concepts illustrated:
;   - Overlaying vectors and filled contours on a map
;   - Manually attaching lat/lon coordinate arrays to a variable
;   - Changing the length of the smallest vector as a fraction of the reference vector
;   - Moving the vector reference annotation to the top right of the plot
;   - Drawing curly vectors
;   - Setting the color for vectors
;   - Making the labelbar be vertical
;   - Increasing the thickness of vectors
;----------------------------------------------------------------------
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

begin

;---Open netCDF file
  ;orig = "jma"
  dev = "png"
  idate = "2019100900"
  a = addfile("/Volumes/dandelion/netcdf/tigge/2019/"+orig+"/glb_"+idate+"_mean.nc","r")
;
; Read in Sea Surface Temperature Anomalies
; Read in U and V at 1000 mb [subscript 0]
; Determine the subscripts corresponding to Jan 88
;
  u    = a->u(:,:,:,:)
  v    = a->v(:,:,:,:)
  t    = a->t(:,:,:,:)
  lat  = a->latitude
  lon  = a->longitude
  plev = a->level
  lev  = plev
  lev  = lev * 1e2 ;hPa=>Pa
  lev@units="Pa"
  time = a->time
  date = cd_string(time, "%Y%N%D%H")
  lon0 = 120.
  lon1 = 160.
  lat0 = 0.
  lat1 = 40.
  
  out = "pv_"+orig+"_"+date(0)
;  wks = gsn_open_wks("png","prtb_mem26")
;  plot = new(4,graphic)
  cmap = read_colormap_file("hotcold_18lev")
  printVarSummary(cmap)
  newcmap = new((/12,4/), float)
  newcmap(:2,:) =  cmap(0:6:3,:)
  newcmap(3:,:) = cmap(10:,:)

  vdate = 2019100900
  index = ind(date.eq.vdate)
  print(index)
  pv = pot_vort_isobaric(lev,u(index,:,:,:),v(index,:,:,:),t(index,:,:,:),lat,1,0)
  printMinMax(pv,0)
  wks = gsn_open_wks(dev,out+"_"+cd_string(time(index),"%D%H"))
  plots = new(3, graphic)
  iplot = 0
  do k = 2, 4
    hgt = plev(k)
;---Create plot
    res                = True               ; plot mods desired
    res@gsnDraw	      = False
    res@gsnFrame	      = False
;     res@cnInfoLabelOn	      = False
     
    res@gsnCenterString	   = hgt+"hPa"
    res@gsnCenterStringOrthogonalPosF = 0.05
    res@gsnCenterStringParallelPosF = 0.6
    ;res@gsnLeftString	    = "K"
    res@gsnLeftStringFontHeightF = 0.03
    res@gsnRightString	    = "PVU"
    res@gsnRightStringFontHeightF = 0.03
    ;res@gsnRightStringParallelPosF = 1.1    
    res@gsnAddCyclic	      = False
    res@cnFillOn             = True               ; turn on color for contours
    res@cnLinesOn            = True               ; turn on contour lines
    res@cnLineLabelsOn       = False              ; turn off contour line labels
    res@gsnContourNegLineDashPattern = 2
    res@gsnContourZeroLineThicknessF = 3.0
    ;res@gsnScalarContour     = True               ; contours desired

    res@cnFillPalette        = newcmap

;     res@mpLandFillColor      = "gray"            ; set land to be gray
    res@mpFillOn	      = False
    res@mpMinLonF            =  lon0               ; select a subregion
    res@mpMaxLonF            =  lon1
    res@mpMinLatF            =  lat0
    res@mpMaxLatF            =  lat1

    res@lbOrientation            = "Vertical"     ; vertical label bar
    res@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer

    res@vpXF = 0.15
    res@vpWidthF = 0.7

; note, when doing a subregion, NCL determines the range of the data from
; the full domain. If you wish to just consider the domain you are plotting,
; you must manually set those levels.

    res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
;     res@cnMaxLevelCount	   = 10
    res@cnMinLevelValF       = -1.               ; set min contour level
    res@cnMaxLevelValF       = 4.                 ; set max contour level
    res@cnLevelSpacingF      = .5               ; set contour spacing

    var = pv(k,:,:)
    var = var * 1e6
    plots(iplot) = gsn_csm_contour_map_ce(wks,var,res)
    iplot = iplot + 1
   end do

   resP	                   = True
   ;resP@gsnPanelMainString	   = "perturbation member 26 T&Q&wind,slp "+vdate
   gsn_panel(wks,plots,(/1,3/),resP)
end
