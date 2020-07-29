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
  f = addfile("../sens/reg-t.nc","r")
  reg      = f->b(0,:,:,:)
  cor      = f->cor(0,:,:,:)
  thr      = f->thr(0,:,:,:)
  lat      = f->lat
  lon      = f->lon
  plev     = f->lev
  
  wks = gsn_open_wks("png","reg-prtb-tmp_100912")
  plot = new(3,graphic)
  cmap = read_colormap_file("hotcold_18lev")
  
  vdate = 2019100912
  do k = 0,2
     hgt = plev(k)
;---Create plot
     res                = True               ; plot mods desired
     res@gsnDraw	      = False
     res@gsnFrame	      = False
;     res@cnInfoLabelOn	      = False
     
     res@gsnCenterString	   = hgt+"hPa"
     res@gsnCenterStringOrthogonalPosF = 0.02
     res@gsnLeftString	    = ""
     res@gsnAddCyclic	      = False
     res@cnFillOn             = True               ; turn on color for contours
     res@cnLinesOn            = False              ; turn off contour lines
     res@cnLineLabelsOn       = False              ; turn off contour line labels
;     res@gsnScalarContour     = True               ; contours desired
     res@cnFillPalette        = cmap

;     res@mpLandFillColor      = "gray"            ; set land to be gray
     res@mpFillOn	      = False
     res@mpMinLonF            =  120.               ; select a subregion
     res@mpMaxLonF            =  160.
     res@mpMinLatF            =  10. 
     res@mpMaxLatF            =  60.

     res@lbOrientation            = "Vertical"     ; vertical label bar
     res@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer

; note, when doing a subregion, NCL determines the range of the data from
; the full domain. If you wish to just consider the domain you are plotting,
; you must manually set those levels.

      res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
;     res@cnMaxLevelCount	   = 10
      res@cnMinLevelValF       = -1.0               ; set min contour level
      res@cnMaxLevelValF       = 1.0                 ; set max contour level
      res@cnLevelSpacingF      = 0.1               ; set contour spacing

      plotc = gsn_csm_contour_map_ce(wks,cor(k,:,:),res)

      res@cnFillOn             = False             ; turn on color for contours
      res@cnLinesOn            = True              ; turn off contour lines
      res@cnLineLabelsOn       = True
      res@gsnContourZeroLineThicknessF = 0. 	  ; doubles thickness of zero contour
      res@gsnContourNegLineDashPattern = 1 	  ; sets negative contours to dash pattern 1
      res@cnLevelSelectionMode = "ManualLevels"    ; set manual contour levels
      res@cnMinLevelValF       = -1.0              ; set min contour level
      res@cnMaxLevelValF       = 1.0                 ; set max contour level
      res@cnLevelSpacingF      = 0.1               ; set contour spacing
      
; when using the gsn_csm plot templates and overlay, it is necessary to
; turn off the automatic label strings so that you don't get text on 
; top of text
      res@gsnRightString   = ""    
      res@gsnLeftString    = ""
      res@gsnCenterString  = ""
      
      plotb = gsn_csm_contour(wks,reg(k,:,:),res)
      overlay(plotc,plotb)

      res@cnFillOn             = False             ; turn on color for contours
      res@cnLinesOn            = False             ; turn off contour lines
      res@cnLineLabelsOn       = False
      ;res@cnMonoFillPattern    = True 	  ; doubles thickness of zero contour
      ;res@cnFillPattern        = 6 	  ; sets negative contours to dash pattern 1
      ;res@cnMonoFillColor      = True
      ;res@cnFillColor          = "gray"
      ;res@gsnContourZeroLineThicknessF = 3.
      res@cnLevelSelectionMode = "ManualLevels"    ; set manual contour levels
      res@cnMinLevelValF       = -1.0              ; set min contour level
      res@cnMaxLevelValF       = 1.0                 ; set max contour level
      res@cnLevelSpacingF      = 0.1               ; set contour spacing
      res@pmLabelBarDisplayMode = "NoCreate"
      res@pmLegendDisplayMode = "NoCreate"
; when using the gsn_csm plot templates and overlay, it is necessary to
; turn off the automatic label strings so that you don't get text on 
; top of text
      res@gsnRightString   = ""    
      res@gsnLeftString    = ""
      res@gsnCenterString  = ""
      

      plott = gsn_csm_contour(wks,thr(k,:,:),res)
      opt                  = True
      opt@gsnShadeFillType = "Pattern"
      opt@gsnShadeHigh     = 3
      opt@gsnShadeLow      = 3
      opt@gsnShadeFillScaleF = 0.5
      plott = gsn_contour_shade(plott, -0.1, 0.1, opt)
      overlay(plotc,plott)

      plot(k) = plotc
      delete(plotc)
      delete(plotb)
      delete(plott)
      delete(res)
   end do

   resP	                   = True
   resP@gsnPanelMainString	   = "regression of perturbation T against position error "+vdate
   gsn_panel(wks,plot,(/1,3/),resP)

  delete(reg)
  delete(cor)
  delete(thr)
  f = addfile("../sens/reg-slp.nc","r")
  g = addfile("../../netcdf/tigge/2019/jma/100912_mean.nc","r")
;
; Read in Sea Surface Temperature Anomalies
; Read in U and V at 1000 mb [subscript 0]
; Determine the subscripts corresponding to Jan 88
;
  reg      = f->b(0,:,:)
  cor      = f->cor(0,:,:)
  thr      = f->thr(0,:,:)
  lat      = f->lat
  lon      = f->lon
  slp      = g->PRES_meansealevel(0,:,:)
  
  wks = gsn_open_wks("png","reg-prtb-slp_100912")
  cmap = read_colormap_file("hotcold_18lev")
  
  vdate = 2019100912
;---Create plot
  res                = True               ; plot mods desired
  res@gsnDraw	      = False
  res@gsnFrame	      = False
;     res@cnInfoLabelOn	      = False
     
  res@tiMainString	   = "regression of perturbation SLP against position error "+vdate 
  res@tiMainFontHeightF = 0.02
  res@gsnCenterString	   = ""
  res@gsnCenterStringOrthogonalPosF = 0.02
  res@gsnLeftString	    = ""
  res@gsnAddCyclic	      = False
  res@cnFillOn             = True               ; turn on color for contours
  res@cnLinesOn            = False              ; turn off contour lines
  res@cnLineLabelsOn       = False              ; turn off contour line labels
;     res@gsnScalarContour     = True               ; contours desired
  res@cnFillPalette        = cmap

;     res@mpLandFillColor      = "gray"            ; set land to be gray
  res@mpFillOn	      = False
  res@mpMinLonF            =  120.               ; select a subregion
  res@mpMaxLonF            =  160.
  res@mpMinLatF            =  10. 
  res@mpMaxLatF            =  60.

  res@lbOrientation            = "Vertical"     ; vertical label bar
  res@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer

; note, when doing a subregion, NCL determines the range of the data from
; the full domain. If you wish to just consider the domain you are plotting,
; you must manually set those levels.

  res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
;     res@cnMaxLevelCount	   = 10
  res@cnMinLevelValF       = -1.0               ; set min contour level
  res@cnMaxLevelValF       = 1.0                 ; set max contour level
  res@cnLevelSpacingF      = 0.1               ; set contour spacing

  plotc = gsn_csm_contour_map_ce(wks,cor,res)

  res@cnFillOn             = False             ; turn on color for contours
  res@cnLinesOn            = True              ; turn off contour lines
  res@cnLineLabelsOn       = True
  res@gsnContourZeroLineThicknessF = 0. 	  ; doubles thickness of zero contour
  res@gsnContourNegLineDashPattern = 1 	  ; sets negative contours to dash pattern 1
  res@cnLevelSelectionMode = "ManualLevels"    ; set manual contour levels
  res@cnMinLevelValF       = -1.0              ; set min contour level
  res@cnMaxLevelValF       = 1.0                 ; set max contour level
  res@cnLevelSpacingF      = 0.1               ; set contour spacing
      
; when using the gsn_csm plot templates and overlay, it is necessary to
; turn off the automatic label strings so that you don't get text on 
; top of text
  res@gsnRightString   = ""    
  res@gsnLeftString    = ""
  res@gsnCenterString  = ""
      
  plotb = gsn_csm_contour(wks,reg,res)
  overlay(plotc,plotb)

  res@cnFillOn             = False             ; turn on color for contours
  res@cnLinesOn            = False             ; turn off contour lines
  res@cnLineLabelsOn       = False
  ;res@cnMonoFillPattern    = True 	  ; doubles thickness of zero contour
  ;res@cnFillPattern        = 6 	  ; sets negative contours to dash pattern 1
  ;res@cnMonoFillColor      = True
  ;res@cnFillColor          = "gray"
  ;res@gsnContourZeroLineThicknessF = 3.
  res@cnLevelSelectionMode = "ManualLevels"    ; set manual contour levels
  res@cnMinLevelValF       = -1.0              ; set min contour level
  res@cnMaxLevelValF       = 1.0                 ; set max contour level
  res@cnLevelSpacingF      = 0.1               ; set contour spacing
  res@pmLabelBarDisplayMode = "NoCreate"
  res@pmLegendDisplayMode = "NoCreate"
; when using the gsn_csm plot templates and overlay, it is necessary to
; turn off the automatic label strings so that you don't get text on 
; top of text
  res@gsnRightString   = ""    
  res@gsnLeftString    = ""
  res@gsnCenterString  = ""

  plott = gsn_csm_contour(wks,thr,res)
  opt                  = True
  opt@gsnShadeFillType = "Pattern"
  opt@gsnShadeHigh     = 3
  opt@gsnShadeLow      = 3
  opt@gsnShadeFillScaleF = 0.5
  plott = gsn_contour_shade(plott, -0.1, 0.1, opt)
  overlay(plotc,plott)

; note, when doing a subregion, NCL determines the range of the data from
; the full domain. If you wish to just consider the domain you are plotting,
; you must manually set those levels.
  res@cnFillOn             = False               ; turn on color for contours
  res@cnLinesOn            = True              ; turn off contour lines
  res@cnLineLabelsOn       = False              ; turn off contour line labels
  res@cnLineThicknessF     = 2.

  res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
;     res@cnMaxLevelCount	   = 10
  res@cnMinLevelValF       = 96000.0               ; set min contour level
  res@cnMaxLevelValF       = 108000.0                 ; set max contour level
  res@cnLevelSpacingF      = 400.0               ; set contour spacing

  plots = gsn_csm_contour(wks,slp,res)
  overlay(plotc,plots)
  draw(plotc)
  frame(wks)
end
