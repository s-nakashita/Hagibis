;*************************************************
; NCL tutorial script: vort_1.ncl
;*************************************************
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
begin
;************************************************
; variable and file handling
;************************************************
  in  = addfile("../../netcdf/gsm/gl/2019/10/init_vort.nc","r")                             ; open netcdf file
  vrt = in->vort                                         ; pull u off file
;************************************************
; calculate vorticity on a Gaussian Grid
; scale for visual convenience
;************************************************
  scale = 1.e05                   
  vrt   = vrt * scale               
;  vrt   = u                                      ; retain coordinates
;  vrt   = uv2vrG_Wrap(u,v) * scale
;  vrt@long_name = "vorticity"
;  vrt@units     = "scaled"
;************************************************
; create plot
;************************************************
  wks  = gsn_open_wks("png","vort_2")              ; send graphics to PNG file

  res                   = True                   ; plot mods desired
  res@cnFillOn          = True                   ; turn on color
  res@cnFillPalette     = "BlWhRe"               ; set color map
  res@cnLinesOn         = False                  ; no contour lines

  res@gsnAddCyclic	      = False 
  res@mpLandFillColor      = "gray"            ; set land to be gray
  res@mpMinLonF            =  120.               ; select a subregion
  res@mpMaxLonF            =  160.
  res@mpMinLatF            =  10. 
  res@mpMaxLatF            =  60.

  res@cnLevelSelectionMode = "ManualLevels"       ; set manual contour levels
  res@cnMinLevelValF       = -4                   ; set min contour level
  res@cnMaxLevelValF       =  4                   ; set max contour level
  res@cnLevelSpacingF      =  0.5                 ; set contour spacing


  plot = gsn_csm_contour_map(wks,vrt(0,4,:,:),res) ; create plot  
;************************************************
end
