;*************************************************
; NCL tutorial script: vort_1.ncl
;*************************************************
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "write_nc.ncl"
begin
;************************************************
; variable and file handling
;************************************************
  in  = addfile("../../netcdf/tigge/2019/jma/glb_100912_mean.nc","r")  ; open netcdf file
  u   = in->UGRD;(:,:,180:360,180:360)                       ; pull u off file
  v   = in->VGRD;(:,:,180:360,180:360)                       ; pull v off file
  t0  = in->time
  lev0 = in->level
  lat0 = in->lat
  lon0 = in->lon
  print(t0)
  printVarSummary(u)
  ofile = "../../netcdf/tigge/2019/jma/glb_100912_mean_v.nc"
;  setfileoption("nc", "format", "netcdf4")
;  setfileoption("nc", "headerReserveSpace", 64000)
;  setfileoption("nc", "preFill", False)
;  if (isfilepresent(ofile)) then
;    system("rm " + ofile)
;  end if
;  out = addfile(ofile,"c") 
;************************************************
; calculate vorticity on a Gaussian Grid
; scale for visual convenience
;************************************************
  scale = 1.e05                                  
  vrt   = u                                      ; retain coordinates
  ;vrt   = uv2dv_cfd(v, -u, lat0, lon0, 2)
  vrt = uv2vrF_Wrap(u,v) ;* scale
  vrt@long_name = "Vorticity"
  vrt@units = "/s"
  printVarSummary(vrt)
;************************************************
; create netcdf4 file
;************************************************
  write_nc(ofile,vrt,t0,lev0,lat0,lon0)

;************************************************
; create plot
;************************************************
  wks  = gsn_open_wks("png","vort")              ; send graphics to PNG file

  res                   = True                   ; plot mods desired
  res@cnFillOn          = True                   ; turn on color
  res@cnFillPalette     = "BlWhRe"               ; set color map
  res@cnLinesOn         = False                  ; no contour lines

  res@gsnAddCyclic	      = False 
  res@mpLandFillColor      = "gray"            ; set land to be gray
  ;res@mpMinLonF            =  120.               ; select a subregion
  ;res@mpMaxLonF            =  160.
  ;res@mpMinLatF            =  10. 
  ;res@mpMaxLatF            =  60.

  res@cnLevelSelectionMode = "ManualLevels"       ; set manual contour levels
  res@cnMinLevelValF       = -4                   ; set min contour level
  res@cnMaxLevelValF       =  4                   ; set max contour level
  res@cnLevelSpacingF      =  0.5                 ; set contour spacing

  vrt = vrt*scale
  plot = gsn_csm_contour_map(wks,vrt(0,4,:,:),res) ; create plot  
;************************************************
end
