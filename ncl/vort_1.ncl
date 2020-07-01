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
  in  = addfile("../../netcdf/gsm/gl/2019/10/init.nc","r")  ; open netcdf file
  u   = in->UGRD(:,:,180:360,180:360)                       ; pull u off file
  v   = in->VGRD(:,:,180:360,180:360)                       ; pull v off file
  t0  = in->time
  lev0 = in->level
  lat0 = in->lat
  lon0 = in->lon
  printVarSummary(u)
  ofile = "../../netcdf/gsm/gl/2019/10/init_vort.nc"
;  setfileoption("nc", "format", "netcdf4")
;  setfileoption("nc", "headerReserveSpace", 64000)
;  setfileoption("nc", "preFill", False)
  if (isfilepresent(ofile)) then
    system("rm " + ofile)
  end if
  out = addfile(ofile,"c") 
;************************************************
; calculate vorticity on a Gaussian Grid
; scale for visual convenience
;************************************************
  scale = 1.e05                                  
  vrt   = u                                      ; retain coordinates
  vrt   = uv2vrF_Wrap(u,v) ;* scale
  printVarSummary(vrt)
;************************************************
; create netcdf4 file
;************************************************
  ndim = dimsizes(vrt)
  ctim = ndim(0)
  nlev = ndim(1)
  nlat = ndim(2)
  nlon = ndim(3)
  print(ctim)
  print(nlev)
  print(nlat)
  print(nlon)
  dim_names = (/"time","level","lat","lon"/)
  dim_sizes = (/ctim,nlev,nlat,nlon/)
  dim_unlimited = (/True,False,False,False/)
  filedimdef(out,dim_names,dim_sizes,dim_unlimited)

  atts = True
  atts@description = "compute Vorticity"
  atts@timestamp = systemfunc("date")
  print("Defining file attributes")
  fileattdef(out,atts)

  print("Defining file variables")
  filevardef(out,"time","double","time")
  filevardef(out,"level","double","level")
  filevardef(out,"lat","double","lat")
  filevardef(out,"lon","double","lon")
  filevardef(out,"vort","float",(/"time","level","lat","lon"/))

  print("Defining file variable attributes")
  tatts = True
  tatts@standard_name = t0@standard_name
  tatts@units = t0@units
  filevarattdef(out,"time",tatts)

  levatts = True
  levatts@long_name = lev0@long_name
  levatts@units = lev0@units
  filevarattdef(out,"level",levatts)

  latatts = True
  latatts@long_name = lat0@long_name
  latatts@units = lat0@units
  filevarattdef(out,"lat",latatts)

  lonatts = True
  lonatts@long_name = lon0@long_name
  lonatts@units = lon0@units
  filevarattdef(out,"lon",lonatts)

  vatts = 0.0
  vatts@long_name = "vorticity"
  vatts@units     = "/s"
  vatts@_FillValue = 1e32
  filevarattdef(out,"vort",vatts)

  time  = t0
  level = lev0
  lat   = lat0(180:360)
  lon   = lon0(180:360)
  out->level = (/level/)
  out->lat   = (/lat/)
  out->lon   = (/lon/)

  do n = 0, ctim-1
    out->time(n) = (/time(n)/)
    out->vort(n,:,:,:) = (/vrt(n,:,:,:)/)
  end do
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
  res@mpMinLonF            =  120.               ; select a subregion
  res@mpMaxLonF            =  160.
  res@mpMinLatF            =  10. 
  res@mpMaxLatF            =  60.

  res@cnLevelSelectionMode = "ManualLevels"       ; set manual contour levels
  res@cnMinLevelValF       = -4                   ; set min contour level
  res@cnMaxLevelValF       =  4                   ; set max contour level
  res@cnLevelSpacingF      =  0.5                 ; set contour spacing

  vrt = vrt*scale
  plot = gsn_csm_contour_map(wks,vrt(0,4,:,:),res) ; create plot  
;************************************************
end
