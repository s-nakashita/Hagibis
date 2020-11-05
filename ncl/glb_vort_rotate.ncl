;*************************************************
; NCL tutorial script: vort_1.ncl
;*************************************************
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "../rotate/librotate.ncl"
load "write_nc.ncl"

procedure rotate_pole(u, v, lon, ur, vr, opt)
; u,v dimension is (time,lev,{89.5:90},lon) or (time,lev,{-90:-89.5},lon)
; North pole; x-axis=>90E, y-axis=>180  South pole; x-axis=>90E, y-axis=>0
local coslon, sinlon
begin
  coslon = cos(lon*deg2rad)
  sinlon = sin(lon*deg2rad)
  if(opt.eq.0) then
  ;North pole
    ur =  u * conform_dims(dimsizes(u), coslon, 3) + v * conform_dims(dimsizes(v), sinlon, 3)
    vr = -u * conform_dims(dimsizes(u), sinlon, 3) + v * conform_dims(dimsizes(v), coslon, 3)
  else if(opt.eq.1) then
  ;South pole
    ur = u * conform_dims(dimsizes(u), coslon, 3) - v * conform_dims(dimsizes(v), sinlon, 3)
    vr = u * conform_dims(dimsizes(u), sinlon, 3) + v * conform_dims(dimsizes(v), coslon, 3)
  end if
  end if
end

begin
;************************************************
; variable and file handling
;************************************************
  init = "2019100900"
  yyyy = str_get_cols(init, 0, 3)
  mmddhh = str_get_cols(init, 4, 9)
  in  = addfile("../../netcdf/tigge/"+yyyy+"/jma/glb_"+mmddhh+"_mean.nc","r")  ; open netcdf file
  ;fu  = addfile("../../netcdf/nc-reanl/2019/uwnd.2019.nc","r")  ; open netcdf file
  ;fv  = addfile("../../netcdf/nc-reanl/2019/vwnd.2019.nc","r")  ; open netcdf file
  u   = in->UGRD;(:,:,180:360,180:360)                       ; pull u off file
  v   = in->VGRD;(:,:,180:360,180:360)                       ; pull v off file
  t0  = in->time
  lev0 = in->level
  lat0 = dble2flt(in->lat)
  lon0 = dble2flt(in->lon)
  ;u   = fu->uwnd(0,:,:,:)                       ; pull u off file
  ;v   = fv->vwnd(0,:,:,:)                       ; pull v off file
  ;t = fu->time
  ;lev0 = fu->level
  ;lat0 = fu->lat ;N->S
  ;lon0 = fu->lon
  ;print(lon0)
  uout = "../../netcdf/tigge/"+yyyy+"/jma/glb_"+mmddhh+"_mean_u.nc"
  vout = "../../netcdf/tigge/"+yyyy+"/jma/glb_"+mmddhh+"_mean_vg.nc"
  ofile = "../../netcdf/tigge/"+yyyy+"/jma/rotate/np_glb_"+mmddhh+"_mean_v.nc"
;  setfileoption("nc", "format", "netcdf4")
;  setfileoption("nc", "headerReserveSpace", 64000)
;  setfileoption("nc", "preFill", False)
  if (isfilepresent(uout)) then
    system("rm " + uout)
  end if
  if (isfilepresent(vout)) then
    system("rm " + vout)
  end if
  if (isfilepresent(ofile)) then
    system("rm " + ofile)
  end if
;  out = addfile(ofile,"c") 
;************************************************
; rotate pole vector
;************************************************
  print(dim_avg_n(u(0,0,{89.5},:),0))
  print(dim_avg_n(v(0,0,{89.5},:),0))
  print(dim_avg_n(u(0,0,{90},:),0))
  print(dim_avg_n(v(0,0,{90},:),0))
  unp = u(:,:,{89.5:90},:)
  vnp = v(:,:,{89.5:90},:)
  ur = new(dimsizes(unp),typeof(unp))
  vr = new(dimsizes(vnp),typeof(vnp))
  rotate_pole(unp, vnp, lon0, ur, vr, 0)
  print(dim_avg_n(ur(0,0,0,:),0))
  u(:,:,{89.5:90},:) = (/ur/)
  v(:,:,{89.5:90},:) = (/vr/)
  usp = u(:,:,{-90:-89.5},:)
  vsp = v(:,:,{-90:-89.5},:)
  rotate_pole(usp, vsp, lon0, ur, vr, 1)
  print(dim_avg_n(ur(0,0,0,:),0))
  u(:,:,{-90:-89.5},:) = (/ur/)
  v(:,:,{-90:-89.5},:) = (/vr/)
  ;print(u(0,0,{-90},360))
  ;print(v(0,0,{-90},360))
  write_nc(uout, u, t0, lev0, lat0, lon0)
  write_nc(vout, v, t0, lev0, lat0, lon0)
  ;exit
;************************************************
; calculate vorticity on a fixed Grid
; scale for visual convenience
;************************************************
  scale = 1.e05                                  
  vrt   = u                                      ; retain coordinates
  vrt = uv2vrF_Wrap(u, v) ;* scale
  vrt@long_name = "Vorticity"
  vrt@units = "/s"
  if(any(ismissing(vrt)))
      print("Missing values in vorticity")
      exit
  end if  
  printVarSummary(vrt)
;************************************************
; rotate coordinate
;************************************************
  ndim = dimsizes(vrt)
  ntim = ndim(0)
  nlev = ndim(1)
  nlat = ndim(2)
  nlon = ndim(3)
  lonout = new(nlon * nlat, "float")
  latout = new(nlon * nlat, "float")
  trackfile = "../pytrack/jma/track"+init+"_mean.txt"
  track = asciiread(trackfile, (/21, 7/), "float")
  vrt_r = vrt
  do l = 0, 20
    dd = sprinti("%0.2d",toint(track(l, 2)))
    hh = sprinti("%0.2d",toint(track(l, 3)))
    print(dd+hh)
    lonc = track(l, 4) * deg2rad
    latc = track(l, 5) * deg2rad
    print("lonc = "+sprintf("%5.3f",lonc*rad2deg)+", latc = "+sprintf("%5.3f",latc*rad2deg))
    lonin = lon0 * deg2rad
    latin = lat0 * deg2rad
    rotate_lonlat(lonc, latc, lonin, latin, lonout, latout)
    printMinMax(lonout, 0)
    printMinMax(latout, 0)
    datain = vrt(l,:,:,:)
    dataout = linint2_points(lon0, lat0, datain(:,:,:), True, lonout*rad2deg, latout*rad2deg, 0)
    if(any(ismissing(dataout(:,:))))
      print("Missing values in interpolation")
      ;lmiss = new(nlev, "logical")
      do k = 0, nlev-1
        nsize = dimsizes(lonout)
        do ij = 0, nsize(0)-1
          if(ismissing(dataout(k,ij)))
            ;lmiss(k) = True
            print(sprinti("%d",ij))
            print("level "+sprintf("%4.1f",lev0(k))+", lon = "+sprintf("%5.3f",lonout(ij)*rad2deg)+", lat = "+sprintf("%5.3f",latout(ij)*rad2deg))
          end if
        end do
      end do
      datamsg = linmsg(dataout(:,:),0)
      dataout = (/datamsg/)
      delete(datamsg)
      ;exit
    end if
    ;printVarSummary(dataout)
    do k = 0, nlev-1
      ;ij = 0
      ;do j = 0, nlat-1
      ;  do i = 0, nlon-1
      vrt_r(l,k,:,:) = (/reshape(dataout(k,:),(/nlat,nlon/))/)
      ;    ij = ij + 1
      ;  end do
      ;end do
    end do
  end do

;************************************************
; create netcdf4 file
;************************************************
  write_nc(ofile,vrt_r,t0,lev0,lat0,lon0)
  ;exit()

;************************************************
; create plot
;************************************************
  wks  = gsn_open_wks("png","vort")              ; send graphics to PNG file

  res                   = True                   ; plot mods desired
  res@gsnPolar          = "SH"
  ;res@mpMinLatF         = 80
  res@mpMaxLatF         = -80
  res@cnFillOn          = True                   ; turn on color
  res@cnFillPalette     = "BlWhRe"               ; set color map
  res@cnLinesOn         = False                  ; no contour lines

  ;res@gsnAddCyclic	      = False 
  res@mpLandFillColor      = "gray"            ; set land to be gray
  ;res@mpMinLonF            =  120.               ; select a subregion
  ;res@mpMaxLonF            =  160.
  ;res@mpMinLatF            =  10. 
  ;res@mpMaxLatF            =  60.

  res@cnLevelSelectionMode = "ManualLevels"       ; set manual contour levels
  res@cnMinLevelValF       = -20                   ; set min contour level
  res@cnMaxLevelValF       =  20                   ; set max contour level
  res@cnLevelSpacingF      =  2.                 ; set contour spacing

  res@gsnScalarContour      = True              ; vectors over scalars
  res@vcRefMagnitudeF       = 20.0              ; vec magnitude
  res@vcRefLengthF          = 0.05              ; size of ref vector
  res@vcLineArrowColor      = "black"           ; white vectors
  res@vcMinDistanceF        = 0.01              ; thins arrows near pole  
  res@vcGlyphStyle          = "LineArrow"     ; turns on curly vectors
  res@vcLineArrowThicknessF = 2.
  res@vcMinMagnitudeF       = 5.0

  vrt = vrt*scale
  ;plot = gsn_csm_contour_map_polar(wks,vrt(3,:,:),res)
  plot = gsn_csm_vector_scalar_map_polar(wks,u(0,3,:,:),v(0,3,:,:),vrt(0,3,:,:),res)
  ;plot = gsn_csm_contour_map(wks,vrt(0,4,:,:),res) ; create plot  
;************************************************
end
