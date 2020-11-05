; ***********************************************
; barb_4.ncl
;
; Concepts illustrated:
;   - Drawing wind barbs
;   - Drawing a wind profile plot
;   - Reversing the Y axis in a wind profile plot
;   - Changing the thickness of a wind barb tick
;   - Turning off the vector reference annotation label
;   - Explicitly setting wind barb levels
;   - Changing the length of the wind barbs
;   - Making the labelbar be vertical
;   - Linearizing the Y axis
;   - Creating a LogLin object
;   - Retrieving the axes mins and maxs of a plot
;
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "../rotate/librotate.ncl"
function rotate_uv(lonc, latc, lonnp, latnp, unp, vnp)
local lonout, latout, xd, yd, zd, lon3d, lat3d, xdtc, ydtc, zdtc, lonc2d, latc2d, lontc, lattc, utc, vtc
begin
  Olist = NewList("lifo")
  ndim = dimsizes(unp)
  nlon = ndim(1)
  nlat = ndim(2)
  lonout = new(nlon*nlat, "float")
  latout = new(nlon*nlat, "float") 
  inv_rotate_lonlat(lonc,latc,lonnp,latnp,lonout,latout)
  ;(u, v) => (xdot, ydot, zdot)
  xd = unp
  yd = unp
  zd = unp
  lon3d = conform(unp, lonnp, 2)
  lat3d = conform(unp, latnp, 1)
  uv2xyzd(unp, vnp, lon3d, lat3d, xd, yd, zd)
  ; bilinear interpolation of winds
  xdnp = linint2_points(lon, lat(::-1), xd(:,::-1,:), True, lonout, latout, 0)
  ydnp = linint2_points(lon, lat(::-1), yd(:,::-1,:), True, lonout, latout, 0)
  zdnp = linint2_points(lon, lat(::-1), zd(:,::-1,:), True, lonout, latout, 0)
  if(any(ismissing(xdnp)).and.any(ismissing(ydnp)).and.any(ismissing(zdnp)))
    print("Missing values in interpolation")
    exit
  end if
  ; rotate coordinates to TC grid
  xdtc = xd
  ydtc = yd
  zdtc = zd
  lonc2d = conform(xd, lonc, -1)
  latc2d = conform(xd, latc, -1)
  np2tc(lonc2d, latc2d, xd, yd, zd, xdtc, ydtc, zdtc)
  if(any(ismissing(xdtc)).and.any(ismissing(ydtc)).and.any(ismissing(zdtc)))
    print("Missing values in rotation")
    exit
  end if
  ; (xdot, ydot, zdot) => (u, v)
  lontc = lonout
  lattc = latout
  utc = xdtc
  vtc = xdtc
  xyzd2uv(xdtc, ydtc, zdtc, lonout, utc, vtc)
  if(any(ismissing(utc)).and.any(ismissing(vtc)))
    print("Missing values in tc")
    exit
  end if
  ListPush(Olist,vtc)
  ListPush(Olist,utc)
  ListPush(Olist,lattc)
  ListPush(Olist,lontc)
  return(Olist)
end

begin
  year = "2019"
  init = "2019100912"
  a = addfile("../../netcdf/tigge/"+year+"/jma/rotate/np_"+init+"_mean.nc","r")
  time = a->time
  date = cd_string(time, "%Y%N%D%H")
  sdate = 2019100912
  edate = 2019101212
  ds = ind(date.eq.sdate)
  de = ind(date.eq.edate)
  unp = a->UGRD(ds:de,:,:,:)			   ; (time, lev)
  vnp = a->VGRD(ds:de,:,:,:)		     ; (time, lev)
  lonnp = dble2flt(a->lon)
  latnp = dble2flt(a->lat)
  lev  = a->level
  t = time(ds:de)
  ;lonc = 120.
  ;latc = 30.
  print(typeof(lonnp))
  ;Olist = rotate_uv(lonc,latc,lonnp,latnp,unp(0,:,:,:),vnp(0,:,:,:))
  ;lontc = Olist[0]
  ;print(lontc)

  ; (lon,lat) => (x,y)
  pi = atan(1.0) * 4.0
  deg2rad = pi / 180.0
  u = unp
  v = vnp
  u = unp * conform_dims(dimsizes(unp), cos(lonnp*deg2rad), 3) - vnp * conform_dims(dimsizes(vnp), sin(lonnp*deg2rad), 3)
  v = unp * conform_dims(dimsizes(unp), sin(lonnp*deg2rad), 3) + vnp * conform_dims(dimsizes(vnp), cos(lonnp*deg2rad), 3)

  ; annulus average
  maxlatF = 87.
  print(latnp({maxlatF:}))
  uav = dim_avg_n(u(:,:,{maxlatF:},:), (/2,3/))
  uav!0 = "time"
  uav!1 = "lev"
  uav&time = (t - min(time)) / 3600
  uav&lev = lev
  vav = uav
  vav = dim_avg_n(v(:,:,{maxlatF:},:), (/2,3/))
  
  wspd = uav
  wspd = (/ sqrt(uav^2+vav^2) /)
  printMinMax(wspd, 0)
;===========================================================
  wks = gsn_open_wks("png","steer_"+init)    ; send graphics to PNG file
  gsn_define_colormap(wks,"WhBlGrYeRe")

  res                = True                    
  res@gsnDraw        = False          ; do not draw
  res@gsnFrame       = False          ; do not advance the frame
  res@tiMainString   = "Steering flow (3-8 degree annulus mean) "+init ; add title
  res@vcGlyphStyle   = "WindBarb"     ; select wind barbs 
  res@vcRefLengthF   = 0.045          ; define length of vec ref
  res@vcRefAnnoOn    = False          ; turn off ref wind barb
  res@trYReverse     = True           ; reverse the Y-axis
  res@vcMonoWindBarbColor     = False            ; color barbs by scalar
  res@lbOrientation           = "vertical"       ; vertical labelbar
  res@vcLevelSelectionMode    = "ExplicitLevels" ; will explicitly set the levels
  res@vcWindBarbLineThicknessF= 1.5              ; set the wind barb thickness
  res@vcLevels                = ispan(5,65,5)    ; set the levels 
  res@vcLevelColors           = (/20,30,40,45,50,55,60,65,70,75,80,85,90,95/)   ; set the colors associated with the wind speed
  
  res@vcMapDirection = False			; needed to disassociate the map from the vectors
  res@gsnAddCyclic = False              ; data is not cyclic
  res@gsnLeftString = ""
  res@tiXAxisString = "forecast time (hour)"
  res@tiYAxisString = "pressure (hPa)"
  
  plot = gsn_csm_vector_scalar(wks,uav(lev|:,time|:),vav(lev|:,time|:),wspd(lev|:,time|:),res)

  getvalues plot                        ; get the X/Y axis min/max for use in the loglin plot
    "trXMinF"    : trxmin
    "trXMaxF"    : trxmax
    "trYMinF"    : trymin
    "trYMaxF"    : trymax
  end getvalues

  loglin = create "logling" logLinPlotClass wks     ; draw a loglin plot, with expanded X/Y axis
    "trXMinF"        : trxmin-6
    "trXMaxF"        : trxmax+6
    "trYMinF"        : trymin-40
    "trYMaxF"        : trymax+50
    "trYReverse"     : True
    "vpXF"           : .15            ; set the X-axis NDC starting point
    "vpYF"           : .8             ; set the Y-axis NDC starting point
    "vpWidthF"       : .7             ; set the width of the plot in NDC units
    "vpHeightF"      : .45            ; set the height of the plot in NDC units
  end create
  overlay(loglin,plot)              ; overlay plot with the loglin plot
  draw(loglin)                      ; draw the plot
  frame(wks)                        ; advance the frame
end
