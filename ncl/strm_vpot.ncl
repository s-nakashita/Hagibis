;*************************************************
; wind_3.ncl
;
; Concepts illustrated:
;   - Use spherical harmonics on a global grid
;   - Use uv2dvG_Wrap and uv2vr_Wrap to compute divergence and vortivity
;     and maintain meta data
;   - Use ilapsG_Wrap to compute velocity potential and streamfunction
;     and maintain meta data
;   - Plot vector winds and scalr quantities
;*************************************************
; wind_4.ncl
load "write_nc.ncl"
;*************************************************
; open file and read in data: data are on a gaussian grid
;*************************************************
  yyyymmddhh = "2019100900"
  yyyy = str_get_cols(yyyymmddhh, 0, 3)
  mm   = str_get_cols(yyyymmddhh, 4, 5)
  mmddhh = str_get_cols(yyyymmddhh, 4, 9)
  f    = addfile ("../../netcdf/tigge/"+yyyy+"/jma/glb_"+mmddhh+"_mean.nc", "r")
  ;f    = addfile ("../../netcdf/gsm/gl/"+yyyy+"/"+mm+"/init.nc", "r")
  ug   = f->UGRD(:,:,:,:)                      ; (time,lat,lon)
  vg   = f->VGRD(:,:,:,:)
  t    = f->time
  lev  = f->level
  lat  = f->lat
  lon  = f->lon

  strm = new(dimsizes(ug),typeof(ug)) 
  strm@long_name = "streamfunction"
  strm@units     = "m^2/s^2"
  vpot = new(dimsizes(ug),typeof(ug))
  vpot@long_name = "velocity potential"
  vpot@units     = "m^2/s^2"
  nlev = dimsizes(lev)
  do l = 0, nlev(0)-1
  u = ug(:,l,:,:)
  v = vg(:,l,:,:)
;*************************************************
; Calculate divergence and vortivity via spherical harmonics
; Use Wrap version to maintain meta data
;*************************************************
  div = uv2dvG_Wrap(u,v)           ; u,v ==> divergence
  vrt = uv2vrG_Wrap(u,v)           ; u,v ==> vorticity

;*************************************************
; Calculate velocity potential and stream function  
; Use Wrap version to maintain meta data
;*************************************************
  chi = ilapsG_Wrap ( div , 0)
  chi@long_name = "velocity potential"
  chi@units     = "m/s" 

  vpot(:,l,:,:) = chi

  psi = ilapsG_Wrap ( vrt , 0)
  psi@long_name = "streamfunction"
  psi@units     = "m/s" 

  strm(:,l,:,:) = psi
  end do
  ofile = "../../netcdf/tigge/"+yyyy+"/jma/strm_"+mmddhh+"_mean.nc"
  ;ofile = "../../netcdf/gsm/gl/"+yyyy+"/"+mm+"/strm_init.nc"
  if (isfilepresent(ofile)) then
  system("rm " + ofile)
  end if
  write_nc(ofile,strm,t,lev,lat,lon)
  ofile = "../../netcdf/tigge/"+yyyy+"/jma/vpot_"+mmddhh+"_mean.nc"
  ;ofile = "../../netcdf/gsm/gl/"+yyyy+"/"+mm+"/vpot_init.nc"
  if (isfilepresent(ofile)) then
  system("rm " + ofile)
  end if
  write_nc(ofile,vpot,t,lev,lat,lon)
