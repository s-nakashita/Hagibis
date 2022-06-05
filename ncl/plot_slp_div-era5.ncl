;*************************************************
; plot_slp_div.ncl
;
; Concepts illustrated:
;   - Use spherical harmonics on a global grid
;   - Use uv2dvF_Wrap and uv2vrF_Wrap to compute divergence and vortivity
;     and maintain meta data
;   - Use ilapsG_Wrap to compute velocity potential and streamfunction
;     and maintain meta data
;   - Plot vector winds and scalar quantities
;*************************************************
; wind_4.ncl
  load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
;*************************************************
; open file and read in data: data are on a gaussian grid
;*************************************************
;  yyyymmddhh = "2019100912"
;  idate = "2019101000"
;  edate = "2019101118"
;  dt = 1
;  orig = "jma"
;  yyyymmddhh = idate
;  orig = "ecmwf"
  yyyy = str_get_cols(yyyymmddhh, 0, 3)
  mm   = str_get_cols(yyyymmddhh, 4, 5)
  dd   = str_get_cols(yyyymmddhh, 6, 7)
  hh   = str_get_cols(yyyymmddhh, 8, 9)
  mmddhh = str_get_cols(yyyymmddhh, 4, 9)
  f    = addfile ("../../netcdf/tigge/"+yyyy+"/"+orig+"/glb_"+yyyymmddhh+"_mean.nc", "r")
  tf   = f->time
  lev  = f->level
  lat  = f->latitude
  lon  = f->longitude
  ntim = dimsizes(tf)
  nlev = dimsizes(lev)
  date  = cd_string(tf, "%Y%N%D%H")
  print(date)
  ds = ind(date.eq.idate) ;start
  de = ind(date.eq.edate) ;end
  nd = de - ds + 1
  nplot = nd
  if ( nplot%dt .ne. 0 ) then
    nplot = (nplot + 1) / dt
  else
    nplot = nplot / dt
  end if
  print(nplot)
  print(tf@units)
  hour  = cd_inv_calendar(toint(yyyy),toint(mm),toint(dd),toint(hh),0,0,tf@units,0)
  shour = hour
  shour = hour + 6*ds
  ehour = hour
  ehour	= hour + 6*de
  sday 	= cd_string(shour, "%N%D%H")
  eday 	= cd_string(ehour, "%N%D%H")
  print("start = "+sday+", end = "+eday)
  delete(date)
  uf   = f->u(ds:de,{200},:,:)      ; (time,lat,lon) 
  vf   = f->v(ds:de,{200},:,:)
  psf  = f->msl(ds:de,:,:)
  print(dimsizes(uf))
  
  ; era5
  f1 = addfile("../../netcdf/era5/"+yyyy+"/"+mm+"/msl+lwr.nc","r")
  f2 = addfile("../../netcdf/era5/"+yyyy+"/"+mm+"/uv200.nc","r")
  ta  = f2->time
  datea = cd_string(ta, "%Y%N%D%H")
  i = ind(datea.eq.idate)
  j = ind(datea.eq.edate)
  psa = f1->msl(i:j,:,:)
  ua  = f2->u(i:j,:,:)
  va  = f2->v(i:j,:,:)
;*************************************************
; interpolation
;*************************************************
  opt = False
  newdim = dimsizes(psf)
  newnlat = newdim(1)
  newnlon = newdim(2)
  psnew = f2fsh_Wrap(psa, (/newnlat, newnlon/))
  psnew@long_name = psa@long_name
  delete(psa)
  psa = psnew
  unew = f2fsh_Wrap(ua, (/newnlat, newnlon/))
  unew@long_name = ua@long_name
  delete(ua)
  ua = unew
  vnew = f2fsh_Wrap(va, (/newnlat, newnlon/))
  vnew@long_name = va@long_name
  delete(va)
  va = vnew  


  date = cd_string(tf, "%H UTC, %D")
  ps = psf
  ps = (psf - psa) * 0.01
;*************************************************
; calc divergent wind
;*************************************************
  divf = uv2dvF_Wrap(uf, vf)
  printMinMax(divf,0)

  diva = uv2dvF_Wrap(ua, va)
  printMinMax(diva,0)

  div = divf
  div = div - diva

  ud = new(dimsizes(uf), typeof(uf), uf@_FillValue)
  vd = new(dimsizes(vf), typeof(vf), vf@_FillValue)

  dv2uvf(div, ud, vd)
  copy_VarCoords(uf, ud)
  copy_VarCoords(vf, vd)
  ud@long_name = "Zonal Divergent Wind"
  ud@units = uf@units
  vd@long_name = "Meridional Divergent Wind"
  vd@units = vf@units
  ;div = div * 1e6
  ;copy_VarCoords(slp, div)
  ;div@long_name = "Divergence"
  ;div@units = "10^(-6)/s"

;*************************************************
; hoskins filtering
;*************************************************
  ndim = dimsizes(div)
  ntim = ndim(0)
  nlat = ndim(1)
  nlon = ndim(2)

  N = nlat
  if (nlon%2 .eq. 0)then
    N = min((/ nlat, (nlon+2)/2 /))
  else
    N = min((/ nlat, (nlon+1)/2 /))
  end if

  T = 106
  NT = 63
  a = new((/ ntim, nlat, N /),float)
  b = new((/ ntim, nlat, N /),float)
  shaec(div, a, b)
  tri_trunc(a, b, T)
  exp_tapersh(a, b, NT, 2)
  shsec(a, b, div)
  div = div * 1e6
  copy_VarCoords(psf, div)
  div@long_name = "Divergence"
  div@units = "10^(-6)/s"
;*************************************************
; plot results
;*************************************************    
  wks  = gsn_open_wks("png","slp+div200-era5_"+mmddhh+"_"+sday+"_to_"+eday+"_"+orig) ; send graphics to PNG file
  plot = new(8,graphic)
  
  res                      = True

  res@gsnDraw              = False
  res@gsnFrame             = False

  res@cnFillOn             = True                 ; color on
  res@cnLinesOn            = False               ; turn off contour lines
  res@cnLineColor          = "blue"
  res@gsnContourZeroLineThicknessF = 0.0
  res@gsnContourNegLineDashPattern = 1
  res@gsnScalarContour     = True                 ; vectors over contours
  ;res@cnFillPalette        = "BlAqGrYeOrReVi200"
  res@cnFillPalette        = "GMT_polar"

  res@mpFillOn             = False              ; turn off map fill
  res@mpMinLonF            = 110.
  res@mpMaxLonF            = 180.
  res@mpMinLatF            = -10.
  res@mpMaxLatF            = 30.
    
  res@gsnRightString       = div@units
  res@gsnRightStringFontHeightF = 0.015
  res@gsnRightStringOrthogonalPosF = 0.05
  res@gsnLeftString        = ""

  res@cnLevelSelectionMode = "ManualLevels"       ; set manual contour levels
  ;res@cnMinLevelValF       = 1000.0                 ; set min contour level
  ;res@cnMaxLevelValF       = 1020.0                 ; set max contour level
  ;res@cnLevelSpacingF      =  1.0                 ; set contour spacing
  res@cnMinLevelValF       = -50.0                 ; set min contour level
  res@cnMaxLevelValF       = 50.0                 ; set max contour level
  res@cnLevelSpacingF      =  10.0                 ; set contour spacing
    
  res@vcRefMagnitudeF           = 10.0             ; define vector ref mag
  res@vcRefLengthF              = 0.045           ; define length of vec ref
  res@vcRefAnnoOrthogonalPosF   = -1.0            ; move ref vector
  res@vcRefAnnoArrowLineColor   = "black"         ; change ref vector color
  res@vcRefAnnoArrowUseVecColor = False           ; don't use vec color for ref

  res@vcGlyphStyle = "CurlyVector"
  res@vcLineArrowColor        = "black"           ; change vector color
  res@vcLineArrowThicknessF   = 2.0               ; change vector thickness
  res@vcMinDistanceF	      = 0.03
  res@vcVectorDrawOrder       = "PostDraw"        ; draw vectors last
  
  res2                      = True

  res2@gsnDraw              = False
  res2@gsnFrame             = False

  res2@cnFillOn             = False                 ; color on
  res2@cnLinesOn            = True               ; turn off contour lines
  res2@cnLineColor          = "black"
  res2@gsnContourZeroLineThicknessF = 0.0
  res2@gsnContourNegLineDashPattern = 1
  res2@gsnRightString       = ""
  res2@gsnRightStringFontHeightF = 0.015
  res2@gsnRightStringOrthogonalPosF = 0.05
  res2@gsnLeftString        = ""
  res2@gsnCenterString      = ""

  res2@cnLevelSelectionMode = "ManualLevels"       ; set manual contour levels
  res2@cnMinLevelValF       = -10.0                 ; set min contour level
  res2@cnMaxLevelValF       = 10.0                 ; set max contour level
  res2@cnLevelSpacingF      =  1.0                 ; set contour spacing
  
  i = 8 - nplot
  print(dt)
  do nt = 0, nd-1, dt
    ;res@tiMainString         = "Velocity Potential via Spherical Harmonics"
    res@gsnCenterString      = date(ds+nt)
    
    plotd=gsn_csm_vector_scalar_map(wks, ud(nt,:,:), vd(nt,:,:), div(nt,:,:), res)

    plots=gsn_csm_contour(wks, ps(nt,:,:), res2)

    overlay(plotd, plots)

    plot(i) = plotd
  
    i = i + 1
    
  end do

  resp = True
  resp@gsnPanelMainString  =  "slp + divergence at 200hPa itime="+mmddhh+" "+orig
  gsn_panel(wks,plot,(/2,4/),resp)
  