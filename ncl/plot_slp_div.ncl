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
  yyyymmddhh = "2019100912"
  orig = "ecmwf"
  yyyy = str_get_cols(yyyymmddhh, 0, 3)
  mm   = str_get_cols(yyyymmddhh, 4, 5)
  dd   = str_get_cols(yyyymmddhh, 6, 7)
  hh   = str_get_cols(yyyymmddhh, 8, 9)
  mmddhh = str_get_cols(yyyymmddhh, 4, 9)
  f   = addfile ("../../netcdf/tigge/"+yyyy+"/"+orig+"/glb_"+yyyymmddhh+"_mean.nc", "r")
  u  = f->u(:,{200},:,:)      ; (time,lat,lon) 
  v  = f->v(:,{200},:,:)
  slp = f->msl(:,:,:)
  t    = f->time
  lev  = f->level
  lat  = f->latitude
  lon  = f->longitude
  ntim = dimsizes(t)
  nlev = dimsizes(lev)

  date  = cd_string(t, "%Y%N%D%H")
  print(date)
  ds = ind(date.eq.yyyymmddhh) ;start
  de = ind(date.eq."2019101212") ;end
  nd = de - ds + 1
  if ( nd%2 .ne. 0 ) then
    nd = (nd + 1) / 2
  else
    nd = nd / 2
  end if
  print(nd)
  print(t@units)
  hour  = cd_inv_calendar(toint(yyyy),toint(mm),toint(dd),toint(hh),0,0,t@units,0)
  shour = hour
  shour = hour + 6*ds
  ehour = hour
  ehour	= hour + 6*de
  sday 	= cd_string(shour, "%N%D%H")
  eday 	= cd_string(ehour, "%N%D%H")
  print("start = "+sday+", end = "+eday)
  delete(date)

  date = cd_string(t, "%H UTC, %D")
      
  slp = slp * 0.01
;*************************************************
; calc divergent wind
;*************************************************
  div = uv2dvF_Wrap(u, v)
  printMinMax(div,0)

  ud = new(dimsizes(u), typeof(u), u@_FillValue)
  vd = new(dimsizes(v), typeof(v), v@_FillValue)

  dv2uvf(div, ud, vd)
  copy_VarCoords(u, ud)
  copy_VarCoords(v, vd)
  ud@long_name = "Zonal Divergent Wind"
  ud@units = u@units
  vd@long_name = "Meridional Divergent Wind"
  vd@units = v@units
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
  a = new((/ ntim, nlat, N /),float)
  b = new((/ ntim, nlat, N /),float)
  shaec(div, a, b)
  tri_trunc(a, b, T)
  exp_tapersh(a, b, 45, 2)
  shsec(a, b, div)
  div = div * 1e6
  copy_VarCoords(slp, div)
  div@long_name = "Divergence"
  div@units = "10^(-6)/s"
;*************************************************
; plot results
;*************************************************    
  wks  = gsn_open_wks("png","slp+div200_"+mmddhh+"_"+sday+"_to_"+eday+"_"+orig) ; send graphics to PNG file
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
  
  res2@gsnRightString       = ""
  res2@gsnRightStringFontHeightF = 0.015
  res2@gsnRightStringOrthogonalPosF = 0.05
  res2@gsnLeftString        = ""
  res2@gsnCenterString      = ""

  res2@cnLevelSelectionMode = "ManualLevels"       ; set manual contour levels
  res2@cnMinLevelValF       = 1000.0                 ; set min contour level
  res2@cnMaxLevelValF       = 1020.0                 ; set max contour level
  res2@cnLevelSpacingF      =  1.0                 ; set contour spacing
  
  i = 1
  do nt = ds, de, 2
    ;res@tiMainString         = "Velocity Potential via Spherical Harmonics"
    res@gsnCenterString      = date(nt)
    
    plotd=gsn_csm_vector_scalar_map(wks, ud(nt,:,:), vd(nt,:,:), div(nt,:,:), res)

    plots=gsn_csm_contour(wks, slp(nt,:,:), res2)

    overlay(plotd, plots)

    plot(i) = plotd
  
    i = i + 1
  end do

  resp = True
  resp@gsnPanelMainString  =  "slp + divergence at 200hPa itime="+mmddhh+" "+orig
  gsn_panel(wks,plot,(/2,4/),resp)
  