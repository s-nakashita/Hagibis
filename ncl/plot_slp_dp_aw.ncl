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
  orig = "jma"
  mem = "mean"
  yyyy = str_get_cols(yyyymmddhh, 0, 3)
  mm   = str_get_cols(yyyymmddhh, 4, 5)
  dd   = str_get_cols(yyyymmddhh, 6, 7)
  hh   = str_get_cols(yyyymmddhh, 8, 9)
  mmddhh = str_get_cols(yyyymmddhh, 4, 9)
  f   = addfile ("../../netcdf/tigge/"+yyyy+"/"+orig+"/glb_"+yyyymmddhh+"_"+mem+".nc", "r")
  ;u  = f->u(:,{200},:,:)      ; (time,lat,lon) 
  ;v  = f->v(:,{200},:,:)
  ;slp = f->msl(:,:,:)
  ; asymmetric wind
  fu = addfile("../../netcdf/tigge/"+yyyy+"/"+orig+"/rotate/inv_glb_"+yyyymmddhh+"_"+mem+"_ug.nc", "r")
  fv = addfile("../../netcdf/tigge/"+yyyy+"/"+orig+"/rotate/inv_glb_"+yyyymmddhh+"_"+mem+"_vg.nc", "r")
  ; asymmetric slp
  fp = addfile("../../netcdf/tigge/"+yyyy+"/"+orig+"/rotate/inv_glb_"+yyyymmddhh+"_"+mem+"_slp.nc", "r")
  ua = fu->var(:,{850:300},:,:)
  va = fv->var(:,{850:300},:,:)
  slp = fp->var(:,:,:)
  t    = fu->time
  lev  = fu->level
  lat  = fu->lat
  lon  = fu->lon
  ntim = dimsizes(t)
  nlev = dimsizes(lev)
  
  date  = cd_string(t, "%Y%N%D%H")
  print(date)
  ds = ind(date.eq."2019101012") ;start
  de = ind(date.eq."2019101012") ;end
  nd = de - ds + 1
  ;if ( nd%2 .ne. 0 ) then
  ;  nd = (nd + 1) / 2
  ;else
  ;  nd = nd / 2
  ;end if
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
; trackfile
  trackfile = "../pytrack/jma/gtrack"+yyyymmddhh+"_"+mem+".txt"
  track = asciiread(trackfile, (/21, 7/), "float")    
;*************************************************
; calc pressure gradient
;*************************************************
  dpx = slp
  dpy = slp
  gradsf(slp, dpx, dpy)
  printMinMax(dpy,0)

  ;div = div * 1e6
  ;copy_VarCoords(slp, div)
  ;div@long_name = "Divergence"
  ;div@units = "10^(-6)/s"
  slp = slp * 0.01
;*************************************************
; hoskins filtering
;*************************************************
  ndim = dimsizes(dpy)
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
  shaec(dpy, a, b)
  tri_trunc(a, b, T)
  exp_tapersh(a, b, NT, 2)
  shsec(a, b, dpy)
  dpy = dpy * 1e3
  copy_VarCoords(slp, dpx)
  copy_VarCoords(slp, dpy)
  dpx@long_name = "longitudinal gradient of slp"
  dpx@units = "Pa/m"
  dpy@long_name = "latitudinal gradient of slp"
  dpy@units = "Pa/m"
;*************************************************
; vertical integration
;*************************************************
  sig = new(3,float)
  sig(0) = 7./22.
  sig(1) = 1./2.
  sig(2) = 2./11.
  print(sig)
  ua = ua * conform_dims(dimsizes(ua), sig, 1)
  va = va * conform_dims(dimsizes(va), sig, 1)
  up = dim_sum_n(ua, 1)
  up!0 = "time"
  up!1 = "lat"
  up!2 = "lon"
  up&time = t
  up&lat = lat
  up&lon = lon
  up@_FillValue = default_fillvalue("float")
  vp = up
  vp = dim_sum_n(va, 1)
;*************************************************
; plot results
;*************************************************    
  wks  = gsn_open_wks("pdf","slp+dpdy+aw_"+mmddhh+"_"+sday+"_to_"+eday+"_"+orig+mem) ; send graphics to PNG file
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
  res@mpMaxLatF            = 40.
    
  res@gsnRightString       = "10^-3 "+dpy@units
  res@gsnRightStringFontHeightF = 0.015
  res@gsnRightStringOrthogonalPosF = 0.05
  res@gsnLeftString        = ""

  res@cnLevelSelectionMode = "ManualLevels"       ; set manual contour levels
  ;res@cnMinLevelValF       = 1000.0                 ; set min contour level
  ;res@cnMaxLevelValF       = 1020.0                 ; set max contour level
  ;res@cnLevelSpacingF      =  1.0                 ; set contour spacing
  res@cnMinLevelValF       = -2.5                 ; set min contour level
  res@cnMaxLevelValF       = 2.5                 ; set max contour level
  res@cnLevelSpacingF      =  .25                 ; set contour spacing
    
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

  res2@cnLevelSelectionMode = "ExplicitLevels"       ; set manual contour levels
  res2@cnLevels = (/1000.,1002.,1004.,1006.,1008.,1010.,1011.,1012.,1013.,1014.,1015.,1016.,1017.,1018.,1019.,1020./)
  ;res2@cnLevelSelectionMode = "ManualLevels"       ; set manual contour levels
  ;res2@cnMinLevelValF       = 1008.0                 ; set min contour level
  ;res2@cnMaxLevelValF       = 1020.0                 ; set max contour level
  ;res2@cnLevelSpacingF      =  1.0                 ; set contour spacing

  ; TC center
  mres = True
  mres@gsMarkerIndex = 12
  mres@gsMarkerSizeF = 0.025
  mres@gsMarkerColor = "black"
  mres@gsMarkerThicknessF = 5.0
  
  ;i = 8 - nd
  i = 0
  do nt = ds, de, 1
    ;res@tiMainString         = "Velocity Potential via Spherical Harmonics"
    res@gsnCenterString      = date(nt)
    
    ;plotd=gsn_csm_vector_scalar_map(wks, ud(nt,:,:), vd(nt,:,:), div(nt,:,:), res)
    plotd=gsn_csm_vector_scalar_map(wks, up(nt,:,:), vp(nt,:,:), dpy(nt,:,:), res)

    plots=gsn_csm_contour(wks, slp(nt,:,:), res2)

    overlay(plotd, plots)

    lonc = track(nt, 4)
    latc = track(nt, 5)
    dum = gsn_add_polymarker(wks, plotd, lonc, latc, mres)

    plot(i) = plotd
  
    i = i + 1
  end do

  ;resp = True
  ;;resp@gsnPanelMainString  =  "slp + divergence at 200hPa + asymmetric wind itime="+mmddhh+" "+orig
  ;gsn_panel(wks,plot,(/2,4/),resp)
  draw(plot(0))
  frame(wks)
  