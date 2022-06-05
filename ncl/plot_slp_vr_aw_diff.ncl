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
  ; asymmetric wind
  fu = addfile("../../netcdf/tigge/"+yyyy+"/"+orig+"/rotate/inv_r8_glb_"+yyyymmddhh+"_"+mem+"_ug.nc", "r")
  fv = addfile("../../netcdf/tigge/"+yyyy+"/"+orig+"/rotate/inv_r8_glb_"+yyyymmddhh+"_"+mem+"_vg.nc", "r")
  fp = addfile("../../netcdf/tigge/"+yyyy+"/"+orig+"/rotate/inv_r8_glb_"+yyyymmddhh+"_"+mem+"_slp.nc", "r")
  ;fh = addfile("../../netcdf/tigge/"+yyyy+"/"+orig+"/rotate/inv_glb_"+yyyymmddhh+"_"+mem+"_gh.nc", "r")
  fua = addfile("../../netcdf/gsm/gl/"+yyyy+"/"+mm+"/rotate/inv_r8_init_ug.nc", "r")
  fva = addfile("../../netcdf/gsm/gl/"+yyyy+"/"+mm+"/rotate/inv_r8_init_vg.nc", "r")
  fpa = addfile("../../netcdf/gsm/gl/"+yyyy+"/"+mm+"/rotate/inv_r8_init_slp.nc", "r")
  ;;fh = addfile("../../netcdf/gsm/gl/"+yyyy+"/"+mm+"/rotate/inv_init_gh.nc", "r")
  uf   = fu->var(:,{850:300},:,:)
  vf   = fv->var(:,{850:300},:,:)
  slpf = fp->var(:,:,:)
  ua   = fua->var(:,{850:300},:,:)
  va   = fva->var(:,{850:300},:,:)
  slpa = fpa->var(:,:,:)
  ;gh = fh->var(:,{500},:,:)
  tf    = fu->time
  ta    = fua->time
  lev  = fu->level
  lat  = fu->lat
  lon  = fu->lon
  
  datef  = cd_string(tf, "%Y%N%D%H")
  print(datef)
  datea  = cd_string(ta, "%Y%N%D%H")
  print(datea)
  indf = ind(datef.eq."2019101012") ;start
  inda = ind(datea.eq."2019101012") ;end
  hour  = cd_inv_calendar(toint(yyyy),toint(mm),toint(dd),toint(hh),0,0,tf@units,0)
  shour = hour
  shour = hour + 6*indf
  sday 	= cd_string(shour, "%N%D%H")

  print(indf)
  print(inda)
  print(sday)
  ;exit
  date = cd_string(tf, "%H UTC, %D")
      
  slp = slpf(indf,:,:)
  slp = slpf(indf,:,:) - slpa(inda,:,:)
  slp = slp * 0.01
  delete(slp@long_name)
  slp@long_name = "mean sea level pressure"
  printMinMax(slp, 0)
;*************************************************
; calc vorticity
;*************************************************
  u = uf(indf,{500},:,:)
  ;u = uf(indf,{500},:,:) - ua(inda,{500},:,:)
  v = vf(indf,{500},:,:)
  ;v = vf(indf,{500},:,:) - va(inda,{500},:,:)
  delete(u@long_name)
  delete(v@long_name)
  u@long_name = "U component of wind"
  v@long_name = "V component of wind"
  vrf = uv2vrF_Wrap(u, v)
  printMinMax(vrf,0)

;*************************************************
; hoskins filtering
;*************************************************
  ndim = dimsizes(vrf)
  ;ntim = ndim(0)
  ;nlev = ndim(1)
  nlat = ndim(0)
  nlon = ndim(1)

  N = nlat
  if (nlon%2 .eq. 0)then
    N = min((/ nlat, (nlon+2)/2 /))
  else
    N = min((/ nlat, (nlon+1)/2 /))
  end if

  T = 106
  NT = 63
  a = new((/ nlat, N /),float)
  b = new((/ nlat, N /),float)
  shaec(vrf, a, b)
  tri_trunc(a, b, T)
  exp_tapersh(a, b, NT, 2)
  shsec(a, b, vrf)
  delete(u)
  delete(v)
;*************************************************
; calc vorticity
;*************************************************
  u = ua(inda,{500},:,:)
  ;u = uf(indf,{500},:,:) - ua(inda,{500},:,:)
  v = va(inda,{500},:,:)
  ;v = vf(indf,{500},:,:) - va(inda,{500},:,:)
  delete(u@long_name)
  delete(v@long_name)
  u@long_name = "U component of wind"
  v@long_name = "V component of wind"
  vra = uv2vrF_Wrap(u, v)
  printMinMax(vra,0)

;*************************************************
; hoskins filtering
;*************************************************
  ndim = dimsizes(vra)
  ;ntim = ndim(0)
  ;nlev = ndim(1)
  nlat = ndim(0)
  nlon = ndim(1)

  N = nlat
  if (nlon%2 .eq. 0)then
    N = min((/ nlat, (nlon+2)/2 /))
  else
    N = min((/ nlat, (nlon+1)/2 /))
  end if

  T = 106
  NT = 63
  a = new((/ nlat, N /),float)
  b = new((/ nlat, N /),float)
  shaec(vra, a, b)
  tri_trunc(a, b, T)
  exp_tapersh(a, b, NT, 2)
  shsec(a, b, vra)

  vr = vrf
  vr = vrf - vra
  vr = vr * 1e6
  copy_VarCoords(u, vr)
  vr@long_name = "Vorticity"
  vr@units = "10^(-6)/s"
  printMinMax(vr, 0)
  ;exit
;*************************************************
; vertical integration
;*************************************************
  sigf = new(3,float)
  sigf(0) = 7./22.
  sigf(1) = 1./2.
  sigf(2) = 2./11.
  print(sigf)
  uf = uf * conform_dims(dimsizes(uf), sigf, 1)
  vf = vf * conform_dims(dimsizes(vf), sigf, 1)
  upf = dim_sum_n(uf, 1)
  upf!0 = "time"
  upf!1 = "lat"
  upf!2 = "lon"
  upf&time = tf
  upf&lat = lat
  upf&lon = lon
  upf@_FillValue = default_fillvalue("float")
  vpf = upf
  vpf = dim_sum_n(vf, 1)

  siga = new(6, float)
  siga(0) = 3./22.
  siga(1) = 5./22.
  siga(2:4) = 2./11.
  siga(5) = 1./11.
  print(siga)
  ua = ua * conform_dims(dimsizes(ua), siga, 1)
  va = va * conform_dims(dimsizes(va), siga, 1)
  upa = dim_sum_n(ua, 1)
  upa!0 = "time"
  upa!1 = "lat"
  upa!2 = "lon"
  upa&time = ta
  upa&lat = lat
  upa&lon = lon
  upa@_FillValue = default_fillvalue("float")
  vpa = upa
  vpa = dim_sum_n(va, 1)

  up = upf(indf,:,:)
  up = upf(indf,:,:) - upa(inda,:,:)
  vp = vpf(indf,:,:)
  vp = vpf(indf,:,:) - vpa(inda,:,:)
;*************************************************
; plot results
;*************************************************    
  wks  = gsn_open_wks("png","slp+vor500+aw_"+mmddhh+"-gsm_"+sday+"_"+orig+mem+"_r8") ; send graphics to PNG file
  ;wks  = gsn_open_wks("png","slp+vor500+aw_gsm_"+sday+"_to_"+eday+"_r8") ; send graphics to PNG file
  ;plot = new(8,graphic)
  
  res                      = True

  res@gsnDraw              = False
  res@gsnFrame             = False

  res@cnFillOn             = True                 ; color on
  res@cnLinesOn            = False               ; turn off contour lines
  res@cnLineColor          = "black"
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
    
  res@gsnRightString       = slp@units
  res@gsnRightStringFontHeightF = 0.015
  res@gsnRightStringOrthogonalPosF = 0.05
  res@gsnLeftString        = ""

  res@cnLevelSelectionMode = "ManualLevels"       ; set manual contour levels
  ;res@cnMinLevelValF       = 1000.0                 ; set min contour level
  ;res@cnMaxLevelValF       = 1020.0                 ; set max contour level
  ;res@cnLevelSpacingF      =  1.0                 ; set contour spacing
  res@cnMinLevelValF       = -20.0                 ; set min contour level
  res@cnMaxLevelValF       =  20.0                 ; set max contour level
  res@cnLevelSpacingF      =  2.0                 ; set contour spacing
    
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
  res2@cnLineLabelsOn       = False
  res2@cnLineColor          = "black"
  
  res2@gsnRightString       = ""
  res2@gsnRightStringFontHeightF = 0.015
  res2@gsnRightStringOrthogonalPosF = 0.05
  res2@gsnLeftString        = ""
  res2@gsnCenterString      = ""

  res2@cnLevelSelectionMode = "ExplicitLevels"       ; set manual contour levels
  res2@cnLevels = (/1000.,1002.,1004.,1006.,1008.,1010.,1011.,1012.,1013.,1014.,1015.,1016.,1017.,1018.,1019.,1020./)
  ;res2@cnLevelSelectionMode = "ManualLevels"       ; set manual contour levels
  ;res2@cnMinLevelValF       = -10.0                ; set min contour level
  ;res2@cnMaxLevelValF       =  10.0                 ; set max contour level
  ;res2@cnLevelSpacingF      =  1.0                 ; set contour spacing
  
  res@gsnCenterString      = date(indf)
    
    ;plotd=gsn_csm_vector_scalar_map(wks, ud(nt,:,:), vd(nt,:,:), div(nt,:,:), res)
  plotd=gsn_csm_vector_scalar_map(wks, up, vp, vr, res)

  ;plots=gsn_csm_contour(wks, slp, res2)
  slp = slpf(indf,:,:) * 0.01
  plots=gsn_csm_contour(wks, slp, res2)

  overlay(plotd, plots)

  draw(plotd)
  frame(wks)