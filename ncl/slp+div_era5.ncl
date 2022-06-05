;----------------------------------------------------------------------
; vector_1.ncl
;
; Concepts illustrated:
;   - Overlaying vectors and filled contours on a map
;   - Manually attaching lat/lon coordinate arrays to a variable
;   - Changing the length of the smallest vector as a fraction of the reference vector
;   - Moving the vector reference annotation to the top right of the plot
;   - Drawing curly vectors
;   - Setting the color for vectors
;   - Making the labelbar be vertical
;   - Increasing the thickness of vectors
;----------------------------------------------------------------------
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

begin

;---Open netCDF file
  ;hh = 12
  id = 0
  idate = 2019100900
  edate = 2019101018
  sday  = str_get_cols(sprinti("%0.10d",idate), 4, 7)
  eday  = str_get_cols(sprinti("%0.10d",edate), 4, 7)
  cmap  = read_colormap_file("GMT_polar")
  wks   = gsn_open_wks("png","slp+div200_era5_"+sday+"-"+eday)
  plot  = new(4*2,graphic)
  yyyymmdd = str_get_cols(sprinti("%0.10d",idate), 0, 7)
  yyyy     = str_get_cols(sprinti("%0.10d",idate), 0, 3)
  mm       = str_get_cols(sprinti("%0.10d",idate), 4, 5)
  f1 = addfile("../../netcdf/era5/"+yyyy+"/"+mm+"/msl+lwr.nc","r")
  f2 = addfile("../../netcdf/era5/"+yyyy+"/"+mm+"/uv200.nc","r")
  ;b = addfile("../../netcdf/nasa/merra2/MERRA2_400.tavg1_2d_rad_Nx."+yyyymmdd+".nc4","r")
  lat     = f1->latitude
  lon     = f1->longitude
  time    = f1->time
  date = cd_string(time, "%Y%N%D%H")
  ds = ind(date.eq.sprinti("%0.10d",idate))
  de = ind(date.eq.sprinti("%0.10d",edate))
  print(ds)
  print(de)
  delete(date)
  date    = cd_string(time, "%N%D %H UTC")
  print(date)
  
  ps = f1->msl(:,:,:)
  u = f2->u(:,:,:)
  v = f2->v(:,:,:) 
  ps = ps * 0.01
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
  NT = 63
  a = new((/ ntim, nlat, N /),float)
  b = new((/ ntim, nlat, N /),float)
  shaec(div, a, b)
  tri_trunc(a, b, T)
  exp_tapersh(a, b, NT, 2)
  shsec(a, b, div)
  div = div * 1e6
  copy_VarCoords(ps, div)
  div@long_name = "Divergence"
  div@units = "10^(-6)/s"
;*************************************************
; plot results
;************************************************* 
     res                  = True               ; plot mods desired
     res@gsnDraw	      = False
     res@gsnFrame	      = False
;     res@cnInfoLabelOn	      = False
     
     res@gsnRightString	   = div@units
     res@gsnLeftString	   = ""
     res@gsnCenterStringOrthogonalPosF = 0.02
     res@gsnAddCyclic	      = False
     res@cnFillOn             = True               ; turn on color for contours
     res@cnLinesOn            = False               ; turn off contour lines
     res@cnLineLabelsOn       = False              ; turn off contour line labels
     res@gsnScalarContour     = True               ; (for gsn_csm_vector_scalar_map)draw contours
     res@cnFillPalette        = cmap

     res@mpFillOn             = False 
;     res@mpLandFillColor      = "gray"            ; set land to be gray
     res@mpMinLonF            =  110.               ; select a subregion
     res@mpMaxLonF            =  180.
     res@mpMinLatF            =  -10. 
     res@mpMaxLatF            =  30.

     res@lbLabelBarOn         = False
     ;res@lbOrientation            = "Vertical"     ; vertical label bar
     ;res@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer

; note, when doing a subregion, NCL determines the range of the data from
; the full domain. If you wish to just consider the domain you are plotting,
; you must manually set those levels.

     res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
     res@cnMaxLevelCount	   = 20
     res@cnMinLevelValF       = -50.0               ; set min contour level
     res@cnMaxLevelValF       = 50.0                 ; set max contour level
     res@cnLevelSpacingF      = 10.0               ; set contour spacing

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
     
     resp = True
     resp@gsnDraw	      = False
     resp@gsnFrame	      = False
     resp@cnLevelSelectionMode = "ManualLevels"
     resp@cnMinLevelValF       = 1000.0               ; set min contour level
     resp@cnMaxLevelValF       = 1080.0                 ; set max contour level
     resp@cnLevelSpacingF      = 1.0               ; set contour spacing
     resp@cnFillOn             = False               ; turn on color for contours
     resp@cnLinesOn            = True                ; turn off contour lines
     resp@cnLineLabelsOn        = True
     resp@cnLineColor          = "black"
     resp@gsnRightString   = ""    
     resp@gsnLeftString    = ""
     resp@gsnCenterString  = ""
  do d = ds, de
     res@gsnCenterString	   = date(d)
     plotdz = gsn_csm_vector_scalar_map_ce(wks, ud(d,:,:), vd(d,:,:), div(d,:,:), res)

     plotz = gsn_csm_contour(wks,ps(d,:,:),resp)
     overlay(plotdz,plotz)
     plot(id) = plotdz
     id = id + 1
     delete(plotdz)
     delete(plotz)

   end do
   resP	                   = True
   resP@gsnPanelMainString	   = "SLP + Divergence ad 200hPa"
   resP@gsnPanelLabelBar     = True
   gsn_panel(wks,plot,(/2,4/),resP)
   delete(plot)
end
