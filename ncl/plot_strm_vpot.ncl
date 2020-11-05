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
  load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
;*************************************************
; open file and read in data: data are on a gaussian grid
;*************************************************
  yyyymmddhh = "2019100900"
  yyyy = str_get_cols(yyyymmddhh, 0, 3)
  mm   = str_get_cols(yyyymmddhh, 4, 5)
  mmddhh = str_get_cols(yyyymmddhh, 4, 9)
  ;fs   = addfile ("../../netcdf/tigge/"+yyyy+"/jma/strm_"+mmddhh+"_mean.nc", "r")
  ;fv   = addfile ("../../netcdf/tigge/"+yyyy+"/jma/vpot_"+mmddhh+"_mean.nc", "r")
  fs    = addfile ("../../netcdf/gsm/gl/"+yyyy+"/"+mm+"/strm_init.nc", "r")
  fv    = addfile ("../../netcdf/gsm/gl/"+yyyy+"/"+mm+"/vpot_init.nc", "r")
  psi  = fs->var(:,:,:,:)      ; (time,lat,lon)
  chi  = fv->var(:,:,:,:)
  t    = fs->time
  lev  = fs->level
  lat  = fs->lat
  lon  = fs->lon
  ntim = dimsizes(t)
  nlev = dimsizes(lev)

  date = cd_string(t, "%H UTC, %D")
      
  scale = 1e6                                      ; arbitrary: used for nicer plot values
  psi = psi/scale
  chi = chi/scale
;*************************************************
; plot results
;*************************************************    
  do l = 0, nlev(0)-1
    
    ;wks1  = gsn_open_wks("png","vpot_"+mmddhh+"_"+sprinti("%0.3d",toint(lev(l)))) ; send graphics to PNG file
    wks1  = gsn_open_wks("png","vpot_anal_"+sprinti("%0.3d",toint(lev(l)))) ; send graphics to PNG file
    plot1 = new(8,graphic)
    ;wks2  = gsn_open_wks("png","strm_"+mmddhh+"_"+sprinti("%0.3d",toint(lev(l)))) ; send graphics to PNG file
    wks2  = gsn_open_wks("png","strm_anal_"+sprinti("%0.3d",toint(lev(l)))) ; send graphics to PNG file
    plot2 = new(8,graphic)

    res                      = True

    res@gsnDraw              = False
    res@gsnFrame             = False

    res@cnFillOn             = False                 ; color on
    res@cnLinesOn            = True                ; turn off contour lines
    res@cnLineColor          = "blue"
    res@gsnContourZeroLineThicknessF = 0.0
    res@gsnContourNegLineDashPattern = 1
    ;res@gsnScalarContour     = True                 ; vectors over contours
    ;res@cnFillPalette        = "BlAqGrYeOrReVi200"
    res@cnFillPalette        = "amwg256"

    res@mpFillOn             = False              ; turn off map fill
    res@mpMinLonF            = 110.
    res@mpMaxLonF            = 180.
    res@mpMinLatF            = -10.
    res@mpMaxLatF            = 60.
    
    res@gsnRightString       = "10^6 m^2/s^2"
    res@gsnRightStringFontHeightF = 0.015
    res@gsnRightStringOrthogonalPosF = 0.05
    res@gsnLeftString        = ""

    res@cnLevelSelectionMode = "ManualLevels"       ; set manual contour levels
     
    i = 0
    do nt = 0, 7
    ;res@tiMainString         = "Velocity Potential via Spherical Harmonics"
      res@gsnCenterString      = date(nt)
      ;if (l .lt. 2) then
      if (l .lt. 6) then
      res@cnMinLevelValF       = -10.0                 ; set min contour level
      res@cnMaxLevelValF       =  10.0                 ; set max contour level
      res@cnLevelSpacingF      =  1.0                 ; set contour spacing
      else 
      res@cnMinLevelValF       = -20.0                 ; set min contour level
      res@cnMaxLevelValF       =  20.0                 ; set max contour level
      res@cnLevelSpacingF      =  2.0                 ; set contour spacing
      end if
      plot1(i)=gsn_csm_contour_map(wks1,chi(nt,l,:,:),res)
  
      ;res@tiMainString         = "Streamfunction via Spherical Harmonics"     ; title
      ;res@gsnCenterString      = "Psi scaled by 1e6"  ; center string
      ;res@gsnLeftString        = "10^6 m^2/s^2"
      
      ;res@cnMinLevelValF       = -40.0               ; set min contour level
      ;res@cnMaxLevelValF       =  40.0               ; set max contour level
      ;res@cnLevelSpacingF      =  8.0               ; set contour spacing
      ;if (l .lt. 2) then
      if (l .lt. 6) then
      res@cnMinLevelValF       = -40.0                 ; set min contour level
      res@cnMaxLevelValF       =  40.0                 ; set max contour level
      res@cnLevelSpacingF      =  5.0                 ; set contour spacing
      else 
      res@cnMinLevelValF       = -100.0                 ; set min contour level
      res@cnMaxLevelValF       =  100.0                 ; set max contour level
      res@cnLevelSpacingF      =  10.0                 ; set contour spacing
      end if
      plot2(i)=gsn_csm_contour_map(wks2,psi(nt,l,:,:),res)
      i = i + 1
    end do

    resp = True
    resp@gsnPanelMainString  =  "Velocity potential "+sprinti("%0.3d",toint(lev(l)))+"hPa"
    gsn_panel(wks1,plot1,(/3,3/),resp)
    resp@gsnPanelMainString  =  "Stream function "+sprinti("%0.3d",toint(lev(l)))+"hPa"
    gsn_panel(wks2,plot2,(/3,3/),resp)
    delete(plot1)
    delete(plot2)
    delete(wks1)
    delete(wks2)
  end do