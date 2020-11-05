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
  yyyymmdd = "20191009"
  yyyy = str_get_cols(yyyymmdd, 0, 3)
  mm   = str_get_cols(yyyymmdd, 4, 5)
  mmdd = str_get_cols(yyyymmdd, 4, 7)
  fs00 = addfile ("../../netcdf/tigge/"+yyyy+"/jma/strm_"+mmdd+"00_mean.nc", "r")
  fv00 = addfile ("../../netcdf/tigge/"+yyyy+"/jma/vpot_"+mmdd+"00_mean.nc", "r")
  fs12 = addfile ("../../netcdf/tigge/"+yyyy+"/jma/strm_"+mmdd+"12_mean.nc", "r")
  fv12 = addfile ("../../netcdf/tigge/"+yyyy+"/jma/vpot_"+mmdd+"12_mean.nc", "r")
  fs    = addfile ("../../netcdf/gsm/gl/"+yyyy+"/"+mm+"/strm_init.nc", "r")
  fv    = addfile ("../../netcdf/gsm/gl/"+yyyy+"/"+mm+"/vpot_init.nc", "r")
  psi00  = fs00->var(:,:,:,:)      ; (time,lat,lon)
  chi00  = fv00->var(:,:,:,:)
  psi12  = fs12->var(:,:,:,:)      ; (time,lat,lon)
  chi12  = fv12->var(:,:,:,:)
  psia   = fs->var(:,:,:,:)      ; (time,lat,lon)
  chia   = fv->var(:,:,:,:)
  t00    = fs00->time
  t12    = fs12->time
  ta     = fs->time
  lev  = fs00->level
  leva = fs->level
  lat  = fs00->lat
  lon  = fs00->lon
  nlev = dimsizes(lev)

  date = cd_string(t12, "%H UTC, %D")
      
  scale = 1e6                                      ; arbitrary: used for nicer plot values
  psi00 = psi00/scale
  chi00 = chi00/scale
  psi12 = psi12/scale
  chi12 = chi12/scale
  psia  = psia/scale
  chia  = chia/scale
;*************************************************
; plot results
;*************************************************    
  do l = 0, nlev(0)-1
    la = ind(leva.eq.lev(l))
    wks1  = gsn_open_wks("png","vpot_0912-0900_"+sprinti("%0.3d",toint(lev(l)))) ; send graphics to PNG file
    plot1 = new(12,graphic)
    wks2  = gsn_open_wks("png","strm_0912-0900_"+sprinti("%0.3d",toint(lev(l)))) ; send graphics to PNG file
    plot2 = new(12,graphic)

    resd                      = True

    resd@gsnDraw              = False
    resd@gsnFrame             = False

    resd@cnFillOn             = True                 ; color on
    resd@cnLinesOn            = False                ; turn off contour lines
    ;res@cnLineColor          = "blue"
    ;res@gsnScalarContour     = True                 ; vectors over contours
    ;res@cnFillPalette        = "BlAqGrYeOrReVi200"
    resd@cnFillPalette        = "GMT_polar"

    resd@mpFillOn             = False              ; turn off map fill
    resd@mpMinLonF            = 110.
    resd@mpMaxLonF            = 180.
    resd@mpMinLatF            = 10.
    resd@mpMaxLatF            = 80.
    resd@gsnRightString       = "10^6 m^2/s^2"
    resd@gsnRightStringFontHeightF = 0.015
    resd@gsnRightStringOrthogonalPosF = 0.05
    resd@gsnLeftString        = ""
    resd@cnLevelSelectionMode = "ManualLevels"       ; set manual contour levels
    
    resc                     = True
    resc@gsnDraw             = False
    resc@gsnFrame            = False
    resc@cnFillOn            = False
    resc@cnLinesOn           = True
    resc@gsnContourZeroLineThicknessF = 0.0
    resc@gsnContourNegLineDashPattern = 1
    resc@cnLevelSelectionMode = "ManualLevels"       ; set manual contour levels

    resc@gsnRightString       = ""
    resc@gsnLeftString        = ""
     
    i = 0
    do nt = 0, 11
      ind00 = ind(t00.eq.t12(nt))
      inda = ind(ta.eq.t12(nt))
      chi = chi12(nt,l,:,:)
      ;chi = chi00(ind00,l,:,:)
      chi = chi - chi00(ind00,l,:,:)
      ;chi = chi - chia(inda,la,:,:)
      psi = psi12(nt,l,:,:)
      ;psi = psi00(ind00,l,:,:)
      psi = psi - psi00(ind00,l,:,:)
      ;psi = psi - psia(inda,la,:,:)
    ;res@tiMainString         = "Velocity Potential via Spherical Harmonics"
      resd@gsnCenterString      = date(nt)
      if (l .lt. 2) then
      resd@cnMinLevelValF       = -2.0                 ; set min contour level
      resd@cnMaxLevelValF       =  2.0                 ; set max contour level
      resd@cnLevelSpacingF      =  .2                 ; set contour spacing
      resc@cnMinLevelValF       = -10.0                 ; set min contour level
      resc@cnMaxLevelValF       =  10.0                 ; set max contour level
      resc@cnLevelSpacingF      =  1.0                 ; set contour spacing
      else 
      resd@cnMinLevelValF       = -5.0                 ; set min contour level
      resd@cnMaxLevelValF       =  5.0                 ; set max contour level
      resd@cnLevelSpacingF      =  .5                 ; set contour spacing
      resc@cnMinLevelValF       = -20.0                 ; set min contour level
      resc@cnMaxLevelValF       =  20.0                 ; set max contour level
      resc@cnLevelSpacingF      =  2.0                 ; set contour spacing
      end if
      plotd = gsn_csm_contour_map(wks1,chi,resd)
      plotc = gsn_csm_contour(wks1,chi12(nt,l,:,:),resc)
      ;plotc = gsn_csm_contour(wks1,chi00(ind00,l,:,:),resc)
      overlay(plotd,plotc)
      plot1(i)=plotd
      delete(plotd)
      delete(plotc)

      ;res@tiMainString         = "Streamfunction via Spherical Harmonics"     ; title
      ;res@gsnCenterString      = "Psi scaled by 1e6"  ; center string
      ;res@gsnLeftString        = "10^6 m^2/s^2"
      
      ;res@cnMinLevelValF       = -40.0               ; set min contour level
      ;res@cnMaxLevelValF       =  40.0               ; set max contour level
      ;res@cnLevelSpacingF      =  8.0               ; set contour spacing
      if (l .lt. 2) then
      resd@cnMinLevelValF       = -5.0                 ; set min contour level
      resd@cnMaxLevelValF       =  5.0                 ; set max contour level
      resd@cnLevelSpacingF      =  .5                 ; set contour spacing
      resc@cnMinLevelValF       = -40.0                 ; set min contour level
      resc@cnMaxLevelValF       =  40.0                 ; set max contour level
      resc@cnLevelSpacingF      =  5.0                 ; set contour spacing
      else 
      resd@cnMinLevelValF       = -10.0                 ; set min contour level
      resd@cnMaxLevelValF       =  10.0                 ; set max contour level
      resd@cnLevelSpacingF      =  1.0                 ; set contour spacing
      resc@cnMinLevelValF       = -100.0                 ; set min contour level
      resc@cnMaxLevelValF       =  100.0                 ; set max contour level
      resc@cnLevelSpacingF      =  10.0                 ; set contour spacing
      end if
      plotd = gsn_csm_contour_map(wks2,psi,resd)
      plotc = gsn_csm_contour(wks2,psi12(nt,l,:,:),resc)
      ;plotc = gsn_csm_contour(wks2,psi00(ind00,l,:,:),resc)
      overlay(plotd,plotc)
      plot2(i)=plotd
      delete(plotd)
      delete(plotc)
      i = i + 1
    end do

    resp = True
    resp@gsnPanelMainString  =  "Velocity potential 0912-0900 "+sprinti("%0.3d",toint(lev(l)))+"hPa"
    gsn_panel(wks1,plot1,(/3,4/),resp)
    resp@gsnPanelMainString  =  "Stream function 0912-0900 "+sprinti("%0.3d",toint(lev(l)))+"hPa"
    gsn_panel(wks2,plot2,(/3,4/),resp)
    delete(plot1)
    delete(plot2)
    delete(wks1)
    delete(wks2)
  end do