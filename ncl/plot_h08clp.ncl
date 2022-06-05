load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
begin 
   ;Reading Historical data
   yyyymmddhhnn = "201910090000"
   yyyy = str_get_cols(yyyymmddhhnn, 0, 3)
   mm   = str_get_cols(yyyymmddhhnn, 4, 5)
   dd   = str_get_cols(yyyymmddhhnn, 6, 7)
   hh   = str_get_cols(yyyymmddhhnn, 8, 9)
   nn   = str_get_cols(yyyymmddhhnn, 10, 11)
   fh1     = addfile("/Volumes/dandelion/netcdf/jaxa_himawari/L2/CLP/010/"+yyyy+mm+"/"+dd+"/"+hh+"/NC_H08_"+yyyy+mm+dd+"_"+hh+nn+"_L2CLP010_FLDK.02401_02401.nc","r")
   
   var = fh1->CLTH(:, :) ; cloud top height (km)
   lon = fh1->longitude
   lat = fh1->latitude
   minlon = min(lon)
   maxlon = max(lon)
   minlat = min(lat)
   maxlat = max(lat)

   pltDir = "./"                ; where output will be written
   pltName= "H08CLP_clth_"+yyyymmddhhnn
   pltType= "png"               ; ps, eps, png, pdf, x11 .....

;****************************************************************************
; Graphic resources
;****************************************************************************
  pltPath = pltDir+pltName
  wks     = gsn_open_wks(pltType, pltPath) 

  res                            = True     ; plot mods desired

  res@cnFillOn                   = True     ; turn on color fill
  res@cnLinesOn                  = False    ; turn off contour lines
  res@cnLineLabelsOn             = False    ; turn off contour line labels
  res@cnFillPalette              = "WhBlGrYeRe"

  ; Set resources for map
  res@mpMinLatF                  = -10.
  res@mpMaxLatF                  = 55.
  res@mpMinLonF                  = 120.
  res@mpMaxLonF                  = 180.

  res@mpDataBaseVersion          = "MediumRes"
  res@mpOutlineOn                = True
  res@mpFillOn                   = False
  ;res@mpOutlineBoundarySets      = "National"

  res@gsnAddCyclic               = False            ; regional data, don't add cyclic point
 
  ; Set resources for tickmarks
  res@gsnStringFont              = "helvetica-bold"  ; will be used for gsnCenterString later
  res@gsnCenterStringFontHeightF = 0.02              ; ditto
  ;res@gsnLeftString              = ""                ; remove the left and right substrings
  ;res@gsnRightString             = ""
  res@gsnCenterString            = yyyy+mm+dd+" "+hh+":"+nn

  ; Set resources specifically for plots in column 1
  res@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
  res@cnMinLevelValF            = 0.           ; min level
  res@cnMaxLevelValF            = 18.           ; max level
  res@cnLevelSpacingF           = 1.            ; contour interval

;****************************************************************************
; Thinning data
;****************************************************************************
  lon_out = lon(5::10)
  lat_out = lat(5::10)
  varin   = new(dimsizes(var), float)
  varin!0 = "longitude"
  varin&longitude = lon
  varin!1 = "latitude"
  varin&latitude = lat
  varin = (/short2flt(var)/)
  printMinMax(varin, 0)
  deg2rad = acos(-1.) / 180.
  varout = area_hi2lores_Wrap(lon, lat, varin, False, cos(lat * deg2rad), lon_out, lat_out, False)
  printMinMax(varout, 0)
  plot = gsn_csm_contour_map_ce(wks, varout, res)
end
