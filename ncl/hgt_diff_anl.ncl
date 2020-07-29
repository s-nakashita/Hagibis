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
  hh = "00"
  a = addfile("../../netcdf/gsm/gl/2019/10/init.nc","r")
  b = addfile("../../netcdf/tigge/2019/jma/1009"+hh+"_mean.nc","r")
;
; Read in Sea Surface Temperature Anomalies
; Read in U and V at 1000 mb [subscript 0]
; Determine the subscripts corresponding to Jan 88
;
  za      = a->HGT(:,:,180:340,180:360)
  zf      = b->HGT(:,:,:,:)
  lat_uv   = b->lat
  lon_uv   = b->lon
  leva    = a->level
  levf    = b->level
  time_a  = a->time
  time_f  = b->time
  date_a  = cd_string(time_a, "%Y%N%D%H")
  date_f  = cd_string(time_f, "%Y%N%D%H")
  print(date_a)
  print(date_f)

  sig = new(3,float)
  sig(0) = 3./7. + 4./49.
  sig(1) = 10./49. + 1./7.
  sig(2) = 1./7.
  print(sig)

  do l = 0,4 ; level
  ds = 2
  de = 12
  vdate = 2019100900 + toint(hh)
  year  = toint(str_get_cols(sprinti("%0.10d",vdate), 0, 3))
  month = toint(str_get_cols(sprinti("%0.10d",vdate), 4, 5))
  day   = toint(str_get_cols(sprinti("%0.10d",vdate), 6, 7))
  hour  = toint(str_get_cols(sprinti("%0.10d",vdate), 8, 9))
  sec = cd_inv_calendar(year,month,day,hour,0,0, time_f@units, 0)
  newsec = sec
  newsec = sec + 6*60*60*ds
  vdate  = toint(cd_string(newsec, "%Y%N%D%H"))
  wks = gsn_open_wks("png","hgt_09"+hh+"-anl_"+levf(l)+"hPa_"+date_f(ds)+"to"+date_f(de))
  plot = new(6,graphic)
  cmap = read_colormap_file("BlueDarkRed18")
    
  ilev = ind(leva.eq.levf(l))
  do d = 0,5
     print(vdate)
     ind_a   = ind(date_a.eq.vdate)
     ind_f   = ind(date_f.eq.vdate)
     ;print(ind_00)
     ;print(ind_12)

     dz = zf(ind_f,l,:,:)
     dz = zf(ind_f,l,:,:) - za(ind_a,ilev,:,:)
     dz@_FillValue = default_fillvalue("float")
     printVarSummary(dz)
     ;printMinMax(pres,1)
;---Create plot
     res                      = True               ; plot mods desired
     res@gsnDraw	      = False
     res@gsnFrame	      = False
;     res@cnInfoLabelOn	      = False
     
     res@gsnCenterString	   = sprinti("%0.10d",vdate)
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
     res@mpMinLonF            =  120.               ; select a subregion
     res@mpMaxLonF            =  160.
     res@mpMinLatF            =  10. 
     res@mpMaxLatF            =  60.

     res@lbOrientation            = "Vertical"     ; vertical label bar
     res@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer

; note, when doing a subregion, NCL determines the range of the data from
; the full domain. If you wish to just consider the domain you are plotting,
; you must manually set those levels.

     res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
     res@cnMaxLevelCount	   = 20
     res@cnMinLevelValF       = -50.0               ; set min contour level
     res@cnMaxLevelValF       = 50.0                 ; set max contour level
     res@cnLevelSpacingF      = 5.0               ; set contour spacing

     plotdz = gsn_csm_contour_map_ce(wks,dz,res)
     z = zf(ind_f,l,:,:)
     if (l.eq.0) then
      res@cnMinLevelValF       = 1200.0               ; set min contour level
      res@cnMaxLevelValF       = 1500.0                 ; set max contour level
      res@cnLevelSpacingF      = 30.0               ; set contour spacing
     else if (l.eq.1) then 
      res@cnMinLevelValF       = 5200.0               ; set min contour level
      res@cnMaxLevelValF       = 5800.0                 ; set max contour level
      res@cnLevelSpacingF      = 50.0               ; set contour spacing
     else if (l.eq.2) then
      res@cnMinLevelValF       = 8800.0               ; set min contour level
      res@cnMaxLevelValF       = 9600.0                 ; set max contour level
      res@cnLevelSpacingF      = 80.0               ; set contour spacing
     else if (l.eq.3) then
      res@cnMinLevelValF       = 9600.0               ; set min contour level
      res@cnMaxLevelValF       = 11000.0                 ; set max contour level
      res@cnLevelSpacingF      = 100.0               ; set contour spacing
     else if (l.eq.4) then
      res@cnMinLevelValF       = 11000.0               ; set min contour level
      res@cnMaxLevelValF       = 13000.0                 ; set max contour level
      res@cnLevelSpacingF      = 100.0               ; set contour spacing
     end if
     end if
     end if
     end if
     end if
     res@cnFillOn             = False               ; turn on color for contours
     res@cnLinesOn            = True                ; turn off contour lines
     res@cnLineLabelsOn        = True
     res@cnLineColor          = "blue"
     res@gsnRightString   = ""    
     res@gsnLeftString    = ""
     res@gsnCenterString  = ""
     plotz = gsn_csm_contour(wks,z,res)
     overlay(plotdz,plotz)
     plot(d) = plotdz
     delete(plotdz)
     delete(plotz)

     year  = toint(str_get_cols(sprinti("%0.10d",vdate), 0, 3))
     month = toint(str_get_cols(sprinti("%0.10d",vdate), 4, 5))
     day   = toint(str_get_cols(sprinti("%0.10d",vdate), 6, 7))
     hour  = toint(str_get_cols(sprinti("%0.10d",vdate), 8, 9))
     sec = cd_inv_calendar(year,month,day,hour,0,0, time_f@units, 0)
     newsec = sec
     newsec = sec + 12*60*60
     vdate  = toint(cd_string(newsec, "%Y%N%D%H"))
     delete(dz)
     delete(z)
   end do
   resP	                   = True
   resP@gsnPanelMainString	   = "hgt difference 1009"+hh+"-analysis @"+levf(l)+"hPa"
   gsn_panel(wks,plot,(/2,3/),resP)
   delete(plot)
   end do
end
