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
  hh = "12"
  a = addfile("../../netcdf/gsm/gl/2019/10/init.nc","r")
  b = addfile("../../netcdf/tigge/2019/jma/1009"+hh+"_mean.nc","r")
;
; Read in Sea Surface Temperature Anomalies
; Read in U and V at 1000 mb [subscript 0]
; Determine the subscripts corresponding to Jan 88
;
;  pres00   = a->HGT(2,:,:,:)
;  z12      = b->HGT(:,:,:,:)
  ua      = a->UGRD(:,:,180:340,180:360)
  va      = a->VGRD(:,:,180:340,180:360)
  za      = a->HGT(:,:,180:340,180:360)
  uf      = b->UGRD(:,:,:,:)
  vf      = b->VGRD(:,:,:,:)
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

  do l = 0,4 ; level
  ds = ind(date_f.eq.2019101100)
  de = ind(date_f.eq.2019101212)
  vdate = 2019100900 + toint(hh)
  year  = toint(str_get_cols(sprinti("%0.10d",vdate), 0, 3))
  month = toint(str_get_cols(sprinti("%0.10d",vdate), 4, 5))
  day   = toint(str_get_cols(sprinti("%0.10d",vdate), 6, 7))
  hour  = toint(str_get_cols(sprinti("%0.10d",vdate), 8, 9))
  sec = cd_inv_calendar(year,month,day,hour,0,0, time_f@units, 0)
  newsec = sec
  newsec = sec + 6*60*60*ds
  vdate  = toint(cd_string(newsec, "%Y%N%D%H"))
  ;print(vdate)
  wks = gsn_open_wks("png","wind_hgt_09"+hh+"-anl_"+levf(l)+"hPa_"+date_f(ds)+"to"+date_f(de))
  plot = new(4,graphic)
  cmap = read_colormap_file("GMT_polar")
    
  ilev = ind(leva.eq.levf(l))
  ;print(ilev)
  do d = 0,3
     print(vdate)
     dd = str_get_cols(sprinti("%0.10d",vdate), 6, 7)
     hh2 = str_get_cols(sprinti("%0.10d",vdate), 8, 9)
     ind_a   = ind(date_a.eq.vdate)
     ind_f   = ind(date_f.eq.vdate)
     ;print(ind_a)
     ;print(ind_f)

     u = uf(ind_f,l,:,:)
     v = vf(ind_f,l,:,:)
     wf = u
     wf = sqrt(u^2+v^2)
     delete(u)
     delete(v)
     u = ua(ind_a,ilev,:,:)
     v = va(ind_a,ilev,:,:)
     u@_FillValue = default_fillvalue("float")
     v@_FillValue = default_fillvalue("float")
     wa = u
     wa = sqrt(u^2+v^2)
     w = wf
     w = w - wa
     z = zf(ind_f,l,:,:)
     z = z - za(ind_a,ilev,:,:)
     ;printVarSummary(w)
     ;printMinMax(pres,1)
;---Create plot
     res                      = True               ; plot mods desired
     res@gsnDraw	      = False
     res@gsnFrame	      = False
;     res@cnInfoLabelOn	      = False
     
     res@gsnCenterString	   = hh2 + "UTC, " + dd
     res@gsnRightString	   = z@units
     res@gsnRightStringFontHeightF = 0.025
     res@gsnRightStringOrthogonalPosF = 0.05
     res@gsnLeftString	   = ""
     res@gsnCenterStringOrthogonalPosF = 0.02
     res@gsnAddCyclic	      = False
     res@cnFillOn             = False               ; turn on color for contours
     res@cnLinesOn            = True               ; turn off contour lines
     res@cnLineLabelsOn       = False              ; turn off contour line labels
     ;res@gsnScalarContour     = True               ; (for gsn_csm_vector_scalar_map)draw contours
     res@gsnContourZeroLineThicknessF = 0.0
     res@gsnContourNegLineDashPattern = 1
     
     
     res@lbOrientation            = "Vertical"     ; vertical label bar
     res@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer

     plotz = gsn_csm_contour(wks, z, res)
; note, when doing a subregion, NCL determines the range of the data from
; the full domain. If you wish to just consider the domain you are plotting,
; you must manually set those levels.

     res@mpFillOn             = False 
;     res@mpLandFillColor      = "gray"            ; set land to be gray
     res@mpMinLonF            =  120.               ; select a subregion
     res@mpMaxLonF            =  180.
     res@mpMinLatF            =  10. 
     res@mpMaxLatF            =  60.

     res@gsnCenterString	   = ""
     res@gsnLeftString	   = w@units
     res@gsnLeftStringFontHeightF = 0.025
     res@gsnLeftStringOrthogonalPosF = 0.05
     res@gsnRightString	   = ""
     res@cnFillOn             = True               ; turn on color for contours
     res@cnLinesOn            = False               ; turn off contour lines
     res@cnLineLabelsOn       = False              ; turn off contour line labels
     ;res@gsnScalarContour     = True               ; (for gsn_csm_vector_scalar_map)draw contours
     res@cnFillPalette        = cmap
     res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
     ;res@cnLevels  = (/-20.,-18.,-16.,-14.,-12.,-10.,-8.,-6.,-4.,-2.,2.,4.,6.,8.,10.,12.,14.,16.,18.,20./)
     ;res@cnMaxLevelCount	     = 10
     res@cnMinLevelValF       = -20.0               ; set min contour level
     res@cnMaxLevelValF       = 20.0                 ; set max contour level
     res@cnLevelSpacingF      = 2.0               ; set contour spacing

     
     plotw = gsn_csm_contour_map_ce(wks,w,res)
     overlay(plotw,plotz)
     plot(d) = plotw

     year  = toint(str_get_cols(sprinti("%0.10d",vdate), 0, 3))
     month = toint(str_get_cols(sprinti("%0.10d",vdate), 4, 5))
     day   = toint(str_get_cols(sprinti("%0.10d",vdate), 6, 7))
     hour  = toint(str_get_cols(sprinti("%0.10d",vdate), 8, 9))
     sec = cd_inv_calendar(year,month,day,hour,0,0, time_f@units, 0)
     newsec = sec
     newsec = sec + 12*60*60
     vdate  = toint(cd_string(newsec, "%Y%N%D%H"))
     delete(u)
     delete(v)
     delete(res)
   end do
   resP	                   = True
   resP@gsnPanelMainString	   = "wind hgt difference 1009"+hh+"-analysis @"+levf(l)+"hPa"
   gsn_panel(wks,plot,(/2,2/),resP)
   delete(plot)
   end do
end
