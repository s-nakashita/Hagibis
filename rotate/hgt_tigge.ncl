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
; for regrid
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
  
begin
  yyyymmddhh = 2019100912
  MDH = str_get_cols(sprinti("%0.10d",yyyymmddhh),4,9)
;---Open netCDF file
  a = addfile("../netcdf/jma/"+MDH+"_mean.nc","r")
; z -> 850,500,300,250,200
  b = addfile("../netcdf/reanl/hgt.mon.1981-2010.ltm.nc","r")
; climatology data

; Read in Sea Surface Temperature Anomalies
; Read in U and V at 1000 mb [subscript 0]
; Determine the subscripts corresponding to Jan 88
;
  l = 0
  hgt   = a->HGT(:,l,:,:)
  lat   = a->lat
  lon   = a->lon
  plev	= a->level
  time  = a->time
  date  = cd_string(time, "%Y%N%D%H")
  print(date)
  printVarSummary(hgt)
  
  levc	= b->level
  latc	= b->lat
  lonc	= b->lon

  lind	= ind(levc.eq.plev(l))
  print(lind)
  hgtc	= b->hgt(9,lind,:,:)
  
  level = sprintf("%.0f", plev(l))
  wks = gsn_open_wks("png","../png/"+level+"hgt_JMA_epsw_"+MDH)
  plot = new(8,graphic)
;  cmap = read_colormap_file("BlAqGrYeOrReVi200")
    
  vdate = yyyymmddhh
  do d = 1,7
     print(vdate)
     index   = ind(date.eq.vdate)
     print(index)
;---Regrid
     rhgt = linint2_Wrap(lon,lat,hgt(index,:,:),True,lonc,latc(::-1),0)
     rhgt!0 = "lat"
     rhgt!1 = "lon"
     rhgt&lat = latc(::-1)
     rhgt&lon = lonc
     printVarSummary(rhgt)
;---Calcurate Anomaly
     hgta = rhgt
     hgta = rhgt - hgtc(::-1,:)
;---Create plot
     res                      = True               ; plot mods desired
     res@gsnDraw	      = False
     res@gsnFrame	      = False
;     res@cnInfoLabelOn	      = False

     res@gsnCenterString      = "ft="+vdate
     res@gsnCenterStringOrthogonalPosF = 0.05
     res@gsnLeftString        = " "
     res@gsnAddCyclic	      = False
     res@cnFillOn             = True               ; turn on color for contours
     res@cnFillMode	      = "RasterFill"
     res@cnRasterSmoothingOn  = True
     res@cnLinesOn            = False              ; turn off contour lines
     res@cnLineLabelsOn       = True              ; turn off contour line labels
;     res@gsnScalarContour     = True               ; contours desired
;     res@cnFillPalette        = cmap


     res@mpLandFillColor      = "gray"            ; set land to be gray
     res@mpMinLonF            =  120.               ; select a subregion
     res@mpMaxLonF            =  150.
     res@mpMinLatF            =  20. 
     res@mpMaxLatF            =  50.

     res@lbOrientation            = "Vertical"     ; vertical label bar
     res@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer

; note, when doing a subregion, NCL determines the range of the data from
; the full domain. If you wish to just consider the domain you are plotting,
; you must manually set those levels.

      res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
;     res@cnMaxLevelCount	   = 10
      res@cnMinLevelValF       = -300.0               ; set min contour level
      res@cnMaxLevelValF       = 300                 ; set max contour level
      res@cnLevelSpacingF      = 50.0               ; set contour spacing

      plot(d) = gsn_csm_contour_map(wks,hgta,res)

      year  = toint(str_get_cols(sprinti("%0.10d",vdate), 0, 3))
      month = toint(str_get_cols(sprinti("%0.10d",vdate), 4, 5))
      day   = toint(str_get_cols(sprinti("%0.10d",vdate), 6, 7))
      hour  = toint(str_get_cols(sprinti("%0.10d",vdate), 8, 9))
      sec = cd_inv_calendar(year,month,day,hour,0,0, time@units, 0)
      newsec = sec
      newsec = sec + 12*60*60
      vdate  = toint(cd_string(newsec, "%Y%N%D%H"))
   end do
   resP	                   = True
   resP@gsnPanelMainString	   = "hgt anomaly at "+level+"mb JMA ensemble forecast from "+sprinti("%0.10d",yyyymmddhh)
   gsn_panel(wks,plot,(/3,3/),resP)
end
