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
  ;cmap = read_colormap_file("precip3_16lev") ; all level
  cmap = read_colormap_file("GMT_polar") ; low level
  ;cmap = cmap(::-1,:) ; for dyn
  wks = gsn_open_wks("png","slp_era5_0900-1018")
  plot = new(8,graphic)

  id = 0
  
  do idy = 9, 10
  idate = 2019100000 + idy*100
  yyyy = str_get_cols(sprinti("%0.10d",idate), 0, 3)
  mm   = str_get_cols(sprinti("%0.10d",idate), 4, 5)
  dd   = str_get_cols(sprinti("%0.10d",idate), 6, 7)
  yyyymmdd = str_get_cols(sprinti("%0.10d",idate), 0, 7)
  print(yyyymmdd)
  a1 = addfile("/Volumes/dandelion/netcdf/era5/"+yyyy+"/"+mm+"/"+dd+"/slev_"+dd+".nc", "r")
  ;b1 = addfile("../../netcdf/nasa/merra2/MERRA2_400.tavg3_3d_tdt_Np."+yyyymmdd+".nc4","r")
  ;b2 = addfile("../../netcdf/nasa/merra2/MERRA2_400.tavg3_3d_tdt_Np."+yyyymmdd2+".nc4","r")
  lat     = a1->latitude
  lon     = a1->longitude
  ;print(lon)
  ;lev     = b1->lev({600:200})
  ;print(lev)
  ;dp      = dpres_plevel_Wrap(lev, lev({600}), lev({200}), 0)
  ;dpall = lev({600})-lev({200})
  ;dp = dp/dpall
  ;print(dp)
  ;time2d  = a->time
  time3d  = a1->time
  date3d    = cd_string(time3d, "%N%D %H:%M")
  print(date3d)
  
  do d = 0, 6, 2
;
; Read in Sea Surface Temperature Anomalies
; Read in U and V at 1000 mb [subscript 0]
; Determine the subscripts corresponding to Jan 88
;
     ps = a1->msl(d,:,:)
     ps = ps/100
     ps@units = "hPa"
;---Create plot
     res                  = True               ; plot mods desired
     res@gsnDraw	      = False
     res@gsnFrame	      = False
;     res@cnInfoLabelOn	      = False
     
     ; Set resources for tickmarks
     res@gsnStringFont              = "helvetica-bold"  ; will be used for gsnCenterString later
     res@gsnCenterStringFontHeightF = 0.04              ; ditto
     res@gsnRightString             = ""
     res@gsnCenterString	   = date3d(d)
     ;res@gsnRightString     = dTdt@units;"1.0e-4 /s"
     res@gsnRightString     = ""
     res@gsnLeftString	   = ps@units
     res@gsnCenterStringOrthogonalPosF = 0.04
     res@gsnCenterStringFontHeightF = 0.03
     res@gsnAddCyclic	      = False
     res@cnFillOn             = True               ; turn on color for contours
     res@cnLinesOn            = False              ; turn off contour lines
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

     res@cnLevelSelectionMode = "ManualLevels"
     res@cnMinLevelValF       = 1000.0               ; set min contour level
     res@cnMaxLevelValF       = 1020.0                 ; set max contour level
     res@cnLevelSpacingF      = 1.0               ; set contour spacing
     ;resp@cnFillOn             = False               ; turn on color for contours
     res@cnLinesOn            = True                ; turn off contour lines
     res@cnLineLabelsOn        = True
     ;resp@cnLineColor          = "blue"
     ;resp@gsnRightString   = ""    
     ;resp@gsnLeftString    = ""
     ;resp@gsnCenterString  = ""
     plot(id) = gsn_csm_contour_map_ce(wks,ps,res)
     id = id + 1
     delete(ps)
   end do
   end do
   resP	                   = True
   ;resP@gsnPanelMainString	   = "MERRA-2 slp + divergence at 200hPa"
   resP@gsnPanelLabelBar     = True
   gsn_panel(wks,plot,(/2,4/),resP)
   delete(plot)
end
