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
  cmap = read_colormap_file("BlueDarkRed18")
  wks = gsn_open_wks("png","slp+lwr")
  plot = new(4*4,graphic)
  
  id = 0
  do dd=9, 12
  idate = 2019100000 + dd*100
  yyyymmdd = str_get_cols(sprinti("%0.10d",idate), 0, 7)
  a = addfile("../../netcdf/nasa/merra2/MERRA2_400.tavg1_2d_slv_Nx."+yyyymmdd+".nc4","r")
  b = addfile("../../netcdf/nasa/merra2/MERRA2_400.tavg1_2d_rad_Nx."+yyyymmdd+".nc4","r")
  lat     = a->lat
  lon     = a->lon
  time    = a->time
  date    = cd_string(time, "%N%D %H:%M")
  print(date)
  
  do d = 1, 19, 6
;
; Read in Sea Surface Temperature Anomalies
; Read in U and V at 1000 mb [subscript 0]
; Determine the subscripts corresponding to Jan 88
;
     ps      = a->SLP(d,:,:);PRMSL_meansealevel(:,180:340,180:360)
     ps = ps/100
     ps@units = "hPa"
     lwr     = b->LWTUP(d,:,:)
   
;---Create plot
     res                  = True               ; plot mods desired
     res@gsnDraw	      = False
     res@gsnFrame	      = False
;     res@cnInfoLabelOn	      = False
     
     res@gsnCenterString	   = date(d)
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
     res@mpMinLonF            =  130.               ; select a subregion
     res@mpMaxLonF            =  180.
     res@mpMinLatF            =  -10. 
     res@mpMaxLatF            =  20.

     res@lbLabelBarOn         = False
     ;res@lbOrientation            = "Vertical"     ; vertical label bar
     ;res@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer

; note, when doing a subregion, NCL determines the range of the data from
; the full domain. If you wish to just consider the domain you are plotting,
; you must manually set those levels.

     ;res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
     ;res@cnMaxLevelCount	   = 20
     ;res@cnMinLevelValF       = -500.0               ; set min contour level
     ;res@cnMaxLevelValF       = 500.0                 ; set max contour level
     ;res@cnLevelSpacingF      = 50.0               ; set contour spacing

     plotdz = gsn_csm_contour_map_ce(wks,lwr,res)
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
     resp@cnLineColor          = "blue"
     resp@gsnRightString   = ""    
     resp@gsnLeftString    = ""
     resp@gsnCenterString  = ""
     plotz = gsn_csm_contour(wks,ps,resp)
     overlay(plotdz,plotz)
     plot(id) = plotdz
     id = id + 1
     delete(plotdz)
     delete(plotz)

   end do
   end do
   resP	                   = True
   resP@gsnPanelMainString	   = "SLP + Longwave Flux at top of atmosphere"
   resP@gsnPanelLabelBar     = True
   gsn_panel(wks,plot,(/4,4/),resP)
   delete(plot)
end
