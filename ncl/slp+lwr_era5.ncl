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
  
  id = 0
  idate = 2019100900
  edate = 2019101218
  sday  = str_get_cols(sprinti("%0.10d",idate), 4, 9)
  eday  = str_get_cols(sprinti("%0.10d",edate), 4, 9)
  cmap  = read_colormap_file("spread_15lev")
  wks   = gsn_open_wks("png","slp+lwr_era5_"+sday+"-"+eday)
  plot  = new(4*4,graphic)
  yyyymmdd = str_get_cols(sprinti("%0.10d",idate), 0, 7)
  yyyy     = str_get_cols(sprinti("%0.10d",idate), 0, 3)
  mm       = str_get_cols(sprinti("%0.10d",idate), 4, 5)
  a = addfile("/Volumes/dandelion/netcdf/era5/"+yyyy+"/"+mm+"/msl+lwr.nc","r")
  ;b = addfile("../../netcdf/nasa/merra2/MERRA2_400.tavg1_2d_rad_Nx."+yyyymmdd+".nc4","r")
  lat     = a->latitude
  lon     = a->longitude
  time    = a->time
  date = cd_string(time, "%Y%N%D%H")
  ds = ind(date.eq.sprinti("%0.10d",idate))
  de = ind(date.eq.sprinti("%0.10d",edate))
  print(ds)
  print(de)
  delete(date)
  date    = cd_string(time, "%N%D %H UTC")
  print(date)
  
  do d = ds, de
;
; Read in Sea Surface Temperature Anomalies
; Read in U and V at 1000 mb [subscript 0]
; Determine the subscripts corresponding to Jan 88
;
     ps      = a->msl(d,:,:);PRMSL_meansealevel(:,180:340,180:360)
     ps = ps/100
     ps@units = "hPa"
     lwr     = a->mtnlwrf(d,:,:)
     lwr = -1.0*lwr
;---Create plot
     res                  = True               ; plot mods desired
     res@gsnDraw	      = False
     res@gsnFrame	      = False
;     res@cnInfoLabelOn	      = False
     
     res@gsnCenterString	   = date(d)
     res@gsnCenterStringFontHeightF = 0.02
     res@gsnLeftString	   = ""
     res@gsnCenterStringOrthogonalPosF = 0.01
     res@gsnAddCyclic	      = False
     res@cnFillOn             = True               ; turn on color for contours
     res@cnLinesOn            = False               ; turn off contour lines
     res@cnLineLabelsOn       = False              ; turn off contour line labels
     res@cnFillPalette        = cmap(::-1,:)

     res@mpFillOn             = False 
;     res@mpLandFillColor      = "gray"            ; set land to be gray
     res@mpMinLonF            =  120.               ; select a subregion
     res@mpMaxLonF            =  180.
     res@mpMinLatF            =  -10. 
     res@mpMaxLatF            =  55.

     res@lbLabelBarOn         = False
     ;res@lbOrientation            = "Vertical"     ; vertical label bar
     ;res@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer

; note, when doing a subregion, NCL determines the range of the data from
; the full domain. If you wish to just consider the domain you are plotting,
; you must manually set those levels.

     res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
     ;res@cnMaxLevelCount	   = 20
     res@cnMinLevelValF       = 100.0               ; set min contour level
     res@cnMaxLevelValF       = 280.0                 ; set max contour level
     res@cnLevelSpacingF      = 20.0               ; set contour spacing

     plotdz = gsn_csm_contour_map_ce(wks,lwr,res)
     resp = True
     resp@gsnDraw	      = False
     resp@gsnFrame	      = False
     resp@cnLevelSelectionMode = "ManualLevels"
     resp@cnMinLevelValF       = 1000.0               ; set min contour level
     resp@cnMaxLevelValF       = 1020.0                 ; set max contour level
     resp@cnLevelSpacingF      = 4.0               ; set contour spacing
     resp@cnFillOn             = False               ; turn on color for contours
     resp@cnLinesOn            = True                ; turn off contour lines
     resp@cnLineLabelsOn        = True
     resp@cnLineColor          = "black"
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
   resP	                   = True
   ;resP@gsnPanelMainString	   = "SLP + "+lwr@long_name
   resP@gsnPanelLabelBar     = True
   gsn_panel(wks,plot,(/4,4/),resP)
   delete(plot)
end
