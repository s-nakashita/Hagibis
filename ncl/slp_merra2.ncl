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
  wks = gsn_open_wks("png","slp_0700-1118")
  plot = new(20,graphic)

  id = 0
  
  do dd = 7, 11
  idate = 2019100000 + dd*100
  yyyymmdd = str_get_cols(sprinti("%0.10d",idate), 0, 7)
  yyyymmdd2 = str_get_cols(sprinti("%0.10d",idate-100), 0, 7)
  print(yyyymmdd)
  print(yyyymmdd2)
  a1 = addfile("/Volumes/dandelion/netcdf/nasa/merra2/MERRA2_400.tavg3_3d_asm_Nv."+yyyymmdd+".nc4","r")
  a2 = addfile("/Volumes/dandelion/netcdf/nasa/merra2/MERRA2_400.tavg3_3d_asm_Nv."+yyyymmdd2+".nc4","r")
  ;b1 = addfile("../../netcdf/nasa/merra2/MERRA2_400.tavg3_3d_tdt_Np."+yyyymmdd+".nc4","r")
  ;b2 = addfile("../../netcdf/nasa/merra2/MERRA2_400.tavg3_3d_tdt_Np."+yyyymmdd2+".nc4","r")
  lat     = a1->lat
  lon     = a1->lon
  ;print(lon)
  ;lev     = b1->lev({600:200})
  ;print(lev)
  ;dp      = dpres_plevel_Wrap(lev, lev({600}), lev({200}), 0)
  ;dpall = lev({600})-lev({200})
  ;dp = dp/dpall
  ;print(dp)
  ;time2d  = a->time
  time3d  = a1->time
  time3d  = time3d - conform(time3d, 90, -1)
  ;date2d    = cd_string(time2d, "%N%D %H:%M")
  date3d    = cd_string(time3d, "%N%D %H:%M")
  print(date3d)
  
  do d = 0, 6, 2
;
; Read in Sea Surface Temperature Anomalies
; Read in U and V at 1000 mb [subscript 0]
; Determine the subscripts corresponding to Jan 88
;
     ;d2      = ind(date2d.eq.date3d(d))
     if(d.eq.0)then
          ps1   = a1->SLP(d,:,:)
          ps2   = a2->SLP(7,:,:)
          ps    = ps1
          ps    = (ps1+ps2)/2.0
     else  
          ps1   = a1->SLP(d,:,:)
          ps2   = a1->SLP(d-1,:,:)
          ps    = ps1
          ps    = (ps1+ps2)/2.0
     end if
;*************************************************
; interpolation
;*************************************************
     opt = False
     newnlat = 361
     newnlon = 720
     psnew = f2fsh_Wrap(ps, (/newnlat, newnlon/))
     psnew@long_name = ps@long_name
     delete(ps)
     ps = psnew

     ;ps      = a->SLP(d,:,:);PRMSL_meansealevel(:,180:340,180:360)
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
     res@cnFillOn             = False              ; turn on color for contours
     res@cnFillPalette        = cmap
     ;res@gsnScalarContour     = True               ; (for gsn_csm_vector_scalar_map)draw contours

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

     ;res@cnLevelSelectionMode = "ManualLevels"
     ;res@cnMinLevelValF       = 1000.0               ; set min contour level
     ;res@cnMaxLevelValF       = 1020.0                 ; set max contour level
     ;res@cnLevelSpacingF      = 1.0               ; set contour spacing
     ;;resp@cnFillOn             = False               ; turn on color for contours
     res@cnLinesOn            = True                ; turn off contour lines
     res@cnLevelSelectionMode = "ExplicitLevels"       ; set manual contour levels
     res@cnLevels = (/1008.,1009.,1010.,1011.,1012./)
     res@cnMonoLineColor = False
     res@cnLineColors = (/(/0.0,0.0,1.0/),(/.25,0.0,.75/),\
                        (/0.5,0.0,0.5/),(/.75,0.0,.25/),\
                        (/1.0,0.0,0.0/)/)
     res@cnMonoLineThickness = False
     res@cnLineThicknesses = (/2.,1.,2.,1.,2./)
     res@cnInfoLabelOn         = False             ; don't draw the info label
     res@cnLineLabelsOn        = True             ; don't draw the info label
     res@cnLineLabelInterval = 1
     res@cnMonoLineLabelFontColor = False
     res@cnLineLabelFontColors = (/(/0.0,0.0,1.0/),(/.25,0.0,.75/),\
                        (/0.5,0.0,0.5/),(/.75,0.0,.25/),\
                        (/1.0,0.0,0.0/)/)
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
   ;resP@gsnPanelLabelBar     = True
   gsn_panel(wks,plot,(/5,4/),resP)
   delete(plot)
end
