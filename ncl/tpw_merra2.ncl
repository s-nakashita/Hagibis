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
  cmap = read_colormap_file("NCV_bright")
  printVarSummary(cmap)
  newcmap = new((/207,4/), float)
  newcmap(:94,:) = cmap(17:111,:)
  newcmap(95:,:) = cmap(144:,:)
  
  plev = 850
  idd  = 9
  edd  = 11
  date_range = "10"+sprinti("%0.2d",idd)+"-10"+sprinti("%0.2d",edd)
  wks = gsn_open_wks("png","tpw_merra2_"+date_range+"_p"+sprinti("%d",plev))
  plot = new(12,graphic)

  id = 0
  do dd = idd, edd
  idate = 2019100000 + dd*100
  yyyymmdd = str_get_cols(sprinti("%0.10d",idate), 0, 7)
  a = addfile("/Volumes/dandelion/netcdf/nasa/merra2/MERRA2_400.tavg1_2d_slv_Nx."+yyyymmdd+".nc4","r")
  b = addfile("/Volumes/dandelion/netcdf/nasa/merra2/MERRA2_400.inst3_3d_asm_Np."+yyyymmdd+".nc4","r")
  lat     = b->lat
  lon     = b->lon
  lev     = b->lev
  ;print(lev)
  time2d  = a->time
  time3d  = b->time
  date1h  = cd_string(time2d, "%N%D %H:%M")
  date    = cd_string(time3d, "%N%D %H:%M")
  print(date)
  plotnum = (/0,3,6,9,1,4,7,10,2,5,8,11/)
  do d = 0, 6, 2
     u = b->U(d,{plev},:,:)
     v = b->V(d,{plev},:,:)
     tpw = a->TQV(d,:,:)

;---Create plot
     res                  = True               ; plot mods desired
     res@gsnDraw	      = False
     res@gsnFrame	      = False
;     res@cnInfoLabelOn	      = False
     
     res@gsnCenterString	   = date(d)
     res@gsnRightString     = u@units
     res@gsnLeftString	   = tpw@units
     res@gsnCenterStringOrthogonalPosF = 0.02
     res@gsnAddCyclic	      = True
     res@gsnScalarContour     = True               ; (for gsn_csm_vector_scalar_map)draw contours

     res@mpFillOn             = False 
;     res@mpLandFillColor      = "gray"            ; set land to be gray
     res@mpMinLonF            =  120.               ; select a subregion
     res@mpMaxLonF            =  200.
     res@mpCenterLonF         =  160.
     res@mpMinLatF            =  -10. 
     res@mpMaxLatF            =  40.

     ;res@lbLabelBarOn         = False
     ;res@lbOrientation            = "Vertical"     ; vertical label bar
     ;res@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer

; note, when doing a subregion, NCL determines the range of the data from
; the full domain. If you wish to just consider the domain you are plotting,
; you must manually set those levels.

     res@cnFillOn             = True               ; turn on color for contours
     res@cnLinesOn            = False              ; turn off contour lines
     res@cnLineLabelsOn       = False              ; turn off contour line labels
     res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
     res@cnMinLevelValF            = 7.5           ; min level
     res@cnMaxLevelValF            = 72.5          ; max level
     res@cnLevelSpacingF           = 2.5           ; contour interval
     res@cnFillPalette             = newcmap

     res@vcRefLengthF              = 0.045           ; define length of vec ref
     res@vcRefAnnoOrthogonalPosF   = -1.0            ; move ref vector
     res@vcFillArrowsOn       = True
     res@vcMinFracLengthF     = 0.2
     res@vcFillArrowEdgeColor = "black"
     res@vcFillArrowFillColor = "white"
     res@vcMinDistanceF	      = 0.02
     res@vcVectorDrawOrder       = "PostDraw"        ; draw vectors last
     wind = sqrt(u^2+v^2)
     u = where(wind.lt.0.5, u@_FillValue, u)
     plot(plotnum(id)) = gsn_csm_vector_scalar_map_ce(wks, u, v, tpw, res)
     id = id + 1
   end do
   end do
   resP	                   = True
   ;resP@gsnPanelMainString	   = "SLP + "+dT@long_name+" 600-200 hPa"
   ;resP@gsnPanelLabelBar     = True
   gsn_panel(wks,plot,(/4,3/),resP)
   delete(plot)
end
