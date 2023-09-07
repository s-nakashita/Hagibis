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
  dev="ps"
  ;hh = 12
  cmap = read_colormap_file("hotcold_18lev") ; low level
  ;cmap = cmap(::-1,:) ; for dyn
  plev = 500
  idd  = 9
  edd  = 11
  date_range = "10"+sprinti("%0.2d",idd)+"-10"+sprinti("%0.2d",edd)
  wks = gsn_open_wks(dev,"strm_merra2_"+date_range+"_p"+sprinti("%d",plev))
  plot = new(12,graphic)
  lon0  =  120.               ; select a subregion
  lon1  =  200.
  lonc  =  (lon0+lon1)/2
  lat0  =  -10. 
  lat1  =  40.
  id = 0
  do dd = idd, edd
  idate = 2019100000 + dd*100
  yyyymmdd = str_get_cols(sprinti("%0.10d",idate), 0, 7)
  ;a = addfile("/Volumes/dandelion/netcdf/nasa/merra2/MERRA2_400.tavg1_2d_slv_Nx."+yyyymmdd+".nc4","r")
  if (plev.le.500) then
  b = addfile("/Volumes/dandelion/netcdf/nasa/merra2/MERRA2_400.inst3_3d_asm_Np."+yyyymmdd+".nc4","r")
  v_int = False
  else
  b = addfile("/Volumes/dandelion/netcdf/nasa/merra2/MERRA2_400.tavg3_3d_asm_Nv."+yyyymmdd+".nc4","r")
  v_int = True
  end if
  lat     = b->lat
  lon     = b->lon
  lev     = b->lev
  ;print(lev)
  ;time2d  = a->time
  time3d  = b->time
  ;date2d    = cd_string(time2d, "%N%D %H:%M")
  date    = cd_string(time3d, "%N%D %H:%M")
  print(date)
  plotnum = (/0,3,6,9,1,4,7,10,2,5,8,11/)
  do d = 0, 6, 2
     if(v_int) then
     ;; vertical interpolation
          u = b->U(d,:,:,:)
          v = b->V(d,:,:,:)
          sfvp = uv2sfvpF(u,v)
          sf3d = sfvp(0,:,:,:)
          copy_VarCoords(u,sf3d)
          sf3d@long_name = "stream function"
          sf3d@units     = "m^2/s"
          sf_zm = zonalAve(sf3d)
          sf3d = sf3d - conform(sf3d, sf_zm, (/0,1/))
          p   = b->PL(d,:,:,:)
          p = p * 1.0e-2 ; hPa
          p@units = "hPa"
          printVarSummary(p)
          printMinMax(p,0)
          sf = sf3d(0,:,:)
          do j=0,dimsizes(lat)-1
               do i=0,dimsizes(lon)-1
                    ptmp = p(:,j,i)
                    ptmp2 = abs(ptmp - plev)
                    k1 = minind(ptmp2)
                    if (ptmp(k1) .lt. plev) then
                         k2 = k1 + 1
                    else
                         k1 = k1 - 1
                         k2 = k1 + 1
                    end if
                    if ((i.eq.0) .and. (j.eq.0)) then
                         print(ptmp)
                         print(ptmp2)
                         print(k1)
                         print(k2)
                         print(ptmp(k1))
                    end if
                    if (k2.eq.dimsizes(lev)) then
                         sf(j,i) = (/sf3d(k1,j,i)/)
                    else
                         wgt = (ptmp(k2)-plev)/(ptmp(k2)-ptmp(k1))
                         stmp = sf3d(k1,j,i)*wgt + sf3d(k2,j,i)*(1.0-wgt)
                         sf(j,i) = (/stmp/)
                    end if
               end do
          end do
          sf!0 = "lat"
          sf&lat = lat
          sf!1 = "lon"
          sf&lon = lon
          printVarSummary(sf)
          printMinMax(sf,0)
     else
          u = b->U(d,{plev},:,:)
          v = b->V(d,{plev},:,:)
          sfvp = uv2sfvpF(u,v)
          sf = sfvp(0,:,:)
          copy_VarCoords(u,sf)
          sf@long_name = "stream function"
          sf@units     = "m^2/s"
          sf_zm = zonalAve(sf)
          sf = sf - conform(sf, sf_zm, 0)
          printVarSummary(sf)
          printMinMax(sf,True)
          ;vp = sfvp(1,:,:)
          ;copy_VarCoords(u,vp)
          ;vp@long_name = "velocity potential"
          ;vp@units     = "m^2/s"
     end if
;---Create plot
     res                  = True               ; plot mods desired
     res@gsnDraw	      = False
     res@gsnFrame	      = False
;     res@cnInfoLabelOn	      = False
     
     res@gsnCenterString	   = date(d)
     res@gsnRightString     = "10~S~6~N~"+sf@units
     res@gsnLeftString	   = ""
     res@gsnCenterStringOrthogonalPosF = 0.02
     res@gsnAddCyclic	      = True
     ;res@gsnScalarContour     = True               ; (for gsn_csm_vector_scalar_map)draw contours

     res@mpFillOn             = False 
;     res@mpLandFillColor      = "gray"            ; set land to be gray
     res@mpMinLonF            =  lon0               ; select a subregion
     res@mpMaxLonF            =  lon1
     res@mpCenterLonF         =  lonc
     res@mpMinLatF            =  lat0
     res@mpMaxLatF            =  lat1

     ;res@lbLabelBarOn         = False
     ;res@lbOrientation            = "Vertical"     ; vertical label bar
     ;res@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer

; note, when doing a subregion, NCL determines the range of the data from
; the full domain. If you wish to just consider the domain you are plotting,
; you must manually set those levels.

     ;res@cnFillOn             = True               ; turn on color for contours
     res@cnFillOn             = False               ; turn on color for contours
     res@cnLinesOn            = True               ; turn off contour lines
     res@cnLineLabelsOn       = False              ; turn off contour line labels
     res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
     res@cnMinLevelValF       = -16.               ; set min contour level
     res@cnMaxLevelValF       = 16.                ; set max contour level
     res@cnLevelSpacingF      = 2.                 ; set contour spacing
     res@gsnContourNegLineDashPattern = 2
     res@gsnContourZeroLineThicknessF = 4.
     res@cnLineThicknessF = 2.
     res@cnFillPalette        = cmap

     sf = sf * 1.0e-6
     plot(plotnum(id)) = gsn_csm_contour_map_ce(wks,sf,res)
     id = id + 1
     delete(sf)
     ;delete(vp)
   end do
   end do
   resP	                   = True
   ;resP@gsnPanelMainString	   = "SLP + "+dT@long_name+" 600-200 hPa"
   ;resP@gsnPanelLabelBar     = True
   gsn_panel(wks,plot,(/4,3/),resP)
   delete(plot)
end
