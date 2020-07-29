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
  a = addfile("../sens/ensvsa-KE-m1_a.nc","r")
;
; Read in Sea Surface Temperature Anomalies
; Read in U and V at 1000 mb [subscript 0]
; Determine the subscripts corresponding to Jan 88
;
  d = 0
  ug       = a->ugrd(:,:,:,:)
  vg       = a->vgrd(:,:,:,:)
  ke       = a->ke(:,:,:)
  lat_uv   = a->lat
  lon_uv   = a->lon
  plev     = a->lev
  time     = a->time
  
  wks = gsn_open_wks("png","ensvsa_ke+wind_100912_a")
  plot = new(7,graphic)
  cmap = read_colormap_file("spread_15lev")
    
  sig = new(3,float)
  sig(0) = 3./7. + 4./49.
  sig(1) = 10./49. + 1./7.
  sig(2) = 1./7.
  print(sig)
  do d = 0,6
     vdate  = cd_string(time(d), "%Y%N%D%H")
     u = ug(d,:,:,:)
     v = vg(d,:,:,:)
     k = ke(d,:,:)
     ndim = dimsizes(u)
     do ilev=0,ndim(0)-1
      do ilat=0,ndim(1)-1
        do ilon=0,ndim(2)-1
          u(ilev,ilat,ilon) = u(ilev,ilat,ilon)*sig(ilev)
          v(ilev,ilat,ilon) = v(ilev,ilat,ilon)*sig(ilev)
        end do
      end do
     end do
     up = dim_sum_n(u,0)
     up!0 = "lat"
     up!1 = "lon"
     up&lat = lat_uv
     up&lon = lon_uv
     up@_FillValue = default_fillvalue("float")
     vp = up
     vp = dim_sum_n(v,0)
;---Create plot
     res                = True               ; plot mods desired
     res@gsnDraw	      = False
     res@gsnFrame	      = False
;     res@cnInfoLabelOn	      = False
     
     res@gsnCenterString	   = "KE+wind "+vdate
     res@gsnCenterStringOrthogonalPosF = 0.02
     res@gsnLeftString	    = "J/kg/m^2"
     res@gsnRightString   = ""    
     res@gsnAddCyclic	      = False
     res@cnFillOn             = True               ; turn on color for contours
     res@cnLinesOn            = False              ; turn off contour lines
     res@cnLineLabelsOn       = False              ; turn off contour line labels
     res@gsnScalarContour     = True               ; contours desired
     res@cnFillPalette        = cmap

     ;res@mpLandFillColor      = "gray"            ; set land to be gray
     res@mpFillOn             = False
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
;     res@cnMaxLevelCount	   = 10
   if (d.le.3) then
;  d=0
      res@cnMinLevelValF       = 0.               ; set min contour level
      res@cnMaxLevelValF       = 1.                 ; set max contour level
      res@cnLevelSpacingF      = 0.1               ; set contour spacing
      res@vcRefMagnitudeF      = 5.0
   else ;if (d.le.5)
      res@cnMinLevelValF       = 0.               ; set min contour level
      res@cnMaxLevelValF       = 5.                 ; set max contour level
      res@cnLevelSpacingF      = 0.5               ; set contour spacing
      res@vcRefMagnitudeF      = 10.0
;   else
;  d=6
;      res@cnMinLevelValF       = 0.               ; set min contour level
;      res@cnMaxLevelValF       = 10.                 ; set max contour level
;      res@cnLevelSpacingF      = 1.               ; set contour spacing
;      res@vcRefMagnitudeF      = 10.0
;   end if
   end if
      res@vcRefLengthF              = 0.045           ; define length of vec ref
      res@vcRefAnnoOrthogonalPosF   = -1.0            ; move ref vector
      res@vcRefAnnoArrowLineColor   = "black"         ; change ref vector color
      res@vcRefAnnoArrowUseVecColor = False           ; don't use vec color for ref

;     res@vcMinMagnitudeF	      = 1.0		  ; no effect with CurlyVector
      res@vcGlyphStyle            = "CurlyVector"     ; turn on curly vectors
      res@vcLineArrowColor        = "black"           ; change vector color
      res@vcLineArrowThicknessF   = 2.0               ; change vector thickness
      res@vcMinDistanceF	      = 0.02
      res@vcVectorDrawOrder       = "PostDraw"        ; draw vectors last

      wind = sqrt(up^2+vp^2)
      up = where(wind.lt.0.1, up@_FillValue, up)
  
      plot(d) = gsn_csm_vector_scalar_map_ce(wks,up,vp,k,res)

   end do

   resP	                   = True
   resP@gsnPanelMainString	   = "ensvsa-KE KE&wind"
   gsn_panel(wks,plot,(/2,4/),resP)
end
