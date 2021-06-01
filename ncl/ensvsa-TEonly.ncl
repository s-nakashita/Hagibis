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
  a = addfile("../sens/"+ncfile,"r")
;
; Read in Sea Surface Temperature Anomalies
; Read in U and V at 1000 mb [subscript 0]
; Determine the subscripts corresponding to Jan 88
;
  te       = a->te(:,:,:)
  lat_uv   = a->lat
  lon_uv   = a->lon
  plev     = a->lev
  time     = a->time
  
  wks = gsn_open_wks(dev,out)
  plot = new(nd+1,graphic)
  cmap = read_colormap_file("spread_15lev")
  
  ypts = (/ 37.0, 37.0, 33.0, 33.0, 37.0/)
  xpts = (/137.0,142.0,142.0,137.0,137.0/)
  resp                  = True                      ; polyline mods desired
  resp@gsLineColor      = "black"                     ; color of lines
  resp@gsLineThicknessF = 3.0                       ; thickness of lines
  i = 0
  do d = 0,nd-1
     vdate  = cd_string(time(d), "%H UTC, %D")
     t = te(d,:,:)
;---Create plot
     res                = True               ; plot mods desired
     res@gsnDraw	      = False
     res@gsnFrame	      = False
;     res@cnInfoLabelOn	      = False
     
     res@gsnCenterString  = vdate
     res@gsnCenterStringOrthogonalPosF = 0.02
     res@gsnCenterStringFontHeightF = 0.03
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
     res@mpMinLatF            =  0. 
     res@mpMaxLatF            =  60.

     res@lbOrientation            = "Vertical"     ; vertical label bar
     res@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer

; note, when doing a subregion, NCL determines the range of the data from
; the full domain. If you wish to just consider the domain you are plotting,
; you must manually set those levels.

      ;res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
      res@cnLevelSelectionMode = "ExplicitLevels"     ; set manual contour levels
;     res@cnMaxLevelCount	   = 10
   ;if (d.eq.0) then
;  d=0
   ;   res@cnMinLevelValF       = 1.               ; set min contour level
   ;   res@cnMaxLevelValF       = 10.                 ; set max contour level
   ;   res@cnLevelSpacingF      = 1.               ; set contour spacing
   ;   res@vcRefMagnitudeF      = 5.0
   ;else if (d.le.2) 
   if (d.le.2) then
   ;if (d.le.6) then
;  d=0
      res@cnLevels = (/4.,8.,16.,32.,64./)
      ;res@cnMinLevelValF       = 5.               ; set min contour level
      ;res@cnMaxLevelValF       = 50.                 ; set max contour level
      ;res@cnLevelSpacingF      = 5.               ; set contour spacing
      ;res@vcRefMagnitudeF      = 5.0
   else ;if (d.le.5)
      res@cnLevels = (/8.,16.,32.,64.,128.,256.,512.,1024./)
      ;res@cnMinLevelValF       = 20.               ; set min contour level
      ;res@cnMaxLevelValF       = 200.                 ; set max contour level
      ;res@cnLevelSpacingF      = 20.               ; set contour spacing
      ;res@vcRefMagnitudeF      = 10.0
   ;else
;  d=6
   ;   res@cnMinLevelValF       = 50.               ; set min contour level
   ;   res@cnMaxLevelValF       = 500.                 ; set max contour level
   ;   res@cnLevelSpacingF      = 50.               ; set contour spacing
   ;   res@vcRefMagnitudeF      = 10.0
   ;end if
   end if
   ;end if
      plot(i) = gsn_csm_contour_map_ce(wks,t,res)
      dum = new(4,graphic)  
; draw each line separately. Each line must contain two points.
      do j = 0 , 3
         dum(j)=gsn_add_polyline(wks,plot(i),xpts(j:j+1),ypts(j:j+1),resp)      
      end do
      i = i + 1
      delete(res)
   end do

   resP	                   = True
   ;resP@gsnPanelMainString	   = "ensvsa moist TE"
   gsn_panel(wks,plot,(/2,(nd+1)/2/),resP)
end
