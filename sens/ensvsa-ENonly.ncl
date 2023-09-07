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
  load "config.ncl"
begin

;---Open netCDF file
  a = addfile(ncfile,"r")
;
; Read in Sea Surface Temperature Anomalies
; Read in U and V at 1000 mb [subscript 0]
; Determine the subscripts corresponding to Jan 88
;
  te3d     = a->$EN$(:,:,:,:)
  print(te3d&lev)
  lat_uv   = a->lat
  lon_uv   = a->lon
  plev     = a->lev
  time     = a->time
  lon0 = 120.0
  lon1 = 160.0
  lat0 = 0.0
  lat1 = 60.0
  dp   = dpres_plevel_Wrap(plev, plev(nlev-1), plev(0), 0)
  dpall= 1000.0
  sig  = dble2flt( dp / dpall )
  print(sig)
  te = te3d(:,0,:,:)
  te3d = te3d * conform(te3d,sig,1)
  te = dim_sum_n(te3d,1)
  if ((EN.eq."pe").or.(EN.eq."te")) then 
    peps   = a->peps(:,:,:)
    printMinMax(peps(0,{lat0:lat1},{lon0:lon1}),0)
    te     = te + peps
  end if
  
  wks = gsn_open_wks(dev,out)
  nplots = nd / 4 + 1
  plot = new(nplots,graphic)
  ;cmap = read_colormap_file("spread_15lev")
  cmap = read_colormap_file("gsltod")
  
  ypts = (/ 37.0, 37.0, 33.0, 33.0, 37.0/)
  xpts = (/137.0,142.0,142.0,137.0,137.0/)
  printMinMax(te(0,{33:37},{137:142}),0)
  resp                  = True                      ; polyline mods desired
  resp@gsLineColor      = "black"                     ; color of lines
  resp@gsLineThicknessF = 3.0                       ; thickness of lines
  i = 0;1
  dum = new((/nplots,4/),graphic)  
  cnLevels = (/\
  (/2.,4.,8.,16.,32.,64.,128./),\
  (/4.,8.,16.,32.,64.,128.,256./),\
  (/8.,16.,32.,64.,128.,256.,512./),\
  (/16.,32.,64.,128.,256.,512.,1024./)\
  /)
  do d = 0,nd-1,4
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
     res@gsnLeftString	    = "J/kg"
     res@gsnRightString   = ""    
     res@gsnAddCyclic	      = False
     res@cnFillOn             = True               ; turn on color for contours
     res@cnLinesOn            = False              ; turn off contour lines
     res@cnLineLabelsOn       = False              ; turn off contour line labels
     res@cnFillPalette        = cmap

     ;res@mpLandFillColor      = "gray"            ; set land to be gray
     res@mpFillOn             = False
     res@mpMinLonF            = lon0               ; select a subregion
     res@mpMaxLonF            = lon1
     res@mpMinLatF            = lat0
     res@mpMaxLatF            = lat1

     res@lbOrientation            = "Vertical"     ; vertical label bar
     res@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer

; note, when doing a subregion, NCL determines the range of the data from
; the full domain. If you wish to just consider the domain you are plotting,
; you must manually set those levels.

      ;res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
      res@cnLevelSelectionMode = "ExplicitLevels"     ; set manual contour levels
      res@cnLevels = cnLevels(i,:)
;     res@cnMaxLevelCount	   = 10
;   if (d.eq.0) then
;         res@cnLevels = (/0.5,1.,2.,4.,6.,8.,10./)
;  d=0
   ;   res@cnMinLevelValF       = 1.               ; set min contour level
   ;   res@cnMaxLevelValF       = 10.                 ; set max contour level
   ;   res@cnLevelSpacingF      = 1.               ; set contour spacing
   ;   res@vcRefMagnitudeF      = 5.0
   ;else if (d.le.2) 
;   if (d.le.4) then
;;  d=0
;      ;if (((EN.eq."ke").and.d.eq.0).or.(EN.eq."pe")) then
;      ;   res@cnLevels = (/0.5,1.,2.,4.,8.,16.,32./)
;      ;else ;if ((EN.eq."le").or.(EN.eq."ke")) then
;         res@cnLevels = (/1.,2.,4.,8.,16.,32.,64.,128./)
;      ;;else if (EN.eq."ke") then
;      ;;   res@cnLevels = (/0.25,0.5,1.,2.,4.,8.,16.,32./)
;      ;;else
;      ;;   res@cnLevels = (/0.125,0.25,0.5,1.,2.,4.,8.,16./)
;      ;end if
;      ;end if
;      ;end if
;      ;res@cnMinLevelValF       = 5.               ; set min contour level
;      ;res@cnMaxLevelValF       = 50.                 ; set max contour level
;      ;res@cnLevelSpacingF      = 5.               ; set contour spacing
;      ;res@vcRefMagnitudeF      = 5.0
;   else ;if (d.le.5)
;      ;if ((EN.eq."te").or.(EN.eq."le")) then
;         res@cnLevels = (/8.,16.,32.,64.,128.,256.,512.,1024./)
;      ;else if (EN.eq."ke") then
;      ;   res@cnLevels = (/4.,8.,16.,32.,64.,128.,256.,512./)
;      ;else
;      ;   res@cnLevels = (/2.,4.,8.,16.,32.,64.,128.,256./)
;      ;end if
;      ;end if
;      ;res@cnMinLevelValF       = 20.               ; set min contour level
;      ;res@cnMaxLevelValF       = 200.                 ; set max contour level
;      ;res@cnLevelSpacingF      = 20.               ; set contour spacing
;      ;res@vcRefMagnitudeF      = 10.0
;   ;else
;;  d=6
;   ;   res@cnMinLevelValF       = 50.               ; set min contour level
;   ;   res@cnMaxLevelValF       = 500.                 ; set max contour level
;   ;   res@cnLevelSpacingF      = 50.               ; set contour spacing
;   ;   res@vcRefMagnitudeF      = 10.0
;   ;end if
;   ;end if
;   end if
      plot(i) = gsn_csm_contour_map_ce(wks,t,res)
; draw each line separately. Each line must contain two points.
      do j = 0 , 3
         dum(i,j)=gsn_add_polyline(wks,plot(i),xpts(j:j+1),ypts(j:j+1),resp)      
      end do
      i = i + 1
      delete(res)
   end do

   resP	                   = True
   ;resP@gsnPanelMainString	   = "ensvsa moist TE"
   ;gsn_panel(wks,plot,(/2,nplots/2+1/),resP)
   gsn_panel(wks,plot,(/1,nplots/),resP)
end
