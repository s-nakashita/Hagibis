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
  cmap = read_colormap_file("BlueDarkRed18") ; low level
  ;cmap = cmap(::-1,:) ; for dyn
  wks = gsn_open_wks("png","slp+divD250_1100-1218")
  plot = new(8,graphic)

  id = 0
  
  do dd = 11, 12
  idate = 2019100000 + dd*100
  yyyymmdd = str_get_cols(sprinti("%0.10d",idate), 0, 7)
  yyyymmdd2 = str_get_cols(sprinti("%0.10d",idate-100), 0, 7)
  print(yyyymmdd)
  print(yyyymmdd2)
  a1 = addfile("../../netcdf/nasa/merra2/MERRA2_400.tavg3_3d_asm_Nv."+yyyymmdd+".nc4","r")
  a2 = addfile("../../netcdf/nasa/merra2/MERRA2_400.tavg3_3d_asm_Nv."+yyyymmdd2+".nc4","r")
  b1 = addfile("../../netcdf/nasa/merra2/MERRA2_400.tavg3_3d_tdt_Np."+yyyymmdd+".nc4","r")
  b2 = addfile("../../netcdf/nasa/merra2/MERRA2_400.tavg3_3d_tdt_Np."+yyyymmdd2+".nc4","r")
  lat     = a1->lat
  lon     = a1->lon
  lev     = b1->lev({600:200})
  print(lev)
  dp      = dpres_plevel_Wrap(lev, lev({600}), lev({200}), 0)
  dpall = lev({600})-lev({200})
  dp = dp/dpall
  print(dp)
  ;time2d  = a->time
  time3d  = a1->time
  time3d  = time3d - conform(time3d, 90, -1)
  ;date2d    = cd_string(time2d, "%N%D %H:%M")
  date3d    = cd_string(time3d, "%N%D %H:%M")
  print(date3d)
  
  ;delp = a1->PL(0,:,0,0)
  ;print(delp)

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
          u1   = a1->U(d,:,:,:)
          u2   = a2->U(7,:,:,:)
          v1   = a1->V(d,:,:,:)
          v2   = a2->V(7,:,:,:)
          pl1  = a1->PL(d,:,:,:)
          pl2  = a2->PL(7,:,:,:)
          pnew = (/25000.0/)
          up1  = int2p_n(pl1, u1, pnew, 1, 0)
          up2  = int2p_n(pl2, u2, pnew, 1, 0)
          u    = up1(0,:,:)
          u    = (up1(0,:,:)+up2(0,:,:))/2.0
          vp1  = int2p_n(pl1, v1, pnew, 1, 0)
          vp2  = int2p_n(pl2, v2, pnew, 1, 0)
          v    = vp1(0,:,:)
          v    = (vp1(0,:,:)+vp2(0,:,:))/2.0
          ;u1   = a1->U(d,42,:,:)
          ;u2   = a2->U(7,42,:,:)
          ;u    = u1
          ;u    = (u1+u2)/2.0
          ;v1   = a1->V(d,42,:,:)
          ;v2   = a2->V(7,42,:,:)
          ;v    = v1
          ;v    = (v1+v2)/2.0
     else  
          ps1   = a1->SLP(d,:,:)
          ps2   = a1->SLP(d-1,:,:)
          ps    = ps1
          ps    = (ps1+ps2)/2.0
          u1   = a1->U(d,:,:,:)
          u2   = a1->U(d-1,:,:,:)
          v1   = a1->V(d,:,:,:)
          v2   = a1->V(d-1,:,:,:)
          pl1  = a1->PL(d,:,:,:)
          pl2  = a2->PL(d-1,:,:,:)
          pnew = (/25000.0/)
          up1  = int2p_n(pl1, u1, pnew, 1, 0)
          up2  = int2p_n(pl2, u2, pnew, 1, 0)
          u    = up1(0,:,:)
          u    = (up1(0,:,:)+up2(0,:,:))/2.0
          vp1  = int2p_n(pl1, v1, pnew, 1, 0)
          vp2  = int2p_n(pl2, v2, pnew, 1, 0)
          v    = vp1(0,:,:)
          v    = (vp1(0,:,:)+vp2(0,:,:))/2.0
          ;u1   = a1->U(d,42,:,:)
          ;u2   = a1->U(d-1,42,:,:)
          ;u    = u1
          ;u    = (u1+u2)/2.0
          ;v1   = a1->V(d,42,:,:)
          ;v2   = a1->V(d-1,42,:,:)
          ;v    = v1
          ;v    = (v1+v2)/2.0
     end if
     copy_VarCoords(ps, u)
     copy_VarCoords(ps, v)
     print(dimsizes(u))
     ;ps      = a->SLP(d,:,:);PRMSL_meansealevel(:,180:340,180:360)
     ps = ps/100
     ps@units = "hPa"
     ;u      = a->U(d,{200},:,:)
     ;v      = a->V(d,{200},:,:)
     dv     = u
     dv     = uv2dvF(u, v)
     
     ud     = u
     vd     = v
     uvd    = dv2uvF(dv)
     ud     = (/uvd(0,:,:)/)
     vd     = (/uvd(1,:,:)/)
     copy_VarCoords(u, ud)
     copy_VarCoords(v, vd)
     ud@long_name = "Zonal Divergent Wind"
     ud@units     = "m/s" ;u@units
     vd@long_name = "Meridional Divergent Wind"
     vd@units     = "m/s" ;v@units

     if(d.eq.0)then
          dT1   = b1->DTDTMST(d,{600:200},:,:)
          dT2   = b2->DTDTMST(7,{600:200},:,:)
          dTm   = dT1
          dTm   = (dT1+dT2)/2.0
     else  
          dT1   = b1->DTDTMST(d,{600:200},:,:)
          dT2   = b1->DTDTMST(d-1,{600:200},:,:)
          dTm   = dT1
          dTm   = (dT1+dT2)/2.0
     end if 
     dTm   = dTm*conform(dTm, dble2flt(dp), 0)
     dTdt  = ps
     dTdt@standard_name = dTm@standard_name
     dTdt@units = "K/h"
     dTdt  = dim_sum_n(dTm, 0) * 3600.0
;---Create plot
     res                  = True               ; plot mods desired
     res@gsnDraw	      = False
     res@gsnFrame	      = False
;     res@cnInfoLabelOn	      = False
     
     res@gsnCenterString	   = date3d(d)
     res@gsnRightString     = dTdt@units;"1.0e-4 /s"
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

     dv     = dv * 10^4
     res@cnLevelSelectionMode = "ExplicitLevels"     ; set manual contour levels
     ;res@cnLevels = (/0.1,0.2,0.3,0.4,0.5,0.75,1.0,1.5,2.0,2.5,3.0,4.0,5.0/)
     res@cnLevels = (/-2.5,-2.0,-1.5,-1.0,-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5,1.0,1.5,2.0,2.5/)
     ;res@cnMaxLevelCount	   = 20
     ;res@cnMinLevelValF       = 0.5               ; set min contour level
     ;res@cnMaxLevelValF       = 5.0                ; set max contour level
     ;res@cnLevelSpacingF      = 0.5               ; set contour spacing

     res@vcRefMagnitudeF           = 10.0             ; define vector ref mag
     res@vcRefLengthF              = 0.045           ; define length of vec ref
     res@vcRefAnnoOrthogonalPosF   = -1.0            ; move ref vector
     res@vcRefAnnoArrowLineColor   = "black"         ; change ref vector color
     res@vcRefAnnoArrowUseVecColor = False           ; don't use vec color for ref

     res@vcGlyphStyle = "CurlyVector"
     res@vcLineArrowColor        = "black"           ; change vector color
     res@vcLineArrowThicknessF   = 2.0               ; change vector thickness
     res@vcMinDistanceF	      = 0.03
     res@vcVectorDrawOrder       = "PostDraw"        ; draw vectors last

     plotdz = gsn_csm_vector_scalar_map_ce(wks, ud, vd, dTdt, res)
     
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
     ;resp@cnLineColor          = "blue"
     resp@gsnRightString   = ""    
     resp@gsnLeftString    = ""
     resp@gsnCenterString  = ""
     plotz = gsn_csm_contour(wks,ps,resp)
     overlay(plotdz,plotz)
     plot(id) = plotdz
     id = id + 1
     delete(plotdz)
     delete(plotz)
     delete(u)
     delete(v)
   end do
   end do
   resP	                   = True
   ;resP@gsnPanelMainString	   = "SLP + "+dT@long_name+" 600-200 hPa"
   resP@gsnPanelLabelBar     = True
   gsn_panel(wks,plot,(/2,4/),resP)
   delete(plot)
end
