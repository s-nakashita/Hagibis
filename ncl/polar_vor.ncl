;*************************************************
; polar_8.ncl
;************************************************
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
;************************************************

;************************************************
; Function to retrieve the location of a plot so we can calcurate
; a location for a title on the left side of the plot.
;************************************************
function get_title_position(plot)
local vpy, vpx, vph
begin
  getvalues plot
      "vpYF"		: vpy
      "vpXF"		: vpx
      "vpHeightF"	: vph
  end getvalues
  pos = (/vpx*0.2,vpy-(vph/2.)/)
  return(pos)
end

;************************************************
; Main code
;************************************************
begin
;************************************************
; read in atmospheric data
;************************************************
  yyyymmddhh = "2019100900"
  ;b = addfile("../../netcdf/tigge/2019/jma/rotate/np_"+yyyymmddhh+"_mean.nc","r")
  ;a = addfile("../../netcdf/tigge/2019/jma/rotate/np_"+str_get_cols(yyyymmddhh,4,9)+"_mean_v.nc","r")
  a = addfile("../../netcdf/tigge/2019/jma/rotate/np_glb_"+str_get_cols(yyyymmddhh,4,9)+"_mean_v.nc","r")
  t = a->time
  lat = a->lat
  lon = a->lon
  ;t = todouble(time)
  ;t@units = time@units
  ;print(typeof(t))
  date  = cd_string(t, "%Y%N%D%H")
  ds = ind(date.eq.2019100912) ;start
  de = ind(date.eq.2019101212) ;end
  nd = de - ds + 1
  print(nd)
  year 	     = str_get_cols(yyyymmddhh,0,3)
  month	     = str_get_cols(yyyymmddhh,4,5)
  day	     = str_get_cols(yyyymmddhh,6,7)
  hour	     = str_get_cols(yyyymmddhh,8,9)
  sec = cd_inv_calendar(toint(year),toint(month),toint(day),toint(hour),0,0,t@units,0)
  ssec 	= sec
  ssec  = sec + 6*60*60*ds
  esec 	= sec
  esec	= sec + 6*60*60*de
  sday 	= cd_string(ssec, "%N%D%H")
  eday 	= cd_string(esec, "%N%D%H")
  wks = gsn_open_wks("png","asym_ini"+yyyymmddhh+"_"+sday+"_to_"+eday+"_vor500")             ; send graphics to PNG file
  ;wks = gsn_open_wks("png","full_vor")             ; send graphics to PNG file
  gsn_define_colormap(wks, "GMT_polar")
  plot = new(nd,"graphic")

  sig = new(3,float)
  ;sig(0) = 3./7. + 4./49.
  ;sig(1) = 10./49. + 1./7.
  ;sig(2) = 1./7.
  sig(0) = 7./22.
  sig(1) = 1./2.
  sig(2) = 2./11.
  print(sig)
  i = 0
  do d = ds,de;,2
    dsec = sec
    dsec = sec + 6*60*60*d
    index = ind(t.eq.dsec)
    day   = cd_string(dsec, "%H UTC,%D")
    print(day)
    ;vor = a->var(index,:2,:,:)
    vor = a->var(index,1,:,:)
    ;vor = vor * 10^4
    ;vor!0 = "lat"
    ;vor!1 = "lon"
    ;vor&lat = lat
    ;vor&lon = lon
;************************************************
; this data only has an missing_value, so we are
; copying this to _FillValue for plotting purposes
;************************************************
    
    ;vrav = dim_avg_n(vor,2)
    vrav = dim_avg_n(vor,1)
    vor_ = vor
    ;vor_ = vor - conform_dims(dimsizes(vor),vrav,(/0,1/))
    vor_ = vor - conform_dims(dimsizes(vor),vrav,0)
    ;vor_!1 = "lat"
    ;vor_!2 = "lon"
    ;vor_&lat = lat
    ;vor_&lon = lon
    vrpav = vor_
    ;vor_ = vor_ * conform_dims(dimsizes(vor), sig, 0)
    ;vrpav = dim_sum_n(vor_,0)
    vrpav!0 = "lat"
    vrpav!1 = "lon"
    vrpav&lat = lat
    vrpav&lon = lon

    ur = new(dimsizes(vrpav), "float")
    vr = new(dimsizes(vrpav), "float")
    vr2uvf(vrpav, ur, vr)
    printMinMax(ur, 0)
    ur@long_name = "tangential wind"
    ur@units     = "m/s"
    vr@long_name = "radial wind"
    vr@units     = "m/s"
    copy_VarCoords(vrpav, ur)
    copy_VarCoords(vrpav, vr)
    ;ndim = dimsizes(u_)
    ;print(ndim)
    ;u_ = u_ * conform_dims(dimsizes(u_), sig, 0)
    ;v_ = v_ * conform_dims(dimsizes(v_), sig, 0)
    ;vr_ = vr_ * conform_dims(dimsizes(vr), sig, 0)
    ;do lev=0,ndim(0)-1
    ;  do ilat=0,ndim(1)-1
    ;    do ilon=0,ndim(2)-1
    ;      u_(lev,ilat,ilon) = u_(lev,ilat,ilon)*sig(lev)
    ;      v_(lev,ilat,ilon) = v_(lev,ilat,ilon)*sig(lev)
    ;      w(lev,ilat,ilon) = w(lev,ilat,ilon)*sig(lev)
    ;    end do
    ;  end do
    ;end do
    ;upav = dim_sum_n(u_,0)
    ;vrpav = dim_sum_n(vr,0)
    ;vrpav!0 = "lat"
    ;vrpav!1 = "lon"
    ;vrpav&lat = lat
    ;vrpav&lon = lon
    ;vpav = upav
    ;vpav = dim_sum_n(v_,0)
    ;wpav = upav
    ;wpav = dim_sum_n(w,0)
    ;wpav = sqrt(upav^2 + vpav^2)
    ;printVarSummary(vrpav) 
;      print(uzav)
;************************************************
; create plot
;************************************************

    res            = True                         ; plot mods desired
    res@gsnDraw    = False
    res@gsnFrame   = False
    res@gsnPolar   = "NH"                         ; specify the hemisphere
    res@mpMinLatF  = 82                           ; minimum lat to plot
    res@mpFillOn   = False
    res@mpOutlineOn = False
    ;res@mpGridAndLimbOn = False
    res@mpGridLatSpacingF = 2.
    res@lbLabelBarOn = False

    res@cnFillOn	       = True                  ; color fill
    ;res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLevelSelectionMode = "ManualLevels"
;      if (lev.eq.0) then
;      	res@cnLevels = (/1.,2.,4.,6.,8.,10.,12.,14.,16./)
;	res@cnMinLevelValF = 1.
;	res@cnMaxLevelValF = 16.
;      elseif (lev.eq.1) then
;    res@cnLevels = (/-10.,-8.,-6.,-4.,-2.,0.,2.,4.,6.,8.,10./)
	  res@cnLevelSpacingF   = .25
    res@cnMinLevelValF    = -2.5
	  res@cnMaxLevelValF    = 2.5
;      elseif (lev.eq.2) then
;      	res@cnLevels = (/3.,6.,9.,12.,15.,18.,21.,24.,27.,30./)
;	res@cnMinLevelValF = 3.
;	res@cnMaxLevelValF = 30.
;      else	
;        res@cnLevels = (/6.,8.,10.,12.,14.,16.,18.,20.,22.,24.,26.,28.,30.,32.,34./)
;	res@cnMinLevelValF = 6.
;	res@cnMaxLevelValF = 34.
;      end if
    res@cnLinesOn            = True                 ; no contour lines
;    res@gsnContourZeroLineThicknessF = 2.0
;      res@lbLabelBarOn	       = False
    
    res@gsnSpreadColorEnd     = -3
    res@gsnScalarContour      = True              ; vectors over scalars
    res@vcRefMagnitudeF       = 10.0              ; vec magnitude
    res@vcRefLengthF          = 0.05              ; size of ref vector
    res@vcLineArrowColor      = "black"           ; white vectors
    res@vcMinDistanceF        = 0.04              ; thins arrows near pole  
    res@vcGlyphStyle          = "LineArrow"     ; turns on curly vectors
    res@vcLineArrowThicknessF = 2.
    res@vcMinMagnitudeF       = 2.0

    ;res@vcRefAnnoArrowLineColor   = "black"       ; black ref vector
    ;res@vcRefAnnoArrowUseVecColor = False
    res@gsnCenterStringFontHeightF = 0.05
    res@gsnLeftStringFontHeightF = 0.02
    res@gsnRightStringFontHeightF = 0.02
    ;if(d.eq.ds) then
    ;  res@gsnCenterString     = u_@level+"mb"
    res@gsnLeftString       = "Vorticity"; at "+v@level+"hPa"
    ;else
	  res@gsnCenterString	= day
	  res@gsnRightString	= "10^(-4)/s"
    ;end if
; because this is ice data, which has a gap in the tropics, we need to
; explicitly pass the range of the data to plot.  Since we are coloring 
; the vectors, this range should also match the MinLatF above, since the
; range for the colors is chosen over the full data passed, and not the
; map limits.
    scale = 1e4
    vrpav = vrpav * scale
    plot(i) = gsn_csm_vector_scalar_map_polar(wks,ur,vr,vrpav,res)
    ;plot(i) = gsn_csm_contour_map_polar(wks,vrpav,res)
    i = i + 1
    ;end do
  end do
  resP	= True
  resP@gsnFrame	= False
  resP@gsnPanelMainString = str_get_cols(yyyymmddhh,8,9)+"00 UTC 500hPa"
;  resP@gsnPanelMainString = "vertical integrated asymmetric wind "+str_get_cols(yyyymmddhh,4,9)
  resP@gsnPanelLabelBar	  = True
  resP@lbOrientation      = "Vertical"
  resP@lbLabelFontHeightF = 0.007
  resP@gsnPanelLeft	  = 0.005
  resP@gsnPanelBottom	  = 0.08
  resP@gsnPanelSave	  = True
  gsn_panel(wks,plot,(/3,nd/3/),resP)
  ;gsn_panel(wks,plot,(/1,nd/),resP)

; Add rotated text strings on left side of paneled plots.
; Use the resized plots to determine the position of the titles
;  pos = new(2,float)
;  pos(0:1) = get_title_position(plot(0))
;  pos(2:3) = get_title_position(plot(4))
;  pos(4:5) = get_title_position(plot(8))
;  pos(6:7) = get_title_position(plot(12))
;  txres		      = True
;  txres@txFontHeightF = 0.012
;  txres@txAngleF      = 90.
;  gsn_text_ndc(wks, str_get_cols(yyyymmddhh,8,9)+" UTC", pos(0), pos(1), txres)
;  year = str_get_cols(yyyymmddhh,0,3)
;  month = str_get_cols(yyyymmddhh,4,5)
;  day = str_get_cols(yyyymmddhh,6,7)	
;  hour = str_get_cols(yyyymmddhh,8,9)
;;  print(toint(year))
;  sec = cd_inv_calendar(toint(year),toint(month),toint(day),toint(hour),0,0,t@units,0)
;  do d = ds,de
;     newsec = sec
;     newsec = sec + 6*60*60*d
;     date = cd_string(newsec,"%Y%N%D%H")
;     gsn_text_ndc(wks, date, pos(2*(d-ds)), pos(2*(d-ds)+1), txres)
;  end do    
;;  gsn_text_ndc(wks, "2019100912", pos1(0), pos1(1), txres)
;;  gsn_text_ndc(wks, "2019101000", pos2(0), pos2(1), txres)
;;  gsn_text_ndc(wks, "2019101012", pos3(0), pos3(1), txres)
;;  gsn_text_ndc(wks, "2019101100", pos4(0), pos4(1), txres)

  frame(wks)
end
