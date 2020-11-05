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
  year 	     = str_get_cols(yyyymmddhh,0,3)
  month	     = str_get_cols(yyyymmddhh,4,5)
  day	     = str_get_cols(yyyymmddhh,6,7)
  hour	     = str_get_cols(yyyymmddhh,8,9)
  fu = addfile("../../netcdf/tigge/"+year+"/jma/rotate/np_glb_"+str_get_cols(yyyymmddhh,4,9)+"_mean_ug.nc","r")
  fv = addfile("../../netcdf/tigge/"+year+"/jma/rotate/np_glb_"+str_get_cols(yyyymmddhh,4,9)+"_mean_vg.nc","r")
  b = addfile("../../netcdf/tigge/"+year+"/jma/rotate/np_glb_"+str_get_cols(yyyymmddhh,4,9)+"_mean_v.nc","r")
  t = b->time
  lat = b->lat
  lon = b->lon
  ;t = todouble(time)
  ;t@units = time@units
  ;print(typeof(t))
  date  = cd_string(t, "%Y%N%D%H")
  ds = ind(date.eq.2019100912) ;start
  de = ind(date.eq.2019101206) ;end
  nd = de - ds + 1
  print(nd)
  sec = cd_inv_calendar(toint(year),toint(month),toint(day),toint(hour),0,0,t@units,0)
  ssec 	= sec
  ssec  = sec + 6*60*60*ds
  esec 	= sec
  esec	= sec + 6*60*60*de
  sday 	= cd_string(ssec, "%N%D%H")
  eday 	= cd_string(esec, "%N%D%H")
  wks = gsn_open_wks("png","env_ini"+yyyymmddhh+"_"+sday+"_to_"+eday)             ; send graphics to PNG file
  gsn_define_colormap(wks, "sunshine_9lev")
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
    ;do lev = 0,3
;      print(i)
    u = fu->var(index,:2,:,:)
    v = fv->var(index,:2,:,:)
    vor = b->var(index,:2,:,:)
    ;u = fu->var(index,1,:,:)
    ;v = fv->var(index,1,:,:)
    ;vor = b->var(index,1,:,:)
    ;vr = vor * 10^4
    ;vr!0 = "lat"
    ;vr!1 = "lon"
    ;vr&lat = lat
    ;vr&lon = lon
;      print(u@level)
;  t = a->PRES_meansealevel(0,:,:)
;************************************************
; this data only has an missing_value, so we are
; copying this to _FillValue for plotting purposes
;************************************************
    assignFillValue(v,u)                           ; in contributed.ncl   
    assignFillValue(u,v)

    uzav = dim_avg_n(u,2)
    vzav = dim_avg_n(v,2)
    vrav = dim_avg_n(vor,2)
    u_ = u
    u_ = u - conform_dims(dimsizes(u),uzav,(/0,1/))
    v_ = v
    v_ = v - conform_dims(dimsizes(v),vzav,(/0,1/))
    vor_ = vor
    vor_ = vor - conform_dims(dimsizes(vor),vrav,(/0,1/))
    u_ = u_ * conform_dims(dimsizes(u_), sig, 0)
    v_ = v_ * conform_dims(dimsizes(v_), sig, 0)
    vor_ = vor_ * conform_dims(dimsizes(vor_),sig,0)
    
    upav = dim_sum_n(u_,0)
    upav!0 = "lat"
    upav!1 = "lon"
    upav&lat = lat
    upav&lon = lon
    vpav = upav
    vpav = dim_sum_n(v_,0)
    vrpav = dim_sum_n(vor_,0)
    print(dimsizes(vrpav))
    vrpav!0 = "lat"
    vrpav!1 = "lon"
    vrpav&lat = lat
    vrpav&lon = lon
    vrpav@long_name = "vorticity"
    vrpav@units     = "/s"

    ;uzav = dim_avg_n(u,1)
    ;vzav = dim_avg_n(v,1)
    ;vrav = dim_avg_n(vor,1)
    ;;printVarSummary(uzav)
    ;u_ = u
    ;u_ = u - conform_dims(dimsizes(u),uzav,0)
    ;v_ = v
    ;v_ = v - conform_dims(dimsizes(v),vzav,0)
    ;vor_ = vor
    ;vor_ = vor - conform_dims(dimsizes(vor),vrav,0)
    ;printVarSummary(u_)

    ;upav = new(dimsizes(u_),"float")
    ;upav = (/u_/)
    ;;upav@long_name = u@long_name
    ;upav!0 = "lat"
    ;upav!1 = "lon"
    ;upav&lat = lat
    ;upav&lon = lon
    ;printVarSummary(upav)
    ;vpav = upav
    ;vpav = (/v_/)
    ;vrpav = vor_
    ;print(dimsizes(vrpav))
    ;vrpav!0 = "lat"
    ;vrpav!1 = "lon"
    ;vrpav&lat = lat
    ;vrpav&lon = lon
    ;vrpav@long_name = "vorticity"
    ;vrpav@units     = "/s"

    ; calcurate rotational wind
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

    ; subtract asymmetrical wind - rotational wind
    ua = upav
    ua = upav - ur
    va = vpav
    va = vpav - vr
    wa = ua
    wa = sqrt(ua^2 + va^2)
    printMinMax(wa,0)
    ;ndim = dimsizes(u_)
    ;print(ndim)
    
    ;wpav = upav
    ;wpav = dim_sum_n(w,0)
    ;wpav = sqrt(upav^2 + vpav^2)
    ;printVarSummary(upav) 
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
    res@mpGridAndLimbOn = False
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
	  res@cnLevelSpacingF   = .5
    res@cnMinLevelValF    = 1.
	  res@cnMaxLevelValF    = 5.
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
;      res@lbLabelBarOn	       = False
    
    res@gsnSpreadColorEnd     = -3
    res@gsnScalarContour      = True              ; vectors over scalars
    res@vcRefMagnitudeF       = 1.0              ; vec magnitude
    res@vcRefLengthF          = 0.025             ; size of ref vector
    res@vcLineArrowColor      = "black"           ; white vectors
    res@vcMinDistanceF        = 0.04              ; thins arrows near pole  
    res@vcGlyphStyle          = "LineArrow"     ; turns on curly vectors
    res@vcLineArrowThicknessF = 2.
    res@vcMinMagnitudeF       = .5

    ;res@vcRefAnnoArrowLineColor   = "black"       ; black ref vector
    ;res@vcRefAnnoArrowUseVecColor = False
    res@gsnCenterStringFontHeightF = 0.05
    res@gsnLeftStringFontHeightF = 0.02
    res@gsnRightStringFontHeightF = 0.02
    ;if(d.eq.ds) then
    ;  res@gsnCenterString     = u_@level+"mb"
    res@gsnLeftString       = "Wind speed"; at "+v@level+"hPa"
    ;else
	  res@gsnCenterString	= day
	  res@gsnRightString	= "m/s"
    ;end if
; because this is ice data, which has a gap in the tropics, we need to
; explicitly pass the range of the data to plot.  Since we are coloring 
; the vectors, this range should also match the MinLatF above, since the
; range for the colors is chosen over the full data passed, and not the
; map limits.
    plot(i) = gsn_csm_vector_scalar_map_polar(wks,ua,va,wa,res)
    i = i + 1
    ;end do
  end do
  resP	= True
  resP@gsnFrame	= False
  resP@gsnPanelMainString = str_get_cols(yyyymmddhh,8,9)+"00 UTC, vertical mean"
;  resP@gsnPanelMainString = "vertical integrated asymmetric wind "+str_get_cols(yyyymmddhh,4,9)
  resP@gsnPanelLabelBar	  = True
  resP@lbOrientation      = "Vertical"
  resP@lbLabelFontHeightF = 0.007
  resP@gsnPanelLeft	  = 0.005
  resP@gsnPanelBottom	  = 0.08
  resP@gsnPanelSave	  = True
  gsn_panel(wks,plot,(/3,nd/3/),resP)

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
