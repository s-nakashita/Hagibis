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
  pos = (/vpx*0.75,vpy-(vph/2.)/)
  return(pos)
end

;************************************************
; Main code
;************************************************
begin
;************************************************
; read in atmospheric data
;************************************************
  yyyymmddhh = "2019100912"
  ds = 2 ;start
  de = 11 ;end
  av = addfile("../../netcdf/tigge/2019/jma/rotate/np_glb_100912_mean_v.nc","r")
  bv = addfile("../../netcdf/tigge/2019/jma/rotate/np_glb_100900_mean_v.nc","r")
  t12 = av->time
  t00 = bv->time
  lat = av->lat
  lon = av->lon
  year 	     = str_get_cols(yyyymmddhh,0,3)
  month	     = str_get_cols(yyyymmddhh,4,5)
  day	     = str_get_cols(yyyymmddhh,6,7)
  hour	     = str_get_cols(yyyymmddhh,8,9)
  sec = cd_inv_calendar(toint(year),toint(month),toint(day),toint(hour),0,0,t12@units,0)
  ssec 	= sec
  ssec  = sec + 6*60*60*ds
  esec 	= sec
  esec	= sec + 6*60*60*de
  sday 	= cd_string(ssec, "%N%D%H")
  eday 	= cd_string(esec, "%N%D%H")
  wks = gsn_open_wks("png","asym_0912-0900_"+sday+"_to_"+eday+"_vor500")             ; send graphics to PNG file
  gsn_define_colormap(wks, "GMT_polar")
  nd = de-ds+1
  plot = new(nd,"graphic")

  sig = new(3,float)
  sig(0) = 3./7. + 4./49.
  sig(1) = 10./49. + 1./7.
  sig(2) = 1./7.
  print(sig)
  i = 0
  do d = ds,de
    dsec = sec
    dsec = sec + 6*60*60*d
    ind12 = ind(t12.eq.dsec)
    ind00 = ind(t00.eq.dsec)
    day   = cd_string(dsec, "%H UTC,%D")
    print(day)
    ;do lev = 0,3
;      print(i)
    l = 1
    v12 = av->var(ind12,l,:,:)
    v00 = bv->var(ind00,l,:,:)
    ;v12 = av->var(ind12,:2,:,:)
    ;v00 = bv->var(ind00,:2,:,:)
;      print(u@level)
;  t = a->PRES_meansealevel(0,:,:)
;************************************************
; this data only has an missing_value, so we are
; copying this to _FillValue for plotting purposes
;************************************************
    ;12UTC
    vrav = dim_avg_n(v12,1)
   ; vrav = dim_avg_n(v12,2)
    printVarSummary(vrav)
    vrz12 = v12
    vrz12 = v12 - conform_dims(dimsizes(v12),vrav,0)
    ;vrz12 = v12 - conform_dims(dimsizes(v12),vrav,(/0,1/))
    vrz12!0 = "lat"
    vrz12!1 = "lon"
    vrz12&lat = lat
    vrz12&lon = lon
    printVarSummary(vrz12)
    vrp12 = vrz12
    ;vrz12 = vrz12 * conform_dims(dimsizes(vrz12), sig, 0)
    ;vrp12 = dim_sum_n(vrz12,0) 
    ;vrp12!0 = "lat"
    ;vrp12!1 = "lon"
    ;vrp12&lat = lat
    ;vrp12&lon = lon
    delete(vrav)
    
    ur12 = new(dimsizes(vrp12), "float")
    vr12 = new(dimsizes(vrp12), "float")
    vr2uvf(vrp12, ur12, vr12)
    printMinMax(ur12, 0)
    ur12@long_name = "tangential wind"
    ur12@units     = "m/s"
    vr12@long_name = "radial wind"
    vr12@units     = "m/s"
    copy_VarCoords(vrp12, ur12)
    copy_VarCoords(vrp12, vr12)
  
    ;0UTC
    vrav = dim_avg_n(v00,1)
    ;vrav = dim_avg_n(v00,2)
    printVarSummary(vrav)
    vrz00 = v00
    vrz00 = v00 - conform_dims(dimsizes(v00),vrav,0)
    ;vrz00 = v00 - conform_dims(dimsizes(v00),vrav,(/0,1/))
    vrz00!0 = "lat"
    vrz00!1 = "lon"
    vrz00&lat = lat
    vrz00&lon = lon
    printVarSummary(vrz00)
    vrp00 = vrz00
    ;vrz00 = vrz00 * conform_dims(dimsizes(vrz00), sig, 0)
    ;vrp00 = dim_sum_n(vrz00,0) 
    ;vrp00!0 = "lat"
    ;vrp00!1 = "lon"
    ;vrp00&lat = lat
    ;vrp00&lon = lon
    delete(vrav)
    
    ur00 = new(dimsizes(vrp00), "float")
    vr00 = new(dimsizes(vrp00), "float")
    vr2uvf(vrp00, ur00, vr00)
    printMinMax(ur00, 0)
    ur00@long_name = "tangential wind"
    ur00@units     = "m/s"
    vr00@long_name = "radial wind"
    vr00@units     = "m/s"
    copy_VarCoords(vrp00, ur00)
    copy_VarCoords(vrp00, vr00) 
    
    ;12UTC-0UTC
    ur = ur00
    ur = ur12 - ur00
    vr = vr00
    vr = vr12 - vr00
    vrp = vrp00
    vrp = vrp12 - vrp00
    vrp = vrp *10^4
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
    res@lbLabelBarOn = False

    res@cnFillOn	       = True                  ; color fill
    res@cnLevelSelectionMode = "ManualLevels"
;      if (lev.eq.0) then
;      	res@cnLevels = (/1.,2.,4.,6.,8.,10.,12.,14.,16./)
;	res@cnMinLevelValF = 1.
;	res@cnMaxLevelValF = 16.
;      elseif (lev.eq.1) then
;       res@cnLevels = (/2.,4.,6.,8.,10.,12.,14.,16.,18./)
    res@cnLevelSpacingF   = .25
	  res@cnMinLevelValF = -2.5
	  res@cnMaxLevelValF = 2.5
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
    res@vcRefMagnitudeF       = 10.0              ; vec magnitude
    res@vcRefLengthF          = 0.05              ; size of ref vector
    res@vcLineArrowColor      = "black"           ; white vectors
    res@vcMinDistanceF        = 0.03              ; thins arrows near pole  
    res@vcGlyphStyle          = "LineArrow"     ; turns on curly vectors
    res@vcLineArrowThicknessF = 3.
    res@vcMinMagnitudeF       = 3.0

    res@vcRefAnnoArrowLineColor   = "black"       ; black ref vector
    res@vcRefAnnoArrowUseVecColor = False
    res@gsnCenterStringFontHeightF = 0.05
    res@gsnLeftStringFontHeightF = 0.02
    res@gsnRightStringFontHeightF = 0.02
    ;if(d.eq.ds) then
    ;  res@gsnCenterString     = u_@level+"mb"
    res@gsnLeftString       = "Vorticity"
    ;else
	  res@gsnCenterString	= day
	  res@gsnRightString	= "10^(-4)/s"
    ;end if
; because this is ice data, which has a gap in the tropics, we need to
; explicitly pass the range of the data to plot.  Since we are coloring 
; the vectors, this range should also match the MinLatF above, since the
; range for the colors is chosen over the full data passed, and not the
; map limits.
    plot(i) = gsn_csm_vector_scalar_map_polar(wks,ur,vr,vrp,res)
    i = i + 1
    delete(ur)
    delete(vr)
    delete(vrp)
  end do
  resP	= True
  resP@gsnFrame	= False
  resP@gsnPanelMainString = "500hPa asymmetric vorticity diff. 0912-0900"
  resP@gsnPanelLabelBar	  = True
  resP@lbOrientation      = "Vertical"
  resP@lbLabelFontHeightF = 0.007
  resP@gsnPanelLeft	  = 0.005
  resP@gsnPanelBottom	  = 0.08
  resP@gsnPanelSave	  = True
  gsn_panel(wks,plot,(/4,nd/4/),resP)

; Add rotated text strings on left side of paneled plots.
; Use the resized plots to determine the position of the titles
;  pos = new(4*2,float)
;  pos(0:1) = get_title_position(plot(0))
;  pos(2:3) = get_title_position(plot(4))
;  pos(4:5) = get_title_position(plot(8))
;  pos(6:7) = get_title_position(plot(12))
;  txres		      = True
;  txres@txFontHeightF = 0.012
;  txres@txAngleF      = 90.
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
