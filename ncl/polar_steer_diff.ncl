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
  ds = 0 ;start
  de = 11 ;end
  a = addfile("../../netcdf/tigge/2019/jma/rotate/np_2019100912_mean.nc","r")
  b = addfile("../../netcdf/tigge/2019/jma/rotate/np_2019100900_mean.nc","r")
  t12 = a->time
  t00 = b->time
  lat = a->lat
  lon = a->lon
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
  wks = gsn_open_wks("png","steer_diff_"+sday+"_to_"+eday)             ; send graphics to PNG file
  gsn_define_colormap(wks, "sunshine_9lev")
  plot = new(4*3,"graphic")

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
    day   = cd_string(dsec, "%Y%N%D%H")
    print(day)
    ;do lev = 0,3
;      print(i)
    u12 = a->UGRD(ind12,:2,:,:)
    v12 = a->VGRD(ind12,:2,:,:)
    u00 = b->UGRD(ind00,:2,:,:)
    v00 = b->VGRD(ind00,:2,:,:)
;      print(u@level)
;  t = a->PRES_meansealevel(0,:,:)
;************************************************
; this data only has an missing_value, so we are
; copying this to _FillValue for plotting purposes
;************************************************
    ;12UTC
    assignFillValue(v12,u12)                           ; in contributed.ncl   
    assignFillValue(u12,v12)

    uzav = dim_avg_n(u12,2)
    vzav = dim_avg_n(v12,2)
    printVarSummary(uzav)
    u_ = u12
    u_ = u12 - conform_dims(dimsizes(u12),uzav,(/0,1/))
    v_ = v12
    v_ = v12 - conform_dims(dimsizes(v12),vzav,(/0,1/))
;      printVarSummary(u_)
;      printVarSummary(v_)
    ndim = dimsizes(u_)
    print(ndim)
    do lev=0,ndim(0)-1
      do ilat=0,ndim(1)-1
        do ilon=0,ndim(2)-1
          u_(lev,ilat,ilon) = u_(lev,ilat,ilon)*sig(lev)
          v_(lev,ilat,ilon) = v_(lev,ilat,ilon)*sig(lev)
        end do
      end do
    end do
    upav12 = dim_sum_n(u_,0)
    upav12!0 = "lat"
    upav12!1 = "lon"
    upav12&lat = lat
    upav12&lon = lon
    vpav12 = upav12
    vpav12 = dim_sum_n(v_,0)
    printVarSummary(upav12) 
    delete(uzav)
    delete(vzav)
    delete(u_)
    delete(v_)
    
    ;0UTC
    assignFillValue(v00,u00)                           ; in contributed.ncl   
    assignFillValue(u00,v00)

    uzav = dim_avg_n(u00,2)
    vzav = dim_avg_n(v00,2)
    printVarSummary(uzav)
    u_ = u00
    u_ = u00 - conform_dims(dimsizes(u00),uzav,(/0,1/))
    v_ = v00
    v_ = v00 - conform_dims(dimsizes(v00),vzav,(/0,1/))
;      printVarSummary(u_)
;      printVarSummary(v_)
    ndim = dimsizes(u_)
    print(ndim)
    do lev=0,ndim(0)-1
      do ilat=0,ndim(1)-1
        do ilon=0,ndim(2)-1
          u_(lev,ilat,ilon) = u_(lev,ilat,ilon)*sig(lev)
          v_(lev,ilat,ilon) = v_(lev,ilat,ilon)*sig(lev)
        end do
      end do
    end do
    upav00 = dim_sum_n(u_,0)
    upav00!0 = "lat"
    upav00!1 = "lon"
    upav00&lat = lat
    upav00&lon = lon
    vpav00 = upav00
    vpav00 = dim_sum_n(v_,0)
    printVarSummary(upav00)
    ;12UTC-0UTC
    upav = upav00
    upav = upav12 - upav00
    vpav = vpav00
    vpav = vpav12 - vpav00
    wpav = upav00
    wpav = sqrt(upav^2 + vpav^2)
    ;mask out center
    minlatF = 87.
    ustr = upav
    ustr@_FillValue = -99999
    ustr({:minlatF},:) = ustr@_FillValue
    vstr = vpav
    vstr({:minlatF},:) = ustr@_FillValue
    wstr = wpav
    wstr({:minlatF},:) = ustr@_FillValue
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
      res@mpGridLatSpacingF = 2.

      res@cnFillOn	       = True                  ; color fill
      res@cnLevelSelectionMode = "ManualLevels"
;      if (lev.eq.0) then
;      	res@cnLevels = (/1.,2.,4.,6.,8.,10.,12.,14.,16./)
;	res@cnMinLevelValF = 1.
;	res@cnMaxLevelValF = 16.
;      elseif (lev.eq.1) then
;       res@cnLevels = (/2.,4.,6.,8.,10.,12.,14.,16.,18./)
      res@cnLevelSpacingF   = 1.
	    res@cnMinLevelValF = 1.
	    res@cnMaxLevelValF = 10.
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
    	;res@gsnLeftString       = "Wind Velocity"
      ;else
	    res@gsnCenterString	= day
	    res@gsnLeftString	= ""
      ;end if
; because this is ice data, which has a gap in the tropics, we need to
; explicitly pass the range of the data to plot.  Since we are coloring 
; the vectors, this range should also match the MinLatF above, since the
; range for the colors is chosen over the full data passed, and not the
; map limits.
      plot(i) = gsn_csm_vector_scalar_map_polar(wks,ustr,vstr,wstr,res)
      i = i + 1
    ;end do
  end do
  resP	= True
  resP@gsnFrame	= False
  resP@gsnPanelMainString = "vertical integrated asymmetric wind diff. 0912-0900"
;  resP@gsnPanelLabelBar	  = True
  resP@lbLabelFontHeightF = 0.007
  resP@gsnPanelLeft	  = 0.005
  resP@gsnPanelBottom	  = 0.08
  resP@gsnPanelSave	  = True
  gsn_panel(wks,plot,(/3,4/),resP)

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
