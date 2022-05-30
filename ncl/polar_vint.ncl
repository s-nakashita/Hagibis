;*************************************************
; polar_8.ncl
;************************************************
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
;  load "calc_e.ncl"
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
  yyyymmddhh = "2019100912"
  ;a = addfile("../../netcdf/tigge/2019/jma/rotate/np_"+yyyymmddhh+"_mean.nc","r")
  a = addfile("../rotate/np_glb_"+yyyymmddhh+"_mean.nc","r")
  ;a = addfile("../../netcdf/tigge/2019/jma/rotate/np_glb_"+yyyymmddhh+"_mean_ug.nc","r")
  ;b = addfile("../../netcdf/tigge/2019/jma/rotate/np_glb_"+yyyymmddhh+"_mean_vg.nc","r")
  t = a->time
  lat = a->lat
  lon = a->lon
  ;lat@units = "degrees_north"
  ;lon@units = "degrees_east"
  print(lat)
  ;lat = a->lat
  ;lon = a->lon
  lev = a->level
  print(lev(2::-1))
  date  = cd_string(t, "%Y%N%D%H")
  ds = ind(date.eq.2019101000) ;start
  de = ind(date.eq.2019101018) ;end
  nd = de - ds + 1
  year 	     = str_get_cols(yyyymmddhh,0,3)
  month	     = str_get_cols(yyyymmddhh,4,5)
  day	     = str_get_cols(yyyymmddhh,6,7)
  hour	     = str_get_cols(yyyymmddhh,8,9)
  sec = cd_inv_calendar(toint(year),toint(month),toint(day),toint(hour),0,0,t@units,0)
  ssec 	= sec
  ssec  = sec + 6*ds
  esec 	= sec
  esec	= sec + 6*de
  sday 	= cd_string(ssec, "%N%D%H")
  eday 	= cd_string(esec, "%N%D%H")
  ;wks = gsn_open_wks("png","asym_ini"+yyyymmddhh+"_"+sday+"_to_"+eday+"_vint")             ; send graphics to PNG file
  wks = gsn_open_wks("png","asym_glb"+yyyymmddhh)             ; send graphics to PNG file
  gsn_define_colormap(wks, "sunshine_9lev")
  plot = new(4,"graphic")

;************************************************
; prepare filtering function
;************************************************
;  latr0 = 77.5
;  indr0 = ind(lat.eq.latr0)
;  r0 = latr0
;  deg2r(latr0, r0)
;  print(r0)
;  r = lat({latr0:})
;  lenr = dimsizes(r)
;  deg2r(lat({latr0:}), r)
;  E = new(lenr(0), float)
;  calc_e(dble2flt(r), dble2flt(r0), E)
;  E = conform(E, 1.0, -1) - E

;************************************************
; prepare vertical weight
;************************************************
  sig = new(3,float)
  ;sig(0) = 3./7. + 4./49.
  ;sig(1) = 10./49. + 1./7.
  ;sig(2) = 1./7.
  sig(0) = 7./22.
  sig(1) = 1./2.
  sig(2) = 2./11.
  print(sig)
  i = 0
  do d = ds,de
    dsec = sec
    dsec = sec + 6*d
    index = ind(t.eq.dsec)
    day   = cd_string(dsec, "%N-%D %H UTC")
    print(day)
    ;do lev = 0,3
;      print(i)
    ;u = a->u(index,2::-1,:,:)
    ;v = a->v(index,2::-1,:,:)
    u = a->u(index,{850:300},:,:)
    v = a->v(index,{850:300},:,:)
    ;u = a->var(index,{850:300},:,:)
    ;v = b->var(index,{850:300},:,:)
    printVarSummary(u)
    print(u&level)
    ;u = a->UGRD(index,:2,:,:)
    ;v = a->VGRD(index,:2,:,:)
    
;  t = a->PRES_meansealevel(0,:,:)
;************************************************
; this data only has an missing_value, so we are
; copying this to _FillValue for plotting purposes
;************************************************
    assignFillValue(v,u)                           ; in contributed.ncl   
    assignFillValue(u,v)

    uzav = dim_avg_n(u,2)
    vzav = dim_avg_n(v,2)
    ;uzav(:,indr0:) = uzav(:,indr0:) * conform_dims(dimsizes(uzav(:,indr0:)), E, 1)
    ;vzav(:,indr0:) = vzav(:,indr0:) * conform_dims(dimsizes(vzav(:,indr0:)), E, 1)
    printVarSummary(uzav)
    u_ = u
    u_ = u - conform_dims(dimsizes(u),uzav,(/0,1/))
    v_ = v
    v_ = v - conform_dims(dimsizes(v),vzav,(/0,1/))
;      printVarSummary(u_)
;      printVarSummary(v_)
    w = u_
    w = sqrt(u_^2 + v_^2)

    ndim = dimsizes(u_)
    print(ndim)
    u_ = u_ * conform_dims(dimsizes(u_), sig, 0)
    v_ = v_ * conform_dims(dimsizes(v_), sig, 0)
    ;do lev=0,ndim(0)-1
    ;  do ilat=0,ndim(1)-1
    ;    do ilon=0,ndim(2)-1
    ;      u_(lev,ilat,ilon) = u_(lev,ilat,ilon)*sig(lev)
    ;      v_(lev,ilat,ilon) = v_(lev,ilat,ilon)*sig(lev)
    ;      w(lev,ilat,ilon) = w(lev,ilat,ilon)*sig(lev)
    ;    end do
    ;  end do
    ;end do
    upav = dim_sum_n(u_,0)
    upav!0 = "lat"
    upav!1 = "lon"
    upav&lat = lat
    upav&lon = lon
    vpav = upav
    vpav = dim_sum_n(v_,0)
    wpav = upav
    ;wpav = dim_sum_n(w,0)
    wpav = sqrt(upav^2 + vpav^2)
    printVarSummary(upav) 
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
    res@cnLevelSelectionMode = "ManualLevels"
;      if (lev.eq.0) then
;      	res@cnLevels = (/1.,2.,4.,6.,8.,10.,12.,14.,16./)
;	res@cnMinLevelValF = 1.
;	res@cnMaxLevelValF = 16.
;      elseif (lev.eq.1) then
;       res@cnLevels = (/2.,4.,6.,8.,10.,12.,14.,16.,18./)
	  res@cnLevelSpacingF   = 2.
    res@cnMinLevelValF    = 2.
	  res@cnMaxLevelValF    = 20.
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
    res@vcRefAnnoFontHeightF  = 0.025
    res@vcRefMagnitudeF       = 10.0              ; vec magnitude
    res@vcRefLengthF          = 0.08               ; size of ref vector
    res@vcLineArrowColor      = "black"           ; white vectors
    res@vcMinDistanceF        = 0.04              ; thins arrows near pole  
    res@vcGlyphStyle          = "LineArrow"     ; turns on curly vectors
    res@vcLineArrowThicknessF = 2.
    res@vcMinMagnitudeF       = 2.0

    res@vcRefAnnoArrowLineColor   = "black"       ; black ref vector
    res@vcRefAnnoArrowUseVecColor = False
    res@gsnCenterStringFontHeightF = 0.05
    res@gsnLeftStringFontHeightF = 0.02
    res@gsnRightStringFontHeightF = 0.02
      ;if(d.eq.ds) then
      ;  res@gsnCenterString     = u_@level+"mb"
    res@gsnLeftString       = "Wind Velocity"
      ;else
	  res@gsnCenterString	= day
	  ;res@gsnCenterString	= ""
	  res@gsnLeftString	= ""
      ;end if
; because this is ice data, which has a gap in the tropics, we need to
; explicitly pass the range of the data to plot.  Since we are coloring 
; the vectors, this range should also match the MinLatF above, since the
; range for the colors is chosen over the full data passed, and not the
; map limits.
    plot(i) = gsn_csm_vector_scalar_map_polar(wks,upav,vpav,wpav,res)
    i = i + 1
    ;end do
  end do
  resP	= True
  resP@gsnFrame	= False
;  resP@gsnPanelMainString = str_get_cols(yyyymmddhh,8,9)+"00 UTC"
;  resP@gsnPanelMainString = "vertical integrated asymmetric wind "+str_get_cols(yyyymmddhh,4,9)
  resP@gsnPanelLabelBar	  = True
  resP@lbOrientation      = "Horizontal"
  resP@lbLabelFontHeightF = 0.015
  resP@gsnPanelLeft	  = 0.005
  resP@gsnPanelBottom	  = 0.08
  resP@gsnPanelSave	  = True
  gsn_panel(wks,plot,(/1,nd/),resP)

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
