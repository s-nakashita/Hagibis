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
  ;Exp = "cntl"
  if (Exp .eq. "cntl") then
    f = "/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est/np_fcst_"+yyyymmddhh+".nc"
  else if (Exp .eq. "prtb1") then
    f = "/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est_0912+p/np_fcst_"+yyyymmddhh+".nc"
  else if (Exp .eq. "prtb2")
    f = "/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est_0912+p2/np_fcst_"+yyyymmddhh+".nc"
  else if (Exp .eq. "prtbn") then
    f = "/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est_0912+pn/np_fcst_"+yyyymmddhh+".nc"
  else if (Exp .eq. "prtbf")
    f = "/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est_0912+pf/np_fcst_"+yyyymmddhh+".nc"
  end if
  end if
  end if
  end if
  end if
  a = addfile(f,"r")
  t = a->time
  yyyy 	     = str_get_cols(yyyymmddhh,0,3)
  mm  	     = str_get_cols(yyyymmddhh,4,5)
  dd  	     = str_get_cols(yyyymmddhh,6,7)
  hh  	     = str_get_cols(yyyymmddhh,8,9)
  ;;       850    800    700    600    500    400    300
  ;wgt = (/1./22.,3./22.,2./11.,2./11.,2./11.,2./11.,1./11./)
  ;       850    700    600    500    400    300
  wgt = (/3./22.,5./22.,2./11.,2./11.,2./11.,1./11./)
  ;do k = 0, 0
  ;ds = 5*k+2 ;start
  ;de = min((/5*(k+1)+1,14/)) ;end
  sec = cd_inv_calendar(toint(yyyy),toint(mm),toint(dd),toint(hh),0,0,t@units,0)
  ;ssec 	= sec
  ;ssec  = sec + 6*60*60*ds
  ;esec 	= sec
  ;esec	= sec + 6*60*60*de
  ssec = cd_inv_calendar(2019,10,10,18,0,0,t@units,0)
  sday 	= cd_string(ssec, "%N%D%H")
  ds = ind(t.eq.ssec)
  esec = cd_inv_calendar(2019,10,11,18,0,0,t@units,0)
  eday 	= cd_string(esec, "%N%D%H")
  de = ind(t.eq.esec)
  nplots = (de - ds) / 2 + 1
  wks = gsn_open_wks("png","uvasym_vint_"+Exp+"_ini"+yyyymmddhh+"_"+sday+"_to_"+eday)             ; send graphics to PNG file
  gsn_define_colormap(wks, "sunshine_9lev")
  ;gsn_define_colormap(wks, "CBR_wet")
  plot = new(nplots,"graphic")

  i = 0
  do d = ds,de,2
    print(sec)
    dsec = sec
    dsec = sec + 3*60*60*d
    ;dsec = sec + 6*60*60*d
    index = ind(t.eq.dsec)
    day   = cd_string(dsec, "%N%D %H:00")
    print(day)
    u = a->U(index,{850:300},:,:)
    v = a->V(index,{850:300},:,:)
;        rh = a->RH(index,{level_list(lev)},:,:)
;        temp = a->T(index,{level_list(lev)},:,:)
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
    u_ = u 
    u_ = u_ - conform_dims(dimsizes(u), uzav, (/0,1/)) ; level x lat x lon
    v_ = v
    v_ = v_ - conform_dims(dimsizes(v), vzav, (/0,1/)) ; level x lat x lon

    u_ = u_ * conform_dims(dimsizes(u_), wgt, 0)
    v_ = v_ * conform_dims(dimsizes(v_), wgt, 0)
    upav = dim_sum_n(u_, 0)
    upav!0 = "lat"
    upav!1 = "lon"
    upav&lat = u&lat
    upav&lon = u&lon
    vpav = upav
    vpav = dim_sum_n(v_, 0)
    w = upav
    w = sqrt(upav^2 + vpav^2)
;;************************************************
;; relative humidity => specific humidity
;;************************************************
;      q = temp
;      q@standard_name = "specific humidity"
;      q@units = "g/kg"
;      q = mixhum_ptrh(conform(temp,int2flt(level_list(lev)),-1), temp, rh, -2) ; [g/kg]
;      printVarSummary(q)
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

    res@cnFillOn	       = True                  ; color fill
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMinLevelValF = 2.
    res@cnMaxLevelValF = 20.
    res@cnLevelSpacingF = 2.
;      res@cnLevelSelectionMode = "ExplicitLevels"
;      if (lev.eq.0) then
;      	res@cnLevels = (/8.,10.,12.,14.,15.,16.,17.,18.,19.,20./)
;;      	res@cnMinLevelValF = 6.
;;      	res@cnMaxLevelValF = 16.
;;        res@cnLevelSpacingF = 1.
;      else if (lev.eq.1) then
;       res@cnLevels = (/6.,8.,9.,10.,11.,12.,13.,14.,15.,16./)
;;      	res@cnMinLevelValF = 1.
;;	      res@cnMaxLevelValF = 10.
;;        res@cnLevelSpacingF = 1.
;      else if (lev.eq.2) then
;       res@cnLevels = (/5.,6.,7.,8.,9.,10.,11.,12.,13.,14./)
;;      	res@cnLevels = (/3.,6.,9.,12.,15.,18.,21.,24.,27.,30./)
;;	res@cnMinLevelValF = 3.
;;	res@cnMaxLevelValF = 30.
;      else	
;        res@cnLevels = (/1.,2.,3.,4.,5.,6.,7.,8.,9.,10./)
;;	res@cnMinLevelValF = 6.
;;	res@cnMaxLevelValF = 34.
;      end if
;      end if
;      end if
    res@cnLinesOn            = True                ; no contour lines
    res@lbLabelBarOn	       = False
    ;if (d.eq.de) then
    ;  res@lbLabelBarOn = True
    ;end if
    
    res@gsnSpreadColorEnd     = -3
    res@gsnScalarContour      = True              ; vectors over scalars
    res@vcRefMagnitudeF       = 10.0              ; vec magnitude
    res@vcRefLengthF          = 0.05              ; size of ref vector
    res@vcLineArrowColor      = "black"           ; white vectors
    res@vcMinDistanceF        = 0.04              ; thins arrows near pole  
    res@vcGlyphStyle          = "LineArrow"     ; turns on curly vectors

    res@vcRefAnnoArrowLineColor   = "black"       ; black ref vector
    res@vcRefAnnoArrowUseVecColor = False
    res@gsnCenterStringFontHeightF = 0.05
    res@gsnLeftStringFontHeightF = 0.02
    res@gsnRightStringFontHeightF = 0.02
    res@gsnCenterString	= day
	  res@gsnLeftString	= ""
; because this is ice data, which has a gap in the tropics, we need to
; explicitly pass the range of the data to plot.  Since we are coloring 
; the vectors, this range should also match the MinLatF above, since the
; range for the colors is chosen over the full data passed, and not the
; map limits.
    plot(i) = gsn_csm_vector_scalar_map_polar(wks,upav,vpav, \
         w(:,:),res)
        ; q(:,:),res)
    i = i + 1
  end do
  resP	= True
  resP@gsnFrame	= False
  resP@gsnPanelMainString = ""
  resP@gsnPanelLabelBar	  = True
  resP@lbLabelFontHeightF = 0.007
  resP@gsnPanelLeft	  = 0.005
  resP@gsnPanelBottom	  = 0.08
  resP@gsnPanelSave	  = True
  gsn_panel(wks,plot,(/1,nplots/),resP)

; Add rotated text strings on left side of paneled plots.
; Use the resized plots to determine the position of the titles
;  pos = new(4*2,float)
;  pos(0:1) = get_title_position(plot(0))
;  pos(2:3) = get_title_position(plot(4))
;  pos(4:5) = get_title_position(plot(8))
;  if (de - ds .gt. 2) then
;    pos(6:7) = get_title_position(plot(12))
;  end if
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
;     newsec = sec + 3*60*60*d
;     date = cd_string(newsec,"%Y%N%D%H")
;     gsn_text_ndc(wks, date, pos(2*(d-ds)), pos(2*(d-ds)+1), txres)
;  end do    
;;  gsn_text_ndc(wks, "2019100912", pos1(0), pos1(1), txres)
;;  gsn_text_ndc(wks, "2019101000", pos2(0), pos2(1), txres)
;;  gsn_text_ndc(wks, "2019101012", pos3(0), pos3(1), txres)
;;  gsn_text_ndc(wks, "2019101100", pos4(0), pos4(1), txres)

  frame(wks)
  delete(wks)
  ;end do
end
