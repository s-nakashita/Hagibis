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
  yyyymmddhh = "2019100900"
  Exp = "clim"
  if (Exp .eq. "clim") then
    f = "/Volumes/dandelion/GSMJob/Jobwk_Tl959L100/np_fcst_asia_"+yyyymmddhh+".nc"
  else if (Exp .eq. "est") then
    f = "/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_est/np_fcst_asia_"+yyyymmddhh+".nc"
  else if (Exp .eq. "mgdsst")
    f = "/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_mgdsst/np_fcst_asia_"+yyyymmddhh+".nc"
  end if
  end if
  end if
  a = addfile(f,"r")
  t = a->time
  yyyy 	     = str_get_cols(yyyymmddhh,0,3)
  mm  	     = str_get_cols(yyyymmddhh,4,5)
  dd  	     = str_get_cols(yyyymmddhh,6,7)
  hh  	     = str_get_cols(yyyymmddhh,8,9)
  ;       850    800    700    600    500    400    300
  wgt = (/1./22.,3./22.,2./11.,2./11.,2./11.,2./11.,1./11./)
  level = 500
  ;do k = 0, 5
  ;ds = 5*k ;start
  ;de = min((/5*(k+1)-1,28/)) ;end
  sec = cd_inv_calendar(toint(yyyy),toint(mm),toint(dd),toint(hh),0,0,t@units,0)
  ssec = cd_inv_calendar(2019,10,11,6,0,0,t@units,0)
  sday 	= cd_string(ssec, "%N%D%H")
  ds = ind(t.eq.ssec)
  esec = cd_inv_calendar(2019,10,12,3,0,0,t@units,0)
  eday 	= cd_string(esec, "%N%D%H")
  de = ind(t.eq.esec)
  wks = gsn_open_wks("png","vorasym"+sprinti("%3d",level)+"_"+Exp+"_ini"+yyyymmddhh+"_"+sday+"_to_"+eday)             ; send graphics to PNG file
  gsn_define_colormap(wks, "GMT_polar")
  ;gsn_define_colormap(wks, "CBR_wet")
  plot = new(8,"graphic")

  i = 0
  do d = ds,de
    print(sec)
    dsec = sec
    dsec = sec + 3*60*60*d
    index = ind(t.eq.dsec)
    day   = cd_string(dsec, "%N%D %H:00")
    print(day)
    u = a->U(index,{level},:,:)
    v = a->V(index,{level},:,:)
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
    vor = u 
    vor = uv2vr_cfd(u, v, u&lat, u&lon, 3)
    vor@standard_name = "Relative Vorticity"
    vor@units = "10^(-4)/s"
    zav = dim_avg_n(vor,1)
    vor_ = vor
    vor_ = vor_ - conform_dims(dimsizes(vor), zav, 0) ; lat x lon
    vor_ = vor_ * 1.0e4
    printMinMax(vor_,0)
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
    res@cnMinLevelValF = -1.0
    res@cnMaxLevelValF = 1.0
    res@cnLevelSpacingF = 0.1
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
    plot(i) = gsn_csm_contour_map_polar(wks,vor_,res)
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
  gsn_panel(wks,plot,(/2,4/),resP)

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
