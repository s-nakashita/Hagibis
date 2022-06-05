;*************************************************
; polar_rain.ncl
; plot rain distribution. North is set to the moving direction.
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
  ds = 16 ;start
  de = 23 ;end
  a = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100/np_fcst_asia_"+yyyymmddhh+".nc","r")
  t = a->time
  year 	     = str_get_cols(yyyymmddhh,0,3)
  month	     = str_get_cols(yyyymmddhh,4,5)
  day	     = str_get_cols(yyyymmddhh,6,7)
  hour	     = str_get_cols(yyyymmddhh,8,9)
  sec = cd_inv_calendar(toint(year),toint(month),toint(day),toint(hour),0,0,t@units,0)
  ;ssec 	= sec
  ;ssec  = sec + 6*60*60*ds
  ;esec 	= sec
  ;esec	= sec + 6*60*60*de
  ssec = t(ds)
  esec = t(de)
  sday 	= cd_string(ssec, "%N%D%H")
  eday 	= cd_string(esec, "%N%D%H")
  ; direction data
  data = readAsciiTable("../pytrack/velocity"+yyyymmddhh+"_gsm_tl959.txt", 6, "float", 0)
  nslot = dimsizes(data(:,0))
  print(nslot)
  dd = toint(data(:,2))
  hh = toint(data(:,3))
  deg = data(:,5)

  wks = gsn_open_wks("png","rain+vshear_clim_ini"+yyyymmddhh+"_"+sday+"_to_"+eday)             ; send graphics to PNG file
  gsn_define_colormap(wks, "MPL_YlGnBu")
  plot = new(4*2,"graphic")

  i = 0
  do d = ds,de
    ;dsec = sec
    ;dsec = sec + 6*60*60*d
    ;index = ind(t.eq.dsec)
    dsec = t(d)
    day   = cd_string(dsec, "%Y%N%D%H")
    print(day)
    rain = a->RAIN(d,:,:)
    u850 = a->U(d,{850},{82:89},:)
    v850 = a->V(d,{850},{82:89},:)
    u200 = a->U(d,{200},{82:89},:)
    v200 = a->V(d,{200},{82:89},:)
    printMinMax(u850,0)
    printMinMax(v850,0)
    printMinMax(u200,0)
    printMinMax(v200,0)
    ;psea = a->PSEA(d,:,:)
    ;if (d .gt. 1) then
    ;  rain = rain - a->RAIN(d-1,:,:)
    ;end if
    printMinMax(rain, 0)
    ; direction
    do j = 0, nslot-1
      jsec = cd_inv_calendar(toint(year),toint(month),dd(j),hh(j),0,0,t@units,0)
      if (jsec .eq. dsec) then
        break
      end if
    end do
    lonax = -deg(j)
    print(lonax)
;************************************************
; calculate vertical wind shear
;************************************************
    u200_ = u200
    u200_ = u200 - conform_dims(dimsizes(u200), dim_avg_n(u200,1), 0)
    u850_ = u850
    u850_ = u850 - conform_dims(dimsizes(u850), dim_avg_n(u850,1), 0)
    v200_ = v200
    v200_ = v200 - conform_dims(dimsizes(v200), dim_avg_n(v200,1), 0)
    v850_ = v850
    v850_ = v850 - conform_dims(dimsizes(v850), dim_avg_n(v850,1), 0)
    ud = u200_
    ud = u200_ - u850_
    vd = v200_
    vd = v200_ - v850_
    uvs = avg(ud)
    vvs = avg(vd)
    print(uvs)
    print(vvs)
;************************************************
; create plot
;************************************************

    res            = True                         ; plot mods desired
    res@gsnDraw    = False
    res@gsnFrame   = False
    res@gsnPolar   = "NH"                         ; specify the hemisphere
    res@mpMinLatF  = 82                           ; minimum lat to plot
    res@mpCenterLonF = lonax
    res@mpFillOn   = False
    res@mpOutlineOn = False
    res@mpGridAndLimbOn = False

    res@cnFillOn	       = True                  ; color fill
    res@cnLevelSelectionMode = "ManualLevels"
;      if (lev.eq.0) then
;      	res@cnLevels = (/1.,2.,4.,6.,8.,10.,12.,14.,16./)
;	res@cnMinLevelValF = 1.
;	res@cnMaxLevelValF = 16.
;      elseif (lev.eq.1) then
;       res@cnLevels = (/2.,4.,6.,8.,10.,12.,14.,16.,18./)
	  res@cnMinLevelValF = 10.
	  res@cnMaxLevelValF = 90.
    res@cnLevelSpacingF = 10.
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
    
    res@gsnCenterStringFontHeightF = 0.04
    res@gsnCenterString = day
    res@gsnLeftStringFontHeightF = 0.02
    res@gsnRightStringFontHeightF = 0.02
    plot(i) = gsn_csm_contour_map_polar(wks, rain(:,:),res)
    i = i + 1
  end do
  resP	= True
  resP@gsnFrame	= False
  resP@gsnPanelMainString = "3 hour Precipitation"
;  resP@gsnPanelLabelBar	  = True
  resP@lbLabelFontHeightF = 0.007
  resP@gsnPanelLeft	  = 0.005
  resP@gsnPanelBottom	  = 0.08
  resP@gsnPanelSave	  = True
  gsn_panel(wks,plot,(/2,4/),resP)

; Add rotated text strings on left side of paneled plots.
; Use the resized plots to determine the position of the titles
  pos = new(2*2,float)
  pos(0:1) = get_title_position(plot(0))
  pos(2:3) = get_title_position(plot(3))
  ;pos(4:5) = get_title_position(plot(6))
  ;pos(6:7) = get_title_position(plot(12))
  txres		      = True
  txres@txFontHeightF = 0.012
  txres@txAngleF      = 90.
  year = str_get_cols(yyyymmddhh,0,3)
  month = str_get_cols(yyyymmddhh,4,5)
  day = str_get_cols(yyyymmddhh,6,7)	
  hour = str_get_cols(yyyymmddhh,8,9)
;  print(toint(year))
  sec = cd_inv_calendar(toint(year),toint(month),toint(day),toint(hour),0,0,t@units,0)
  do d = ds, de, 4
     newsec = sec
     newsec = sec + 6*60*60*d
     date = cd_string(newsec,"%Y%N%D%H")
     ;gsn_text_ndc(wks, date, pos((d-ds)/3*2), pos((d-ds)/3*2+1), txres)
  end do    
;  gsn_text_ndc(wks, "2019100912", pos1(0), pos1(1), txres)
;  gsn_text_ndc(wks, "2019101000", pos2(0), pos2(1), txres)
;  gsn_text_ndc(wks, "2019101012", pos3(0), pos3(1), txres)
;  gsn_text_ndc(wks, "2019101100", pos4(0), pos4(1), txres)

  frame(wks)
end
