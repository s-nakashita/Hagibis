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

function np2lonlat(unp,vnp,lon)
local pi, deg2rad
begin
  pi = acos(-1.0)
  deg2rad = pi / 180.0
  ull = unp * cos(lon*deg2rad) - vnp * sin(lon*deg2rad)
  vll = unp * sin(lon*deg2rad) + vnp * cos(lon*deg2rad)
  ull!0 = "lon"
  ull&lon = lon
  vll!0 = "lon"
  vll&lon = lon
  return((/ull,vll/))
end
;************************************************
; Main code
;************************************************
begin
;************************************************
; read in atmospheric data
;************************************************
  yyyymmddhh = "2019100900"
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

  wks = gsn_open_wks("png","vshear_clim_ini"+yyyymmddhh+"_"+sday+"_to_"+eday)             ; send graphics to PNG file
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
    u850 = a->U(d,{850},{82:89},:)
    v850 = a->V(d,{850},{82:89},:)
    u200 = a->U(d,{200},{82:89},:)
    v200 = a->V(d,{200},{82:89},:)
    printMinMax(u850,0)
    printMinMax(v850,0)
    printMinMax(u200,0)
    printMinMax(v200,0)
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
    uvs = dim_avg_n(ud, 0)
    uvs!0 = "lon"
    uvs&lon = dble2flt(u200&lon)
    vvs = dim_avg_n(vd, 0)
    vvs!0 = "lon"
    vvs&lon = dble2flt(v200&lon)
    uv = np2lonlat(uvs,vvs,uvs&lon)
;************************************************
; create plot
;************************************************

    res            = True                         ; plot mods desired
    res@gsnDraw    = False
    res@gsnFrame   = False
    res@gsnCenterStringFontHeightF = 0.04
    res@gsnCenterString = day
    res@gsnLeftStringFontHeightF = 0.02
    res@gsnRightStringFontHeightF = 0.02

    res@trXMaxF = max(uvs&lon)
    res@trXMinF = min(uvs&lon)
    res@xyLineColors = (/"blue","red"/)
    plot(i) = gsn_csm_xy(wks, uvs&lon, uv, res)
    ;overlay(plot(i), plotu)
    i = i + 1
  end do
  resP	= True
  resP@gsnFrame	= False
  resP@gsnPanelMainString = "850-200 hPa wind shear"
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
