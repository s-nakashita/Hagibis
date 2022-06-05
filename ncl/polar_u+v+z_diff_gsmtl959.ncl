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
  a = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100/np_fcst_asia_"+yyyymmddhh+".nc","r")
  b = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_mgdsst/np_fcst_asia_"+yyyymmddhh+".nc","r")
  t = a->time
  year 	     = str_get_cols(yyyymmddhh,0,3)
  month	     = str_get_cols(yyyymmddhh,4,5)
  day	     = str_get_cols(yyyymmddhh,6,7)
  hour	     = str_get_cols(yyyymmddhh,8,9)
  do k = 0, 1
  ds = 4*k + 20 ;start (from 10-11 12:00:00)
  de = min((/4*(k+1)+19,28/)) ;end
  sec = cd_inv_calendar(toint(year),toint(month),toint(day),toint(hour),0,0,t@units,0)
  ssec 	= sec
  ssec  = sec + 3*60*60*ds
  esec 	= sec
  esec	= sec + 3*60*60*de
  sday 	= cd_string(ssec, "%N%D%H")
  eday 	= cd_string(esec, "%N%D%H")
  level = 250
  wks = gsn_open_wks("png","u+v+zassym_"+sprinti("%.3d", level)+"_mgdsst-clim_ini"+yyyymmddhh+"_"+sday+"_to_"+eday)             ; send graphics to PNG file
  gsn_define_colormap(wks, "GMT_polar")
  ;gsn_define_colormap(wks, "CBR_wet")
  plot = new(4*3,"graphic")
  
  i = 0
  do d = ds,de
    dsec = sec
    dsec = sec + 3*60*60*d
    index = ind(t.eq.dsec)
    day   = cd_string(dsec, "%Y%N%D%H")
    print(day)
    if (level .eq. 1000) then
      u1 = a->U10m(index,:,:)
      v1 = a->V10m(index,:,:)
      z1 = a->PSEA(index,:,:)
      u2 = b->U10m(index,:,:)
      v2 = b->V10m(index,:,:)
      z2 = b->PSEA(index,:,:)
    else
      u1 = a->U(index,{level},:,:)
      v1 = a->V(index,{level},:,:)
      z1 = a->Z(index,{level},:,:)
      u2 = b->U(index,{level},:,:)
      v2 = b->V(index,{level},:,:)
      z2 = b->Z(index,{level},:,:)
    end if
    u1_ = u1
    u1_ = u1_ - conform_dims(dimsizes(u1), dim_avg_n(u1,1), 0)
    v1_ = v1
    v1_ = v1_ - conform_dims(dimsizes(v1), dim_avg_n(v1,1), 0)
    z1_ = z1
    z1_ = z1_ - conform_dims(dimsizes(z1), dim_avg_n(z1,1), 0)
    u2_ = u2
    u2_ = u2_ - conform_dims(dimsizes(u2), dim_avg_n(u2,1), 0)
    v2_ = v2
    v2_ = v2_ - conform_dims(dimsizes(v2), dim_avg_n(v2,1), 0)
    z2_ = z2
    z2_ = z2_ - conform_dims(dimsizes(z2), dim_avg_n(z2,1), 0)

    ud = u2_
    ud = u2_ - u1_
    vd = v2_
    vd = v2_ - v1_
    zd = z2_
    zd = z2_ - z1_
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
    res@mpMinLatF  = 62                           ; minimum lat to plot
    res@mpFillOn   = False
    res@mpOutlineOn = False
    res@mpGridAndLimbOn = False

    res@cnFillOn	       = True                  ; color fill
    res@cnLinesOn        = True                  ; no contour lines
    res@lbLabelBarOn	   = False
    if (d.eq.de) then
      res@lbLabelBarOn = True
    end if
    
    res@gsnCenterStringFontHeightF = 0.05
    res@gsnLeftStringFontHeightF = 0.02
    res@gsnRightStringFontHeightF = 0.02
    ; U & V
    resuv = res
    resuv@cnLevelSelectionMode = "ManualLevels"
    resuv@cnMinLevelValF = -10.
    resuv@cnMaxLevelValF = 10.
    resuv@cnLevelSpacingF = 1.

    if(d.eq.ds) then
      resuv@gsnCenterString     = "tangential"
      resuv@gsnLeftString       = ""
    else
	    resuv@gsnCenterString	= ""
	    resuv@gsnLeftString	= ""
    end if
    plot(3*i) = gsn_csm_contour_map_polar(wks,ud(:,:),resuv)
    if(d.eq.ds) then
      resuv@gsnCenterString     = "radial"
      resuv@gsnLeftString       = ""
    else
	    resuv@gsnCenterString	= ""
	    resuv@gsnLeftString	= ""
    end if
    plot(3*i+1) = gsn_csm_contour_map_polar(wks,vd(:,:),resuv)

    ; Z
    resz = res
    resz@cnLevelSelectionMode = "ManualLevels"
    resz@cnMinLevelValF = -30.
    resz@cnMaxLevelValF = 30.
    resz@cnLevelSpacingF = 5.
    if(d.eq.ds) then
      resz@gsnCenterString     = "Geopotential height"
      resz@gsnLeftString       = ""
    else
	    resz@gsnCenterString	= ""
	    resz@gsnLeftString	= ""
    end if
    plot(3*i+2) = gsn_csm_contour_map_polar(wks, zd(:,:), resz)
    i = i + 1
  end do
  resP	= True
  resP@gsnFrame	= False
  resP@gsnPanelMainString = ""
;  resP@gsnPanelLabelBar	  = True
  resP@lbLabelFontHeightF = 0.007
  resP@gsnPanelLeft	  = 0.005
  resP@gsnPanelBottom	  = 0.08
  resP@gsnPanelSave	  = True
  gsn_panel(wks,plot,(/4,3/),resP)

; Add rotated text strings on left side of paneled plots.
; Use the resized plots to determine the position of the titles
  pos = new(4*2,float)
  pos(0:1) = get_title_position(plot(0))
  pos(2:3) = get_title_position(plot(3))
  pos(4:5) = get_title_position(plot(6))
  if (de - ds .gt. 2) then
    pos(6:7) = get_title_position(plot(9))
  end if
  txres		      = True
  txres@txFontHeightF = 0.012
  txres@txAngleF      = 90.
  year = str_get_cols(yyyymmddhh,0,3)
  month = str_get_cols(yyyymmddhh,4,5)
  day = str_get_cols(yyyymmddhh,6,7)	
  hour = str_get_cols(yyyymmddhh,8,9)
;  print(toint(year))
  sec = cd_inv_calendar(toint(year),toint(month),toint(day),toint(hour),0,0,t@units,0)
  do d = ds,de
     newsec = sec
     newsec = sec + 3*60*60*d
     date = cd_string(newsec,"%Y%N%D%H")
     gsn_text_ndc(wks, date, pos(2*(d-ds)), pos(2*(d-ds)+1), txres)
  end do    
;  gsn_text_ndc(wks, "2019100912", pos1(0), pos1(1), txres)
;  gsn_text_ndc(wks, "2019101000", pos2(0), pos2(1), txres)
;  gsn_text_ndc(wks, "2019101012", pos3(0), pos3(1), txres)
;  gsn_text_ndc(wks, "2019101100", pos4(0), pos4(1), txres)

  frame(wks)
  delete(wks)
  end do
end
