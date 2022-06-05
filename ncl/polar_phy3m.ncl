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
  yyyymmddhh = "2019100900"
  ds = 0 ;start
  de = 3 ;end
  a = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_est/np_fcst_phy_"+yyyymmddhh+".nc","r")
  t = a->time
  print(t)
  ;exit()
  year 	     = str_get_cols(yyyymmddhh,0,3)
  month	     = str_get_cols(yyyymmddhh,4,5)
  day	     = str_get_cols(yyyymmddhh,6,7)
  hour	     = str_get_cols(yyyymmddhh,8,9)
  sec = cd_inv_calendar(toint(year),toint(month),toint(day),toint(hour),0,0,t@units,0)
  print(sec)
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

  ;plev = 500
  wks = gsn_open_wks("png","hrp600-200_est_ini"+yyyymmddhh+"_"+sday+"_to_"+eday)             ; send graphics to PNG file
  gsn_define_colormap(wks, "hotcold_18lev")
  plot = new(4*4,"graphic")

  lev = a->level({600:200})
  dp  = dpres_plevel_Wrap(lev, lev(0), lev(dimsizes(lev)-1), 0)
  dp  = dp / (lev(0) - lev(dimsizes(lev)-1))
  print(dp)
  ;exit()
  i = 0
  do d = ds,de
    ;dsec = sec
    ;dsec = sec + 6*60*60*d
    ;index = ind(t.eq.dsec)
    dsec = t(d)
    day   = cd_string(dsec, "%Y%N%D%H")
    print(day)
    hr = a->HR(d,{600:200},:,:)
    hrrs = a->HRRS(d,{600:200},:,:)
    hrrl = a->HRRL(d,{600:200},:,:)
    hrcv = a->HRCV(d,{600:200},:,:)
    hrlc = a->HRLC(d,{600:200},:,:)
    hrvd = a->HRVD(d,{600:200},:,:)
    hrad = a->HRAD(d,{600:200},:,:)
    hrr = hrrs
    hrr = hrrs + hrrl
    hrc = hrcv
    hrc = hrcv + hrlc
    hrd = hrvd
    hrd = hrvd + hrad
    ;vertical integration
    vhr = hr(0,:,:)
    hr  = hr * conform(hr, dble2flt(dp), 0)
    vhr = dim_sum_n(hr, 0)
    vhrr = hrr(0,:,:)
    hrr  = hrr * conform(hrr, dble2flt(dp), 0)
    vhrr = dim_sum_n(hrr, 0)
    vhrc = hrc(0,:,:)
    hrc  = hrc * conform(hrc, dble2flt(dp), 0)
    vhrc = dim_sum_n(hrc, 0)
    vhrd = hrd(0,:,:)
    hrd  = hrd * conform(hrd, dble2flt(dp), 0)
    vhrd = dim_sum_n(hrd, 0)
    printMinMax(vhr,0)
    printMinMax(vhrr,0)
    printMinMax(vhrc,0)
    printMinMax(vhrd,0)
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
    res@gsnRightStringFontHeightF = 0.02
    res@gsnLeftString = ""
    res@gsnCenterString = ""
    if (d.eq.ds) then
      res@gsnCenterString = "total"
    end if
	  res@cnMinLevelValF = -.8
	  res@cnMaxLevelValF = .8
    res@cnLevelSpacingF = .1
    plot(4*i)   = gsn_csm_contour_map_polar(wks, vhr(:,:),res)
    if (d.eq.ds) then
    res@gsnCenterString = "radiation"
    end if
    plot(4*i+1) = gsn_csm_contour_map_polar(wks, vhrr(:,:),res)
    if (d.eq.ds) then
    res@gsnCenterString = "cloud"
    end if
	  res@cnMinLevelValF = -8.
	  res@cnMaxLevelValF = 8.
    res@cnLevelSpacingF = 1.
    plot(4*i+2) = gsn_csm_contour_map_polar(wks, vhrc(:,:),res)
    if (d.eq.ds) then
    res@gsnCenterString = "turbulance+advection"
    end if
    plot(4*i+3) = gsn_csm_contour_map_polar(wks, vhrd(:,:),res)
    i = i + 1
  end do
  resP	= True
  resP@gsnFrame	= False
  resP@gsnPanelMainString = "Heating rate integrated from 600 to 200 hPa"
;  resP@gsnPanelLabelBar	  = True
  resP@lbLabelFontHeightF = 0.007
  resP@gsnPanelLeft	  = 0.005
  resP@gsnPanelBottom	  = 0.08
  resP@gsnPanelSave	  = True
  gsn_panel(wks,plot,(/4,4/),resP)

; Add rotated text strings on left side of paneled plots.
; Use the resized plots to determine the position of the titles
  pos = new(2*4,float)
  pos(0:1) = get_title_position(plot(0))
  pos(2:3) = get_title_position(plot(4))
  pos(4:5) = get_title_position(plot(8))
  pos(6:7) = get_title_position(plot(12))
  txres		      = True
  txres@txFontHeightF = 0.012
  txres@txAngleF      = 90.
  year = str_get_cols(yyyymmddhh,0,3)
  month = str_get_cols(yyyymmddhh,4,5)
  day = str_get_cols(yyyymmddhh,6,7)	
  hour = str_get_cols(yyyymmddhh,8,9)
;  print(toint(year))
  sec = cd_inv_calendar(toint(year),toint(month),toint(day),toint(hour),0,0,t@units,0)
  do d = ds, de
     newsec = t(d)
     date = cd_string(newsec,"%Y%N%D%H")
     gsn_text_ndc(wks, date, pos(2*(d-ds)), pos(2*(d-ds)+1), txres)
  end do    
;  gsn_text_ndc(wks, "2019100912", pos1(0), pos1(1), txres)
;  gsn_text_ndc(wks, "2019101000", pos2(0), pos2(1), txres)
;  gsn_text_ndc(wks, "2019101012", pos3(0), pos3(1), txres)
;  gsn_text_ndc(wks, "2019101100", pos4(0), pos4(1), txres)

  frame(wks)
end
