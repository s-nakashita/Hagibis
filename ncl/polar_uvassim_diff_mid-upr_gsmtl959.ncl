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
  wks = gsn_open_wks("png","uvassim_mid-upr_mgdsst-clim_ini"+yyyymmddhh+"_"+sday+"_to_"+eday)             ; send graphics to PNG file
  gsn_define_colormap(wks, "sunshine_9lev")
  ;gsn_define_colormap(wks, "CBR_wet")
  plot = new(4*3,"graphic")
  level_list = (/500, 300, 250/)

  i = 0
  do d = ds,de
     dsec = sec
     dsec = sec + 3*60*60*d
     index = ind(t.eq.dsec)
     day   = cd_string(dsec, "%Y%N%D%H")
     print(day)
    do lev = 0,2
;      print(i)
      if (level_list(lev) .eq. 1000) then
        u1 = a->U10m(index,:,:)
        v1 = a->V10m(index,:,:)
        u2 = b->U10m(index,:,:)
        v2 = b->V10m(index,:,:)
      else
        u1 = a->U(index,{level_list(lev)},:,:)
        v1 = a->V(index,{level_list(lev)},:,:)
        u2 = b->U(index,{level_list(lev)},:,:)
        v2 = b->V(index,{level_list(lev)},:,:)
      end if
;      print(u@level)
;  t = a->PRES_meansealevel(0,:,:)
;************************************************
; this data only has an missing_value, so we are
; copying this to _FillValue for plotting purposes
;************************************************
      assignFillValue(v1,u1)                           ; in contributed.ncl   
      assignFillValue(u1,v1)
      assignFillValue(v2,u2)                           ; in contributed.ncl   
      assignFillValue(u2,v2)
      u1_ = u1 
      u1_ = u1_ - conform_dims(dimsizes(u1), dim_avg_n(u1,1), 0)
      v1_ = v1
      v1_ = v1_ - conform_dims(dimsizes(v1), dim_avg_n(v1,1), 0)
      assignFillValue(v2,u2)                           ; in contributed.ncl   
      assignFillValue(u2,v2)
      u2_ = u2 
      u2_ = u2_ - conform_dims(dimsizes(u2), dim_avg_n(u2,1), 0)
      v2_ = v2
      v2_ = v2_ - conform_dims(dimsizes(v2), dim_avg_n(v2,1), 0)
      ud = u2_
      ud = u2_ - u1_
      vd = v2_
      vd = v2_ - v1_
      w = ud
      w = sqrt(ud^2 + vd^2)
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
      res@cnLevelSelectionMode = "ManualLevels"
      res@cnMinLevelValF = 1.
      res@cnMaxLevelValF = 10.
      res@cnLevelSpacingF = 1.
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
      res@cnLinesOn            = False                 ; no contour lines
      res@lbLabelBarOn	       = False
      if (d.eq.de) then
        res@lbLabelBarOn = True
      end if
    
      res@gsnSpreadColorEnd     = -3
      res@gsnScalarContour      = True              ; vectors over scalars
      res@vcRefMagnitudeF       = 5.0              ; vec magnitude
      res@vcRefLengthF          = 0.04              ; size of ref vector
      res@vcLineArrowColor      = "black"           ; white vectors
      res@vcMinDistanceF        = 0.03              ; thins arrows near pole  
      res@vcGlyphStyle          = "LineArrow"     ; turns on curly vectors

      res@vcRefAnnoArrowLineColor   = "black"       ; black ref vector
      res@vcRefAnnoArrowUseVecColor = False
      res@gsnCenterStringFontHeightF = 0.05
      res@gsnLeftStringFontHeightF = 0.02
      res@gsnRightStringFontHeightF = 0.02
      if(d.eq.ds) then
        res@gsnCenterString     = ud@level+"hPa"
    	  res@gsnLeftString       = "Wind Velocity"
      else
	      res@gsnCenterString	= ""
	      res@gsnLeftString	= ""
      end if
; because this is ice data, which has a gap in the tropics, we need to
; explicitly pass the range of the data to plot.  Since we are coloring 
; the vectors, this range should also match the MinLatF above, since the
; range for the colors is chosen over the full data passed, and not the
; map limits.
        plot(i) = gsn_csm_vector_scalar_map_polar(wks,ud(:,:),vd(:,:), \
         w(:,:),res)
        i = i + 1
    end do
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
