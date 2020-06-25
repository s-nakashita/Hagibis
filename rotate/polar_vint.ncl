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
  ds = 5 ;start
  de = 8 ;end
  a = addfile("../../netcdf/tigge/2019/jma/rotate/np_"+yyyymmddhh+"_mean.nc","r")
  t = a->time
  lat = a->lat
  lon = a->lon
  year 	     = str_get_cols(yyyymmddhh,0,3)
  month	     = str_get_cols(yyyymmddhh,4,5)
  day	     = str_get_cols(yyyymmddhh,6,7)
  hour	     = str_get_cols(yyyymmddhh,8,9)
  sec = cd_inv_calendar(toint(year),toint(month),toint(day),toint(hour),0,0,t@units,0)
  ssec 	= sec
  ssec  = sec + 6*60*60*ds
  esec 	= sec
  esec	= sec + 6*60*60*de
  sday 	= cd_string(ssec, "%N%D%H")
  eday 	= cd_string(esec, "%N%D%H")
  wks = gsn_open_wks("png","../png/asym_ini"+yyyymmddhh+"_"+sday+"_to_"+eday+"_vint")             ; send graphics to PNG file
  gsn_define_colormap(wks, "matlab_jet")
  plot = new(4,"graphic")

  sig = new(3,float)
  sig(0) = 3./7. + 4./49.
  sig(1) = 10./49. + 1./7.
  sig(2) = 1./7.
  print(sig)
  i = 0
  do d = ds,de
    dsec = sec
    dsec = sec + 6*60*60*d
    index = ind(t.eq.dsec)
    day   = cd_string(dsec, "%Y%N%D%H")
    print(day)
    ;do lev = 0,3
;      print(i)
    u = a->UGRD(index,:2,:,:)
    v = a->VGRD(index,:2,:,:)
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
    do lev=0,ndim(0)-1
      do ilat=0,ndim(1)-1
        do ilon=0,ndim(2)-1
          u_(lev,ilat,ilon) = u_(lev,ilat,ilon)*sig(lev)
          v_(lev,ilat,ilon) = v_(lev,ilat,ilon)*sig(lev)
          w(lev,ilat,ilon) = w(lev,ilat,ilon)*sig(lev)
        end do
      end do
    end do
    upav = dim_sum_n(u_,0)
    upav!0 = "lat"
    upav!1 = "lon"
    upav&lat = lat
    upav&lon = lon
    vpav = upav
    vpav = dim_sum_n(v_,0)
    wpav = upav
    wpav = dim_sum_n(w,0)
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

      res@cnFillOn	       = True                  ; color fill
      res@cnLevelSelectionMode = "ManualLevels"
;      if (lev.eq.0) then
;      	res@cnLevels = (/1.,2.,4.,6.,8.,10.,12.,14.,16./)
;	res@cnMinLevelValF = 1.
;	res@cnMaxLevelValF = 16.
;      elseif (lev.eq.1) then
;       res@cnLevels = (/2.,4.,6.,8.,10.,12.,14.,16.,18./)
	res@cnMinLevelValF = 2.
	res@cnMaxLevelValF = 20.
;      elseif (lev.eq.2) then
;      	res@cnLevels = (/3.,6.,9.,12.,15.,18.,21.,24.,27.,30./)
;	res@cnMinLevelValF = 3.
;	res@cnMaxLevelValF = 30.
;      else	
;        res@cnLevels = (/6.,8.,10.,12.,14.,16.,18.,20.,22.,24.,26.,28.,30.,32.,34./)
;	res@cnMinLevelValF = 6.
;	res@cnMaxLevelValF = 34.
;      end if
      res@cnLinesOn            = False                 ; no contour lines
;      res@lbLabelBarOn	       = False
    
      res@gsnSpreadColorEnd     = -3
      res@gsnScalarContour      = True              ; vectors over scalars
      res@vcRefMagnitudeF       = 20.0              ; vec magnitude
      res@vcRefLengthF          = 0.04              ; size of ref vector
      res@vcLineArrowColor      = "black"           ; white vectors
      res@vcMinDistanceF        = 0.02              ; thins arrows near pole  
      res@vcGlyphStyle          = "LineArrow"     ; turns on curly vectors

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
  resP@gsnPanelMainString = "asymmetric wind"
;  resP@gsnPanelLabelBar	  = True
  resP@lbLabelFontHeightF = 0.007
  resP@gsnPanelLeft	  = 0.005
  resP@gsnPanelBottom	  = 0.08
  resP@gsnPanelSave	  = True
  gsn_panel(wks,plot,(/1,4/),resP)

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
