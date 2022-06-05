;----------------------------------------------------
; tigge_3.ncl
;----------------------------------------------------
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"		; needed for dim_avg_Wrap, dim_variance_Wrap
  load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

begin
  yyyymmddhh = "2019100912"
  orig = "jma"
  yyyy = str_get_cols(yyyymmddhh, 0, 3)
  mm   = str_get_cols(yyyymmddhh, 4, 5)
  dd   = str_get_cols(yyyymmddhh, 6, 7)
  hh   = str_get_cols(yyyymmddhh, 8, 9)
  mmddhh = str_get_cols(yyyymmddhh, 4, 9)
  am = addfile("/Volumes/dandelion/netcdf/tigge/"+yyyy+"/"+orig+"/glb_"+yyyymmddhh+"_mean.nc", "r")
  t = am->time
  lon = am->longitude
  lat = am->latitude
  date = cd_string(t,"%N%D%H")
  d = ind(date.eq."101012")
  ;d = ind(date.eq."100912")
  slp = am->msl(d,:,:)
  slp = slp*0.01
  ;res@gsnCenterString = date(d)
  print(date)
  delete(am)

  wks = gsn_open_wks("png","slp_stamp_"+orig+"_i"+mmddhh+"_"+date(d))                      ; send graphics to PNG file
  ;cmap = read_colormap_file("WhiteBlue")
  ;plot_base = gsn_csm_map(wks,res)               ; create the blank base plot
	
  
  mres = True                                    ; create a resource list for contour plots
  mres@mpMinLatF = -10.                                       ; set the minimum latitude
  mres@mpMaxLatF = 40.                                      ; set the maximum latitude
  mres@mpMinLonF = 110.                                     ; set the minimum longitude
  mres@mpMaxLonF = 180.                                     ; set the maximum longitude
  mres@mpCenterLonF = (mres@mpMaxLonF + mres@mpMinLonF)/2.    ; set the longitude at the center of the plot
  ;res@mpOutlineOn = True                                   ; turn the map outlines on
  ;res@mpOutlineBoundarySets = "National"                   ; draw national boundaries
  ;res@mpOutlineSpecifiers = (/"Canada : Provinces","United States : States"/)  ; draw US States, Canadian provinces
  mres@mpLandFillColor = "gray90"			            ; shade land with gray90
  mres@cnLevelSelectionMode = "ExplicitLevels"       ; set manual contour levels
  mres@cnLevels = (/1008.,1009.,1010.,1011.,1012./)
  mres@cnMonoLineColor = False
  mres@cnLineColors = (/(/0.0,0.0,1.0/),(/.25,0.0,.75/),\
                        (/0.5,0.0,0.5/),(/.75,0.0,.25/),\
                        (/1.0,0.0,0.0/)/)
  mres@cnMonoLineThickness = False
  mres@cnLineThicknesses = (/2.,1.,2.,1.,2./)
  ;mres@cnFillOn = True
  ;mres@cnFillPalette = cmap(::-1,:)
  ;mres@cnLevelSelectionMode = "ManualLevels"     ; manually select the levels
  ;mres@cnMinLevelValF  = 1008.                    ; set the minimum contour level
  ;mres@cnMaxLevelValF  = 1016.                    ; set the maximum contour level
  ;mres@cnLevelSpacingF = .5                      ; set the contour interval
  ;mres@cnLineColor = "blue"                      ; contour the lines blue
  mres@gsnFrame = False                          ; don't advance the frame
  mres@gsnDraw  = False                          ; don't draw the plot
  ;mres@cnLineLabelInterval = 4                   ; label every other line
  mres@cnInfoLabelOn         = False             ; don't draw the info label
  mres@cnLineLabelsOn        = True             ; don't draw the info label
  mres@cnLineLabelInterval = 1
  mres@cnMonoLineLabelFontColor = False
  mres@cnLineLabelFontColors = (/(/0.0,0.0,1.0/),(/.25,0.0,.75/),\
                        (/0.5,0.0,0.5/),(/.75,0.0,.25/),\
                        (/1.0,0.0,0.0/)/)
  ;mres@cnLineLabelPerimOn = False                ; turn the line label perimeter off
;  mres@cnLineLabelBackgroundColor = -1
  ;mres@cnLineLabelFontColor = "blue"             ; color the line labels blue
  ;mres@cnLineLabelFontHeightF = 0.005             ; set the line label font height
  ;mres@tiMainString      = date(d)
  ;mres@tiMainFontHeightF = 0.012                 ; set the main title font height
	
  mres@gsnLeftString = "hPa"
  mres@gsnRightString = ""

  polyres = True
  polyres@sMarkerIndex = 16
  polyres@gsMarkerSizeF = 5.0
  polyres@gsMarkerColor = (/"white"/)
	
  plot = new(30,"graphic") 
  dum = new(30,graphic)        
  mres@cnLineColor = "blue"                      ; contour the lines blue
  mres@gsnCenterString = "mean"
  plot_temp = gsn_csm_contour_map(wks,slp,mres)
  llon = lon({150:180})
  ;llon = lon({mres@mpMinLonF:mres@mpMaxLonF})
  llat = lat({-10:20})
  ;llat = lat({mres@mpMinLatF:mres@mpMaxLatF})
  print(llon)
  print(llat)
  ;lslp = slp({mres@mpMinLatF:mres@mpMaxLatF},{mres@mpMinLonF:mres@mpMaxLonF})
  lslp = slp({-10:20},{150:180})
  locmin = local_min(lslp, False, 0.0)
  if .not. any(ismissing(locmin@xi)) then
    x = llon(locmin@xi)
    y = llat(locmin@yi)
    print(x)
    dum(0) = gsn_add_polymarker(wks,plot_temp,x,y,polyres)
    delete(locmin)
    delete(x)
    delete(y)
  end if
  plot(0) = plot_temp
  delete(plot_temp)
  mres@gsnLeftString = ""
  mres@cnLineLabelsOn        = False             ; don't draw the info label
  j = 4
  do i = 1, 26
    if (i .lt. 10) then
      nm = "0"+sprinti("%d",i)
    else 
      nm = sprinti("%2d",i)
    end if
    ;print(nm)
    mres@gsnCenterString = nm
    a = addfile("/Volumes/dandelion/netcdf/tigge/"+yyyy+"/"+orig+"/glb_"+yyyymmddhh+"_"+nm+".nc", "r")
    slp = a->msl(d,:,:)
    slp = slp*0.01
    plot_temp = gsn_csm_contour_map(wks,slp,mres)
    ;lslp = slp({mres@mpMinLatF:mres@mpMaxLatF},{mres@mpMinLonF:mres@mpMaxLonF})
    lslp = slp({-10:20},{150:180})
    locmin = local_min(lslp, False, 0.0)
    if any(ismissing(locmin@xi)) then
      plot(j) = plot_temp
      j = j + 1
      continue
    end if
    x = llon(locmin@xi)
    y = llat(locmin@yi)
    print(x)
    print(y)
    dum(i) = gsn_add_polymarker(wks,plot_temp,x,y,polyres)
    delete(locmin)
    delete(x)
    delete(y)
    ;overlay(plot_base,plot_temp(i-1))
    plot(j) = plot_temp
    j = j + 1
    delete(plot_temp)
    ;delete(slp)
  end do
  resP	                   = True
  resP@gsnPanelMainString	   = "i"+date(0)+"-v"+date(d)
  gsn_panel(wks,plot,(/6,5/),resP)
  ;delete(plot)
end