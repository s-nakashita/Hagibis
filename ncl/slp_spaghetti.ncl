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
  res = True
  ;res@mpDataBaseVersion="MediumRes"
  ;res@mpProjection = "LambertEqualArea"
  ;res@mpLimitMode = "LatLon"
  res@mpMinLatF = -10.                                       ; set the minimum latitude
  res@mpMaxLatF = 40.                                      ; set the maximum latitude
  res@mpMinLonF = 110.                                     ; set the minimum longitude
  res@mpMaxLonF = 180.                                     ; set the maximum longitude
  res@mpCenterLonF = (res@mpMaxLonF + res@mpMinLonF)/2.    ; set the longitude at the center of the plot
  ;res@mpOutlineOn = True                                   ; turn the map outlines on
  ;res@mpOutlineBoundarySets = "National"                   ; draw national boundaries
  ;res@mpOutlineSpecifiers = (/"Canada : Provinces","United States : States"/)  ; draw US States, Canadian provinces
  res@mpLandFillColor = "gray90"			            ; shade land with gray90
  res@gsnFrame = False                                     ; don't advance the frame
  res@gsnDraw  = False                                     ; don't draw the plot
  res@gsnMaximize = True                                   ; maximize the plot
  ;res@gsnPaperOrientation = "portrait"                     ; when maximizing, keep the orientation as portrait
  ;res@gsnLeftString = egrr@initial_time
  ;res@gsnRightString = egrr@units+" *10~S~-1~N~"
  ;res@tiMainString = "500hPa Temperature"	
  ;res@tiMainFontHeightF = 0.012                 ; set the main title font height
;  res@gsnStringFontHeightF = 0.015
  am = addfile("/Volumes/dandelion/netcdf/tigge/"+yyyy+"/"+orig+"/glb_"+yyyymmddhh+"_mean.nc", "r")
  t = am->time
  date = cd_string(t,"%N%D%H")
  d = ind(date.eq."101012")
  res@gsnCenterString = date(d)
  print(date)
  delete(am)

  wks = gsn_open_wks("png","slp_spaghetti_"+orig+"_i"+mmddhh+"_"+date(d))                        ; send graphics to PNG file

  plot_base = gsn_csm_map(wks,res)               ; create the blank base plot
	
  mres = True                                    ; create a resource list for contour plots
  mres@cnLevelSelectionMode = "ExplicitLevels"       ; set manual contour levels
  ;mres@cnLevels = (/1000.,1002.,1004.,1006.,1008.,1010.,1011.,1012.,1013.,1014.,1015.,1016.,1017.,1018.,1019.,1020./)
  mres@cnLevels = (/1008.,1010.,1012./)
  ;mres@cnLevelSelectionMode = "ManualLevels"     ; manually select the levels
  ;mres@cnMinLevelValF  = 1000.                    ; set the minimum contour level
  ;mres@cnMaxLevelValF  = 310.                    ; set the maximum contour level
  ;mres@cnLevelSpacingF = 2.                      ; set the contour interval
  mres@cnMonoLineColor = False
  mres@cnLineColors = (/"blue","purple","red"/)   ; contour the lines blue
  mres@gsnFrame = False                          ; don't advance the frame
  mres@gsnDraw  = False                          ; don't draw the plot
  mres@cnInfoLabelOn         = False             ; don't draw the info label
  ;mres@cnLineLabelPerimOn = False                ; turn the line label perimeter off
;  mres@cnLineLabelBackgroundColor = -1
  mres@cnLineLabelsOn = False
;  mres@cnLineLabelInterval = 5                   ; label every other line
;  mres@cnMonoLineLabelFontColor = False
;  mres@cnLineLabelFontColors = (/"blue","purple","red"/) ; color the line labels blue
;  mres@cnLineLabelFontHeightF = 0.005             ; set the line label font height
  ;mres@tiMainFontHeightF = 0.012                 ; set the main title font height
	
  mres@gsnLeftString = "hPa"
  mres@gsnCenterString = ""
  mres@gsnRightString = ""
	
  plot_temp = new(26,"graphic")         
  ;mres@cnLineColor = "blue"                      ; contour the lines blue
  do i = 1, 26
    if (i .gt. 1) then
      mres@gsnLeftString = ""
    end if
    if (i .lt. 10) then
      nm = "0"+sprinti("%d",i)
    else 
      nm = sprinti("%2d",i)
    end if
    print(nm)
    a = addfile("/Volumes/dandelion/netcdf/tigge/"+yyyy+"/"+orig+"/glb_"+yyyymmddhh+"_"+nm+".nc", "r")
    slp = a->msl(d,:,:)
    slp = slp*0.01
    plot_temp(i-1) = gsn_csm_contour(wks,slp,mres)
    overlay(plot_base,plot_temp(i-1))
    ;delete(slp)
  end do

  draw(plot_base)
  frame(wks)
end