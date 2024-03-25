;************************************
; trans_1.ncl
;************************************
;
; Concepts illustrated:
;   - Calculating and plotting a transect
;   - Using gc_latlon to calculate a great circle path
;   - Using linint2_points to interpolate rectilinear grid values to set of lat/lon points
;   - Attaching polylines to a map plot
;   - Explicitly setting tickmarks and labels on the bottom X axis
;************************************
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
;************************************
begin
  init = "2019100900"
  tl = 959
  expnlist = (/"CLIM","EST","MGD"/)
  do l = 0, 2
  ;l = 0
  expn = expnlist(l)
  if (expn .eq. "CLIM") then
    suffix = ""
  else if (expn .eq. "EST") then 
    suffix = "_est"
  else
    suffix = "_mgdsst"
  end if
  end if
  suffix2 = suffix
  in = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl"+sprinti("%d",tl)+"L100"+suffix+"/fcst_pflxpp_"+init+".nc","r")
  print(getfilevarnames(in))
  do d = 18,18
  wks = gsn_open_wks("pdf","cross_outer_pflx_"+init+"_tl"+sprinti("%d",tl)+"_"+expn+"_"+sprinti("%d",d))
  rflx = in->RFLX(d,{1000:70},::-1,:)
  sflx = in->SFLX(d,{1000:70},::-1,:)
  printVarSummary(rflx)
  printMinMax(rflx,0)
  printVarSummary(sflx)
  printMinMax(sflx,0)
  pflx = rflx
  pflx = rflx + sflx
  pflx = pflx*3.6e3
  pflx@long_name = "Precipitation flux"
  pflx@units = "kg/m~S~2~N~/h"
  printVarSummary(pflx)
  printMinMax(pflx,0)
  
  date = cd_string(in->time(d), "%H00 UTC, %D") ;"%Y/%N/%D %H:00")
;************************************
; calculate great circle along transect
;************************************
  leftlat  =  25.1710
  rightlat =  36.5310

  leftlon  =  131.195 + 1.
  rightlon =  137.929 + 1.

  npts     =   100                    ; number of points in resulting transect

  dist     = gc_latlon(leftlat,leftlon,rightlat,rightlon,npts,2)
  points   = ispan(0,npts-1,1)*1.0 
;********************************
; interpolate data to great circle
;********************************
  trans   = linint2_points(pflx&lon,pflx&lat,pflx,True,dist@gclon,dist@gclat,2)
  copy_VarAtts(pflx,trans)          ; copy attributes

  trans!0      = "level"           ; create named dimension and assign
  trans&level    = pflx&level           ; coordinate variable for 0th dimension only
;********************************
; create plot
;********************************

  res                     = True          ; plot mods desired
; label values
;  res@tmXBMode            = "Explicit"    ; explicitly label x-axis
;  res@tmXBValues          = (/points(0),points(npts-1)/) ; points to label
;  res@tmXBLabels          = (/leftlat +", "+leftlon,rightlat+", "+rightlon/)
  nLabels                  = 10               ; arbitrary  
  res@tmXBLabels           = new(nLabels,"string")
  res@tmXBMode             = "Explicit"

;**************************************************
; Draw Full cross section at a specified blue grid line
;**************************************************
  res@tmXBValues         := toint( fspan(0,npts-1,nLabels) )
  do i=0,nLabels-1
    x = dist@gclon(res@tmXBValues(i))
    y = dist@gclat(res@tmXBValues(i))
    res@tmXBLabels(i) = sprintf("%5.1f", y)+"~C~"+sprintf("%5.1f", x)
  end do
  res@tmYRTickSpacingF    = 2
  res@tmYLOn = False
  res@tmYLLabelsOn = False

  res@cnFillOn            = True         ; turn on color
  res@cnLinesOn           = False        ; turn off countour lines
  res@lbOrientation       = "horizontal"   ; vertical label bar
  res@lbBoxLinesOn        = False
  res@lbLabelBarOn = False
  res@pmLabelBarOrthogonalPosF = -0.005        ; move label bar closer to plot
  res@pmLabelBarParallelPosF = 0.45
  res@pmLabelBarHeightF = 0.05

  res@gsnLeftString = ""
  res@gsnRightString = ""

  res@tiMainString        = "";date + " " + expn ; add title
  res@tiXAxisString       = "";"lat/lon along transect"
  res@trYReverse          = True         ; reverse y axis
  ;res@trXReverse          = True         ; reverse x axis (neg longitudes)
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels = (/0.1,0.5,1.0,2.0,3.0,5.0,10.0,15.0,20.0,25.0/)
  res@cnInfoLabelOn = False
  res@cnLineLabelsOn = False
  cmap_default = read_colormap_file("BkBlAqGrYeOrReViWh200")
  cmap = new((/11,4/), float)
  cmap(0,:) = (/0.0,0.0,0.0,0.0/)
  cmap(1,:) = cmap_default(18,:)
  cmap(2,:) = cmap_default(34,:)
  cmap(3,:) = cmap_default(48,:)
  cmap(4,:) = cmap_default(64,:)
  cmap(5,:) = cmap_default(96,:)
  cmap(6,:) = cmap_default(100,:)
  cmap(7,:) = cmap_default(116,:)
  cmap(8,:) = cmap_default(132,:)
  cmap(9,:) = cmap_default(148,:)
  cmap(10,:) = cmap_default(164,:)
  res@cnFillPalette       = cmap  ; set color map
  
  ;plot = gsn_csm_contour(wks,trans,res)  ; create plot
  plot = gsn_csm_pres_hgt(wks,trans,res)  ; create plot
  ;********************************
; show transect on a map
;********************************
  wks = gsn_open_wks("pdf","cross_outer_map_"+init+"_tl"+sprinti("%d",tl)+"_"+expn+"_"+sprinti("%d",d))
  mres              = True                ; plot mods desired
  mres@gsnFrame     = False               ; don't turn page yet
  mres@gsnDraw      = False               ; don't draw yet
  mres@mpDataBaseVersion = "HighRes"
  mres@mpMaxLatF    = 39.0
  mres@mpMinLatF    = 24.0
  mres@mpMaxLonF    = 141.0
  mres@mpMinLonF    = 126.0
  mres@tiMainString        = date + " " + expn ; add title
  mres@tiMainFontHeightF = 0.04
  mres@mpGeophysicalLineColor = "white"
  mres@mpGeophysicalLineThicknessF = 2.0
  ;mres@tiMainString = "Transect Location + prcp" ; title

  mres@cnFillOn = True
  mres@cnLinesOn = False
  mres@cnLineLabelsOn = False
  mres@cnFillPalette             = cmap(1:,:);_default(16:160,:)
  mres@cnLevelSelectionMode      = "ExplicitLevels"  ; manually specify contour levels
  mres@cnLevels = (/.5, 1., 2., 3., 5., 10., 15., 20., 25./)
  mres@lbTitleString = "mm/hr"
  mres@lbTitleFontHeightF = .02
  mres@lbTitlePosition            = "Bottom"            ; title position
  mres@lbTitleDirection           = "Across"            ; title direction
  mres@pmLabelBarOrthogonalPosF = 0.15        ; move label bar closer to plot
  mres@pmLabelBarHeightF = 0.05
  mres@gsnAddCyclic = True
  ;mres@lbLabelBarOn = False
  
  prcp2 = in->RAIN(d, :, :)
  prcp = prcp2
  if (d.gt.0) then
    prcp1 = in->RAIN(d-1, :, :)
    prcp = prcp - prcp1
  end if
  prcp = prcp / 3. ; per 3 hours -> per 1 hour
  printMinMax(prcp, 0)
  ;map = gsn_csm_map(wks,mres)         ; create map
  map = gsn_csm_contour_map(wks,prcp,mres)

; add polyline to map
  pres                  = True           ; polyline mods desired
  pres@gsLineColor      = "gray"          ; color of lines
  pres@gsLineThicknessF = 5.0            ; line thickness   
  ;id = gsn_add_polyline(wks,map,(/leftlon,rightlon/),(/leftlat,rightlat/),pres)
  pmres = True
  pmres@gsMarkerIndex = 16
  pmres@gsMarkerColor = "gray"
  pmres@gsMarkerSizeF = 10
  ldum = new(nLabels,graphic)
  mdum = new(nLabels,graphic)
  do i=0,nLabels-1
    x = dist@gclon(res@tmXBValues(i))
    y = dist@gclat(res@tmXBValues(i))
    if (i .lt. nLabels-1) then
      x2 = dist@gclon(res@tmXBValues(i+1))
      y2 = dist@gclat(res@tmXBValues(i+1))
      ldum(i) = gsn_add_polyline(wks,map,(/x,x2/),(/y,y2/),pres)
    end if
    mdum(i) = gsn_add_polymarker(wks, map, x, y, pmres)
  end do
  draw(map)    ; draws map and polyline
  frame(wks)
  end do
  end do
end
