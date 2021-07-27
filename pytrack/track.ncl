begin
  ncols = 7
  wks = gsn_open_wks("pdf", "track+prmsl")

  in  = addfile("fcst_p.nc","r")
  time = in ->time 
  ntim = dimsizes(time)
  pres = in ->prmslmsl(ntim-1,:,:)

  res = True 
  res@gsnDraw = False 
  res@gsnFrame = False 
  ;res@tiMainString = "Best Track of Hagibis (2019)"

  res@mpDataBaseVersion = "RANGS"
  res@mpOutlineOn = True
  res@mpMinLatF = 28.0
  res@mpMaxLatF = 38.0
  res@mpMinLonF = 135.0
  res@mpMaxLonF = 145.0
  res@gsnRightString = " "
  res@gsnLeftString = " "
  res@cnInfoLabelOn = False
  res@cnLineThicknessF = 2.0
  res@cnLevelSpacingF = 4.0

  plot = gsn_csm_contour_map_ce(wks, pres, res)

  data = readAsciiTable("bst_hagibis.txt", ncols, "float", 0)
  ;print(data(:,4:6))
  npoints = dimsizes(data(:,0))
  dd = toint(data(:,2))
  hh = toint(data(:,3))
  lon = data(:,4)
  lat = data(:,5)

  dum = new(npoints, graphic)
  dum2 = new(npoints, graphic)
  plres = True 
  plres@gsLineThicknessF = 6.0
  pmres = True 
  pmres@gsMarkerIndex = 4
  pmres@gsMarkerThicknessF = 3.0
  do i = 0, npoints-2
    dum(i) = gsn_add_polyline(wks, plot, (/lon(i),lon(i+1)/),(/lat(i),lat(i+1)/),plres)
    if (dd(i) .eq. 12) .and. (hh(i) .eq. 12) then
      pmres@gsMarkerIndex = 4
      pmres@gsMarkerSizeF = 16.0
      pmres@gsMarkerThicknessF = 5.0
    else 
      pmres@gsMarkerIndex = 16
      pmres@gsMarkerSizeF = 10.0
      pmres@gsMarkerThicknessF = 3.0
    end if
    dum2(i) = gsn_add_polymarker(wks, plot, lon(i), lat(i), pmres)
  end do
  dum2(npoints-1) = gsn_add_polymarker(wks, plot, lon(npoints-1), lat(npoints-1), pmres)
  delete(data)

  files = (/"track2019100900_gsm.txt", "track2019100912_gsm.txt", "track2019100900_gsm_sstclim.txt", "track2019100912_gsm_sstclim.txt"/)
  dum3 = new((/npoints, 4/), graphic)
  dum4 = new((/npoints, 4/), graphic)
  plres@gsLineThicknessF = 4.0
  ;pmres@gsMarkerIndex = 4
  ;pmres@gsMarkerThicknessF = 2.0
  df = new((/npoints, 4/), integer)
  hf = new((/npoints, 4/), integer)
  lonf = new((/npoints, 4/), float)
  latf = new((/npoints, 4/), float)
  colors   = (/"red","blue","red","blue"/)
  patterns = (/    2,     2,    0,     0/)
  do j = 0, 3
    tfile = files(j)
    data = readAsciiTable(tfile, ncols, "float", 0)
    mypoints = dimsizes(data(:,0))
    df(0:mypoints-1,j) = toint(data(:,2))
    hf(0:mypoints-1,j) = toint(data(:,3))
    lonf(0:mypoints-1,j) = data(:,4)
    latf(0:mypoints-1,j) = data(:,5)
    plres@gsLineColor = colors(j)
    plres@gsLineDashPattern = patterns(j)
    pmres@gsMarkerColor = colors(j)
    do i = 0, mypoints-2
      dum3(i,j) = gsn_add_polyline(wks, plot, (/lonf(i,j),lonf(i+1,j)/),(/latf(i,j),latf(i+1,j)/),plres)
      if (df(i,j) .eq. 12) .and. (hf(i,j) .eq. 12) then
        pmres@gsMarkerIndex = 4
        pmres@gsMarkerSizeF = 16.0
        pmres@gsMarkerThicknessF = 5.0
      else 
        pmres@gsMarkerIndex = 16
        pmres@gsMarkerSizeF = 10.0
        pmres@gsMarkerThicknessF = 3.0
      end if
      dum4(i,j) = gsn_add_polymarker(wks, plot, lonf(i,j), latf(i,j), pmres)
    end do
    if (df(i,j) .eq. 12) .and. (hf(i,j) .eq. 12) then
      pmres@gsMarkerIndex = 4
      pmres@gsMarkerSizeF = 16.0
      pmres@gsMarkerThicknessF = 5.0
    else 
      pmres@gsMarkerIndex = 16
      pmres@gsMarkerSizeF = 10.0
      pmres@gsMarkerThicknessF = 3.0
    end if
    dum4(i,j) = gsn_add_polymarker(wks, plot, lonf(i,j), latf(i,j), pmres)
    delete(data)
  end do
;  ;; legends
;  gres = True
;  gres@lgPerimOn = True
;  gres@YPosPercent = 85.5    ; expressed as %, 0->100, sets position of top border of legend
;                            ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
;  gres@XPosPercent = 55.    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)
;;
;  lineres = True
;  lineres@lgLineColors = colors ; line colors
;  lineres@lgDashIndexes = patterns
;  lineres@lgLineThicknesses = 4                          ; line thicknesses
;  lineres@LineLengthPercent = 10.                        ; expressed as %, 0->100, length of line
;;
;  textres = True
;  textres@lgLabels = (/"ope. 00Z","ope. 12Z","exp. 00Z","exp. 12Z"/)  ; legend labels (required)
;
;  simple_legend_ndc(wks,gres,lineres,textres)
  draw(plot)
  frame(wks)

  ;; text
  dum5 = new(5, graphic)
  txres = True 
  txres@txPerimOn = False
  ;txres@txBackgroundFillColor = "white"
  txres@txFontThicknessF = 10.0
  txres@txFontHeightF = 0.02
  txres@txPerimThicknessF = 3.0
  i = 0
  dum5(i) = gsn_add_text(wks, plot, "best", 139.0, 28.5, txres)
  i = i + 1
  txres@txFontHeightF = 0.015
  ;delete(txres@txBackgroundFillColor)
  txres@txPerimColor = "blue"
  dum5(i) = gsn_add_text(wks, plot, "exp. 12UTC", 135.5, 33.0, txres)
  i = i + 1
  txres@txPerimColor = "red"
  dum5(i) = gsn_add_text(wks, plot, "exp. 00UTC", 140.5, 33.5, txres)
  i = i + 1
  txres@txPerimColor = "blue"
  txres@txPerimDashPattern = 2
  dum5(i) = gsn_add_text(wks, plot, "ope. 12UTC", 136.5, 36.0, txres)
  i = i + 1
  dum5(i) = gsn_add_text(wks, plot, "ope. 00UTC", 142.5, 36.5, txres)
  ;draw(plot)
  ;frame(wks)
end