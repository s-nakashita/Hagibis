load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
begin
  ncols = 7
  init = "2019100912"
  ;wks = gsn_open_wks("png", "track+prmsl_tl479_"+init)
  wks = gsn_open_wks("png", "track_tl959_0900_sst")
  wks = gsn_open_wks("png", "track_tl479_0912")
  ;wks = gsn_open_wks("png", "track+prmsl_merra2")

  in  = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_mgdsst/fcst_surf_"+init+".nc","r")
  ;in  = addfile("/Volumes/dandelion/netcdf/era5/2019/10/12/slev_12.nc","r")
  ;in  = addfile("/Volumes/dandelion/netcdf/nasa/merra2/MERRA2_400.tavg3_3d_asm_Nv.20191012.nc4","r")
  time = in ->time 
  ntim = dimsizes(time)
  time_land = cd_inv_calendar( 2019, 10, 12, 12, 00, 00, time@units, 0)
  print(time_land)
  pres = in ->PSEA({time_land},:,:)
  ;pres = in ->msl({time_land},:,:)
  ;time_pre = cd_inv_calendar( 2019, 10, 12, 12, 00, 00, time@units, 0)
  ;time_pos = cd_inv_calendar( 2019, 10, 12, 15, 00, 00, time@units, 0)
  ;pres_pre = in ->SLP({time_pre},:,:)
  ;pres_pos = in ->SLP({time_pos},:,:)
  ;pres = pres_pre
  ;pres = 0.5 * (pres_pre + pres_pos)

  res = True 
  res@gsnDraw = False 
  res@gsnFrame = False 
  ;res@tiMainString = "Best Track of Hagibis (2019)"

  res@mpDataBaseVersion = "RANGS"
  res@mpOutlineOn = True
  res@mpFillOn    = True
  res@mpLandFillColor = "tan"
  res@mpOceanFillColor = "LightCyan"
  res@mpInlandWaterFillColor = res@mpOceanFillColor
  ;res@mpMinLatF = 28.0
  res@mpMinLatF = 18.0
  res@mpMaxLatF = 40.0
  res@mpMinLonF = 133.0
  res@mpMaxLonF = 145.0
  res@gsnRightString = " "
  res@gsnLeftString = " "
  ;res@cnInfoLabelOn = False
  ;res@cnLineThicknessF = 2.0
  ;res@cnLevelSpacingF = 4.0

  ;pres = pres * 0.01
  ;plot = gsn_csm_contour_map_ce(wks, pres, res)
  plot = gsn_csm_map_ce(wks, res)
  data = readAsciiTable("bst_hagibis.txt", ncols, "float", 0)
  ;print(data(:,4:6))
  npoints = dimsizes(data(:,0))
  delete(data)
  plres = True 
  pmres = True 
  
  ;files = (/"track"+init+"_gsm.txt", \
  ;"rjtd/track"+init+".txt", \
  ;"track"+init+"_gsm_tl959.txt", \
  ;"track"+init+"_gsm_tl959_est.txt", \
  ;"track"+init+"_gsm_tl959_mgdsst.txt"/);, \
  ;"track"+init+"_gsm_tl479_est.txt",\
  ;"track"+init+"_gsm_tl479_est+p.txt",\
  ;"track"+init+"_gsm_tl479_est+p2.txt"/)
  ;files = (/;"track2019100900_gsm_tl479_est.txt",\
  ;"track"+init+"_gsm_tl479_est.txt",\
  ;"track"+init+"_gsm_tl479_est+p.txt"/)
  files = (/;"track_era5.txt",\
  "track"+init+"_gsm_tl479_est.txt",\
  "track"+init+"_gsm_tl479_est+p.txt",\
  "track"+init+"_gsm_tl479_est+p2.txt",\
  ;"track"+init+"_gsm_tl479_est+pf.txt",\
  "track"+init+"_gsm_tl479_est+pn.txt"\
  /)
  ;files = (/"track2019100900_gsm.txt", \
  ;"track2019100912_gsm.txt", \
  ;"rjtd/track2019100900.txt", \
  ;"rjtd/track2019100912.txt", \
  ;"track2019100900_gsm_tl959_est.txt",\
  ;"track2019100912_gsm_tl959_est.txt"/)
  ;"track"+init+"_gsm_sstclim.txt", \
  ;"track"+init+"_gsm_tl319l100_sstclim.txt", \
  ;"track"+init+"_gsm_tl319l100_sstclim+mysub.txt"/)
  ;"track"+init+"_gsm_sstclim.txt", \
  ;"track"+init+"_gsm_sstjra55.txt", \
  ;"track"+init+"_gsm_mgdsst.txt"/)
  ;files = (/"track_merra2.txt"/)
  ;files = (/"track2019100900_gsm_tl959.txt", \
  ;"track2019100900_gsm_tl959_est.txt", \
  ;"track2019100900_gsm_tl959_mgdsst.txt"/);, \
  ;;"track2019100912_gsm_tl959.txt", \
  ;;"track2019100912_gsm_tl959_est.txt", \
  ;;"track2019100912_gsm_tl959_mgdsst.txt"/)
  dum3 = new((/npoints, dimsizes(files)/), graphic)
  dum4 = new((/npoints, dimsizes(files)/), graphic)
  plres@gsLineThicknessF = 4.0
  ;pmres@gsMarkerIndex = 4
  ;pmres@gsMarkerThicknessF = 2.0
  df = new((/npoints, dimsizes(files)/), integer)
  hf = new((/npoints, dimsizes(files)/), integer)
  lonf = new((/npoints, dimsizes(files)/), float)
  latf = new((/npoints, dimsizes(files)/), float)
  ;colors   = (/"red","blue","red","blue","red","blue"/)
  ;colors   = (/"magenta","cyan","red","blue","green"/)
  ;colors   = (/"blue","red"/)
  colors   = (/"red","blue","green","blue"/)
  ;colors   = (/"magenta","magenta","cyan","cyan","blue","blue"/)
  ;colors   = (/"magenta","cyan","blue","green"/)
  ;colors   = (/"darkorange"/)
  ;colors   = (/"red","blue","green","red","blue","green"/)
  ;patterns = (/    2,     2,    0,     0,    3,     3/)
  ;patterns = (/    2,    2,     0,    0,    0/)
  ;patterns = (/    2,    0,     0/)
  patterns = (/    0,    0,    0,    2,    0/)
  ;patterns = (/    2,    0,     2,    0,    2,    0/)
  ;patterns = (/     2,    2,    2,    0,    0,    0/)
  ;markers  = (/     4,    4,    4,   16,   16,   16/)
  ;patterns = (/     0,    0,    0/)
  markers  = (/    16,   16,   16,   16,   16,   16/)
  marker_seq = (/16,9,7,6/)
  do j = 0, dimsizes(files)-1
    tfile = files(j)
    data = readAsciiTable(tfile, ncols, "float", 0)
    print(dimsizes(data))
    mypoints = dimsizes(data(:,0))
    df(0:mypoints-1,j) = toint(data(:,2))
    hf(0:mypoints-1,j) = toint(data(:,3))
    lonf(0:mypoints-1,j) = data(:,4)
    latf(0:mypoints-1,j) = data(:,5)
    plres@gsLineColor = colors(j)
    plres@gsLineDashPattern = patterns(j)
    pmres@gsMarkerColor = colors(j)
    ;pmres@gsMarkerIndex = markers(j)
    k = 0
    do i = 0, mypoints-2
      dum3(i,j) = gsn_add_polyline(wks, plot, (/lonf(i,j),lonf(i+1,j)/),(/latf(i,j),latf(i+1,j)/),plres)
      if (df(i,j) .eq. 12) .and. (hf(i,j) .eq. 12) then
        ;pmres@gsMarkerIndex = 4
        pmres@gsMarkerSizeF = 16.0
        pmres@gsMarkerThicknessF = 8.0
      else 
        ;pmres@gsMarkerIndex = 16
        ;pmres@gsMarkerSizeF = 10.0
        pmres@gsMarkerSizeF = 8.0
        pmres@gsMarkerThicknessF = 5.0
      end if
      pmres@gsMarkerIndex = marker_seq(k)
      k = k + 1
      if (k .gt. 3) then
        k = 0
      end if
      dum4(i,j) = gsn_add_polymarker(wks, plot, lonf(i,j), latf(i,j), pmres)
    end do
    if (df(i,j) .eq. 12) .and. (hf(i,j) .eq. 12) then
      ;pmres@gsMarkerIndex = 4
      pmres@gsMarkerSizeF = 16.0
      pmres@gsMarkerThicknessF = 5.0
    else 
      ;pmres@gsMarkerIndex = 16
      pmres@gsMarkerSizeF = 10.0
      pmres@gsMarkerThicknessF = 3.0
    end if
    pmres@gsMarkerIndex = marker_seq(k)
    dum4(i,j) = gsn_add_polymarker(wks, plot, lonf(i,j), latf(i,j), pmres)
    delete(data)
  end do
;; best track
  data = readAsciiTable("bst_hagibis.txt", ncols, "float", 0)
  ;print(data(:,4:6))
  npoints = dimsizes(data(:,0))
  dd = toint(data(:,2))
  hh = toint(data(:,3))
  lon = data(:,4)
  lat = data(:,5)

  dum = new(npoints, graphic)
  dum2 = new(npoints, graphic)
  plres@gsLineThicknessF = 6.0
  pmres@gsMarkerIndex = 4
  pmres@gsMarkerThicknessF = 3.0
  pmres@gsMarkerColor = "black"
  do i = 0, npoints-2
    ;dum(i) = gsn_add_polyline(wks, plot, (/lon(i),lon(i+1)/),(/lat(i),lat(i+1)/),plres)
    if (dd(i) .eq. 12) .and. (hh(i) .eq. 12) then
      pmres@gsMarkerIndex = 4
      pmres@gsMarkerSizeF = 16.0
      pmres@gsMarkerThicknessF = 5.0
    else 
      pmres@gsMarkerIndex = 4
      pmres@gsMarkerSizeF = 8.0
      pmres@gsMarkerThicknessF = 3.0
    end if
    dum2(i) = gsn_add_polymarker(wks, plot, lon(i), lat(i), pmres)
  end do
  dum2(npoints-1) = gsn_add_polymarker(wks, plot, lon(npoints-1), lat(npoints-1), pmres)
  delete(data)

;  ;; legends
  xleg = (/0.71,0.71,0.71,0.71/)   ; Location of
  xtxt = (/0.76,0.76,0.76,0.76/)   ; legend markers
  yleg = (/0.29,0.25,0.21,0.17/)   ; and text
  ytxt = (/0.29,0.25,0.21,0.17/)   ; strings.

  num_distinct_markers = dimsizes(marker_seq)
  labels = (/"00/12Z","03/15Z","06/18Z","09/21Z"/)
  txres = True 
  txres@txFontHeightF = 0.015
  pmres@gsMarkerColor  = "black"
  pmres@gsMarkerSizeF = 10.0
  pmres@gsMarkerThicknessF = 3.0
  do i = 0, num_distinct_markers-1
    pmres@gsMarkerIndex = marker_seq(i)
;
; Add marker and text for the legend.
;
    gsn_polymarker_ndc(wks,          xleg(i),yleg(i),pmres)
    gsn_text_ndc      (wks,labels(i),xtxt(i),ytxt(i),txres)
  end do
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