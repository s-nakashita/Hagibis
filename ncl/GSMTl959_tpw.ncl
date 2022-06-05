load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
begin 
  ;Reading netcdf data
  yyyymmddhh = 2019100912
  sst = "_mgdsst"
  ;fhp1    = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100"+sst+"/fcst_surf_asia_"+yyyymmddhh+".nc", "r")
  fh     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100"+sst+"/fcst_phy2m_snp_"+yyyymmddhh+".nc", "r")

  t6h  = fh->time
  date = cd_string(t6h, "%N%D %H:%M")
  ntim = dimsizes(date)

  ;Reading track data
  data = readAsciiTable("../pytrack/track"+yyyymmddhh+"_gsm_tl959"+sst+".txt", 7, "float", 0)
  npoints = dimsizes(data(:,0))
  dd = toint(data(:,2))
  hh = toint(data(:,3))
  lon1 = data(:,4)
  lat1 = data(:,5)
  delete(data)

  pltDir = "./"                ; where output will be written
  pltName= "GSMtl959_tpw"+sst+"_"+yyyymmddhh                                   
  pltType= "png"               ; ps, eps, png, pdf, x11 .....

;****************************************************************************
; Graphic resources
;****************************************************************************
  do d = 8, ntim-1
  lh = fh->TPW(d, :, :)
  i = 0
  pltPath = pltDir+pltName+"_"+cd_string(t6h(d), "%D%H")
  wks     = gsn_open_wks(pltType, pltPath) 
  
  res1                            = True     ; plot mods desired
  res1@gsnDraw                    = False    ; will panel plots later
  res1@gsnFrame                   = False

  res1@cnFillOn                   = True     ; turn on color fill
  res1@cnLinesOn                  = False    ; turn off contour lines
  res1@cnLineLabelsOn             = False    ; turn off contour line labels

;**********************************************************************
; We are turning the labelbar off for all but the bottom row of plots
;**********************************************************************
  ;res1@lbLabelBarOn               = False
  res1@lbTitleFontHeightF         = .01                 ; make title smaller
  res1@lbTitleString              = "kg/m^2"            ; title string
  res1@lbTitlePosition            = "Bottom"            ; title position
  res1@lbTitleDirection           = "Across"            ; title direction
  ;res1@lbLabelStride              = 2
  res1@lbLabelFontHeightF         = 0.01
  ;res1@pmLabelBarHeightF          = 0.1
  ;res1@lbBoxEndCapStyle           = "TriangleBothEnds"
  ; Set resources for map
  res1@mpMinLatF                  = 15.
  res1@mpMaxLatF                  = 45.
  res1@mpMinLonF                  = 120.
  res1@mpMaxLonF                  = 150.

  res1@mpDataBaseVersion          = "HighRes"
  res1@mpOutlineOn                = True
  res1@mpFillOn                   = False                ; turn off gray fill
  ;res@mpOutlineBoundarySets      = "National"

  ;res1@gsnAddCyclic               = False            ; regional data, don't add cyclic point
 
  ; Set resources for tickmarks
  res1@gsnStringFont              = "helvetica-bold"  ; will be used for gsnCenterString later
  res1@gsnCenterStringFontHeightF = 0.04              ; ditto
  res1@gsnLeftString              = ""                ; remove the left and right substrings
  res1@gsnRightString             = date(d)

  ; Set resources specifically for plots in column 1
  res1@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels

  ; For track
  pmres = True
  pmres@gsMarkerIndex = 16
  pmres@gsMarkerSizeF = 5.0
  pmres@gsMarkerColor = "red"
;****************************************************************************
; Create the plots for the panel
;****************************************************************************
  ; first column : raw values
  ind_t = 0
  day = toint(cd_string(t6h(d), "%D"))
  hour = toint(cd_string(t6h(d), "%H"))
  do ip = 0, npoints-1
    if(dd(ip).eq.day) .and. (hh(ip).eq.hour) then
    break
    end if
    ind_t = ind_t + 1
  end do
  print(ind_t)
  res1@cnMinLevelValF            = 7.5           ; min level
  res1@cnMaxLevelValF            = 72.5          ; max level
  res1@cnLevelSpacingF           = 2.5           ; contour interval
  cmap = read_colormap_file("NCV_bright")
  printVarSummary(cmap)
  newcmap = new((/207,4/), float)
  newcmap(:94,:) = cmap(17:111,:)
  newcmap(95:,:) = cmap(144:,:)
  res1@cnFillPalette             = newcmap
  
  plot = gsn_csm_contour_map_ce(wks, lh, res1)
  ;dum  = gsn_add_polymarker(wks, plot, lon1(ind_t), lat1(ind_t), pmres)
  draw(plot)
  frame(wks)
  delete(plot)
  delete(wks)
  end do
end
