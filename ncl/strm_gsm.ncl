;----------------------------------------------------------------------
; vector_1.ncl
;
; Concepts illustrated:
;   - Overlaying vectors and filled contours on a map
;   - Manually attaching lat/lon coordinate arrays to a variable
;   - Changing the length of the smallest vector as a fraction of the reference vector
;   - Moving the vector reference annotation to the top right of the plot
;   - Drawing curly vectors
;   - Setting the color for vectors
;   - Making the labelbar be vertical
;   - Increasing the thickness of vectors
;----------------------------------------------------------------------
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

begin

;---Open netCDF file
  ;"cntl","p","pn","p2","pf"
  ;prtb = "pn"
  dev = "ps"
  d_offset = 0
  t_offset = 1
  if (prtb .eq. "cntl") then
  a = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est/fcst_p_2019100912.nc","r")
  else
  a = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl479L100_est_0912+"+prtb+"/fcst_p_2019100912.nc","r")
  end if
  um   = a->U(:,:,:,:)
  vm   = a->V(:,:,:,:)
  lat  = a->lat
  lon  = a->lon
  plev = a->level
  time = a->time
  date = cd_string(time, "%Y%N%D%H")
  lon0 = 120.
  lon1 = 180.
  lat0 = 0.
  lat1 = 40.
  
  ;out = "strm_gsmtl479_"+prtb
  out = prtb
  figdir = "strm_gsmtl479"
  system("mkdir -p "+figdir)
  cmap = read_colormap_file("hotcold_18lev")
  nrow=1
  ncol=4
  plotnum = ispan(0, nrow*ncol-1, 1)
  ;plotnum = (/0,3,6,9,1,4,7,10,2,5,8,11/)
  ;plevs = (/850,500,300/)
  plevs = (/100,300,500,850/)
  do k = 0, dimsizes(plevs)-1
    u = um(:,{plevs(k)},::-1,:)
    v = vm(:,{plevs(k)},::-1,:)
    hgt = plevs(k)
    sfvp = uv2sfvpF(u,v)
    sf = sfvp(0,:,:,:)
    copy_VarCoords(u,sf)
    sf@long_name = "stream function"
    sf@units     = "m^2/s"
    sf_zm = zonalAve(sf)
    ;sf = sf - conform(sf,sf_zm,(/0,1/))
    printVarSummary(sf)
    printMinMax(sf,0)
    sf = sf * 1.0e-6

    plot = new(dimsizes(plotnum),graphic)
    wks = gsn_open_wks(dev,figdir+"/"+out+"_"+hgt)
;---Create plot
    res                = True               ; plot mods desired
    res@gsnDraw	      = False
    res@gsnFrame	      = False
;     res@cnInfoLabelOn	      = False
     
    res@gsnCenterStringOrthogonalPosF = 0.02
    res@gsnLeftString	    = "10~S~6~N~"+sf@units
    res@gsnLeftStringFontHeightF = 0.02
    res@gsnRightString	    = ""
    res@gsnRightStringFontHeightF = 0.02
    res@gsnRightStringParallelPosF = 1.1    
    res@gsnAddCyclic	       = True
    res@cnFillOn             = True               ; turn on color for contours
    res@cnLinesOn            = True               ; turn off contour lines
    res@cnLineLabelsOn       = False              ; turn off contour line labels
    res@cnFillPalette        = cmap

;    res@mpLandFillColor      = "gray"            ; set land to be gray
    res@mpFillOn	      = False
    res@mpMinLonF            =  lon0               ; select a subregion
    res@mpMaxLonF            =  lon1
    res@mpCenterLonF         = (lon0+lon1)/2
    res@mpMinLatF            =  lat0
    res@mpMaxLatF            =  lat1

    res@lbLabelBarOn = False
    ;res@lbOrientation            = "Vertical"     ; vertical label bar
    ;res@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer
    ;res@vpXF = 0.15
    ;res@vpWidthF = 0.7

; note, when doing a subregion, NCL determines the range of the data from
; the full domain. If you wish to just consider the domain you are plotting,
; you must manually set those levels.

    res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
    res@cnMinLevelValF       = -24.               ; set min contour level
    res@cnMaxLevelValF       = 24.                ; set max contour level
    res@cnLevelSpacingF      = 4.                 ; set contour spacing
    res@gsnContourNegLineDashPattern = 2
    res@gsnContourZeroLineThicknessF = 3.

    i = d_offset
    d_end = min((/dimsizes(date),t_offset+dimsizes(plotnum)/))
    do d=t_offset,d_end-1
      res@gsnCenterString = date(d)
      plot(plotnum(i)) = gsn_csm_contour_map_ce(wks,sf(d,:,:),res)
      i = i + 1
      if (i.ge.dimsizes(plotnum)) then
        break
      end if
    end do
    delete(sf)
    resp = True
    resp@gsnPanelLabelBar = True
    gsn_panel(wks,plot,(/nrow,ncol/),resp)
    delete(plot)
    delete(wks)
    ;if(dev.eq."x11") then
    ;exit()
    ;end if
   end do
end
