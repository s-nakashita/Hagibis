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
; for regrid
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
  
begin
  yyyymmddhh = 2019100912
  MDH = str_get_cols(sprinti("%0.10d",yyyymmddhh),4,9)
;---Open netCDF file
  a = addfile("./inv1d_"+yyyymmddhh+"_mean.nc","r")

; Read in Sea Surface Temperature Anomalies
; Read in U and V at 1000 mb [subscript 0]
; Determine the subscripts corresponding to Jan 88
;
  l = 4
  plev	= a->level
  time  = a->time
  print(time)
  print(time@units)
  ;date  = cd_string(time, "%Y%N%D%H")
  ;print(date)
  
  level = sprintf("%.0f", plev(l))
  wks = gsn_open_wks("png","inv_"+yyyymmddhh)
  plot = new(8,graphic)
  cmap = read_colormap_file("BlAqGrYeOrReVi200")
  
  ; trackfile
  trackfile = "../pytrack/jma/gtrack"+yyyymmddhh+"_mean.txt"
  track = asciiread(trackfile, (/21, 7/), "float") 

  sig = new(3, float)
  sig(0) = 7./22.
  sig(1) = 1./2.
  sig(2) = 2./11.
  ;vdate = time(0)
  do d = 2, 2
    ngrid = toint(a->GridSize(2*d))
    print(ngrid)
    msl1d = a->msl(2*d,:ngrid-1)
    u1d   = a->u(2*d,{850:300},:ngrid-1)
    v1d   = a->v(2*d,{850:300},:ngrid-1)
    lat1d = a->lat1d(2*d,:ngrid-1)
    lon1d = a->lon1d(2*d,:ngrid-1)
    printVarSummary(lat1d)
    printVarSummary(msl1d)
    printVarSummary(u1d)
    printVarSummary(v1d)

    ; vertical integrate
    u1d = u1d * conform_dims(dimsizes(u1d), sig, 0)
    v1d = v1d * conform_dims(dimsizes(v1d), sig, 0)
    up  = dim_sum_n(u1d, 0)
    vp  = dim_sum_n(v1d, 0)

    ; interpolate 2D grid
    lon0 = 110.
    lon1 = 180.
    lat0 = -10.
    lat1 = 40.
    dlon = 1.0
    dlat = 1.0
    nlon = toint((lon1 - lon0)/dlon)
    nlat = toint((lat1 - lat0)/dlat)
    lon = fspan(lon0, lon1, nlon+1)
    lat = fspan(lat0, lat1, nlat+1)
    ;print(lon)
    ;print(lat)

    ;ndim = dimsizes(msl1d)
    ;print(ndim)
    lonin = new(ngrid,float)
    latin = new(ngrid,float)
    lonin = tofloat((/lon1d/))
    latin = tofloat((/lat1d/))
    printMinMax(lonin, 0)
    printMinMax(latin, 0)
    msl = natgrid_Wrap(lonin, latin, msl1d, lon, lat)
    mslt = transpose(msl)
    mslt = mslt * 0.01
    mslt!0 = "lat"
    mslt!1 = "lon"
    mslt&lat = lat
    mslt&lon = lon
    mslt&lat@units = "degrees_north"
    mslt&lon@units = "degrees_east"
    u = natgrid_Wrap(lonin, latin, up, lon, lat)
    ut = transpose(u)
    ut!0 = "lat"
    ut!1 = "lon"
    ut&lat = lat
    ut&lon = lon
    ut&lat@units = "degrees_north"
    ut&lon@units = "degrees_east"
    printVarSummary(ut)
    v = natgrid_Wrap(lonin, latin, vp, lon, lat)
    vt = transpose(v)
    vt!0 = "lat"
    vt!1 = "lon"
    vt&lat = lat
    vt&lon = lon
    vt&lat@units = "degrees_north"
    vt&lon@units = "degrees_east"
    printVarSummary(vt)
    delete(msl1d)
    delete(u1d)
    delete(v1d)
    delete(up)
    delete(vp)
    delete(lon1d)
    delete(lat1d)
    delete(lonin)
    delete(latin)
    ;exit()

     ;print(vdate)
     ;index   = ind(date.eq.vdate)
     ;print(index)
;---Create plot
    res                      = True               ; plot mods desired
    res@gsnDraw	      = False
    res@gsnFrame	      = False
;   res@cnInfoLabelOn	      = False

    res@gsnCenterString      = "ft="+time(2*d)+"h"
    ;res@gsnCenterStringOrthogonalPosF = 0.05
    res@gsnLeftString        = " "
    res@gsnAddCyclic	       = False
    res@cnFillOn             = False               ; turn on color for contours
    ;res@cnFillMode	      = "RasterFill"
    ;res@cnRasterSmoothingOn  = True
    res@cnLinesOn            = True              ; turn off contour lines
    res@cnLineLabelsOn       = True              ; turn off contour line labels
;     res@gsnScalarContour     = True               ; contours desired
;     res@cnFillPalette        = cmap


    res@mpLandFillColor      = "gray"            ; set land to be gray
    res@mpMinLonF            =  lon0              ; select a subregion
    res@mpMaxLonF            =  lon1
    res@mpMinLatF            =  lat0 
    res@mpMaxLatF            =  lat1

    res@lbOrientation            = "Vertical"     ; vertical label bar
    res@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer

; note, when doing a subregion, NCL determines the range of the data from
; the full domain. If you wish to just consider the domain you are plotting,
; you must manually set those levels.
    res@cnLevelSelectionMode = "ExplicitLevels"       ; set manual contour levels
    res@cnLevels = (/1000.,1002.,1004.,1006.,1008.,1010.,1011.,1012.,1013.,1014.,1015.,1016.,1017.,1018.,1019.,1020./)
    ;res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
;     res@cnMaxLevelCount	   = 10
    ;res@cnMinLevelValF       = 100000.0               ; set min contour level
    ;res@cnMaxLevelValF       = 102000.0                 ; set max contour level
    ;res@cnLevelSpacingF      = 200.0               ; set contour spacing

    res@gsnSpreadColorEnd     = -3
    res@gsnScalarContour      = True              ; vectors over scalars
    res@vcRefMagnitudeF       = 10.0              ; vec magnitude
    res@vcRefLengthF          = 0.045             ; size of ref vector
    res@vcLineArrowColor      = "black"           ; white vectors
    res@vcLineArrowThicknessF = 2.0
    res@vcMinDistanceF        = 0.03            ; thins arrows near pole  
    res@vcGlyphStyle          = "CurlyVector"     ; turns on curly vectors

    res@vcRefAnnoArrowLineColor   = "black"       ; black ref vector
    res@vcRefAnnoArrowUseVecColor = False

    ; TC center
    mres = True
    mres@gsMarkerIndex = 12
    mres@gsMarkerSizeF = 0.025
    mres@gsMarkerColor = "black"
    mres@gsMarkerThicknessF = 5.0
    ;plot(d) = gsn_csm_contour_map(wks,mslt,res)
    plotd = gsn_csm_vector_scalar_map(wks,ut,vt,mslt,res)
    delete(ut)
    delete(vt)
    delete(mslt)

    lonc = track(2*d, 4)
    latc = track(2*d, 5)
    print(track(2*d,:))
    dum = gsn_add_polymarker(wks, plotd, lonc, latc, mres)

    plot(0) = plotd
    delete(plotd)
    ;year  = toint(str_get_cols(sprinti("%0.10d",vdate), 0, 3))
    ;month = toint(str_get_cols(sprinti("%0.10d",vdate), 4, 5))
    ;day   = toint(str_get_cols(sprinti("%0.10d",vdate), 6, 7))
    ;hour  = toint(str_get_cols(sprinti("%0.10d",vdate), 8, 9))
    ;sec = cd_inv_calendar(year,month,day,hour,0,0, time@units, 0)
    ;newsec = sec
    ;newsec = sec + 6*60*60
    ;vdate  = toint(cd_string(newsec, "%Y%N%D%H"))
  end do
  ;resP	                   = True
  ;;resP@gsnPanelMainString	   = "hgt anomaly at "+level+"mb JMA ensemble forecast from "+sprinti("%0.10d",yyyymmddhh)
  ;gsn_panel(wks,plot,(/2,4/),resP)
  draw(plot(0))
  frame(wks)
end
