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
  load "config.ncl"
begin

;---Open netCDF file
  f = addfile(ncfile,"r")
;  if (d.eq.0) then
    fm= addfile(mfile, "r")
;  end if
;
; Read in Sea Surface Temperature Anomalies
; Read in U and V at 1000 mb [subscript 0]
; Determine the subscripts corresponding to Jan 88
;
  ;d = 0
  ;nlev = 6
  mmddhh = str_get_cols(yyyymmddhh, 4, 9)
  lon0 = 120.
  lon1 = 160.
  ;lon1 = 180.
  if (d.eq.0) then
      out = gfile
      lat0 = 0.
      ;lat0 = -10.
      lat1 = 40.
  else if (d.eq.2) then
      out = gfile+"_v12h"
      ;lat0 = 20.
      lat0 = 0.
      lat1 = 60.
  else if (d.eq.4) then
      out = gfile+"_v24h"
      ;lat0 = 20.
      lat0 = 0.
      lat1 = 40.
  else if (d.eq.9) then
      out = gfile+"_v54h"
      ;lat0 = 20.
      lat0 = 0.
      lat1 = 60.
  else
      out = gfile+"_v"
      lat0 = 20.
      ;lat0 = 0.
      lat1 = 60.
  end if
  end if
  end if
  end if
  ug       = f->ugrd(d,:,:,:)
  vg       = f->vgrd(d,:,:,:)
  ps       = f->pres_meansealev(d,:,:)
;  if (d.eq.0) then
  ps_m     = fm->msl(d,:,:)
  ps_m = ps_m * 0.01 ; Pa->hPa
;  end if
;; 00Z -> *(-1)
  ug = (/-1*ug/)
  vg = (/-1*vg/)
  ps = (/-1*ps/)
  te       = f->te(d,:,:,:)
  lat_uv   = f->lat
  lon_uv   = f->lon
  plev     = f->lev
  time     = f->time
  
  ;vl_pl = [/ug, vg, temp, spfh/]
  ;vl_sfc = [/ps, te/]
  ;cnt_pl = ListCount(vl_pl)
  ;cnt_sfc = ListCount(vl_sfc)
  ;print(cnt_pl)
  ;print(cnt_sfc)
;  out = "ensvsa_00Zv";"ensvsa_m1_t+wind_n_100912"
;  dev = "png"
;  wks = gsn_open_wks("png",out)
;  plot = new(6,graphic)
;  cmap = read_colormap_file("hotcold_18lev")
    
  vdate  = cd_string(time(d), "%Y%N%D%H")
;*************************************************
; hoskins filtering
;*************************************************
  if (lfilter) then
     ndim = dimsizes(ug)
     nlev = ndim(0)
     nlat = ndim(1)
     nlon = ndim(2)

     N = nlat
     if (nlon%2 .eq. 0)then
       N = min((/ nlat, (nlon+2)/2 /))
     else
       N = min((/ nlat, (nlon+1)/2 /))
     end if

     ;vl_plt = NewList("lifo")
     ;vl_sfct = NewList("lifo")
     T = 106
     NT = 63
     ;do i = 0, cnt_pl-1
     ;   var = vl_pl[i]
     ;   printVarSummary(var)
     ;   var_t = var
     a = new((/ nlev, nlat, N /),float)
     b = new((/ nlev, nlat, N /),float)
     ;; ug
     shaec(ug, a, b)
     ;tri_trunc(a, b, T)
     exp_tapersh(a, b, NT, 2)
     shsec(a, b, ug)
     copy_VarCoords(vg, ug)
     ug@long_name = "U component of wind"
     ug@units     = "m/s" 
     printMinMax(ug, 0)
     ;; vg
     shaec(vg, a, b)
     ;tri_trunc(a, b, T)
     exp_tapersh(a, b, NT, 2)
     shsec(a, b, vg)
     copy_VarCoords(ug, vg)
     vg@long_name = "V component of wind"
     vg@units     = "m/s" 
     printMinMax(vg, 0)
     ;   ListPush(vl_plt, var_t)
     delete(a)
     delete(b)
        ;delete(var)
        ;delete(var_t)
     ; end do
     ; print(ListCount(vl_plt))
     ; spfh = (/vl_plt[0]/)
     ; temp = (/vl_plt[1]/)
     ; vg = (/vl_plt[2]/)
     ; ug = (/vl_plt[3]/)
     ; delete(var)
     ; delete(var_t)
     ; do i = 0, cnt_sfc-1
     ;   var = vl_sfc[i]
     ;   printVarSummary(var)
     ;   var_t = var
     a = new((/ nlat, N /),float)
     b = new((/ nlat, N /),float)
     ;; ps
     shaec(ps, a, b)
     ;tri_trunc(a, b, T)
     exp_tapersh(a, b, NT, 2)
     shsec(a, b, ps)
     copy_VarCoords(te, ps)
     ps@long_name = "Mean sea level pressure"
     ps@units = "Pa"
     printMinMax(ps, 0)
     ;; te
     shaec(te, a, b)
     ;tri_trunc(a, b, T)
     exp_tapersh(a, b, NT, 2)
     shsec(a, b, te)
     copy_VarCoords(ps, te)
     te@long_name = "Total moist energy"
     te@units = "J/kg"
     printMinMax(te, 0)
     ;   printVarSummary(var_t)
     ;   ListPush(vl_sfct, var_t)
     delete(a)
     delete(b)
        ;delete(var)
        ;delete(var_t)
     ; end do
     ; te = (/vl_sfct[0]/)
     ; ps = (/vl_sfct[1]/)
   end if
;*************************************************
; vertical integration
;*************************************************
  ;sig = new(nlev,float)
  dp   = dpres_plevel_Wrap(plev, plev(nlev-1), plev(0), 0)
  dpall= plev(nlev-1) - plev(0)
  sig  = dble2flt( dp / dpall )
  ;;sig(0) = 1./7.
  ;;sig(1) = 10./49. + 1./7.
  ;;sig(2) = 3./7. + 4./49.
  ;sig(0) = 7./22.
  ;sig(1) = 1./2.
  ;sig(2) = 2./11.
  print(sig)
  ;sig = new(6, float)
  ;sig(0) = 3./22.
  ;sig(1) = 5./22.
  ;sig(2:4) = 2./11.
  ;sig(5) = 1./11.
  ;print(sig)
  ug = ug * conform_dims(dimsizes(ug), sig, 0)
  vg = vg * conform_dims(dimsizes(vg), sig, 0)
  up = dim_sum_n(ug, 0)
  up!0 = "lat"
  up!1 = "lon"
  up&lat = lat_uv
  up&lon = lon_uv
  up@_FillValue = default_fillvalue("float")
  vp = up
  vp = dim_sum_n(vg, 0)
;************************************************
; create points for box
;************************************************
 ypts = (/ 37.0, 37.0, 33.0, 33.0, 37.0/)
 xpts = (/137.0,142.0,142.0,137.0,137.0/)
; xpts = (/135.0,140.0,140.0,135.0,135.0/)
;  ypts = (/ 27.0, 27.0, 24.0, 24.0, 27.0/)
;  xpts = (/136.0,139.0,139.0,136.0,136.0/)
;************************************************
; add the box
;************************************************
  resp                  = True                      ; polyline mods desired
  resp@gsLineColor      = "black"                     ; color of lines
  resp@gsLineThicknessF = 4.0                       ; thickness of lines
;  resp@gsLineLabelString= "test"                    ; adds a line label string

   ;cmap = read_colormap_file("hotcold_18lev")
   cmap1 = read_colormap_file("gsdtol")
   cmap2 = read_colormap_file("gsltod")
   print(dimsizes(cmap1))
   print(dimsizes(cmap2))
   cmap = new((/32,4/),float)
   cmap(:15,:) = cmap1(:15,:)
   cmap(15:,:) = cmap2(:16,:)
   wks = gsn_open_wks(dev,out+"_slp+wind")
   reste                = True               ; plot mods desired
   reste@gsnDraw	      = False
   reste@gsnFrame	      = False
;   res@cnInfoLabelOn	      = False
     
   reste@gsnCenterString	   = "" ;"Total moist energy + slp"
   reste@gsnCenterStringOrthogonalPosF = 0.02
   reste@gsnCenterStringParallelPosF   = 0.55
   reste@gsnCenterStringFontHeightF    = 0.025
   reste@gsnLeftString	    = "" ;"J/kg/m^2"
   reste@gsnLeftStringFontHeightF = 0.02
   reste@gsnLeftStringParallelPosF = -0.1
   reste@gsnRightString	    = "hPa"
   reste@gsnRightStringFontHeightF = 0.02
   reste@gsnRightStringParallelPosF = 1.1
   reste@gsnAddCyclic	      = False
   reste@cnFillOn             = True               ; turn on color for contours
   reste@cnLinesOn            = True             ; turn off contour lines
   reste@cnLineLabelsOn       = False              ; turn off contour line labels
   reste@cnFillPalette        = cmap
   reste@gsnScalarContour     = True               ; contours desired
   reste@gsnContourZeroLineThicknessF = 0.0
   reste@gsnContourNegLineDashPattern = 2
   reste@gsnContourLineThicknessesScale = 2.0

   ;reste@mpLandFillColor      = "gray"            ; set land to be gray
   reste@mpFillOn             =  False            ; not land fill 
   reste@mpMinLonF            =  lon0               ; select a subregion
   reste@mpMaxLonF            =  lon1
   reste@mpMinLatF            =  lat0 
   reste@mpMaxLatF            =  lat1
   reste@lbOrientation            = "Vertical"     ; vertical label bar
   ;reste@lbOrientation            = "Horizontal"     ; vertical label bar
   ;reste@lbLabelStride = 2
   ;reste@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer

   reste@vpXF = 0.15
   reste@vpWidthF = 0.7

; note, when doing a subregion, NCL determines the range of the data from
; the full domain. If you wish to just consider the domain you are plotting,
; you must manually set those levels.

   ;reste@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
   reste@cnLevelSelectionMode = "ExplicitLevels"
;     res@cnMaxLevelCount	   = 10
   if (d.eq.0) then
;  d=0
      ;reste@cnLevels = (/1.,2.,4.,6.,8.,10.,15.,20.,25.,30.,35.,40.,45.,50.,60.,70.,80./)
      ;reste@cnLevels = (/1.,2.,3.,4.,5.,6.,8.,10.,12.5,15.,17.5,20./)
      ;reste@cnMinLevelValF       = -1.4              ; set min contour level
      ;reste@cnMaxLevelValF       = 1.4                ; set max contour level
      ;reste@cnLevelSpacingF      = 0.2               ; set contour spacing
      reste@cnLevels = (/-1.4,-1.0,-0.6,-0.2,0.2,0.6,1.0,1.4/)
      ;reste@cnMinLevelValF       = .1               ; set min contour level
      ;reste@cnMinLevelValF       = 10.               ; set min contour level
      ;reste@cnMaxLevelValF       = 100.0                 ; set max contour level
      ;reste@cnLevelSpacingF      = 10.               ; set contour spacing
      reste@vcRefMagnitudeF      = 10.0
   else if (d.eq.2) then
;  d=1
      ;reste@cnLevels = (/10.,20.,30.,40.,50.,75.,100.,150.,200.,250.,300.,400.,500.,600.,800.,1000./)
      ;reste@cnMinLevelValF       = -10.0               ; set min contour level
      ;reste@cnMaxLevelValF       = 10.                 ; set max contour level
      ;reste@cnLevelSpacingF      = 1.               ; set contour spacing
      reste@cnLevels = (/-8.,-6.,-4.,-2.,2.,4.,6.,8./)
      ;reste@cnMinLevelValF       = 10.               ; set min contour level
      ;reste@cnMaxLevelValF       = 100.0                 ; set max contour level
      ;reste@cnLevelSpacingF      = 10.               ; set contour spacing
      reste@vcRefMagnitudeF      = 20.0
   else if (d.eq.4) then
;  d=2
      ;reste@cnLevels = (/10.,20.,30.,40.,50.,75.,100.,150.,200.,250.,300.,400.,500.,600.,800.,1000./)
      ;reste@cnMinLevelValF       = -10.0               ; set min contour level
      ;reste@cnMaxLevelValF       = 10.                 ; set max contour level
      ;reste@cnLevelSpacingF      = 1.               ; set contour spacing
      reste@cnLevels = (/-8.,-6.,-4.,-2.,2.,4.,6.,8./)
      ;reste@cnMinLevelValF       = 10.               ; set min contour level
      ;reste@cnMaxLevelValF       = 100.0                 ; set max contour level
      ;reste@cnLevelSpacingF      = 10.               ; set contour spacing
      reste@vcRefMagnitudeF      = 20.0
   else 
;  d=6
      ;reste@cnLevels = (/10.,20.,30.,40.,50.,75.,100.,150.,200.,250.,300.,400.,500.,600.,800.,1000./)
      ;reste@cnMinLevelValF       = -50.0               ; set min contour level
      ;reste@cnMaxLevelValF       = 50.                 ; set max contour level
      ;reste@cnLevelSpacingF      = 5.               ; set contour spacing
      reste@cnLevels = (/-48.,-36.,-24.,-12.,12.,24.,36.,48./)
      ;reste@cnMinLevelValF       = 100.               ; set min contour level
      ;reste@cnMaxLevelValF       = 1000.0                 ; set max contour level
      ;reste@cnLevelSpacingF      = 100.               ; set contour spacing
      reste@vcRefMagnitudeF      = 100.0
   end if
   end if
   end if
   reste@vcRefLengthF              = 0.045           ; define length of vec ref
   reste@vcRefAnnoOrthogonalPosF   = -1.0            ; move ref vector
;   reste@vcRefAnnoArrowLineColor   = "black"         ; change ref vector color
;   reste@vcRefAnnoArrowUseVecColor = False           ; don't use vec color for ref

;     res@vcMinMagnitudeF	      = 1.0		  ; no effect with CurlyVector
;   reste@vcGlyphStyle            = "CurlyVector"     ; turn on curly vectors
   reste@vcFillArrowsOn       = True
   reste@vcMinFracLengthF     = 0.5
   reste@vcFillArrowEdgeColor = "black"
   ;reste@vcFillArrowFillColor = "green"
   reste@vcFillArrowFillColor = "gray"
   reste@vcMinDistanceF	      = 0.025
   reste@vcVectorDrawOrder       = "PostDraw"        ; draw vectors last

   wind = sqrt(up^2+vp^2)
   up = where(wind.lt.0.1, up@_FillValue, up)
   ;plote = gsn_csm_contour_map(wks, te, reste)
   plote = gsn_csm_vector_scalar_map(wks, up, vp, ps, reste)

   resps                = True               ; plot mods desired
   resps@gsnDraw	      = False
   resps@gsnFrame	      = False
;   res@cnInfoLabelOn	      = False     
   resps@gsnAddCyclic	      = False
   resps@cnFillOn             = False               ; turn on color for contours
   resps@cnLinesOn            = True              ; turn off contour lines
   resps@cnLineLabelsOn       = False              ; turn off contour line labels
   resps@cnLineColor          = "black"
   ;resps@cnLineColor          = "blue"
   resps@cnLineThicknessF     = 4.
   resps@gsnContourZeroLineThicknessF = 0. 	  ; doubles thickness of zero contour
   resps@gsnContourNegLineDashPattern = 1 	  ; sets negative contours to dash pattern 1
      
   resps@lbOrientation            = "Vertical"     ; vertical label bar
   resps@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer

; note, when doing a subregion, NCL determines the range of the data from
; the full domain. If you wish to just consider the domain you are plotting,
; you must manually set those levels.

   ;resps@cnLevelSelectionMode = "ExplicitLevels"
   resps@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
;     res@cnMaxLevelCount	   = 10
;   if (d.eq.0) then
;;  d=0
;      ;resps@cnLevels = (/1.,2.,3.,4.,5.,6.,8.,10.,12.5,15.,17.5,20./)
;      ;resps@cnMinLevelValF       = 10.               ; set min contour level
;      ;resps@cnMaxLevelValF       = 50.                 ; set max contour level
;      ;resps@cnLevelSpacingF      = 10.               ; set contour spacing
      resps@cnMinLevelValF       = 1003.            ; set min contour level
      resps@cnMaxLevelValF       = 1019.            ; set max contour level
      resps@cnLevelSpacingF      = 4.               ; set contour spacing
;   else if (d.eq.2) then
;;  d=2
;      ;resps@cnLevels = (/10.,20.,30.,40.,50.,75.,100.,150.,200.,250.,300.,400.,500.,600.,800.,1000./)
;      resps@cnMinLevelValF       = 10.               ; set min contour level
;      resps@cnMaxLevelValF       = 50.                 ; set max contour level
;      resps@cnLevelSpacingF      = 10.               ; set contour spacing
;      ;resps@cnMinLevelValF       = -20.0               ; set min contour level
;      ;resps@cnMaxLevelValF       = 20.                 ; set max contour level
;      ;resps@cnLevelSpacingF      = 2.               ; set contour spacing
;   else if (d.eq.4) then
;;  d=2
;      ;resps@cnLevels = (/10.,20.,30.,40.,50.,75.,100.,150.,200.,250.,300.,400.,500.,600.,800.,1000./)
;      resps@cnMinLevelValF       = 10.               ; set min contour level
;      resps@cnMaxLevelValF       = 50.                 ; set max contour level
;      resps@cnLevelSpacingF      = 10.               ; set contour spacing
;      ;resps@cnMinLevelValF       = -20.0               ; set min contour level
;      ;resps@cnMaxLevelValF       = 20.                 ; set max contour level
;      ;resps@cnLevelSpacingF      = 2.               ; set contour spacing
;   else 
;;  d=6
;      ;resps@cnLevels = (/10.,20.,30.,40.,50.,75.,100.,150.,200.,250.,300.,400.,500.,600.,800.,1000./)
;      resps@cnMinLevelValF       = 100.0               ; set min contour level
;      resps@cnMaxLevelValF       = 1000.                 ; set max contour level
;      resps@cnLevelSpacingF      = 100.               ; set contour spacing
;      ;resps@cnMinLevelValF       = -100.0               ; set min contour level
;      ;resps@cnMaxLevelValF       = 100.                 ; set max contour level
;      ;resps@cnLevelSpacingF      = 10.               ; set contour spacing
;   end if
;   end if
;   end if
   resps@gsnRightString   = ""    
   resps@gsnLeftString    = ""
   resps@gsnCenterString  = ""
;   if (d.eq.0) then
   plotp = gsn_csm_contour(wks, ps_m, resps)
   ;plotp = gsn_csm_contour(wks, te, resps)
   overlay(plote,plotp)
;   end if
; create array of dummy graphic variables. This is required, b/c each line
; must be associated with a unique dummy variable.
   dum = new(4,graphic)  
; draw each line separately. Each line must contain two points.
   do i = 0 , 3
      dum(i)=gsn_add_polyline(wks,plote,xpts(i:i+1),ypts(i:i+1),resp)      
   end do   
   ;plot(4) = plote
   draw(plote)
   frame(wks)

   ;resP	                   = True
   ;resP@gsnPanelMainString	   = "ensvsa-TE T&wind,slp,TE "+vdate+" "+orig
   ;gsn_panel(wks,plot,(/2,3/),resP)
end
