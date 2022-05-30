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
  f = addfile("../sens/"+ncfile,"r")
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
  else if (d.lt.8) then
      hh = sprinti("%d",6*d)
      out = gfile+"_v"+hh+"h"
      ;lat0 = 20.
      lat0 = 0.
      lat1 = 40.
;  else if (d.eq.9) then
;      out = gfile+"_v54h"
;      ;lat0 = 20.
;      lat0 = 0.
;      lat1 = 60.
  else
      hh = sprinti("%d",6*d)
      out = gfile+"_v"+hh+"h"
      lat0 = 20.
      ;lat0 = 0.
      lat1 = 60.
;  end if
  end if
  end if
;; 00Z -> *(-1)
  ug       = f->ugrd(d,:,:,:)
  ug = (/-1*ug/)
  vg       = f->vgrd(d,:,:,:)
  vg = (/-1*vg/)
  temp     = f->t(d,:,:,:)
  temp = (/-1*temp/)
  if (ntype .eq. "TE") then
    spfh	  = f->q(d,:,:,:)
    spfh = (/-1*spfh/)
  end if
  ps       = f->pres_meansealev(d,:,:)
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
     ;; temp
     shaec(temp, a, b)
     ;tri_trunc(a, b, T)
     exp_tapersh(a, b, NT, 2)
     shsec(a, b, temp)
     copy_VarCoords(ug, temp)
     temp@long_name = "Temparature"
     temp@units     = "K" 
     printMinMax(temp, 0)
   if (ntype.eq."TE") then
     ;; spfh
     shaec(spfh, a, b)
     ;tri_trunc(a, b, T)
     exp_tapersh(a, b, NT, 2)
     shsec(a, b, spfh)
     copy_VarCoords(ug, spfh)
     spfh@long_name = "Specific humidity"
     spfh@units     = "kg/kg" 
     printMinMax(spfh, 0)
     ;   ListPush(vl_plt, var_t)
   end if
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
  resp@gsLineThicknessF = 3.0                       ; thickness of lines
;  resp@gsLineLabelString= "test"                    ; adds a line label string

  do k = 0, nlev-1
      ilev = nlev-k-1
      wks = gsn_open_wks(dev,out+"_"+sprinti("%0.3d",toint(plev(ilev))))
      cmap = read_colormap_file("hotcold_18lev")
      u = ug(ilev,:,:)
      v = vg(ilev,:,:)
      t = temp(ilev,:,:)
      if (ntype.eq."TE") then
      q = spfh(ilev,:,:)
      q = q*1000
      end if
      hgt = plev(ilev)
      ;wks = gsn_open_wks(dev,out+"_"+hgt)
;---Create plot
      res                = True               ; plot mods desired
      res@gsnDraw	      = False
      res@gsnFrame	      = False
;     res@cnInfoLabelOn	      = False
     
      res@gsnCenterString	   = "T,q+wind_"+hgt+"hPa"
      res@gsnCenterStringFontHeightF = 0.025
      res@gsnCenterStringOrthogonalPosF = 0.04
      ;res@gsnCenterStringParallelPosF   = 0.45
      res@gsnLeftString	    = "K"
      res@gsnLeftStringFontHeightF = 0.02
      res@gsnRightString   = "g/kg"
      res@gsnRightStringFontHeightF = 0.02
      res@gsnRightStringParallelPosF = 1.1    
      res@gsnAddCyclic	      = False
      res@cnFillOn             = True               ; turn on color for contours
      res@cnLinesOn            = False              ; turn off contour lines
      res@cnLineLabelsOn       = False              ; turn off contour line labels
      res@cnFillPalette        = cmap

      res@vpXF = 0.15
      res@vpWidthF = 0.7

      res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels

      resv = res
      ;res@mpLandFillColor      = "gray"            ; set land to be gray
      resv@mpFillOn             = False
      resv@mpMinLonF            =  lon0               ; select a subregion
      resv@mpMaxLonF            =  lon1
      resv@mpMinLatF            =  lat0 
      resv@mpMaxLatF            =  lat1
      resv@lbOrientation            = "Vertical"     ; vertical label bar
      resv@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer
;     res@cnMaxLevelCount	   = 10
   if (d.le.4) then
;  d=0
      resv@cnMinLevelValF       = -2.              ; set min contour level
      resv@cnMaxLevelValF       = 2.                ; set max contour level
      resv@cnLevelSpacingF      = .2               ; set contour spacing
      resv@vcRefMagnitudeF      = 10.0
   else if (d.le.8) then 
;  d=2
      resv@cnMinLevelValF       = -5.               ; set min contour level
      resv@cnMaxLevelValF       = 5.                 ; set max contour level
      resv@cnLevelSpacingF      = .5               ; set contour spacing
      resv@vcRefMagnitudeF      = 20.0
   else 
;  d=6
      resv@cnMinLevelValF       = -10.               ; set min contour level
      resv@cnMaxLevelValF       = 10.                 ; set max contour level
      resv@cnLevelSpacingF      = 1.               ; set contour spacing
      resv@vcRefMagnitudeF      = 100.0
   end if
   end if
   ;end if
      resv@gsnScalarContour     = True               ; contours desired
      resv@vcRefLengthF              = 0.045           ; define length of vec ref
      resv@vcRefAnnoOrthogonalPosF   = -1.0            ; move ref vector
      resv@vcRefAnnoArrowLineColor   = "black"         ; change ref vector color
      resv@vcRefAnnoArrowUseVecColor = False           ; don't use vec color for ref

;     resv@vcMinMagnitudeF	      = 1.0		  ; no effect with CurlyVector
      resv@vcGlyphStyle            = "CurlyVector"     ; turn on curly vectors
      resv@vcLineArrowColor        = "black"           ; change vector color
      resv@vcLineArrowThicknessF   = 2.0               ; change vector thickness
      resv@vcMinDistanceF	      = 0.02
      resv@vcVectorDrawOrder       = "PostDraw"        ; draw vectors last

      wind = sqrt(u^2+v^2)
      u = where(wind.lt.0.1, u@_FillValue, u)
  
      plott = gsn_csm_vector_scalar_map_ce(wks,u,v,t,resv)

   if (ntype.eq."TE") then
      res@cnFillOn             = False             ; turn on color for contours
      res@cnLinesOn            = True              ; turn off contour lines
      res@cnLineLabelsOn       = True
      res@cnLineColor          = "blue"
      res@gsnContourZeroLineThicknessF = 0. 	  ; doubles thickness of zero contour
      res@gsnContourNegLineDashPattern = 1 	  ; sets negative contours to dash pattern 1
      res@cnLevelSelectionMode = "ManualLevels"    ; set manual contour levels
;     res@cnMaxLevelCount	   = 10
      if (d.le.4) then
;  d=0
      res@cnMinLevelValF       = -5.               ; set min contour level
      res@cnMaxLevelValF       = 5.                 ; set max contour level
      res@cnLevelSpacingF      = 1.               ; set contour spacing
   else if (d.le.8) then
;  d=2
      res@cnMinLevelValF       = -10.               ; set min contour level
      res@cnMaxLevelValF       = 10.                 ; set max contour level
      res@cnLevelSpacingF      = 2.               ; set contour spacing
   else 
;  d=6
      res@cnMinLevelValF       = -20.               ; set min contour level
      res@cnMaxLevelValF       = 20.                 ; set max contour level
      res@cnLevelSpacingF      = 4.               ; set contour spacing
   end if
   end if
   ;end if
; when using the gsn_csm plot templates and overlay, it is necessary to
; turn off the automatic label strings so that you don't get text on 
; top of text
      res@gsnRightString   = ""    
      res@gsnLeftString    = ""
      res@gsnCenterString  = ""
      
      plotq = gsn_csm_contour(wks,q,res)
      overlay(plott,plotq)
      delete(plotq)
   end if
; create array of dummy graphic variables. This is required, b/c each line
; must be associated with a unique dummy variable.
      dum = new(4,graphic)  
; draw each line separately. Each line must contain two points.
      do i = 0 , 3
         dum(i)=gsn_add_polyline(wks,plott,xpts(i:i+1),ypts(i:i+1),resp)      
      end do   
      ;plot(k) = plott
      draw(plott)
      frame(wks)
      delete(plott)
      delete(dum)
   end do

   ;delete(cmap)
   ;cmap = read_colormap_file("spread_15lev")
   wks = gsn_open_wks(dev,out+"_sfc")
   reste                = True               ; plot mods desired
   reste@gsnDraw	      = False
   reste@gsnFrame	      = False
;   res@cnInfoLabelOn	      = False
     
   ;reste@gsnCenterString	   = "Total moist energy + slp"
   reste@gsnCenterString	   = "slp"
   reste@gsnCenterStringOrthogonalPosF = 0.02
   reste@gsnCenterStringParallelPosF   = 0.55
   reste@gsnCenterStringFontHeightF    = 0.025
   ;reste@gsnLeftString	    = "J/kg/m^2"
   reste@gsnLeftString	    = ""
   reste@gsnLeftStringFontHeightF = 0.02
   reste@gsnLeftStringParallelPosF = -0.1
   reste@gsnRightString	    = "hPa"
   ;reste@gsnRightString	    = ""
   reste@gsnRightStringFontHeightF = 0.02
   reste@gsnRightStringParallelPosF = 1.1
   reste@gsnAddCyclic	      = False
   reste@cnFillOn             = True               ; turn on color for contours
   reste@cnLinesOn            = False             ; turn off contour lines
   reste@cnLineLabelsOn       = False              ; turn off contour line labels
   reste@cnFillPalette        = cmap

   reste@mpLandFillColor      = "gray"            ; set land to be gray
   reste@mpMinLonF            =  lon0               ; select a subregion
   reste@mpMaxLonF            =  lon1
   reste@mpMinLatF            =  lat0 
   reste@mpMaxLatF            =  lat1
   reste@lbOrientation            = "Vertical"     ; vertical label bar
   ;reste@lbOrientation            = "Horizontal"     ; vertical label bar
   ;reste@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer

   reste@vpXF = 0.15
   reste@vpWidthF = 0.7

; note, when doing a subregion, NCL determines the range of the data from
; the full domain. If you wish to just consider the domain you are plotting,
; you must manually set those levels.

   reste@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
   ;reste@cnLevelSelectionMode = "ExplicitLevels"
;     res@cnMaxLevelCount	   = 10
   if (d.le.4) then
;  d=0
      ;reste@cnLevels = (/1.,2.,4.,6.,8.,10.,15.,20.,25.,30.,35.,40.,45.,50.,60.,70.,80./)
      ;reste@cnLevels = (/1.,2.,3.,4.,5.,6.,8.,10.,12.5,15.,17.5,20./)
      reste@cnMinLevelValF       = -2.5               ; set min contour level
      reste@cnMaxLevelValF       = 2.5                 ; set max contour level
      reste@cnLevelSpacingF      = 0.25               ; set contour spacing
      ;reste@cnMinLevelValF       = .1               ; set min contour level
      ;reste@cnMinLevelValF       = 10.               ; set min contour level
      ;reste@cnMaxLevelValF       = 100.0                 ; set max contour level
      ;reste@cnLevelSpacingF      = 10.               ; set contour spacing
   else if (d.le.8) then
;  d=2
      ;reste@cnLevels = (/10.,20.,30.,40.,50.,75.,100.,150.,200.,250.,300.,400.,500.,600.,800.,1000./)
      reste@cnMinLevelValF       = -10.0               ; set min contour level
      reste@cnMaxLevelValF       = 10.                 ; set max contour level
      reste@cnLevelSpacingF      = 1.               ; set contour spacing
      ;reste@cnMinLevelValF       = 10.               ; set min contour level
      ;reste@cnMaxLevelValF       = 100.0                 ; set max contour level
      ;reste@cnLevelSpacingF      = 10.               ; set contour spacing
   else 
;  d=6
      ;reste@cnLevels = (/10.,20.,30.,40.,50.,75.,100.,150.,200.,250.,300.,400.,500.,600.,800.,1000./)
      reste@cnMinLevelValF       = -50.0               ; set min contour level
      reste@cnMaxLevelValF       = 50.                 ; set max contour level
      reste@cnLevelSpacingF      = 5.               ; set contour spacing
      ;reste@cnMinLevelValF       = 100.               ; set min contour level
      ;reste@cnMaxLevelValF       = 1000.0                 ; set max contour level
      ;reste@cnLevelSpacingF      = 100.               ; set contour spacing
   end if
   end if
   ;plote = gsn_csm_contour_map(wks, te, reste)
   plote = gsn_csm_contour_map(wks, ps, reste)

   resps                = True               ; plot mods desired
   resps@gsnDraw	      = False
   resps@gsnFrame	      = False
;   res@cnInfoLabelOn	      = False
     
   resps@gsnCenterString	   = ""
   resps@gsnCenterStringOrthogonalPosF = 0.02
   resps@gsnLeftString	    = ""
   resps@gsnAddCyclic	      = False
   resps@cnFillOn             = False               ; turn on color for contours
   resps@cnLinesOn            = True              ; turn off contour lines
   resps@cnLineLabelsOn       = False              ; turn off contour line labels
   resps@cnLineColor          = "black"
   ;resps@cnLineColor          = "blue"
   resps@cnLineThicknessF     = 2.
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
   if (d.le.4) then
;  d=0
      ;resps@cnLevels = (/1.,2.,3.,4.,5.,6.,8.,10.,12.5,15.,17.5,20./)
      resps@cnMinLevelValF       = 20.               ; set min contour level
      resps@cnMaxLevelValF       = 200.                 ; set max contour level
      resps@cnLevelSpacingF      = 20.               ; set contour spacing
      ;resps@cnMinLevelValF       = -5.0               ; set min contour level
      ;resps@cnMaxLevelValF       = 5.                 ; set max contour level
      ;resps@cnLevelSpacingF      = 0.5               ; set contour spacing
   else if (d.le.8) then
;  d=2
      ;resps@cnLevels = (/10.,20.,30.,40.,50.,75.,100.,150.,200.,250.,300.,400.,500.,600.,800.,1000./)
      resps@cnMinLevelValF       = 20.0               ; set min contour level
      resps@cnMaxLevelValF       = 200.                 ; set max contour level
      resps@cnLevelSpacingF      = 20.               ; set contour spacing
      ;resps@cnMinLevelValF       = -20.0               ; set min contour level
      ;resps@cnMaxLevelValF       = 20.                 ; set max contour level
      ;resps@cnLevelSpacingF      = 2.               ; set contour spacing
   else 
;  d=6
      ;resps@cnLevels = (/10.,20.,30.,40.,50.,75.,100.,150.,200.,250.,300.,400.,500.,600.,800.,1000./)
      resps@cnMinLevelValF       = 100.0               ; set min contour level
      resps@cnMaxLevelValF       = 1000.                 ; set max contour level
      resps@cnLevelSpacingF      = 100.               ; set contour spacing
      ;resps@cnMinLevelValF       = -100.0               ; set min contour level
      ;resps@cnMaxLevelValF       = 100.                 ; set max contour level
      ;resps@cnLevelSpacingF      = 10.               ; set contour spacing
   end if
   end if
   ;end if
   
   ;plotp = gsn_csm_contour(wks, ps, resps)
   ;plotp = gsn_csm_contour(wks, te, resps)
   ;overlay(plote,plotp)
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
