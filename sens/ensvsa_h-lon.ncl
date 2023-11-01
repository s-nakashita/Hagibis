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
  fm= addfile(mfile,"r")
;
; Read in Sea Surface Temperature Anomalies
; Read in U and V at 1000 mb [subscript 0]
; Determine the subscripts corresponding to Jan 88
;
  ;d = 0
  ;nlev = 6
  mmddhh = str_get_cols(yyyymmddhh, 4, 9)
  lon0 = 135.
  lon1 = 160.
  ;lon1 = 180.
  if (d.eq.0) then
  ;    latc = 15.
      out = gfile+"_h-lon_lat"+sprinti("%0.2d",toint(latc))
  else
  ;    latc = 35.
      out = gfile+"_h-lon_v"+sprinti("%d",6*d)+"h_lat"+sprinti("%0.2d",toint(latc))
  end if
  lat0 = latc - 5.0
  lat1 = latc + 5.0
;; 00Z -> *(-1)
  ug      = f->ugrd(d,:,:,:)
  ug = (/-1*ug/)
  vg      = f->vgrd(d,:,:,:)
  vg = (/-1*vg/)
  temp    = f->t(d,:,{latc},:)
  temp = (/-1*temp/)
  if (ntype.eq."TE") then
  spfh	  = f->q(d,:,{latc},:)
  spfh = (/-1*spfh/)
  end if
  gh      = f->gh(d,:,{latc},:)
  gh   = (/-1*gh/)
  ps      = f->pres_meansealev(d,:,:)
  ps   = (/-1*ps/)
  lon      = f->lon
  plev     = f->lev
  plev@units = "hPa"
  plev_pa = plev
  plev_pa = plev * 1.0e2
  plev_pa@units = "Pa"
  time     = f->time
;;; ensemble mean
  u_m    = fm->u(d,:,{latc},:)
  temp_m = fm->t(d,:,{latc},:)
  spfh_m = fm->q(d,:,{latc},:)
  gh_m   = fm->gh(d,:,{latc},:)
  plev_m = int2flt(fm->level)
  ps_m    = fm->msl(d,:,:) ;Pa
  ;pbase = ps
  ;pbase = pbase * 1.0e2
  pbase = ps_m ;+ ps
  ;pbase@units = ps_m@units
  printVarSummary(pbase)
  printMinMax(pbase,0)
  tbase = temp
  do k=0,nlev-1
    tbase({plev(k)},:) = (/temp_m({plev(k)},:)/) ;+ temp({plev(k)},:)
  end do
  printVarSummary(tbase)
  printMinMax(tbase,0)
;;; calculate equivalent potential temperature
  ;plev_m = plev_m * 1.0e2 ;Pa
  ;;; spfh -> rh
  ;rh = relhum(temp_m, spfh_m, conform(temp_m, plev_m, 0))
  ;;; t, rh -> tlcl
  ;tlcl = tlcl_rh_bolton(temp_m,rh,(/1,1/))
  ;;; theta_e
  ;theta_e = pot_temp_equiv_tlcl(conform(temp_m, plev_m, 0), temp_m, tlcl, spfh_m, (/1,1,0,1/))
  ;printMinMax(theta_e,0)
  ;plev_m = plev_m * 1.0e-2
;; Compute divergence
  div = uv2dvF(ug, vg)
  copy_VarCoords(ug, div)
  printVarSummary(div)
  printMinMax(div(:,{latc},:),0)
;; Compute pressure layer thicknesses
  dp = dpres_plevel(plev_pa, pbase, 0.0, 0)
  printVarSummary(dp)
  printMinMax(dp,0)
;; Compute d(omega) 
;; (Atmospheric science an introductory survey, 2nd ed., Wallace and Hobbs, eq (7.39b))
  omega = div
  omega = div * dp
  omega@long_name = "d(omega)"
  omega@units     = plev_pa@units + "/s"
  omega@method    = "derived via kinematic method: div*dp"
  printVarSummary(omega)
  printMinMax(omega(:,{latc},:),0)
;; integrate d(omega)
  tmp = omega(lat|:,lon|:,lev|::-1)
  oBot = 0.0 ; explicit
  do kl=0,nlev-1
    omega(nlev-1-kl,:,:) = oBot + dim_sum(tmp(:,:,0:kl))
  end do
  omega@long_name = "pressure velocity"
  printVarSummary(omega)
  printMinMax(omega(:,{latc},:),0)
;; omega -> w
  wvel = omega(:,{latc},:)
  wvel = omega_to_w(omega(:,{latc},:), conform(tbase,dble2flt(plev_pa),0), tbase)
  wvel = wvel * 1.0e2
  wvel@long_name = "vertical velocity"
  wvel@units     = "cm/s"
  printVarSummary(wvel)
  printMinMax(wvel,0)
;; T -> theta
  theta = pot_temp(plev, temp, 0, False)

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
     ;; div
     shaec(div, a, b)
     ;tri_trunc(a, b, T)
     exp_tapersh(a, b, NT, 2)
     shsec(a, b, div)
     copy_VarCoords(ug, div)
     div@long_name = "Divergence"
     div@units     = "/s"
     printVarSummary(div)
     printMinMax(div, 0)
     ;; temp
     shaec(temp, a, b)
     ;tri_trunc(a, b, T)
     exp_tapersh(a, b, NT, 2)
     shsec(a, b, temp)
     copy_VarCoords(ug, temp)
     temp@long_name = "Temparature"
     temp@units     = "K" 
     printVarSummary(temp)
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
     printVarSummary(spfh)
     printMinMax(spfh, 0)
     end if
     ;   ListPush(vl_plt, var_t)
     ;; gh
     shaec(gh, a, b)
     ;tri_trunc(a, b, T)
     exp_tapersh(a, b, NT, 2)
     shsec(a, b, gh)
     copy_VarCoords(ug, gh)
     gh@long_name = "Geopotential height"
     gh@units     = "gpm" 
     printMinMax(gh, 0)
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
     ;   delete(var)
     ;   delete(var_t)
     ; end do
     ; te = (/vl_sfct[0]/)
     ; ps = (/vl_sfct[1]/)
   end if
  resp                  = True                      ; polyline mods desired
  resp@gsLineColor      = "black"                     ; color of lines
  resp@gsLineThicknessF = 3.0                       ; thickness of lines
;  resp@gsLineLabelString= "test"                    ; adds a line label string

  wks = gsn_open_wks(dev,out)
  ;cmap = read_colormap_file("temp_19lev")
  ;;cmap = cmap(::-1,:) ;;reverse
  cmap1 = read_colormap_file("gsdtol")
  cmap2 = read_colormap_file("gsltod")
  print(dimsizes(cmap1))
  print(dimsizes(cmap2))
  cmap = new((/32,4/),float)
  cmap(:15,:) = cmap1(:15,:)
  cmap(15:,:) = cmap2(:16,:)
  ;u = ug(:,{lon0:lon1})
  ;v = vg(:,{lon0:lon1})
  ;t = temp(:,{lon0:lon1})
  ;t = dim_avg_n(theta(:,:,{lon0:lon1}),1)
  t = temp(:,{lon0:lon1})
  if (ntype.eq."TE") then
  ;q = dim_avg_n(spfh(:,:,{lon0:lon1}),1)
  q = spfh(:,{lon0:lon1})
  q = q*1000
  q!0 = "level"
  q&level = plev
  q!1 = "lon"
  q&lon = lon({lon0:lon1})
  end if
  ;w = wvel(:,{lon0:lon1})
  dv = div(:,{latc},{lon0:lon1})
  ;z = gh(:,{lon0:lon1})
  t!0 = "level"
  t&level = plev
  t!1 = "lon"
  t&lon = lon({lon0:lon1})
  ;w!0 = "level"
  ;w&level = plev
  ;w!1 = "lon"
  ;w&lon = lon({lon0:lon1})
  dv!0 = "level"
  dv&level = plev
  dv!1 = "lon"
  dv&lon = lon({lon0:lon1})
  ;u = w
  ;u = 0.0 ;dummy
  ;z!0 = "level"
  ;z&level = plev
  ;te = dim_avg_n(theta_e(:,:,{lon0:lon1}),1)
  ;te = theta_e(:,{lon0:lon1})
  ;te!0 = "level"
  ;te&level = plev_m
  ;te!1 = "lon"
  ;te&lon = lon({lon0:lon1})
  ;z = dim_avg_n(gh_m(:,:,{lon0:lon1}),1)
  ;z = gh_m(:,{lon0:lon1})
  ;z!0 = "level"
  ;z&level = plev_m
  ;z!1 = "lon"
  ;z&lon = lon({lon0:lon1})
  ;u = u_m(:,{lon0:lon1})
;---Create plot
  res                = True               ; plot mods desired
  res@gsnDraw	      = False
  res@gsnFrame	      = False
;     res@cnInfoLabelOn	      = False
  ;res@gsnCenterString	   = "lat="+lat0+"N-"+lat1+"N"
  res@gsnCenterString	   = "lat="+latc+"N"
  res@gsnCenterStringFontHeightF = 0.025
  res@gsnCenterStringOrthogonalPosF = 0.04
  ;res@gsnCenterStringParallelPosF   = 0.45
  if (ntype.eq."TE") then
  res@gsnLeftString	    = "g/kg"
  else
  res@gsnLeftString     = ""
  end if
  res@gsnLeftStringFontHeightF = 0.02
  res@gsnRightString   = "K"
  res@gsnRightStringFontHeightF = 0.02
  res@gsnRightStringParallelPosF = 1.1    
  res@gsnAddCyclic	      = False
  ;res@gsnScalarContour     = True               ; contours desired
  
  res@tmYRMode = "Automatic"
  res@lbOrientation            = "Vertical"     ; vertical label bar
  ;res@lbOrientation            = "Horizontal"     ; vertical label bar
  res@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer

  res@vpXF = 0.15
  res@vpWidthF = 0.65
  
  res@cnFillOn             = True               ; turn on color for contours
  res@cnFillPalette        = cmap
  res@cnLinesOn            = True              ; turn off contour lines
  res@cnLineThicknessF     = 2.
  res@gsnContourZeroLineThicknessF = 0. 	  ; doubles thickness of zero contour
  res@gsnContourNegLineDashPattern = 2 	  ; sets negative contours to dash pattern 1
  ;res@cnLevelSelectionMode = "ManualLevels"    ; set manual contour levels
  res@cnLevelSelectionMode = "ExplicitLevels"    ; set manual contour levels
  if (d.lt.4) then
;  d=0
    ;res@cnMinLevelValF       = -1.8              ; set min contour level
    ;res@cnMaxLevelValF       = 1.8              ; set max contour level
    ;res@cnLevelSpacingF      = .2               ; set contour spacing
    res@cnLevels = (/-1.6,-1.2,-0.8,-0.4,0.4,0.8,1.2,1.6/)
  ;  res@vcRefMagnitudeF      = 10.0
  else if (d.lt.8) then 
;  d=2
    ;res@cnMinLevelValF       = -5.               ; set min contour level
    ;res@cnMaxLevelValF       = 5.                 ; set max contour level
    ;res@cnLevelSpacingF      = .5               ; set contour spacing
    res@cnLevels = (/-4.,-3.,-2.,-1.,1.,2.,3.,4./)
  ;  res@vcRefMagnitudeF      = 20.0
  else 
;  d=6
    ;res@cnMinLevelValF       = -9.               ; set min contour level
    ;res@cnMaxLevelValF       = 9.                 ; set max contour level
    ;res@cnLevelSpacingF      = 1.               ; set contour spacing
    res@cnLevels = (/-8.,-6.,-4.,-2.,2.,4.,6.,8./)
  ;  res@vcRefMagnitudeF      = 100.0
  end if
  end if
  ;end if
  ;resv=res
  ;resv@vcRefAnnoOrthogonalPosF = -0.165
  ;resv@pmLabelBarOrthogonalPosF = -0.05
  ;resv@vcRefMagnitudeF = 1.0
  ;resv@vcRefLengthF    = 0.045
  ;resv@vcMapDirection  = False
  ;plott = gsn_csm_vector_scalar_map_ce(wks,u,v,t,res)
  plott = gsn_csm_pres_hgt(wks,t,res)
  ;plott = gsn_csm_pres_hgt_vector(wks,t,u,w,resv)

if (ntype.eq."TE") then
  res@cnFillOn             = False             ; turn on color for contours
  res@cnLinesOn            = True              ; turn off contour lines
  res@cnLineLabelsOn       = True
  res@cnLineColor          = "black"
  res@cnLineThicknessF     = 4.0
  res@cnInfoLabelOrthogonalPosF = -0.085
  res@gsnContourZeroLineThicknessF = 0. 	  ; doubles thickness of zero contour
  res@gsnContourNegLineDashPattern = 1 	  ; sets negative contours to dash pattern 1
;     res@cnMaxLevelCount	   = 10
  res@cnLevelSelectionMode = "ManualLevels"    ; set manual contour levels
  if (d.lt.4) then
;;  d=0
;    res@cnMinLevelValF       = -5.               ; set min contour level
;    res@cnMaxLevelValF       = 5.                 ; set max contour level
;    res@cnLevelSpacingF      = 1.               ; set contour spacing
;  elseif (d.lt.8) then
;  d=2
    res@cnMinLevelValF       = -10.               ; set min contour level
    res@cnMaxLevelValF       = 10.                 ; set max contour level
    res@cnLevelSpacingF      = 2.               ; set contour spacing
  else 
;  d=6
    res@cnMinLevelValF       = -10.               ; set min contour level
    res@cnMaxLevelValF       = 10.                 ; set max contour level
    res@cnLevelSpacingF      = 2.               ; set contour spacing
  end if
; when using the gsn_csm plot templates and overlay, it is necessary to
; turn off the automatic label strings so that you don't get text on 
; top of text
  res@gsnRightString   = ""    
  res@gsnLeftString    = ""
  res@gsnCenterString  = ""
  plotq = gsn_csm_pres_hgt(wks,q,res)
  overlay(plott,plotq)
end if
;;; div
;  res@cnLinesOn            = True
;  res@cnFillOn             = False
;  res@cnLineLabelsOn       = True
;  res@cnInfoLabelOn        = False
;  res@cnLineColor          = "black"
;  res@cnLineThicknessF     = 5.0
;  res@gsnContourZeroLineThicknessF = 0. 	  ; doubles thickness of zero contour
;  res@gsnContourNegLineDashPattern = 1 	  ; sets negative contours to dash pattern 1
;  res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
;  res@cnMinLevelValF = -5.
;  res@cnMaxLevelValF = 5.
;  res@cnLevelSpacingF = 1.
;  dv = -dv*1.0e5
;  plotz = gsn_csm_pres_hgt(wks,dv,res)
;  overlay(plott,plotz)
;  if (d.eq.0) then
;;  d=0
;    res@cnMinLevelValF       = -20.               ; set min contour level
;    res@cnMaxLevelValF       = 20.                 ; set max contour level
;    res@cnLevelSpacingF      = 4.               ; set contour spacing
;  else if (d.eq.4) then
;;  d=2
;    res@cnMinLevelValF       = -20.               ; set min contour level
;    res@cnMaxLevelValF       = 20.                 ; set max contour level
;    res@cnLevelSpacingF      = 4.               ; set contour spacing
;  else 
;;  d=6
;    res@cnMinLevelValF       = -150.               ; set min contour level
;    res@cnMaxLevelValF       = 150.                 ; set max contour level
;    res@cnLevelSpacingF      = 25.               ; set contour spacing
;  end if
;  end if
;  plotz = gsn_csm_pres_hgt(wks,z,res)
;;; theta_e
;  res@cnMinLevelValF = 300.0
;  res@cnMaxLevelValF = 360.0
;  res@cnLevelSpacingF = 4.0
;  plotz = gsn_csm_pres_hgt(wks,te,res)
;;; gh
;  res@cnMinLevelValF = 1000.0
;  res@cnMaxLevelValF = 8000.0
;  res@cnLevelSpacingF = 1000.0
;  plotz = gsn_csm_pres_hgt(wks,z,res)
;;; u
;  res@cnMinLevelValF = -10.0
;  res@cnMaxLevelValF = 10.0
;  res@cnLevelSpacingF = 2.0
;  plotz = gsn_csm_pres_hgt(wks,u,res)
;  overlay(plott,plotz)
  draw(plott)
  frame(wks)
end