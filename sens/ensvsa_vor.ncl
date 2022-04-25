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

  ;d = 0
  ;nlev = 6
  mmddhh = str_get_cols(yyyymmddhh, 4, 9)
  lon0 = 130.
  lon1 = 150.
  ;lon1 = 180.
  if (d.eq.0) then
      out = gfile+"_vor"
      lat0 = 10.
      ;lat0 = -10.
      lat1 = 30.
  else if (d.lt.8) then
      hh = sprinti("%d",6*d)
      out = gfile+"_vor_v"+hh+"h"
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
;---Open netCDF file
  f = addfile("../sens/"+ncfile,"r")
  ug       = f->ugrd(d,:,:,:)
  vg       = f->vgrd(d,:,:,:)
  vor = uv2vrF(ug, vg)
  copy_VarCoords(ug, vor)
  vor@long_name = "Vorticity"
  vor@units     = "/s"
  printVarSummary(vor)
  printMinMax(vor,0)
;; 00Z -> *(-1)
  lat_uv   = f->lat
  lon_uv   = f->lon
  plev     = f->lev
  time     = f->time
;; ensemble mean
  fm = addfile(mfile,"r")
  u_m      = fm->u(d,:,:,:)
  v_m      = fm->v(d,:,:,:)
  wspd_m   = u_m
  wspd_m   = sqrt(u_m*u_m + v_m*v_m)
  wspd_m@long_name = "Wind speed"
  wspd_m@units     = "m/s"
  printVarSummary(wspd_m)
  printMinMax(wspd_m,0)
    
  vdate  = cd_string(time(d), "%Y%N%D%H")
;*************************************************
; hoskins filtering
;*************************************************
  if (lfilter) then
     ndim = dimsizes(vor)
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
     T = 213
     NT = 63
     ;do i = 0, cnt_pl-1
     ;   var = vl_pl[i]
     ;   printVarSummary(var)
     ;   var_t = var
     a = new((/ nlev, nlat, N /),float)
     b = new((/ nlev, nlat, N /),float)
     shaec(vor, a, b)
     tri_trunc(a, b, T)
     ;exp_tapersh(a, b, NT, 2)
     shsec(a, b, vor)
     copy_VarCoords(ug, vor)
     vor@long_name = "Vorticity"
     vor@units     = "/s"
     printVarSummary(vor)
     printMinMax(vor, 0)
     ;   ListPush(vl_plt, var_t)
     delete(a)
     delete(b)
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
     ;cmap = read_colormap_file("hotcold_18lev")
     cmap = read_colormap_file("sunshine_9lev")
     v = vor(ilev,:,:)
     hgt = plev(ilev)
     w = wspd_m({hgt},:,:)
;---Create plot
     res                = True               ; plot mods desired
     res@gsnDraw	      = False
     res@gsnFrame	      = False
;     res@cnInfoLabelOn	      = False
     
     res@gsnCenterString	   = "Wspd+Vor_"+hgt+"hPa"
     res@gsnCenterStringFontHeightF = 0.025
     res@gsnCenterStringOrthogonalPosF = 0.04
     ;res@gsnCenterStringParallelPosF   = 0.45
     res@gsnLeftStringFontHeightF = 0.02
     res@gsnLeftString   = "10~S~-5~N~ "+v@units
     res@gsnRightStringFontHeightF = 0.02
     res@gsnRightStringParallelPosF = 1.1    
     res@gsnRightString	    = w@units
     res@gsnAddCyclic	      = False
     res@cnFillOn             = True               ; turn on color for contours
     res@cnLinesOn            = False              ; turn off contour lines
     res@cnLineLabelsOn       = False              ; turn off contour line labels
     res@cnFillPalette        = cmap

     res@vpXF = 0.15
     res@vpWidthF = 0.7

     res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
     res@cnMaxLevelValF       = 45.0
     res@cnMinLevelValF       = 5.0
     res@cnLevelSpacingF      = 5.0
     
     resm = res
     ;resm@mpLandFillColor      = "gray"            ; set land to be gray
     resm@mpFillOn             = False
     resm@mpMinLonF            =  lon0               ; select a subregion
     resm@mpMaxLonF            =  lon1
     resm@mpMinLatF            =  lat0 
     resm@mpMaxLatF            =  lat1
     resm@lbOrientation            = "Vertical"     ; vertical label bar
     resm@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer
     plott = gsn_csm_contour_map_ce(wks,w,resm)

     res@cnFillOn             = False             ; turn on color for contours
     res@cnLinesOn            = True              ; turn off contour lines
     res@cnLineLabelsOn       = True
     res@cnLineColor          = "black"
     res@cnLineThicknessF     = 5.0
     res@gsnContourZeroLineThicknessF = 0. 	  ; doubles thickness of zero contour
     res@gsnContourNegLineDashPattern = 1 	  ; sets negative contours to dash pattern 1
     res@cnLevelSelectionMode = "ManualLevels"    ; set manual contour levels
     if (hgt.le.300) then
       res@cnMinLevelValF       = -10.               ; set min contour level
       res@cnMaxLevelValF       = 10.                 ; set max contour level
       res@cnLevelSpacingF      = .5               ; set contour spacing
     else 
       res@cnMinLevelValF       = -20.               ; set min contour level
       res@cnMaxLevelValF       = 20.                 ; set max contour level
       res@cnLevelSpacingF      = 1.               ; set contour spacing
     end if
      res@gsnRightString   = ""    
      res@gsnLeftString    = ""
      res@gsnCenterString  = ""
      
      v = v*1.0e5
      plotq = gsn_csm_contour(wks,v,res)
      overlay(plott,plotq)
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
      delete(plotq)
      delete(dum)
   end do
end
