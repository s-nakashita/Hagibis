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
  ;hh = 12
  ;cmap = read_colormap_file("precip3_16lev") ; all level
  ;cmap = read_colormap_file("BlueDarkRed18") ; low level
  cmap = read_colormap_file("MPL_YlOrRd") ; low level
  ;cmap = cmap(::-1,:) ; for dyn
  wks = gsn_open_wks("png","slp+dTdt_mst")
  plot = new(16,graphic)

  id = 0
  pmax = 850
  pmin = 250
  do dd = 9, 12
  idate = 2019100000 + dd*100
  yyyymmdd = str_get_cols(sprinti("%0.10d",idate), 0, 7)
  a = addfile("/Volumes/dandelion/netcdf/nasa/merra2/MERRA2_400.tavg1_2d_slv_Nx."+yyyymmdd+".nc4","r")
  b = addfile("/Volumes/dandelion/netcdf/nasa/merra2/MERRA2_400.tavg3_3d_tdt_Np."+yyyymmdd+".nc4","r")
  lat     = b->lat
  lon     = b->lon
  lev     = b->lev({pmax:pmin})
  print(lev)
  dp      = dpres_plevel_Wrap(lev, lev({pmax}), lev({pmin}), 0)
  dpall = lev({pmax})-lev({pmin})
  dp = dp/dpall
  print(dp)
  time2d  = a->time
  time3d  = b->time
  date2d    = cd_string(time2d, "%N%D %H:%M")
  date3d    = cd_string(time3d, "%N%D %H:%M")
  print(date3d)
  
  do d = 0, 6, 2
;
; Read in Sea Surface Temperature Anomalies
; Read in U and V at 1000 mb [subscript 0]
; Determine the subscripts corresponding to Jan 88
;
     ind_p    = ind(date2d.eq.date3d(d))
     ps       = a->SLP(ind_p,:,:);PRMSL_meansealevel(:,180:340,180:360)
     ind_p2   = ind(date2d.eq.date3d(d+1))
     ps       = ps + a->SLP(ind_p2,:,:)
     ps       = ps*1e-2/2.
     ps@units = "hPa"
     dTp      = b->DTDTTOT(d,{pmax:pmin},:,:)
     dTp      = dTp + b->DTDTTOT(d+1,{pmax:pmin},:,:)
     dTp      = dTp * 0.5
     dTm      = b->DTDTMST(d,{pmax:pmin},:,:)
     dTm      = dTm + b->DTDTMST(d+1,{pmax:pmin},:,:)
     dTm      = dTm * 0.5
     dT       = dTm
     ctime    = time3d(d)
     ctime    = (ctime + time3d(d+1))/2
     date     = cd_string(ctime, "%N%D %H:%M")
     ;dT = dTp - dTm
     dT = dT*conform(dT, dble2flt(dp), 0)
     dTdt = ps
     dTdt@standard_name = dT@standard_name
     dTdt@units = "K/h"
     dTdt = dim_sum_n(dT, 0) * 3600.0

;---Create plot
     res                  = True               ; plot mods desired
     res@gsnDraw	      = False
     res@gsnFrame	      = False
;     res@cnInfoLabelOn	      = False
     
     res@gsnCenterString	   = date
     res@gsnRightString     = dTdt@units
     res@gsnLeftString	   = ps@units
     res@gsnCenterStringOrthogonalPosF = 0.02
     res@gsnAddCyclic	      = False
     res@cnFillOn             = True               ; turn on color for contours
     res@cnLinesOn            = False               ; turn off contour lines
     res@cnLineLabelsOn       = False              ; turn off contour line labels
     res@gsnScalarContour     = True               ; (for gsn_csm_vector_scalar_map)draw contours
     res@cnFillPalette        = cmap

     res@mpFillOn             = False 
;     res@mpLandFillColor      = "gray"            ; set land to be gray
     res@mpMinLonF            =  130.               ; select a subregion
     res@mpMaxLonF            =  160.
     res@mpMinLatF            =  10. 
     res@mpMaxLatF            =  40.

     res@lbLabelBarOn         = False
     ;res@lbOrientation            = "Vertical"     ; vertical label bar
     ;res@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer

; note, when doing a subregion, NCL determines the range of the data from
; the full domain. If you wish to just consider the domain you are plotting,
; you must manually set those levels.

     res@cnLevelSelectionMode = "ExplicitLevels"     ; set manual contour levels
     res@cnLevels = (/0.0,0.1,0.2,0.3,0.4,0.5,0.75,1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0,7.0,8.0/)
     ;res@cnLevels = (/-2.5,-2.0,-1.5,-1.0,-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5,1.0,1.5,2.0,2.5/)
     ;res@cnMaxLevelCount	   = 20
     ;res@cnMinLevelValF       = 0.5               ; set min contour level
     ;res@cnMaxLevelValF       = 5.0                ; set max contour level
     ;res@cnLevelSpacingF      = 0.5               ; set contour spacing

     plotdz = gsn_csm_contour_map_ce(wks,dTdt,res)
     resp = True
     resp@gsnDraw	      = False
     resp@gsnFrame	      = False
     resp@cnLevelSelectionMode = "ManualLevels"
     resp@cnMinLevelValF       = 1000.0               ; set min contour level
     resp@cnMaxLevelValF       = 1020.0                 ; set max contour level
     resp@cnLevelSpacingF      = 4.0               ; set contour spacing
     resp@cnFillOn             = False               ; turn on color for contours
     resp@cnLinesOn            = True                ; turn off contour lines
     resp@cnLineLabelsOn        = True
     ;resp@cnLineColor          = "blue"
     resp@gsnRightString   = ""    
     resp@gsnLeftString    = ""
     resp@gsnCenterString  = ""
     plotz = gsn_csm_contour(wks,ps,resp)
     overlay(plotdz,plotz)
     plot(id) = plotdz
     id = id + 1
     delete(plotdz)
     delete(plotz)

   end do
   end do
   resP	                   = True
   ;resP@gsnPanelMainString	   = "SLP + "+dT@long_name+" 600-200 hPa"
   resP@gsnPanelLabelBar     = True
   gsn_panel(wks,plot,(/4,4/),resP)
   delete(plot)
end
