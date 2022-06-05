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
  yyyymmddhh = 2019100912
  sst = "_mgdsst"
  a = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100"+sst+"/fcst_surf_asia_"+yyyymmddhh+".nc","r")
  b = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100"+sst+"/fcst_phy3mpp_"+yyyymmddhh+".nc","r")
  wks = gsn_open_wks("png","slp+hr_cloud"+yyyymmddhh+sst)
  plot = new(16,graphic)

  id = 2
  pmax = 850
  pmin = 250

  lat     = b->lat
  lon     = b->lon
  lev     = b->level({pmax:pmin})
  print(lev)
  dp      = dpres_plevel_Wrap(lev, lev({pmax}), lev({pmin}), 0)
  dpall = lev({pmax})-lev({pmin})
  dp = dp/dpall
  print(dp)
  time1h  = a->time
  time6h  = b->time
  date1h  = cd_string(time1h, "%N%D %H:%M")
  date6h  = cd_string(time6h, "%N%D %H:%M")
  ntime   = dimsizes(time6h)
  print(date6h)
  factor = 1.0 / 24.0 ; per day -> per hour
  do d = 0, ntime-1
;
; Read in Sea Surface Temperature Anomalies
; Read in U and V at 1000 mb [subscript 0]
; Determine the subscripts corresponding to Jan 88
;
     ind_p    = ind(date1h.eq.date6h(d))
     ps       = a->PSEA(ind_p,:,:)
     ps@units = "hPa"
     dTcv     = b->aHRCV(d,{pmax:pmin},:,:)
     dTlc     = b->aHRLC(d,{pmax:pmin},:,:)
     dT       = dTcv
     dT       = dTcv + dTlc
     ctime    = time6h(d)
     date     = cd_string(ctime, "%N%D %H:%M")
     ;dT = dTp - dTm
     dT = dT*conform(dT, dp, 0)
     dTdt = dT(0,:,:)
     dTdt@standard_name = "Heating Rate (Cloud)"
     dTdt@units         = "K/h"
     dTdt = dim_sum_n(dT, 0) * factor

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
   resP	                   = True
   ;resP@gsnPanelMainString	   = "SLP + "+dT@long_name+" 600-200 hPa"
   resP@gsnPanelLabelBar     = True
   gsn_panel(wks,plot,(/4,4/),resP)
   delete(plot)
end
