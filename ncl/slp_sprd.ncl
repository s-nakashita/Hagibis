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
  hh = "00"
  orig = "jma"
  a = addfile("../../netcdf/tigge/2019/"+orig+"/glb_20191009"+hh+"_mean.nc","r")
  b = addfile("../../netcdf/tigge/2019/"+orig+"/glb_20191009"+hh+"_spr.nc","r")
;
; Read in Sea Surface Temperature Anomalies
; Read in U and V at 1000 mb [subscript 0]
; Determine the subscripts corresponding to Jan 88
;
  ps    = a->msl(:,:,:);PRMSL_meansealevel(:,180:340,180:360)
  spr   = b->msl(:,:,:)
  ;ps12    = b->PRES_meansealevel(:,:,:)
  lat   = a->latitude
  lon   = a->longitude
  time  = a->time
  date  = cd_string(time, "%Y%N%D%H")
  print(date)

  ;nd = 8
  ds = ind(date.eq."2019101012")
  de = ind(date.eq."2019101012")
  ;ds = 0
  ;de = 7
  nd = de - ds + 1
  vdate = 2019100900 + toint(hh)
  year  = toint(str_get_cols(sprinti("%0.10d",vdate), 0, 3))
  month = toint(str_get_cols(sprinti("%0.10d",vdate), 4, 5))
  day   = toint(str_get_cols(sprinti("%0.10d",vdate), 6, 7))
  hour  = toint(str_get_cols(sprinti("%0.10d",vdate), 8, 9))
  sec = cd_inv_calendar(year,month,day,hour,0,0, time@units, 0)
  newsec = sec
  newsec = sec + 6*ds
  vdate  = toint(cd_string(newsec, "%Y%N%D%H"))
  wks = gsn_open_wks("png","slp_sprd_09"+hh+"_"+orig+"_"+date(ds)+"to"+date(de))
  plot = new(nd,graphic)
  cmap = read_colormap_file("spread_15lev")
    
  i = 0
  do d = ds,de
;---Create plot
     res                      = True               ; plot mods desired
     res@gsnDraw	      = False
     res@gsnFrame	      = False
;     res@cnInfoLabelOn	      = False
     
     res@gsnCenterString	   = sprinti("%0.10d",vdate)
     res@gsnLeftString	   = ""
     res@gsnCenterStringOrthogonalPosF = 0.02
     res@gsnAddCyclic	      = False
     res@cnFillOn             = True               ; turn on color for contours
     res@cnLinesOn            = False               ; turn off contour lines
     res@cnLineLabelsOn       = False              ; turn off contour line labels
     ;res@gsnScalarContour     = True               ; (for gsn_csm_vector_scalar_map)draw contours
     res@cnFillPalette        = cmap

     res@mpFillOn             = False 
;     res@mpLandFillColor      = "gray"            ; set land to be gray
     res@mpMinLonF            =  110.               ; select a subregion
     res@mpMaxLonF            =  180.
     res@mpMinLatF            =  -10. 
     res@mpMaxLatF            =  40.

     res@lbOrientation            = "Vertical"     ; vertical label bar
     res@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer

; note, when doing a subregion, NCL determines the range of the data from
; the full domain. If you wish to just consider the domain you are plotting,
; you must manually set those levels.

     res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
     res@cnMaxLevelCount	   = 20
     res@cnMinLevelValF       = 40.0               ; set min contour level
     res@cnMaxLevelValF       = 200.0                 ; set max contour level
     res@cnLevelSpacingF      = 10.0               ; set contour spacing

     plotdz = gsn_csm_contour_map_ce(wks,spr(d,:,:),res)
     resp = True
     resp@gsnDraw	      = False
     resp@gsnFrame	      = False
     resp@cnLevelSelectionMode = "ManualLevels"
     resp@cnMinLevelValF       = 100000.0               ; set min contour level
     resp@cnMaxLevelValF       = 102000.0                 ; set max contour level
     resp@cnLevelSpacingF      = 100.0               ; set contour spacing
     resp@cnFillOn             = False               ; turn on color for contours
     resp@cnLinesOn            = True                ; turn off contour lines
     resp@cnLineLabelsOn        = True
     resp@cnLineColor          = "black"
     resp@gsnRightString   = ""    
     resp@gsnLeftString    = ""
     resp@gsnCenterString  = ""
     plotz = gsn_csm_contour(wks,ps(d,:,:),resp)
     overlay(plotdz,plotz)
     plot(i) = plotdz
     delete(plotdz)
     delete(plotz)
     year  = toint(str_get_cols(sprinti("%0.10d",vdate), 0, 3))
     month = toint(str_get_cols(sprinti("%0.10d",vdate), 4, 5))
     day   = toint(str_get_cols(sprinti("%0.10d",vdate), 6, 7))
     hour  = toint(str_get_cols(sprinti("%0.10d",vdate), 8, 9))
     sec = cd_inv_calendar(year,month,day,hour,0,0, time@units, 0)
     newsec = sec
     newsec = sec + 6 ;*60*60
     vdate  = toint(cd_string(newsec, "%Y%N%D%H"))
     i = i + 1
   end do
   ;resP	                   = True
   ;resP@gsnPanelMainString	   = "slp spread "+orig+" 1010"+hh
   ;gsn_panel(wks,plot,(/2,nd/2/),resP)
   ;delete(plot)
   draw(plot(0))
   frame(wks)
end
