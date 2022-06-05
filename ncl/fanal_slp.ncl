;----------------------------------------------------------------------
; fanal_1.ncl
;
; Concepts illustrated:
;   - Calculating zonal Fourier information of a monthly climatology
;---------------------------------------------------------
; Calculate zonal Fourier information of a monthly climatology           
; .  Amplitude, Phase, Percent Variance at each latitude
;---------------------------------------------------------
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
 
   orig = "ukmo"
   if (orig.eq."ecmwf") then
   member = 50
   nlev=3
   etype="TE"
   else if (orig.eq."jma") then
   member = 26
   nlev=6
   etype="TE"
   else if (orig.eq."ncep") then
   member = 20
   nlev=3
   etype="TE"
   else
   member = 17
   nlev=3
   etype="dTE"
   end if
   end if
   end if
   diri = "/Volumes/dandelion/netcdf/tigge/2019/"+orig+"/rotate/"
   film = "np_2019100912_mean.nc"
   fm   = addfile(diri+film,"r")               ; mean
   lat  =fm->lat
   slpm =fm->PRES_meansealevel(0,{87:90},:)
   filw = "../sens/"+etype+"-"+orig+"-2019100912_nlev"+sprinti("%d",nlev)+".txt" ; EnSVSA weights
   wgts = asciiread(filw,(/10,member/),"float")
   wgt_mem = wgts(0,:)
   if (orig.ne."ncep") then
   wgt_mem = (/-1*wgt_mem/)
   end if
   print(wgt_mem)
   ;exit
   ;wgt_mem = (/-0.1779, 0.4010, 0.4165, -0.09494, 0.1487,\
   ;-0.06161, -0.1471, -0.1457, -0.1755, -0.07958, -0.1726,\
   ;0.1764, -0.1650, -0.03709, -0.1821, 0.1902,  -0.2722,\
   ;0.3560, 0.1550, -0.1771, 0.1133, -0.1444, 0.09761,\
   ;-0.1003, 0.1932, -0.1148/) ; Ensemble members' weight for EnSVSA of jma
   ;; for plot
   res = True
   res@gsnPolar = "NH"
   res@mpMinLatF = 87
   res@mpOutlineOn = False
   res@cnFillOn = True
   res@cnLevelSelectionMode = "ManualLevels"
   res@cnMinLevelValF = -2.5
   res@cnMaxLevelValF = 2.5
   res@cnLevelSpacingF = 0.25
   res@cnFillPalette = "hotcold_18lev"
   
   slp1 = slpm ; EnSVSA
   slp1 = 0.0
do m = 1, member
;   m=1
   mem = sprinti("%0.2d",m)
   fili = "np_2019100912_"+mem+".nc"
   out = "prtb_"+orig+"/finfo_0912_"+orig+"_"+mem
   f    = addfile(diri+fili,"r")               ; member
   slp  = f->PRES_meansealevel(0,{87:90},:)          ; (lat,lon)   
   slp = slp - slpm
   printVarSummary(slp)                        ; (3,360)
   printMinMax(slp,0)
   ;if (orig .eq. "jma") then
   slp1 = slp1 + wgt_mem(m-1)*slp
   ;end if
   wks = gsn_open_wks("png",out)
   slp = slp*0.01
   plot = gsn_csm_contour_map_polar(wks,slp,res)
   delete(wks)
   delete(plot)
   
   dims = dimsizes(slp)
   nlat = dims(0)                              ; 3
   mlon = dims(1)                              ; 360
 
;****************************************
; Calculate harmonics           
;   finfo(0,:,:,nhar) - amplitude of each harmonic
;   finfo(1,:,:,nhar) - phase of each harmonic
;   finfo(2,:,:,nhar) - percent variance explained by each harmonic 
;****************************************
   ; extract wave-number 0 component
   slp_zave = zonalAve(slp)
   slp = slp - conform(slp, slp_zave, 0)

   nhar  = mlon/2                           ; # of harmonics
   finfo = fourier_info(slp, nhar, 1.)      ; fourier_info works on *right* dim => lon
                                            ; (3,nhar)
   copy_VarCoords(slp(:,0), finfo(0,:,0)) ; trick
   finfo!0 = "ftype"
   finfo!2 = "harmonic"
 ;;finfo@long_name = "Z500: amp, phase, % variance""
   finfo_avg = finfo(:,0,:)
   finfo_avg = dim_avg_n(finfo, 1)
   print(finfo_avg(:,0))                   ; [ftype|3]x[lat|3]x[harmonic|3]
   ;recalculate variance including wave-number 0
   amp_0   = abs(avg(slp_zave))
   amp_sum = sum(finfo_avg(0,:))
   amp_sum = amp_sum + amp_0
   var     = finfo_avg(0,:) / amp_sum
   var0    = amp_0 / amp_sum
   finfo_avg(2,:) = (/var/)

   fmt_str = "%s"
   dlist = NewList("lifo")
   ListPush(dlist,"#amplitude(Pa) phase(deg) variance(%)")
   write_table(out,"w",dlist,fmt_str)
   delete(dlist)
   fmt_str = str_concat(conform_dims(3,"%16.8e",-1))
   dlist = NewList("lifo")
   ListPush(dlist,(/var0/)) ;variance
   ListPush(dlist,(/0.0/))  ;phase
   ListPush(dlist,(/avg(slp_zave)/));amplitude
   write_table(out,"a",dlist,fmt_str)
   delete(dlist)
   dlist = NewList("lifo")
   do j=2,0,1
   ListPush(dlist,(/finfo_avg(j,:)/))
   end do
   write_table(out,"a",dlist,fmt_str)
   delete(dlist)
end do
;if (orig .eq. "jma") then
; Calculate harmonics for EnSVSA
   out = "prtb_"+orig+"/finfo_0912_"+orig+"_EnSVSA"
   printVarSummary(slp1)
   printMinMax(slp1,0)
   ; check EnSVSA
   wks = gsn_open_wks("png",out)
   res@cnMinLevelValF = -2.5
   res@cnMaxLevelValF = 2.5
   res@cnLevelSpacingF = 0.25
   slp1 = slp1*0.01
   plot = gsn_csm_contour_map_polar(wks,slp1,res)
   ; extract wave-number 0 component
   slp_zave = zonalAve(slp1)
   slp1 = slp1 - conform(slp1, slp_zave, 0)

   nhar  = mlon/2                           ; # of harmonics
   finfo = fourier_info(slp1, nhar, 1.)      ; fourier_info works on *right* dim => lon
                                            ; (3,nhar)
   copy_VarCoords(slp1(:,0), finfo(0,:,0)) ; trick
   finfo!0 = "ftype"
   finfo!2 = "harmonic"
 ;;finfo@long_name = "Z500: amp, phase, % variance""
   finfo_avg = finfo(:,0,:)
   finfo_avg = dim_avg_n(finfo, 1)
   print(finfo_avg(:,0))                   ; [ftype|3]x[lat|3]x[harmonic|3]
   ;recalculate variance including wave-number 0
   amp_0   = abs(avg(slp_zave))
   amp_sum = sum(finfo_avg(0,:))
   amp_sum = amp_sum + amp_0
   var     = finfo_avg(0,:) / amp_sum
   var0    = amp_0 / amp_sum
   finfo_avg(2,:) = (/var/)

   fmt_str = "%s"
   dlist = NewList("lifo")
   ListPush(dlist,"#amplitude(Pa) phase(deg) variance(%)")
   write_table(out,"w",dlist,fmt_str)
   delete(dlist)
   fmt_str = str_concat(conform_dims(3,"%16.8e",-1))
   dlist = NewList("lifo")
   ListPush(dlist,(/var0/)) ;variance
   ListPush(dlist,(/0.0/))  ;phase
   ListPush(dlist,(/avg(slp_zave)/));amplitude
   write_table(out,"a",dlist,fmt_str)
   delete(dlist)
   dlist = NewList("lifo")
   do j=2,0,1
   ListPush(dlist,(/finfo_avg(j,:)/))
   end do
   write_table(out,"a",dlist,fmt_str)
   delete(dlist)
;end if
exit
;****************************************
; Create plot
;****************************************
   nharPlt = 2
   plot    = new (3,"graphic")
                             
   maxAmp  = max(finfo(0,:,:))
  
   pltType = "png"                               ; send graphics to PNG file
   pltName = "fanal_0900_mem1"
  
   wks     = gsn_open_wks(pltType, pltName)    
  
   res                      = True          ; individual plot
   res@gsnDraw              = False
   res@gsnFrame             = False
   res@xyLineThicknessF     = 2.0
   res@xyLineColors         = (/"blue","red"/)   ; change line colors 
   res@xyDashPattern        = 0                  ; Make curves all solid
  
; add a legend
   res@pmLegendDisplayMode    = "Always"              ; turn on legend
   res@pmLegendSide           = "Top"                 ; Change location of 
   res@pmLegendParallelPosF   = .85                   ; move units right
   res@pmLegendOrthogonalPosF = -0.85                 ; more neg = down
   res@pmLegendWidthF         = 0.13                  ; Change width and
   res@pmLegendHeightF        = 0.075                 ; height of legend.
   res@lgLabelFontHeightF     = .02                   ; change font height
   res@lgPerimOn              = False                 ; no box around
   res@xyExplicitLegendLabels = (/"HAR 1", "HAR 2"/) ; legend labels
  
   np = -1
   np = np+1
   res@gsnLeftString  = "Amplitude Harmonics 1 & 2"
   res@gsnRightString = slp@units
   res@trXMaxF        = maxAmp   
   work               = finfo(ftype|0, harmonic|0:1 ,lat|:)  ; reorder
   plot(np)           = gsn_csm_xy(wks,work,slp&lat,res)
  
   np = np+1
   res@gsnLeftString  = "Phase Harmonics 1 & 2"
   res@gsnRightString = "deg"
   res@trXMaxF        = 360.0      
   work               = finfo(ftype|1, harmonic|0:1 ,lat|:)  ; reorder
   plot(np)           = gsn_csm_xy(wks,work,slp&lat,res)
  
   np = np+1
   res@gsnLeftString  = "% Variances Harmonics 1 & 2"
   res@gsnRightString = "%"
   res@trXMaxF        = 100.0      
   work               = finfo(ftype|2, harmonic|0:1 ,lat|:)  ; reorder
   plot(np)           = gsn_csm_xy(wks,work,slp&lat,res)
  
   resP                     = True          ; panel resources
   resP@gsnMaximize         = True
   resP@gsnPanelMainString  = "SLP perturbation : member 1"
   gsn_panel(wks,plot,(/1,3/),resP)
