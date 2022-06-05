;**********************************************************************
; Function to retrieve the location of a plot so we can calculate
; a location for a title on the left side of the plot.
;**********************************************************************
function get_title_position(plot)
local vpy, vpx, vph
begin
  getvalues plot
      "vpYF"       : vpy
      "vpXF"       : vpx
      "vpHeightF"  : vph
  end getvalues
  pos = (/vpx*0.75,vpy-(vph/2.)/)
  return(pos)
end

;**********************************************************************
; Main code
;**********************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
begin 
   ;Reading Historical data

   sst = "est"
   init = "2019100900"
   ;fh1     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100/fcst_p_asia_2019100900.nc", "r")
   fh     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_"+sst+"/fcst_p_asia_"+init+".nc","r")
   tfile  = "../pytrack/track"+init+"_gsm_tl959_"+sst+".txt"
   ;fh3     = addfile("/Volumes/dandelion/GSMJob/Jobwk_Tl959L100_mgdsst/fcst_surf_2019100900.nc", "r")
   lmax = 100
   u   = fh->U(:, {:lmax}, :, :)
   v   = fh->V(:, {:lmax}, :, :)
   t3h = fh->time
   date = cd_string(t3h, "%N%D %H:%M")
   ntim = dimsizes(date)
   ;track data
   data = readAsciiTable(tfile, 7, "float", 0)
  ;print(data(:,4:6))
   npoints = dimsizes(data(:,0))
   dd   = toint(data(:,2))
   hh   = toint(data(:,3))
   lonc = data(:,4)
   latc = data(:,5)
   yyyy = dd
   mm   = dd
   mn   = hh
   sc   = hh
   yyyy = 2019
   mm   = 10
   mn   = 0
   sc   = 0
   tracktime = cd_inv_calendar(yyyy,mm,dd,hh,mn,sc,t3h@units,0)
   tracktime!0 = "time"
   ;print(tracktime)

   pltDir = "./"                ; where output will be written
   pltName= "GSMtl959_"+sst+"_vor_h-lon_2019100900"
   pltType= "png"               ; ps, eps, png, pdf, x11 .....

;****************************************************************************
; Graphic resources
;****************************************************************************
    ;ver = 4
  ds = 0
  de = min((/ds+15, 28/)) ;max 28
  i = 0
  pltPath = pltDir+pltName+"_"+cd_string(t3h(ds), "%D%H")+"-"+cd_string(t3h(de), "%D%H")
  wks     = gsn_open_wks(pltType, pltPath) 
  gsn_define_colormap(wks, "GMT_polar")
  nrow = 4
  ncol = 4
  plot = new(nrow*ncol , "graphic")
  pi = atan(1.0)*4.0
  do d=ds, de
    indc = ind(tracktime.eq.t3h(d))
    lon0 = lonc(indc)
    lat0 = latc(indc)
    rlon = 200.0 / 6.371e3 / cos(lat0*pi/180.0) * 180.0 / pi
    print(rlon)
    lonmin = lon0 - rlon
    lonmax = lon0 + rlon
;****************************************************************************
; Calculate vorticity
;****************************************************************************
    vor = u(d,:,:,:)
    vor = uv2vr_cfd(u(d,:,:,:), v(d,:,:,:), u&lat, u&lon, 3)
    vor@standard_name = "Relative vorticity"
    vor@units = "/s"
;****************************************************************************
; Set plot resources
;****************************************************************************
    res1                            = True     ; plot mods desired
    res1@gsnDraw                    = False    ; will panel plots later
    res1@gsnFrame                   = False

    res1@cnFillOn                   = True     ; turn on color fill
    res1@cnLinesOn                  = True    ; turn off contour lines
    res1@cnLineLabelsOn             = False    ; turn off contour line labels
    ;res1@cnFillPalette              = "GMT_polar"
    res1@cnLevelSelectionMode      = "ManualLevels"  ; manually specify contour levels
    res1@cnMinLevelValF            = -2.e-3          ; min level
    res1@cnMaxLevelValF            = 2.e-3           ; max level
    res1@cnLevelSpacingF           = 2.e-4           ; contour interval
    res1@gsnContourZeroLineThicknessF = 3.0

;**********************************************************************
; We are turning the labelbar off for all but the bottom row of plots
;**********************************************************************
    res1@lbLabelBarOn               = False

    ; Set resources for tickmarks
    res1@gsnStringFont              = "helvetica-bold"  ; will be used for gsnCenterString later
    res1@gsnCenterStringFontHeightF = 0.04              ; ditto
    res1@gsnLeftString              = ""                ; remove the left and right substrings
    res1@gsnRightString             = ""
    res1@gsnCenterString            = date(d)

;    res1@gsnScalarContour = True
;    res1@vcRefMagnitudeF           = 40.0            ; define vector ref mag
;    res1@vcRefLengthF              = 0.045           ; define length of vec ref
;    res1@vcRefAnnoOrthogonalPosF   = -1.0            ; move ref vector
;    res1@vcRefAnnoArrowLineColor   = "black"         ; change ref vector color
;    res1@vcRefAnnoArrowUseVecColor = False           ; don't use vec color for ref
;
;    res1@vcGlyphStyle            = "LineArrow"
;    res1@vcLineArrowColor        = "black"           ; change vector color
;    res1@vcLineArrowThicknessF   = 2.0               ; change vector thickness
;    res1@vcMinDistanceF	         = 0.05
;    res1@vcMinMagnitudeF	       = 5.0
;    res1@vcVectorDrawOrder       = "PostDraw"        ; draw vectors last

;****************************************************************************
; Create the plots for the panel
;****************************************************************************
    plot(i) = gsn_csm_pres_hgt(wks, vor(:,{lat0},{lonmin:lonmax}), res1)
    i = i + 1
;    delete(ud)
;    delete(vd)
  end do
;****************************************************************************
; Draw the panel and label bars
;****************************************************************************
    resP                = True
    resP@gsnMaximize    = True
    resP@gsnFrame       = False
    resP@gsnPanelLeft   = 0.005
    resP@gsnPanelBottom = 0.08 
    resP@gsnPanelSave   = True      ; Save the state of the paneled plots so we can
                                    ; query for text strings locations.
   ;resP@gsnPanelMainString = "Land Only: Tmax: 2 Colorbars"
    resP@gsnPanelLabelBar = True

    gsn_panel(wks, plot, (/nrow,ncol/), resP)  ; plots
    frame(wks)
    delete(wks)
end
