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
  orig = "jma"
  dev = "png"
  init = "0912"
  if (init.eq."0900") then
    d_offset = 0
  else
    d_offset = 2
  end if
  a = addfile("/Volumes/dandelion/netcdf/tigge/2019/"+orig+"/glb_201910"+init+"_mean.nc","r")

  um   = a->u(:,:,:,:)
  vm   = a->v(:,:,:,:)
  tm   = a->t(:,:,:,:)
  qm   = a->q(:,:,:,:)
  psm  = a->msl(:,:,:)
  psm  = psm*0.01
  zm   = a->gh(:,:,:,:)
  lat  = a->latitude
  lon  = a->longitude
  plev = a->level
  time = a->time
  date = cd_string(time, "%Y%N%D%H")
  lon0 = 120.
  lon1 = 200.
  lat0 = -10.
  lat1 = 40.
  
  out = "strm_ensm_"+orig+"_"+init
  cmap = read_colormap_file("hotcold_18lev")
  plotnum = (/0,3,6,9,1,4,7,10,2,5,8,11/)
  do k = 3, 8
    u = um(:,k,:,:)
    v = vm(:,k,:,:)
    hgt = plev(k)
    sfvp = uv2sfvpF(u,v)
    sf = sfvp(0,:,:,:)
    copy_VarCoords(u,sf)
    sf@long_name = "stream function"
    sf@units     = "m^2/s"
    sf_zm = zonalAve(sf)
    sf = sf - conform(sf,sf_zm,(/0,1/))
    printVarSummary(sf)
    printMinMax(sf,0)
    sf = sf * 1.0e-6

    plot = new(dimsizes(plotnum),graphic)
    wks = gsn_open_wks(dev,out+"_"+hgt)
;---Create plot
    res                = True               ; plot mods desired
    res@gsnDraw	      = False
    res@gsnFrame	      = False
;     res@cnInfoLabelOn	      = False
     
    res@gsnCenterStringOrthogonalPosF = 0.02
    res@gsnLeftString	    = "10~S~6~N~"+sf@units
    res@gsnLeftStringFontHeightF = 0.02
    res@gsnRightString	    = ""
    res@gsnRightStringFontHeightF = 0.02
    res@gsnRightStringParallelPosF = 1.1    
    res@gsnAddCyclic	       = True
    res@cnFillOn             = True               ; turn on color for contours
    res@cnLinesOn            = True               ; turn off contour lines
    res@cnLineLabelsOn       = False              ; turn off contour line labels
    res@cnFillPalette        = cmap

;    res@mpLandFillColor      = "gray"            ; set land to be gray
    res@mpFillOn	      = False
    res@mpMinLonF            =  lon0               ; select a subregion
    res@mpMaxLonF            =  lon1
    res@mpCenterLonF         = (lon0+lon1)/2
    res@mpMinLatF            =  lat0
    res@mpMaxLatF            =  lat1

    ;res@lbOrientation            = "Vertical"     ; vertical label bar
    ;res@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer
    ;res@vpXF = 0.15
    ;res@vpWidthF = 0.7

; note, when doing a subregion, NCL determines the range of the data from
; the full domain. If you wish to just consider the domain you are plotting,
; you must manually set those levels.

    res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
    res@cnMinLevelValF       = -9.               ; set min contour level
    res@cnMaxLevelValF       = 9.                ; set max contour level
    res@cnLevelSpacingF      = 1.                 ; set contour spacing
    res@gsnContourNegLineDashPattern = 2
    res@gsnContourZeroLineThicknessF = 3.

    i = d_offset
    do d=0,dimsizes(plotnum)-1
      res@gsnCenterString = date(d)
      plot(plotnum(i)) = gsn_csm_contour_map_ce(wks,sf(d,:,:),res)
      i = i + 1
      if (i.ge.dimsizes(plotnum)) then
        break
      end if
    end do
    delete(sf)
    gsn_panel(wks,plot,(/4,3/),False)
    delete(plot)
    delete(wks)
   end do
end
