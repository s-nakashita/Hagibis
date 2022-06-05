load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
begin
;; hourly data
;f = addfile("/Volumes/dandelion/netcdf/era5/2019/10/cape_msl.nc","r")
;; monthly averaged hourly data
f = addfile("/Volumes/dandelion/netcdf/era5/2019/cape_msl.nc","r")

cape = f->cape
msl = f->msl
msl = msl*1e-2 ; Pa->hPa
t = f->time
date = cd_string(t, "%N%D%H")
month = cd_string(t, "%N")
hour = cd_string(t, "%H")
;nd = dimsizes(date)
print(date)
hours = (/"00","03","06","09","12","15","18","21"/)
;nd = dimsizes(hours)
ind10 = ind(month.eq."10")
nd = dimsizes(ind10)
print(ind10)
;; time average
;mcape = cape(0,:,:)
;mcape = dim_avg_n(cape, 0)
;mmsl = msl(0,:,:)
;mmsl = dim_avg_n(msl, 0)
;printMinMax(mcape, 0)
;printMinMax(mmsl, 0)
;wks = gsn_open_wks("png", "cape+msl_ave"+date(0)+"-"+date(nd-1))
do d = 0, nd-1
;wks = gsn_open_wks("png", "cape+msl_"+date(d))
;wks = gsn_open_wks("png", "cape_ave1006-1014_"+hours(d)+"Z")
wks = gsn_open_wks("png", "cape_oct"+hour(ind10(d))+"Z")
;indh = ind(hour.eq.hours(d))
;print(date(indh))
;mcape = cape(0,:,:)
;mcape = dim_avg_n(cape(indh,:,:), 0)
;printMinMax(mcape, 0)
res1 = True
res1@gsnDraw = False
res1@gsnFrame = False

res1@cnLevelSelectionMode = "ManualLevels"
res1@cnMinLevelValF = 125.
res1@cnMaxLevelValF = 2000.
res1@cnLevelSpacingF = 125.
res1@cnFillOn = True
res1@cnLinesOn = False
res1@cnFillPalette   = "spread_15lev"

res1@lbOrientation = "Vertical"
;res1@lbLabelStride = 2

res1@mpMinLonF = 100.
res1@mpMaxLonF = 180.
res1@mpMinLatF = -20.
res1@mpMaxLatF = 40.
res1@mpFillOn = False

res1@vpXF = 0.15
res1@vpWidthF = 0.7

res1@cnFillDrawOrder = "Predraw"

;res1@gsnCenterString = date(d)
res1@gsnCenterString = hour(ind10(d))+"Z"
res1@gsnLeftString = "";"hPa"
res1@gsnRightString = "(J kg~S~-1~N~)"
res1@tiMainString = cape@long_name
plot = gsn_csm_contour_map(wks, cape(ind10(d),:,:), res1)
;plot = gsn_csm_contour_map(wks, mcape(:,:), res1)

res2 = True
res2@gsnDraw = False
res2@gsnFrame = False
res2@cnLevelSelectionMode = "ManualLevels"
res2@cnMinLevelValF = 1004.
res2@cnMaxLevelValF = 1012.
res2@cnLevelSpacingF = 2.
res2@cnLineThicknessF = 4.0
res2@cnLineLabelsOn = True
res2@gsnLeftString = ""
res2@gsnRightString = ""

;plotp = gsn_csm_contour(wks, msl(d,:,:), res2)
;overlay(plot,plotp)
draw(plot)
frame(wks)
;delete(mcape)
delete(plot)
;delete(plotp)
delete(wks)
end do
end