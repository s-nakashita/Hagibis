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
;  rank = "worst"
lon0 = 145
lon1 = 155
lat0 = 8
lat1 = 18
;orig = "ecmwf"
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
init = "0912"
out = "prtb_"+orig+"/ridge_"+orig+"_"+init+".txt"
a = addfile("/Volumes/dandelion/netcdf/tigge/2019/"+orig+"/glb_201910"+init+"_mean.nc","r")
psm = a->msl(0,{lat0:lat1},{lon0:lon1})
psm = psm*0.01
slp_ave = new(member+1,float)
do m=1,member
  mem = sprinti("%0.2d", m) ;"26"
  b = addfile("/Volumes/dandelion/netcdf/tigge/2019/"+orig+"/glb_201910"+init+"_"+mem+".nc","r")
  ps = b->msl(0,{lat0:lat1},{lon0:lon1})
  ps = ps*0.01
  ps = ps - psm
  slp_ave(m-1) = avg(ps)
  delete(b)
  delete(ps)
end do
b = addfile("../sens/ensvsa-"+etype+"-m1-"+orig+"-2019100912-2019101212_nlev"+sprinti("%d",nlev)+".nc", "r")
ps = b->pres_meansealev(0,{lat0:lat1},{lon0:lon1})
if(orig.ne."ncep") then
ps = (/-1*ps/)
end if
slp_ave(member) = avg(ps)
asciiwrite(out,slp_ave)
end
