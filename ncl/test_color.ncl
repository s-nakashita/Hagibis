load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
begin 
  wks = gsn_open_wks("x11","test")
  gsn_define_colormap(wks,"mycmap")
  gsn_draw_colormap(wks)
end