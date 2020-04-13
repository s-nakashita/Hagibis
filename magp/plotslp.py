import sys
import subprocess
from pathlib import Path
from datetime import datetime, timedelta
import Magics.macro as mg

path = Path("../netcdf/gsm/gl")
initstr = "2019100712" if len(sys.agv) < 2 else sys.argv[1]
infile = "surface.nc"
t0 = datetime.strptime(initstr, "%Y%m%d%H")
dt = timedelta(hours = 6)
outdir = Path("../png") / initstr
outdir.mkdir(0o755, True, True)

page = mg.page(
    page_id_line = "off")

area = mg.mmap(
    subpage_map_projection = "cylindrical",
    subpage_lower_left_longitude  = 100,
    subpage_lower_left_latitude   =   0,
    subpage_upper_right_longitude = 180,
    subpage_upper_right_latitude  =  60,
)

coast = mg.mcoast(
    map_coastline_land_shade = "on",
    map_coastline_land_shade_colour = "cream",
    map_coastline_resolution = "medium",
    map_grid_line_style = "dash",
    map_grid_colour = "grey",
    map_label= "on",
    map_coastline_colour = "grey")

contour = mg.mcont(
    contour_highlight_colour = "black",
    contour_highlight_thickness = 4,
    contour_hilo = "off",
    contour_interval = 5,
    contour_label = "on",
    contour_label_frequency = 2,
    contour_label_height = 0.4,
    contour_level_selection_type = "interval",
    contour_line_colour = "black",
    contour_line_thickness = 2)

imax = 45 if t0.hour == 12 else 23
for i in range(imax):
    valid = t0 + i * dt
    print(valid.isoformat(" "))
    outname = "slp{}FT{:03d}".format(initstr, i * 6)
    ext = "png"
    outfile = outname + "." + ext

    output = mg.output(
        output_formats = [ext],
        output_name = outname,
        output_name_first_page_number = "off")

    msl = mg.mnetcdf(
        netcdf_filename = str(path / initstr / infile),
        netcdf_field_automatic_scaling = "off",
        netcdf_field_scaling_factor = 0.01,
        netcdf_time_dimension_setting = valid.isoformat(" "),
        netcdf_value_variable = "prmsl")

    title = mg.mtext(
        text_lines = ["init "+initstr, "valid "+valid.strftime("%Y%m%d%H")],
        text_font_size = 0.8,
        text_colour = "charcoal",)

    mg.plot(output, page, area, coast, msl, contour, title)
    subprocess.run(["gm", "convert", "-trim", outfile, str(outdir / outfile)])
    outfile.unlink()
