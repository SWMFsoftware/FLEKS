import math
import numpy as np

plot_unit_planet = {
    "time": "s",
    "t": "s",
    "mass": "amu",
    "rho": "amu/cm**3",
    "u": "km/s",
    "ux": "km/s",
    "uy": "km/s",
    "uz": "km/s",
    "p": "nPa",
    "pxx": "nPa",
    "pxy": "nPa",
    "pxz": "nPa",
    "pyy": "nPa",
    "pyz": "nPa",
    "pzz": "nPa",
    "b": "nT",
    "bx": "nT",
    "by": "nT",
    "bz": "nT",
    "e": "nT*km/s",
    "ex": "nT*km/s",
    "ey": "nT*km/s",
    "ez": "nT*km/s",
    "x": "Planet_Radius",
    "y": "Planet_Radius",
    "z": "Planet_Radius",
    "p_x": "Planet_Radius",
    "p_y": "Planet_Radius",
    "p_z": "Planet_Radius",
    "p_ux": "km/s",
    "p_uy": "km/s",
    "p_uz": "km/s",
    "p_w": "amu"
}

plot_unit_si = {
    "time": "s",
    "t": "s",
    "mass": "kg",
    "rho": "kg/m**3",
    "u": "m/s",
    "ux": "m/s",
    "uy": "m/s",
    "uz": "m/s",
    "p": "Pa",
    "pxx": "Pa",
    "pxy": "Pa",
    "pxz": "Pa",
    "pyy": "Pa",
    "pyz": "Pa",
    "pzz": "Pa",
    "b": "T",
    "bx": "T",
    "by": "T",
    "bz": "T",
    "e": "T*m/s",
    "ex": "T*m/s",
    "ey": "T*m/s",
    "ez": "T*m/s",
    "x": "m",
    "y": "m",
    "z": "m",
    "p_x": "m",
    "p_y": "m",
    "p_z": "m",
    "p_ux": "m/s",
    "p_uy": "m/s",
    "p_uz": "m/s",
    "p_w": "kg"
}


def get_unit(var, unit_type="planet"):
    if var[-1].isdigit():
        # Example: pxxs0 -> pxx
        var = var[0:-2]
    var = var.lower()

    if not (var in plot_unit_planet.keys()):
        return "dimensionless"

    unit_type = unit_type.lower()
    if unit_type == "planet" or unit_type == "planetary":
        return plot_unit_planet[var]
    elif unit_type == "si":
        return plot_unit_si[var]
    else:
        return "dimensionless"


def get_ticks(vmin, vmax):    
    dv = vmax - vmin
    if dv == 0: 
        return [vmin]
    norder = 10**(math.floor(math.log10(dv))-1)

    v0 = int(vmin/norder)
    v1 = int(vmax/norder) 

    dtick = int((v1-v0)/4) 
    dv = dtick*norder 

    tickMin = math.ceil(vmin/dv)*dv
    tickMax = math.floor(vmax/dv)*dv
    nticks = int((tickMax - tickMin)/dv) + 1
    return np.linspace(tickMin, tickMax, nticks)