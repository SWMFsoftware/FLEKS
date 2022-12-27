import yt
import numpy as np
import matplotlib.pyplot as plt
import math
from utilities import get_unit, get_ticks
import streamplot
from copy import deepcopy


def compare(d1, d2):
    header = ("var", "min|d1-d2|", "max|d1-d2|",
              "mean|d1-d2|", "mean|d1|", "mean|d2|")
    s = '{:8}   '+'{:18}'*5
    print(s.format(*header))
    for var in d1.vars:
        dif = abs(d1.data[var] - d2.data[var])
        l = (var, float(dif.min()), float(dif.max()), float(dif.mean()),
             float(d1.data[var].mean()), float(d2.data[var].mean()))
        s = '{:10}'+'{:+.6e},   '*5
        print(s.format(*l))


class dataContainer(object):
    def __init__(self, dataSets, x, y, z, xlabel, ylabel, zlabel, step=-1, time=-1, gencoord=False, filename=""):

        # Re-generate the coordinates to make sure they are equally spaced.
        # x = np.linspace(x[0], x[-1], len(x))
        # y = np.linspace(y[0], y[-1], len(y))
        # z = np.linspace(z[0], z[-1], len(z))

        self.data = dataSets
        self.x = x
        self.y = y
        self.z = z
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.zlabel = zlabel

        self.vars = [x for x in self.data.keys()]
        self.range = [[x[0], x[-1]], [y[0], y[-1]], [z[0], z[-1]]]
        self.dimensions = self.data[self.vars[0]].shape
        self.nstep = step
        self.time = time
        self.filename = filename
        self.gencoord = gencoord

    def __repr__(self):
        print("\n-----------------------------")
        print("Variables       :", self.vars)
        print("\nData range      :", self.range)
        print("\nData dimensions :", self.dimensions)
        print("\n")

        header = ("var", "min", "max", "mean")
        s = '{:8}   '+'{:18}'*3
        print(s.format(*header))
        for v in self.vars:
            d = self.data[v]
            s = '{:10}'+'{:+.6e},   '*3
            print(s.format(v, float(d.min()), float(d.max()), float(d.mean())))
        print("-----------------------------")
        return "\n"

    def add_bottom_line(self, f, verbose):
        if verbose > 0:
            s = "time = " + str(self.time)
        if verbose > 1:
            s += " "*5 + "nstep = "+str(self.nstep)
        if verbose > 2:
            s += " "*20+self.filename

        if verbose > 0:
            f.text(0.01, 0.01, s)

    def analyze_variable_string(self, var):
        r"""
        This method analyzes the plot string and return the plot variable and plot range. 

        Parameters
        ----------------------
        var: String
        Example: var = "{bb}<(-10)>(-9.8)"

        Return: a tuple contains the variable name, variable min and max. 
        Example: return "{bb}", -10, -9.8
        """
        vMin = None
        vMax = None

        varName = var
        if varName.find(">") > 0:
            varName = varName[:varName.find(">")]

        if varName.find("<") > 0:
            varName = varName[:varName.find("<")]

        if var.find(">") > 0:
            tmpVar = var[var.find(">")+2:]
            p1 = tmpVar.find(")")
            vMin = float(tmpVar[:p1])

        if var.find("<") > 0:
            tmpVar = var[var.find("<")+2:]
            p1 = tmpVar.find(")")
            vMax = float(tmpVar[:p1])

        return varName, vMin, vMax

    def evaluate_expression(self, expression, unit='planet'):
        r""" 
        This method calculates the variable expression and return the result, which is 
        a YTArray.  

        Parameters
        -------------------
        expression: String
        Example: expression = "np.log({rhos0}+{rhos1})"
        """

        if expression.find('{') < 0:
            return self.get_variable(expression, unit)

        exp1 = expression.replace('{', "self.get_variable('")

        exp2 = exp1.replace("}", "',unit)")
        exp2 = r"result="+exp2

        print(exp2)
        exec(exp2)

        return locals()['result']

    def add_variable(self, name, val):
        r"""
        This method adds a variable to the dataset for visualization purpose

        Parameters
        ------------------
        name: String
        The name of the variable to be added into self.data

        val: array-like structure
        Values of the variable stored in an array

        """

        if type(val) != yt.units.yt_array.YTArray:
            val = yt.YTArray(val, 'dimensionless')
        
        self.data[name] = val
        self.vars.append(name)

    def get_variable(self, var, unit='planet'):
        r"""
        This method calculate the value of the variable. 

        Parameters
        ------------------
        var: String
        Example: var = "pbeta"

        Return: A YTArray 
        """

        if var in self.data.keys():
            varUnit = get_unit(var, unit)
            ytarr = self.data[var]
        else:
            var = var.lower()
            expression = None
            if var == "b":
                expression = "np.sqrt({Bx}**2+{By}**2+{Bz}**2)"
                varUnit = get_unit('b', unit)
            elif var == 'bb':
                expression = "{Bx}**2+{By}**2+{Bz}**2"
                varUnit = get_unit('b', unit)+"**2"
            elif var[0:2] == 'ps':
                ss = var[2:3]
                expression = "({pxxs"+ss+"}+"+"{pyys"+ss+"}+"+"{pzzs"+ss+"})/3"
                varUnit = get_unit('p', unit)
            elif var == 'pb':
                coef = 0.5/(yt.units.mu_0.value)
                ytarr = coef*self.get_variable('bb', 'si')
                ytarr = yt.YTArray(ytarr.value, 'Pa')
                varUnit = get_unit('p', unit)
            elif var == "pbeta":
                ytarr = (self.get_variable('ps0', unit) +
                         self.get_variable('ps1', unit)) / self.get_variable('pb', unit)
                varUnit = 'dimensionless'
            elif var == "calfven":
                ytarr = self.get_variable(
                    'b', 'si')/np.sqrt(yt.units.mu_0.value * self.get_variable('rhos1', 'si'))
                ytarr = yt.YTArray(ytarr.value, 'm/s')
                varUnit = get_unit('u', unit)

            if expression != None:
                ytarr = self.evaluate_expression(expression, unit)
                if type(ytarr) != yt.units.yt_array.YTArray:
                    varUnit = 'dimensionless'
                    ytarr = yt.YTArray(ytarr, varUnit)

        return ytarr if str(ytarr.units) == 'dimensionless' else ytarr.in_units(varUnit)


class dataContainer3D(dataContainer):
    r"""
    A class handles 3D box data sets.     
    """

    def __init__(self, dataSets, x, y, z, xlabel='X', ylabel='Y', zlabel='Z', *args, **kwargs):
        r"""
        Parameters
        ---------------------
        dataSets: dictonary 
        The key is the variable name, and the dictionary value is usually a YTArray. 

        x/y/z: A 1D YTArray 
        """
        super(dataContainer3D, self).__init__(
            dataSets, x, y, z, xlabel, ylabel, zlabel, *args, **kwargs)

    def get_slice(self, norm, cut_loc):
        r"""
        Get a 2D slice from the 3D box data. 

        Parameters
        ------------------
        norm: String
        The normal direction of the slice. Potential value: 'x', 'y' or 'z'

        cur_loc: Float 
        The position of slicing. 

        Return: A dataContainer2D object 
        """

        axDir = {'X': 0, 'Y': 1, 'Z': 2}
        idir = axDir[norm.upper()]

        axes = [self.x, self.y, self.z]

        axSlice = axes[idir]

        for iSlice in range(axSlice.size):
            if axSlice[iSlice] > cut_loc:
                break

        dataSets = {}
        for varname, val in self.data.items():
            if idir == 0:
                arr = val[iSlice, :, :]
            elif idir == 1:
                arr = val[:, iSlice, :]
            elif idir == 2:
                arr = val[:, :, iSlice]
            dataSets[varname] = np.squeeze(arr)

        axLabes = {0: ('Y', 'Z'), 1: ('X', 'Z'), 2: ('X', 'Y')}
        ax = {0: (1, 2), 1: (0, 2), 2: (0, 1)}

        return dataContainer2D(
            dataSets, axes[ax[idir][0]], axes[ax[idir][1]], axLabes[idir][0], axLabes[idir][1], norm, cut_loc)


class dataContainer2D(dataContainer):
    r""" 
    A class handles 2D Cartesian data. 
    """

    def __init__(self, dataSets, x, y, xlabel, ylabel, cut_norm=None, cut_loc=None, *args, **kwargs):
        r"""
        Parameters
        ---------------------
        dataSets: dictonary 
        The key is the variable name, and the dictionary value is usually a YTArray. 

        x/y: A 1D YTArray 

        xlabel/ylabel: String 

        cut_norm: String
        'x', 'y' or 'z'

        cut_loc: Float 
        cut_norm and cut_loc are used to record the position of slice if this 2D 
        data set is obtained from a 3D box. 
        """

        zlabel = None
        z = (0, 0)
        super(dataContainer2D, self).__init__(
            dataSets, x, y, z, xlabel, ylabel, zlabel, *args, **kwargs)

        self.cut_norm = cut_norm
        self.cut_loc = cut_loc

    def __sub__(self, other):
        dif = deepcopy(self)
        for var in dif.vars:
            dif.data[var] -= other.data[var]
        return dif

    def contour(self, vars, xlim=None, ylim=None, unit="planet", nlevels=200,
                cmap="rainbow", figsize=(12, 8), pcolor=False, log=False,
                addgrid=False, bottomline=10, *args, **kwargs):
        r""" 
        Contour plots. 

        Parameters
        ----------------------
        vars: String
        Ploting variables and ploting range. 
        Example: vars = "Bx<(50)>(-50) By (np.log(2*{rhos0}))>(-5)"

        xlim/ylim: A list/tuple contains the x- y-axis range

        unit: String
        'planet' or 'si' 

        nlevels: Integer 
        Number of the countour plot color levels

        cmap: String 
        Color map type 

        figsize: A tuple         

        log: Bool 
        Using log plot or not. 

        Examples
        ----------------
        >>> f,axes = dc.contour("Bx<(50)>(-50) By (np.log(2*{rhos0}))>(-5)",xlim=[-40,-5])
        """

        if type(vars) == str:
            vars = vars.split()

        nvar = len(vars)
        nRow = int(round(np.sqrt(nvar)))
        nCol = math.ceil(nvar/nRow)

        varNames = []
        varMin = []
        varMax = []
        for var in vars:
            vname, vmin, vmax = self.analyze_variable_string(var)
            varNames.append(vname)
            varMin.append(vmin)
            varMax.append(vmax)

        f, axes = plt.subplots(nRow, nCol, figsize=figsize)
        axes = np.array(axes)  # in case nRow = nCol = 1

        axes = axes.reshape(-1)

        for isub, ax in zip(range(nvar), axes):
            ytVar = self.evaluate_expression(varNames[isub], unit)
            v = ytVar
            varUnit = 'dimensionless'
            if type(ytVar) == yt.units.yt_array.YTArray:
                v = ytVar.value
                varUnit = str(ytVar.units)

            vmin = v.min() if varMin[isub] == None else varMin[isub]
            vmax = v.max() if varMax[isub] == None else varMax[isub]

            logplot = log and vmin > 0
            if logplot:
                v = np.log10(v)

            levels = np.linspace(vmin, vmax, nlevels)
            if self.gencoord:
                if pcolor or abs(vmin-vmax) < 1e-20*abs(vmax):
                    cs = ax.tripcolor(self.x.value, self.y.value, v.T,
                                      cmap=cmap, *args, **kwargs)
                else:
                    cs = ax.tricontourf(self.x.value, self.y.value, v.T, levels=levels,
                                        cmap=cmap, extend="both", *args, **kwargs)
            else:
                if pcolor or abs(vmin-vmax) < 1e-20*abs(vmax):
                    cs = ax.pcolormesh(self.x.value, self.y.value, v.T,
                                       cmap=cmap, *args, **kwargs)
                else:
                    cs = ax.contourf(self.x.value, self.y.value, v.T, levels=levels,
                                     cmap=cmap, extend="both", *args, **kwargs)
            if addgrid:
                if self.gencoord:
                    gx = self.x.value
                    gy = self.y.value
                else:
                    gg = np.meshgrid(self.x.value, self.y.value)
                    gx = np.reshape(gg[0], -1)
                    gy = np.reshape(gg[1], -1)

                ax.plot(gx, gy, 'x')

            ticks = get_ticks(vmin, vmax)
            cb = f.colorbar(cs, ax=ax, ticks=ticks)
            cb.formatter.set_powerlimits((0, 0))

            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.set_xlabel(self.xlabel)
            ax.set_ylabel(self.ylabel)
            title = varNames[isub]
            if varUnit != 'dimensionless':
                title = title + " ["+varUnit+"]"
            if logplot:
                title = "$log_{10}$("+title+")"
            ax.set_title(title)

        # Delete axes for empty subplots
        for ax in axes[nvar:nRow*nCol]:
            f.delaxes(ax)

        plt.tight_layout()
        if self.cut_norm != None and self.cut_loc != None:
            print("Plots at "+self.cut_norm+' = ', self.cut_loc)

        self.add_bottom_line(f, bottomline)
        return f, axes.reshape(nRow, nCol)

    def add_contour(self, ax, var, unit='planet', *args, **kwargs):
        r""" 
        Adding contour lines to an axis. 

        Parameters
        --------------------
        ax: matplotlib axis
        The aixs to plot contour lines. 

        var: String 

        Examples
        ----------------------
        >>> f,axes = dc.contour("Bx<(50)>(-50) By (np.log(2*{rhos0}))>(-5)",xlim=[-40,-5])
        >>> dc.add_contour(axes[0,0], "rhos1>1")        
        """

        vname, vmin, vmax = self.analyze_variable_string(var)

        ytVar = self.evaluate_expression(vname, unit)
        v = ytVar
        if type(ytVar) == yt.units.yt_array.YTArray:
            v = ytVar.value

        vmin = v.min() if vmin == None else vmin
        vmax = v.max() if vmin == None else vmax
        v = np.clip(v, vmin, vmax)

        ax.contour(self.x, self.y, v.T, *args, **kwargs)

    def add_stream(self, ax, var1, var2, density=1, *args, **kwargs):
        r""" 
        Adding streamlines to an axis. 

        Parameters
        --------------------
        ax: matplotlib axis
        The aixs to add streamlines 

        var1/var2: String 

        density: Int
        It controls the number of streamlines. 

        Examples
        ----------------------
        >>> f,axes = dc.contour("Bx<(50)>(-50) By (np.log(2*{rhos0}))>(-5)",xlim=[-40,-5])
        >>> dc.add_stream(axes[1,0], "Bx", "Bz", density=2)
        """

        v1 = self.evaluate_expression(var1).value
        v2 = self.evaluate_expression(var2).value
        if type(v1) == yt.units.yt_array.YTArray:
            v1 = v1.value
        if type(v2) == yt.units.yt_array.YTArray:
            v2 = v2.value
        streamplot.streamplot(ax, self.x.value, self.y.value,
                              v1.T, v2.T, density=density, *args, **kwargs)


class dataContainer1D(dataContainer):
    r""" 
    A class handles 1D Cartesian data. 
    """

    def __init__(self, dataSets, x, xlabel, *args, **kwargs):
        r"""
        Parameters
        ---------------------
        dataSets: dictonary 
        The key is the variable name, and the dictionary value is usually a YTArray. 

        x: A 1D YTArray 

        xlabel: String 
        """

        ylabel = None
        y = (0, 0)

        zlabel = None
        z = (0, 0)

        super(dataContainer1D, self).__init__(
            dataSets, x, y, z, xlabel, ylabel, zlabel, *args, **kwargs)

    def plot(self, vars, xlim=None, ylim=None, unit="planet",
             figsize=(12, 8), log=False, bottomline=10, *args, **kwargs):
        """
        Examples
        ----------------
        >>> f,axes=dc.plot("absdivb bx",xlim=[-5,5],color='k',linestyle='solid',marker='x')
        """

        if type(vars) == str:
            vars = vars.split()

        nvar = len(vars)
        nRow = int(round(np.sqrt(nvar)))
        nCol = math.ceil(nvar/nRow)

        varNames = []
        varMin = []
        varMax = []
        for var in vars:
            vname, vmin, vmax = self.analyze_variable_string(var)
            varNames.append(vname)
            varMin.append(vmin)
            varMax.append(vmax)

        f, axes = plt.subplots(nRow, nCol, figsize=figsize)
        axes = np.array(axes)  # in case nRow = nCol = 1

        axes = axes.reshape(-1)

        for isub, ax in zip(range(nvar), axes):
            ytVar = self.evaluate_expression(varNames[isub], unit)
            v = ytVar
            varUnit = 'dimensionless'
            if type(ytVar) == yt.units.yt_array.YTArray:
                v = ytVar.value
                varUnit = str(ytVar.units)

            vmin = v.min() if varMin[isub] == None else varMin[isub]
            vmax = v.max() if varMax[isub] == None else varMax[isub]
            # v = np.clip(v, vmin, vmax)
            ax.plot(self.x.value, v, label=varNames[isub], *args, **kwargs)

            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.set_xlabel(self.xlabel)
            ax.legend()

        self.add_bottom_line(f, bottomline)
        return f, axes
