
import os
from os.path import dirname, abspath, join
import numpy as np
from ctypes import cdll

libdir = abspath(dirname(__file__))
libfile = join(libdir, "libsme.so.5")

try:
    os.add_dll_directory(libdir)
except AttributeError:
    cdll.LoadLibrary(libfile)
    # os.chdir(libdir)
    # if "LD_LIBRARY_PATH" in os.environ:
    #     newpath = libdir + os.pathsep + os.environ["LD_LIBRARY_PATH"]
    # else:
    #     newpath = libdir
    # os.environ["LD_LIBRARY_PATH"] = newpath

import _smelib

# print(_smelib.LibraryVersion())
# print(_smelib.GetDataFiles())
# print(_smelib.GetLibraryPath())

# print(_smelib.SetLibraryPath("/bla/"))
# print(_smelib.GetLibraryPath())
# print(_smelib.InputWaveRange(1000, 2000))

species = np.asarray(["Na 1", "Li 1"], dtype="S8")
linelist = np.zeros((8, 2))
linelist[3, :] = 1
linelist[2, 0] = 100
linelist[2, 1] = 1000

_smelib.InputLineList(species, linelist)
# _smelib.UpdateLineList(["Na 1"], linelist[:, 0:1], [0])
# data = _smelib.OutputLineList()

# print(linelist[2:])
# print(data.T.shape)
# print(data.T)

opflag = np.zeros(20, dtype=np.int16)
depth = temp = xna = xne = rho = vt = height = np.zeros(10)
_smelib.InputModel(5000, 4.4, 5000, "RHOX", opflag, depth, temp, xna, xne, rho, vt, radius=-10)

bmat = np.ones((2, 10))
_smelib.InputDepartureCoefficients(0, bmat)

data = _smelib.GetDepartureCoefficients(0)
print(data)

# _smelib.InputLineList(species, "bla")


# print(_smelib.InputWaveRange(1000, "2000"))

