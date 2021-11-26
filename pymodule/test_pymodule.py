
import os
from os.path import dirname, abspath
libdir = abspath(dirname(__file__))
try:
    os.add_dll_directory(libdir)
except AttributeError:
    os.chdir(libdir)
import _smelib

print(_smelib.LibraryVersion())
print(_smelib.GetDataFiles())
print(_smelib.GetLibraryPath())

print(_smelib.SetLibraryPath("/bla/"))
print(_smelib.GetLibraryPath())
print(_smelib.InputWaveRange(1000, 2000))

# print(_smelib.InputWaveRange(1000, "2000"))

