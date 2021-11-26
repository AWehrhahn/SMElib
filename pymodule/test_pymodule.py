import _smelib

print(_smelib.LibraryVersion())
print(_smelib.GetDataFiles())
print(_smelib.GetLibraryPath())

print(_smelib.SetLibraryPath("/bla/"))
print(_smelib.GetLibraryPath())
print(_smelib.InputWaveRange(1000, 2000))
