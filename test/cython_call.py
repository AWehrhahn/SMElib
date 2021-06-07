import numpy as np
import cywrapper

version = cywrapper.SMELibraryVersion()
print(version)

files = cywrapper.GetDataFiles()
print(files)

path = cywrapper.GetLibraryPath()
print(path)

resp = cywrapper.SetLibraryPath(path)
print(resp)

resp = cywrapper.InputWaveRange(2000, 4000)
print(resp)

resp = cywrapper.SetVWscale(1)
print(resp)

resp = cywrapper.SetH2broad()
print(resp)

resp = cywrapper.ClearH2broad()
print(resp)

species = np.array(["Fe 1", "S 1"])
atomic = np.zeros((2, 8))
resp = cywrapper.InputLineList(species, atomic)
print(resp)

atomic = np.zeros((2, 8))
resp = cywrapper.OutputLineList()
print(resp)

species = np.array(["S 1"])
atomic = np.zeros((1, 8))
index = np.array([1], dtype=np.short)
resp = cywrapper.UpdateLineList(species, atomic, index)
print(resp)

opflag = np.zeros(15, dtype="short")
depth = np.zeros(2)
temp = np.zeros(2)
xne = np.zeros(2)
xna = np.zeros(2)
rho = np.zeros(2)
vt = np.zeros(2)
radius = 1000
height = np.zeros(2)

resp = cywrapper.InputModel(
    5000,
    4.0,
    "SPH",
    5000,
    opflag,
    depth,
    temp,
    xne,
    xna,
    rho,
    vt,
    radius, 
    height,
    )
print(resp)

bmat = np.zeros((2, 2))
resp = cywrapper.InputDepartureCoefficients(bmat, 1)
print(resp)
