# distutils: language = c++
# cython: infer_types=True

import numpy as np
cimport numpy as np
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free

np_int = np.int
np_double = np.double
np_str = np.str_
np_short = np.short

ctypedef np.int_t np_int_t
ctypedef np.double_t np_double_t
ctypedef str np_str_t
ctypedef short np_short_t
ctypedef float np_float_t


cdef extern from "../src/sme/sme_python_bridge.h":
    cdef struct s_IDL_STRING:
        int slen
        short stype
        char * s

    ctypedef s_IDL_STRING IDL_STRING

    const char * Python_SMELibraryVersion()
    const char * Python_GetDataFiles()
    const char * Python_GetLibraryPath()
    const char * Python_SetLibraryPath(IDL_STRING * path)
    const char * Python_InputWaveRange(double wmin, double wmax)
    const char * Python_SetVWscale(double gamma6)
    const char * Python_SetH2broad()
    const char * Python_ClearH2broad()
    const char * Python_InputLineList(int nlines, IDL_STRING * species, double * atomic)
    const char * Python_OutputLineList(int nlines, double * atomic)
    const char * Python_UpdateLineList(short nlines, IDL_STRING * species, double * atomic, short * index)
    const char * Python_InputModel(short ndepth, double teff, double grav, double wlstd, IDL_STRING * motype, short * opflag, double * depth, double * temp, double * xne, double * xna, double * rho, double * vt, double radius, double * height)
    const char * Python_InputDepartureCoefficients(double * bmat, int lineindex)
    const char * Python_GetDepartureCoefficients(double * bmat, int nrhox, int line)
    const char * Python_GetNLTEflags(short * nlte_flags, int nlines)
    const char * Python_ResetDepartureCoefficients()
    const char * Python_InputAbund(double * abund)
    const char * Python_Opacity()
    const char * Python_GetOpacity(short ifop, short length, double * result, IDL_STRING * species, IDL_STRING * key)
    const char * Python_Ionization(short ion)
    const char * Python_GetDensity(short length, double * result)
    const char * Python_GetNatom(short length, double * result)
    const char * Python_GetNelec(short length, double * result)
    const char * Python_Transf(short nmu, double * mu, double * cint_seg, double * cintr_seg, int nwmax, int nw, double * wint_seg, double * sint_seg, double accrt, double accwi, short keep_lineop, short long_continuum)
    const char * Python_CentralDepth(int nmu, double * mu, int nwsize, float * table, double accrt)
    const char * Python_GetLineOpacity(double wave, short nmu, double * lop, double * cop, double * scr, double * tsf, double * csf)
    const char * Python_GetLineRange(double * linerange, int nlines)

cdef char * _chars(s):
    if isinstance(s, unicode):
        # encode to the specific encoding used inside of the module
        s = (<unicode>s).encode('utf8')
    return s

cdef IDL_STRING _idl(s):
    cdef IDL_STRING idl_str
    idl_str.slen = len(s)
    idl_str.stype = 1
    idl_str.s = _chars(s)
    return idl_str

cdef IDL_STRING * to_strarr(np.ndarray arr):
    cdef Py_ssize_t i
    cdef Py_ssize_t size
    cdef IDL_STRING * strarr
    
    size = arr.shape[0]
    strarr = <IDL_STRING*> PyMem_Malloc(size * sizeof(IDL_STRING *))
    for i in range(size):
        strarr[i] = _idl(arr[i])

    return strarr

cdef double * to_arr_double_2d(np_double_t[:,::1] arr):
    return &arr[0, 0]

cdef double * to_arr_double(np_double_t[::1] arr):
    return &arr[0]

cdef short * to_arr_short(np_short_t[::1] arr):
    return &arr[0]

cdef float * to_arr_float(np_float_t[::1] arr):
    return &arr[0]

def SMELibraryVersion() -> str:
    cdef const char * byte
    byte = Python_SMELibraryVersion()
    return byte.decode("utf8")

def GetDataFiles() -> list:
    cdef const char * byte
    byte = Python_GetDataFiles()
    files = byte.decode("utf8")
    files = files.split(";")
    return files

def GetLibraryPath() -> str:
    cdef const char * byte
    byte = Python_GetLibraryPath()
    return byte.decode("utf8")

def SetLibraryPath(path:str) -> str:
    cdef IDL_STRING path_idl
    cdef const char * byte
    path_idl = _idl(path)
    byte = Python_SetLibraryPath(&path_idl)
    return byte.decode("utf8")

def InputWaveRange(double wmin, double wmax) -> str:
    cdef const char * byte
    byte = Python_InputWaveRange(wmin, wmax)
    return byte.decode("utf8")

def SetVWscale(double gamma6) -> str:
    cdef const char * byte
    byte = Python_SetVWscale(gamma6)
    return byte.decode("utf8")

def SetH2broad(flag=True) -> str:
    cdef const char * byte
    if flag:
        byte = Python_SetH2broad()
    else:
        byte = Python_ClearH2broad()
    return byte.decode("utf8")

def ClearH2broad() -> str:
    return SetH2broad(flag=False)

def InputLineList(np.ndarray species, np_double_t[:, ::1] atomic) -> str:
    cdef int nlines
    cdef double * atomic_data
    cdef IDL_STRING * species_data
    cdef const char * byte

    # TODO: Is there a way to avoid creating a new species array everytime?
    nlines = species.shape[0]
    species_data = to_strarr(species)
    atomic_data = to_arr_double_2d(atomic)
    byte = Python_InputLineList(nlines, species_data, atomic_data)
    PyMem_Free(species_data)
    return byte.decode("utf8")

def OutputLineList(np_double_t[:, ::1] atomic) -> str:
    cdef int nlines
    cdef const char * byte

    nlines = atomic.shape[0]
    byte = Python_OutputLineList(nlines, to_arr_double_2d(atomic))
    return byte.decode("utf8")

def UpdateLineList(np.ndarray species, np_double_t[:, ::1] atomic, np_short_t[::1] index) -> str:
    cdef short nlines
    cdef double * atomic_data
    cdef short * index_data
    cdef IDL_STRING * species_data
    cdef const char * byte

    # TODO: Is there a way to avoid creating a new species array everytime?
    nlines = species.shape[0]
    species_data = to_strarr(species)
    atomic_data = to_arr_double_2d(atomic)
    index_data = &index[0]
    byte = Python_UpdateLineList(nlines, species_data, atomic_data, index_data)
    PyMem_Free(species_data)
    return byte.decode("utf8")

def InputModel(
    np_double_t teff,
    np_double_t grav,
    motype:str,
    np_double_t wlstd,
    np_short_t[::1] opflag,
    np_double_t[::1] depth,
    np_double_t[::1] temp,
    np_double_t[::1] xne,
    np_double_t[::1] xna,
    np_double_t[::1] rho,
    np_double_t[::1] vt,
    np_double_t radius, 
    np_double_t[::1] height,
    )->str:
    cdef short ndepth;
    cdef IDL_STRING motype_str;
    cdef const char * byte

    motype_str = _idl(motype)
    ndepth = depth.shape[0]
    byte = Python_InputModel(
        ndepth, 
        teff, 
        grav, 
        wlstd, 
        &motype_str, 
        to_arr_short(opflag),
        to_arr_double(depth),
        to_arr_double(temp),
        to_arr_double(xne),
        to_arr_double(xna),
        to_arr_double(rho), 
        to_arr_double(vt), 
        radius, 
        to_arr_double(height),
        )
    return byte.decode("utf8")

def InputDepartureCoefficients(np_double_t[:, ::1] bmat, int lineindex) -> str:
    cdef const char * byte
    byte = Python_InputDepartureCoefficients(to_arr_double_2d(bmat), lineindex)
    return byte.decode("utf8")

def GetDepartureCoefficients(np_double_t[:, ::1] bmat, int ndepth, int line) -> str:
    cdef const char * byte
    byte = Python_GetDepartureCoefficients(to_arr_double_2d(bmat), ndepth, line)
    return byte.decode("utf8")

def GetNLTEflags(np_short_t[::1] nlte_flags) -> str:
    cdef const char * byte
    cdef int nlines
    nlines = nlte_flags.shape[0]
    byte = Python_GetNLTEflags(to_arr_short(nlte_flags), nlines)
    return byte.decode("utf8")

def ResetDepartureCoefficients() -> str:
    cdef const char * byte
    byte = Python_ResetDepartureCoefficients()
    return byte.decode("utf8")

def InputAbund(np_double_t[::1] abund) -> str:
    cdef const char * byte
    byte = Python_InputAbund(to_arr_double(abund))
    return byte.decode()

def Opacity() -> str:
    cdef const char * byte
    byte = Python_Opacity()
    return byte.decode()

def GetOpacity(short ifop, short length, np_double_t[::1] result, species:str = "", key:str = "") -> str:
    cdef const char * byte
    cdef IDL_STRING species_idl
    cdef IDL_STRING key_idl
    species_idl = _idl(species)
    key_idl = _idl(key)
    byte = Python_GetOpacity(ifop, length, to_arr_double(result), &species_idl, &key_idl)
    return byte.decode()

def Ionization(short ion = 0) -> str:
    cdef const char * byte
    byte = Python_Ionization(ion)
    return byte.decode()

def GetDensity(np_double_t[::1] result) -> str:
    cdef const char * byte
    cdef short length
    length = result.shape[0]
    byte = Python_GetDensity(length, to_arr_double(result))
    return byte.decode("utf8")

def GetNatom(np_double_t[::1] result) -> str:
    cdef const char * byte
    cdef short length
    length = result.shape[0]
    byte = Python_GetNatom(length, to_arr_double(result))
    return byte.decode("utf8")

def GetNelec(np_double_t[::1] result) -> str:
    cdef const char * byte
    cdef short length
    length = result.shape[0]
    byte = Python_GetNelec(length, to_arr_double(result))
    return byte.decode("utf8")

def Transf(
    np_double_t[::1] mu, 
    np_double_t[:, ::1] cint_seg, 
    np_double_t[::1] cintr_seg, 
    int nwmax,
    int nw, 
    np_double_t[::1] wint_seg, 
    np_double_t[:, ::1] sint_seg, 
    double accrt = 0.01, 
    double accwi = 0.03, 
    short keep_lineop = 0, 
    short long_continuum = 1
    ) -> str:
    cdef const char * byte
    cdef short nmu
    cdef short keep_lineop_c
    cdef short long_continuum_c
    keep_lineop_c = 1 if keep_lineop else 0
    long_continuum_c = 1 if long_continuum else 0
    nmu = mu.shape[0]
    byte = Python_Transf(
        nmu, 
        to_arr_double(mu), 
        to_arr_double_2d(cint_seg), 
        to_arr_double(cintr_seg), 
        nwmax, 
        nw,
        to_arr_double(wint_seg), 
        to_arr_double_2d(sint_seg), 
        accrt, 
        accwi, 
        keep_lineop_c, 
        long_continuum_c,
        )
    return byte.decode("utf8")

def CentralDepth(np_double_t[::1] mu, np_float_t[::1] table, double accrt) -> str:
    cdef const char * byte
    cdef int nmu

    nmu = mu.shape[0]
    nwsize = table.shape[0]

    byte = Python_CentralDepth(nmu, to_arr_double(mu), nwsize, to_arr_float(table), accrt)
    return byte.decode("utf8")

def GetLineOpacity(double wave, short nmu, np_double_t[::1] lop, np_double_t[::1] cop, np_double_t[::1] scr, np_double_t[::1] tsf, np_double_t[::1] csf) -> str:
    cdef const char * byte
    byte = Python_GetLineOpacity(wave, nmu, to_arr_double(lop), to_arr_double(cop), to_arr_double(scr), to_arr_double(tsf), to_arr_double(csf))
    return byte.decode("utf8")

def GetLineRange(np_double_t[:,::1] linerange) -> str:
    cdef const char * byte
    cdef int nlines
    nlines = linerange.shape[0]
    byte = Python_GetLineRange(to_arr_double_2d(linerange), nlines)
    return byte.decode("utf8")
