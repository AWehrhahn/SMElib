#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>

// Header of the SME library
#include "sme_synth_faster.h"

// Everything else is considered an error
const char OK_response = '\0';


static char module_docstring[] =
    "This module provides a Python interface to the SME libary";

static char smelib_LibraryVersion_docstring[] = "Return SME library version";
static PyObject * smelib_LibraryVersion(PyObject *self, PyObject *args)
{
    const char * version;
    PyObject *ret;
    
    version = SMELibraryVersion(0, NULL);
    ret = Py_BuildValue("s", version);

    return ret;
}

static char smelib_GetDataFiles_docstring[] = "Return the required data files";
static PyObject * smelib_GetDataFiles(PyObject *self, PyObject *args)
{
    const char * datafiles;
    PyObject *ret;
    
    datafiles = GetDataFiles(0, NULL);
    ret = Py_BuildValue("s", datafiles);

    return ret;   
}

static char smelib_GetLibraryPath_docstring[] = "Return the current data file directory";
static PyObject * smelib_GetLibraryPath(PyObject* self, PyObject *args)
{
    const char * path;
    PyObject *ret;
    
    path = GetLibraryPath(0, NULL);
    ret = Py_BuildValue("s", path);

    return ret;  
}

static char smelib_SetLibraryPath_docstring[] = "Set the data file directory";
static PyObject * smelib_SetLibraryPath(PyObject * self, PyObject *args)
{
    const char * path;
    IDL_STRING idl_path;
    const char * result = NULL;
    void * args_c[1];

    if (!PyArg_ParseTuple(args, "s", &path))
        return NULL;
    

    // Create IDL String
    idl_path.slen = strlen(path);
    idl_path.s = (char *) malloc(sizeof(char) * idl_path.slen);
    if (idl_path.s == NULL) {
        PyErr_SetString(PyExc_MemoryError, "Could not assign memory");
        return NULL;
    }
    strcpy(idl_path.s, path);

    // Pass to SMELIB
    args_c[0] = &idl_path;
    result = SetLibraryPath(1, args_c);

    // Clean pointers
    free(idl_path.s);

    // Check for errors
    if (result != NULL && result[0] != OK_response)
    {
        PyErr_SetString(PyExc_RuntimeError, result);
        return NULL;
    }

    // Return the Ok
    Py_RETURN_NONE;
}



static char smelib_InputWaveRange_docstring[] = "Read in Wavelength range";
static PyObject * smelib_InputWaveRange(PyObject * self, PyObject *args)
{
    int n = 2;
    const char * result = NULL;
    void * args_c[n];
    double wmin, wmax;

    if (!PyArg_ParseTuple(args, "dd", &wmin, &wmax))
        return NULL;

    args_c[0] = &wmin;
    args_c[1] = &wmax;
    result = InputWaveRange(n, args_c);

    if (result != NULL && result[0] != OK_response)
    {
        PyErr_SetString(PyExc_RuntimeError, result);
        return NULL;
    }

    Py_RETURN_NONE;
}


static char smelib_SetVWscale_docstring[] = "Set van der Waals scaling factor";
static PyObject * smelib_SetVWscale(PyObject * self, PyObject *args)
{
    int n = 1;
    const char * result = NULL;
    void * args_c[n];
    double vwscale;

    if (!PyArg_ParseTuple(args, "d", &vwscale))
        return NULL;

    args_c[0] = &vwscale;
    result = SetVWscale(n, args_c);

    if (result != NULL && result[0] != OK_response)
    {
        PyErr_SetString(PyExc_RuntimeError, result);
        return NULL;
    }

    Py_RETURN_NONE;
}

static char smelib_SetH2broad_docstring[] = "Set flag for H2 molecule";
static PyObject * smelib_SetH2broad(PyObject * self, PyObject *args)
{
    int n = 0;
    const char * result = NULL;
    void ** args_c = NULL;

    result = SetH2broad(n, args_c);

    if (result != NULL && result[0] != OK_response)
    {
        PyErr_SetString(PyExc_RuntimeError, result);
        return NULL;
    }

    Py_RETURN_NONE;
}

static char smelib_ClearH2broad_docstring[] = "Clear flag for H2 molecule";
static PyObject * smelib_ClearH2broad(PyObject * self, PyObject *args)
{
    int n = 0;
    const char * result = NULL;
    void ** args_c = NULL;

    result = ClearH2broad(n, args_c);

    if (result != NULL && result[0] != OK_response)
    {
        PyErr_SetString(PyExc_RuntimeError, result);
        return NULL;
    }

    Py_RETURN_NONE;
}

static char smelib_InputLineList_docstring[] = "Read in line list";
static PyObject * smelib_InputLineList(PyObject * self, PyObject *args)
{
    const int n = 3;
    const char * result = NULL;
    void * args_c[n];
    int nlines = 0, nchar = 0;
    PyObject * species_obj = NULL, * linelist_obj = NULL;
    PyArrayObject * species_array = NULL, *linelist_array = NULL;
    IDL_STRING * species = NULL;
    char * species_data = NULL;
    double * linelist = NULL;
    PyArray_Descr * dtype = NULL;


    if (!PyArg_ParseTuple(args, "OO", &species_obj, &linelist_obj))
        return NULL;

    // Convert to Numpy arrays
    species_array = (PyArrayObject *) PyArray_FROM_OTF(species_obj, NPY_STRING, NPY_ARRAY_IN_ARRAY);
    linelist_array = (PyArrayObject *) PyArray_FROM_OTF(linelist_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);

    
    // Check if that worked
    if (species_array == NULL || linelist_array == NULL){
        goto fail;
    }

    // Check dimensions
    if (PyArray_NDIM(species_array) != 1){
        PyErr_SetString(PyExc_ValueError, "Expected species array of ndim == 1");
        goto fail;
    }
    if (PyArray_NDIM(linelist_array) != 2){
        PyErr_SetString(PyExc_ValueError, "Expected linelist of ndim == 2");
        goto fail;
    }
    // Get sizes
    nlines = PyArray_DIM(species_array, 0);
    dtype = PyArray_DESCR(species_array);
    nchar = dtype->elsize;

    if (PyArray_DIM(linelist_array, 0) != 8){
        PyErr_SetString(PyExc_ValueError, "Expected linelist to have 8 values");
        goto fail;
    }
    if (PyArray_DIM(linelist_array, 1) != nlines){
        PyErr_SetString(PyExc_ValueError, "Expected both arrays to have the same length");
        goto fail;
    }

    // Get the pointers to pass to SMELIB
    linelist = (double*) PyArray_DATA(linelist_array);
    species_data = (char*) PyArray_DATA(species_array);
    species = (IDL_STRING *) malloc(nlines * sizeof(IDL_STRING));

    for (int i = 0; i < nlines; i++)
    {
        species[i].slen = nchar;
        species[i].s = &species_data[nchar*i];
        species[i].stype = 0;
    }
  
    // printf("0: %s\n", species[0].s);
    args_c[0] = &nlines;
    args_c[1] = species;
    args_c[2] = linelist;
    result = InputLineList(n, args_c);


    free(species);
    Py_DECREF(species_array);
    Py_DECREF(linelist_array);

    if (result != NULL && result[0] != OK_response)
    {
        PyErr_SetString(PyExc_RuntimeError, result);
        return NULL;
    }

    Py_RETURN_NONE;

fail:
    Py_XDECREF(species_array);
    Py_XDECREF(linelist_array);
    return NULL;
}


static char smelib_OutputLineList_docstring[] = "Return line list";
static PyObject * smelib_OutputLineList(PyObject * self, PyObject *args)
{
    const int n = 2;
    const char * result = NULL;
    void * args_c[n];
    int nlines;
    double * linelist = NULL;
    PyArrayObject * linelist_array;

    nlines = GetNLINES();

    npy_intp dims[2] = {nlines, 6};
    linelist_array = (PyArrayObject*) PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    linelist = (double*) PyArray_DATA(linelist_array);

    args_c[0] = &nlines;
    args_c[1] = linelist;
    result = OutputLineList(n, args_c);

    if (result != NULL && result[0] != OK_response)
    {
        Py_DecRef((PyObject*)linelist_array);
        PyErr_SetString(PyExc_RuntimeError, result);
        return NULL;
    }

    return (PyObject*) linelist_array;
}

static char smelib_UpdateLineList_docstring[] = "Change line list parameters";
static PyObject * smelib_UpdateLineList(PyObject * self, PyObject *args)
{
    const int n = 4;
    const char * result = NULL;
    void * args_c[n];
    int nlines, nchar;
    double * linelist = NULL;
    IDL_STRING * species = NULL;
    short * index = NULL;
    PyObject * linelist_obj = NULL, * species_obj = NULL, * index_obj = NULL;
    PyArrayObject * linelist_array = NULL, * species_array = NULL, * index_array = NULL;
    char * species_data = NULL;
    PyArray_Descr * dtype = NULL;
    
    
    if (!PyArg_ParseTuple(args, "OOO", &species_obj, &linelist_obj, &index_obj))
        return NULL;

    // Convert to Numpy arrays
    species_array = (PyArrayObject *) PyArray_FROM_OTF(species_obj, NPY_STRING, NPY_ARRAY_IN_ARRAY);
    linelist_array = (PyArrayObject *) PyArray_FROM_OTF(linelist_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    index_array = (PyArrayObject *) PyArray_FROM_OTF(index_obj, NPY_SHORT, NPY_ARRAY_IN_ARRAY);

    if (species_array == NULL || linelist_array == NULL || index_array == NULL){
        goto fail;
    }

     // Check dimensions
    if (PyArray_NDIM(species_array) != 1){
        PyErr_SetString(PyExc_ValueError, "Expected species array of ndim == 1");
        goto fail;
    }
    if (PyArray_NDIM(linelist_array) != 2){
        PyErr_SetString(PyExc_ValueError, "Expected linelist of ndim == 2");
        goto fail;
    }
    if (PyArray_NDIM(index_array) != 1){
        PyErr_SetString(PyExc_ValueError, "Expected index array of ndim == 1");
        goto fail;
    }

    // Get sizes
    nlines = PyArray_DIM(species_array, 0);
    dtype = PyArray_DESCR(species_array);
    nchar = dtype->elsize;

    // Check sizes
    if (PyArray_DIM(linelist_array, 0) != 8){
        PyErr_SetString(PyExc_ValueError, "Expected linelist to have 8 values");
        goto fail;
    }
    if (PyArray_DIM(linelist_array, 1) != nlines){
        PyErr_SetString(PyExc_ValueError, "Expected all arrays to have the same length");
        goto fail;
    }
    if (PyArray_DIM(index_array, 0) != nlines){
        PyErr_SetString(PyExc_ValueError, "Expected all arrays to have the same length");
        goto fail;
    }

    // Get the pointers to pass to SMELIB
    linelist = (double*) PyArray_DATA(linelist_array);
    species_data = (char*) PyArray_DATA(species_array);
    index = (short *) PyArray_DATA(index_array);
    species = (IDL_STRING *) malloc(nlines * sizeof(IDL_STRING));

    for (int i = 0; i < nlines; i++)
    {
        species[i].slen = nchar;
        species[i].s = &species_data[nchar*i];
        species[i].stype = 0;
    }

    args_c[0] = &nlines;
    args_c[1] = species;
    args_c[2] = linelist;
    args_c[3] = index;
    result = UpdateLineList(n, args_c);

    free(species);
    Py_XDECREF(species_array);
    Py_XDECREF(linelist_array);
    Py_XDECREF(index_array);

    if (result != NULL && result[0] != OK_response)
    {
        PyErr_SetString(PyExc_RuntimeError, result);
        return NULL;
    }

    Py_RETURN_NONE;

fail:
    Py_XDECREF(species_array);
    Py_XDECREF(linelist_array);
    Py_XDECREF(index_array);
    return NULL;
}

static char smelib_InputModel_docstring[] = "Read in model atmosphere";
static PyObject * smelib_InputModel(PyObject * self, PyObject *args, PyObject *kwds)
{
    int n = 14;
    void * args_c[14];
    const char * result = NULL;
    double teff, grav, wlstd, radius = NAN;
    char * motype;
    IDL_STRING motype_idl;
    int nrhox;
    PyObject * opflag_obj = NULL, *depth_obj = NULL, *temp_obj = NULL;
    PyObject * xne_obj = NULL, *xna_obj = NULL, *rho_obj = NULL, *vt_obj = NULL;
    PyObject * height_obj = NULL;

    PyArrayObject * opflag_arr = NULL, *depth_arr = NULL, *temp_arr = NULL;
    PyArrayObject * xne_arr = NULL, *xna_arr = NULL, *rho_arr = NULL, *vt_arr = NULL;
    PyArrayObject * height_arr = NULL;

    // Need to make this constant because C++ needs it
    // but we cast to non constant and trust Python
    static const char * keywords[] = {"teff", "grav", "wlstd", "motype", "opflag", "depth", "temp", "xna", "xne", "rho", "vt", "radius", "height", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "dddsOOOOOOO|dO", const_cast<char **>(keywords), &teff, &grav, &wlstd, &motype, 
            &opflag_obj, &depth_obj, &temp_obj, &xna_obj, &xne_obj, &rho_obj, &vt_obj, &radius, &height_obj))
        return NULL;

    // Convert to Numpy arrays
    opflag_arr = (PyArrayObject *) PyArray_FROM_OTF(opflag_obj, NPY_SHORT, NPY_ARRAY_IN_ARRAY);
    depth_arr = (PyArrayObject *) PyArray_FROM_OTF(depth_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    temp_arr = (PyArrayObject *) PyArray_FROM_OTF(temp_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    xna_arr = (PyArrayObject *) PyArray_FROM_OTF(xna_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    xne_arr = (PyArrayObject *) PyArray_FROM_OTF(xne_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    rho_arr = (PyArrayObject *) PyArray_FROM_OTF(rho_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    vt_arr = (PyArrayObject *) PyArray_FROM_OTF(vt_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);

    if (opflag_arr == NULL || depth_arr == NULL || temp_arr == NULL || xna_arr == NULL ||
        xne_arr == NULL || rho_arr == NULL || vt_arr == NULL){
            goto fail;
        }

    // Only if given
    if (height_obj != NULL){
        height_arr = (PyArrayObject *) PyArray_FROM_OTF(height_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
        if (height_arr == NULL){
            goto fail;
        }
    }

    if (PyArray_NDIM(opflag_arr) != 1){
        PyErr_SetString(PyExc_ValueError, "Expected opflag array of ndim == 1");
        goto fail;
    }
    if (PyArray_NDIM(depth_arr) != 1){
        PyErr_SetString(PyExc_ValueError, "Expected depth of ndim == 1");
        goto fail;
    }
    if (PyArray_NDIM(temp_arr) != 1){
        PyErr_SetString(PyExc_ValueError, "Expected temp array of ndim == 1");
        goto fail;
    }
    if (PyArray_NDIM(xna_arr) != 1){
        PyErr_SetString(PyExc_ValueError, "Expected XNA array of ndim == 1");
        goto fail;
    }
    if (PyArray_NDIM(xne_arr) != 1){
        PyErr_SetString(PyExc_ValueError, "Expected XNE array of ndim == 1");
        goto fail;
    }
    if (PyArray_NDIM(rho_arr) != 1){
        PyErr_SetString(PyExc_ValueError, "Expected rho array of ndim == 1");
        goto fail;
    }
    if (PyArray_NDIM(vt_arr) != 1){
        PyErr_SetString(PyExc_ValueError, "Expected v_turb array of ndim == 1");
        goto fail;
    }
    if (height_arr != NULL){
        if (PyArray_NDIM(height_arr) != 1){
            PyErr_SetString(PyExc_ValueError, "Expected height array of ndim == 1");
            goto fail;
        }
    }

    // Get size
    nrhox = PyArray_DIM(depth_arr, 0);

    // Check size
    if (PyArray_DIM(opflag_arr, 0) != 20){
        PyErr_SetString(PyExc_ValueError, "Expected opflag array of size 20");
        goto fail;
    }
    if (PyArray_DIM(depth_arr, 0) != nrhox){
        PyErr_SetString(PyExc_ValueError, "Expected depth of the same size");
        goto fail;
    }
    if (PyArray_DIM(temp_arr, 0) != nrhox){
        PyErr_SetString(PyExc_ValueError, "Expected temp array of the same size");
        goto fail;
    }
    if (PyArray_DIM(xna_arr, 0) != nrhox){
        PyErr_SetString(PyExc_ValueError, "Expected XNA array of the same size");
        goto fail;
    }
    if (PyArray_DIM(xne_arr, 0) != nrhox){
        PyErr_SetString(PyExc_ValueError, "Expected XNE array of the same size");
        goto fail;
    }
    if (PyArray_DIM(rho_arr, 0) != nrhox){
        PyErr_SetString(PyExc_ValueError, "Expected rho array of the same size");
        goto fail;
    }
    if (PyArray_DIM(vt_arr, 0) != nrhox){
        PyErr_SetString(PyExc_ValueError, "Expected v_turb array of the same size");
        goto fail;
    }
    if (height_arr != NULL){
        if (PyArray_DIM(height_arr, 0) != nrhox){
            PyErr_SetString(PyExc_ValueError, "Expected height array of the same size");
            goto fail;
        }
    }

    // Convert to IDL string
    motype_idl.s = motype;
    motype_idl.slen = strlen(motype);
    motype_idl.stype = 0;

    if (strncmp(motype, "SPH", 3) == 0){
        if (isnan(radius) || height_arr == NULL){
            PyErr_SetString(PyExc_ValueError, "Model type is SPH but no height and/or radius are given");
            goto fail;
        }
        args_c[0] = &nrhox;
        args_c[1] = &teff;
        args_c[2] = &grav;
        args_c[3] = &wlstd;
        args_c[4] = &motype_idl;
        args_c[5] = &radius;
        args_c[6] = (short*) PyArray_DATA(opflag_arr);
        args_c[7] = (double*) PyArray_DATA(depth_arr);
        args_c[8] = (double*) PyArray_DATA(temp_arr);
        args_c[9] = (double*) PyArray_DATA(xne_arr);
        args_c[10] = (double*) PyArray_DATA(xna_arr);
        args_c[11] = (double*) PyArray_DATA(rho_arr);
        args_c[12] = (double*) PyArray_DATA(vt_arr);
        args_c[13] = (double*) PyArray_DATA(height_arr);
        n = 14;
    } else {
        args_c[0] = &nrhox;
        args_c[1] = &teff;
        args_c[2] = &grav;
        args_c[3] = &wlstd;
        args_c[4] = &motype_idl;
        args_c[5] = (short*) PyArray_DATA(opflag_arr);
        args_c[6] = (double*) PyArray_DATA(depth_arr);
        args_c[7] = (double*) PyArray_DATA(temp_arr);
        args_c[8] = (double*) PyArray_DATA(xne_arr);
        args_c[9] = (double*) PyArray_DATA(xna_arr);
        args_c[10] = (double*) PyArray_DATA(rho_arr);
        args_c[11] = (double*) PyArray_DATA(vt_arr);
        n = 12;
    }

    result = InputModel(n, args_c);

    Py_XDECREF(opflag_arr);
    Py_XDECREF(depth_arr);
    Py_XDECREF(temp_arr);
    Py_XDECREF(xna_arr);
    Py_XDECREF(xne_arr);
    Py_XDECREF(rho_arr);
    Py_XDECREF(vt_arr);
    Py_XDECREF(height_arr);

    if (result != NULL && result[0] != OK_response)
    {
        PyErr_SetString(PyExc_RuntimeError, result);
        return NULL;
    }

    Py_RETURN_NONE;

fail:
    Py_XDECREF(opflag_arr);
    Py_XDECREF(depth_arr);
    Py_XDECREF(temp_arr);
    Py_XDECREF(xna_arr);
    Py_XDECREF(xne_arr);
    Py_XDECREF(rho_arr);
    Py_XDECREF(vt_arr);
    Py_XDECREF(height_arr);
    return NULL;
}

static char smelib_InputDepartureCoefficients_docstring[] = "Input Departure coefficients";
static PyObject * smelib_InputDepartureCoefficients(PyObject * self, PyObject *args)
{
    int n = 2;
    void * args_c[n];
    const char * result = NULL;
    int linenumber;
    short nrhox;
    PyObject * bmatrix_obj = NULL;
    PyArrayObject * bmatrix_arr = NULL;

    if (!PyArg_ParseTuple(args, "iO", &linenumber, &bmatrix_obj))
        return NULL;

    bmatrix_arr = (PyArrayObject *) PyArray_FROM_OTF(bmatrix_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);

    if (bmatrix_arr == NULL){
        return NULL;
    }

    nrhox = GetNRHOX();

    if (PyArray_NDIM(bmatrix_arr) != 2){
        Py_XDECREF(bmatrix_arr);
        PyErr_SetString(PyExc_ValueError, "Expected bmatrix with ndim == 2");
        return NULL;
    }
    if (PyArray_DIM(bmatrix_arr, 0) != 2){
        Py_XDECREF(bmatrix_arr);
        PyErr_SetString(PyExc_ValueError, "Expected bmatrix with shape (2, nrhox)");
        return NULL;
    }
    if (PyArray_DIM(bmatrix_arr, 1) != nrhox){
        Py_XDECREF(bmatrix_arr);
        PyErr_SetString(PyExc_ValueError, "Expected bmatrix with shape (2, nrhox)");
        return NULL;
    }

    args_c[0] = (double*) PyArray_DATA(bmatrix_arr);
    args_c[1] = &linenumber;
    result = InputDepartureCoefficients(n, args_c);

    Py_DECREF(bmatrix_arr);

    if (result != NULL && result[0] != OK_response)
    {
        PyErr_SetString(PyExc_RuntimeError, result);
        return NULL;
    }

    Py_RETURN_NONE;
}

static char smelib_GetDepartureCoefficients_docstring[] = "Get NLTE b's for specific line";
static PyObject * smelib_GetDepartureCoefficients(PyObject * self, PyObject *args)
{
    const int n = 3;
    void * args_c[n];
    const char * result = NULL;
    int nrhox, linenumber;

    PyArrayObject * bmatrix_arr;
    double * bmatrix_data;

        
    if (!PyArg_ParseTuple(args, "i", &linenumber))
        return NULL;

    nrhox = GetNRHOX();
    npy_intp dims[2] = {2, nrhox};    
    bmatrix_arr = (PyArrayObject*) PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    bmatrix_data = (double*) PyArray_DATA(bmatrix_arr);

    args_c[0] = bmatrix_data;
    args_c[1] = &nrhox;
    args_c[2] = &linenumber;
    result = GetDepartureCoefficients(n, args_c);

    if (result != NULL && result[0] != OK_response)
    {
        Py_DecRef((PyObject*)bmatrix_arr);
        PyErr_SetString(PyExc_RuntimeError, result);
        return NULL;
    }

    return (PyObject*) bmatrix_arr;
}


static PyMethodDef module_methods[] = {
    {"LibraryVersion", smelib_LibraryVersion, METH_VARARGS, smelib_LibraryVersion_docstring},
    {"GetDataFiles", smelib_GetDataFiles, METH_VARARGS, smelib_GetDataFiles_docstring},
    {"GetLibraryPath", smelib_GetLibraryPath, METH_VARARGS, smelib_GetLibraryPath_docstring},
    {"SetLibraryPath", smelib_SetLibraryPath, METH_VARARGS, smelib_SetLibraryPath_docstring},
    {"InputWaveRange", smelib_InputWaveRange, METH_VARARGS, smelib_InputWaveRange_docstring},
    {"SetVWscale", smelib_SetVWscale, METH_VARARGS, smelib_SetVWscale_docstring},
    {"SetH2broad", smelib_SetH2broad, METH_VARARGS, smelib_SetH2broad_docstring},
    {"ClearH2broad", smelib_ClearH2broad, METH_VARARGS, smelib_ClearH2broad_docstring},
    {"InputLineList", smelib_InputLineList, METH_VARARGS, smelib_InputLineList_docstring},
    {"OutputLineList", smelib_OutputLineList, METH_VARARGS, smelib_OutputLineList_docstring},
    {"UpdateLineList", smelib_UpdateLineList, METH_VARARGS, smelib_UpdateLineList_docstring},
    {"InputModel", (PyCFunction)(void(*)(void))smelib_InputModel, METH_VARARGS|METH_KEYWORDS, smelib_InputModel_docstring},
    {"InputDepartureCoefficients", smelib_InputDepartureCoefficients, METH_VARARGS, smelib_InputDepartureCoefficients_docstring},
    {"GetDepartureCoefficients", smelib_GetDepartureCoefficients, METH_VARARGS, smelib_GetDepartureCoefficients_docstring},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC PyInit__smelib(void)
{
    static struct PyModuleDef smelibmodule = {
    PyModuleDef_HEAD_INIT,
    "smelib",          /* name of module */
    module_docstring,  /* module documentation, may be NULL */
    -1,                /* size of per-interpreter state of the module,
                          or -1 if the module keeps state in global variables. */
    module_methods
    };

    PyObject *m = PyModule_Create(&smelibmodule);
    if (m == NULL)
        return NULL;

    /* Load `numpy` functionality. */
    import_array();

    return m;
}
