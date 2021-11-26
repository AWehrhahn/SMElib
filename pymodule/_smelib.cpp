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




static PyMethodDef module_methods[] = {
    {"LibraryVersion", smelib_LibraryVersion, METH_VARARGS, smelib_LibraryVersion_docstring},
    {"GetDataFiles", smelib_GetDataFiles, METH_VARARGS, smelib_GetDataFiles_docstring},
    {"GetLibraryPath", smelib_GetLibraryPath, METH_VARARGS, smelib_GetLibraryPath_docstring},
    {"SetLibraryPath", smelib_SetLibraryPath, METH_VARARGS, smelib_SetLibraryPath_docstring},
    {"InputWaveRange", smelib_InputWaveRange, METH_VARARGS, smelib_InputWaveRange_docstring},
    {"SetVWscale", smelib_SetVWscale, METH_VARARGS, smelib_SetVWscale_docstring},
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
