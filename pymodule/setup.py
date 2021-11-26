from os.path import join, dirname, abspath
from distutils.core import setup, Extension
import numpy.distutils.misc_util

libdir = abspath(join(dirname(__file__), "../lib"))
include_dirs = numpy.distutils.misc_util.get_numpy_include_dirs()
# include_dirs += ["/usr/include/python3.8"]
include_dirs += [libdir]

module = Extension(
    "_smelib",
    sources=["_smelib.cpp"],
    language="c++",
    include_dirs=include_dirs,
    libraries=["sme"],
    library_dirs=[libdir],
)

setup(ext_modules=[module])

