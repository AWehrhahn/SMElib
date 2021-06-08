from os.path import join, realpath, dirname
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

lib_path = realpath(join(dirname(__file__), "../lib"))
file_path = realpath(join(dirname(__file__), "smelib.pyx"))

examples_extension = Extension(
    name="smelib",
    sources=[file_path],
    libraries=["gfortran", "sme"],
    library_dirs=[lib_path],
    runtime_library_dirs=[lib_path],
)
setup(
    name="smelib",
    ext_modules=cythonize([examples_extension], language_level=3)
)