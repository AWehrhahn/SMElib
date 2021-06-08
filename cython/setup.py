from distutils.core import setup
from distutils.extension import Extension
from os.path import dirname, join, realpath

import numpy as np
from Cython.Build import cythonize

lib_path = realpath(join(dirname(__file__), "lib"))
file_path = realpath(join(dirname(__file__), "smelib.pyx"))

print("lib_path: " + lib_path)
print("file_path: " + file_path)

examples_extension = Extension(
    name="smelib",
    sources=[file_path],
    libraries=["gfortran", "sme"],
    library_dirs=[lib_path],
    runtime_library_dirs=[lib_path],
    include_dirs=[np.get_include()],
)
setup(
    name="smelib",
    ext_modules=cythonize([examples_extension], language_level=3)
)
