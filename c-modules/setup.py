from distutils.core import setup
from distutils.extension import Extension
 
setup(name="STATIUM C++ Query Distances Module", ext_modules=[Extension("statium_cpp", ["extensions.cpp"], libraries = ["boost_python"], extra_compile_args=['-std=c++11'])])
