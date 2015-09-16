from distutils.core import setup, Extension
from Cython.Build import cythonize


setup(name='cycuba', version='0.0.1',
      description=
      'A Cython wrapper to the Cuba library (http://www.feynarts.de/cuba/)',
      long_description=None, author='C. Nathan Woods',
      author_email='woodscn@lanl.gov', license='LGPL',
      install_requires=['cython'],
      ext_modules=cythonize(
          [Extension('cycuba', ['cycuba.pyx'], library_dirs=['./Cuba-4.2/'],
                     libraries=['cuba'], include_dirs=['./Cuba-4.2/'])],
          gdb_debug=True)
      )
