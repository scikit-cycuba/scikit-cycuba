from distutils.core import setup, Extension
from Cython.Build import cythonize


setup(name='cycuba', version='0.1.0',
      description=
      'A Cython wrapper to the Cuba library (http://www.feynarts.de/cuba/)',
      long_description=None, author='C. Nathan Woods',
      author_email='woodscn@lanl.gov', license='BSD',
      ext_modules=cythonize(
          [Extension('_cycuba', ['_cycuba.pyx'], libraries=['cuba'])])
      )
