from setuptools import setup, Extension
import os
import sys
import subprocess
CYTHONIZE = True
cwd = os.path.abspath(os.path.dirname(__file__))
cython_directory = os.path.join(cwd, 'cycuba')
extension_filenames =  [os.path.join(cython_directory, '_cycuba')]
extension_sources = []
for filename in extension_filenames:
    file_name_c = filename + '.c'
    if not os.path.isfile(filename) or CYTHONIZE:
        # Call Cython on the associated *.pyx file
        p = subprocess.call(['cython', filename + '.pyx'], cwd=cython_directory)
        if p != 0:
            raise RuntimeError("Cythonize failed!")
    extension_sources.append(file_name_c)


extensions = [Extension('_cycuba', extension_sources, libraries=['cuba'])]

setup(name='cycuba', version='0.1.0',
      description=
      'A Cython wrapper to the Cuba library (http://www.feynarts.de/cuba/)',
      long_description=None, author='C. Nathan Woods',
      author_email='woodscn@lanl.gov', license='BSD',
      ext_modules=extensions
      )
