======
CyCuba
======

A Cython wrapper around the Cuba multidimensional integration library.

What is Cuba?
-------------

The Cuba integration library (http://www.feynarts.de/cuba/) provides both 
Monte Carlo and deterministic rules for the evaluation of multidimensional 
integrals.

What does CyCuba provide?
-------------------------

CyCuba provides an interface between the Python interpreter and the Cuba
integration library. It allows the user to quickly and easily integrate their
Python functions, including for situations where the SciPy quadrature routines
(http://docs.scipy.org/doc/scipy/reference/integrate.html) are inappropriate due
to excessive dimensionality or non-smooth integrand functions. CyCuba is
intended to greatly simplify the process of evaluating these integrals, while
simultaneously providing transparent access to all of the different options of
the Cuba library for those users who need it.

Notably, CyCuba does NOT yet support concurrent evaluation of integrals, and
testing for some integration options remains incomplete. When you find bugs or
incomplete features that you need, please post an issue to GitHub
(https://github.com/woodscn/CyCuba), and we will try to resolve it right away.

Requirements
------------

- Cuba
- cython
- numpy
- pytest


Installation
------------

CyCuba relies on a working copy of the Cuba library. The most current version of
this library can be obtained from http://www.feynarts.de/cuba/, along with
papers which describe in detail the algorithms used in Cuba, along with the
various options and capabilities the library provides. This library and its
accompanying header files should be placed in a "normal" location for your
system, or else environment variables can be modified in order to point to their
location. For instance, if the library is placed in ``~/lib/`` and the library
header file ``cuba.h`` is placed in ``~/include/``, then one could write:

    export CFLAGS="-I ~/include/ $CFLAGS"
    export LDFLAGS="-L ~/lib/ $LDFLAGS"

Once this library is available, the CyCuba package can be
installed using the normal command: ``python setup.py install``. You will need
to have the ``cython``, ``numpy``, and ``pytest`` packages installed.  The
package tests may be run by executing ``py.test`` in the installation
directory.



References
----------
- Cuba - a library for multidimensional numerical integration
(http://arxiv.org/abs/hep-ph/0404043)

- Concurrent Cuba (http://arxiv.org/abs/1408.6373)


