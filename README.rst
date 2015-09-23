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
header file ``cuba.h`` is placed in ``~/include/``, then one could write
(OS X or Linux): 

    export CFLAGS="-I ~/include/ $CFLAGS"
    export LDFLAGS="-L ~/lib/ $LDFLAGS"

Once this library is available, the CyCuba package can be
installed using the normal command: ``python setup.py install``. You will need
to have the ``cython``, ``numpy``, and ``pytest`` packages installed.  The
package tests may be run as ``cycuba.test()``.


Usage
-----
CyCuba provides wrappers to the four Cuba routines: Vegas, Suave, Divonne, and
Cuhre. Detailed descriptions of these routines are available in the Cuba
documentation (Ref. 1), and only the necessary Python equivalents are given
here.

Simplest possible usage:

``Vegas/Suave/Divonne/Cuhre(integrand, ranges=None)``

- ``integrand``: Python callable object. The signature is y = f(*args), where y
  is an iterable. The Python wrapper handles the necessary conversion to the
  required Cuba form. The ``userdata`` pointer is not available; if context must
  be provided to ``integrand``, then it should be defined as a callable class,
  and the context stored as state. The additional information Cuba provides,
  e.g. ``nvec``, or ``weight`` for Vegas, are not supported at this time.

- ``ranges``: Python iterable of len(2) iterables. Defines hypercubic
  integration domain, e.g. [[x0min, x0max], [x1min, x1max], ... [xnmin, xnmax]]
  If unspecified, defaults to [[0, 1], [0, 1], ...]



Common arguments:
- ``epsrel=1e-3``, ``epsabs=1e-12``: The requested relative and absolute accuracies.

- ``mineval=0``, ``maxeval=1e6``: The minimum and maximum allowed number of
  ``integrand`` evaluations.

- ``verbosity=0``: The requested verbosity level, ranging from 0 (lowest) to 3
  (highest).

- ``last_samples_only=False``: Use only the last, largest set of samples to
  compute the final integral value.

- ``do_not_smooth=False`` (Vegas and Suave only): Do not smooth the importance
  function. Recommended for non-smooth integrands.

- ``retain_state_file=False``: Do not delete the state file upon successful
  termination of the integration. Note that the use of state files is currently
  unsupported.

- ``file_grid_only=False``: Ignore integrator state, even if state file is
  present. Allows the reuse of one grid for a different integrand. Note that the
  use of state files is currently unsupported.

- ``seed=0``: Seed the random number generator. Uses Sobol quasi-random numbers
  if set to 0. If non-zero, the specific random number generator is determined
  by ``level``.

- ``level=0``: Further specify a random number generator, provided ``seed`` is
  non-zero. If ``level`` is set to zero, use Mersenne Twister. If non-zero, use
  Ranlux with a generation period ``p`` defined by ``level`` as follows:
  - 1: ``p`` = 48: Very long period; passes gap test but fails spectral test. 
  - 2: ``p`` = 97: Passes all known tests, but theoretically still defective.
  - 3: ``p`` = 223: Very small chance of observing theoretical correlations.
  - 4: ``p`` = 389: Highest possible luxury.
  - 5-23: Default to 3.
  - 24+: Specify ``p`` directly.

``def Vegas(integrand, ranges=None, **kwargs)``

Vegas-specific arguments:

  - ``nstart=1e3``: Starting number of integrand evaluations per iteration.
    
  - ``nincrease=5e2``: Increase in number of integrand evaluations per
    iteration.
    
  - ``nbatch=1e3``: Batch size for sampling.
    
  - ``gridno=0``: Slot in internal grid table. See Cuba documentation.

``def Suave(integrand, ranges=None, **kwargs)``

Suave-specific arguments:

 - ``nnew=1e3``: Number of new integrand evaluations in each subdivision.

 - ``nmin=2``: Minimum number of samples a former pass must congtribute to a
   subregion to be considered in that region's compound integration
   value. Increasing ``nmin`` may reduce jumps in the chi-squared value.

 - ``flatness=25``: The type of norm used to compute the fluctuation of a
   sample. Choose a large value for flat integrands, and a smaller value for
   volatile integrands. Very large values (> ~ 200) may cause double-precision
   overflow.

``def Divonne(integrand, ranges=None, **kwargs)``

The user is cautioned that Divonne is the most complex of the routines in Cuba,
and testing for this part of the wrapper does not provide complete coverage of
all the available options. In particular, specification of ``xgiven`` and
``peakfinder`` is untested, and may contain bugs. The developers welcome any
test routines that may be contributed to extend test coverage for Divonne (and
the other routines). 

Divonne-specific arguments (See Ref. 1 for description):

 - ``key1=47``

 - ``key2=1``

 - ``key3=1``
   
 - ``maxpass=5``

 - ``border=0``

 - ``maxchisq=10``

 - ``mindeviation=0.25``

 - ``xgiven=[]``

 - ``nextra=0``

 - ``peakfinder=None``

``def Cuhre(integrand, ranges=None, **kwargs)``

Cuhre-specific arguments:

 - ``key=1``: Select the cubature rule of degree ``key``. Available choices are
   7, 9, 11 (3-dimensions only), 13 (2-dimensions only). For other values, the
   highest available rule for the dimensionality is used.


References
----------
1. Cuba - a library for multidimensional numerical integration
(http://arxiv.org/abs/hep-ph/0404043)

2. Concurrent Cuba (http://arxiv.org/abs/1408.6373)


