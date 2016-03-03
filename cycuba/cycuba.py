import sys
import inspect
from functools import reduce
from os import linesep, environ, path
from warnings import warn, simplefilter
from operator import mul

import pytest

from _cycuba import _cuba

simplefilter('always')
# Parallelization within Cuba is not yet supported. Testing procedures are
# required in order to make certain that the Python wrapper is thread-safe,
# since the PyMultiNest project found that there were problems.
environ['CUBACORES'] = '1'

ibf_v = 0
ibf_lso = False
ibf_dns = False
ibf_rsf = False
ibf_fgo = False
ibf_l = 0


def get_ndim(integrand):
    if sys.version_info[0] < 3: # Python 2 Compatibility
        if inspect.isfunction(integrand):
            out = len(inspect.getargspec(integrand)[0])
        else:
            out = len(inspect.getargspec(integrand.__call__)[0]) - 1
    else:
        out = len(inspect.signature(integrand).parameters)
    return out

    
class ScaledIntegrand(object):
    def __init__(self, f, ranges):
        self.f = f
        self.ranges = ranges
        self.diffs = [range_[1] - range_[0] for range_ in ranges]
        self.jacobian = reduce(mul, self.diffs, 1)

    def __call__(self, *args):
        if len(args) != len(self.ranges):
            raise CyCubaError("different lengths of args and ranges!")
        scaled_args = [range_[0] + x * diff for (range_, diff, x) in
                       zip(self.ranges, self.diffs, args)]
        return [item * self.jacobian for item in self.f(*scaled_args)]


def integer_bit_flags(verbosity, last_samples_only, do_not_smooth_vs,
                      retain_state_file, file_grid_only_v, level):
    level_key = 3 if level > 4 and level < 24 else level
    flag_string = (
        format(level_key, 'b') +
        {True: '1', False: '0'}[file_grid_only_v] +
        {True: '1', False: '0'}[retain_state_file] +
        {True: '1', False: '0'}[do_not_smooth_vs] +
        {True: '1', False: '0'}[last_samples_only] +
        {0: '00', 1: '01', 2: '10', 3: '11'}[verbosity]
    )
    return int(flag_string, base=2)


class CyCubaError(Exception):
    pass


class CyCubaWarning(UserWarning):
    pass


class CyCubaIntegration(object):
    def __init__(self, integrand, ranges, epsrel=1e-3, epsabs=1e-12, seed=0,
                 mineval=0, maxeval=1e6, statefile=""):
        self.ndim = get_ndim(integrand)
        self.ranges = ranges
        # Set things up for scaling integrands
        if self.ranges:
            self.integrand = ScaledIntegrand(integrand, ranges)
        else:
            self.integrand = integrand
        self.test_args = [0.5 for ind in range(self.ndim)]
        self.ncomp = len(self.integrand(*self.test_args))
        self.nvec = 1  # Vectorized integrands not yet supported.
        self.epsrel = epsrel
        self.epsabs = epsabs
        self.seed = seed
        self.mineval = mineval
        self.maxeval = maxeval
        self.statefile = statefile
        self.spin = 0  # Spin variable is not yet supported.
        self.empty_vegas_args = (0, 0, 0, 0)
        self.empty_suave_args = (0, 0, 0)
        self.empty_divonne_args = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        self.empty_cuhre_args = (0,)

    def base_args(self, flags):
        return (self.ndim, self.ncomp, self.integrand, self.nvec,
                self.epsrel, self.epsabs, flags, self.seed, self.mineval,
                self.maxeval, self.statefile, self.spin)

    def _vegas(self, flags, *args):
        vegas_args = (self.base_args(flags) +
                      args +
                      self.empty_suave_args +
                      self.empty_divonne_args +
                      self.empty_cuhre_args)
        out = _cuba('vegas', *vegas_args)
        [self.neval, self.fail, self.integral, self.error, self.prob] = out
        self._check_fail()
        return [self.integral, self.error, self.prob, self.neval]

    def _suave(self, flags, *args):
        suave_args = (self.base_args(flags) +
                      self.empty_vegas_args +
                      args +
                      self.empty_divonne_args +
                      self.empty_cuhre_args)
        out = _cuba('suave', *suave_args)
        [self.neval, self.nregions, self.fail,
         self.integral, self.error, self.prob] = out
        self._check_fail()
        return [self.integral, self.error, self.prob, self.neval]

    def _divonne(self, flags, *args):
        divonne_args = (self.base_args(flags) +
                        self.empty_vegas_args +
                        self.empty_suave_args +
                        args +
                        self.empty_cuhre_args)
        out = _cuba('divonne', *divonne_args)
        [self.neval, self.nregions, self.fail,
         self.integral, self.error, self.prob] = out
        self._check_fail()
        return [self.integral, self.error, self.prob, self.neval]

    def _cuhre(self, flags, *args):
        cuhre_args = (self.base_args(flags) +
                      self.empty_vegas_args +
                      self.empty_suave_args +
                      self.empty_divonne_args +
                      args)
        out = _cuba('cuhre', *cuhre_args)
        [self.neval, self.nregions, self.fail,
         self.integral, self.error, self.prob] = out
        try:
            self._check_fail()
        finally:
            return [self.integral, self.error, self.prob, self.neval]

    def _check_fail(self):
        if self.fail == -1:
            raise CyCubaError("Dimension out of range!")
        elif self.fail > 0:
            if self.fail == 1:
                warn("Accuracy goal not achieved!", CyCubaWarning)
            else:
                warn(
                    "Accuracy goal not achieved! " + linesep +
                    "Approximately "+str(self.fail)+" additional points" +
                    " required to reach desired accuracy.", CyCubaWarning)



def Vegas(integrand, ranges=None, nstart=1e3, nincrease=5e2, nbatch=1e3,
          gridno=0, verbosity=ibf_v, last_samples_only=ibf_lso,
          do_not_smooth=ibf_dns, retain_state_file=ibf_rsf,
          file_grid_only=ibf_fgo, level=ibf_l, **kwargs):
    cycuba_integration = CyCubaIntegration(integrand, ranges, **kwargs)
    flags = integer_bit_flags(
        verbosity=verbosity, last_samples_only=last_samples_only,
        do_not_smooth_vs=do_not_smooth, retain_state_file=retain_state_file,
        file_grid_only_v=file_grid_only, level=level)
    out = cycuba_integration._vegas(
        flags, nstart, nincrease, nbatch, gridno)
    return out


def Suave(integrand, ranges=None, nnew=1000, nmin=2, flatness=25,
          verbosity=ibf_v, last_samples_only=ibf_lso, do_not_smooth=ibf_dns,
          retain_state_file=ibf_rsf, level=ibf_l,
          **kwargs):
    cycuba_integration = CyCubaIntegration(integrand, ranges, **kwargs)
    flags = integer_bit_flags(
        verbosity=verbosity, last_samples_only=last_samples_only,
        do_not_smooth_vs=do_not_smooth, retain_state_file=retain_state_file,
        file_grid_only_v=False, level=level)
    out = cycuba_integration._suave(flags, nnew, nmin, flatness)
    return out


def Divonne(integrand, ranges=None, key1=47, key2=1, key3=1, maxpass=5,
            border=0, maxchisq=10, mindeviation=.25, xgiven=[], nextra=0,
            peakfinder=None,
            verbosity=ibf_v,
            last_samples_only=ibf_lso, retain_state_file=ibf_rsf,
            level=ibf_l,
            **kwargs):
    cycuba_integration = CyCubaIntegration(integrand, ranges, **kwargs)
    flags = integer_bit_flags(
        verbosity=verbosity, last_samples_only=last_samples_only,
        do_not_smooth_vs=False, retain_state_file=retain_state_file,
        file_grid_only_v=False, level=level)
    ngiven = len(xgiven)
    try:
        ldxgiven = len(xgiven[0])
    except IndexError:
        ldxgiven = cycuba_integration.ndim
    out = cycuba_integration._divonne(
        flags, key1, key2, key3, maxpass, border, maxchisq, mindeviation,
        ngiven, ldxgiven, xgiven, nextra, peakfinder)
    return out


def Cuhre(integrand, ranges=None, key=0, verbosity=ibf_v,
          last_samples_only=ibf_lso, retain_state_file=ibf_rsf,
          **kwargs):
    if 'seed' in kwargs:
        warn("Cuhre is a deterministic routine. Ignoring value for 'seed'",
             CyCubaWarning)
    cycuba_integration = CyCubaIntegration(integrand, ranges, **kwargs)
    flags = integer_bit_flags(
        verbosity=verbosity, last_samples_only=last_samples_only,
        do_not_smooth_vs=False, retain_state_file=retain_state_file,
        file_grid_only_v=False, level=0)
    out = cycuba_integration._cuhre(flags, key)
    return out

def test():
    pytest.main(path.dirname(path.abspath(__file__)))
