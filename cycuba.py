from functools import reduce
from inspect import signature

from _cycuba import _cuba


ibf_v = 0
ibf_lso = False
ibf_dns = False
ibf_rsf = False
ibf_fgo = False
ibf_l = 0


class ScaledIntegrand(object):
    def __init__(self, f, ranges):
        self.f = f
        self.ranges = ranges
        self.diffs = [range_[1] - range_[0] for range_ in ranges]
        self.jacobian = reduce(operator.mul, self.diffs, 1)

    def __call__(self, *args):
        if len(args) == len(ranges):
            raise CyCubaError("different lengths of args and ranges!")
        scaled_args = [range_[0] + x * diff for (range_, diff, x) in
                       zip(self.ranges, self.diffs, args)]
        return self.f(*scaled_args) * self.jacobian
    
    
def integer_bit_flags(verbosity, last_samples_only, do_not_smooth_vs,
                      retain_state_file, file_grid_only_v, level):
    level_key = 3 if level > 4 and level < 24 else level
    flag_string = (
        format(level_key, 'b')
        + {True: '0', False: '1'}[file_grid_only_v]
        + {True: '0', False: '1'}[retain_state_file]
        + {True: '0', False: '1'}[do_not_smooth_vs]
        + {True: '0', False: '1'}[last_samples_only]
        + {0: '00', 1: '01', 2: '10', 3: '11'}[verbosity]
    )
    return int(flag_string,base=2)


class CyCubaError(Exception):
    pass


class CyCubaIntegration(object):
    def __init__(self, integrand, ranges, epsrel=1e-3, epsabs=1e-12, seed=0,
                 mineval=0, maxeval=1e6, statefile=""):
        self.integrand = integrand
        self.ndim = len(signature(integrand).parameters)
        self.ranges = ranges
        # Set things up for scaling integrands
        if self.ranges:
            self.diffs = [range_[1] - range_[0] for range_ in ranges]
            self.jacobian = reduce(operator.mul, self.diffs, 1)
            self.test_args = [range_[0] + 0.5 * diff
                              for range_, diff in zip(ranges, diffs)]
        else:
            self.test_args = [0.5 for ind in range(self.ndim)]
        self.ncomp = len(self.integrand(*self.test_args))
        self.nvec = 1 # Vectorized integrands not yet supported.
        self.epsrel = epsrel
        self.epsabs = epsabs
        self.seed = seed
        self.mineval = mineval
        self.maxeval = maxeval
        self.statefile = statefile
        self.spin = 0 # Spin variable is not yet supported.
        self.empty_vegas_args = (0, 0, 0, 0)
        self.empty_suave_args = (0, 0, 0)
        self.empty_divonne_args = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        self.empty_cuhre_args = (0,)

    def base_args(self, flags):
        return (self.ndim, self.ncomp, self.integrand, self.nvec,
                self.epsrel, self.epsabs, flags, self.seed, self.mineval,
                self.maxeval, self.statefile, self.spin)

    def _vegas(self, flags, nstart, nincrease, nbatch, gridno):
        args = (self.base_args(flags) +
                (nstart, nincrease, nbatch, gridno) +
                self.empty_suave_args + 
                self.empty_divonne_args +
                self.empty_cuhre_args)
        return _cuba('vegas', *args)

    def _cuhre(self, flags, key):
        args = (self.base_args(flags) +
                self.empty_vegas_args +
                self.empty_suave_args + 
                self.empty_divonne_args +
                (key,))
        return _cuba('cuhre', *args)

                     
def Vegas(integrand, ranges=None, nstart=1e3, nincrease=5e2, nbatch=1e3,
          gridno=0, verbosity=ibf_v, last_samples_only=ibf_lso,
          do_not_smooth=ibf_dns, retain_state_file=ibf_rsf,
          file_grid_only=ibf_fgo, level=ibf_l, **kwargs):
    cycuba_integration = CyCubaIntegration(integrand, ranges, **kwargs)
    flags = integer_bit_flags(verbosity, last_samples_only, do_not_smooth,
                              retain_state_file, file_grid_only, level)
    return cycuba_integration._vegas(flags, nstart, nincrease, nbatch, gridno)


def Cuhre(integrand, ranges=None, key=0, verbosity=ibf_v,
          last_samples_only=ibf_lso, do_not_smooth=ibf_dns,
          retain_state_file=ibf_rsf, file_grid_only=ibf_fgo, level=ibf_l,
          **kwargs):
    cycuba_integration = CyCubaIntegration(integrand, ranges, **kwargs)
    flags = integer_bit_flags(verbosity, last_samples_only, do_not_smooth,
                              retain_state_file, file_grid_only, level)
    return cycuba_integration._cuhre(flags, key)


if __name__ == "__main__":
    test_function = lambda x, y: [1 if x**2 + y**2 < 1 else 0]
    out = Vegas(test_function, verbosity=0)
    print('Vegas results: ', out)
    out = Cuhre(test_function, verbosity=0)
    print('Cuhre results: ', out)
    
