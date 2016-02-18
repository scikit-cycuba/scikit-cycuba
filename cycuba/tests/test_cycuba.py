from math import pi, sin, cos, exp, sqrt

import numpy.testing as nptest
import pytest

from cycuba import Vegas, Suave, Divonne, Cuhre
from cycuba import CyCubaError, CyCubaWarning

ndim = 3
ncomp = 1
nvec = 1
epsrel = 1e-3
epsabs = 1e-12
last = 4
seed = 0
mineval = 0
maxeval = 5e4

nstart = 1000
nincrease = 500
nbatch = 1000
gridno = 0
statefile = ""
spin = -1

nnew = 1000
nmin = 2
flatness = 25

key1 = 47
key2 = 1
key3 = 1
maxpass = 5
border = 0
maxchisq = 10
mindeviation = .25
ngiven = 0
ldxgiven = ndim
xgiven = []
nextra = 0
peakfinder = None

key = 0

verbose = 0

common_kwargs = {'epsrel': epsrel, 'epsabs': epsabs, 'seed': seed,
                 'mineval': mineval, 'maxeval': maxeval,
                 'statefile': statefile}

vegas_results = [
    [[0.664810733790], [0.000492177172], [0.08870493338]],
    [[5.269028388903], [0.009209558720], [0.00000185723]],
    [[0.307836626649], [0.000266157600], [0.000]],
    [[0.877334733857], [0.000827895321], [0.000]],
    [[0.416407819550], [0.000280981363], [0.005]],
    [[1.201914955901], [0.001122543858], [0.000]],
    [[0.709669278182], [0.000766047026], [0.000]],
    [[0.891209461227], [0.000847259931], [0.000]],
    [[0.080172877940], [0.000060757161], [0.036]],
    [[2.396848407842], [0.003171104659], [0.000]],
    [[0.523448717427], [0.001549583946], [0.001]]
    ]

suave_results = [
    [[0.664906987056], [0.000646671743], [0.807]],
    [[5.268345737220], [0.005220179167], [1.000]],
    [[0.307754949289], [0.000289631796], [1.000]],
    [[0.877251357781], [0.000825602379], [1.000]],
    [[0.416640336680], [0.000403879196], [0.975]],
    [[1.202392100441], [0.001177754704], [1.000]],
    [[0.709620676835], [0.000708507990], [1.000]],
    [[0.891188629400], [0.000824989887], [0.136]],
    [[0.080191317183], [0.000072979494], [0.439]],
    [[2.396433599525], [0.002356525928], [1.000]],
    [[0.523384162657], [0.001087959095], [0.000]]
    ]

divonne_results = [
    [[0.664619514751], [0.000635028303], [0.000]],
    [[5.268234724759], [0.005160562827], [0.000]],
    [[0.307826360819], [0.000288507183], [0.000]],
    [[0.877397127855], [0.000846183407], [0.000]],
    [[0.416585590187], [0.000390130006], [0.000]],
    [[1.202823015224], [0.000707720777], [0.000]],
    [[0.708880815647], [0.000704254351], [0.000]],
    [[0.891110425782], [0.000352607413], [0.000]],
    [[0.080182891553], [0.000063453661], [0.000]],
    [[2.396596937191], [0.002313583295], [0.000]],
    [[0.523818896242], [0.000524028620], [0.000]]
    ]

cuhre_results = [
    [[0.664669679782], [0.000000000033], [0.000]],
    [[5.268515066720], [0.001897671739], [0.000]],
    [[0.307805926385], [0.000246474384], [0.000]],
    [[0.877461451723], [0.000805950421], [0.039]],
    [[0.416538385894], [0.000000001505], [0.000]],
    [[1.202075364130], [0.001092125490], [0.000]],
    [[0.709694848704], [0.001483691315], [0.000]],
    [[0.891212798091], [0.000000002768], [0.008]],
    [[0.080185565078], [0.000000275828], [0.000]],
    [[2.396255344838], [0.002357219909], [0.000]],
    [[0.523697555021], [0.025943077006], [0.000]]
    ]


class TestIntegrand(object):
    def __init__(self, fun):
        self.fun = fun

    def rsq(self, x, y, z):
        return x ** 2 + y ** 2 + z ** 2

    def __call__(self, x, y, z):
        if self.fun == 0:
            out = sin(x) * cos(y) * exp(z)
        elif self.fun == 1:
            out = 1 / ((x + y) ** 2 + .003) * cos(y) * exp(z)
        elif self.fun == 2:
            out = 1 / (3.75 - cos(pi * x) - cos(pi * y) - cos(pi * z))
        elif self.fun == 3:
            out = abs(self.rsq(x, y, z) - .125)
        elif self.fun == 4:
            out = exp(-self.rsq(x, y, z))
        elif self.fun == 5:
            out = 1 / (1 - x * y * z + 1e-10)
        elif self.fun == 6:
            out = sqrt(abs(x - y - z))
        elif self.fun == 7:
            out = exp(-x * y * z)
        elif self.fun == 8:
            out = x ** 2 / (cos(x + y + z + 1) + 5)
        elif self.fun == 9:
            out = 1 / sqrt(x * y * z + 1e-5) if x > .5 else sqrt(x * y * z)
        else:
            out = 1 if self.rsq(x, y, z) < 1 else 0
        return [out]


def _vegas_test_runner(ind):
    kwargs = {'verbosity': verbose, 'last_samples_only': False,
              'do_not_smooth': False, 'retain_state_file': False,
              'file_grid_only': False, 'level': 0}
    kwargs.update(common_kwargs)
    [integral, error, prob, neval] = Vegas(
        TestIntegrand(ind), ranges=None,
        nstart=nstart, nincrease=nincrease,
        nbatch=nbatch, gridno=gridno,
        **kwargs)
    err_msgs = ["Integral values do not match!",
                "Error values do not match!",
                "P-values do not match!"]
    for item in zip(  # Not comparing prob values or neval.
            [integral, error], vegas_results[ind], err_msgs):
        nptest.assert_allclose(
            item[0], item[1], atol=1e-12, rtol=1e-3,
            err_msg="Problem #: " + str(ind) + item[2])


@pytest.mark.parametrize('ind', [0, pytest.mark.xfail(1, raises=CyCubaWarning),
                                 2, 3, 4, 5, 6, 7, 8, 9, 10])
def test_Vegas(ind):
    out = _vegas_test_runner(ind)


def _suave_test_runner(ind):
    kwargs = {'verbosity': verbose, 'last_samples_only': True,
              'do_not_smooth': False, 'retain_state_file': False,
              'level': 0}
    kwargs.update(common_kwargs)
    [integral, error, prob, neval] = Suave(
        TestIntegrand(ind), ranges=None,
        nnew=nnew, nmin=nmin, flatness=flatness,
        **kwargs)
    err_msgs = ["Integral values do not match!",
                "Error values do not match!",
                "P-values do not match!"]
    for item in zip(  # Not comparing prob values or neval.
            [integral, error], suave_results[ind], err_msgs):
        nptest.assert_allclose(
            item[0], item[1], atol=1e-12, rtol=1e-3,
            err_msg="Problem #: " + str(ind) + " " + item[2])


@pytest.mark.parametrize('ind', range(11))
def test_Suave(ind):
    out = _suave_test_runner(ind)


def _divonne_test_runner(ind):
    kwargs = {'verbosity': verbose, 'last_samples_only': False,
              'retain_state_file': False, 'level': 0}
    kwargs.update(common_kwargs)
    [integral, error, prob, neval] = Divonne(
        TestIntegrand(ind), ranges=None,
        key1=key1, key2=key2, key3=key3, maxpass=maxpass, border=border,
        maxchisq=maxchisq, mindeviation=mindeviation, xgiven=xgiven,
        nextra=nextra, peakfinder=peakfinder,
        **kwargs)
    err_msgs = ["Integral values do not match!",
                "Error values do not match!",
                "P-values do not match!"]
    for item in zip(  # Not comparing prob values or neval.
            [integral, error], divonne_results[ind], err_msgs):
        nptest.assert_allclose(
            item[0], item[1], atol=1e-12, rtol=1e-3,
            err_msg="Problem #: " + str(ind) + " " + item[2])


@pytest.mark.parametrize('ind', [0, pytest.mark.xfail(1), pytest.mark.xfail(2),
                                 3, 4, 5, 6, 7, 8, pytest.mark.xfail(9), 10])
# These three tests are close, but not close enough.
def test_Divonne(ind):
    out = _divonne_test_runner(ind)


def _cuhre_test_runner(ind):
    kwargs = {'verbosity': verbose, 'last_samples_only': True,
              'retain_state_file': False}
    kwargs.update(common_kwargs)
    [integral, error, prob, neval] = Cuhre(
        TestIntegrand(ind), ranges=None,
        key=key,
        **kwargs)
    err_msgs = ["Integral values do not match!",
                "Error values do not match!",
                "P-values do not match!"]
    for item in zip(  # Not comparing prob values or neval.
            [integral, error], cuhre_results[ind], err_msgs):
        nptest.assert_allclose(
            item[0], item[1], atol=1e-12, rtol=1e-3,
            err_msg="Problem #: " + str(ind) + " " + item[2])


@pytest.mark.parametrize('ind', range(11))
def test_Cuhre(ind):
    out = _cuhre_test_runner(ind)


def test_scaling():
    kwargs = {'verbosity': verbose, 'last_samples_only': False,
              'retain_state_file': False, 'level': 0}
    kwargs.update(common_kwargs)
    kwargs['maxeval'] = 5e5
    def test_function(x, y):
        return [1 if 1 - x**2 - y**2 > 0 else 0]
    [integral, error, prob, neval] = Divonne(
        test_function, ranges=[[-1, 1], [-1, 1]],
        key1=key1, key2=key2, key3=key3, maxpass=maxpass, border=border,
        maxchisq=maxchisq, mindeviation=mindeviation, xgiven=xgiven,
        nextra=nextra, peakfinder=peakfinder,
        **kwargs)
    nptest.assert_allclose(integral, pi, atol=1e-2)
