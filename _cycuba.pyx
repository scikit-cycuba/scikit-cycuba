# file: cycuba.pyx
from cpython.mem cimport PyMem_Malloc, PyMem_Free

from ccuba cimport cubareal
cimport ccuba

cdef:
    # Common arguments
    int ndim_c
    int ncomp_c
    void *userdata_c
    int nvec_c
    ccuba.cubareal epsrel_c
    ccuba.cubareal epsabs_c
    int flags_c
    int seed_c
    int mineval_c
    int maxeval_c
    char *statefile_c
    void *spin_c

    # Vegas-specific arguments
    int nstart_c
    int nincrease_c
    int nbatch_c
    int gridno_c

    # Suave-specific arguments
    int nnew_c
    int nmin_c
    double flatness_c

    # Divonne-specific arguments
    int key1_c
    int key2_c
    int key3_c
    int maxpass_c
    double border_c
    double maxchisq_c
    double mindeviation_c
    int ngiven_c
    int ldxgiven_c
    int nextra_c

    # Cuhre-specific arguments
    int key_c

    # Common returns
    int nregions_value = 0  # Not used by Vegas
    int *nregions_c = &nregions_value
    int neval_value = 0
    int *neval_c = &neval_value
    int fail_value = 0
    int *fail_c = &fail_value

    # Peakfinder declarations, required for Divonne
    int peakfinder_n_value
    int *peakfinder_n_c


def _cuba(integrator, ndim, ncomp, integrand, nvec, epsrel, epsabs,
          flags, seed, mineval, maxeval, statefile, spin, nstart, nincrease,
          nbatch, gridno, nnew, nmin, flatness, key1, key2, key3, maxpass,
          border, maxchisq, mindeviation, ngiven, ldxgiven, xgiven, nextra,
          peakfinder, key):
    """Convert all necessary arguments to their C equivalents and call the
appropriate Cuba routine.

Integrator values:
 - 'vegas'
 - 'suave'
 - 'divonne'
 - 'cuhre'

Function signatures:

int (*integrand_t)(const int* ndim, const double x[],
                   const int* ncomp, double f[], void* userdata)

void (*peakfinder_t)(const int* ndim, const double b[],
                     int* n, double x[], void* userdata)

"""
    # Common arguments
    ndim_c = <int> ndim
    ncomp_c = <int> ncomp
    userdata_c = <void*>integrand
    nvec_c = <int> nvec
    epsrel_c = <double> epsrel
    epsabs_c = <double> epsabs
    flags_c = <int> flags
    seed_c = <int> seed
    mineval_c = <int> mineval
    maxeval_c = <int> maxeval
    if statefile:
        raise Exception("Statefile recording not yet supported!")
    statefile_c = NULL  # &statefile.encode('utf-8')
    if spin:
        raise Exception("Spin specification not yet supported!")
    spin_c = NULL  # Not sure what to do here.

    # Common returns
    integral_c = <ccuba.cubareal*> PyMem_Malloc(ncomp * sizeof(ccuba.cubareal))
    error_c = <ccuba.cubareal*> PyMem_Malloc(ncomp * sizeof(ccuba.cubareal))
    prob_c = <ccuba.cubareal*> PyMem_Malloc(ncomp * sizeof(ccuba.cubareal))

    out = []

    if integrator == 'vegas':
        nstart_c = <int> nstart
        nincrease_c = <int> nincrease
        nbatch_c = <int> nbatch
        gridno_c = <int> gridno

        ccuba.Vegas(
            ndim_c, ncomp_c, _integrand_c, userdata_c, nvec_c, epsrel_c,
            epsabs_c, flags_c, seed_c, mineval_c, maxeval_c, nstart_c,
            nincrease_c, nbatch_c, gridno_c, statefile_c, spin_c, neval_c,
            fail_c, integral_c, error_c, prob_c)

    elif integrator == 'suave':
        nnew_c = <int> nnew
        nmin_c = <int> nmin
        flatness_c = <double> flatness

        ccuba.Suave(
            ndim_c, ncomp_c, _integrand_c, userdata_c, nvec_c, epsrel_c,
            epsabs_c, flags_c, seed_c, mineval_c, maxeval_c, nnew_c, nmin_c,
            flatness_c, statefile_c, spin_c, nregions_c, neval_c, fail_c,
            integral_c, error_c, prob_c)

    elif integrator == 'divonne':
        key1_c = <int> key1
        key2_c = <int> key2
        key3_c = <int> key3
        maxpass_c = <int> maxpass
        border_c = <double> border
        maxchisq_c = <double> maxchisq
        mindeviation_c = <double> mindeviation
        ngiven_c = <int> ngiven
        ldxgiven_c = <int> ldxgiven
        xgiven_c = <double*> PyMem_Malloc(
            ldxgiven_c * ngiven_c * sizeof(double))
        for inda in range(ngiven):
            for indb in range(ldxgiven):
                xgiven_c[indb + ldxgiven * inda] = xgiven[indb][inda]
        if peakfinder:
            raise Exception(
                "The use of peakfinder functions in CyCuba's Divonne " +
                "integrator has not been fully implemented. If you need " +
                "peakfinder, then please contact the developers " +
                "in order to help them set up test problems and APIs.")
            nextra_c = <int> nextra
            _peakfinder_c = _peakfinder_c_true
        else:
            nextra_c = <int> 0

        peakfinder_b_c = <double*> PyMem_Malloc(ndim * 2 * sizeof(double))

        ccuba.Divonne(
            ndim_c, ncomp_c, _integrand_c, userdata_c, nvec_c, epsrel_c,
            epsabs_c, flags_c, seed_c, mineval_c, maxeval_c, key1_c,
            key2_c, key3_c, maxpass_c, border_c, maxchisq_c, mindeviation_c,
            ngiven_c, ldxgiven_c, xgiven_c, nextra_c, _peakfinder_c,
            statefile_c, spin_c, nregions_c, neval_c, fail_c, integral_c,
            error_c, prob_c)
        PyMem_Free(xgiven_c)
        PyMem_Free(peakfinder_b_c)
    elif integrator == 'cuhre':
        key_c = <int> key

        ccuba.Cuhre(
            ndim_c, ncomp_c, _integrand_c, userdata_c, nvec_c, epsrel_c,
            epsabs_c, flags_c, mineval_c, maxeval_c, key_c, statefile_c,
            spin_c, nregions_c, neval_c, fail_c, integral_c, error_c, prob_c)
    else:
        raise Exception("Bad value for integrator")

    if integrator in ['suave', 'divonne', 'cuhre']:
        out.append(<object> nregions_c[0])
    integral = [integral_c[i] for i in range(ncomp)]
    error = [error_c[i] for i in range(ncomp)]
    prob = [prob_c[i] for i in range(ncomp)]
    PyMem_Free(integral_c)
    PyMem_Free(error_c)
    PyMem_Free(prob_c)
    print(fail_c[0])
    out.extend([<object> neval_c[0], <object> fail_c[0], integral, error, prob])
    return out

cdef int _integrand_c(const int *ndim, const ccuba.cubareal x[],
                      const int *ncomp, ccuba.cubareal f[],
                      void *userdata):
    try:
        func = <object> userdata
        x_py = [<double> entry for entry in x[:ndim[0]]]
        out_py = func(*x_py)
        for ind in range(ncomp[0]):
            f[ind] = <ccuba.cubareal> <double> out_py[ind]
    except:
        return <int> -999
    return 0

cdef void _peakfinder_c(const int *ndim, const double b[],
                        int *n, double x[], void *userdata):
    obj = <object> userdata
    func = obj.peakfinder
    b_py = [[<double> b[ind * ndim[0]], <double> b[ind * ndim[0] + 1]]
            for ind in ndim[0]]
    n_py = n[0]
    out_py_n, out_py_x = func(b_py, n_py)
    for inda in range(out_py_n):
        for indb in range(ndim[0]):
            x[inda * ndim[0] + indb] = <double> out_py_x[inda][indb]
    n[0] = <int> out_py_n
