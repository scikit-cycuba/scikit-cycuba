# file: cycuba.pyx
from cpython.mem cimport PyMem_Malloc, PyMem_Free
cimport ccuba
from ccuba cimport cubareal

#cdef int _integrand_c(int ndim, int ncomp, integrand_py):
#    return ccuba.integrand_t integrand(integrand_py)

#def void _vegas(ndim, ncomp, ):
#
#    ccuba.Vegas(ndim, ncomp, integrand, userdata, nvec, epsrel, epsabs,
#                flags, seed, mineval, maxeval, nstart, nincrease, nbatch,
#                gridno, statefile, spin, neval, fail, integral, error,
#                prob)

cdef:
    # Common arguments
    int ndim_c
    int ncomp_c
    void* userdata_c
    int nvec_c
    ccuba.cubareal epsrel_c
    ccuba.cubareal epsabs_c
    int flags_c
    int seed_c
    int mineval_c
    int maxeval_c
    char* statefile_c
    void* spin_c

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
    double xgiven_c
    int nextra_c

    #Cuhre-specific arguments
    int key_c

    #Common returns
    int nregions_value = 0 # Not used by Vegas
    int* nregions_c = &nregions_value
    int neval_value = 0
    int* neval_c = &neval_value
    int fail_value = 0
    int* fail_c = &fail_value
    

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
    statefile_c = NULL #&statefile.encode('utf-8')
    spin_c = NULL # Not sure what to do here.
        
    
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
        raise Exception("Suave integration not yet implemented!")
    elif integrator == 'divonne':
        raise Exception("Divonne integration not yet implemented!")
    elif integrator == 'cuhre':
        key_c = <int> key

        ccuba.Cuhre(
            ndim_c, ncomp_c, _integrand_c, userdata_c, nvec_c, epsrel_c,
            epsabs_c, flags_c, mineval_c, maxeval_c, key_c, statefile_c,
            spin_c, nregions_c, neval_c, fail_c, integral_c, error_c, prob_c)
    else:
        raise Exception("Bad value for integrator")
    
    if integrator in ['suave', 'divonne', 'cuhre']:
        out.append(nregions_c[0])
    integral = [integral_c[i] for i in range(ncomp)]
    error = [error_c[i] for i in range(ncomp)]
    prob = [prob_c[i] for i in range(ncomp)]
    PyMem_Free(integral_c)
    PyMem_Free(error_c)
    PyMem_Free(prob_c)
    
    out.extend([neval_c[0], fail_c[0], integral, error, prob])
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
        return <int>-999
    return 0

#cdef void _peakfinder_c(const int *ndim, const double b[],
#                        int *n, double f[], void *userdata):
#    try:
#        func = <object> userdata
#        x_py = [<double> entry for entry in x[:ndim[0]]]
#        b_py = [[

def integer_bit_flags(verbosity=3, all_samples=True,
                      vegas_suave_smoothing=False, delete_state_file=True,
                      vegas_use_file_state=True, level=0):
    level_key = 3 if level > 4 and level < 24 else level
    flag_string = (
        format(level_key, 'b')
        + {True: '0', False: '1'}[vegas_use_file_state]
        + {True: '0', False: '1'}[delete_state_file]
        + {True: '0', False: '1'}[vegas_suave_smoothing]
        + {True: '0', False: '1'}[all_samples]
        + {0: '00', 1: '01', 2: '10', 3: '11'}[verbosity]
    )
    return int(flag_string,base=2)


def test_run():
    def test_function(x, y):
        return [1 if x**2 + y**2 < 1 else 0]
    args = ('cuhre', # integrator
            #Common arguments
            2, # ndim
            1, # ncomp
            test_function, # integrand
            1, # nvec
            1e-3, # epsrel
            1e-12, # epsabs
            11, # flags
            0, # seed
            1e2, # mineval
            1e5, # maxeval
            "", # statefile
            1, # spin
            # Vegas-specific
            1000, # nstart
            500, # nincrease
            1000, # nbatch
            0, # gridno
            # Suave-specific
            0, # nnew
            0, # nmin
            0, # flatness
            #Divonne-specific
            0, # key1
            0, # key2
            0, # key3
            0, # maxpass
            0, # border
            0, # maxchisq
            0, # mindeviation
            0, # ngiven
            0, # ldxgiven
            0, # xgiven
            0, # nextra
            0, # peakfinder
            # Cuhre-specific
            0 # key
    )
    return _cuba(*args)


    
