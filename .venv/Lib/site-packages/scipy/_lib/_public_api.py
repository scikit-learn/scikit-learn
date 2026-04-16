"""PUBLIC_MODULES was once included in scipy._lib.tests.test_public_api.

It has been separated into this file so that this list of public modules
could be used when generating tables showing support for alternative
array API backends across modules in
scipy/doc/source/array_api_capabilities.py.
"""

# Historically SciPy has not used leading underscores for private submodules
# much.  This has resulted in lots of things that look like public modules
# (i.e. things that can be imported as `import scipy.somesubmodule.somefile`),
# but were never intended to be public.  The PUBLIC_MODULES list contains
# modules that are either public because they were meant to be, or because they
# contain public functions/objects that aren't present in any other namespace
# for whatever reason and therefore should be treated as public.
PUBLIC_MODULES = ["scipy." + s for s in [
    "cluster",
    "cluster.vq",
    "cluster.hierarchy",
    "constants",
    "datasets",
    "differentiate",
    "fft",
    "fftpack",
    "integrate",
    "interpolate",
    "io",
    "io.arff",
    "io.matlab",
    "io.wavfile",
    "linalg",
    "linalg.blas",
    "linalg.cython_blas",
    "linalg.lapack",
    "linalg.cython_lapack",
    "linalg.interpolative",
    "ndimage",
    "odr",
    "optimize",
    "optimize.elementwise",
    "signal",
    "signal.windows",
    "sparse",
    "sparse.linalg",
    "sparse.csgraph",
    "spatial",
    "spatial.distance",
    "spatial.transform",
    "special",
    "stats",
    "stats.contingency",
    "stats.distributions",
    "stats.mstats",
    "stats.qmc",
    "stats.sampling"
]]
