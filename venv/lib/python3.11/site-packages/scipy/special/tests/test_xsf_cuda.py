import os
import pytest
import scipy.special as sc
import shutil
import tempfile

from uuid import uuid4

from scipy.special._testutils import check_version
from scipy.special._testutils import MissingModule

try:
    import cupy  # type: ignore
except (ImportError, AttributeError):
    cupy = MissingModule('cupy')


def get_test_cases():
    cases_source = [
        (sc.beta, "cephes/beta.h", "out0 = xsf::cephes::beta(in0, in1)"),
        (sc.binom, "binom.h", "out0 = xsf::binom(in0, in1)"),
        (sc.digamma, "digamma.h", "xsf::digamma(in0)"),
        (sc.expn, "cephes/expn.h", "out0 = xsf::cephes::expn(in0, in1)"),
        (sc.hyp2f1, "hyp2f1.h", "out0 = xsf::hyp2f1(in0, in1, in2, in3)"),
        (sc._ufuncs._lambertw, "lambertw.h", "out0 = xsf::lambertw(in0, in1, in2)"),
        (sc.ellipkinc, "cephes/ellik.h", "out0 = xsf::cephes::ellik(in0, in1)"),
        (sc.ellipeinc, "cephes/ellie.h", "out0 = xsf::cephes::ellie(in0, in1)"),
        (sc.gdtrib, "cdflib.h", "out0 = xsf::gdtrib(in0, in1, in2)"),
        (sc.sici, "sici.h", "xsf::sici(in0, &out0, &out1)"),
        (sc.shichi, "sici.h", "xsf::shichi(in0, &out0, &out1)"),
    ]

    cases = []
    for ufunc, header, routine in cases_source:
        preamble = f"#include <xsf/{header}>"
        for signature in ufunc.types:
            cases.append((signature, preamble, routine))
    return cases


dtype_map = {
    "f": "float32",
    "d": "float64",
    "F": "complex64",
    "D": "complex128",
    "i": "int32",
    "l": "int64",
}


def get_params(signature):
    in_, out = signature.split("->")
    in_params = []
    out_params = []
    for i, typecode in enumerate(in_):
        in_params.append(f"{dtype_map[typecode]} in{i}")
    for i, typecode in enumerate(out):
        out_params.append(f"{dtype_map[typecode]} out{i}")
    in_params = ", ".join(in_params)
    out_params = ", ".join(out_params)
    return in_params, out_params


def get_sample_input(signature, xp):
    dtype_map = {
        "f": xp.float32,
        "d": xp.float64,
        "F": xp.complex64,
        "D": xp.complex128,
        "i": xp.int32,
        "l": xp.int64,
    }

    in_, _ = signature.split("->")
    args = []
    for typecode in in_:
        args.append(xp.zeros(2, dtype=dtype_map[typecode]))
    return args


@pytest.fixture(scope="module", autouse=True)
def manage_cupy_cache():
    # Temporarily change cupy kernel cache location so kernel cache will not be polluted
    # by these tests. Remove temporary cache in teardown.
    temp_cache_dir = tempfile.mkdtemp()
    original_cache_dir = os.environ.get('CUPY_CACHE_DIR', None)
    os.environ['CUPY_CACHE_DIR'] = temp_cache_dir

    yield

    if original_cache_dir is not None:
        os.environ['CUPY_CACHE_DIR'] = original_cache_dir
    else:
        del os.environ['CUPY_CACHE_DIR']
    shutil.rmtree(temp_cache_dir)


@check_version(cupy, "13.0.0")
@pytest.mark.parametrize("signature,preamble,routine", get_test_cases())
@pytest.mark.xslow
def test_compiles_in_cupy(signature, preamble, routine, manage_cupy_cache):
    name = f"x{uuid4().hex}"
    in_params, out_params = get_params(signature)

    func = cupy.ElementwiseKernel(
        in_params,
        out_params,
        routine,
        name,
        preamble=preamble,
        options=(f"--include-path={sc._get_include()}", "-std=c++17")
    )

    _ = func(*get_sample_input(signature, cupy))
