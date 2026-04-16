"""Generate flat tables showing Array API capabilities for use in docs.

These tables are intended for presenting Array API capabilities across
a wide number of functions at once. Rows correspond to functions and
columns correspond to library/device/option combinations.
"""

from collections import defaultdict
from enum import auto, Enum
from importlib import import_module
from types import ModuleType

from scipy._lib._array_api import xp_capabilities_table
from scipy._lib._array_api import _make_sphinx_capabilities

# For undocumented aliases of public functions which are kept around for
# backwards compatibility reasons. These should be excluded from the
# tables since they would be redundant. There are also no docs pages to
# link entries to.
ALIASES = {
    "scipy.linalg": {
        # Alias of scipy.linalg.solve_continuous_lyapunov
        "solve_lyapunov",
    },
    "scipy.ndimage": {
        # Alias of scipy.ndimage.sum_labels
        "sum",
    },
    "scipy.special": {
        # Alias of scipy.special.jv
        "jn",
        # Alias of scipy.special.roots_legendre
        "p_roots",
        # Alias of scipy.special.roots_chebyt
        "t_roots",
        # Alias of scipy.special.roots_chebyu
        "u_roots",
        # Alias of scipy.special.roots_chebyc
        "c_roots",
        # Alias of scipy.special.roots_chebys
        "s_roots",
        # Alias of scipy.special.roots_jacobi
        "j_roots",
        # Alias of scipy.special.roots_laguerre
        "l_roots",
        # Alias of scipy.special.roots_genlaguerre
        "la_roots",
        # Alias of scipy.special.roots_hermite
        "h_roots",
        # Alias of scipy.special.roots_hermitenorm
        "he_roots",
        # Alias of scipy.special.roots_gegenbauer
        "cg_roots",
        # Alias of scipy.special.roots_sh_legendre
        "ps_roots",
        # Alias of scipy.special.roots_sh_chebyt
        "ts_roots",
        # Alias of scipy.special.roots_chebyu
        "us_roots",
        # Alias of scipy.special.roots_sh_jacobi
        "js_roots",
    }
}

# Shortened names for use in table.
BACKEND_NAMES_MAP = {
    "jax.numpy": "jax",
    "dask.array": "dask",
}


class BackendSupportStatus(Enum):
    YES = auto()
    NO = auto()
    OUT_OF_SCOPE = auto()
    UNKNOWN = auto()


def _process_capabilities_table_entry(entry: dict | None) -> dict[str, dict[str, bool]]:
    """Returns dict showing alternative backend support in easy to consume form.

    Parameters
    ----------
    entry : Optional[dict]
       A dict with the structure of the values of the dict
       scipy._lib._array_api.xp_capabilities_table. If None, it is
       assumped that no alternative backends are supported.
       Default: None.

    Returns
    -------
    dict[str, dict[str, bool]]
        The output dict currently has keys "cpu", "gpu", "jit" and "lazy".
        The value associated to each key is itself a dict. The keys of
        the inner dicts correspond to backends, with bool values stating
        whether or not the backend is supported with a given device or
        mode. Inapplicable backends do not appear in the inner dicts
        (e.g. since cupy is gpu-only, it does not appear in the inner
        dict keyed on "cpu"). Only alternative backends to NumPy are
        included since NumPY support should be guaranteed.

    """
    # This is a template for the output format. If more backends and
    # backend options are added, it will need to be updated manually.
    # Entries start as boolean, but upon returning, will take values
    # from the BackendSupportStatus Enum.
    output = {
        "cpu": {"torch": False, "jax": False, "dask": False},
        "gpu": {"cupy": False, "torch": False, "jax": False},
        "jit": {"jax": False},
        "lazy": {"dask": False},
    }
    S = BackendSupportStatus
    if entry is None:
        # If there is no entry, assume no alternative backends are supported.
        # If the list of supported backends will grows, this hard-coded dict
        # will need to be updated.
        return {
            outer_key: {inner_key: S.UNKNOWN for inner_key in outer_value}
            for outer_key, outer_value in output.items()
        }

    if entry["out_of_scope"]:
        # None is used to signify out-of-scope functions.
        return {
            outer_key: {inner_key: S.OUT_OF_SCOPE for inner_key in outer_value}
            for outer_key, outer_value in output.items()
        }

    # For now, use _make_sphinx_capabilities because that's where
    # the relevant logic for determining what is and isn't
    # supported based on xp_capabilities_table entries lives.
    # Perhaps this logic should be decoupled from sphinx.
    for backend, capabilities in _make_sphinx_capabilities(**entry).items():
        if backend in {"array_api_strict", "numpy"}:
            continue
        backend = BACKEND_NAMES_MAP.get(backend, backend)
        cpu, gpu = capabilities.cpu, capabilities.gpu
        if cpu is not None:
            if backend not in output["cpu"]:
                raise ValueError(
                    "Input capabilities table entry contains unhandled"
                    f" backend {backend} on cpu."
                )
            output["cpu"][backend] = cpu
        if gpu is not None:
            if backend not in output["gpu"]:
                raise ValueError(
                    "Input capabilities table entry contains unhandled"
                    f" backend {backend} on gpu."
                )
            output["gpu"][backend] = gpu
        if backend == "jax":
            output["jit"]["jax"] = entry["jax_jit"] and output["cpu"]["jax"]
        if backend == "dask.array":
            support_lazy = not entry["allow_dask_compute"] and output["dask"]
            output["lazy"]["dask"] = support_lazy
    return {
        outer_key: {
            inner_key: S.YES if inner_value else S.NO
            for inner_key, inner_value in outer_value.items()
        }
        for outer_key, outer_value in output.items()
    }


def is_named_function_like_object(obj):
    return (
        not isinstance(obj, ModuleType | type)
        and callable(obj) and hasattr(obj, "__name__")
    )


def make_flat_capabilities_table(
        modules: str | list[str],
        backend_type: str,
        /,
        *,
        capabilities_table: list[str] | None = None,
) -> list[dict[str, str]]:
    """Generate full table of array api capabilities across public functions.

    Parameters
    ----------
    modules : str | list[str]
        A string containing single SciPy module, (e.g `scipy.stats`, `scipy.fft`)
        or a list of such strings.

    backend_type : {'cpu', 'gpu', 'jit', 'lazy'}

    capabilities_table : Optional[list[str]]
        Table in the form of `scipy._lib._array_api.xp_capabilities_table`.
        If None, uses `scipy._lib._array_api.xp_capabilities_table`.
        Default: None.

    Returns
    -------
    output : list[dict[str, str]]
        `output` is a table in dict format
        (keys corresponding to column names). The first column is "module".
        The other columns correspond to supported backends for the given
        `backend_type`, e.g. jax.numpy, torch, and dask on cpu.
         numpy is excluded because it should always be supported.
         See the helper function
        `_process_capabilities_table_entry` above).

    """
    if backend_type not in {"cpu", "gpu", "jit", "lazy"}:
        raise ValueError(f"Received unhandled backend type {backend_type}")

    if isinstance(modules, str):
        modules = [modules]

    if capabilities_table is None:
        capabilities_table = xp_capabilities_table

    output = []

    for module_name in modules:
        module = import_module(module_name)
        public_things = module.__all__
        for name in public_things:
            if name in ALIASES.get(module_name, {}):
                # Skip undocumented aliases that are kept
                # for backwards compatibility reasons.
                continue
            thing = getattr(module, name)
            if not is_named_function_like_object(thing):
                continue
            entry = xp_capabilities_table.get(thing, None)
            capabilities = _process_capabilities_table_entry(entry)[backend_type]
            row = {"module": module_name}
            row.update({"function": name})
            row.update(capabilities)
            output.append(row)
    return output


def calculate_table_statistics(
    flat_table: list[dict[str, str]]
) -> dict[str, tuple[dict[str, str], bool]]:
    """Get counts of what is supported per module.

    Parameters
    ----------
    flat_table : list[dict[str, str]]
        A table as returned by `make_flat_capabilities_table`

    Returns
    -------
    dict[str, tuple[dict[str, str], bool]]
        dict mapping module names to 2-tuples containing an inner dict and a
        bool. The inner dicts have a key "total" along with keys for each
        backend column of the supplied flat capabilities table. The value
        corresponding to total is the total count of functions in the given
        module, and the value associated to the other keys is the count of
        functions that support that particular backend. The bool is False if
        the calculation may be innacurate due to missing xp_capabilities
        decorators, and True if all functions for that particular module have
        been decorated with xp_capabilities.
    """
    if not flat_table:
        return []

    counter = defaultdict(lambda: defaultdict(int))

    S = BackendSupportStatus
    # Keep track of which modules have functions with missing xp_capabilities
    # decorators so this information can be passed back to the caller.
    missing_xp_capabilities = set()
    for entry in flat_table:
        entry = entry.copy()
        entry.pop("function")
        module = entry.pop("module")
        current_counter = counter[module]

        # By design, all backends and options must be considered out-of-scope
        # if one is, so just pick an arbitrary entry here to test if function is
        # in-scope.
        if next(iter(entry.values())) != S.OUT_OF_SCOPE:
            current_counter["total"] += 1
            for key, value in entry.items():
                # Functions missing xp_capabilities will be tabulated as
                # unsupported, but may actually be supported. There is a
                # note about this in the documentation and this function is
                # set up to return information needed to put asterisks next
                # to percentages impacted by missing xp_capabilities decorators.
                current_counter[key] += 1 if value == S.YES else 0
                if value == S.UNKNOWN:
                    missing_xp_capabilities.add(module)
    return {
        key: (dict(value), key not in missing_xp_capabilities)
        for key, value in counter.items()
    }
