# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from itertools import product


def number_of_kernels(
    X,
    kernels,
    kernels_scopes,
    kernels_param_grids,
    precomputed_kernels,
):
    """
    Returns the number of kernels based on the input data and kernels configuration.
    """
    if precomputed_kernels or kernels is None:
        return len(X)

    number = 0
    for scope, params in zip(kernels_scopes, kernels_param_grids):
        if scope == "all":
            number += len(list(product(*params.values())))
        else:
            number += (X.shape[1] if len(X.shape) > 1 else 1) * len(
                list(product(*params.values()))
            )

    return number


def kernel_generator(
    X,
    kernels,
    kernels_scopes,
    kernels_param_grids,
    precomputed_kernels,
    Y=None,
):
    """
    Generates kernels based on the input data and kernels configuration, yielding
    one kernel at a time.
    """
    if precomputed_kernels or kernels is None:
        for kernel in X:
            yield kernel
    else:
        if Y is None:
            Y = X
        for kernel, scope, params in zip(kernels, kernels_scopes, kernels_param_grids):
            for params_prod in (
                dict(zip(params.keys(), values)) for values in product(*params.values())
            ):
                if scope == "all":
                    yield kernel(
                        X,
                        Y,
                        **params_prod,
                    )
                else:
                    if len(X.shape) == 1:
                        yield kernel(
                            X.reshape(-1, 1),
                            Y.reshape(-1, 1),
                            **params_prod,
                        )
                    else:
                        for j in range(X.shape[1]):
                            yield kernel(
                                X[:, j].reshape(-1, 1),
                                Y[:, j].reshape(-1, 1),
                                **params_prod,
                            )
