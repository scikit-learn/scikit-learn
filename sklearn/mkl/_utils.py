# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause


def number_of_kernels(
    X,
    kernels,
    kernels_scopes,
    kernels_param_grids,
    precomputed_kernels,
):
    if precomputed_kernels or kernels is None:
        return len(X)

    number = 0
    for scope, params in zip(kernels_scopes, kernels_param_grids):
        for _ in range(len(params[next(iter(params))]) if len(params) > 0 else 1):
            if scope == "all":
                number += 1
            else:
                number += X.shape[1]

    return number


def kernel_generator(
    X,
    kernels,
    kernels_scopes,
    kernels_param_grids,
    precomputed_kernels,
    Y=None,
):
    if precomputed_kernels or kernels is None:
        for kernel in X:
            yield kernel
    else:
        if Y is None:
            Y = X
        for kernel, scope, params in zip(kernels, kernels_scopes, kernels_param_grids):
            for i in range(len(params[next(iter(params))]) if len(params) > 0 else 1):
                if scope == "all":
                    yield kernel(X, Y, **{k: v[i] for k, v in params.items()})
                else:
                    for j in range(X.shape[1]):
                        yield kernel(
                            X[:, j].reshape(-1, 1),
                            Y[:, j].reshape(-1, 1),
                            **{k: v[i] for k, v in params.items()},
                        )
