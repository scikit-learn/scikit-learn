# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause


def number_of_kernels(X, kernels, kernels_params):
    if kernels == "precomputed":
        return len(X)

    number = 0
    for kernel, kernel_params in zip(kernels, kernels_params):
        dct = kernel_params[1]
        for i in range(len(dct[next(iter(dct))]) if len(dct) > 0 else 1):
            if kernel_params[0] == "all":
                number += 1
            else:
                number += X.shape[1]

    return number


def kernel_generator(X, kernels, kernels_params):
    if kernels == "precomputed":
        for kernel in X:
            yield kernel
    else:
        for kernel, kernel_params in zip(kernels, kernels_params):
            dct = kernel_params[1]
            for i in range(len(dct[next(iter(dct))]) if len(dct) > 0 else 1):
                if kernel_params[0] == "all":
                    yield kernel(X, **{k: v[i] for k, v in dct.items()})
                else:
                    for j in range(X.shape[1]):
                        yield kernel(
                            X[:, j].reshape(-1, 1), **{k: v[i] for k, v in dct.items()}
                        )
