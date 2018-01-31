import numpy as np


def _two_point_interp(x0, x1, y0, y1, y):
    """Finds x such that the point (x, y) lies on the line through (x0, y0)
    and (x1, y1). Calling this function should be slightly faster than
    np.interp(y, [y0, y1], [x0, x1]).
    """
    dy = y1 - y0
    if dy == 0.:
        raise ValueError("Line is horizontal; can't interpolate for x.")
    dx = x1 - x0
    inv_slope = dx / dy
    return x0 + (y - y0) * inv_slope


def prox_box(w, lambd, a=0., b=1., c=1., eps=1e-15, verbose=0):
    """Computes the proximal operator of the squared box-norm.
    with parameters a, b and c

        ||w||_box^2 := \supp_{\theta \in \Theta} \sqrt{w^T diag(1 / \theta) w}

    where

        \Theta := \{\theta | a < \theta_i \le b,\; \sum_i \theta_i \le c\}.

    The case a = 0, b = 1, c = k (integer) corresponds to the well-known
    k-support norm.

    The time-complexity of the algorithm is O(p log p).
    """
    # misc
    p = len(w)
    absw = np.abs(w)
    mask = absw > eps
    absw = absw[mask]
    w = w[mask]
    min_alpha = (lambd + a) / absw
    max_alpha = (lambd + b) / absw
    alphas = np.sort(np.concatenate((min_alpha, max_alpha)))

    def _left_condition(alpha):
        return alpha <= min_alpha

    def _middle_condition(alpha):
        return np.logical_and(min_alpha < alpha, alpha < max_alpha)

    def _right_condition(alpha):
        return alpha > max_alpha

    def _compute_s(alpha):
        s = 0.
        if b != 0.:
            s += b * _right_condition(alpha).sum()
        s += np.sum(alpha * absw[_middle_condition(alpha)] - lambd)
        if a != 0.:
            s += a * _left_condition(alpha).sum()
        return s

    # Binary-search for an index i such that s(alpha_i) <= c <= s(alpha_i_+_1).
    # The loops lasts for at most O(log p) iterations each of complexity s,
    # if w is s-sparse, for a total complexity of O(s log p).
    start = 0
    stop = 2 * p - 1
    hit = False
    curr = p
    it = 0
    while start < stop:
        it += 1

        # check s(alpha[start] <= c <= s(alphas[start + 1])
        if start + 1 == stop:
            break

        alpha = alphas[curr]
        s = _compute_s(alpha)
        if verbose:
            print(" start: % 7i stop: % 7i curr: % 7i alpha: % 7f s: % 7f" % (
                start, stop, curr, alpha, s))
        if s > c:
            stop = curr
        elif s < c:
            start = curr
        else:
            hit = True
            break

        # no, luck, half the interval and continue
        curr = (start + stop) // 2
    if verbose:
        print("Loop exited after %i iterations." % it)

    # 1D interpolation: runs in constant time O(1)
    if not hit:
        low_alpha = alphas[start]
        high_alpha = alphas[stop]
        low_s = _compute_s(low_alpha)
        high_s = _compute_s(high_alpha)
        alpha = _two_point_interp(low_alpha, high_alpha, low_s, high_s, c)

    # compute the prox proper
    theta = np.full(p, a, dtype=w.dtype)
    theta_masked = theta[mask]
    theta_masked[_right_condition(alpha)] = b
    middle_condition = _middle_condition(alpha)
    theta_masked[middle_condition] = alpha * absw[middle_condition] - lambd
    theta_masked = (theta_masked * w) / (theta_masked + lambd)
    theta[mask] = theta_masked
    return theta


def prox_ksupp(w, lambd, k=1, **kwargs):
    return prox_box(w, lambd, a=0., b=1., c=k, **kwargs)
