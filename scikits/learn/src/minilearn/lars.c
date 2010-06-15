/*
 * Least-angle regression (LARS) is a regression algorithm for
 * high-dimensional data, developed by Bradley Efron, Trevor Hastie,
 * Iain Johnstone and Robert Tibshirani.
 *
 * The advantages of LARS are:
 *
 * - It is computationally just as fast as forward selection and has
 *   the same order of complexity as an ordinary least squares.
 *
 * - It produces a full piecewise linear solution path, which is
 *   useful in cross-validation or similar attempts to tune the model.
 *
 * - If two variables are almost equally correlated with the response,
 *   then their coefficients should increase at approximately the same
 *   rate. The algorithm thus behaves as intuition would expect, and
 *   also is more stable.
 *
 * - It is easily modified to produce solutions for other estimators,
 *   like the Lasso. We implement lasso in method lars_lasso_TODO
 *
 * - It is effective in contexts where p >> n (IE, when the number of
 *   dimensions is significantly greater than the number of points)
 *
 * The disadvantages of the LARS method include:
 *
 * - Because LARS is based upon an iterative refitting of the
 *   residuals, it would appear to be especially sensitive to the
 *   effects of noise. This problem is discussed in detail by Weisberg
 *   in the discussion section of the Efron et al. (2004) Annals of
 *   Statistics article.
 * 
 *
 * Dependencies
 * ------------
 * cblas
 * 
 * Author
 * ------
 * Fabian Pedregosa <fabian.pedregosa@inria.fr>.
 *
 * License
 * -------
 * Public Domain.
 *
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <cblas.h>
#include "minilearn.h"

/* 
 * Doubly linked list structure, we will use this to iterate over
 * indices.
 */
struct dllist {
    double *ptr;
    double cov; /* current covariance */
    struct dllist *next;
    struct dllist *prev;
};


/*
 * Fit LARS algorithm.
 *
 * Parameters
 * ----------
 * predictors are expected to be normalized
 *
 *  X is (features, nsamples), b is (nsamples,1) 
 *
 * beta is the returned array with the parameters of the regression
 * problem.
 * 
 * b will be modified to store the current residual.
 *
 * L will hold the cholesky decomposition.
 * 
 * TODO
 * ----
 * The case when two vectors have equal correlation.
 *
 * Notes
 * -----
 * Most iterations need only be done with indices in the active set
 * (or in its complement), this limits the use of BLAS leves 2, 3.
 *
 * TODO
 * ----
 * Consider swapping columns/vector in the problem formulation, it
 * should make vector multiplication faster.
 *
 */
void lars_fit(int nfeatures, int nsamples, double *X, double *res, 
              double *beta, double *lambdas, int *row, int *col, 
              double *L, int niter)
{
    /* temp variables */
    /* TODO: pass uu by reference */
    double *uu  = (double *) calloc(nsamples, sizeof(double));

    double *v, *dir, C, temp, gamma=0, tgamma, aj, Aa = 0;
    int i, k, sum_k=0;

    struct dllist *active_set, *head, *cur, *top_active;
    struct dllist *pmax=NULL; /* current working variable */

    active_set = (struct dllist *) malloc((nfeatures+1) * sizeof(struct dllist));

    /* create index list as a circular doubly linked list */
    for (i=0; i<=nfeatures; ++i) {
        active_set[i].ptr = X + i;
        active_set[i].next = active_set + i + 1;
        active_set[i].prev = active_set + i - 1;
    }

    head = active_set + nfeatures;
    head->next = active_set;
    active_set->prev = head;

    /* set sentinels */
    top_active = NULL;
    head->ptr = NULL;

    /* main loop, we iterate over the user-suplied number of iterations */
    for (k=0; k < niter; ++k) {

        sum_k = k * (k+1) / 2;
        v = L + sum_k;
        dir = beta + sum_k;

        /* 
         * Update residual.
         */
        cblas_daxpy (nsamples, -gamma, uu, 1, res, 1);

        /* calculate covariances (c_hat), and get maximum covariance (C_hat)*/
        for (pmax=head->next, cur=head->next; cur->ptr; cur=cur->next) {
            cur->cov = cblas_ddot (nsamples, cur->ptr, nfeatures, res, 1);
            pmax = (fabs(cur->cov) > fabs(pmax->cov)) ? cur : pmax;
        }

        /* remove it from the unused set */
        pmax->prev->next = pmax->next;
        pmax->next->prev = pmax->prev;

        /* push pmax into the active set */
        pmax->prev = top_active;
        top_active = pmax;

        /* 
         * Compute the least-squares direction of the coefficients.
         */
        dir[k] = temp = pmax->cov;

        for (i=k-1, cur=top_active->prev; cur; cur=cur->prev, --i) {
            /*
             * To update the cholesky decomposition, we need to compute
             *
             *                    v   = Xa' * cur
             *                    dir = Xa' * res
             *
             * where res is the current residual and cur is the last
             * vector to enter the active set.
             *
             * Covariances in the active set are kept tied and decreasing.
             */
            v[i] = cblas_ddot (nsamples, cur->ptr, nfeatures, 
                               pmax->ptr, nfeatures);

            temp = copysign(temp, cur->cov);
            dir[i] = temp;
        }

        /*
         * Compute least squares solution by the method of normal equations
         * (Golub & Van Loan, 1996)
         * 
         * Update the cholesky decomposition of (Xa * Xa')
         *
         *          ( L   0 )
         *   L  ->  (       )  , where L * w = v , z = ||u||
         *          ( w   z )
         *
         * and u is the last vector added to the active set.
         *
         */
        cblas_dtpsv (CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit, 
                     k, L, v, 1);

        v[k] = sqrt(1 - cblas_ddot(k, v, 1, v, 1));

        cblas_dtpsv (CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit,
                     k + 1, L, dir, 1);

        cblas_dtpsv (CblasRowMajor, CblasLower, CblasTrans, CblasNonUnit,
                     k + 1, L, dir, 1);

        /* 
         * Update uu with the current direction 
         * uu = Xa' * dir
         */
        cblas_dscal (nsamples, 0., uu, 1);
        for (i=k, cur = top_active; cur; cur = cur->prev, --i)
            cblas_daxpy (nsamples, dir[i], cur->ptr, nfeatures, uu, 1);
        Aa = cblas_ddot (nsamples, top_active->ptr, nfeatures, uu, 1);
        Aa = fabs (Aa)        + DBL_EPSILON;
        C  = fabs (pmax->cov) + DBL_EPSILON;

        /*
         * Compute gamma.
         */
        gamma = DBL_MAX;
        for (cur = head->next; cur->ptr; cur = cur->next) {
            aj = cblas_ddot (nsamples, cur->ptr, nfeatures, uu, 1);
            if (cur->cov > 0) {
                tgamma = (C - cur->cov) / (Aa - aj);
            }
            else {
                tgamma = (C + cur->cov) / (Aa + aj);
            }
            gamma = fmin (tgamma, gamma);
        }

        /* 
         * Set up return values.
         * TODO: we should iterate over used_indices (for Lasso variant).
         */
        cblas_daxpy (k,  1./gamma, dir-k, 1, dir, 1);
        cblas_dscal (k + 1, gamma, dir, 1);
        lambdas[k] = pmax->cov;

        /* TODO: this is proper of LAR */
        memcpy (row + sum_k, row + sum_k - k, k * sizeof(int));
        row[sum_k + k] = (pmax->ptr - X);
        for (i=0; i<=k; ++i) col[sum_k + i] = k;
    }
 
    free (active_set);
    free (uu);
}
