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

/* some compatibility for MS Visual C compiler */
#ifdef _MSC_VER
    #define copysign _copysign
    #define fmin __min
#endif

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
void lars_fit(int itype, int m, int nsamples, double *X, double *res, 
              double *beta, double *lambdas, int *row, int *col, 
              double *L, int niter)
{
    /* temp variables */
    /* TODO: pass uu by reference */
    double *uu  = (double *) calloc(nsamples, sizeof(double));

    double *v, *dir, C, temp, gamma, tgamma, aj, Aa = 0;
    double cj, ga, gb;
    int i, k=0, drop=0, sum_k=0;

    struct dllist *active_set, *head, *cur, *top_active;
    struct dllist *pmax=NULL; /* current working variable */

    active_set = (struct dllist *) malloc((m+1) * sizeof(struct dllist));

    /* create index list as a circular doubly linked list */
    for (i=0; i<=m; ++i) {
        active_set[i].ptr = X + i;
        active_set[i].next = active_set + i + 1;
        active_set[i].prev = active_set + i - 1;
    }

    head = active_set + m;
    head->next = active_set;
    active_set->prev = head;

    /* set sentinels */
    top_active = NULL;
    head->ptr = NULL;

    /* calculate covariances (c_hat), and get maximum covariance (C_hat)*/
    for (pmax=head->next, cur=head->next; cur->ptr; cur=cur->next) {
        cur->cov = cblas_ddot (nsamples, cur->ptr, m, res, 1);
        pmax = (fabs(cur->cov) > fabs(pmax->cov)) ? cur : pmax;
    }


    while (k+drop < niter) {

        sum_k += k;
        v = L + sum_k;
        dir = beta + sum_k;

        /* remove it from the unused set */
        pmax->prev->next = pmax->next;
        pmax->next->prev = pmax->prev;

        /* push pmax into the active set */
        pmax->prev = top_active;
        top_active = pmax;

        /* 
         * Compute the least-squares direction of the coefficients.
         */
        dir[k-drop] = temp = pmax->cov;

        for (i=k-drop-1, cur=top_active->prev; cur; cur=cur->prev, --i) {
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
            v[i] = cblas_ddot (nsamples, cur->ptr, m, 
                               pmax->ptr, m);

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
                     k-drop, L, v, 1);

        v[k] = sqrt(1 - cblas_ddot(k-drop, v, 1, v, 1));

        cblas_dtpsv (CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit,
                     k + 1, L, dir, 1);

        cblas_dtpsv (CblasRowMajor, CblasLower, CblasTrans, CblasNonUnit,
                     k + 1, L, dir, 1);

        /* 
         * Update uu with the current direction 
         * uu = Xa' * dir
         * 
         */
        cblas_dscal (nsamples, 0., uu, 1);
        for (i=k-drop, cur = top_active; cur; cur = cur->prev, --i)
            cblas_daxpy (nsamples, dir[i], cur->ptr, m, uu, 1);


        /*
         * Compute gamma
         */
        gamma = 1.;
        Aa = C  = fabs (pmax->cov) + DBL_EPSILON;
        for (cur = head->next; cur->ptr; cur = cur->next) {
            aj = cblas_ddot (nsamples, cur->ptr, m, uu, 1);

            cj = cur->cov;

            ga = (C - cj) / (Aa - aj);
            gb = (C + cj) / (Aa + aj);

            ga = ga > 0 ? ga : gamma;
            gb = gb > 0 ? gb : gamma;

            tgamma = fmin(ga, gb);
            gamma  = fmin(tgamma, gamma);
        }

        /* 
         * Set up return values.
         * TODO: we should iterate over used_indices (for Lasso variant).
         */
        cblas_daxpy (k,  1./gamma, dir-k, 1, dir, 1);
        cblas_dscal (k + 1, gamma, dir, 1);
        lambdas[k] = C;

        memcpy (row + sum_k, row + sum_k - k, k * sizeof(int));
        row[sum_k + k] = (pmax->ptr - X);
        for (i=0; i<=k; ++i) col[sum_k + i] = k + 1;

        switch (itype) {

        case 0:
            k++;
            break;

        case 1: /* Lasso modification */
            tgamma = gamma;
            for (i=0; i<k; ++i) {
                if (dir[i - k]*dir[i] < 0) {
                    printf("%d - %d\n", k , i);
                    /* it's possible that more than 1 variables
                     *  have changed sign
                     */
                    temp = - gamma * dir[i-k] / (dir[k] - dir[i-k]);
                    tgamma = fmin (tgamma, temp);
                    printf ("Change of sign!! %f \n", tgamma / gamma);
                }
            }

            if (tgamma < gamma) {

                for (i=0; i<k; ++i) {
                    dir[i] = dir[i-k] + (dir[i] - dir[i-k]) \
                        * tgamma / gamma ;

                }
                dir[k] = dir[k] * tgamma / gamma;
                gamma = tgamma;

                /* remove pmax from the active set */
                top_active = pmax->prev;

                pmax->prev = head;
                pmax->next = head->next;
                head->next->prev = pmax;
                head->next = pmax;

                drop ++;
            }
            break;
        }

        /* 
         * Update residual.
         */
        cblas_daxpy (nsamples, -gamma, uu, 1, res, 1);

        /* calculate next active set */
        for (pmax=head->next, cur=head->next; cur->ptr; cur=cur->next) {
            cur->cov = cblas_ddot (nsamples, cur->ptr, m, res, 1);
            pmax = (fabs(cur->cov) > fabs(pmax->cov)) ? cur : pmax;
        }

    }

    lambdas[k] = fabs(pmax->cov);
 
    free (active_set);
    free (uu);
}
