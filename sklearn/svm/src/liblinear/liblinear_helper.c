#include <stdlib.h>
#include <numpy/arrayobject.h>
#include "linear.h"

/*
 * Convert matrix to sparse representation suitable for liblinear. x is
 * expected to be an array of length n_samples*n_features.
 *
 * Whether the matrix is densely or sparsely populated, the fastest way to
 * convert it to liblinear's sparse format is to calculate the amount of memory
 * needed and allocate a single big block.
 *
 * Special care must be taken with indices, since liblinear indices start at 1
 * and not at 0.
 *
 * If bias is > 0, we append an item at the end.
 */
static struct feature_node **dense_to_sparse(char *x, int double_precision,
        int n_samples, int n_features, int n_nonzero, double bias)
{
    float *x32 = (float *)x;
    double *x64 = (double *)x;
    struct feature_node **sparse;
    int i, j;                           /* number of nonzero elements in row i */
    struct feature_node *T;             /* pointer to the top of the stack */
    int have_bias = (bias > 0);

    sparse = malloc (n_samples * sizeof(struct feature_node *));
    if (sparse == NULL)
        return NULL;

    n_nonzero += (have_bias+1) * n_samples;
    T = malloc (n_nonzero * sizeof(struct feature_node));
    if (T == NULL) {
        free(sparse);
        return NULL;
    }

    for (i=0; i<n_samples; ++i) {
        sparse[i] = T;

        for (j=1; j<=n_features; ++j) {
            if (double_precision) {
                if (*x64 != 0) {
                    T->value = *x64;
                    T->index = j;
                    ++ T;
                }
                ++ x64; /* go to next element */
            } else {
                if (*x32 != 0) {
                    T->value = *x32;
                    T->index = j;
                    ++ T;
                }
                ++ x32; /* go to next element */
            }
        }

        /* set bias element */
        if (have_bias) {
                T->value = bias;
                T->index = j;
                ++ T;
            }

        /* set sentinel */
        T->index = -1;
        ++ T;
    }

    return sparse;
}


/*
 * Convert scipy.sparse.csr to liblinear's sparse data structure
 */
static struct feature_node **csr_to_sparse(char *x, int double_precision,
        int *indices, int *indptr, int n_samples, int n_features, int n_nonzero,
        double bias)
{
    float *x32 = (float *)x;
    double *x64 = (double *)x;
    struct feature_node **sparse;
    int i, j=0, k=0, n;
    struct feature_node *T;
    int have_bias = (bias > 0);

    sparse = malloc (n_samples * sizeof(struct feature_node *));
    if (sparse == NULL)
        return NULL;

    n_nonzero += (have_bias+1) * n_samples;
    T = malloc (n_nonzero * sizeof(struct feature_node));
    if (T == NULL) {
        free(sparse);
        return NULL;
    }

    for (i=0; i<n_samples; ++i) {
        sparse[i] = T;
        n = indptr[i+1] - indptr[i]; /* count elements in row i */

        for (j=0; j<n; ++j) {
            T->value = double_precision ? x64[k] : x32[k];
            T->index = indices[k] + 1; /* liblinear uses 1-based indexing */
            ++T;
            ++k;
        }

        if (have_bias) {
            T->value = bias;
            T->index = n_features + 1;
            ++T;
            ++j;
        }

        /* set sentinel */
        T->index = -1;
        ++T;
    }

    return sparse;
}

struct problem * set_problem(char *X, int double_precision_X, int n_samples,
        int n_features, int n_nonzero, double bias, char* sample_weight,
        char *Y)
{
    struct problem *problem;
    /* not performant but simple */
    problem = malloc(sizeof(struct problem));
    if (problem == NULL) return NULL;
    problem->l = n_samples;
    problem->n = n_features + (bias > 0);
    problem->y = (double *) Y;
    problem->sample_weight = (double *) sample_weight;
    problem->x = dense_to_sparse(X, double_precision_X, n_samples, n_features,
                        n_nonzero, bias);
    problem->bias = bias;

    if (problem->x == NULL) { 
        free(problem);
        return NULL;
    }

    return problem;
}

struct problem * csr_set_problem (char *X, int double_precision_X,
        char *indices, char *indptr, int n_samples, int n_features,
        int n_nonzero, double bias, char *sample_weight, char *Y)
{
    struct problem *problem;
    problem = malloc (sizeof (struct problem));
    if (problem == NULL) return NULL;
    problem->l = n_samples;
    problem->n = n_features + (bias > 0);
    problem->y = (double *) Y;
    problem->sample_weight = (double *) sample_weight;
    problem->x = csr_to_sparse(X, double_precision_X, (int *) indices,
                        (int *) indptr, n_samples, n_features, n_nonzero, bias);
    problem->bias = bias;

    if (problem->x == NULL) {
        free(problem);
        return NULL;
    }

    return problem;
}


/* Create a paramater struct with and return it */
struct parameter *set_parameter(int solver_type, double eps, double C,
                                npy_intp nr_weight, char *weight_label,
                                char *weight, int max_iter, unsigned seed, 
                                double epsilon)
{
    struct parameter *param = malloc(sizeof(struct parameter));
    if (param == NULL)
        return NULL;

    srand(seed);
    param->solver_type = solver_type;
    param->eps = eps;
    param->C = C;
    param->p = epsilon;  // epsilon for epsilon-SVR
    param->nr_weight = (int) nr_weight;
    param->weight_label = (int *) weight_label;
    param->weight = (double *) weight;
    param->max_iter = max_iter;
    return param;
}

void copy_w(void *data, struct model *model, int len)
{
    memcpy(data, model->w, len * sizeof(double)); 
}

double get_bias(struct model *model)
{
    return model->bias;
}

void free_problem(struct problem *problem)
{
    free(problem->x[0]);
    free(problem->x);
    free(problem);
}

void free_parameter(struct parameter *param)
{
    free(param);
}

/* rely on built-in facility to control verbose output */
static void print_null(const char *s) {}

static void print_string_stdout(const char *s)
{
    fputs(s ,stdout);
    fflush(stdout);
}

/* provide convenience wrapper */
void set_verbosity(int verbosity_flag){
    if (verbosity_flag)
        set_print_string_function(&print_string_stdout);
    else
        set_print_string_function(&print_null);
}
