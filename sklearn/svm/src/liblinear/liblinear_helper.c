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
static struct feature_node **dense_to_sparse(double *x, int n_samples,
        int n_features, int n_nonzero, double bias)
{
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
            if (*x != 0) {
                T->value = *x;
                T->index = j;
                ++ T;
            }
            ++ x; /* go to next element */
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
 * Convert scipy.sparse.csr to libsvm's sparse data structure
 */
static struct feature_node **csr_to_sparse(double *values, int *indices,
        int *indptr, int n_samples, int n_features, int n_nonzero, double bias)
{
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
            T->value = values[k];
            T->index = indices[k] + 1; /* libsvm uses 1-based indexing */
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

struct problem * set_problem(char *X, char *Y, int n_samples, int n_features,
        int n_nonzero, double bias, char* sample_weight)
{
    struct problem *problem;
    /* not performant but simple */
    problem = malloc(sizeof(struct problem));
    if (problem == NULL) return NULL;
    problem->l = n_samples;

    if (bias > 0) {
        problem->n = (int) n_features + 1;
    } else {
        problem->n = (int) n_features;
    }

    problem->y = (double *) Y;
    problem->sample_weight = (double *) sample_weight;
    problem->x = dense_to_sparse((double *) X, n_samples, n_features, n_nonzero, bias);
    problem->bias = bias;
    problem->sample_weight = sample_weight;
    if (problem->x == NULL) { 
        free(problem);
        return NULL;
    }

    return problem;
}

struct problem * csr_set_problem (char *values, char *indices, char *indptr,
        char *Y, int n_samples, int n_features, int n_nonzero, double bias,
        char *sample_weight) {

    struct problem *problem;
    problem = malloc (sizeof (struct problem));
    if (problem == NULL) return NULL;
    problem->l = n_samples;
    problem->sample_weight = (double *) sample_weight;

    if (bias > 0){
        problem->n = (int) n_features + 1;
    } else {
        problem->n = (int) n_features;
    }

    problem->y = (double *) Y;
    problem->x = csr_to_sparse((double *) values, (int *) indices,
                        (int *) indptr, n_samples, n_features, n_nonzero, bias);
    problem->bias = bias;
    problem->sample_weight = sample_weight;

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
