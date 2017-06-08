#include <stdlib.h>
#include <numpy/arrayobject.h>
#include "linear.h"

/*
 * Convert matrix to sparse representation suitable for libsvm. x is
 * expected to be an array of length nrow*ncol.
 *
 * Typically the matrix will be dense, so we speed up the routine for
 * this case. We create a temporary array temp that collects non-zero
 * elements and after we just memcpy that to the proper array.
 *
 * Special care must be taken with indices, since libsvm indices start
 * at 1 and not at 0.
 *
 * If bias is > 0, we append an item at the end.
 */
static struct feature_node **dense_to_sparse(double *x, npy_intp *dims,
                                             double bias)
{
    struct feature_node **sparse;
    int i, j;                           /* number of nonzero elements in row i */
    struct feature_node *temp;          /* stack for nonzero elements */
    struct feature_node *T;             /* pointer to the top of the stack */
    int count;

    sparse = malloc (dims[0] * sizeof(struct feature_node *));
    if (sparse == NULL)
        goto sparse_error;

    temp = malloc ((dims[1]+2) * sizeof(struct feature_node));
    if (temp == NULL)
        goto temp_error;

    for (i=0; i<dims[0]; ++i) {
        T = temp; /* reset stack pointer */

        for (j=1; j<=dims[1]; ++j) {
            if (*x != 0) {
                T->value = *x;
                T->index = j;
                ++ T;
            }
            ++ x; /* go to next element */
        }

        /* set bias element */
        if (bias > 0) {
                T->value = bias;
                T->index = j;
                ++ T;
            }

        /* set sentinel */
        T->index = -1;
        ++ T;

        /* allocate memory and copy collected items*/
        count = T - temp;
        sparse[i] = malloc(count * sizeof(struct feature_node));
        if (sparse[i] == NULL) {
            int k;
            for (k=0; k<i; k++)
                free(sparse[k]);
            goto sparse_i_error;
        }
        memcpy(sparse[i], temp, count * sizeof(struct feature_node));
    }

    free(temp);
    return sparse;

sparse_i_error:
    free(temp);
temp_error:
    free(sparse);
sparse_error:
    return NULL;
}


/*
 * Convert scipy.sparse.csr to libsvm's sparse data structure
 */
static struct feature_node **csr_to_sparse(double *values,
        npy_intp *shape_indices, int *indices, npy_intp *shape_indptr,
        int *indptr, double bias, int n_features)
{
    struct feature_node **sparse, *temp;
    int i, j=0, k=0, n;

    sparse = malloc ((shape_indptr[0]-1)* sizeof(struct feature_node *));
    if (sparse == NULL)
        return NULL;

    for (i=0; i<shape_indptr[0]-1; ++i) {
        n = indptr[i+1] - indptr[i]; /* count elements in row i */

        sparse[i] = malloc ((n+2) * sizeof(struct feature_node));
        if (sparse[i] == NULL) {
            int l;
            for (l=0; l<i; l++)
                free(sparse[l]);
            break;
        }

        temp = sparse[i];
        for (j=0; j<n; ++j) {
            temp[j].value = values[k];
            temp[j].index = indices[k] + 1; /* libsvm uses 1-based indexing */
            ++k;
        }

        if (bias > 0) {
            temp[j].value = bias;
            temp[j].index = n_features + 1;
            ++j;
        }

        /* set sentinel */
        temp[j].index = -1;
    }

    return sparse;
}

struct problem * set_problem(char *X,char *Y, npy_intp *dims, double bias, char* sample_weight)
{
    struct problem *problem;
    /* not performant but simple */
    problem = malloc(sizeof(struct problem));
    if (problem == NULL) return NULL;
    problem->l = (int) dims[0];

    if (bias > 0) {
        problem->n = (int) dims[1] + 1;
    } else {
        problem->n = (int) dims[1];
    }

    problem->y = (double *) Y;
    problem->sample_weight = (double *) sample_weight;
    problem->x = dense_to_sparse((double *) X, dims, bias);
    problem->bias = bias;
    problem->sample_weight = sample_weight;
    if (problem->x == NULL) { 
        free(problem);
        return NULL;
    }

    return problem;
}

struct problem * csr_set_problem (char *values, npy_intp *n_indices,
	char *indices, npy_intp *n_indptr, char *indptr, char *Y,
        npy_intp n_features, double bias, char *sample_weight) {

    struct problem *problem;
    problem = malloc (sizeof (struct problem));
    if (problem == NULL) return NULL;
    problem->l = (int) n_indptr[0] -1;
    problem->sample_weight = (double *) sample_weight;

    if (bias > 0){
        problem->n = (int) n_features + 1;
    } else {
        problem->n = (int) n_features;
    }

    problem->y = (double *) Y;
    problem->x = csr_to_sparse((double *) values, n_indices, (int *) indices,
			n_indptr, (int *) indptr, bias, n_features);
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
    int i;
    for(i=problem->l-1; i>=0; --i) free(problem->x[i]);
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
