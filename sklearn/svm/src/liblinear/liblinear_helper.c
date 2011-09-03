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
struct feature_node **dense_to_sparse (double *x, npy_intp *dims, double bias)
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
                free(sparse[i]);
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
c * Convert scipy.sparse.csr to libsvm's sparse data structure
 */
struct feature_node **csr_to_sparse (double *values, npy_intp *shape_indices,
		int *indices, npy_intp *shape_indptr, int *indptr, double bias,
                int n_features)
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
                free(sparse[i]);
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

struct problem * set_problem(char *X,char *Y, npy_intp *dims, double bias)
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

    problem->y = (int *) Y;
    problem->x = dense_to_sparse((double *) X, dims, bias);
    problem->bias = bias;
    if (problem->x == NULL) { 
        free(problem);
        return NULL;
    }

    return problem;
}

struct problem * csr_set_problem (char *values, npy_intp *n_indices,
	char *indices, npy_intp *n_indptr, char *indptr, char *Y,
        npy_intp n_features, double bias) {

    struct problem *problem;
    problem = malloc (sizeof (struct problem));
    if (problem == NULL) return NULL;
    problem->l = (int) n_indptr[0] -1;

    if (bias > 0){
        problem->n = (int) n_features + 1;
    } else {
        problem->n = (int) n_features;
    }

    problem->y = (int *) Y;
    problem->x = csr_to_sparse((double *) values, n_indices, (int *) indices,
			n_indptr, (int *) indptr, bias, n_features);
    problem->bias = bias;

    if (problem->x == NULL) {
        free(problem);
        return NULL;
    }

    return problem;
}


/* Create a paramater struct with and return it */
struct parameter * set_parameter(int solver_type, double eps, double C, npy_intp nr_weight, char *weight_label, char *weight)
{
    struct parameter *param;
    param = malloc(sizeof(struct parameter));
    if (param == NULL) return NULL;
    param->solver_type = solver_type;
    param->eps = eps;
    param->C = C;
    param->nr_weight = (int) nr_weight;
    param->weight_label = (int *) weight_label;
    param->weight = (double *) weight;
    return param;
}

struct model * set_model(struct parameter *param, char *coef, npy_intp *dims, 
                         char *label, double bias)
{
    npy_intp len_w = dims[0] * dims[1];
    int m = (int) dims[0], k = (int) dims[1];
    struct model *model;

    if (m == 1) m = 2; /* liblinear collapses the weight vector in the case of two classes */
    if ((model = malloc(sizeof(struct model))) == NULL)
        goto model_error;
    if ((model->w = malloc( len_w * sizeof(double))) == NULL)
        goto w_error;
    if ((model->label = malloc( m * sizeof(int))) == NULL)
        goto label_error;

    memcpy(model->label, label, m * sizeof(int));
    memcpy(model->w, coef, len_w * sizeof(double));

    model->nr_feature = bias > 0 ? k - 1 : k;
    model->nr_class = m;
	
    model->param = *param;
    model->bias = bias;

    return model;

label_error:
    free(model->w);
w_error:
    free(model);
model_error:
    return NULL;
}


void copy_w(char *data, struct model *model, int len)
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
}

void free_parameter(struct parameter *param)
{
    free(param);
}

int copy_predict(char *train, struct model *model_, npy_intp *train_dims,
                 char *dec_values)
{
    int *t = (int *) dec_values;
    int i, n;
    struct feature_node **train_nodes;
    n = train_dims[0];
    train_nodes = dense_to_sparse((double *) train, train_dims, model_->bias);
    if (train_nodes == NULL)
        return -1;
    for(i=0; i<n; ++i) {
        *t = predict(model_, train_nodes[i]);
        free(train_nodes[i]);
        ++t;
    }
    free(train_nodes);
    return 0;
}
/*
 * Predict using a model, where data is expected to be enconded into a csr matrix.
 */
int csr_copy_predict(npy_intp n_features, npy_intp *data_size, char *data,
        npy_intp *index_size,
        char *index, npy_intp *indptr_shape, char *intptr, struct model *model_,
        char *dec_values)
{
    int *t = (int *) dec_values;
    struct feature_node **predict_nodes;
    npy_intp i;

    predict_nodes = csr_to_sparse((double *) data, index_size,
            (int *) index, indptr_shape, (int *) intptr, model_->bias, n_features);

    if (predict_nodes == NULL)
        return -1;
    for (i = 0; i < indptr_shape[0] - 1; ++i) {
        *t = predict(model_, predict_nodes[i]);
        free(predict_nodes[i]);
        ++t;
    }
    free(predict_nodes);
    return 0;
}

int copy_predict_values (char *predict, struct model *model_, 
                         npy_intp *predict_dims, char *dec_values, int nr_class)
{
    npy_intp i;
    struct feature_node **predict_nodes;
    predict_nodes = dense_to_sparse((double *) predict, predict_dims, model_->bias);
    if (predict_nodes == NULL)
        return -1;
    for(i=0; i<predict_dims[0]; ++i) {
        predict_values(model_, predict_nodes[i], 
                       ((double *) dec_values) + i*nr_class);
        free(predict_nodes[i]);
    }

    free(predict_nodes);
    return 0;
}

int csr_copy_predict_values(npy_intp n_features, npy_intp *data_size,
                            char *data, npy_intp *index_size, char
                            *index, npy_intp *indptr_shape, char
                            *intptr, struct model *model_, char
                            *dec_values, int nr_class)
{
    struct feature_node **predict_nodes;
    npy_intp i;

    predict_nodes = csr_to_sparse((double *) data, index_size,
                                  (int *) index, indptr_shape, (int *) intptr, model_->bias, n_features);

    if (predict_nodes == NULL)
        return -1;
    for (i = 0; i < indptr_shape[0] - 1; ++i) {
        predict_values(model_, predict_nodes[i],
                       ((double *) dec_values) + i*nr_class);
        free(predict_nodes[i]);
    }
    free(predict_nodes);
    return 0;
}


int copy_prob_predict(char *predict, struct model *model_, npy_intp *predict_dims,
                 char *dec_values)
{
    struct feature_node **predict_nodes;
    int i;
    int n, m;
    n = predict_dims[0];
    m = model_->nr_class;
    predict_nodes = dense_to_sparse((double *) predict, predict_dims, model_->bias);
    if (predict_nodes == NULL)
        return -1;
    for(i=0; i<n; ++i) {
        predict_probability(model_, predict_nodes[i],
                            ((double *) dec_values) + i*m);
        free(predict_nodes[i]);
    }
    free(predict_nodes);
    return 0;
}


int csr_copy_predict_proba(npy_intp n_features, npy_intp *data_size, 
                           char *data, npy_intp *index_size,
                           char *index, npy_intp *indptr_shape, 
                           char *indptr, struct model *model_,
                           char *dec_values)
{
    struct feature_node **predict_nodes;
    int i;
    double *tx = (double *) dec_values;

    predict_nodes = csr_to_sparse((double *) data, index_size,
                                  (int *) index, indptr_shape, 
                                  (int *) indptr, model_->bias, 
                                  n_features);
    if (predict_nodes == NULL)
        return -1;
    for(i=0; i<indptr_shape[0] - 1; ++i) {
        predict_probability(model_, predict_nodes[i], tx);
        tx += model_->nr_class;
        free(predict_nodes[i]);
    }
    free(predict_nodes);
    return 0;
}


int copy_label(char *data, struct model *model_, int nr_class)
{
    memcpy(data, model_->label, nr_class * sizeof(int));
    return 0;
}
