#include <stdlib.h>
#include "svm.h"
#include <numpy/arrayobject.h>

/*
 * Some helper methods for libsvm bindings.
 *
 * We need to access from python some parameters stored in svm_model
 * but libsvm does not expose this structure, so we define it here
 * along some utilities to convert from numpy arrays.
 *
 * License: New BSD.
 *
 * Author: 2010 Fabian Pedregosa <fabian.pedregosa@inria.fr>
 */

struct svm_model
{
    struct svm_parameter param;   // parameter
    int nr_class;                 // number of classes, = 2 in
                                  // regression/one class svm
    int l;                        // total #SV
    struct svm_node **SV;         // SVs (SV[l])
    double **sv_coef;             // coefficients for SVs in decision
                                  // functions (sv_coef[k-1][l])
    double *rho;                  // constants in decision functions
                                  // (rho[k*(k-1)/2])
    double *probA;                // pairwise probability information
    double *probB;

    // for classification only

    int *label;     // label of each class (label[k])
    int *nSV;       // number of SVs for each class (nSV[k])
                    // nSV[0] + nSV[1] + ... + nSV[k-1] = l
    // XXX
    int free_sv;    // 1 if svm_model is created by svm_load_model
                    // 0 if svm_model is created by svm_train
};



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
 */
struct svm_node **dense_to_sparse (double *x, npy_intp *dims)
{
    struct svm_node **sparse;
    register int i, j;              /* number of nonzero elements in row i */
    struct svm_node *temp;          /* stack for nonzero elements */
    struct svm_node *T;             /* pointer to the top of the stack */
    int count;

    sparse = (struct svm_node **) malloc (dims[0] * sizeof(struct svm_node *));
    temp = (struct svm_node *) malloc ((dims[1]+1) * sizeof(struct svm_node));

    if (sparse == NULL || temp == NULL) return NULL;

    for (i=0; i<dims[0]; ++i) {
        T = temp; /* reset stack pointer */

        for (j=1; j<=dims[1]; ++j) {
            T->value = *x;
            if (T->value != 0) {
                T->index = j;
                ++T;
            }
            ++x; /* go to next element */
        }

        /* set sentinel */
        T->index = -1;
        ++T;

        /* allocate memory and copy collected items*/
        count = T - temp;
        sparse[i] = (struct svm_node *) malloc(count * sizeof(struct svm_node));
        if (sparse[i] == NULL) return NULL;
        memcpy(sparse[i], temp, count * sizeof(struct svm_node));
    }

    free(temp);
    return sparse;
}

/* Create a svm_paramater struct with and return it */
struct svm_parameter * set_parameter(int svm_type, int kernel_type, int degree,  double gamma,
                                     double coef0, double nu, double cache_size, 
                                     double C, double eps, double p, int shrinking, 
                                     int probability, int nr_weight, char *weight_label, 
                                     char *weight)
{
    struct svm_parameter *param;
    param = (struct svm_parameter *) malloc(sizeof(struct svm_parameter));
    if (param == NULL) return NULL;
    param->svm_type = svm_type;
    param->kernel_type = kernel_type;
    param->degree = degree;
    param->coef0 = coef0;
    param->nu = nu;
    param->cache_size = cache_size;
    param->C = C;
    param->eps = eps;
    param->p = p;
    param->shrinking = shrinking;
    param->probability = probability;
    param->nr_weight = nr_weight;
    param->weight_label = (int *) weight_label;
    param->weight = (double *) weight;
    param->gamma = gamma;
    return param;
}


struct svm_problem * set_problem(char *X,char *Y, npy_intp *dims)
{
    struct svm_problem *problem;
    problem = (struct svm_problem *) malloc(sizeof(struct svm_problem));
    if (problem == NULL) return NULL;
    problem->l = (int) dims[0];
    problem->y = (double *) Y;
    problem->x = dense_to_sparse((double *) X, dims);
    if (problem->x == NULL) { 
        free(problem);
        return NULL;
    }
    return problem;
}

/*
 * Create and return an instance of svm_model.
 *
 * The copy of model->sv_coef should be straightforward, but
 * unfortunately to represent a matrix numpy and libsvm use different
 * approaches, so it requires some iteration.  
 *
 * Possible issue: on 64 bits, the number of columns that numpy can 
 * store is a long, but libsvm enforces this number (model->l) to be
 * an int.
 *
 * TODO: check that malloc was successful.
 */
struct svm_model *set_model(struct svm_parameter *param, int nr_class,
                            char *SV, npy_intp *SV_dims, npy_intp *sv_coef_strides, 
                            char *sv_coef, char *rho, char *nSV, char *label, 
                            char *probA, char *probB)
{
    int i;
    char *t = sv_coef;
    int m = nr_class*(nr_class-1)/2;
    struct svm_model *model;

    model = (struct svm_model *)  malloc(sizeof(struct svm_model));
    model->nSV =     (int *)      malloc(nr_class * sizeof(int));
    model->label =   (int *)      malloc(nr_class * sizeof(int));;
    model->sv_coef = (double **)  malloc((nr_class-1)*sizeof(double *));
    model->rho =     (double *)   malloc( m * sizeof(double));
    model->SV = dense_to_sparse((double *) SV, SV_dims);

    model->nr_class = nr_class;
    model->param = *param;
    model->l = (int) SV_dims[0];

    /* 
     * regression and one-class does not use nSV, label.
     * TODO: does this provoque memory leaks (we just malloc'ed them)?
     */
    if (param->svm_type < 2) {
        memcpy(model->nSV, nSV, model->nr_class * sizeof(int));
        memcpy(model->label, label, model->nr_class * sizeof(int));
    }

    for (i=0; i < model->nr_class-1; i++) {
        /*
         * We cannot squash all this mallocs in a single call since
         * svm_destroy_model will free each element of the array.
         */
        model->sv_coef[i] = (double *) malloc((model->l) * sizeof(double));
        memcpy(model->sv_coef[i], t, (model->l) * sizeof(double));
        t += sv_coef_strides[0]; 
    }

    memcpy(model->rho, rho, m * sizeof(double));
    /* 
     * just to avoid segfaults, these features are not wrapped but
     * svm_destroy_model will try to free them.
     */

    if (param->probability) {
        model->probA = (double *) malloc(m * sizeof(double));
        memcpy(model->probA, probA, m * sizeof(double));
        model->probB = (double *) malloc(m * sizeof(double));
        memcpy(model->probB, probB, m * sizeof(double));
    } else {
        model->probA = NULL;
        model->probB = NULL;
    }

    /* We'll free SV ourselves */
    model->free_sv = 0;
    return model;
}

/*
 * Get the number of support vectors in a model.
 */
npy_intp get_l(struct svm_model *model)
{
    return model->l;
}

/*
 * Get the number of classes in a model, = 2 in regression/one class
 * svm.
 */
npy_intp get_nr(struct svm_model *model)
{
    return model->nr_class;
}

/*
 * Some helpers to convert from libsvm sparse data structures 
 * model->sv_coef is a double **, whereas data is just a double *,
 * so we have to do some stupid copying.
 */
void copy_sv_coef(char *data, struct svm_model *model, npy_intp *strides)
{
    register int i;
    char *temp = data;
    npy_intp step = strides[0];
    int len = model->nr_class-1;
    for(i=0; i<len; ++i) {
        memcpy(temp, model->sv_coef[i], step);
        temp += step;
    }
}

void copy_rho(char *data, struct svm_model *model, npy_intp *dims)
{
    memcpy(data, model->rho, dims[0] * sizeof(double));
}

/* 
 * This is a bit more complex since SV are stored as sparse
 * structures, so we have to do the conversion on the fly and also
 * iterate fast over data.
 */
void copy_SV(char *data, struct svm_model *model, npy_intp *strides)
{
    register int i, j;
    char *t = data;
    int k, n = model->l;
    npy_intp step = strides[1];
    for (i=0; i<n; ++i) {
        k = model->SV[i][0].index;
        for(j=0; k >= 0; ++j) {
            * ((double *) (t + (k-1)*step)) = model->SV[i][j].value;
            k = model->SV[i][j+1].index;
        }
        t += strides[0];
    }
}

/* 
 * copy svm_model.nSV, an array with the number of SV for each class 
 * will be NULL in the case of SVR, OneClass
 */
void copy_nSV(char *data, struct svm_model *model)
{
    if (model->label == NULL) return;
    memcpy(data, model->nSV, model->nr_class * sizeof(int));
}

/* 
 * same as above with model->label
 * TODO: maybe merge into the previous?
 */
void copy_label(char *data, struct svm_model *model)
{
    if (model->label == NULL) return;
    memcpy(data, model->label, model->nr_class * sizeof(int));
}

void copy_probA(char *data, struct svm_model *model, npy_intp * dims)
{
    memcpy(data, model->probA, dims[0] * sizeof(double));
}


void copy_probB(char *data, struct svm_model *model, npy_intp * dims)
{
    memcpy(data, model->probB, dims[0] * sizeof(double));
}

/* 
 * Use train to predict using model. Train is expected to be an array
 * in the numpy sense.
 *
 *  It will return -1 if we run out of memory.
 */
int copy_predict(char *predict, struct svm_model *model, npy_intp *predict_dims,
                 char *dec_values)
{
    double *t = (double *) dec_values;
    register int i, n;
    n = predict_dims[0];
    struct svm_node **predict_nodes;
    predict_nodes = dense_to_sparse((double *) predict, predict_dims);
    if (predict_nodes == NULL)
        return -1;
    for(i=0; i<n; ++i) {
        *t = svm_predict(model, predict_nodes[i]);
        free(predict_nodes[i]);
        ++t;
    }
    return 0;
}

int copy_prob_predict(char *predict, struct svm_model *model, npy_intp *predict_dims,
                 char *dec_values)
{
    double *t = (double *) dec_values;
    register int i, n;
    n = predict_dims[0];
    struct svm_node **predict_nodes;
    predict_nodes = dense_to_sparse((double *) predict, predict_dims);
    if (predict_nodes == NULL)
        return -1;
    for(i=0; i<n; ++i) {
        svm_predict_probability(model, predict_nodes[i], &t[i*predict_dims[1]]);
        free(predict_nodes[i]);
    }
    return 0;
}


/* 
 * Some free routines. Some of them are nontrivial since a lot of
 * sharing happens across objects (they *must* be called in the
 * correct order)
 */
int free_problem(struct svm_problem *problem)
{
    if (problem == NULL) return -1;
    register int i;
    for (i=problem->l-1; i>=0; --i) free (problem->x[i]);
    free (problem->x);
    return 0;
}

int free_model(struct svm_model *model)
{
    if (model == NULL) return -1;
    svm_destroy_model(model);
    return 0;
}

int free_model_SV(struct svm_model *model)
{
    register int i;
    for (i=model->l-1; i>=0; --i) free(model->SV[i]);
    /* svn_destroy_model frees model->SV */
    return 0;
}

int free_param(struct svm_parameter *param)
{
    if (param == NULL) return -1;
    free(param);
    return 0;
}
       
                            
