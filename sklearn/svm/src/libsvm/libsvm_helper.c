#include <stdlib.h>
#include <numpy/arrayobject.h>
#include "svm.h"

/*
 * Some helper methods for libsvm bindings.
 *
 * We need to access from python some parameters stored in svm_model
 * but libsvm does not expose this structure, so we define it here
 * along some utilities to convert from numpy arrays.
 *
 * License: BSD 3 clause
 *
 * Author: 2010 Fabian Pedregosa <fabian.pedregosa@inria.fr>
 */


/*
 * Convert matrix to sparse representation suitable for libsvm. x is
 * expected to be an array of length nrow*ncol.
 *
 * Typically the matrix will be dense, so we speed up the routine for
 * this case. We create a temporary array temp that collects non-zero
 * elements and after we just memcpy that to the proper array.
 *
 * Special care must be taken with indinces, since libsvm indices start
 * at 1 and not at 0.
 *
 * Strictly speaking, the C standard does not require that structs are
 * contiguous, but in practice its a reasonable assumption.
 *
 */
struct svm_node *dense_to_libsvm (double *x, npy_intp *dims)
{
    struct svm_node *node;
    npy_intp len_row = dims[1];
    double *tx = x;
    int i;

    node = malloc (dims[0] * sizeof(struct svm_node));

    if (node == NULL) return NULL;
    for (i=0; i<dims[0]; ++i) {
        node[i].values = tx;
        node[i].dim = (int) len_row;
        node[i].ind = i; /* only used if kernel=precomputed, but not
                            too much overhead */
        tx += len_row;
    }

    return node;
}


/*
 * Fill an svm_parameter struct.
 */
void set_parameter(struct svm_parameter *param, int svm_type, int kernel_type, int degree,
		double gamma, double coef0, double nu, double cache_size, double C,
		double eps, double p, int shrinking, int probability, int nr_weight,
		char *weight_label, char *weight, int max_iter, int random_seed)
{
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
    param->max_iter = max_iter;
    param->random_seed = random_seed;
}

/*
 * Fill an svm_problem struct. problem->x will be malloc'd.
 */
void set_problem(struct svm_problem *problem, char *X, char *Y, char *sample_weight, npy_intp *dims, int kernel_type)
{
    if (problem == NULL) return;
    problem->l = (int) dims[0]; /* number of samples */
    problem->y = (double *) Y;
    problem->x = dense_to_libsvm((double *) X, dims); /* implicit call to malloc */
    problem->W = (double *) sample_weight;
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
 * an int, so we might have numpy matrices that do not fit into libsvm's
 * data structure.
 *
 */
struct svm_model *set_model(struct svm_parameter *param, int nr_class,
                            char *SV, npy_intp *SV_dims,
                            char *support, npy_intp *support_dims,
                            npy_intp *sv_coef_strides,
                            char *sv_coef, char *rho, char *nSV,
                            char *probA, char *probB)
{
    struct svm_model *model;
    double *dsv_coef = (double *) sv_coef;
    int i, m;

    m = nr_class * (nr_class-1)/2;

    if ((model = malloc(sizeof(struct svm_model))) == NULL)
        goto model_error;
    if ((model->nSV = malloc(nr_class * sizeof(int))) == NULL)
        goto nsv_error;
    if ((model->label = malloc(nr_class * sizeof(int))) == NULL)
        goto label_error;
    if ((model->sv_coef = malloc((nr_class-1)*sizeof(double *))) == NULL)
        goto sv_coef_error;
    if ((model->rho = malloc( m * sizeof(double))) == NULL)
        goto rho_error;

    model->nr_class = nr_class;
    model->param = *param;
    model->l = (int) support_dims[0];

    if (param->kernel_type == PRECOMPUTED) {
        if ((model->SV = malloc ((model->l) * sizeof(struct svm_node))) == NULL)
            goto SV_error;
        for (i=0; i<model->l; ++i) {
            model->SV[i].ind = ((int *) support)[i];
            model->SV[i].values = NULL;
        }
    } else {
        model->SV = dense_to_libsvm((double *) SV, SV_dims);
    }
    /*
     * regression and one-class does not use nSV, label.
     * TODO: does this provoke memory leaks (we just malloc'ed them)?
     */
    if (param->svm_type < 2) {
        memcpy(model->nSV, nSV,     model->nr_class * sizeof(int));
        for(i=0; i < model->nr_class; i++)
            model->label[i] = i;
    }

    for (i=0; i < model->nr_class-1; i++) {
        model->sv_coef[i] = dsv_coef + i*(model->l);
    }

    for (i=0; i<m; ++i) {
        (model->rho)[i] = -((double *) rho)[i];
    }

    /*
     * just to avoid segfaults, these features are not wrapped but
     * svm_destroy_model will try to free them.
     */

    if (param->probability) {
        if ((model->probA = malloc(m * sizeof(double))) == NULL)
            goto probA_error;
        memcpy(model->probA, probA, m * sizeof(double));
        if ((model->probB = malloc(m * sizeof(double))) == NULL)
            goto probB_error;
        memcpy(model->probB, probB, m * sizeof(double));
    } else {
        model->probA = NULL;
        model->probB = NULL;
    }

    /* We'll free SV ourselves */
    model->free_sv = 0;
    return model;

probB_error:
    free(model->probA);
probA_error:
    free(model->SV);
SV_error:
    free(model->rho);
rho_error:
    free(model->sv_coef);
sv_coef_error:
    free(model->label);
label_error:
    free(model->nSV);
nsv_error:
    free(model);
model_error:
    return NULL;
}



/*
 * Get the number of support vectors in a model.
 */
npy_intp get_l(struct svm_model *model)
{
    return (npy_intp) model->l;
}

/*
 * Get the number of classes in a model, = 2 in regression/one class
 * svm.
 */
npy_intp get_nr(struct svm_model *model)
{
    return (npy_intp) model->nr_class;
}

/*
 * Some helpers to convert from libsvm sparse data structures
 * model->sv_coef is a double **, whereas data is just a double *,
 * so we have to do some stupid copying.
 */
void copy_sv_coef(char *data, struct svm_model *model)
{
    int i, len = model->nr_class-1;
    double *temp = (double *) data;
    for(i=0; i<len; ++i) {
        memcpy(temp, model->sv_coef[i], sizeof(double) * model->l);
        temp += model->l;
    }
}

void copy_intercept(char *data, struct svm_model *model, npy_intp *dims)
{
    /* intercept = -rho */
    npy_intp i, n = dims[0];
    double t, *ddata = (double *) data;
    for (i=0; i<n; ++i) {
        t = model->rho[i];
        /* we do this to avoid ugly -0.0 */
        *ddata = (t != 0) ? -t : 0;
        ++ddata;
    }
}

/*
 * This is a bit more complex since SV are stored as sparse
 * structures, so we have to do the conversion on the fly and also
 * iterate fast over data.
 */
void copy_SV(char *data, struct svm_model *model, npy_intp *dims)
{
    int i, n = model->l;
    double *tdata = (double *) data;
    int dim = model->SV[0].dim;
    for (i=0; i<n; ++i) {
        memcpy (tdata, model->SV[i].values, dim * sizeof(double));
        tdata += dim;
    }
}

void copy_support (char *data, struct svm_model *model)
{
    memcpy (data, model->sv_ind, (model->l) * sizeof(int));
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

void copy_probA(char *data, struct svm_model *model, npy_intp * dims)
{
    memcpy(data, model->probA, dims[0] * sizeof(double));
}

void copy_probB(char *data, struct svm_model *model, npy_intp * dims)
{
    memcpy(data, model->probB, dims[0] * sizeof(double));
}

/*
 * Predict using model.
 *
 *  It will return -1 if we run out of memory.
 */
int copy_predict(char *predict, struct svm_model *model, npy_intp *predict_dims,
                 char *dec_values)
{
    double *t = (double *) dec_values;
    struct svm_node *predict_nodes;
    npy_intp i;

    predict_nodes = dense_to_libsvm((double *) predict, predict_dims);

    if (predict_nodes == NULL)
        return -1;
    for(i=0; i<predict_dims[0]; ++i) {
        *t = svm_predict(model, &predict_nodes[i]);
        ++t;
    }
    free(predict_nodes);
    return 0;
}

int copy_predict_values(char *predict, struct svm_model *model,
                        npy_intp *predict_dims, char *dec_values, int nr_class)
{
    npy_intp i;
    struct svm_node *predict_nodes;
    predict_nodes = dense_to_libsvm((double *) predict, predict_dims);
    if (predict_nodes == NULL)
        return -1;
    for(i=0; i<predict_dims[0]; ++i) {
        svm_predict_values(model, &predict_nodes[i],
                                ((double *) dec_values) + i*nr_class);
    }

    free(predict_nodes);
    return 0;
}



int copy_predict_proba(char *predict, struct svm_model *model, npy_intp *predict_dims,
                 char *dec_values)
{
    npy_intp i, n, m;
    struct svm_node *predict_nodes;
    n = predict_dims[0];
    m = (npy_intp) model->nr_class;
    predict_nodes = dense_to_libsvm((double *) predict, predict_dims);
    if (predict_nodes == NULL)
        return -1;
    for(i=0; i<n; ++i) {
        svm_predict_probability(model, &predict_nodes[i],
                                ((double *) dec_values) + i*m);
    }
    free(predict_nodes);
    return 0;
}


/*
 * Some free routines. Some of them are nontrivial since a lot of
 * sharing happens across objects (they *must* be called in the
 * correct order)
 */

int free_model(struct svm_model *model)
{
    /* like svm_free_and_destroy_model, but does not free sv_coef[i] */
    if (model == NULL) return -1;
    free(model->SV);

    /* We don't free sv_ind, since we did not create them in
       set_model */
    /* free(model->sv_ind); */
    free(model->sv_coef);
    free(model->rho);
    free(model->label);
    free(model->probA);
    free(model->probB);
    free(model->nSV);
    free(model);

    return 0;
}

int free_param(struct svm_parameter *param)
{
    if (param == NULL) return -1;
    free(param);
    return 0;
}


/* borrowed from original libsvm code */
static void print_null(const char *s) {}

static void print_string_stdout(const char *s)
{
	fputs(s,stdout);
	fflush(stdout);
}

/* provide convenience wrapper */
void set_verbosity(int verbosity_flag){
	if (verbosity_flag)
		svm_set_print_string_function(&print_string_stdout);
	else
		svm_set_print_string_function(&print_null);
}
