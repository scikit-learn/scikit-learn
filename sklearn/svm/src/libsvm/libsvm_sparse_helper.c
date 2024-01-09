#include <stdlib.h>
#include <numpy/arrayobject.h>
#include "svm.h"
#include "_svm_cython_blas_helpers.h"


#ifndef MAX
    #define MAX(x, y) (((x) > (y)) ? (x) : (y))
#endif


/*
 * Convert scipy.sparse.csr to libsvm's sparse data structure
 */
struct svm_csr_node **csr_to_libsvm (double *values, int* indices, int* indptr, npy_int n_samples)
{
    struct svm_csr_node **sparse, *temp;
    int i, j=0, k=0, n;
    sparse = malloc (n_samples * sizeof(struct svm_csr_node *));

    if (sparse == NULL)
        return NULL;

    for (i=0; i<n_samples; ++i) {
        n = indptr[i+1] - indptr[i]; /* count elements in row i */
        temp = malloc ((n+1) * sizeof(struct svm_csr_node));

        if (temp == NULL) {
            for (j=0; j<i; j++)
                free(sparse[j]);
            free(sparse);
            return NULL;
        }

        for (j=0; j<n; ++j) {
            temp[j].value = values[k];
            temp[j].index = indices[k] + 1; /* libsvm uses 1-based indexing */
            ++k;
        }
        /* set sentinel */
        temp[n].index = -1;
        sparse[i] = temp;
    }

    return sparse;
}



struct svm_parameter * set_parameter(int svm_type, int kernel_type, int degree,
		double gamma, double coef0, double nu, double cache_size, double C,
		double eps, double p, int shrinking, int probability, int nr_weight,
		char *weight_label, char *weight, int max_iter, int random_seed)
{
    struct svm_parameter *param;
    param = malloc(sizeof(struct svm_parameter));
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
    param->max_iter = max_iter;
    param->random_seed = random_seed;
    return param;
}


/*
 * Create and return a svm_csr_problem struct from a scipy.sparse.csr matrix. It is
 * up to the user to free resulting structure.
 *
 * TODO: precomputed kernel.
 */
struct svm_csr_problem * csr_set_problem (char *values, npy_intp *n_indices,
		char *indices, npy_intp *n_indptr, char *indptr, char *Y,
                char *sample_weight, int kernel_type) {

    struct svm_csr_problem *problem;
    problem = malloc (sizeof (struct svm_csr_problem));
    if (problem == NULL) return NULL;
    problem->l = (int) n_indptr[0] - 1;
    problem->y = (double *) Y;
    problem->x = csr_to_libsvm((double *) values, (int *) indices,
                               (int *) indptr, problem->l);
    /* should be removed once we implement weighted samples */
    problem->W = (double *) sample_weight;

    if (problem->x == NULL) {
        free(problem);
        return NULL;
    }
    return problem;
}


struct svm_csr_model *csr_set_model(struct svm_parameter *param, int nr_class,
                            char *SV_data, npy_intp *SV_indices_dims,
                            char *SV_indices, npy_intp *SV_indptr_dims,
                            char *SV_intptr,
                            char *sv_coef, char *rho, char *nSV,
                            char *probA, char *probB)
{
    struct svm_csr_model *model;
    double *dsv_coef = (double *) sv_coef;
    int i, m;

    m = nr_class * (nr_class-1)/2;

    if ((model = malloc(sizeof(struct svm_csr_model))) == NULL)
        goto model_error;
    if ((model->nSV = malloc(nr_class * sizeof(int))) == NULL)
        goto nsv_error;
    if ((model->label = malloc(nr_class * sizeof(int))) == NULL)
        goto label_error;
    if ((model->sv_coef = malloc((nr_class-1)*sizeof(double *))) == NULL)
        goto sv_coef_error;
    if ((model->rho = malloc( m * sizeof(double))) == NULL)
        goto rho_error;

    // This is only allocated in dynamic memory while training.
    model->n_iter = NULL;

    /* in the case of precomputed kernels we do not use
       dense_to_precomputed because we don't want the leading 0. As
       indices start at 1 (not at 0) this will work */
    model->l = (int) SV_indptr_dims[0] - 1;
    model->SV = csr_to_libsvm((double *) SV_data, (int *) SV_indices,
                              (int *) SV_intptr, model->l);
    model->nr_class = nr_class;
    model->param = *param;

    /*
     * regression and one-class does not use nSV, label.
     */
    if (param->svm_type < 2) {
        memcpy(model->nSV,   nSV,   model->nr_class * sizeof(int));
        for(i=0; i < model->nr_class; i++)
            model->label[i] = i;
    }

    for (i=0; i < model->nr_class-1; i++) {
        /*
         * We cannot squash all this mallocs in a single call since
         * svm_destroy_model will free each element of the array.
         */
        if ((model->sv_coef[i] = malloc((model->l) * sizeof(double))) == NULL) {
            int j;
            for (j=0; j<i; j++)
                free(model->sv_coef[j]);
            goto sv_coef_i_error;
        }
        memcpy(model->sv_coef[i], dsv_coef, (model->l) * sizeof(double));
        dsv_coef += model->l;
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
    for (i=0; i < model->nr_class-1; i++)
        free(model->sv_coef[i]);
sv_coef_i_error:
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
 * Copy support vectors into a scipy.sparse.csr matrix
 */
int csr_copy_SV (char *data, npy_intp *n_indices,
		char *indices, npy_intp *n_indptr, char *indptr,
		struct svm_csr_model *model, int n_features)
{
	int i, j, k=0, index;
	double *dvalues = (double *) data;
	int *iindices = (int *) indices;
	int *iindptr  = (int *) indptr;
	iindptr[0] = 0;
	for (i=0; i<model->l; ++i) { /* iterate over support vectors */
		index = model->SV[i][0].index;
        for(j=0; index >=0 ; ++j) {
        	iindices[k] = index - 1;
            dvalues[k] = model->SV[i][j].value;
            index = model->SV[i][j+1].index;
            ++k;
        }
        iindptr[i+1] = k;
	}

	return 0;
}

/* get number of nonzero coefficients in support vectors */
npy_intp get_nonzero_SV (struct svm_csr_model *model) {
	int i, j;
	npy_intp count=0;
	for (i=0; i<model->l; ++i) {
		j = 0;
		while (model->SV[i][j].index != -1) {
			++j;
			++count;
		}
	}
	return count;
}


/*
 * Predict using a model, where data is expected to be encoded into a csr matrix.
 */
int csr_copy_predict (npy_intp *data_size, char *data, npy_intp *index_size,
		char *index, npy_intp *intptr_size, char *intptr, struct svm_csr_model *model,
		char *dec_values, BlasFunctions *blas_functions) {
    double *t = (double *) dec_values;
    struct svm_csr_node **predict_nodes;
    npy_intp i;

    predict_nodes = csr_to_libsvm((double *) data, (int *) index,
                                  (int *) intptr, intptr_size[0]-1);

    if (predict_nodes == NULL)
        return -1;
    for(i=0; i < intptr_size[0] - 1; ++i) {
        *t = svm_csr_predict(model, predict_nodes[i], blas_functions);
        free(predict_nodes[i]);
        ++t;
    }
    free(predict_nodes);
    return 0;
}

int csr_copy_predict_values (npy_intp *data_size, char *data, npy_intp *index_size,
                char *index, npy_intp *intptr_size, char *intptr, struct svm_csr_model *model,
                char *dec_values, int nr_class, BlasFunctions *blas_functions) {
    struct svm_csr_node **predict_nodes;
    npy_intp i;

    predict_nodes = csr_to_libsvm((double *) data, (int *) index,
                                  (int *) intptr, intptr_size[0]-1);

    if (predict_nodes == NULL)
        return -1;
    for(i=0; i < intptr_size[0] - 1; ++i) {
        svm_csr_predict_values(model, predict_nodes[i],
                               ((double *) dec_values) + i*nr_class,
			       blas_functions);
        free(predict_nodes[i]);
    }
    free(predict_nodes);

    return 0;
}

int csr_copy_predict_proba (npy_intp *data_size, char *data, npy_intp *index_size,
		char *index, npy_intp *intptr_size, char *intptr, struct svm_csr_model *model,
		char *dec_values, BlasFunctions *blas_functions) {

    struct svm_csr_node **predict_nodes;
    npy_intp i;
    int m = model->nr_class;

    predict_nodes = csr_to_libsvm((double *) data, (int *) index,
                                  (int *) intptr, intptr_size[0]-1);

    if (predict_nodes == NULL)
        return -1;
    for(i=0; i < intptr_size[0] - 1; ++i) {
        svm_csr_predict_probability(
		model, predict_nodes[i], ((double *) dec_values) + i*m, blas_functions);
        free(predict_nodes[i]);
    }
    free(predict_nodes);
    return 0;
}


npy_intp get_nr(struct svm_csr_model *model)
{
    return (npy_intp) model->nr_class;
}

void copy_intercept(char *data, struct svm_csr_model *model, npy_intp *dims)
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

void copy_support (char *data, struct svm_csr_model *model)
{
    memcpy (data, model->sv_ind, (model->l) * sizeof(int));
}

/*
 * Some helpers to convert from libsvm sparse data structures
 * model->sv_coef is a double **, whereas data is just a double *,
 * so we have to do some stupid copying.
 */
void copy_sv_coef(char *data, struct svm_csr_model *model)
{
    int i, len = model->nr_class-1;
    double *temp = (double *) data;
    for(i=0; i<len; ++i) {
        memcpy(temp, model->sv_coef[i], sizeof(double) * model->l);
        temp += model->l;
    }
}

/*
 * Get the number of iterations run in optimization
 */
void copy_n_iter(char *data, struct svm_csr_model *model)
{
    const int n_models = MAX(1, model->nr_class * (model->nr_class-1) / 2);
    memcpy(data, model->n_iter, n_models * sizeof(int));
}

/*
 * Get the number of support vectors in a model.
 */
npy_intp get_l(struct svm_csr_model *model)
{
    return (npy_intp) model->l;
}

void copy_nSV(char *data, struct svm_csr_model *model)
{
    if (model->label == NULL) return;
    memcpy(data, model->nSV, model->nr_class * sizeof(int));
}

/*
 * same as above with model->label
 * TODO: merge in the cython layer
 */
void copy_label(char *data, struct svm_csr_model *model)
{
    if (model->label == NULL) return;
    memcpy(data, model->label, model->nr_class * sizeof(int));
}

void copy_probA(char *data, struct svm_csr_model *model, npy_intp * dims)
{
    memcpy(data, model->probA, dims[0] * sizeof(double));
}

void copy_probB(char *data, struct svm_csr_model *model, npy_intp * dims)
{
    memcpy(data, model->probB, dims[0] * sizeof(double));
}


/*
 * Some free routines. Some of them are nontrivial since a lot of
 * sharing happens across objects (they *must* be called in the
 * correct order)
 */
int free_problem(struct svm_csr_problem *problem)
{
    int i;
    if (problem == NULL) return -1;
    for (i=0; i<problem->l; ++i)
        free (problem->x[i]);
    free (problem->x);
    free (problem);
    return 0;
}

int free_model(struct svm_csr_model *model)
{
    /* like svm_free_and_destroy_model, but does not free sv_coef[i] */
    /* We don't free n_iter, since we did not create them in set_model. */
    if (model == NULL) return -1;
    free(model->SV);
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


int free_model_SV(struct svm_csr_model *model)
{
    int i;
    for (i=model->l-1; i>=0; --i) free(model->SV[i]);
    /* svn_destroy_model frees model->SV */
    for (i=0; i < model->nr_class-1 ; ++i) free(model->sv_coef[i]);
    /* svn_destroy_model frees model->sv_coef */
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
