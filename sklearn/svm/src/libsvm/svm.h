#ifndef _LIBSVM_H
#define _LIBSVM_H

#define LIBSVM_VERSION 310

#ifdef __cplusplus
extern "C" {
#endif
#include "_svm_cython_blas_helpers.h"

struct svm_node
{
	int dim;
	int ind; /* index. A bit redundant, but needed if using a
                    precomputed kernel */
	double *values;
};

struct svm_problem
{
	int l;
	double *y;
	struct svm_node *x;
	double *W; /* instance weights */
};


struct svm_csr_node
{
	int index;
	double value;
};

struct svm_csr_problem
{
	int l;
	double *y;
	struct svm_csr_node **x;
        double *W; /* instance weights */
};


enum { C_SVC, NU_SVC, ONE_CLASS, EPSILON_SVR, NU_SVR };	/* svm_type */
enum { LINEAR, POLY, RBF, SIGMOID, PRECOMPUTED }; /* kernel_type */

struct svm_parameter
{
	int svm_type;
	int kernel_type;
	int degree;	/* for poly */
	double gamma;	/* for poly/rbf/sigmoid */
	double coef0;	/* for poly/sigmoid */

	/* these are for training only */
	double cache_size; /* in MB */
	double eps;	/* stopping criteria */
	double C;	/* for C_SVC, EPSILON_SVR and NU_SVR */
	int nr_weight;		/* for C_SVC */
	int *weight_label;	/* for C_SVC */
	double* weight;		/* for C_SVC */
	double nu;	/* for NU_SVC, ONE_CLASS, and NU_SVR */
	double p;	/* for EPSILON_SVR */
	int shrinking;	/* use the shrinking heuristics */
	int probability; /* do probability estimates */
	int max_iter; /* ceiling on Solver runtime */
    int random_seed; /* seed for random number generator */
};

//
// svm_model
//
struct svm_model
{
	struct svm_parameter param;	/* parameter */
	int nr_class;		/* number of classes, = 2 in regression/one class svm */
	int l;			/* total #SV */
	struct svm_node *SV;		/* SVs (SV[l]) */
	double **sv_coef;	/* coefficients for SVs in decision functions (sv_coef[k-1][l]) */

	int *sv_ind;            /* index of support vectors */

	double *rho;		/* constants in decision functions (rho[k*(k-1)/2]) */
	double *probA;		/* pairwise probability information */
	double *probB;

	/* for classification only */

	int *label;		/* label of each class (label[k]) */
	int *nSV;		/* number of SVs for each class (nSV[k]) */
				/* nSV[0] + nSV[1] + ... + nSV[k-1] = l */
	/* XXX */
	int free_sv;		/* 1 if svm_model is created by svm_load_model*/
				/* 0 if svm_model is created by svm_train */
};


struct svm_csr_model
{
	struct svm_parameter param;	/* parameter */
	int nr_class;		/* number of classes, = 2 in regression/one class svm */
	int l;			/* total #SV */
	struct svm_csr_node **SV;		/* SVs (SV[l]) */
	double **sv_coef;	/* coefficients for SVs in decision functions (sv_coef[k-1][l]) */

        int *sv_ind;            /* index of support vectors */

	double *rho;		/* constants in decision functions (rho[k*(k-1)/2]) */
	double *probA;		/* pairwise probability information */
	double *probB;

	/* for classification only */

	int *label;		/* label of each class (label[k]) */
	int *nSV;		/* number of SVs for each class (nSV[k]) */
				/* nSV[0] + nSV[1] + ... + nSV[k-1] = l */
	/* XXX */
	int free_sv;		/* 1 if svm_model is created by svm_load_model*/
				/* 0 if svm_model is created by svm_train */
};


struct svm_model *svm_train(const struct svm_problem *prob, const struct svm_parameter *param, int *status, BlasFunctions *blas_functions);
void svm_cross_validation(const struct svm_problem *prob, const struct svm_parameter *param, int nr_fold, double *target, BlasFunctions *blas_functions);

int svm_save_model(const char *model_file_name, const struct svm_model *model);
struct svm_model *svm_load_model(const char *model_file_name);

int svm_get_svm_type(const struct svm_model *model);
int svm_get_nr_class(const struct svm_model *model);
void svm_get_labels(const struct svm_model *model, int *label);
double svm_get_svr_probability(const struct svm_model *model);

double svm_predict_values(const struct svm_model *model, const struct svm_node *x, double* dec_values, BlasFunctions *blas_functions);
double svm_predict(const struct svm_model *model, const struct svm_node *x, BlasFunctions *blas_functions);
double svm_predict_probability(const struct svm_model *model, const struct svm_node *x, double* prob_estimates, BlasFunctions *blas_functions);

void svm_free_model_content(struct svm_model *model_ptr);
void svm_free_and_destroy_model(struct svm_model **model_ptr_ptr);
void svm_destroy_param(struct svm_parameter *param);

const char *svm_check_parameter(const struct svm_problem *prob, const struct svm_parameter *param);

void svm_set_print_string_function(void (*print_func)(const char *));


/* sparse version */

struct svm_csr_model *svm_csr_train(const struct svm_csr_problem *prob, const struct svm_parameter *param, int *status, BlasFunctions *blas_functions);
void svm_csr_cross_validation(const struct svm_csr_problem *prob, const struct svm_parameter *param, int nr_fold, double *target, BlasFunctions *blas_functions);

int svm_csr_get_svm_type(const struct svm_csr_model *model);
int svm_csr_get_nr_class(const struct svm_csr_model *model);
void svm_csr_get_labels(const struct svm_csr_model *model, int *label);
double svm_csr_get_svr_probability(const struct svm_csr_model *model);

double svm_csr_predict_values(const struct svm_csr_model *model, const struct svm_csr_node *x, double* dec_values, BlasFunctions *blas_functions);
double svm_csr_predict(const struct svm_csr_model *model, const struct svm_csr_node *x, BlasFunctions *blas_functions);
double svm_csr_predict_probability(const struct svm_csr_model *model, const struct svm_csr_node *x, double* prob_estimates, BlasFunctions *blas_functions);

void svm_csr_free_model_content(struct svm_csr_model *model_ptr);
void svm_csr_free_and_destroy_model(struct svm_csr_model **model_ptr_ptr);
void svm_csr_destroy_param(struct svm_parameter *param);

const char *svm_csr_check_parameter(const struct svm_csr_problem *prob, const struct svm_parameter *param);

/* end sparse version */


#ifdef __cplusplus
}
#endif

#endif /* _LIBSVM_H */
