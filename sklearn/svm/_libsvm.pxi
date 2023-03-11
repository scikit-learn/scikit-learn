################################################################################
# Includes
cdef extern from "_svm_cython_blas_helpers.h":
    ctypedef double (*dot_func)(int, double*, int, double*, int)
    cdef struct BlasFunctions:
        dot_func dot


cdef extern from "svm.h":
    cdef struct svm_node
    cdef struct svm_model
    cdef struct svm_parameter:
        int svm_type
        int kernel_type
        int degree	# for poly
        double gamma	# for poly/rbf/sigmoid
        double coef0	# for poly/sigmoid

        # these are for training only
        double cache_size # in MB
        double eps	# stopping criteria
        double C	# for C_SVC, EPSILON_SVR and NU_SVR
        int nr_weight		# for C_SVC
        int *weight_label	# for C_SVC
        double* weight		# for C_SVC
        double nu	# for NU_SVC, ONE_CLASS, and NU_SVR
        double p	# for EPSILON_SVR
        int shrinking	# use the shrinking heuristics
        int probability # do probability estimates
        int max_iter  # ceiling on Solver runtime
        int random_seed  # seed for random generator in probability estimation

    cdef struct svm_problem:
        int l
        double *y
        svm_node *x
        double *W # instance weights

    char *svm_check_parameter(svm_problem *prob, svm_parameter *param)
    svm_model *svm_train(
        svm_problem *prob,
        svm_parameter *param,
        int *status,
        BlasFunctions *blas_functions,
    ) nogil
    void svm_free_and_destroy_model(svm_model** model_ptr_ptr)
    void svm_cross_validation(
        svm_problem *prob,
        svm_parameter *param,
        int nr_fold,
        double *target,
        BlasFunctions *blas_functions,
    ) nogil


cdef extern from "libsvm_helper.c":
    # this file contains methods for accessing libsvm 'hidden' fields
    svm_node **dense_to_sparse(char *, cnp.npy_intp *)
    void set_parameter(
        svm_parameter *param,
        int svm_type,
        int kernel_type,
        int degree,
        double gamma,
        double coef0,
        double nu,
        double cache_size,
        double C,
        double eps,
        double p,
        int shrinking,
        int probability,
        int nr_weight,
        char *weight_label,
        char *weight,
        int max_iter,
        int random_seed,
    )
    void set_problem(
        svm_problem *problem,
        char *X,
        char *Y,
        char *sample_weight,
        cnp.npy_intp *dims,
        int kernel_type,
    )
    svm_model *set_model(
        svm_parameter *param,
        int nr_class,
        char *SV,
        cnp.npy_intp *SV_dims,
        char *support,
        cnp.npy_intp *support_dims,
        cnp.npy_intp *sv_coef_strides,
        char *sv_coef,
        char *rho,
        char *nSV,
        char *probA,
        char *probB,
    )
    void copy_sv_coef(char *data, svm_model *model)
    void copy_n_iter(char *data, svm_model *model)
    void copy_intercept(char *data, svm_model *model, cnp.npy_intp *dims)
    void copy_SV(char *data, svm_model *model, cnp.npy_intp *dims)
    int copy_support(char *data, svm_model *model)
    int copy_predict(
        char *predict,
        svm_model *model,
        cnp.npy_intp *predict_dims,
        char *dec_values,
        BlasFunctions *blas_functions,
    ) nogil
    int copy_predict_proba(
        char *predict,
        svm_model *model,
        cnp.npy_intp *predict_dims,
        char *dec_values,
        BlasFunctions *blas_functions,
    ) nogil
    int copy_predict_values(
        char *predict,
        svm_model *model,
        cnp.npy_intp *predict_dims,
        char *dec_values,
        int nr_class,
        BlasFunctions *blas_functions,
    ) nogil
    void copy_nSV(char *data, svm_model *model)
    void copy_probA(char *data, svm_model *model, cnp.npy_intp *dims)
    void copy_probB(char *data, svm_model *model, cnp.npy_intp *dims)
    cnp.npy_intp get_l(svm_model *model)
    cnp.npy_intp get_nr(svm_model *model)
    int free_problem(svm_problem *)
    int free_model(svm_model *model)
    void set_verbosity(int)
