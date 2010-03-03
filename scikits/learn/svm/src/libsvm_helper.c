#include <stdlib.h>
#include "svm.h"
#include <stdio.h>
#include <numpy/arrayobject.h>

/*

   libsvm does not expose this structure, thus it is not in the .h. We add it 
   here so that cython can use it.

   There are some redundant calls to malloc in set_problem and set_param since we do
   not want to export the structs in libsvm.pyx. But still we could do both in a single
   call.

*/

#define SD sizeof(double) /* just a shortcut */

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
   Convert matrix to sparse representation suitable for libsvm. x is
   expected to be an array of length nrow*ncol.

   Typically the matrix will be dense, so we speed up the routine for
   this case. We create a temporary array temp that collects non-zero
   elements and after we just memcpy that to the proper array.

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

    for (i=0; i<dims[0]; i++) {
        T = temp; /* reset stack to start of array */

        for (j=0; j<dims[1]; j++) {
            if ((T->value = *x) != 0) {
                T->index = j+1;
                ++T; /* go to to next struct*/
            }
            ++x; /* whatever happens, go to the next element of the array*/
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
    problem->l = dims[0];
    problem->y = (double *) Y;
    problem->x = dense_to_sparse((double *) X, dims);
    if (problem->x == NULL) { 
        free(problem);
        return NULL;
    }
    return problem;
}

npy_intp get_l(struct svm_model *model)
{
    return model->l;
}

npy_intp get_nr(struct svm_model *model)
{
    return model->nr_class;
}

/* some helpers to convert from libsvm sparse data structures */
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

void copy_rho(char *data, struct svm_model *model, npy_intp *strides)
{
    memcpy(data, model->rho, strides[0]);
}

/* this is a bit more complex since SV are stored as sparse
   structures, so we have to do the conversion on the fly and also
   iterate fast over data             

   XXX: assigning test = strides[1] fails because these
   strides happen to be zero on 64 bits (I have no clue why)

                          */
void copy_SV(char *data, struct svm_model *model, npy_int *strides)
{
    register int i, j, k;
    char *t = data;
    npy_intp n = (npy_intp) model->l;
    int step = sizeof(double);
    for (i=0; i<n; i++) {
        j = 0;
        while ((k = model->SV[i][j].index) != -1) {
            * ((double *) (t + (k-1)*step)) = model->SV[i][j].value;
            j++;
        }
        t += strides[0];
    }
}

/* we must also free the nodes that we created in the call to dense_to_sparse */
int copy_predict(char *train, struct svm_model *model, npy_intp *train_dims,
                 char *dec_values)
{
    double *t = (double *) dec_values;
    register int i, n;
    n = train_dims[0];
    struct svm_node **train_nodes;
    if ((train_nodes = dense_to_sparse((double *) train, train_dims)) == NULL)
        return -1;
    for(i=0; i<n; ++i) {
        *t = svm_predict(model, train_nodes[i]);
        free(train_nodes[i]);
        ++t;
    }
    return 0;
}

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

int free_param(struct svm_parameter *param)
{
    if (param == NULL) return -1;
    free(param);
    return 0;
}
       
                            
