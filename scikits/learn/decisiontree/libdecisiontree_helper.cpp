//#include <stdlib.h>
//#include <numpy/arrayobject.h>
#include "Node.h"
#include "Object.h"

/*
 * Helper methods for libdecisiontree bindings.
 *
 * License: New BSD.
 *
 * Author: 2011 Noel Dawe <Noel.Dawe@cern.ch>
 */

/* 
 * Initialize the root node prior to fitting
 */
void init_root(char* X, char* Y, char* weights, npy_intp* X_dims, Node* root)
{
    double* x = (double *) X;
    double* y = (double *) Y;
    double* w = (double *) weights;
    npy_intp i;
    Object* obj;
    for(i=0; i<X_dims[0]; ++i)
    {
        obj = new Object();
        obj->attrs = &x[X_dims[0] * i];
        obj->weight = w[i];
        obj->dim =X_dims[1]; 
        if (y[i]==1)
        {
            obj->label = SIGNAL;
            root->add_signal(obj);
        }
        else
        {
            obj->label = BACKGROUND;
            root->add_background(obj);
        }
    }
}

/* 
 * Predict using model.
 *
 *  It will return -1 if we run out of memory.
 */
void copy_predict(char* predict, Node* root, npy_intp* predict_dims, char* dec_values)
{
    double* t = (double *) dec_values;
    double* p = (double *) predict;
    npy_intp i;
    if (root == NULL)
        return;
    for(i=0; i<predict_dims[0]; ++i)
        *(t++) = root->predict(&p[i*predict_dims[1]]);
}
