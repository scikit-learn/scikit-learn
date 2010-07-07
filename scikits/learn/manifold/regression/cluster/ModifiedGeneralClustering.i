// -*-c++-*-

// Matthieu Brucher
// Last Change : 2008-04-07 14:44

%define DOCSTRING
"
Some clustering algorithms
"
%enddef

%include "std_vector.i"
// Instantiate templates used by example
namespace std {
  %template(LongVector) vector<long>;
}

%{
#define SWIG_FILE_WITH_INIT
%}
%include "numpy.i"
%init %{
import_array();
%}
%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(double* matrix, int rows, int cols)};

%module (package="modified_general_clustering", docstring=DOCSTRING) modified_general_clustering
%{
#include "ModifiedGeneralClustering.h"
#include <matrix/matrix_lib.h>
#include <matrix/pointer_matrix.h>

class GeneralCluster
{
  public:
  GeneralCluster(float dissimilarity_min, float dissimilarity_max)
  {
    clustering = new Graph::Clustering::ModifiedGeneralClustering<std::binder2nd<std::less<float> >, std::binder2nd<std::greater<float> > >(std::bind2nd(std::less<float>(), dissimilarity_min), std::bind2nd(std::greater<float>(), dissimilarity_max));
  }

  ~GeneralCluster()
  {
    delete clustering;
  }

  std::vector<long> process(double* matrix, int rows, int cols)
  {
    return clustering->process(Matrix::PointerMatrix<double>(rows, cols, matrix));
  }

  Graph::Clustering::ModifiedGeneralClustering<std::binder2nd<std::less<float> >, std::binder2nd<std::greater<float> > >* clustering;
};
%}

class GeneralCluster
{
public:
  GeneralCluster(float dissimilarity_min, float dissimilarity_max);
  ~GeneralCluster();
  std::vector<long> process(double* matrix, int rows, int cols);
};
