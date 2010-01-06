
#include <matrix/matrix_lib.h>
#include "neighboorsFromMatrix.h"
#include "neighboors.h"

Neighboor::NeighboorsFromMatrix<Matrix::Matrix<double> >* allocate_neighbooorer(double* data, unsigned long dim1, unsigned long dim2, unsigned long levels)
{
  Neighboor::NeighboorsFromMatrix<Matrix::Matrix<double> >* temp = new Neighboor::NeighboorsFromMatrix<Matrix::Matrix<double> >(Matrix::Matrix<double>(dim2, dim1, data), levels);
  return temp;
}

int delete_neighboorer(Neighboor::NeighboorsFromMatrix<Matrix::Matrix<double> >* tree)
{
  delete tree;
  return 0;
}

unsigned long find_kneighboors(Neighboor::NeighboorsFromMatrix<Matrix::Matrix<double> >* tree, double* data, unsigned long dim1, unsigned long neighboors, numpy_callback callback)
{
  std::multimap<double, unsigned long> results = tree->kneighboors(Matrix::Matrix<double>(dim1, 1, data), neighboors);

  for(std::multimap<double, unsigned long>::const_iterator it = results.begin(); it != results.end(); ++it)
    callback(it->first, it->second);

  return results.size();
}

unsigned long find_parzen(Neighboor::NeighboorsFromMatrix<Matrix::Matrix<double> >* tree, double* data, unsigned long dim1, double windowSize, numpy_callback callback)
{
  std::multimap<double, unsigned long> results = tree->parzen(Matrix::Matrix<double>(dim1, 1, data), windowSize);

  for(std::multimap<double, unsigned long>::const_iterator it = results.begin(); it != results.end(); ++it)
    if(it->first <= windowSize * windowSize)
      callback(it->first, it->second);

  return results.size();
}
