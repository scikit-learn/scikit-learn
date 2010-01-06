
#include <matrix/matrix_lib.h>
#include "neighborsFromMatrix.h"
#include "neighbors.h"

Neighbor::NeighborsFromMatrix<Matrix::Matrix<double> >* allocate_neighborer(double* data, unsigned long dim1, unsigned long dim2, unsigned long levels)
{
  Neighbor::NeighborsFromMatrix<Matrix::Matrix<double> >* temp = new Neighbor::NeighborsFromMatrix<Matrix::Matrix<double> >(Matrix::Matrix<double>(dim2, dim1, data), levels);
  return temp;
}

int delete_neighborer(Neighbor::NeighborsFromMatrix<Matrix::Matrix<double> >* tree)
{
  delete tree;
  return 0;
}

unsigned long find_kneighbors(Neighbor::NeighborsFromMatrix<Matrix::Matrix<double> >* tree, double* data, unsigned long dim1, unsigned long neighbors, numpy_callback callback)
{
  std::multimap<double, unsigned long> results = tree->kneighbors(Matrix::Matrix<double>(dim1, 1, data), neighbors);

  for(std::multimap<double, unsigned long>::const_iterator it = results.begin(); it != results.end(); ++it)
    callback(it->first, it->second);

  return results.size();
}

unsigned long find_parzen(Neighbor::NeighborsFromMatrix<Matrix::Matrix<double> >* tree, double* data, unsigned long dim1, double windowSize, numpy_callback callback)
{
  std::multimap<double, unsigned long> results = tree->parzen(Matrix::Matrix<double>(dim1, 1, data), windowSize);

  for(std::multimap<double, unsigned long>::const_iterator it = results.begin(); it != results.end(); ++it)
    if(it->first <= windowSize * windowSize)
      callback(it->first, it->second);

  return results.size();
}
