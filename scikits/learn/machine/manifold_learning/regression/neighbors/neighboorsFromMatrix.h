/**
 * \file neighboorsFromMatrix.h
 * Finds the nearest neighboors with the help of a tree
 */

// Matthieu Brucher
// Last Change : 2007-07-17 17:08

#ifndef NEIGHBOORSFROMMATRIX
#define NEIGHBOORSFROMMATRIX

#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>
#include <matrix/matrix_lib.h>

namespace Neighboor
{
  /// This class creates a set of elements in a K-neighboor or Parzen window of an other point with the help of a tree
  template<class MatrixType>
  class NeighboorsFromMatrix
  {
    /// Inner implementation of the tree
    struct Node
    {
      /// The center of the node
      typename MatrixType::Column center;
      /// The dimension of the region described by the node
      std::vector<typename MatrixType::Data_Type> size;
      /// Distance to the corners
      typename MatrixType::Data_Type squareDistance;
      /// Level in the hierarchy
      unsigned long level;
      /// Sub nodes
      std::vector<boost::shared_ptr<Node> > children;
      /// points in the node
      std::vector<unsigned long> points;
    };

    /// Predicate for removing test
    struct NoChildrenPredicate
    {
      /**
        * Tests if the node has children
        * @param node is the node to test
        * @return true if node->children is empty
        */
      bool operator()(const boost::shared_ptr<Node>& node) const
      {
        return node->points.empty();
      }
    };

  public:
    typedef typename MatrixType::Data_Type DataType;
    /**
     * Constructor of a window factory
     * @param points is the set of points used in the neighbooring
     * @param levels indicates the number of levels in the tree
     */
    NeighboorsFromMatrix(const MatrixType& points, unsigned long levels)
      :levels(levels), points(points), tree(new Node)
    {
      Matrix::Matrix<DataType, 0, 2U> minmax(points.height(), 2U);
      for(unsigned long i = 0; i < points.width(); ++i)
      {
        for(unsigned long j = 0; j < points.height(); ++j)
        {
          if(points(j, i) < minmax(j, 0))
            minmax(j, 0) = points(j, i);
          if(points(j, i) > minmax(j, 1))
            minmax(j, 1) = points(j, i);
        }
      }
      tree->size.resize(points.height());
      for(unsigned long i = 0; i < points.height(); ++i)
      {
        tree->size[i] = (minmax(i, 1) - minmax(i, 0)) * 3. / 2.;
      }
      tree->squareDistance = computeSquareDistanceToCorner(tree->size);

      tree->center = (minmax[1] + minmax[0]) * (1. / 2.);
      tree->level = 0U;

      tree->points.resize(points.width());
      for(unsigned long i = 0; i < points.width(); ++i)
      {
        tree->points[i] = i;
      }

      createTree(tree, levels);
    }

    /**
     * Finds the k-neighboors of a point in a set of points
     * @param points is the set of points
     * @param neighboor is the number of neighboors to consider
     * @param newPoint is the point neighboors are searched for
     * @return a multimap of the neighboors
     */
    template<class ColumnVector>
    std::multimap<DataType, unsigned long> kneighboors(const ColumnVector& newPoint, unsigned long neighboor) const
    {
      assert(newPoint.width() == 1U);
      assert(newPoint.height() == points.height());

      std::multimap<DataType, unsigned long> neighboors;
      std::multimap<DataType, Node*> openNode;
      openNode.insert(std::make_pair(norm2(newPoint - tree->center), tree.get()));

      while(!(openNode.empty()))
      {
        typename std::multimap<DataType, Node*>::iterator it = openNode.begin();
        Node* front = it->second;
        openNode.erase(it);

        if(neighboors.size() < neighboor)
          processNode(newPoint, front, openNode, neighboors);
        else
        {
          typename std::multimap<DataType, unsigned long>::iterator last = neighboors.begin();
          std::advance(last, neighboor - 1);

          if((DataTypeTraits<DataType>::sqrt(norm2(newPoint - front->center)) - DataTypeTraits<DataType>::sqrt(front->squareDistance)) < last->first)
          {
            processNode(newPoint, front, openNode, neighboors);
          }
        }
      }

      typename std::multimap<DataType, unsigned long>::iterator iterator = neighboors.begin();
      std::advance(iterator, neighboor);
      neighboors.erase(iterator, neighboors.end());

      return neighboors;
    }

    /**
     * Finds the neighboors in a Parzen window of a point in a set of points
     * @param points is the set of points
     * @param windowSize is the size of the Parzen window
     * @param newPoint is the point neighboors are searched for
     * @return a multimap of the neighboors
     */
    template<class ColumnVector>
    std::multimap<DataType, unsigned long> parzen(const ColumnVector& newPoint, DataType windowSize) const
    {
      assert(newPoint.width() == 1U);
      assert(newPoint.height() == points.height());

      std::multimap<DataType, unsigned long> neighboors;
      std::multimap<DataType, Node*> openNode;
      openNode.insert(std::make_pair(norm2(newPoint - tree->center), tree.get()));

      DataType squareDistance = windowSize * windowSize;

       while(!(openNode.empty()))
      {
        typename std::multimap<DataType, Node*>::iterator it = openNode.begin();
        Node* front = it->second;
        openNode.erase(it);

        if((DataTypeTraits<DataType>::sqrt(norm2(newPoint - front->center)) - DataTypeTraits<DataType>::sqrt(front->squareDistance)) < squareDistance)
        {
          processNode(newPoint, front, openNode, neighboors);
        }
      }
      typename std::multimap<DataType, unsigned long>::iterator it = neighboors.upper_bound(windowSize);
      neighboors.erase(it, neighboors.end());
      return neighboors;
    }

  public:
    /**
      * Creates the tree
      * @param tree is the subtree that will be processed
      * @param levels is the number of levels that relains to be done
      */
    template<class Container>
    void createTree(Container& tree, unsigned long levels)
    {
      DataType meanSize = std::accumulate(tree->size.begin(), tree->size.end(), DataTypeTraits<DataType>::zero(points(0,0)), std::plus<DataType>()) / tree->size.size();

      std::vector<unsigned long> dimensions;
      for(unsigned long i = 0; i < tree->size.size(); ++i)
        if(tree->size[i] > (meanSize * 2.f / 3.f))
          dimensions.push_back(i);

      std::vector<DataType> newSize = tree->size;
      for(unsigned long i = 0; i < dimensions.size(); ++i)
        newSize[dimensions[i]] *= 1.f / 2.f;
      DataType distanceToCorner = computeSquareDistanceToCorner(newSize);

      if(levels != 0U && !tree->points.empty())
      {
        tree->children.resize(1 << dimensions.size());
        for(unsigned long i = 0; i < (1U << dimensions.size()); ++i)
        {
          tree->children[i].reset(new Node);

          tree->children[i]->center = tree->center;
          for(unsigned long j = 0; j < dimensions.size(); ++j)
          {
            tree->children[i]->center(dimensions[j]) += ((i & (1U << j))==(1U << j)? newSize[dimensions[j]]:-newSize[dimensions[j]]) * (1.f / 2.f);
          }

          tree->children[i]->size = newSize;
          tree->children[i]->squareDistance = distanceToCorner;
          tree->children[i]->level = tree->level + 1;

          populatePoints(tree->points, tree->children[i]->points, tree->children[i]->center, newSize);

          createTree(tree->children[i], levels - 1);
        }
        tree->children.erase(std::remove_if(tree->children.begin(), tree->children.end(), NoChildrenPredicate()), tree->children.end());
      }
    }

    /**
      * Computes the square distance from the center to a corner
      * @param size is the size of the parallelepipede
      * @return the square of the half size norm2
      */
    DataType computeSquareDistanceToCorner(const std::vector<DataType>& size) const
    {
      DataType distance = DataTypeTraits<DataType>::one(points(0, 0));
      for(unsigned long i = 0; i < points.height(); ++i)
      {
        distance += size[i] * size[i];
      }
      return distance * (1.f / 4.f);
    }

    /**
      * Populates a vector based on a set of points that can be in the vector, the center of the space and the size of the space in which the populated points must be
      * @param indicePoints is the indices of the points to consider
      * @param toPopulatePoints is a vector of indices of points that must be populated
      * @param center is the center of the space to consider
      * @param size is the size of the considered space
      */
    void populatePoints(const std::vector<unsigned long>& indicePoints, std::vector<unsigned long>& toPopulatePoints, const typename MatrixType::Column& center, const std::vector<DataType>& size)
    {
      toPopulatePoints.clear();
      typename MatrixType::Column sizes(size.size(), 1U);
      std::copy(size.begin(), size.end(), sizes.begin());
      sizes *= 1.f / 2.f;

      for(unsigned long i = 0; i < indicePoints.size(); ++i)
      {
        if(MatrixType::Column::OwnTraits::absolute(points[indicePoints[i]] - center) < sizes)
        {
          toPopulatePoints.push_back(indicePoints[i]);
        }
      }
    }

    /**
      * Processes a node and populates nodeMap or neighboors
      * @param newPoint is the point for which neighboors must be found
      * @param node is the node to process
      * @param nodesMap is the list of nodes taht remains to be processed
      * @param neighboors is the list of neighboors that are near enough the new point
      */
    template<class ColumnVector>
    void processNode(const ColumnVector& newPoint, Node* node, std::multimap<DataType, Node*>& nodesMap, std::multimap<DataType, unsigned long>& neighboors) const
    {
      if(node->children.empty())
      {
        for(std::vector<unsigned long>::const_iterator it = node->points.begin(); it != node->points.end(); ++it)
        {
          neighboors.insert(std::make_pair(DataTypeTraits<DataType>::sqrt(norm2(newPoint - points[*it])), *it));
        }
      }
      else
      {
        for(typename std::vector<boost::shared_ptr<Node> >::const_iterator it = node->children.begin(); it != node->children.end(); ++it)
        {
          nodesMap.insert(std::make_pair(norm2(newPoint - (*it)->center), it->get()));
        }
      }
    }

    /// The number of levels in the tree
    unsigned long levels;
    /// The points used for the neighbooring
    MatrixType points;
    /// Pointer to the inner tree
    boost::shared_ptr<Node> tree;
  };
}
#endif
