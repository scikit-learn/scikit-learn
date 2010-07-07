/**
 * \file ModifiedGeneralClustering.h
 * Clustering with the General Algorithm described in doi:10.1016/j.jcss.2004.10.012
 * Uses a general matrix type containing the complete correlation between polongs in the graph
 * Once the general algorithm is complete, a last pass is made without checking distances between clusters
 */

#ifndef MODIFIEDGENERALCLUSTERING
#define MODIFIEDGENERALCLUSTERING

#include "clusterInfos.h"

#include <vector>
#include <cmath>


namespace Graph
{
  namespace Clustering
  {
    // std::isnan is a C99 extension
    // some compileres do not implement it
    inline bool isnan(long double element)
    {
      return element != element;
    }

    /// Structure encapsulating the General Cluster algorithm, modified for our purpose
    template<class EnoughDistance, class LittleDistance>
    class ModifiedGeneralClustering
    {
    public:
      /**
       * Constructs the clustering object
       * @param enoughDistance indicates if a distance is big enough
       * @param littleDistance indicates if a distance is small enough
       */
      ModifiedGeneralClustering(EnoughDistance enoughDistance, LittleDistance littleDistance)
        :enoughDistance(enoughDistance), littleDistance(littleDistance)
      {
      }
      /**
       * Make a cluster with each point and the best correlated other polongs
       * @param correlations is a matrix or a pseudo matrix of correlation between points
       * @return a vector with cluster associated to each polongs or -1 for singletons
       * \pre correlations must be symmetric and its diagonal must be one
       */
      template<class CorrelationMatrix>
      std::vector<long> process(const CorrelationMatrix& correlations)
      {
        std::vector<long> clusters(correlations.width(), noCluster);
        long cluster = 0;

        for(unsigned long i = 0; i < correlations.width(); ++i)
        {
          if(clusters[i] == noCluster)
          {
            for(unsigned long j = 0; (j < correlations.width())/* && (clusters[i] == noCluster)*/; ++j)
            {
              if(enoughDistance(correlations(i, j)) && clusters[j] == noCluster)
              {
                makeCluster(i, clusters, correlations, cluster);
                ++cluster;
                break;
              }
            }
          }
        }
        for(unsigned long i = 0; i < correlations.width(); ++i)
        {
          if(clusters[i] == noCluster)
          {
            makeCluster(i, clusters, correlations, cluster);
            ++cluster;
          }
        }
        return clusters;
      }

      /**
       * Associates a point to a given cluster if the correlation is good
       * @param polong is the point to which the correlation is taken
       * @param clusters is the vector of associated clusters that is modified
       * @param correlations is a matrix of correlations
       * @param cluster is the cluster that is analyzed
       * \pre correlations(point, point) should be one if point is in its cluster
       */
      template<class CorrelationMatrix>
      void makeCluster(unsigned long point, std::vector<long>& clusters, const CorrelationMatrix& correlations, long cluster)
      {
        for(unsigned long nbCluster = 0; nbCluster < clusters.size(); ++nbCluster)
        {
          typename CorrelationMatrix::Data_Type value = correlations(point, nbCluster);
          if(clusters[nbCluster] == noCluster && (littleDistance(value) || isnan(value)))
            clusters[nbCluster] = cluster;
        }
      }

    private:
      /// Indicates if the distance between two points is big enough
      EnoughDistance enoughDistance;
      /// Indicates if the distance between two points is small enough
      LittleDistance littleDistance;
    };

    /**
     * Returns a general cluster algorithm ready to be used
     * @param enoughDistance is a functor indicating if the distance between two points is big enough
     * @param littleDistance is a functor indicating if the distance between two points is small enough
     * @return the generated instance of General Clustering
     */
    template<class EnoughDistance, class LittleDistance>
    ModifiedGeneralClustering<EnoughDistance, LittleDistance> makeModifiedGeneralClusterAlgorithm(EnoughDistance enoughDistance, LittleDistance littleDistance)
    {
      return ModifiedGeneralClustering<EnoughDistance, LittleDistance>(enoughDistance, littleDistance);
    }
  }
}
#endif
