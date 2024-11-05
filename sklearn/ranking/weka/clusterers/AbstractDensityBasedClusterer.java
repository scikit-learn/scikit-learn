/*
 *    This program is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program; if not, write to the Free Software
 *    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/*
 *    AbstractDensityBasedClusterer.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.clusterers;

import weka.core.Instance;
import weka.core.SerializedObject;
import weka.core.Utils;

/** 
 * Abstract clustering model that produces (for each test instance)
 * an estimate of the membership in each cluster 
 * (ie. a probability distribution).
 *
 * @author   Mark Hall (mhall@cs.waikato.ac.nz)
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version  $Revision: 1.1 $
 */
public abstract class AbstractDensityBasedClusterer
  extends AbstractClusterer implements DensityBasedClusterer {

  /** for serialization. */
  private static final long serialVersionUID = -5950728041704213845L;

  // ===============
  // Public methods.
  // ===============

  /**
   * Returns the prior probability of each cluster.
   *
   * @return the prior probability for each cluster
   * @exception Exception if priors could not be 
   * returned successfully
   */
  public abstract double[] clusterPriors() 
    throws Exception;

  /**
   * Computes the log of the conditional density (per cluster) for a given instance.
   * 
   * @param instance the instance to compute the density for
   * @return an array containing the estimated densities
   * @exception Exception if the density could not be computed
   * successfully
   */
  public abstract double[] logDensityPerClusterForInstance(Instance instance) 
    throws Exception;

  /**
   * Computes the density for a given instance.
   * 
   * @param instance the instance to compute the density for
   * @return the density.
   * @exception Exception if the density could not be computed successfully
   */
  public double logDensityForInstance(Instance instance) throws Exception {

    double[] a = logJointDensitiesForInstance(instance);
    double max = a[Utils.maxIndex(a)];
    double sum = 0.0;

    for(int i = 0; i < a.length; i++) {
      sum += Math.exp(a[i] - max);
    }

    return max + Math.log(sum);
  }

  /**
   * Returns the cluster probability distribution for an instance.
   *
   * @param instance the instance to be clustered
   * @return the probability distribution
   * @throws Exception if computation fails
   */  
  public double[] distributionForInstance(Instance instance) throws Exception {
    
    return Utils.logs2probs(logJointDensitiesForInstance(instance));
  }

  /** 
   * Returns the logs of the joint densities for a given instance.
   *
   * @param inst the instance 
   * @return the array of values
   * @exception Exception if values could not be computed
   */
  public double[] logJointDensitiesForInstance(Instance inst)
    throws Exception {

    double[] weights = logDensityPerClusterForInstance(inst);
    double[] priors = clusterPriors();

    for (int i = 0; i < weights.length; i++) {
      if (priors[i] > 0) {
	weights[i] += Math.log(priors[i]);
      } else {
	throw new IllegalArgumentException("Cluster empty!");
      }
    }
    return weights;
  }

  /**
   * Creates copies of the current clusterer. Note that this method
   * now uses Serialization to perform a deep copy, so the Clusterer
   * object must be fully Serializable. Any currently built model will
   * now be copied as well.
   *
   * @param model an example clusterer to copy
   * @param num the number of clusterer copies to create.
   * @return an array of clusterers.
   * @exception Exception if an error occurs 
   */
  public static DensityBasedClusterer [] makeCopies(DensityBasedClusterer model,
						    int num) throws Exception {
     if (model == null) {
      throw new Exception("No model clusterer set");
    }
    DensityBasedClusterer [] clusterers = new DensityBasedClusterer [num];
    SerializedObject so = new SerializedObject(model);
    for(int i = 0; i < clusterers.length; i++) {
      clusterers[i] = (DensityBasedClusterer) so.getObject();
    }
    return clusterers;
  }
}
