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
 *    Clusterer.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.clusterers;

import weka.core.Capabilities;
import weka.core.Instance;
import weka.core.Instances;

/**
 * Interface for clusterers. Clients will typically extend either
 * AbstractClusterer or AbstractDensityBasedClusterer.
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @version $Revision: 5987 $
 */
public interface Clusterer {

  /**
   * Generates a clusterer. Has to initialize all fields of the clusterer
   * that are not being set via options.
   *
   * @param data set of instances serving as training data 
   * @exception Exception if the clusterer has not been 
   * generated successfully
   */
  void buildClusterer(Instances data) throws Exception;

  /**
   * Classifies a given instance. Either this or distributionForInstance()
   * needs to be implemented by subclasses.
   *
   * @param instance the instance to be assigned to a cluster
   * @return the number of the assigned cluster as an integer
   * @exception Exception if instance could not be clustered
   * successfully
   */
  int clusterInstance(Instance instance) throws Exception;

  /**
   * Predicts the cluster memberships for a given instance.  Either
   * this or clusterInstance() needs to be implemented by subclasses.
   *
   * @param instance the instance to be assigned a cluster.
   * @return an array containing the estimated membership 
   * probabilities of the test instance in each cluster (this 
   * should sum to at most 1)
   * @exception Exception if distribution could not be 
   * computed successfully 
   */
  public double[] distributionForInstance(Instance instance) throws Exception;

  /**
   * Returns the number of clusters.
   *
   * @return the number of clusters generated for a training dataset.
   * @exception Exception if number of clusters could not be returned
   * successfully
   */
  int numberOfClusters() throws Exception;

  /** 
   * Returns the Capabilities of this clusterer. Derived classifiers have to
   * override this method to enable capabilities.
   *
   * @return            the capabilities of this object
   * @see               Capabilities
   */
  public Capabilities getCapabilities();
  
}
