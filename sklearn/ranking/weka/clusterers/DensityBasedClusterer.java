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
 *    DensityBasedClusterer.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.clusterers;

import weka.core.Instance;

/**
 * Interface for clusterers that can estimate the density for a given instance.
 * Implementations will typically extend AbstractDensityBasedClusterer.
 *
 * @author   Mark Hall (mhall@cs.waikato.ac.nz)
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version  $Revision: 5987 $
 */
public interface DensityBasedClusterer extends Clusterer {

  /**
   * Returns the prior probability of each cluster.
   *
   * @return the prior probability for each cluster
   * @exception Exception if priors could not be 
   * returned successfully
   */
  double[] clusterPriors() throws Exception;

  /**
   * Computes the log of the conditional density (per cluster) for a given instance.
   * 
   * @param instance the instance to compute the density for
   * @return an array containing the estimated densities
   * @exception Exception if the density could not be computed
   * successfully
   */
  double[] logDensityPerClusterForInstance(Instance instance) throws Exception;

  /**
   * Computes the density for a given instance.
   * 
   * @param instance the instance to compute the density for
   * @return the density.
   * @exception Exception if the density could not be computed successfully
   */
  double logDensityForInstance(Instance instance) throws Exception;

  /** 
   * Returns the logs of the joint densities for a given instance.
   *
   * @param inst the instance 
   * @return the array of values
   * @exception Exception if values could not be computed
   */
  double[] logJointDensitiesForInstance(Instance inst) throws Exception;
}
