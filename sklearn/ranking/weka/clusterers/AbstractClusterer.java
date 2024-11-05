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
 *    AbstractClusterer.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.clusterers;

import weka.core.Capabilities;
import weka.core.CapabilitiesHandler;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.core.SerializedObject;
import weka.core.Utils;
import weka.core.Capabilities.Capability;

import java.io.Serializable;

/** 
 * Abstract clusterer.
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @version $Revision: 6624 $
 */
public abstract class AbstractClusterer
  implements Clusterer, Cloneable, Serializable, CapabilitiesHandler, RevisionHandler {

  /** for serialization */
  private static final long serialVersionUID = -6099962589663877632L;

  // ===============
  // Public methods.
  // ===============
 
  /**
   * Generates a clusterer. Has to initialize all fields of the clusterer
   * that are not being set via options.
   *
   * @param data set of instances serving as training data 
   * @exception Exception if the clusterer has not been 
   * generated successfully
   */
  public abstract void buildClusterer(Instances data) throws Exception;

  /**
   * Classifies a given instance. Either this or distributionForInstance()
   * needs to be implemented by subclasses.
   *
   * @param instance the instance to be assigned to a cluster
   * @return the number of the assigned cluster as an integer
   * @exception Exception if instance could not be clustered
   * successfully
   */
  public int clusterInstance(Instance instance) throws Exception {

    double [] dist = distributionForInstance(instance);

    if (dist == null) {
      throw new Exception("Null distribution predicted");
    }

    if (Utils.sum(dist) <= 0) {
      throw new Exception("Unable to cluster instance");
    }
    return Utils.maxIndex(dist);
  }

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
  public double[] distributionForInstance(Instance instance) 
    throws Exception {

    double[] d = new double[numberOfClusters()];

    d[clusterInstance(instance)] = 1.0;
    
    return d;
  }

  /**
   * Returns the number of clusters.
   *
   * @return the number of clusters generated for a training dataset.
   * @exception Exception if number of clusters could not be returned
   * successfully
   */
  public abstract int numberOfClusters() throws Exception;

  /**
   * Creates a new instance of a clusterer given it's class name and
   * (optional) arguments to pass to it's setOptions method. If the
   * clusterer implements OptionHandler and the options parameter is
   * non-null, the clusterer will have it's options set.
   *
   * @param clustererName the fully qualified class name of the clusterer
   * @param options an array of options suitable for passing to setOptions. May
   * be null.
   * @return the newly created search object, ready for use.
   * @exception Exception if the clusterer class name is invalid, or the 
   * options supplied are not acceptable to the clusterer.
   */
  public static Clusterer forName(String clustererName,
				  String [] options) throws Exception {
    return (Clusterer)Utils.forName(Clusterer.class,
				    clustererName,
				    options);
  }

  /**
   * Creates a deep copy of the given clusterer using serialization.
   *
   * @param model the clusterer to copy
   * @return a deep copy of the clusterer
   * @exception Exception if an error occurs
   */
  public static Clusterer makeCopy(Clusterer model) throws Exception {
    return (Clusterer) new SerializedObject(model).getObject();
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
  public static Clusterer [] makeCopies(Clusterer model,
					int num) throws Exception {
     if (model == null) {
      throw new Exception("No model clusterer set");
    }
    Clusterer [] clusterers = new Clusterer [num];
    SerializedObject so = new SerializedObject(model);
    for(int i = 0; i < clusterers.length; i++) {
      clusterers[i] = (Clusterer) so.getObject();
    }
    return clusterers;
  }

  /** 
   * Returns the Capabilities of this clusterer. Derived classifiers have to
   * override this method to enable capabilities.
   *
   * @return            the capabilities of this object
   * @see               Capabilities
   */
  public Capabilities getCapabilities() {
    Capabilities 	result;
    
    result = new Capabilities(this);
    result.enableAll();
//    result.enable(Capability.NO_CLASS);
    
    return result;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return            the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6624 $");
  }
  
  /**
   * runs the clusterer instance with the given options.
   * 
   * @param clusterer		the clusterer to run
   * @param options	the commandline options
   */
  public static void runClusterer(Clusterer clusterer, String[] options) {
    try {
      System.out.println(ClusterEvaluation.evaluateClusterer(clusterer, options));
    } 
    catch (Exception e) {
      if (    (e.getMessage() == null)
	   || (    (e.getMessage() != null)
	        && (e.getMessage().indexOf("General options") == -1) ) )
	e.printStackTrace();
      else
	System.err.println(e.getMessage());
    }
  }
}
