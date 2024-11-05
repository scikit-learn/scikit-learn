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
 *    Classifier.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers;

import weka.core.Attribute;
import weka.core.Capabilities;
import weka.core.CapabilitiesHandler;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.core.SerializedObject;
import weka.core.Utils;

import java.io.Serializable;
import java.util.Enumeration;
import java.util.Set;
import java.util.Vector;

/**
 * Classifier interface. All schemes for numeric or nominal prediction in
 * Weka implement this interface. Note that a classifier MUST either implement
 * distributionForInstance() or classifyInstance().
 *
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 6041 $
 */
public interface Classifier {

  /**
   * Generates a classifier. Must initialize all fields of the classifier
   * that are not being set via options (ie. multiple calls of buildClassifier
   * must always lead to the same result). Must not change the dataset
   * in any way.
   *
   * @param data set of instances serving as training data
   * @exception Exception if the classifier has not been
   * generated successfully
   */
  public abstract void buildClassifier(Instances data) throws Exception;

  /**
   * Classifies the given test instance. The instance has to belong to a
   * dataset when it's being classified. Note that a classifier MUST
   * implement either this or distributionForInstance().
   *
   * @param instance the instance to be classified
   * @return the predicted most likely class for the instance or
   * Utils.missingValue() if no prediction is made
   * @exception Exception if an error occurred during the prediction
   */
  public double classifyInstance(Instance instance) throws Exception;

  /**
   * Predicts the class memberships for a given instance. If
   * an instance is unclassified, the returned array elements
   * must be all zero. If the class is numeric, the array
   * must consist of only one element, which contains the
   * predicted value. Note that a classifier MUST implement
   * either this or classifyInstance().
   *
   * @param instance the instance to be classified
   * @return an array containing the estimated membership
   * probabilities of the test instance in each class
   * or the numeric prediction
   * @exception Exception if distribution could not be
   * computed successfully
   */
  public double[] distributionForInstance(Instance instance) throws Exception;

  /**
   * Returns the Capabilities of this classifier. Maximally permissive
   * capabilities are allowed by default. Derived classifiers should
   * override this method and first disable all capabilities and then
   * enable just those capabilities that make sense for the scheme.
   *
   * @return            the capabilities of this object
   * @see               Capabilities
   */
  public Capabilities getCapabilities();
}

