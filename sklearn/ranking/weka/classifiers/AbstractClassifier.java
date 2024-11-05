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
 *    AbstractClassifier.java
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
import weka.core.labelranking.PreferenceAttribute;

import java.io.Serializable;
import java.util.Enumeration;
import java.util.Vector;

/**
 * Abstract classifier. All schemes for numeric or nominal prediction in
 * Weka extend this class. Note that a classifier MUST either implement
 * distributionForInstance() or classifyInstance().
 *
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 6624 $
 */
public abstract class AbstractClassifier
  implements Classifier, Cloneable, Serializable, OptionHandler,
             CapabilitiesHandler, RevisionHandler {

  /** for serialization */
  private static final long serialVersionUID = 6502780192411755341L;

  /** Whether the classifier is run in debug mode. */
  protected boolean m_Debug = false;

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
  public double classifyInstance(Instance instance) throws Exception {

    double [] dist = distributionForInstance(instance);
    if (dist == null) {
      throw new Exception("Null distribution predicted");
    }
    switch (instance.classAttribute().type()) {
    case Attribute.NOMINAL:
    case PreferenceAttribute.RANKING:
      double max = 0;
      int maxIndex = 0;

      for (int i = 0; i < dist.length; i++) {
        if (dist[i] > max) {
          maxIndex = i;
          max = dist[i];
        }
      }
      if (max > 0) {
        return maxIndex;
      } else {
        return Utils.missingValue();
      }
    case Attribute.NUMERIC:
      return dist[0];
    default:
      return Utils.missingValue();
    }
  }

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
  public double[] distributionForInstance(Instance instance) throws Exception {

    double[] dist = new double[instance.numClasses()];
    switch (instance.classAttribute().type()) {
      case Attribute.NOMINAL:
      case PreferenceAttribute.RANKING:
        double classification = classifyInstance(instance);
        if (Utils.isMissingValue(classification)) {
          return dist;
        } else {
          dist[(int)classification] = 1.0;
        }
        return dist;
      case Attribute.NUMERIC:
        dist[0] = classifyInstance(instance);
        return dist;
      default:
        return dist;
    }
  }

  /**
   * Creates a new instance of a classifier given it's class name and
   * (optional) arguments to pass to it's setOptions method. If the
   * classifier implements OptionHandler and the options parameter is
   * non-null, the classifier will have it's options set.
   *
   * @param classifierName the fully qualified class name of the classifier
   * @param options an array of options suitable for passing to setOptions. May
   * be null.
   * @return the newly created classifier, ready for use.
   * @exception Exception if the classifier name is invalid, or the options
   * supplied are not acceptable to the classifier
   */
  public static Classifier forName(String classifierName,
      String [] options) throws Exception {

    return ((AbstractClassifier)Utils.forName(Classifier.class,
                                              classifierName,
                                              options));
  }

  /**
   * Creates a deep copy of the given classifier using serialization.
   *
   * @param model the classifier to copy
   * @return a deep copy of the classifier
   * @exception Exception if an error occurs
   */
  public static Classifier makeCopy(Classifier model) throws Exception {

    return (Classifier)new SerializedObject(model).getObject();
  }

  /**
   * Creates a given number of deep copies of the given classifier using serialization.
   *
   * @param model the classifier to copy
   * @param num the number of classifier copies to create.
   * @return an array of classifiers.
   * @exception Exception if an error occurs
   */
  public static Classifier [] makeCopies(Classifier model, int num) throws Exception {

    if (model == null) {
      throw new Exception("No model classifier set");
    }
    Classifier [] classifiers = new Classifier [num];
    SerializedObject so = new SerializedObject(model);
    for(int i = 0; i < classifiers.length; i++) {
      classifiers[i] = (Classifier) so.getObject();
    }
    return classifiers;
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {

    Vector newVector = new Vector(1);

    newVector.addElement(new Option(
          "\tIf set, classifier is run in debug mode and\n"
          + "\tmay output additional info to the console",
          "D", 0, "-D"));
    return newVector.elements();
  }

  /**
   * Parses a given list of options. Valid options are:<p>
   *
   * -D  <br>
   * If set, classifier is run in debug mode and
   * may output additional info to the console.<p>
   *
   * @param options the list of options as an array of strings
   * @exception Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {

    setDebug(Utils.getFlag('D', options));
  }

  /**
   * Gets the current settings of the Classifier.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String [] getOptions() {

    String [] options;
    if (getDebug()) {
      options = new String[1];
      options[0] = "-D";
    } else {
      options = new String[0];
    }
    return options;
  }

  /**
   * Set debugging mode.
   *
   * @param debug true if debug output should be printed
   */
  public void setDebug(boolean debug) {

    m_Debug = debug;
  }

  /**
   * Get whether debugging is turned on.
   *
   * @return true if debugging output is on
   */
  public boolean getDebug() {

    return m_Debug;
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String debugTipText() {
    return "If set to true, classifier may output additional info to " +
      "the console.";
  }

  /**
   * Returns the Capabilities of this classifier. Maximally permissive
   * capabilities are allowed by default. Derived classifiers should
   * override this method and first disable all capabilities and then
   * enable just those capabilities that make sense for the scheme.
   *
   * @return            the capabilities of this object
   * @see               Capabilities
   */
  public Capabilities getCapabilities() {
    Capabilities result = new Capabilities(this);
    result.enableAll();

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
   * runs the classifier instance with the given options.
   *
   * @param classifier		the classifier to run
   * @param options	the commandline options
   */
  public static void runClassifier(Classifier classifier, String[] options) {
    try {
      System.out.println(Evaluation.evaluateModel(classifier, options));
    }
    catch (Exception e) {
      if (    ((e.getMessage() != null) && (e.getMessage().indexOf("General options") == -1))
          || (e.getMessage() == null) )
        e.printStackTrace();
      else
        System.err.println(e.getMessage());
    }
  }
}

