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
 *    NaiveBayes.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.bayes;

import weka.classifiers.Classifier;
import weka.classifiers.AbstractClassifier;
import weka.core.Attribute;
import weka.core.Capabilities;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformationHandler;
import weka.core.Utils;
import weka.core.WeightedInstancesHandler;
import weka.core.Capabilities.Capability;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;
import weka.core.labelranking.PreferenceAttribute;
import weka.estimators.DiscreteEstimator;
import weka.estimators.Estimator;
import weka.estimators.KernelEstimator;
import weka.estimators.NormalEstimator;

import java.util.Enumeration;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * Class for a Naive Bayes classifier using estimator classes. Numeric estimator precision values are chosen based on analysis of the  training data. For this reason, the classifier is not an UpdateableClassifier (which in typical usage are initialized with zero training instances) -- if you need the UpdateableClassifier functionality, use the NaiveBayesUpdateable classifier. The NaiveBayesUpdateable classifier will  use a default precision of 0.1 for numeric attributes when buildClassifier is called with zero training instances.<br/>
 * <br/>
 * For more information on Naive Bayes classifiers, see<br/>
 * <br/>
 * George H. John, Pat Langley: Estimating Continuous Distributions in Bayesian Classifiers. In: Eleventh Conference on Uncertainty in Artificial Intelligence, San Mateo, 338-345, 1995.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;inproceedings{John1995,
 *    address = {San Mateo},
 *    author = {George H. John and Pat Langley},
 *    booktitle = {Eleventh Conference on Uncertainty in Artificial Intelligence},
 *    pages = {338-345},
 *    publisher = {Morgan Kaufmann},
 *    title = {Estimating Continuous Distributions in Bayesian Classifiers},
 *    year = {1995}
 * }
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -K
 *  Use kernel density estimator rather than normal
 *  distribution for numeric attributes</pre>
 * 
 * <pre> -D
 *  Use supervised discretization to process numeric attributes
 * </pre>
 *
 * <pre> -O
 *  Display model in old format (good when there are many classes)
 * </pre>
 * 
 <!-- options-end -->
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version $Revision: 5928 $
 */
public class NaiveBayes extends AbstractClassifier 
implements OptionHandler, WeightedInstancesHandler, 
           TechnicalInformationHandler {

  /** for serialization */
  static final long serialVersionUID = 5995231201785697655L;

  /** The attribute estimators. */
  protected Estimator [][] m_Distributions;

  /** The class estimator. */
  protected Estimator m_ClassDistribution;

  /**
   * Whether to use kernel density estimator rather than normal distribution
   * for numeric attributes
   */
  protected boolean m_UseKernelEstimator = false;

  /**
   * Whether to use discretization than normal distribution
   * for numeric attributes
   */
  protected boolean m_UseDiscretization = false;

  /** The number of classes (or 1 for numeric class) */
  protected int m_NumClasses;

  /**
   * The dataset header for the purposes of printing out a semi-intelligible 
   * model 
   */
  protected Instances m_Instances;

  /*** The precision parameter used for numeric attributes */
  protected static final double DEFAULT_NUM_PRECISION = 0.01;

  /**
   * The discretization filter.
   */
  protected weka.filters.supervised.attribute.Discretize m_Disc = null;

  protected boolean m_displayModelInOldFormat = false;

  /**
   * Returns a string describing this classifier
   * @return a description of the classifier suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return "Class for a Naive Bayes classifier using estimator classes. Numeric"
      +" estimator precision values are chosen based on analysis of the "
      +" training data. For this reason, the classifier is not an"
      +" UpdateableClassifier (which in typical usage are initialized with zero"
      +" training instances) -- if you need the UpdateableClassifier functionality,"
      +" use the NaiveBayesUpdateable classifier. The NaiveBayesUpdateable"
      +" classifier will  use a default precision of 0.1 for numeric attributes"
      +" when buildClassifier is called with zero training instances.\n\n"
      +"For more information on Naive Bayes classifiers, see\n\n"
      + getTechnicalInformation().toString();
  }

  /**
   * Returns an instance of a TechnicalInformation object, containing 
   * detailed information about the technical background of this class,
   * e.g., paper reference or book this class is based on.
   * 
   * @return the technical information about this class
   */
  public TechnicalInformation getTechnicalInformation() {
    TechnicalInformation 	result;

    result = new TechnicalInformation(Type.INPROCEEDINGS);
    result.setValue(Field.AUTHOR, "George H. John and Pat Langley");
    result.setValue(Field.TITLE, "Estimating Continuous Distributions in Bayesian Classifiers");
    result.setValue(Field.BOOKTITLE, "Eleventh Conference on Uncertainty in Artificial Intelligence");
    result.setValue(Field.YEAR, "1995");
    result.setValue(Field.PAGES, "338-345");
    result.setValue(Field.PUBLISHER, "Morgan Kaufmann");
    result.setValue(Field.ADDRESS, "San Mateo");

    return result;
  }

  /**
   * Returns default capabilities of the classifier.
   *
   * @return      the capabilities of this classifier
   */
  public Capabilities getCapabilities() {
    Capabilities result = super.getCapabilities();
    result.disableAll();

    // attributes
    result.enable(Capability.NOMINAL_ATTRIBUTES);
    result.enable(Capability.NUMERIC_ATTRIBUTES);
    result.enable(Capability.MISSING_VALUES);

    // class
    result.enable(Capability.NOMINAL_CLASS);
    result.enable(Capability.MISSING_CLASS_VALUES);

    // instances
    result.setMinimumNumberInstances(0);

    return result;
  }

  /**
   * Generates the classifier.
   *
   * @param instances set of instances serving as training data 
   * @exception Exception if the classifier has not been generated 
   * successfully
   */
  public void buildClassifier(Instances instances) throws Exception {

    // can classifier handle the data?
    getCapabilities().testWithFail(instances);

    // remove instances with missing class
    instances = new Instances(instances);
    instances.deleteWithMissingClass();

    m_NumClasses = instances.numClasses();

    // Copy the instances
    m_Instances = new Instances(instances);

    // Discretize instances if required
    if (m_UseDiscretization) {
      m_Disc = new weka.filters.supervised.attribute.Discretize();
      m_Disc.setInputFormat(m_Instances);
      m_Instances = weka.filters.Filter.useFilter(m_Instances, m_Disc);
    } else {
      m_Disc = null;
    }

    // Reserve space for the distributions
    m_Distributions = new Estimator[m_Instances.numAttributes() - 1]
      [m_Instances.numClasses()];
    m_ClassDistribution = new DiscreteEstimator(m_Instances.numClasses(), 
                                                true);
    int attIndex = 0;
    Enumeration enu = m_Instances.enumerateAttributes();
    while (enu.hasMoreElements()) {
      Attribute attribute = (Attribute) enu.nextElement();

      // If the attribute is numeric, determine the estimator 
      // numeric precision from differences between adjacent values
      double numPrecision = DEFAULT_NUM_PRECISION;
      if (attribute.type() == Attribute.NUMERIC) {
	m_Instances.sort(attribute);
	if ((m_Instances.numInstances() > 0)
	    && !m_Instances.instance(0).isMissing(attribute)) {
	  double lastVal = m_Instances.instance(0).value(attribute);
	  double currentVal, deltaSum = 0;
	  int distinct = 0;
	  for (int i = 1; i < m_Instances.numInstances(); i++) {
	    Instance currentInst = m_Instances.instance(i);
	    if (currentInst.isMissing(attribute)) {
	      break;
	    }
	    currentVal = currentInst.value(attribute);
	    if (currentVal != lastVal) {
	      deltaSum += currentVal - lastVal;
	      lastVal = currentVal;
	      distinct++;
	    }
	  }
	  if (distinct > 0) {
	    numPrecision = deltaSum / distinct;
	  }
	}
      }


      for (int j = 0; j < m_Instances.numClasses(); j++) {
	switch (attribute.type()) {
	case Attribute.NUMERIC: 
	  if (m_UseKernelEstimator) {
	    m_Distributions[attIndex][j] = 
	      new KernelEstimator(numPrecision);
	  } else {
	    m_Distributions[attIndex][j] = 
	      new NormalEstimator(numPrecision);
	  }
	  break;
	case PreferenceAttribute.RANKING:
	case Attribute.NOMINAL:
	  m_Distributions[attIndex][j] = 
	    new DiscreteEstimator(attribute.numValues(), true);
	  break;
	default:
	  throw new Exception("Attribute type unknown to NaiveBayes");
	}
      }
      attIndex++;
    }

    // Compute counts
    Enumeration enumInsts = m_Instances.enumerateInstances();
    while (enumInsts.hasMoreElements()) {
      Instance instance = 
	(Instance) enumInsts.nextElement();
      updateClassifier(instance);
    }

    // Save space
    m_Instances = new Instances(m_Instances, 0);
  }


  /**
   * Updates the classifier with the given instance.
   *
   * @param instance the new training instance to include in the model 
   * @exception Exception if the instance could not be incorporated in
   * the model.
   */
  public void updateClassifier(Instance instance) throws Exception {

    if (!instance.classIsMissing()) {
      Enumeration enumAtts = m_Instances.enumerateAttributes();
      int attIndex = 0;
      while (enumAtts.hasMoreElements()) {
	Attribute attribute = (Attribute) enumAtts.nextElement();
	if (!instance.isMissing(attribute)) {
	  m_Distributions[attIndex][(int)instance.classValue()].
            addValue(instance.value(attribute), instance.weight());
	}
	attIndex++;
      }
      m_ClassDistribution.addValue(instance.classValue(),
                                   instance.weight());
    }
  }


  /**
   * Calculates the class membership probabilities for the given test 
   * instance.
   *
   * @param instance the instance to be classified
   * @return predicted class probability distribution
   * @exception Exception if there is a problem generating the prediction
   */
  public double [] distributionForInstance(Instance instance) 
    throws Exception { 

    if (m_UseDiscretization) {
      m_Disc.input(instance);
      instance = m_Disc.output();
    }
    double [] probs = new double[m_NumClasses];
    for (int j = 0; j < m_NumClasses; j++) {
      probs[j] = m_ClassDistribution.getProbability(j);
    }
    Enumeration enumAtts = instance.enumerateAttributes();
    int attIndex = 0;
    while (enumAtts.hasMoreElements()) {
      Attribute attribute = (Attribute) enumAtts.nextElement();
      if (!instance.isMissing(attribute)) {
	double temp, max = 0;
	for (int j = 0; j < m_NumClasses; j++) {
	  temp = Math.max(1e-75, Math.pow(m_Distributions[attIndex][j].
                                          getProbability(instance.value(attribute)), 
                                          m_Instances.attribute(attIndex).weight()));
	  probs[j] *= temp;
	  if (probs[j] > max) {
	    max = probs[j];
	  }
	  if (Double.isNaN(probs[j])) {
	    throw new Exception("NaN returned from estimator for attribute "
                                + attribute.name() + ":\n"
                                + m_Distributions[attIndex][j].toString());
	  }
	}
	if ((max > 0) && (max < 1e-75)) { // Danger of probability underflow
	  for (int j = 0; j < m_NumClasses; j++) {
	    probs[j] *= 1e75;
	  }
	}
      }
      attIndex++;
    }

    // Display probabilities
    Utils.normalize(probs);
    return probs;
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {

    Vector newVector = new Vector(3);

    newVector.addElement(
              new Option("\tUse kernel density estimator rather than normal\n"
                         +"\tdistribution for numeric attributes",
                         "K", 0,"-K"));
    newVector.addElement(
              new Option("\tUse supervised discretization to process numeric attributes\n",
                         "D", 0,"-D"));
    
    newVector.addElement(
              new Option("\tDisplay model in old format (good when there are "
                         + "many classes)\n",
                         "O", 0, "-O"));
    
    return newVector.elements();
  }

  /**
   * Parses a given list of options. <p/>
   *
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -K
   *  Use kernel density estimator rather than normal
   *  distribution for numeric attributes</pre>
   * 
   * <pre> -D
   *  Use supervised discretization to process numeric attributes
   * </pre>
   *
   * <pre> -O
   *  Display model in old format (good when there are many classes)
   * </pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @exception Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {

    boolean k = Utils.getFlag('K', options);
    boolean d = Utils.getFlag('D', options);
    if (k && d) {
      throw new IllegalArgumentException("Can't use both kernel density " +
                                         "estimation and discretization!");
    }
    setUseSupervisedDiscretization(d);
    setUseKernelEstimator(k);
    setDisplayModelInOldFormat(Utils.getFlag('O', options));
    Utils.checkForRemainingOptions(options);
  }

  /**
   * Gets the current settings of the classifier.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String [] getOptions() {

    String [] options = new String [3];
    int current = 0;

    if (m_UseKernelEstimator) {
      options[current++] = "-K";
    }

    if (m_UseDiscretization) {
      options[current++] = "-D";
    }

    if (m_displayModelInOldFormat) {
      options[current++] = "-O";
    }

    while (current < options.length) {
      options[current++] = "";
    }
    return options;
  }

  /**
   * Returns a description of the classifier.
   *
   * @return a description of the classifier as a string.
   */
  public String toString() {
    if (m_displayModelInOldFormat) {
      return toStringOriginal();
    }

    StringBuffer temp = new StringBuffer();
    temp.append("Naive Bayes Classifier");
    if (m_Instances == null) {
      temp.append(": No model built yet.");
    } else {

      int maxWidth = 0;
      int maxAttWidth = 0;
      boolean containsKernel = false;

      // set up max widths
      // class values
      for (int i = 0; i < m_Instances.numClasses(); i++) {
        if (m_Instances.classAttribute().value(i).length() > maxWidth) {
          maxWidth = m_Instances.classAttribute().value(i).length();
        }
      }
      // attributes
      for (int i = 0; i < m_Instances.numAttributes(); i++) {
        if (i != m_Instances.classIndex()) {
          Attribute a = m_Instances.attribute(i);
          if (a.name().length() > maxAttWidth) {
            maxAttWidth = m_Instances.attribute(i).name().length();
          }
          //RANKING BEGIN
          if (a.isNominal() || a.isRanking()) {
        	  //RANKING END
            // check values
            for (int j = 0; j < a.numValues(); j++) {
              String val = a.value(j) + "  ";
              if (val.length() > maxAttWidth) {
                maxAttWidth = val.length();
              }
            }
          }
        }
      }

      for (int i = 0; i < m_Distributions.length; i++) {
        for (int j = 0; j < m_Instances.numClasses(); j++) {
          if (m_Distributions[i][0] instanceof NormalEstimator) {
            // check mean/precision dev against maxWidth
            NormalEstimator n = (NormalEstimator)m_Distributions[i][j];
            double mean = Math.log(Math.abs(n.getMean())) / Math.log(10.0);
            double precision = Math.log(Math.abs(n.getPrecision())) / Math.log(10.0);
            double width = (mean > precision)
              ? mean
              : precision;
            if (width < 0) {
              width = 1;
            }
            // decimal + # decimal places + 1
            width += 6.0;
            if ((int)width > maxWidth) {
              maxWidth = (int)width;
            }
          } else if (m_Distributions[i][0] instanceof KernelEstimator) {
            containsKernel = true;
            KernelEstimator ke = (KernelEstimator)m_Distributions[i][j];
            int numK = ke.getNumKernels();
            String temps = "K" + numK + ": mean (weight)";
            if (maxAttWidth < temps.length()) {
              maxAttWidth = temps.length();
            }
            // check means + weights against maxWidth
            if (ke.getNumKernels() > 0) {
              double[] means = ke.getMeans();
              double[] weights = ke.getWeights();
              for (int k = 0; k < ke.getNumKernels(); k++) {
                String m = Utils.doubleToString(means[k], maxWidth, 4).trim();
                m += " (" + Utils.doubleToString(weights[k], maxWidth, 1).trim() + ")";
                if (maxWidth < m.length()) {
                  maxWidth = m.length();
                }
              }
            }
          } else if (m_Distributions[i][0] instanceof DiscreteEstimator) {
            DiscreteEstimator d = (DiscreteEstimator)m_Distributions[i][j];
            for (int k = 0; k < d.getNumSymbols(); k++) {
              String size = "" + d.getCount(k);
              if (size.length() > maxWidth) {
                maxWidth = size.length();
              }
            }
            int sum = ("" + d.getSumOfCounts()).length();
            if (sum > maxWidth) {
              maxWidth = sum;
            }
          }
        }
      }

      // Check width of class labels
      for (int i = 0; i < m_Instances.numClasses(); i++) {
        String cSize = m_Instances.classAttribute().value(i);
        if (cSize.length() > maxWidth) {
          maxWidth = cSize.length();
        }
      }

      // Check width of class priors
      for (int i = 0; i < m_Instances.numClasses(); i++) {
        String priorP = 
          Utils.doubleToString(((DiscreteEstimator)m_ClassDistribution).getProbability(i),
                               maxWidth, 2).trim();
        priorP = "(" + priorP + ")";
        if (priorP.length() > maxWidth) {
          maxWidth = priorP.length();
        }
      }
    
      if (maxAttWidth < "Attribute".length()) {
        maxAttWidth = "Attribute".length();
      }

      if (maxAttWidth < "  weight sum".length()) {
        maxAttWidth = "  weight sum".length();
      }

      if (containsKernel) {
        if (maxAttWidth < "  [precision]".length()) {
          maxAttWidth = "  [precision]".length();
        }
      }

      maxAttWidth += 2;
    


      temp.append("\n\n");
      temp.append(pad("Class", " ", 
                      (maxAttWidth + maxWidth + 1) - "Class".length(), 
                      true));

      temp.append("\n");
      temp.append(pad("Attribute", " ", maxAttWidth - "Attribute".length(), false));
      // class labels
      for (int i = 0; i < m_Instances.numClasses(); i++) {
        String classL = m_Instances.classAttribute().value(i);
        temp.append(pad(classL, " ", maxWidth + 1 - classL.length(), true));
      }
      temp.append("\n");
      // class priors
      temp.append(pad("", " ", maxAttWidth, true));
      for (int i = 0; i < m_Instances.numClasses(); i++) {
        String priorP = 
          Utils.doubleToString(((DiscreteEstimator)m_ClassDistribution).getProbability(i),
                               maxWidth, 2).trim();
        priorP = "(" + priorP + ")";
        temp.append(pad(priorP, " ", maxWidth + 1 - priorP.length(), true));
      }
      temp.append("\n");
      temp.append(pad("", "=", maxAttWidth + 
                      (maxWidth * m_Instances.numClasses()) 
                      + m_Instances.numClasses() + 1, true));
      temp.append("\n");

      // loop over the attributes
      int counter = 0;
      for (int i = 0; i < m_Instances.numAttributes(); i++) {
        if (i == m_Instances.classIndex()) {
          continue;
        }
        String attName = m_Instances.attribute(i).name();
        temp.append(attName + "\n");
          
        if (m_Distributions[counter][0] instanceof NormalEstimator) {
          String meanL = "  mean";
          temp.append(pad(meanL, " ", maxAttWidth + 1 - meanL.length(), false));
          for (int j = 0; j < m_Instances.numClasses(); j++) {            
            // means
            NormalEstimator n = (NormalEstimator)m_Distributions[counter][j];
            String mean = 
              Utils.doubleToString(n.getMean(), maxWidth, 4).trim();
            temp.append(pad(mean, " ", maxWidth + 1 - mean.length(), true));
          }
          temp.append("\n");            
          // now do std deviations
          String stdDevL = "  std. dev.";
          temp.append(pad(stdDevL, " ", maxAttWidth + 1 - stdDevL.length(), false));
          for (int j = 0; j < m_Instances.numClasses(); j++) {
            NormalEstimator n = (NormalEstimator)m_Distributions[counter][j];
            String stdDev = 
              Utils.doubleToString(n.getStdDev(), maxWidth, 4).trim();
            temp.append(pad(stdDev, " ", maxWidth + 1 - stdDev.length(), true));
          }
          temp.append("\n");
          // now the weight sums
          String weightL = "  weight sum";
          temp.append(pad(weightL, " ", maxAttWidth + 1 - weightL.length(), false));
          for (int j = 0; j < m_Instances.numClasses(); j++) {
            NormalEstimator n = (NormalEstimator)m_Distributions[counter][j];
            String weight = 
              Utils.doubleToString(n.getSumOfWeights(), maxWidth, 4).trim();
            temp.append(pad(weight, " ", maxWidth + 1 - weight.length(), true));
          }
          temp.append("\n");
          // now the precisions
          String precisionL = "  precision";
          temp.append(pad(precisionL, " ", maxAttWidth + 1 - precisionL.length(), false));
          for (int j = 0; j < m_Instances.numClasses(); j++) {
            NormalEstimator n = (NormalEstimator)m_Distributions[counter][j];
            String precision = 
              Utils.doubleToString(n.getPrecision(), maxWidth, 4).trim();
            temp.append(pad(precision, " ", maxWidth + 1 - precision.length(), true));
          }
          temp.append("\n\n");
            
        } else if (m_Distributions[counter][0] instanceof DiscreteEstimator) {
          Attribute a = m_Instances.attribute(i);
          for (int j = 0; j < a.numValues(); j++) {
            String val = "  " + a.value(j);
            temp.append(pad(val, " ", maxAttWidth + 1 - val.length(), false));
            for (int k = 0; k < m_Instances.numClasses(); k++) {
              DiscreteEstimator d = (DiscreteEstimator)m_Distributions[counter][k];
              String count = "" + d.getCount(j);
              temp.append(pad(count, " ", maxWidth + 1 - count.length(), true));
            }
            temp.append("\n");
          }
          // do the totals
          String total = "  [total]";
          temp.append(pad(total, " ", maxAttWidth + 1 - total.length(), false));
          for (int k = 0; k < m_Instances.numClasses(); k++) {
            DiscreteEstimator d = (DiscreteEstimator)m_Distributions[counter][k];
            String count = "" + d.getSumOfCounts();
            temp.append(pad(count, " ", maxWidth + 1 - count.length(), true));
          }
          temp.append("\n\n");
        } else if (m_Distributions[counter][0] instanceof KernelEstimator) {
          String kL = "  [# kernels]";
          temp.append(pad(kL, " ", maxAttWidth + 1 - kL.length(), false));
          for (int k = 0; k < m_Instances.numClasses(); k++) {
            KernelEstimator ke = (KernelEstimator)m_Distributions[counter][k];
            String nk = "" + ke.getNumKernels();
            temp.append(pad(nk, " ", maxWidth + 1 - nk.length(), true));
          }
          temp.append("\n");
          // do num kernels, std. devs and precisions
          String stdDevL = "  [std. dev]";
          temp.append(pad(stdDevL, " ", maxAttWidth + 1 - stdDevL.length(), false));
          for (int k = 0; k < m_Instances.numClasses(); k++) {
            KernelEstimator ke = (KernelEstimator)m_Distributions[counter][k];
            String stdD = Utils.doubleToString(ke.getStdDev(), maxWidth, 4).trim(); 
            temp.append(pad(stdD, " ", maxWidth + 1 - stdD.length(), true));
          }
          temp.append("\n");
          String precL = "  [precision]";
          temp.append(pad(precL, " ", maxAttWidth + 1 - precL.length(), false));
          for (int k = 0; k < m_Instances.numClasses(); k++) {
            KernelEstimator ke = (KernelEstimator)m_Distributions[counter][k];
            String prec = Utils.doubleToString(ke.getPrecision(), maxWidth, 4).trim(); 
            temp.append(pad(prec, " ", maxWidth + 1 - prec.length(), true));
          }
          temp.append("\n");
          // first determine max number of kernels accross the classes
          int maxK = 0;
          for (int k = 0; k < m_Instances.numClasses(); k++) {
            KernelEstimator ke = (KernelEstimator)m_Distributions[counter][k];
            if (ke.getNumKernels() > maxK) {
              maxK = ke.getNumKernels();
            }
          }
          for (int j = 0; j < maxK; j++) {
            // means first
            String meanL = "  K" + (j+1) + ": mean (weight)";
            temp.append(pad(meanL, " ", maxAttWidth + 1 - meanL.length(), false));
            for (int k = 0; k < m_Instances.numClasses(); k++) {
              KernelEstimator ke = (KernelEstimator)m_Distributions[counter][k];
              double[] means = ke.getMeans();
              double[] weights = ke.getWeights();
              String m = "--";
              if (ke.getNumKernels() == 0) {
                m = "" + 0;
              } else if (j < ke.getNumKernels()) {
                m = Utils.doubleToString(means[j], maxWidth, 4).trim();
                m += " (" + Utils.doubleToString(weights[j], maxWidth, 1).trim() + ")";
              }
              temp.append(pad(m, " ", maxWidth + 1 - m.length(), true));
            }
            temp.append("\n");              
          }
          temp.append("\n");
        }


        counter++;
      }
    }
      
    return temp.toString();
  }

  /**
   * Returns a description of the classifier in the old format.
   *
   * @return a description of the classifier as a string.
   */
  protected String toStringOriginal() {
    
    StringBuffer text = new StringBuffer();

    text.append("Naive Bayes Classifier");
    if (m_Instances == null) {
      text.append(": No model built yet.");
    } else {
      try {
	for (int i = 0; i < m_Distributions[0].length; i++) {
	  text.append("\n\nClass " + m_Instances.classAttribute().value(i) +
                      ": Prior probability = " + Utils.
                      doubleToString(m_ClassDistribution.getProbability(i),
                                     4, 2) + "\n\n");
	  Enumeration enumAtts = m_Instances.enumerateAttributes();
	  int attIndex = 0;
	  while (enumAtts.hasMoreElements()) {
	    Attribute attribute = (Attribute) enumAtts.nextElement();
	    if (attribute.weight() > 0) {
	      text.append(attribute.name() + ":  " 
                          + m_Distributions[attIndex][i]);
	    }
	    attIndex++;
	  }
	}
      } catch (Exception ex) {
	text.append(ex.getMessage());
      }
    }

    return text.toString();
  } 

  private String pad(String source, String padChar, 
                     int length, boolean leftPad) {
    StringBuffer temp = new StringBuffer();

    if (leftPad) {
      for (int i = 0; i< length; i++) {
        temp.append(padChar);
      }
      temp.append(source);
    } else {
      temp.append(source);
      for (int i = 0; i< length; i++) {
        temp.append(padChar);
      }
    }
    return temp.toString();
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String useKernelEstimatorTipText() {
    return "Use a kernel estimator for numeric attributes rather than a "
      +"normal distribution.";
  }
  /**
   * Gets if kernel estimator is being used.
   *
   * @return Value of m_UseKernelEstimatory.
   */
  public boolean getUseKernelEstimator() {

    return m_UseKernelEstimator;
  }

  /**
   * Sets if kernel estimator is to be used.
   *
   * @param v  Value to assign to m_UseKernelEstimatory.
   */
  public void setUseKernelEstimator(boolean v) {

    m_UseKernelEstimator = v;
    if (v) {
      setUseSupervisedDiscretization(false);
    }
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String useSupervisedDiscretizationTipText() {
    return "Use supervised discretization to convert numeric attributes to nominal "
      +"ones.";
  }

  /**
   * Get whether supervised discretization is to be used.
   *
   * @return true if supervised discretization is to be used.
   */
  public boolean getUseSupervisedDiscretization() {

    return m_UseDiscretization;
  }

  /**
   * Set whether supervised discretization is to be used.
   *
   * @param newblah true if supervised discretization is to be used.
   */
  public void setUseSupervisedDiscretization(boolean newblah) {

    m_UseDiscretization = newblah;
    if (newblah) {
      setUseKernelEstimator(false);
    }
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String displayModelInOldFormatTipText() {
    return "Use old format for model output. The old format is "
      + "better when there are many class values. The new format "
      + "is better when there are fewer classes and many attributes.";
  }

  /**
   * Set whether to display model output in the old, original
   * format.
   *
   * @param d true if model ouput is to be shown in the old format
   */
  public void setDisplayModelInOldFormat(boolean d) {
    m_displayModelInOldFormat = d;
  }

  /**
   * Get whether to display model output in the old, original
   * format.
   *
   * @return true if model ouput is to be shown in the old format
   */
  public boolean getDisplayModelInOldFormat() {
    return m_displayModelInOldFormat;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5928 $");
  }

  /**
   * Main method for testing this class.
   *
   * @param argv the options
   */
  public static void main(String [] argv) {
    runClassifier(new NaiveBayes(), argv);
  }
}

