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
 *    MultiScheme.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.meta;

import weka.classifiers.Classifier;
import weka.classifiers.AbstractClassifier;
import weka.classifiers.Evaluation;
import weka.classifiers.RandomizableMultipleClassifiersCombiner;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.Utils;

import java.util.Enumeration;
import java.util.Random;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * Class for selecting a classifier from among several using cross validation on the training data or the performance on the training data. Performance is measured based on percent correct (classification) or mean-squared error (regression).
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -X &lt;number of folds&gt;
 *  Use cross validation for model selection using the
 *  given number of folds. (default 0, is to
 *  use training error)</pre>
 * 
 * <pre> -S &lt;num&gt;
 *  Random number seed.
 *  (default 1)</pre>
 * 
 * <pre> -B &lt;classifier specification&gt;
 *  Full class name of classifier to include, followed
 *  by scheme options. May be specified multiple times.
 *  (default: "weka.classifiers.rules.ZeroR")</pre>
 * 
 * <pre> -D
 *  If set, classifier is run in debug mode and
 *  may output additional info to the console</pre>
 * 
 <!-- options-end -->
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 5928 $
 */
public class MultiScheme 
  extends RandomizableMultipleClassifiersCombiner {

  /** for serialization */
  static final long serialVersionUID = 5710744346128957520L;
  
  /** The classifier that had the best performance on training data. */
  protected Classifier m_Classifier;
 
  /** The index into the vector for the selected scheme */
  protected int m_ClassifierIndex;

  /**
   * Number of folds to use for cross validation (0 means use training
   * error for selection)
   */
  protected int m_NumXValFolds;
    
  /**
   * Returns a string describing classifier
   * @return a description suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {

    return  "Class for selecting a classifier from among several using cross "
      + "validation on the training data or the performance on the "
      + "training data. Performance is measured based on percent correct "
      + "(classification) or mean-squared error (regression).";
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {

    Vector newVector = new Vector(1);
    newVector.addElement(new Option(
	      "\tUse cross validation for model selection using the\n"
	      + "\tgiven number of folds. (default 0, is to\n"
	      + "\tuse training error)",
	      "X", 1, "-X <number of folds>"));

    Enumeration enu = super.listOptions();
    while (enu.hasMoreElements()) {
      newVector.addElement(enu.nextElement());
    }
    return newVector.elements();
  }

  /**
   * Parses a given list of options. <p/>
   *
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -X &lt;number of folds&gt;
   *  Use cross validation for model selection using the
   *  given number of folds. (default 0, is to
   *  use training error)</pre>
   * 
   * <pre> -S &lt;num&gt;
   *  Random number seed.
   *  (default 1)</pre>
   * 
   * <pre> -B &lt;classifier specification&gt;
   *  Full class name of classifier to include, followed
   *  by scheme options. May be specified multiple times.
   *  (default: "weka.classifiers.rules.ZeroR")</pre>
   * 
   * <pre> -D
   *  If set, classifier is run in debug mode and
   *  may output additional info to the console</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    
    String numFoldsString = Utils.getOption('X', options);
    if (numFoldsString.length() != 0) {
      setNumFolds(Integer.parseInt(numFoldsString));
    } else {
      setNumFolds(0);
    }
    super.setOptions(options);
  }

  /**
   * Gets the current settings of the Classifier.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String [] getOptions() {


    String [] superOptions = super.getOptions();
    String [] options = new String [superOptions.length + 2];

    int current = 0;
    options[current++] = "-X"; options[current++] = "" + getNumFolds();

    System.arraycopy(superOptions, 0, options, current, 
		     superOptions.length);

    return options;
  }
  
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String classifiersTipText() {
    return "The classifiers to be chosen from.";
  }

  /**
   * Sets the list of possible classifers to choose from.
   *
   * @param classifiers an array of classifiers with all options set.
   */
  public void setClassifiers(Classifier [] classifiers) {

    m_Classifiers = classifiers;
  }

  /**
   * Gets the list of possible classifers to choose from.
   *
   * @return the array of Classifiers
   */
  public Classifier [] getClassifiers() {

    return m_Classifiers;
  }
  
  /**
   * Gets a single classifier from the set of available classifiers.
   *
   * @param index the index of the classifier wanted
   * @return the Classifier
   */
  public Classifier getClassifier(int index) {

    return m_Classifiers[index];
  }
  
  /**
   * Gets the classifier specification string, which contains the class name of
   * the classifier and any options to the classifier
   *
   * @param index the index of the classifier string to retrieve, starting from
   * 0.
   * @return the classifier string, or the empty string if no classifier
   * has been assigned (or the index given is out of range).
   */
  protected String getClassifierSpec(int index) {
    
    if (m_Classifiers.length < index) {
      return "";
    }
    Classifier c = getClassifier(index);
    if (c instanceof OptionHandler) {
      return c.getClass().getName() + " "
	+ Utils.joinOptions(((OptionHandler)c).getOptions());
    }
    return c.getClass().getName();
  }
  
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String seedTipText() {
    return "The seed used for randomizing the data " +
      "for cross-validation.";
  }

  /**
   * Sets the seed for random number generation.
   *
   * @param seed the random number seed
   */
  public void setSeed(int seed) {
    
    m_Seed = seed;;
  }

  /**
   * Gets the random number seed.
   * 
   * @return the random number seed
   */
  public int getSeed() {

    return m_Seed;
  }
  
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String numFoldsTipText() {
    return "The number of folds used for cross-validation (if 0, " +
      "performance on training data will be used).";
  }

  /** 
   * Gets the number of folds for cross-validation. A number less
   * than 2 specifies using training error rather than cross-validation.
   *
   * @return the number of folds for cross-validation
   */
  public int getNumFolds() {

    return m_NumXValFolds;
  }

  /**
   * Sets the number of folds for cross-validation. A number less
   * than 2 specifies using training error rather than cross-validation.
   *
   * @param numFolds the number of folds for cross-validation
   */
  public void setNumFolds(int numFolds) {
    
    m_NumXValFolds = numFolds;
  }
  
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String debugTipText() {
    return "Whether debug information is output to console.";
  }

  /**
   * Set debugging mode
   *
   * @param debug true if debug output should be printed
   */
  public void setDebug(boolean debug) {

    m_Debug = debug;
  }

  /**
   * Get whether debugging is turned on
   *
   * @return true if debugging output is on
   */
  public boolean getDebug() {

    return m_Debug;
  }
  
  /**
   * Get the index of the classifier that was determined as best during 
   * cross-validation.
   * 
   * @return the index in the classifier array
   */
  public int getBestClassifierIndex() {
    return m_ClassifierIndex;
  }

  /**
   * Buildclassifier selects a classifier from the set of classifiers
   * by minimising error on the training data.
   *
   * @param data the training data to be used for generating the
   * boosted classifier.
   * @throws Exception if the classifier could not be built successfully
   */
  public void buildClassifier(Instances data) throws Exception {

    if (m_Classifiers.length == 0) {
      throw new Exception("No base classifiers have been set!");
    }

    // can classifier handle the data?
    getCapabilities().testWithFail(data);

    // remove instances with missing class
    Instances newData = new Instances(data);
    newData.deleteWithMissingClass();
    
    Random random = new Random(m_Seed);
    newData.randomize(random);
    if ((newData.classAttribute().isNominal() || newData.classAttribute().isRanking()) && (m_NumXValFolds > 1)) {
      newData.stratify(m_NumXValFolds);
    }
    Instances train = newData;               // train on all data by default
    Instances test = newData;               // test on training data by default
    Classifier bestClassifier = null;
    int bestIndex = -1;
    double bestPerformance = Double.NaN;
    int numClassifiers = m_Classifiers.length;
    for (int i = 0; i < numClassifiers; i++) {
      Classifier currentClassifier = getClassifier(i);
      Evaluation evaluation;
      if (m_NumXValFolds > 1) {
	evaluation = new Evaluation(newData);
	for (int j = 0; j < m_NumXValFolds; j++) {

          // We want to randomize the data the same way for every 
          // learning scheme.
	  train = newData.trainCV(m_NumXValFolds, j, new Random (1));
	  test = newData.testCV(m_NumXValFolds, j);
	  currentClassifier.buildClassifier(train);
	  evaluation.setPriors(train);
	  evaluation.evaluateModel(currentClassifier, test);
	}
      } else {
	currentClassifier.buildClassifier(train);
	evaluation = new Evaluation(train);
	evaluation.evaluateModel(currentClassifier, test);
      }

      double error = evaluation.errorRate();
      if (m_Debug) {
	System.err.println("Error rate: " + Utils.doubleToString(error, 6, 4)
			   + " for classifier "
			   + currentClassifier.getClass().getName());
      }

      if ((i == 0) || (error < bestPerformance)) {
	bestClassifier = currentClassifier;
	bestPerformance = error;
	bestIndex = i;
      }
    }
    m_ClassifierIndex = bestIndex;
    if (m_NumXValFolds > 1) {
      bestClassifier.buildClassifier(newData);
    }
    m_Classifier = bestClassifier;
  }

  /**
   * Returns class probabilities.
   *
   * @param instance the instance to be classified
   * @return the distribution for the instance
   * @throws Exception if instance could not be classified
   * successfully
   */
  public double[] distributionForInstance(Instance instance) throws Exception {

    return m_Classifier.distributionForInstance(instance);
  }

  /**
   * Output a representation of this classifier
   * @return a string representation of the classifier
   */
  public String toString() {

    if (m_Classifier == null) {
      return "MultiScheme: No model built yet.";
    }

    String result = "MultiScheme selection using";
    if (m_NumXValFolds > 1) {
      result += " cross validation error";
    } else {
      result += " error on training data";
    }
    result += " from the following:\n";
    for (int i = 0; i < m_Classifiers.length; i++) {
      result += '\t' + getClassifierSpec(i) + '\n';
    }

    result += "Selected scheme: "
      + getClassifierSpec(m_ClassifierIndex)
      + "\n\n"
      + m_Classifier.toString();
    return result;
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
   * @param argv should contain the following arguments:
   * -t training file [-T test file] [-c class index]
   */
  public static void main(String [] argv) {
    runClassifier(new MultiScheme(), argv);
  }
}
