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
 *    AdditiveRegression.java
 *    Copyright (C) 2000 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.meta;

import weka.classifiers.Classifier;
import weka.classifiers.AbstractClassifier;
import weka.classifiers.IteratedSingleClassifierEnhancer;
import weka.classifiers.rules.ZeroR;
import weka.core.AdditionalMeasureProducer;
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

import java.util.Enumeration;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * Meta classifier that enhances the performance of a regression base classifier. Each iteration fits a model to the residuals left by the classifier on the previous iteration. Prediction is accomplished by adding the predictions of each classifier. Reducing the shrinkage (learning rate) parameter helps prevent overfitting and has a smoothing effect but increases the learning time.<br/>
 * <br/>
 * For more information see:<br/>
 * <br/>
 * J.H. Friedman (1999). Stochastic Gradient Boosting.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;techreport{Friedman1999,
 *    author = {J.H. Friedman},
 *    institution = {Stanford University},
 *    title = {Stochastic Gradient Boosting},
 *    year = {1999},
 *    PS = {http://www-stat.stanford.edu/\~jhf/ftp/stobst.ps}
 * }
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -S
 *  Specify shrinkage rate. (default = 1.0, ie. no shrinkage)
 * </pre>
 * 
 * <pre> -I &lt;num&gt;
 *  Number of iterations.
 *  (default 10)</pre>
 * 
 * <pre> -D
 *  If set, classifier is run in debug mode and
 *  may output additional info to the console</pre>
 * 
 * <pre> -W
 *  Full name of base classifier.
 *  (default: weka.classifiers.trees.DecisionStump)</pre>
 * 
 * <pre> 
 * Options specific to classifier weka.classifiers.trees.DecisionStump:
 * </pre>
 * 
 * <pre> -D
 *  If set, classifier is run in debug mode and
 *  may output additional info to the console</pre>
 * 
 <!-- options-end -->
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @version $Revision: 5928 $
 */
public class AdditiveRegression 
  extends IteratedSingleClassifierEnhancer 
  implements OptionHandler,
	     AdditionalMeasureProducer,
	     WeightedInstancesHandler,
	     TechnicalInformationHandler {

  /** for serialization */
  static final long serialVersionUID = -2368937577670527151L;
  
  /**
   * Shrinkage (Learning rate). Default = no shrinkage.
   */
  protected double m_shrinkage = 1.0;

  /** The number of successfully generated base classifiers. */
  protected int m_NumIterationsPerformed;

  /** The model for the mean */
  protected ZeroR m_zeroR;

  /** whether we have suitable data or nor (if not, ZeroR model is used) */
  protected boolean m_SuitableData = true;
  
  /**
   * Returns a string describing this attribute evaluator
   * @return a description of the evaluator suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return " Meta classifier that enhances the performance of a regression "
      +"base classifier. Each iteration fits a model to the residuals left "
      +"by the classifier on the previous iteration. Prediction is "
      +"accomplished by adding the predictions of each classifier. "
      +"Reducing the shrinkage (learning rate) parameter helps prevent "
      +"overfitting and has a smoothing effect but increases the learning "
      +"time.\n\n"
      +"For more information see:\n\n"
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
    
    result = new TechnicalInformation(Type.TECHREPORT);
    result.setValue(Field.AUTHOR, "J.H. Friedman");
    result.setValue(Field.YEAR, "1999");
    result.setValue(Field.TITLE, "Stochastic Gradient Boosting");
    result.setValue(Field.INSTITUTION, "Stanford University");
    result.setValue(Field.PS, "http://www-stat.stanford.edu/~jhf/ftp/stobst.ps");
    
    return result;
  }

  /**
   * Default constructor specifying DecisionStump as the classifier
   */
  public AdditiveRegression() {

    this(new weka.classifiers.trees.DecisionStump());
  }

  /**
   * Constructor which takes base classifier as argument.
   *
   * @param classifier the base classifier to use
   */
  public AdditiveRegression(Classifier classifier) {

    m_Classifier = classifier;
  }

  /**
   * String describing default classifier.
   * 
   * @return the default classifier classname
   */
  protected String defaultClassifierString() {
    
    return "weka.classifiers.trees.DecisionStump";
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {

    Vector newVector = new Vector(4);

    newVector.addElement(new Option(
	      "\tSpecify shrinkage rate. "
	      +"(default = 1.0, ie. no shrinkage)\n", 
	      "S", 1, "-S"));

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
   * <pre> -S
   *  Specify shrinkage rate. (default = 1.0, ie. no shrinkage)
   * </pre>
   * 
   * <pre> -I &lt;num&gt;
   *  Number of iterations.
   *  (default 10)</pre>
   * 
   * <pre> -D
   *  If set, classifier is run in debug mode and
   *  may output additional info to the console</pre>
   * 
   * <pre> -W
   *  Full name of base classifier.
   *  (default: weka.classifiers.trees.DecisionStump)</pre>
   * 
   * <pre> 
   * Options specific to classifier weka.classifiers.trees.DecisionStump:
   * </pre>
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

    String optionString = Utils.getOption('S', options);
    if (optionString.length() != 0) {
      Double temp = Double.valueOf(optionString);
      setShrinkage(temp.doubleValue());
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

    options[current++] = "-S"; options[current++] = "" + getShrinkage();

    System.arraycopy(superOptions, 0, options, current, 
		     superOptions.length);

    current += superOptions.length;
    while (current < options.length) {
      options[current++] = "";
    }
    return options;
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String shrinkageTipText() {
    return "Shrinkage rate. Smaller values help prevent overfitting and "
      + "have a smoothing effect (but increase learning time). "
      +"Default = 1.0, ie. no shrinkage."; 
  }

  /**
   * Set the shrinkage parameter
   *
   * @param l the shrinkage rate.
   */
  public void setShrinkage(double l) {
    m_shrinkage = l;
  }

  /**
   * Get the shrinkage rate.
   *
   * @return the value of the learning rate
   */
  public double getShrinkage() {
    return m_shrinkage;
  }

  /**
   * Returns default capabilities of the classifier.
   *
   * @return      the capabilities of this classifier
   */
  public Capabilities getCapabilities() {
    Capabilities result = super.getCapabilities();

    // class
    result.disableAllClasses();
    result.disableAllClassDependencies();
    result.enable(Capability.NUMERIC_CLASS);
    result.enable(Capability.DATE_CLASS);
    
    return result;
  }

  /**
   * Build the classifier on the supplied data
   *
   * @param data the training data
   * @throws Exception if the classifier could not be built successfully
   */
  public void buildClassifier(Instances data) throws Exception {

    super.buildClassifier(data);

    // can classifier handle the data?
    getCapabilities().testWithFail(data);

    // remove instances with missing class
    Instances newData = new Instances(data);
    newData.deleteWithMissingClass();

    double sum = 0;
    double temp_sum = 0;
    // Add the model for the mean first
    m_zeroR = new ZeroR();
    m_zeroR.buildClassifier(newData);
    
    // only class? -> use only ZeroR model
    if (newData.numAttributes() == 1) {
      System.err.println(
	  "Cannot build model (only class attribute present in data!), "
	  + "using ZeroR model instead!");
      m_SuitableData = false;
      return;
    }
    else {
      m_SuitableData = true;
    }
    
    newData = residualReplace(newData, m_zeroR, false);
    for (int i = 0; i < newData.numInstances(); i++) {
      sum += newData.instance(i).weight() *
	newData.instance(i).classValue() * newData.instance(i).classValue();
    }
    if (m_Debug) {
      System.err.println("Sum of squared residuals "
			 +"(predicting the mean) : " + sum);
    }

    m_NumIterationsPerformed = 0;
    do {
      temp_sum = sum;

      // Build the classifier
      m_Classifiers[m_NumIterationsPerformed].buildClassifier(newData);

      newData = residualReplace(newData, m_Classifiers[m_NumIterationsPerformed], true);
      sum = 0;
      for (int i = 0; i < newData.numInstances(); i++) {
	sum += newData.instance(i).weight() *
	  newData.instance(i).classValue() * newData.instance(i).classValue();
      }
      if (m_Debug) {
	System.err.println("Sum of squared residuals : "+sum);
      }
      m_NumIterationsPerformed++;
    } while (((temp_sum - sum) > Utils.SMALL) && 
	     (m_NumIterationsPerformed < m_Classifiers.length));
  }

  /**
   * Classify an instance.
   *
   * @param inst the instance to predict
   * @return a prediction for the instance
   * @throws Exception if an error occurs
   */
  public double classifyInstance(Instance inst) throws Exception {

    double prediction = m_zeroR.classifyInstance(inst);

    // default model?
    if (!m_SuitableData) {
      return prediction;
    }
    
    for (int i = 0; i < m_NumIterationsPerformed; i++) {
      double toAdd = m_Classifiers[i].classifyInstance(inst);
      toAdd *= getShrinkage();
      prediction += toAdd;
    }

    return prediction;
  }

  /**
   * Replace the class values of the instances from the current iteration
   * with residuals ater predicting with the supplied classifier.
   *
   * @param data the instances to predict
   * @param c the classifier to use
   * @param useShrinkage whether shrinkage is to be applied to the model's output
   * @return a new set of instances with class values replaced by residuals
   * @throws Exception if something goes wrong
   */
  private Instances residualReplace(Instances data, Classifier c, 
				    boolean useShrinkage) throws Exception {
    double pred,residual;
    Instances newInst = new Instances(data);

    for (int i = 0; i < newInst.numInstances(); i++) {
      pred = c.classifyInstance(newInst.instance(i));
      if (useShrinkage) {
	pred *= getShrinkage();
      }
      residual = newInst.instance(i).classValue() - pred;
      newInst.instance(i).setClassValue(residual);
    }
    //    System.err.print(newInst);
    return newInst;
  }

  /**
   * Returns an enumeration of the additional measure names
   * @return an enumeration of the measure names
   */
  public Enumeration enumerateMeasures() {
    Vector newVector = new Vector(1);
    newVector.addElement("measureNumIterations");
    return newVector.elements();
  }

  /**
   * Returns the value of the named measure
   * @param additionalMeasureName the name of the measure to query for its value
   * @return the value of the named measure
   * @throws IllegalArgumentException if the named measure is not supported
   */
  public double getMeasure(String additionalMeasureName) {
    if (additionalMeasureName.compareToIgnoreCase("measureNumIterations") == 0) {
      return measureNumIterations();
    } else {
      throw new IllegalArgumentException(additionalMeasureName 
			  + " not supported (AdditiveRegression)");
    }
  }

  /**
   * return the number of iterations (base classifiers) completed
   * @return the number of iterations (same as number of base classifier
   * models)
   */
  public double measureNumIterations() {
    return m_NumIterationsPerformed;
  }

  /**
   * Returns textual description of the classifier.
   *
   * @return a description of the classifier as a string
   */
  public String toString() {
    StringBuffer text = new StringBuffer();

    // only ZeroR model?
    if (!m_SuitableData) {
      StringBuffer buf = new StringBuffer();
      buf.append(this.getClass().getName().replaceAll(".*\\.", "") + "\n");
      buf.append(this.getClass().getName().replaceAll(".*\\.", "").replaceAll(".", "=") + "\n\n");
      buf.append("Warning: No model could be built, hence ZeroR model is used:\n\n");
      buf.append(m_zeroR.toString());
      return buf.toString();
    }
    
    if (m_NumIterations == 0) {
      return "Classifier hasn't been built yet!";
    }

    text.append("Additive Regression\n\n");

    text.append("ZeroR model\n\n" + m_zeroR + "\n\n");

    text.append("Base classifier " 
		+ getClassifier().getClass().getName()
		+ "\n\n");
    text.append("" + m_NumIterationsPerformed + " models generated.\n");

    for (int i = 0; i < m_NumIterationsPerformed; i++) {
      text.append("\nModel number " + i + "\n\n" +
		  m_Classifiers[i] + "\n");
    }

    return text.toString();
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
    runClassifier(new AdditiveRegression(), argv);
  }
}
