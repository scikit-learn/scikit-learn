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
 *    LinearRegression.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.functions;

import weka.classifiers.Classifier;
import weka.classifiers.AbstractClassifier;
import weka.core.Capabilities;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Matrix;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.SelectedTag;
import weka.core.Tag;
import weka.core.Utils;
import weka.core.WeightedInstancesHandler;
import weka.core.Capabilities.Capability;
import weka.filters.Filter;
import weka.filters.supervised.attribute.NominalToBinary;
import weka.filters.unsupervised.attribute.ReplaceMissingValues;

import java.util.Enumeration;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * Class for using linear regression for prediction. Uses the Akaike criterion for model selection, and is able to deal with weighted instances.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -D
 *  Produce debugging output.
 *  (default no debugging output)</pre>
 * 
 * <pre> -S &lt;number of selection method&gt;
 *  Set the attribute selection method to use. 1 = None, 2 = Greedy.
 *  (default 0 = M5' method)</pre>
 * 
 * <pre> -C
 *  Do not try to eliminate colinear attributes.
 * </pre>
 * 
 * <pre> -R &lt;double&gt;
 *  Set ridge parameter (default 1.0e-8).
 * </pre>
 * 
 <!-- options-end -->
 *
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 5928 $
 */
public class LinearRegression extends AbstractClassifier implements OptionHandler,
  WeightedInstancesHandler {
  
  /** for serialization */
  static final long serialVersionUID = -3364580862046573747L;

  /** Array for storing coefficients of linear regression. */
  private double[] m_Coefficients;

  /** Which attributes are relevant? */
  private boolean[] m_SelectedAttributes;

  /** Variable for storing transformed training data. */
  private Instances m_TransformedData;

  /** The filter for removing missing values. */
  private ReplaceMissingValues m_MissingFilter;

  /** The filter storing the transformation from nominal to 
      binary attributes. */
  private NominalToBinary m_TransformFilter;

  /** The standard deviations of the class attribute */
  private double m_ClassStdDev;

  /** The mean of the class attribute */
  private double m_ClassMean;

  /** The index of the class attribute */
  private int m_ClassIndex;

  /** The attributes means */
  private double[] m_Means;

  /** The attribute standard deviations */
  private double[] m_StdDevs;

  /** True if debug output will be printed */
  private boolean b_Debug;

  /** The current attribute selection method */
  private int m_AttributeSelection;

  /** Attribute selection method: M5 method */
  public static final int SELECTION_M5 = 0;
  /** Attribute selection method: No attribute selection */
  public static final int SELECTION_NONE = 1;
  /** Attribute selection method: Greedy method */
  public static final int SELECTION_GREEDY = 2;
  /** Attribute selection methods */
  public static final Tag [] TAGS_SELECTION = {
    new Tag(SELECTION_NONE, "No attribute selection"),
    new Tag(SELECTION_M5, "M5 method"),
    new Tag(SELECTION_GREEDY, "Greedy method")
  };

  /** Try to eliminate correlated attributes? */
  private boolean m_EliminateColinearAttributes = true;

  /** Turn off all checks and conversions? */
  private boolean m_checksTurnedOff = false;

  /** The ridge parameter */
  private double m_Ridge = 1.0e-8;

  /**
   * Turns off checks for missing values, etc. Use with caution.
   * Also turns off scaling.
   */
  public void turnChecksOff() {

    m_checksTurnedOff = true;
  }

  /**
   * Turns on checks for missing values, etc. Also turns
   * on scaling.
   */
  public void turnChecksOn() {

    m_checksTurnedOff = false;
  }

  /**
   * Returns a string describing this classifier
   * @return a description of the classifier suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return "Class for using linear regression for prediction. Uses the Akaike "
      +"criterion for model selection, and is able to deal with weighted "
      +"instances.";
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
    result.enable(Capability.DATE_ATTRIBUTES);
    result.enable(Capability.MISSING_VALUES);

    // class
    result.enable(Capability.NUMERIC_CLASS);
    result.enable(Capability.DATE_CLASS);
    result.enable(Capability.MISSING_CLASS_VALUES);
    
    return result;
  }

  /**
   * Builds a regression model for the given data.
   *
   * @param data the training data to be used for generating the
   * linear regression function
   * @throws Exception if the classifier could not be built successfully
   */
  public void buildClassifier(Instances data) throws Exception {
  
    if (!m_checksTurnedOff) {
      // can classifier handle the data?
      getCapabilities().testWithFail(data);

      // remove instances with missing class
      data = new Instances(data);
      data.deleteWithMissingClass();
    }

    // Preprocess instances
    if (!m_checksTurnedOff) {
      m_TransformFilter = new NominalToBinary();
      m_TransformFilter.setInputFormat(data);
      data = Filter.useFilter(data, m_TransformFilter);
      m_MissingFilter = new ReplaceMissingValues();
      m_MissingFilter.setInputFormat(data);
      data = Filter.useFilter(data, m_MissingFilter);
      data.deleteWithMissingClass();
    } else {
      m_TransformFilter = null;
      m_MissingFilter = null;
    }

    m_ClassIndex = data.classIndex();
    m_TransformedData = data;

    // Turn all attributes on for a start
    m_SelectedAttributes = new boolean[data.numAttributes()];
    for (int i = 0; i < data.numAttributes(); i++) {
      if (i != m_ClassIndex) {
	m_SelectedAttributes[i] = true;
      }
    }
    m_Coefficients = null;

    // Compute means and standard deviations
    m_Means = new double[data.numAttributes()];
    m_StdDevs = new double[data.numAttributes()];
    for (int j = 0; j < data.numAttributes(); j++) {
      if (j != data.classIndex()) {
	m_Means[j] = data.meanOrMode(j);
	m_StdDevs[j] = Math.sqrt(data.variance(j));
	if (m_StdDevs[j] == 0) {
	  m_SelectedAttributes[j] = false;
	} 
      }
    }

    m_ClassStdDev = Math.sqrt(data.variance(m_TransformedData.classIndex()));
    m_ClassMean = data.meanOrMode(m_TransformedData.classIndex());

    // Perform the regression
    findBestModel();

    // Save memory
    m_TransformedData = new Instances(data, 0);
  }

  /**
   * Classifies the given instance using the linear regression function.
   *
   * @param instance the test instance
   * @return the classification
   * @throws Exception if classification can't be done successfully
   */
  public double classifyInstance(Instance instance) throws Exception {

    // Transform the input instance
    Instance transformedInstance = instance;
    if (!m_checksTurnedOff) {
      m_TransformFilter.input(transformedInstance);
      m_TransformFilter.batchFinished();
      transformedInstance = m_TransformFilter.output();
      m_MissingFilter.input(transformedInstance);
      m_MissingFilter.batchFinished();
      transformedInstance = m_MissingFilter.output();
    }

    // Calculate the dependent variable from the regression model
    return regressionPrediction(transformedInstance,
				m_SelectedAttributes,
				m_Coefficients);
  }

  /**
   * Outputs the linear regression model as a string.
   * 
   * @return the model as string
   */
  public String toString() {

    if (m_TransformedData == null) {
      return "Linear Regression: No model built yet.";
    }
    try {
      StringBuffer text = new StringBuffer();
      int column = 0;
      boolean first = true;
      
      text.append("\nLinear Regression Model\n\n");
      
      text.append(m_TransformedData.classAttribute().name()+" =\n\n");
      for (int i = 0; i < m_TransformedData.numAttributes(); i++) {
	if ((i != m_ClassIndex) 
	    && (m_SelectedAttributes[i])) {
	  if (!first) 
	    text.append(" +\n");
	  else
	    first = false;
	  text.append(Utils.doubleToString(m_Coefficients[column], 12, 4)
		      + " * ");
	  text.append(m_TransformedData.attribute(i).name());
	  column++;
	}
      }
      text.append(" +\n" + 
		  Utils.doubleToString(m_Coefficients[column], 12, 4));
      return text.toString();
    } catch (Exception e) {
      return "Can't print Linear Regression!";
    }
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    
    Vector newVector = new Vector(4);
    newVector.addElement(new Option("\tProduce debugging output.\n"
				    + "\t(default no debugging output)",
				    "D", 0, "-D"));
    newVector.addElement(new Option("\tSet the attribute selection method"
				    + " to use. 1 = None, 2 = Greedy.\n"
				    + "\t(default 0 = M5' method)",
				    "S", 1, "-S <number of selection method>"));
    newVector.addElement(new Option("\tDo not try to eliminate colinear"
				    + " attributes.\n",
				    "C", 0, "-C"));
    newVector.addElement(new Option("\tSet ridge parameter (default 1.0e-8).\n",
				    "R", 1, "-R <double>"));
    return newVector.elements();
  }

  /**
   * Parses a given list of options. <p/>
   *
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -D
   *  Produce debugging output.
   *  (default no debugging output)</pre>
   * 
   * <pre> -S &lt;number of selection method&gt;
   *  Set the attribute selection method to use. 1 = None, 2 = Greedy.
   *  (default 0 = M5' method)</pre>
   * 
   * <pre> -C
   *  Do not try to eliminate colinear attributes.
   * </pre>
   * 
   * <pre> -R &lt;double&gt;
   *  Set ridge parameter (default 1.0e-8).
   * </pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {

    String selectionString = Utils.getOption('S', options);
    if (selectionString.length() != 0) {
      setAttributeSelectionMethod(new SelectedTag(Integer
						  .parseInt(selectionString),
						  TAGS_SELECTION));
    } else {
      setAttributeSelectionMethod(new SelectedTag(SELECTION_M5,
						  TAGS_SELECTION));
    }
    String ridgeString = Utils.getOption('R', options);
    if (ridgeString.length() != 0) {
      setRidge(new Double(ridgeString).doubleValue());
    } else {
      setRidge(1.0e-8);
    }
    setDebug(Utils.getFlag('D', options));
    setEliminateColinearAttributes(!Utils.getFlag('C', options));
  }

  /**
   * Returns the coefficients for this linear model.
   * 
   * @return the coefficients for this linear model
   */
  public double[] coefficients() {

    double[] coefficients = new double[m_SelectedAttributes.length + 1];
    int counter = 0;
    for (int i = 0; i < m_SelectedAttributes.length; i++) {
      if ((m_SelectedAttributes[i]) && ((i != m_ClassIndex))) {
	coefficients[i] = m_Coefficients[counter++];
      }
    }
    coefficients[m_SelectedAttributes.length] = m_Coefficients[counter];
    return coefficients;
  }

  /**
   * Gets the current settings of the classifier.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String [] getOptions() {

    String [] options = new String [6];
    int current = 0;

    options[current++] = "-S";
    options[current++] = "" + getAttributeSelectionMethod()
      .getSelectedTag().getID();
    if (getDebug()) {
      options[current++] = "-D";
    }
    if (!getEliminateColinearAttributes()) {
      options[current++] = "-C";
    }
    options[current++] = "-R";
    options[current++] = "" + getRidge();

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
  public String ridgeTipText() {
    return "The value of the Ridge parameter.";
  }

  /**
   * Get the value of Ridge.
   *
   * @return Value of Ridge.
   */
  public double getRidge() {
    
    return m_Ridge;
  }
  
  /**
   * Set the value of Ridge.
   *
   * @param newRidge Value to assign to Ridge.
   */
  public void setRidge(double newRidge) {
    
    m_Ridge = newRidge;
  }
  
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String eliminateColinearAttributesTipText() {
    return "Eliminate colinear attributes.";
  }

  /**
   * Get the value of EliminateColinearAttributes.
   *
   * @return Value of EliminateColinearAttributes.
   */
  public boolean getEliminateColinearAttributes() {
    
    return m_EliminateColinearAttributes;
  }
  
  /**
   * Set the value of EliminateColinearAttributes.
   *
   * @param newEliminateColinearAttributes Value to assign to EliminateColinearAttributes.
   */
  public void setEliminateColinearAttributes(boolean newEliminateColinearAttributes) {
    
    m_EliminateColinearAttributes = newEliminateColinearAttributes;
  }
  
  /**
   * Get the number of coefficients used in the model
   *
   * @return the number of coefficients
   */
  public int numParameters()
  {
    return m_Coefficients.length-1;
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String attributeSelectionMethodTipText() {
    return "Set the method used to select attributes for use in the linear "
      +"regression. Available methods are: no attribute selection, attribute "
      +"selection using M5's method (step through the attributes removing the one "
      +"with the smallest standardised coefficient until no improvement is observed "
      +"in the estimate of the error given by the Akaike "
      +"information criterion), and a greedy selection using the Akaike information "
      +"metric.";
  }

  /**
   * Sets the method used to select attributes for use in the
   * linear regression. 
   *
   * @param method the attribute selection method to use.
   */
  public void setAttributeSelectionMethod(SelectedTag method) {
    
    if (method.getTags() == TAGS_SELECTION) {
      m_AttributeSelection = method.getSelectedTag().getID();
    }
  }

  /**
   * Gets the method used to select attributes for use in the
   * linear regression. 
   *
   * @return the method to use.
   */
  public SelectedTag getAttributeSelectionMethod() {
    
    return new SelectedTag(m_AttributeSelection, TAGS_SELECTION);
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String debugTipText() {
    return "Outputs debug information to the console.";
  }

  /**
   * Controls whether debugging output will be printed
   *
   * @param debug true if debugging output should be printed
   */
  public void setDebug(boolean debug) {

    b_Debug = debug;
  }

  /**
   * Controls whether debugging output will be printed
   *
   * @return true if debugging output is printed
   */
  public boolean getDebug() {

    return b_Debug;
  }

  /**
   * Removes the attribute with the highest standardised coefficient
   * greater than 1.5 from the selected attributes.
   *
   * @param selectedAttributes an array of flags indicating which 
   * attributes are included in the regression model
   * @param coefficients an array of coefficients for the regression
   * model
   * @return true if an attribute was removed
   */
  private boolean deselectColinearAttributes(boolean [] selectedAttributes,
					     double [] coefficients) {

    double maxSC = 1.5;
    int maxAttr = -1, coeff = 0;
    for (int i = 0; i < selectedAttributes.length; i++) {
      if (selectedAttributes[i]) {
	double SC = Math.abs(coefficients[coeff] * m_StdDevs[i] 
			     / m_ClassStdDev);
	if (SC > maxSC) {
	  maxSC = SC;
	  maxAttr = i;
	}
	coeff++;
      }
    }
    if (maxAttr >= 0) {
      selectedAttributes[maxAttr] = false;
      if (b_Debug) {
	System.out.println("Deselected colinear attribute:" + (maxAttr + 1)
			   + " with standardised coefficient: " + maxSC);
      }
      return true;
    }
    return false;
  }

  /**
   * Performs a greedy search for the best regression model using
   * Akaike's criterion.
   *
   * @throws Exception if regression can't be done
   */
  private void findBestModel() throws Exception {

    // For the weighted case we still use numInstances in
    // the calculation of the Akaike criterion. 
    int numInstances = m_TransformedData.numInstances();

    if (b_Debug) {
      System.out.println((new Instances(m_TransformedData, 0)).toString());
    }

    // Perform a regression for the full model, and remove colinear attributes
    do {
      m_Coefficients = doRegression(m_SelectedAttributes);
    } while (m_EliminateColinearAttributes && 
	     deselectColinearAttributes(m_SelectedAttributes, m_Coefficients));

    // Figure out current number of attributes + 1. (We treat this model
    // as the full model for the Akaike-based methods.)
    int numAttributes = 1;
    for (int i = 0; i < m_SelectedAttributes.length; i++) {
      if (m_SelectedAttributes[i]) {
	numAttributes++;
      }
    }

    double fullMSE = calculateSE(m_SelectedAttributes, m_Coefficients);
    double akaike = (numInstances - numAttributes) + 2 * numAttributes;
    if (b_Debug) {
      System.out.println("Initial Akaike value: " + akaike);
    }

    boolean improved;
    int currentNumAttributes = numAttributes;
    switch (m_AttributeSelection) {

    case SELECTION_GREEDY:

      // Greedy attribute removal
      do {
	boolean [] currentSelected = (boolean []) m_SelectedAttributes.clone();
	improved = false;
	currentNumAttributes--;

	for (int i = 0; i < m_SelectedAttributes.length; i++) {
	  if (currentSelected[i]) {

	    // Calculate the akaike rating without this attribute
	    currentSelected[i] = false;
	    double [] currentCoeffs = doRegression(currentSelected);
	    double currentMSE = calculateSE(currentSelected, currentCoeffs);
	    double currentAkaike = currentMSE / fullMSE 
	      * (numInstances - numAttributes)
	      + 2 * currentNumAttributes;
	    if (b_Debug) {
	      System.out.println("(akaike: " + currentAkaike);
	    }

	    // If it is better than the current best
	    if (currentAkaike < akaike) {
	      if (b_Debug) {
		System.err.println("Removing attribute " + (i + 1)
				   + " improved Akaike: " + currentAkaike);
	      }
	      improved = true;
	      akaike = currentAkaike;
	      System.arraycopy(currentSelected, 0,
			       m_SelectedAttributes, 0,
			       m_SelectedAttributes.length);
	      m_Coefficients = currentCoeffs;
	    }
	    currentSelected[i] = true;
	  }
	}
      } while (improved);
      break;

    case SELECTION_M5:

      // Step through the attributes removing the one with the smallest 
      // standardised coefficient until no improvement in Akaike
      do {
	improved = false;
	currentNumAttributes--;

	// Find attribute with smallest SC
	double minSC = 0;
	int minAttr = -1, coeff = 0;
	for (int i = 0; i < m_SelectedAttributes.length; i++) {
	  if (m_SelectedAttributes[i]) {
	    double SC = Math.abs(m_Coefficients[coeff] * m_StdDevs[i] 
				 / m_ClassStdDev);
	    if ((coeff == 0) || (SC < minSC)) {
	      minSC = SC;
	      minAttr = i;
	    }
	    coeff++;
	  }
	}

	// See whether removing it improves the Akaike score
	if (minAttr >= 0) {
	  m_SelectedAttributes[minAttr] = false;
	  double [] currentCoeffs = doRegression(m_SelectedAttributes);
	  double currentMSE = calculateSE(m_SelectedAttributes, currentCoeffs);
	  double currentAkaike = currentMSE / fullMSE 
	    * (numInstances - numAttributes)
	    + 2 * currentNumAttributes;
	  if (b_Debug) {
	    System.out.println("(akaike: " + currentAkaike);
	  }

	  // If it is better than the current best
	  if (currentAkaike < akaike) {
	    if (b_Debug) {
	      System.err.println("Removing attribute " + (minAttr + 1)
				 + " improved Akaike: " + currentAkaike);
	    }
	    improved = true;
	    akaike = currentAkaike;
	    m_Coefficients = currentCoeffs;
	  } else {
	    m_SelectedAttributes[minAttr] = true;
	  }
	}
      } while (improved);
      break;

    case SELECTION_NONE:
      break;
    }
  }

  /**
   * Calculate the squared error of a regression model on the 
   * training data
   *
   * @param selectedAttributes an array of flags indicating which 
   * attributes are included in the regression model
   * @param coefficients an array of coefficients for the regression
   * model
   * @return the mean squared error on the training data
   * @throws Exception if there is a missing class value in the training
   * data
   */
  private double calculateSE(boolean [] selectedAttributes, 
			      double [] coefficients) throws Exception {

    double mse = 0;
    for (int i = 0; i < m_TransformedData.numInstances(); i++) {
      double prediction = regressionPrediction(m_TransformedData.instance(i),
					       selectedAttributes,
					       coefficients);
      double error = prediction - m_TransformedData.instance(i).classValue();
      mse += error * error;
    }
    return mse;
  }

  /**
   * Calculate the dependent value for a given instance for a
   * given regression model.
   *
   * @param transformedInstance the input instance
   * @param selectedAttributes an array of flags indicating which 
   * attributes are included in the regression model
   * @param coefficients an array of coefficients for the regression
   * model
   * @return the regression value for the instance.
   * @throws Exception if the class attribute of the input instance
   * is not assigned
   */
  private double regressionPrediction(Instance transformedInstance,
				      boolean [] selectedAttributes,
				      double [] coefficients) 
  throws Exception {
    
    double result = 0;
    int column = 0;
    for (int j = 0; j < transformedInstance.numAttributes(); j++) {
      if ((m_ClassIndex != j) 
	  && (selectedAttributes[j])) {
	result += coefficients[column] * transformedInstance.value(j);
	column++;
      }
    }
    result += coefficients[column];
    
    return result;
  }

  /**
   * Calculate a linear regression using the selected attributes
   *
   * @param selectedAttributes an array of booleans where each element
   * is true if the corresponding attribute should be included in the
   * regression.
   * @return an array of coefficients for the linear regression model.
   * @throws Exception if an error occurred during the regression.
   */
  private double [] doRegression(boolean [] selectedAttributes) 
  throws Exception {

    if (b_Debug) {
      System.out.print("doRegression(");
      for (int i = 0; i < selectedAttributes.length; i++) {
	System.out.print(" " + selectedAttributes[i]);
      }
      System.out.println(" )");
    }
    int numAttributes = 0;
    for (int i = 0; i < selectedAttributes.length; i++) {
      if (selectedAttributes[i]) {
	numAttributes++;
      }
    }

    // Check whether there are still attributes left
    Matrix independent = null, dependent = null;
    double[] weights = null;
    if (numAttributes > 0) {
      independent = new Matrix(m_TransformedData.numInstances(), 
			       numAttributes);
      dependent = new Matrix(m_TransformedData.numInstances(), 1);
      for (int i = 0; i < m_TransformedData.numInstances(); i ++) {
	Instance inst = m_TransformedData.instance(i);
	int column = 0;
	for (int j = 0; j < m_TransformedData.numAttributes(); j++) {
	  if (j == m_ClassIndex) {
	    dependent.setElement(i, 0, inst.classValue());
	  } else {
	    if (selectedAttributes[j]) {
	      double value = inst.value(j) - m_Means[j];
	      
	      // We only need to do this if we want to
	      // scale the input
	      if (!m_checksTurnedOff) {
		value /= m_StdDevs[j];
	      }
	      independent.setElement(i, column, value);
	      column++;
	    }
	  }
	}
      }
      
      // Grab instance weights
      weights = new double [m_TransformedData.numInstances()];
      for (int i = 0; i < weights.length; i++) {
	weights[i] = m_TransformedData.instance(i).weight();
      }
    }

    // Compute coefficients (note that we have to treat the
    // intercept separately so that it doesn't get affected
    // by the ridge constant.)
    double[] coefficients = new double[numAttributes + 1];
    if (numAttributes > 0) {
      double[] coeffsWithoutIntercept  =
	independent.regression(dependent, weights, m_Ridge);
      System.arraycopy(coeffsWithoutIntercept, 0, coefficients, 0,
		       numAttributes);
    }
    coefficients[numAttributes] = m_ClassMean;
	   
    // Convert coefficients into original scale
    int column = 0;
    for(int i = 0; i < m_TransformedData.numAttributes(); i++) {
      if ((i != m_TransformedData.classIndex()) &&
	  (selectedAttributes[i])) {

	// We only need to do this if we have scaled the
	// input.
	if (!m_checksTurnedOff) {
	  coefficients[column] /= m_StdDevs[i];
	}

	// We have centred the input
	coefficients[coefficients.length - 1] -= 
	  coefficients[column] * m_Means[i];
	column++;
      }
    }

    return coefficients;
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
   * Generates a linear regression function predictor.
   *
   * @param argv the options
   */
  public static void main(String argv[]) {
    runClassifier(new LinearRegression(), argv);
  }
}
