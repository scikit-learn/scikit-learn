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
 *    SGD.java
 *    Copyright (C) 2009 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.functions;

import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Random;
import java.util.Vector;

import weka.classifiers.AbstractClassifier;
import weka.classifiers.RandomizableClassifier;
import weka.classifiers.UpdateableClassifier;
import weka.core.Capabilities;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.SelectedTag;
import weka.core.Tag;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.filters.Filter;
import weka.filters.unsupervised.attribute.ReplaceMissingValues;
import weka.filters.unsupervised.attribute.Normalize;

/**
 <!-- globalinfo-start -->
 * Implements stochastic gradient descent for learning various linear models (binary class SVM, binary class logistic regression and linear  regression). globally replaces all missing values and transforms nominal attributes into binary ones. It also normalizes all attributes, so the coefficients in the output are based on the normalized data.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -F
 *  Set the loss function to minimize. 0 = hinge loss (SVM), 1 = log loss (logistic regression),
 *  2 = squared loss (regression).
 *  (default = 0)</pre>
 * 
 * <pre> -L
 *  The learning rate. If normalization is
 *  turned off (as it is automatically for streaming data), then the
 *  default learning rate will need to be reduced (try 0.0001).
 *  (default = 0.01).</pre>
 * 
 * <pre> -R &lt;double&gt;
 *  The lambda regularization constant (default = 0.0001)</pre>
 * 
 * <pre> -E &lt;integer&gt;
 *  The number of epochs to perform (batch learning only, default = 500)</pre>
 * 
 * <pre> -N
 *  Don't normalize the data</pre>
 * 
 * <pre> -M
 *  Don't replace missing values</pre>
 * 
 <!-- options-end -->
 *
 * @author Eibe Frank (eibe{[at]}cs{[dot]}waikato{[dot]}ac{[dot]}nz)
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com)
 * @version $Revision: 6584 $
 *
 */
public class SGD extends RandomizableClassifier
  implements UpdateableClassifier,
    OptionHandler {
  
  /** For serialization */
  private static final long serialVersionUID = -3732968666673530290L;
  
  /** Replace missing values */
  protected ReplaceMissingValues m_replaceMissing;
  
  /** Convert nominal attributes to numerically coded binary ones.
   * Uses supervised NominalToBinary in the batch learning case */
  protected Filter m_nominalToBinary;
  
  /** Normalize the training data */
  protected Normalize m_normalize;
  
  /** The regularization parameter */
  protected double m_lambda = 0.0001;
  
  /** The learning rate */
  protected double m_learningRate = 0.01;
  
  /** Stores the weights (+ bias in the last element) */
  protected double[] m_weights;
  
  /** Holds the current iteration number */
  protected double m_t;
  
  /** The number of training instances */
  protected double m_numInstances;

  /**
   *  The number of epochs to perform (batch learning). Total iterations is
   *  m_epochs * num instances 
   */
  protected int m_epochs = 500;
  
  /** 
   * Turn off normalization of the input data. This option gets
   * forced for incremental training.
   */
  protected boolean m_dontNormalize = false;
  
  /**
   *  Turn off global replacement of missing values. Missing values
   *  will be ignored instead. This option gets forced for
   *  incremental training.
   */
  protected boolean m_dontReplaceMissing = false;
  
  /** Holds the header of the training data */
  protected Instances m_data;
  
  /**
   * Returns default capabilities of the classifier.
   *
   * @return      the capabilities of this classifier
   */
  public Capabilities getCapabilities() {
    Capabilities result = super.getCapabilities();
    result.disableAll();
    
    //attributes
    result.enable(Capability.NOMINAL_ATTRIBUTES);
    result.enable(Capability.NUMERIC_ATTRIBUTES);
    result.enable(Capability.MISSING_VALUES);
    
    // class
    result.enable(Capability.BINARY_CLASS);
    result.enable(Capability.NUMERIC_CLASS);
    result.enable(Capability.MISSING_CLASS_VALUES);
    
    // instances
    result.setMinimumNumberInstances(0);
    
    return result;
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String lambdaTipText() {
    return "The regularization constant. (default = 0.0001)";
  }
  
  /**
   * Set the value of lambda to use
   * 
   * @param lambda the value of lambda to use
   */
  public void setLambda(double lambda) {
    m_lambda = lambda;
  }
  
  /**
   * Get the current value of lambda
   * 
   * @return the current value of lambda
   */
  public double getLambda() {
    return m_lambda;
  }
  
  /**
   * Set the learning rate.
   * 
   * @param lr the learning rate to use.
   */
  public void setLearningRate(double lr) {
    m_learningRate = lr;
  }
  
  /**
   * Get the learning rate.
   * 
   * @return the learning rate
   */
  public double getLearningRate() {
    return m_learningRate;
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String learningRateTipText() {
    return "The learning rate. If normalization is turned off " +
    		"(as it is automatically for streaming data), then" +
    		"the default learning rate will need to be reduced (" +
    		"try 0.0001).";
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String epochsTipText() {
    return "The number of epochs to perform (batch learning). " +
    		"The total number of iterations is epochs * num" +
    		" instances.";
  }
  
  /**
   * Set the number of epochs to use
   * 
   * @param e the number of epochs to use
   */
  public void setEpochs(int e) {
    m_epochs = e;
  }
  
  /**
   * Get current number of epochs
   * 
   * @return the current number of epochs
   */
  public int getEpochs() {
    return m_epochs;
  }
  
  /**
   * Turn normalization off/on.
   * 
   * @param m true if normalization is to be disabled.
   */
  public void setDontNormalize(boolean m) {
    m_dontNormalize = m;
  }
  
  /**
   * Get whether normalization has been turned off.
   * 
   * @return true if normalization has been disabled.
   */
  public boolean getDontNormalize() {
    return m_dontNormalize;
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String dontNormalizeTipText() {
    return "Turn normalization off";
  }
  
  /**
   * Turn global replacement of missing values off/on. If turned off,
   * then missing values are effectively ignored.
   * 
   * @param m true if global replacement of missing values is to be
   * turned off.
   */
  public void setDontReplaceMissing(boolean m) {
    m_dontReplaceMissing = m;
  }
  
  /**
   * Get whether global replacement of missing values has been
   * disabled.
   * 
   * @return true if global replacement of missing values has been turned
   * off
   */
  public boolean getDontReplaceMissing() {
    return m_dontReplaceMissing;
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String dontReplaceMissingTipText() {
    return "Turn off global replacement of missing values";
  }
  
  /**
   * Set the loss function to use.
   * 
   * @param function the loss function to use.
   */
  public void setLossFunction(SelectedTag function) {
    if (function.getTags() == TAGS_SELECTION) {
      m_loss = function.getSelectedTag().getID();
    }
  }
  
  /**
   * Get the current loss function.
   * 
   * @return the current loss function.
   */
  public SelectedTag getLossFunction() {
    return new SelectedTag(m_loss, TAGS_SELECTION);
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String lossFunctionTipText() {
    return "The loss function to use. Hinge loss (SVM), " +
                "log loss (logistic regression) or " +
                "squared loss (regression).";
  }
  
  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration<Option> listOptions() {

    Vector<Option> newVector = new Vector<Option>();
    newVector.add(new Option("\tSet the loss function to minimize. 0 = " +
        "hinge loss (SVM), 1 = log loss (logistic regression),\n\t" +
        "2 = squared loss (regression).\n\t(default = 0)", "F", 1, "-F"));
    newVector.add(new Option("\tThe learning rate. If normalization is\n" +
    		"\tturned off (as it is automatically for streaming data), then the\n\t" +
    		"default learning rate will need to be reduced " +
    		"(try 0.0001).\n\t(default = 0.01).", "L", 1, "-L"));
    newVector.add(new Option("\tThe lambda regularization constant " +
    		"(default = 0.0001)",
    		"R", 1, "-R <double>"));
    newVector.add(new Option("\tThe number of epochs to perform (" +
    		"batch learning only, default = 500)", "E", 1,
    		"-E <integer>"));
    newVector.add(new Option("\tDon't normalize the data", "N", 0, "-N"));
    newVector.add(new Option("\tDon't replace missing values", "M", 0, "-M"));
    
    return newVector.elements();
  }
  
  /**
   * 
   * Parses a given list of options. <p/>
   * 
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -F
   *  Set the loss function to minimize. 0 = hinge loss (SVM), 1 = log loss (logistic regression),
   *  2 = squared loss (regression).
   *  (default = 0)</pre>
   * 
   * <pre> -L
   *  The learning rate. If normalization is
   *  turned off (as it is automatically for streaming data), then the
   *  default learning rate will need to be reduced (try 0.0001).
   *  (default = 0.01).</pre>
   * 
   * <pre> -R &lt;double&gt;
   *  The lambda regularization constant (default = 0.0001)</pre>
   * 
   * <pre> -E &lt;integer&gt;
   *  The number of epochs to perform (batch learning only, default = 500)</pre>
   * 
   * <pre> -N
   *  Don't normalize the data</pre>
   * 
   * <pre> -M
   *  Don't replace missing values</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    reset();
    
    super.setOptions(options);
    
    String lossString = Utils.getOption('F', options);
    if (lossString.length() != 0) {
      setLossFunction(new SelectedTag(Integer.parseInt(lossString), 
          TAGS_SELECTION));
    }
    
    String lambdaString = Utils.getOption('R', options);
    if (lambdaString.length() > 0) {
      setLambda(Double.parseDouble(lambdaString));
    }
    
    String learningRateString = Utils.getOption('L', options);
    if (learningRateString.length() > 0) {
      setLearningRate(Double.parseDouble(learningRateString));
    }
    
    String epochsString = Utils.getOption("E", options);
    if (epochsString.length() > 0) {
      setEpochs(Integer.parseInt(epochsString));
    }
    
    setDontNormalize(Utils.getFlag("N", options));
    setDontReplaceMissing(Utils.getFlag('M', options));
  }
  
  /**
   * Gets the current settings of the classifier.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    ArrayList<String> options = new ArrayList<String>();
    
    options.add("-F"); options.add("" + getLossFunction().getSelectedTag().getID());
    options.add("-L"); options.add("" + getLearningRate());
    options.add("-R"); options.add("" + getLambda());
    options.add("-E"); options.add("" + getEpochs());
    if (getDontNormalize()) {
      options.add("-N");
    }
    if (getDontReplaceMissing()) {
      options.add("-M");
    }
    
    return options.toArray(new String[1]);
  }
  
  /**
   * Returns a string describing classifier
   * @return a description suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return "Implements stochastic gradient descent for learning" +
    		" various linear models (binary class SVM, binary class" +
    		" logistic regression and linear" +
    		" regression)." +
    		" globally replaces all missing values and transforms nominal" +
    		" attributes into binary ones. It also normalizes all attributes," +
    		" so the coefficients in the output are based on the normalized" +
    		" data.";  
  }
  
  /**
   * Reset the classifier.
   */
  public void reset() {
    m_t = 1;
    m_weights = null;
  }

  /**
   * Method for building the classifier.
   * 
   * @param data the set of training instances.
   * @throws Exception if the classifier can't be built successfully.
   */
  public void buildClassifier(Instances data) throws Exception {
    reset();
    
    // can classifier handle the data?
    getCapabilities().testWithFail(data);
    
    data = new Instances(data);
    data.deleteWithMissingClass();
    
    // check loss function with respect to class attribute
    if (m_loss != SQUAREDLOSS && data.classAttribute().isNumeric()) {
      throw new Exception("Must use squared loss with a numeric class!");
    }
    
    if (m_loss == SQUAREDLOSS && (data.classAttribute().isNominal() || data.classAttribute().isRanking())) {
      throw new Exception("Must use either the hinge or log loss " +
      		"with a discrete class!");
    }
    
    if (data.numInstances() > 0 && !m_dontReplaceMissing) {
      m_replaceMissing = new ReplaceMissingValues();
      m_replaceMissing.setInputFormat(data);
      data = Filter.useFilter(data, m_replaceMissing);
    }
    
    // check for only numeric attributes
    boolean onlyNumeric = true;
    for (int i = 0; i < data.numAttributes(); i++) {
      if (i != data.classIndex()) {
        if (!data.attribute(i).isNumeric()) {
          onlyNumeric = false;
          break;
        }
      }
    }
    
    if (!onlyNumeric) {
      if (data.numInstances() > 0) {
        m_nominalToBinary = new weka.filters.supervised.attribute.NominalToBinary();
      } else {
        m_nominalToBinary = new weka.filters.unsupervised.attribute.NominalToBinary();
      }
      m_nominalToBinary.setInputFormat(data);
      data = Filter.useFilter(data, m_nominalToBinary);
    }
    
    if (!m_dontNormalize && data.numInstances() > 0) {

      m_normalize = new Normalize();
      m_normalize.setInputFormat(data);
      data = Filter.useFilter(data, m_normalize);
    }
    
    m_numInstances = data.numInstances();
    
    m_weights = new double[data.numAttributes() + 1];
    m_data = new Instances(data, 0);

    if (data.numInstances() > 0) {
      data.randomize(new Random(getSeed())); // randomize the data
      train(data);    
    }
  }
  
  protected static final int HINGE = 0;
  protected static final int LOGLOSS = 1;
  protected static final int SQUAREDLOSS = 2;
  
  /** The current loss function to minimize */
  protected int m_loss = HINGE;
  
  /** Loss functions to choose from */
  public static final Tag [] TAGS_SELECTION = {
    new Tag(HINGE, "Hinge loss (SVM)"),
    new Tag(LOGLOSS, "Log loss (logistic regression)"),
    new Tag(SQUAREDLOSS, "Squared loss (regression)")
  };
  
  protected double dloss(double z) {
    if (m_loss == HINGE) {
      return (z < 1) ? 1 : 0;
    }
    
    if (m_loss == LOGLOSS) {
      // log loss
      if (z < 0) {
        return 1.0 / (Math.exp(z) + 1.0);  
      } else {
        double t = Math.exp(-z);
        return t / (t + 1);
      }
    }
    
    // squared loss
    return z;
  }
  
  private void train(Instances data) throws Exception {
    for (int e = 0; e < m_epochs; e++) {
      for (int i = 0; i < data.numInstances(); i++) {
        updateClassifier(data.instance(i));
      }
    }
  }
  
  protected static double dotProd(Instance inst1, double[] weights, int classIndex) {
    double result = 0;

    int n1 = inst1.numValues();
    int n2 = weights.length - 1; 

    for (int p1 = 0, p2 = 0; p1 < n1 && p2 < n2;) {
      int ind1 = inst1.index(p1);
      int ind2 = p2;
      if (ind1 == ind2) {
        if (ind1 != classIndex && !inst1.isMissingSparse(p1)) {
          result += inst1.valueSparse(p1) * weights[p2];
        }
        p1++;
        p2++;
      } else if (ind1 > ind2) {
        p2++;
      } else {
        p1++;
      }
    }
    return (result);
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

      double wx = dotProd(instance, m_weights, instance.classIndex());
      
      double y;
      double z;
      if (instance.classAttribute().isNominal() || instance.classAttribute().isRanking()) {
        y = (instance.classValue() == 0) ? -1 : 1;
        z = y * (wx + m_weights[m_weights.length - 1]);        
      } else {
        y = instance.classValue();
        z = y - (wx + m_weights[m_weights.length - 1]);
        y = 1;
      }
      
      // Compute multiplier for weight decay
      double multiplier = 1.0;
      if (m_numInstances == 0) {
        multiplier = 1.0 - (m_learningRate * m_lambda) / m_t;
      } else {
        multiplier = 1.0 - (m_learningRate * m_lambda) / m_numInstances;
      }
      for (int i = 0; i < m_weights.length - 1; i++) {
        m_weights[i] *= multiplier;
      }

      // Only need to do the following if the loss is non-zero
      if (m_loss != HINGE || (z < 1)) {          
        
        // Compute Factor for updates
        double factor = m_learningRate * y * dloss(z);

        // Update coefficients for attributes
        int n1 = instance.numValues();
        for (int p1 = 0; p1 < n1; p1++) {
          int indS = instance.index(p1);
          if (indS != instance.classIndex() &&  !instance.isMissingSparse(p1)) {
            m_weights[indS] += factor * instance.valueSparse(p1);
          }
        }
        
        // update the bias
        m_weights[m_weights.length - 1] += factor;
      }
      m_t++;
    }
  }
    
  /**
   * Computes the distribution for a given instance
   *
   * @param instance the instance for which distribution is computed
   * @return the distribution
   * @throws Exception if the distribution can't be computed successfully
   */
  public double[] distributionForInstance(Instance inst) throws Exception {
    double[] result = (inst.classAttribute().isNominal() || inst.classAttribute().isRanking()) 
      ? new double[2]
      : new double[1];

    if (m_replaceMissing != null) {
      m_replaceMissing.input(inst);
      inst = m_replaceMissing.output();
    }

    if (m_nominalToBinary != null) {
      m_nominalToBinary.input(inst);
      inst = m_nominalToBinary.output();
    }

    if (m_normalize != null){
      m_normalize.input(inst);
      inst = m_normalize.output();
    }

    double wx = dotProd(inst, m_weights, inst.classIndex());// * m_wScale;
    double z = (wx + m_weights[m_weights.length - 1]);

    if (inst.classAttribute().isNumeric()) {
      result[0] = z;
      return result;
    }
    
    if (z <= 0) {
      //  z = 0;
      if (m_loss == LOGLOSS) {
        result[0] = 1.0 / (1.0 + Math.exp(z));
        result[1] = 1.0 - result[0];
      } else {
        result[0] = 1;
      }
    } else {
      if (m_loss == LOGLOSS) {
        result[1] = 1.0 / (1.0 + Math.exp(-z));
        result[0] = 1.0 - result[1];
      } else {
        result[1] = 1;
      }
    }
    return result;
  }
  
  /**
   * Prints out the classifier.
   *
   * @return a description of the classifier as a string
   */
  public String toString() {
    if (m_weights == null) {
      return "SGD: No model built yet.\n";
    }
    StringBuffer buff = new StringBuffer();
    buff.append("Loss function: ");
    if (m_loss == HINGE) {
      buff.append("Hinge loss (SVM)\n\n");
    } else if (m_loss == LOGLOSS){
      buff.append("Log loss (logistic regression)\n\n");
    } else {
      buff.append("Squared loss (linear regression)\n\n");
    }
    
    buff.append(m_data.classAttribute().name() + " = \n\n");
    int printed = 0;
    
    for (int i = 0 ; i < m_weights.length - 1; i++) {
      if (i != m_data.classIndex()) {
        if (printed > 0) {
          buff.append(" + ");
        } else {
          buff.append("   ");
        }

        buff.append(Utils.doubleToString(m_weights[i], 12, 4) +
            " " + ((m_normalize != null) ? "(normalized) " : "") 
            + m_data.attribute(i).name() + "\n");

        printed++;
      }
    }
    
    if (m_weights[m_weights.length - 1] > 0) {
      buff.append(" + " + Utils.doubleToString(m_weights[m_weights.length - 1], 12, 4));
    } else {
      buff.append(" - " + Utils.doubleToString(-m_weights[m_weights.length - 1], 12, 4));
    }
    
    return buff.toString();
  }

  /**
   * Returns the revision string.
   * 
   * @return            the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6584 $");
  }
  
  /**
   * Main method for testing this class.
   */
  public static void main(String[] args) {
    runClassifier(new SGD(), args);
  }
}

