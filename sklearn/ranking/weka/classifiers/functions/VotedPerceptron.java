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
 *    VotedPerceptron.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.classifiers.functions;

import weka.classifiers.Classifier;
import weka.classifiers.AbstractClassifier;
import weka.core.Capabilities;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformationHandler;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;
import weka.filters.Filter;
import weka.filters.unsupervised.attribute.NominalToBinary;
import weka.filters.unsupervised.attribute.ReplaceMissingValues;

import java.util.Enumeration;
import java.util.Random;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * Implementation of the voted perceptron algorithm by Freund and Schapire. Globally replaces all missing values, and transforms nominal attributes into binary ones.<br/>
 * <br/>
 * For more information, see:<br/>
 * <br/>
 * Y. Freund, R. E. Schapire: Large margin classification using the perceptron algorithm. In: 11th Annual Conference on Computational Learning Theory, New York, NY, 209-217, 1998.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;inproceedings{Freund1998,
 *    address = {New York, NY},
 *    author = {Y. Freund and R. E. Schapire},
 *    booktitle = {11th Annual Conference on Computational Learning Theory},
 *    pages = {209-217},
 *    publisher = {ACM Press},
 *    title = {Large margin classification using the perceptron algorithm},
 *    year = {1998}
 * }
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -I &lt;int&gt;
 *  The number of iterations to be performed.
 *  (default 1)</pre>
 * 
 * <pre> -E &lt;double&gt;
 *  The exponent for the polynomial kernel.
 *  (default 1)</pre>
 * 
 * <pre> -S &lt;int&gt;
 *  The seed for the random number generation.
 *  (default 1)</pre>
 * 
 * <pre> -M &lt;int&gt;
 *  The maximum number of alterations allowed.
 *  (default 10000)</pre>
 * 
 <!-- options-end -->
 *
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version $Revision: 5928 $ 
 */
public class VotedPerceptron 
  extends AbstractClassifier 
  implements OptionHandler, TechnicalInformationHandler {
  
  /** for serialization */
  static final long serialVersionUID = -1072429260104568698L;
  
  /** The maximum number of alterations to the perceptron */
  private int m_MaxK = 10000;

  /** The number of iterations */
  private int m_NumIterations = 1;

  /** The exponent */
  private double m_Exponent = 1.0;

  /** The actual number of alterations */
  private int m_K = 0;

  /** The training instances added to the perceptron */
  private int[] m_Additions = null;

  /** Addition or subtraction? */
  private boolean[] m_IsAddition = null;

  /** The weights for each perceptron */
  private int[] m_Weights = null;
  
  /** The training instances */
  private Instances m_Train = null;

  /** Seed used for shuffling the dataset */
  private int m_Seed = 1;

  /** The filter used to make attributes numeric. */
  private NominalToBinary m_NominalToBinary;

  /** The filter used to get rid of missing values. */
  private ReplaceMissingValues m_ReplaceMissingValues;

  /**
   * Returns a string describing this classifier
   * @return a description of the classifier suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "Implementation of the voted perceptron algorithm by Freund and "
      + "Schapire. Globally replaces all missing values, and transforms "
      + "nominal attributes into binary ones.\n\n"
      + "For more information, see:\n\n"
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
    result.setValue(Field.AUTHOR, "Y. Freund and R. E. Schapire");
    result.setValue(Field.TITLE, "Large margin classification using the perceptron algorithm");
    result.setValue(Field.BOOKTITLE, "11th Annual Conference on Computational Learning Theory");
    result.setValue(Field.YEAR, "1998");
    result.setValue(Field.PAGES, "209-217");
    result.setValue(Field.PUBLISHER, "ACM Press");
    result.setValue(Field.ADDRESS, "New York, NY");
    
    return result;
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {

    Vector newVector = new Vector(4);

    newVector.addElement(new Option("\tThe number of iterations to be performed.\n"
				    + "\t(default 1)",
				    "I", 1, "-I <int>"));
    newVector.addElement(new Option("\tThe exponent for the polynomial kernel.\n"
				    + "\t(default 1)",
				    "E", 1, "-E <double>"));
    newVector.addElement(new Option("\tThe seed for the random number generation.\n"
				    + "\t(default 1)",
				    "S", 1, "-S <int>"));
    newVector.addElement(new Option("\tThe maximum number of alterations allowed.\n"
				    + "\t(default 10000)",
				    "M", 1, "-M <int>"));

    return newVector.elements();
  }

  /**
   * Parses a given list of options. <p/>
   *
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -I &lt;int&gt;
   *  The number of iterations to be performed.
   *  (default 1)</pre>
   * 
   * <pre> -E &lt;double&gt;
   *  The exponent for the polynomial kernel.
   *  (default 1)</pre>
   * 
   * <pre> -S &lt;int&gt;
   *  The seed for the random number generation.
   *  (default 1)</pre>
   * 
   * <pre> -M &lt;int&gt;
   *  The maximum number of alterations allowed.
   *  (default 10000)</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    
    String iterationsString = Utils.getOption('I', options);
    if (iterationsString.length() != 0) {
      m_NumIterations = Integer.parseInt(iterationsString);
    } else {
      m_NumIterations = 1;
    }
    String exponentsString = Utils.getOption('E', options);
    if (exponentsString.length() != 0) {
      m_Exponent = (new Double(exponentsString)).doubleValue();
    } else {
      m_Exponent = 1.0;
    }
    String seedString = Utils.getOption('S', options);
    if (seedString.length() != 0) {
      m_Seed = Integer.parseInt(seedString);
    } else {
      m_Seed = 1;
    }
    String alterationsString = Utils.getOption('M', options);
    if (alterationsString.length() != 0) {
      m_MaxK = Integer.parseInt(alterationsString);
    } else {
      m_MaxK = 10000;
    }
  }

  /**
   * Gets the current settings of the classifier.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {

    String[] options = new String [8];
    int current = 0;

    options[current++] = "-I"; options[current++] = "" + m_NumIterations;
    options[current++] = "-E"; options[current++] = "" + m_Exponent;
    options[current++] = "-S"; options[current++] = "" + m_Seed;
    options[current++] = "-M"; options[current++] = "" + m_MaxK;
    while (current < options.length) {
      options[current++] = "";
    }
    return options;
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
    result.enable(Capability.BINARY_CLASS);
    result.enable(Capability.MISSING_CLASS_VALUES);

    // instances
    result.setMinimumNumberInstances(0);
    
    return result;
  }

  /**
   * Builds the ensemble of perceptrons.
   *
   * @param insts the data to train the classifier with
   * @throws Exception if something goes wrong during building
   */
  public void buildClassifier(Instances insts) throws Exception {
 
    // can classifier handle the data?
    getCapabilities().testWithFail(insts);

    // remove instances with missing class
    insts = new Instances(insts);
    insts.deleteWithMissingClass();
    
    // Filter data
    m_Train = new Instances(insts);
    m_ReplaceMissingValues = new ReplaceMissingValues();
    m_ReplaceMissingValues.setInputFormat(m_Train);
    m_Train = Filter.useFilter(m_Train, m_ReplaceMissingValues);
    
    m_NominalToBinary = new NominalToBinary();
    m_NominalToBinary.setInputFormat(m_Train);
    m_Train = Filter.useFilter(m_Train, m_NominalToBinary);

    /** Randomize training data */
    m_Train.randomize(new Random(m_Seed));

    /** Make space to store perceptrons */
    m_Additions = new int[m_MaxK + 1];
    m_IsAddition = new boolean[m_MaxK + 1];
    m_Weights = new int[m_MaxK + 1];

    /** Compute perceptrons */
    m_K = 0;
  out:
    for (int it = 0; it < m_NumIterations; it++) {
      for (int i = 0; i < m_Train.numInstances(); i++) {
	Instance inst = m_Train.instance(i);
	if (!inst.classIsMissing()) {
	  int prediction = makePrediction(m_K, inst);
	  int classValue = (int) inst.classValue();
	  if (prediction == classValue) {
	    m_Weights[m_K]++;
	  } else {
	    m_IsAddition[m_K] = (classValue == 1);
	    m_Additions[m_K] = i;
	    m_K++;
	    m_Weights[m_K]++;
	  }
	  if (m_K == m_MaxK) {
	    break out;
	  }
	}
      }
    }
  }

  /**
   * Outputs the distribution for the given output.
   *
   * Pipes output of SVM through sigmoid function.
   * @param inst the instance for which distribution is to be computed
   * @return the distribution
   * @throws Exception if something goes wrong
   */
  public double[] distributionForInstance(Instance inst) throws Exception {

    // Filter instance
    m_ReplaceMissingValues.input(inst);
    m_ReplaceMissingValues.batchFinished();
    inst = m_ReplaceMissingValues.output();

    m_NominalToBinary.input(inst);
    m_NominalToBinary.batchFinished();
    inst = m_NominalToBinary.output();
    
    // Get probabilities
    double output = 0, sumSoFar = 0;
    if (m_K > 0) {
      for (int i = 0; i <= m_K; i++) {
	if (sumSoFar < 0) {
	  output -= m_Weights[i];
	} else {
	  output += m_Weights[i];
	}
	if (m_IsAddition[i]) {
	  sumSoFar += innerProduct(m_Train.instance(m_Additions[i]), inst);
	} else {
	  sumSoFar -= innerProduct(m_Train.instance(m_Additions[i]), inst);
	}
      }
    }
    double[] result = new double[2];
    result[1] = 1 / (1 + Math.exp(-output));
    result[0] = 1 - result[1];

    return result;
  }

  /**
   * Returns textual description of classifier.
   * 
   * @return the model as string
   */
  public String toString() {

    return "VotedPerceptron: Number of perceptrons=" + m_K;
  }
 
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String maxKTipText() {
    return "The maximum number of alterations to the perceptron.";
  }

  /**
   * Get the value of maxK.
   *
   * @return Value of maxK.
   */
  public int getMaxK() {
    
    return m_MaxK;
  }
  
  /**
   * Set the value of maxK.
   *
   * @param v  Value to assign to maxK.
   */
  public void setMaxK(int v) {
    
    m_MaxK = v;
  }
  
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String numIterationsTipText() {
    return "Number of iterations to be performed.";
  }

  /**
   * Get the value of NumIterations.
   *
   * @return Value of NumIterations.
   */
  public int getNumIterations() {
    
    return m_NumIterations;
  }
  
  /**
   * Set the value of NumIterations.
   *
   * @param v  Value to assign to NumIterations.
   */
  public void setNumIterations(int v) {
    
    m_NumIterations = v;
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String exponentTipText() {
    return "Exponent for the polynomial kernel.";
  }

  /**
   * Get the value of exponent.
   *
   * @return Value of exponent.
   */
  public double getExponent() {
    
    return m_Exponent;
  }
  
  /**
   * Set the value of exponent.
   *
   * @param v  Value to assign to exponent.
   */
  public void setExponent(double v) {
    
    m_Exponent = v;
  }
  
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String seedTipText() {
    return "Seed for the random number generator.";
  }

  /**
   * Get the value of Seed.
   *
   * @return Value of Seed.
   */
  public int getSeed() {
    
    return m_Seed;
  }
  
  /**
   * Set the value of Seed.
   *
   * @param v  Value to assign to Seed.
   */
  public void setSeed(int v) {
    
    m_Seed = v;
  }

  /** 
   * Computes the inner product of two instances
   * 
   * @param i1 first instance
   * @param i2 second instance
   * @return the inner product
   * @throws Exception if computation fails
   */
  private double innerProduct(Instance i1, Instance i2) throws Exception {

    // we can do a fast dot product
    double result = 0;
    int n1 = i1.numValues(); int n2 = i2.numValues();
    int classIndex = m_Train.classIndex();
    for (int p1 = 0, p2 = 0; p1 < n1 && p2 < n2;) {
        int ind1 = i1.index(p1);
        int ind2 = i2.index(p2);
        if (ind1 == ind2) {
            if (ind1 != classIndex) {
                result += i1.valueSparse(p1) *
                          i2.valueSparse(p2);
            }
            p1++; p2++;
        } else if (ind1 > ind2) {
            p2++;
        } else {
            p1++;
        }
    }
    result += 1.0;
    
    if (m_Exponent != 1) {
      return Math.pow(result, m_Exponent);
    } else {
      return result;
    }
  }

  /** 
   * Compute a prediction from a perceptron
   * 
   * @param k
   * @param inst the instance to make a prediction for
   * @return the prediction
   * @throws Exception if computation fails
   */
  private int makePrediction(int k, Instance inst) throws Exception {

    double result = 0;
    for (int i = 0; i < k; i++) {
      if (m_IsAddition[i]) {
	result += innerProduct(m_Train.instance(m_Additions[i]), inst);
      } else {
	result -= innerProduct(m_Train.instance(m_Additions[i]), inst);
      }
    }
    if (result < 0) {
      return 0;
    } else {
      return 1;
    }
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
   * Main method.
   * 
   * @param argv the commandline options
   */
  public static void main(String[] argv) {
    runClassifier(new VotedPerceptron(), argv);
  }
}
