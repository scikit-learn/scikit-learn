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
 *    ReliefFAttributeEval.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.attributeSelection;

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
import weka.core.Capabilities.Capability;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;
import weka.core.labelranking.PreferenceAttribute;

import java.util.Enumeration;
import java.util.Random;
import java.util.Vector;

/** 
 <!-- globalinfo-start -->
 * ReliefFAttributeEval :<br/>
 * <br/>
 * Evaluates the worth of an attribute by repeatedly sampling an instance and considering the value of the given attribute for the nearest instance of the same and different class. Can operate on both discrete and continuous class data.<br/>
 * <br/>
 * For more information see:<br/>
 * <br/>
 * Kenji Kira, Larry A. Rendell: A Practical Approach to Feature Selection. In: Ninth International Workshop on Machine Learning, 249-256, 1992.<br/>
 * <br/>
 * Igor Kononenko: Estimating Attributes: Analysis and Extensions of RELIEF. In: European Conference on Machine Learning, 171-182, 1994.<br/>
 * <br/>
 * Marko Robnik-Sikonja, Igor Kononenko: An adaptation of Relief for attribute estimation in regression. In: Fourteenth International Conference on Machine Learning, 296-304, 1997.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;inproceedings{Kira1992,
 *    author = {Kenji Kira and Larry A. Rendell},
 *    booktitle = {Ninth International Workshop on Machine Learning},
 *    editor = {Derek H. Sleeman and Peter Edwards},
 *    pages = {249-256},
 *    publisher = {Morgan Kaufmann},
 *    title = {A Practical Approach to Feature Selection},
 *    year = {1992}
 * }
 * 
 * &#64;inproceedings{Kononenko1994,
 *    author = {Igor Kononenko},
 *    booktitle = {European Conference on Machine Learning},
 *    editor = {Francesco Bergadano and Luc De Raedt},
 *    pages = {171-182},
 *    publisher = {Springer},
 *    title = {Estimating Attributes: Analysis and Extensions of RELIEF},
 *    year = {1994}
 * }
 * 
 * &#64;inproceedings{Robnik-Sikonja1997,
 *    author = {Marko Robnik-Sikonja and Igor Kononenko},
 *    booktitle = {Fourteenth International Conference on Machine Learning},
 *    editor = {Douglas H. Fisher},
 *    pages = {296-304},
 *    publisher = {Morgan Kaufmann},
 *    title = {An adaptation of Relief for attribute estimation in regression},
 *    year = {1997}
 * }
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -M &lt;num instances&gt;
 *  Specify the number of instances to
 *  sample when estimating attributes.
 *  If not specified, then all instances
 *  will be used.</pre>
 * 
 * <pre> -D &lt;seed&gt;
 *  Seed for randomly sampling instances.
 *  (Default = 1)</pre>
 * 
 * <pre> -K &lt;number of neighbours&gt;
 *  Number of nearest neighbours (k) used
 *  to estimate attribute relevances
 *  (Default = 10).</pre>
 * 
 * <pre> -W
 *  Weight nearest neighbours by distance</pre>
 * 
 * <pre> -A &lt;num&gt;
 *  Specify sigma value (used in an exp
 *  function to control how quickly
 *  weights for more distant instances
 *  decrease. Use in conjunction with -W.
 *  Sensible value=1/5 to 1/10 of the
 *  number of nearest neighbours.
 *  (Default = 2)</pre>
 * 
 <!-- options-end -->
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @version $Revision: 5987 $
 */
public class ReliefFAttributeEval
  extends ASEvaluation
  implements AttributeEvaluator,
             OptionHandler, 
             TechnicalInformationHandler {
  
  /** for serialization */
  static final long serialVersionUID = -8422186665795839379L;

  /** The training instances */
  private Instances m_trainInstances;

  /** The class index */
  private int m_classIndex;

  /** The number of attributes */
  private int m_numAttribs;

  /** The number of instances */
  private int m_numInstances;

  /** Numeric class */
  private boolean m_numericClass;

  /** The number of classes if class is nominal */
  private int m_numClasses;

  /** 
   * Used to hold the probability of a different class val given nearest
   * instances (numeric class)
   */
  private double m_ndc;

  /** 
   * Used to hold the prob of different value of an attribute given
   * nearest instances (numeric class case)
   */
  private double[] m_nda;

  /**
   * Used to hold the prob of a different class val and different att
   * val given nearest instances (numeric class case)
   */
  private double[] m_ndcda;

  /** Holds the weights that relief assigns to attributes */
  private double[] m_weights;

  /** Prior class probabilities (discrete class case) */
  private double[] m_classProbs;

  /** 
   * The number of instances to sample when estimating attributes
   * default == -1, use all instances
   */
  private int m_sampleM;

  /** The number of nearest hits/misses */
  private int m_Knn;

  /** k nearest scores + instance indexes for n classes */
  private double[][][] m_karray;

  /** Upper bound for numeric attributes */
  private double[] m_maxArray;

  /** Lower bound for numeric attributes */
  private double[] m_minArray;

  /** Keep track of the farthest instance for each class */
  private double[] m_worst;

  /** Index in the m_karray of the farthest instance for each class */
  private int[] m_index;

  /** Number of nearest neighbours stored of each class */
  private int[] m_stored;
 
  /** Random number seed used for sampling instances */
  private int m_seed;

  /**
   *  used to (optionally) weight nearest neighbours by their distance
   *  from the instance in question. Each entry holds 
   *  exp(-((rank(r_i, i_j)/sigma)^2)) where rank(r_i,i_j) is the rank of
   *  instance i_j in a sequence of instances ordered by the distance
   *  from r_i. sigma is a user defined parameter, default=20
   **/
  private double[] m_weightsByRank;
  private int m_sigma;
  
  /** Weight by distance rather than equal weights */
  private boolean m_weightByDistance;

  /**
   * Constructor
   */
  public ReliefFAttributeEval () {
    resetOptions();
  }

  /**
   * Returns a string describing this attribute evaluator
   * @return a description of the evaluator suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return "ReliefFAttributeEval :\n\nEvaluates the worth of an attribute by "
      +"repeatedly sampling an instance and considering the value of the "
      +"given attribute for the nearest instance of the same and different "
      +"class. Can operate on both discrete and continuous class data.\n\n"
      + "For more information see:\n\n"
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
    TechnicalInformation        result;
    TechnicalInformation        additional;
    
    result = new TechnicalInformation(Type.INPROCEEDINGS);
    result.setValue(Field.AUTHOR, "Kenji Kira and Larry A. Rendell");
    result.setValue(Field.TITLE, "A Practical Approach to Feature Selection");
    result.setValue(Field.BOOKTITLE, "Ninth International Workshop on Machine Learning");
    result.setValue(Field.EDITOR, "Derek H. Sleeman and Peter Edwards");
    result.setValue(Field.YEAR, "1992");
    result.setValue(Field.PAGES, "249-256");
    result.setValue(Field.PUBLISHER, "Morgan Kaufmann");
    
    additional = result.add(Type.INPROCEEDINGS);
    additional.setValue(Field.AUTHOR, "Igor Kononenko");
    additional.setValue(Field.TITLE, "Estimating Attributes: Analysis and Extensions of RELIEF");
    additional.setValue(Field.BOOKTITLE, "European Conference on Machine Learning");
    additional.setValue(Field.EDITOR, "Francesco Bergadano and Luc De Raedt");
    additional.setValue(Field.YEAR, "1994");
    additional.setValue(Field.PAGES, "171-182");
    additional.setValue(Field.PUBLISHER, "Springer");
    
    additional = result.add(Type.INPROCEEDINGS);
    additional.setValue(Field.AUTHOR, "Marko Robnik-Sikonja and Igor Kononenko");
    additional.setValue(Field.TITLE, "An adaptation of Relief for attribute estimation in regression");
    additional.setValue(Field.BOOKTITLE, "Fourteenth International Conference on Machine Learning");
    additional.setValue(Field.EDITOR, "Douglas H. Fisher");
    additional.setValue(Field.YEAR, "1997");
    additional.setValue(Field.PAGES, "296-304");
    additional.setValue(Field.PUBLISHER, "Morgan Kaufmann");
    
    return result;
  }

  /**
   * Returns an enumeration describing the available options.
   * @return an enumeration of all the available options.
   **/
  public Enumeration listOptions () {
    Vector newVector = new Vector(4);
    newVector
      .addElement(new Option("\tSpecify the number of instances to\n" 
                             + "\tsample when estimating attributes.\n" 
                             + "\tIf not specified, then all instances\n" 
                             + "\twill be used.", "M", 1
                             , "-M <num instances>"));
    newVector.
      addElement(new Option("\tSeed for randomly sampling instances.\n" 
                            + "\t(Default = 1)", "D", 1
                            , "-D <seed>"));
    newVector.
      addElement(new Option("\tNumber of nearest neighbours (k) used\n" 
                            + "\tto estimate attribute relevances\n" 
                            + "\t(Default = 10).", "K", 1
                            , "-K <number of neighbours>"));
    newVector.
      addElement(new Option("\tWeight nearest neighbours by distance", "W"
                            , 0, "-W"));
    newVector.
      addElement(new Option("\tSpecify sigma value (used in an exp\n" 
                            + "\tfunction to control how quickly\n" 
                            + "\tweights for more distant instances\n" 
                            + "\tdecrease. Use in conjunction with -W.\n" 
                            + "\tSensible value=1/5 to 1/10 of the\n" 
                            + "\tnumber of nearest neighbours.\n" 
                            + "\t(Default = 2)", "A", 1, "-A <num>"));
    return  newVector.elements();
  }


  /**
   * Parses a given list of options. <p/>
   *
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -M &lt;num instances&gt;
   *  Specify the number of instances to
   *  sample when estimating attributes.
   *  If not specified, then all instances
   *  will be used.</pre>
   * 
   * <pre> -D &lt;seed&gt;
   *  Seed for randomly sampling instances.
   *  (Default = 1)</pre>
   * 
   * <pre> -K &lt;number of neighbours&gt;
   *  Number of nearest neighbours (k) used
   *  to estimate attribute relevances
   *  (Default = 10).</pre>
   * 
   * <pre> -W
   *  Weight nearest neighbours by distance</pre>
   * 
   * <pre> -A &lt;num&gt;
   *  Specify sigma value (used in an exp
   *  function to control how quickly
   *  weights for more distant instances
   *  decrease. Use in conjunction with -W.
   *  Sensible value=1/5 to 1/10 of the
   *  number of nearest neighbours.
   *  (Default = 2)</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions (String[] options)
    throws Exception {
    String optionString;
    resetOptions();
    setWeightByDistance(Utils.getFlag('W', options));
    optionString = Utils.getOption('M', options);

    if (optionString.length() != 0) {
      setSampleSize(Integer.parseInt(optionString));
    }

    optionString = Utils.getOption('D', options);

    if (optionString.length() != 0) {
      setSeed(Integer.parseInt(optionString));
    }

    optionString = Utils.getOption('K', options);

    if (optionString.length() != 0) {
      setNumNeighbours(Integer.parseInt(optionString));
    }

    optionString = Utils.getOption('A', options);

    if (optionString.length() != 0) {
      setWeightByDistance(true); // turn on weighting by distance
      setSigma(Integer.parseInt(optionString));
    }
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String sigmaTipText() {
    return "Set influence of nearest neighbours. Used in an exp function to "
      +"control how quickly weights decrease for more distant instances. "
      +"Use in conjunction with weightByDistance. Sensible values = 1/5 to "
      +"1/10 the number of nearest neighbours.";
  }

  /**
   * Sets the sigma value.
   *
   * @param s the value of sigma (> 0)
   * @throws Exception if s is not positive
   */
  public void setSigma (int s)
    throws Exception {
    if (s <= 0) {
      throw  new Exception("value of sigma must be > 0!");
    }

    m_sigma = s;
  }


  /**
   * Get the value of sigma.
   *
   * @return the sigma value.
   */
  public int getSigma () {
    return  m_sigma;
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String numNeighboursTipText() {
    return "Number of nearest neighbours for attribute estimation.";
  }

  /**
   * Set the number of nearest neighbours
   *
   * @param n the number of nearest neighbours.
   */
  public void setNumNeighbours (int n) {
    m_Knn = n;
  }


  /**
   * Get the number of nearest neighbours
   *
   * @return the number of nearest neighbours
   */
  public int getNumNeighbours () {
    return  m_Knn;
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String seedTipText() {
    return "Random seed for sampling instances.";
  }

  /**
   * Set the random number seed for randomly sampling instances.
   *
   * @param s the random number seed.
   */
  public void setSeed (int s) {
    m_seed = s;
  }


  /**
   * Get the seed used for randomly sampling instances.
   *
   * @return the random number seed.
   */
  public int getSeed () {
    return  m_seed;
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String sampleSizeTipText() {
    return "Number of instances to sample. Default (-1) indicates that all "
      +"instances will be used for attribute estimation.";
  }

  /**
   * Set the number of instances to sample for attribute estimation
   *
   * @param s the number of instances to sample.
   */
  public void setSampleSize (int s) {
    m_sampleM = s;
  }


  /**
   * Get the number of instances used for estimating attributes
   *
   * @return the number of instances.
   */
  public int getSampleSize () {
    return  m_sampleM;
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String weightByDistanceTipText() {
    return "Weight nearest neighbours by their distance.";
  }

  /**
   * Set the nearest neighbour weighting method
   *
   * @param b true nearest neighbours are to be weighted by distance.
   */
  public void setWeightByDistance (boolean b) {
    m_weightByDistance = b;
  }


  /**
   * Get whether nearest neighbours are being weighted by distance
   *
   * @return m_weightByDiffernce
   */
  public boolean getWeightByDistance () {
    return  m_weightByDistance;
  }


  /**
   * Gets the current settings of ReliefFAttributeEval.
   *
   * @return an array of strings suitable for passing to setOptions()
   */
  public String[] getOptions () {
    String[] options = new String[9];
    int current = 0;

    if (getWeightByDistance()) {
      options[current++] = "-W";
    }

    options[current++] = "-M";
    options[current++] = "" + getSampleSize();
    options[current++] = "-D";
    options[current++] = "" + getSeed();
    options[current++] = "-K";
    options[current++] = "" + getNumNeighbours();
    
    if (getWeightByDistance()) {
      options[current++] = "-A";
      options[current++] = "" + getSigma();
    }

    while (current < options.length) {
      options[current++] = "";
    }

    return  options;
  }


  /**
   * Return a description of the ReliefF attribute evaluator.
   *
   * @return a description of the evaluator as a String.
   */
  public String toString () {
    StringBuffer text = new StringBuffer();

    if (m_trainInstances == null) {
      text.append("ReliefF feature evaluator has not been built yet\n");
    }
    else {
      text.append("\tReliefF Ranking Filter");
      text.append("\n\tInstances sampled: ");

      if (m_sampleM == -1) {
        text.append("all\n");
      }
      else {
        text.append(m_sampleM + "\n");
      }

      text.append("\tNumber of nearest neighbours (k): " + m_Knn + "\n");

      if (m_weightByDistance) {
        text.append("\tExponentially decreasing (with distance) " 
                    + "influence for\n" 
                    + "\tnearest neighbours. Sigma: " 
                    + m_sigma + "\n");
      }
      else {
        text.append("\tEqual influence nearest neighbours\n");
      }
    }

    return  text.toString();
  }

  /**
   * Returns the capabilities of this evaluator.
   *
   * @return            the capabilities of this evaluator
   * @see               Capabilities
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
    result.enable(Capability.NOMINAL_CLASS);
    result.enable(Capability.NUMERIC_CLASS);
    result.enable(Capability.DATE_CLASS);
    result.enable(Capability.MISSING_CLASS_VALUES);
    
    return result;
  }

  /**
   * Initializes a ReliefF attribute evaluator. 
   *
   * @param data set of instances serving as training data 
   * @throws Exception if the evaluator has not been 
   * generated successfully
   */
  public void buildEvaluator (Instances data)
    throws Exception {
    
    int z, totalInstances;
    Random r = new Random(m_seed);

    // can evaluator handle data?
    getCapabilities().testWithFail(data);

    m_trainInstances = data;
    m_classIndex = m_trainInstances.classIndex();
    m_numAttribs = m_trainInstances.numAttributes();
    m_numInstances = m_trainInstances.numInstances();

    if (m_trainInstances.attribute(m_classIndex).isNumeric()) {
      m_numericClass = true;
    }
    else {
      m_numericClass = false;
    }

    if (!m_numericClass) {
      m_numClasses = m_trainInstances.attribute(m_classIndex).numValues();
    }
    else {
      m_ndc = 0;
      m_numClasses = 1;
      m_nda = new double[m_numAttribs];
      m_ndcda = new double[m_numAttribs];
    }

    if (m_weightByDistance) // set up the rank based weights
      {
        m_weightsByRank = new double[m_Knn];

        for (int i = 0; i < m_Knn; i++) {
          m_weightsByRank[i] = 
            Math.exp(-((i/(double)m_sigma)*(i/(double)m_sigma)));
        }
      }

    // the final attribute weights
    m_weights = new double[m_numAttribs];
    // num classes (1 for numeric class) knn neighbours, 
    // and 0 = distance, 1 = instance index
    m_karray = new double[m_numClasses][m_Knn][2];

    if (!m_numericClass) {
      m_classProbs = new double[m_numClasses];

      for (int i = 0; i < m_numInstances; i++) {
        m_classProbs[(int)m_trainInstances.instance(i).value(m_classIndex)]++;
      }

      for (int i = 0; i < m_numClasses; i++) {
        m_classProbs[i] /= m_numInstances;
      }
    }

    m_worst = new double[m_numClasses];
    m_index = new int[m_numClasses];
    m_stored = new int[m_numClasses];
    m_minArray = new double[m_numAttribs];
    m_maxArray = new double[m_numAttribs];

    for (int i = 0; i < m_numAttribs; i++) {
      m_minArray[i] = m_maxArray[i] = Double.NaN;
    }

    for (int i = 0; i < m_numInstances; i++) {
      updateMinMax(m_trainInstances.instance(i));
    }
    
    if ((m_sampleM > m_numInstances) || (m_sampleM < 0)) {
      totalInstances = m_numInstances;
    }
    else {
      totalInstances = m_sampleM;
    }

    // process each instance, updating attribute weights
    for (int i = 0; i < totalInstances; i++) {
      if (totalInstances == m_numInstances) {
        z = i;
      }
      else {
        z = r.nextInt()%m_numInstances;
      }

      if (z < 0) {
        z *= -1;
      }

      if (!(m_trainInstances.instance(z).isMissing(m_classIndex))) {
        // first clear the knn and worst index stuff for the classes
        for (int j = 0; j < m_numClasses; j++) {
          m_index[j] = m_stored[j] = 0;

          for (int k = 0; k < m_Knn; k++) {
            m_karray[j][k][0] = m_karray[j][k][1] = 0;
          }
        }

        findKHitMiss(z);

        if (m_numericClass) {
          updateWeightsNumericClass(z);
        }
        else {
          updateWeightsDiscreteClass(z);
        }
      }
    }

    // now scale weights by 1/m_numInstances (nominal class) or
    // calculate weights numeric class
    // System.out.println("num inst:"+m_numInstances+" r_ndc:"+r_ndc);
    for (int i = 0; i < m_numAttribs; i++) {if (i != m_classIndex) {
      if (m_numericClass) {
        m_weights[i] = m_ndcda[i]/m_ndc - 
          ((m_nda[i] - m_ndcda[i])/((double)totalInstances - m_ndc));
      }
      else {
        m_weights[i] *= (1.0/(double)totalInstances);
      }

      //          System.out.println(r_weights[i]);
    }
    }
  }


  /**
   * Evaluates an individual attribute using ReliefF's instance based approach.
   * The actual work is done by buildEvaluator which evaluates all features.
   *
   * @param attribute the index of the attribute to be evaluated
   * @throws Exception if the attribute could not be evaluated
   */
  public double evaluateAttribute (int attribute)
    throws Exception {
    return  m_weights[attribute];
  }


  /**
   * Reset options to their default values
   */
  protected void resetOptions () {
    m_trainInstances = null;
    m_sampleM = -1;
    m_Knn = 10;
    m_sigma = 2;
    m_weightByDistance = false;
    m_seed = 1;
  }


  /**
   * Normalizes a given value of a numeric attribute.
   *
   * @param x the value to be normalized
   * @param i the attribute's index
   * @return the normalized value
   */
  private double norm (double x, int i) {
    if (Double.isNaN(m_minArray[i]) || 
        Utils.eq(m_maxArray[i], m_minArray[i])) {
      return  0;
    }
    else {
      return  (x - m_minArray[i])/(m_maxArray[i] - m_minArray[i]);
    }
  }


  /**
   * Updates the minimum and maximum values for all the attributes
   * based on a new instance.
   *
   * @param instance the new instance
   */
  private void updateMinMax (Instance instance) {
    //    for (int j = 0; j < m_numAttribs; j++) {
    try {
      for (int j = 0; j < instance.numValues(); j++) {
        if ((instance.attributeSparse(j).isNumeric()) && 
            (!instance.isMissingSparse(j))) {
          if (Double.isNaN(m_minArray[instance.index(j)])) {
            m_minArray[instance.index(j)] = instance.valueSparse(j);
            m_maxArray[instance.index(j)] = instance.valueSparse(j);
          }
        else {
          if (instance.valueSparse(j) < m_minArray[instance.index(j)]) {
            m_minArray[instance.index(j)] = instance.valueSparse(j);
          }
          else {
            if (instance.valueSparse(j) > m_maxArray[instance.index(j)]) {
              m_maxArray[instance.index(j)] = instance.valueSparse(j);
            }
          }
        }
        }
      }
    } catch (Exception ex) {
      System.err.println(ex);
      ex.printStackTrace();
    }
  }

  /**
   * Computes the difference between two given attribute
   * values.
   */
  private double difference(int index, double val1, double val2) {

    switch (m_trainInstances.attribute(index).type()) {
    case PreferenceAttribute.RANKING:
    case Attribute.NOMINAL:
      
      // If attribute is nominal
      if (Utils.isMissingValue(val1) || 
          Utils.isMissingValue(val2)) {
        return (1.0 - (1.0/((double)m_trainInstances.
                            attribute(index).numValues())));
      } else if ((int)val1 != (int)val2) {
        return 1;
      } else {
        return 0;
      }
    case Attribute.NUMERIC:

      // If attribute is numeric
      if (Utils.isMissingValue(val1) || 
          Utils.isMissingValue(val2)) {
        if (Utils.isMissingValue(val1) && 
            Utils.isMissingValue(val2)) {
          return 1;
        } else {
          double diff;
          if (Utils.isMissingValue(val2)) {
            diff = norm(val1, index);
          } else {
            diff = norm(val2, index);
          }
          if (diff < 0.5) {
            diff = 1.0 - diff;
          }
          return diff;
        }
      } else {
        return Math.abs(norm(val1, index) - norm(val2, index));
      }
    default:
      return 0;
    }
  }

  /**
   * Calculates the distance between two instances
   *
   * @param first the first instance
   * @param second the second instance
   * @return the distance between the two given instances, between 0 and 1
   */          
  private double distance(Instance first, Instance second) {  

    double distance = 0;
    int firstI, secondI;

    for (int p1 = 0, p2 = 0; 
         p1 < first.numValues() || p2 < second.numValues();) {
      if (p1 >= first.numValues()) {
        firstI = m_trainInstances.numAttributes();
      } else {
        firstI = first.index(p1); 
      }
      if (p2 >= second.numValues()) {
        secondI = m_trainInstances.numAttributes();
      } else {
        secondI = second.index(p2);
      }
      if (firstI == m_trainInstances.classIndex()) {
        p1++; continue;
      } 
      if (secondI == m_trainInstances.classIndex()) {
        p2++; continue;
      } 
      double diff;
      if (firstI == secondI) {
        diff = difference(firstI, 
                          first.valueSparse(p1),
                          second.valueSparse(p2));
        p1++; p2++;
      } else if (firstI > secondI) {
        diff = difference(secondI, 
                          0, second.valueSparse(p2));
        p2++;
      } else {
        diff = difference(firstI, 
                          first.valueSparse(p1), 0);
        p1++;
      }
      //      distance += diff * diff;
      distance += diff;
    }
    
    //    return Math.sqrt(distance / m_NumAttributesUsed);
    return distance;
  }


  /**
   * update attribute weights given an instance when the class is numeric
   *
   * @param instNum the index of the instance to use when updating weights
   */
  private void updateWeightsNumericClass (int instNum) {
    int i, j;
    double temp,temp2;
    int[] tempSorted = null;
    double[] tempDist = null;
    double distNorm = 1.0;
    int firstI, secondI;

    Instance inst = m_trainInstances.instance(instNum);
   
    // sort nearest neighbours and set up normalization variable
    if (m_weightByDistance) {
      tempDist = new double[m_stored[0]];

      for (j = 0, distNorm = 0; j < m_stored[0]; j++) {
        // copy the distances
        tempDist[j] = m_karray[0][j][0];
        // sum normalizer
        distNorm += m_weightsByRank[j];
      }

      tempSorted = Utils.sort(tempDist);
    }

    for (i = 0; i < m_stored[0]; i++) {
      // P diff prediction (class) given nearest instances
      if (m_weightByDistance) {
        temp = difference(m_classIndex, 
                          inst.value(m_classIndex),
                          m_trainInstances.
                          instance((int)m_karray[0][tempSorted[i]][1]).
                          value(m_classIndex));
        temp *= (m_weightsByRank[i]/distNorm);
      }
      else {
        temp = difference(m_classIndex, 
                          inst.value(m_classIndex), 
                          m_trainInstances.
                          instance((int)m_karray[0][i][1]).
                          value(m_classIndex));
        temp *= (1.0/(double)m_stored[0]); // equal influence
      }

      m_ndc += temp;

      Instance cmp;
      cmp = (m_weightByDistance) 
        ? m_trainInstances.instance((int)m_karray[0][tempSorted[i]][1])
        : m_trainInstances.instance((int)m_karray[0][i][1]);
 
      double temp_diffP_diffA_givNearest = 
        difference(m_classIndex, inst.value(m_classIndex),
                   cmp.value(m_classIndex));
      // now the attributes
      for (int p1 = 0, p2 = 0; 
           p1 < inst.numValues() || p2 < cmp.numValues();) {
        if (p1 >= inst.numValues()) {
          firstI = m_trainInstances.numAttributes();
        } else {
          firstI = inst.index(p1); 
        }
        if (p2 >= cmp.numValues()) {
          secondI = m_trainInstances.numAttributes();
        } else {
          secondI = cmp.index(p2);
        }
        if (firstI == m_trainInstances.classIndex()) {
          p1++; continue;
        } 
        if (secondI == m_trainInstances.classIndex()) {
          p2++; continue;
        } 
        temp = 0.0;
        temp2 = 0.0;
      
        if (firstI == secondI) {
          j = firstI;
          temp = difference(j, inst.valueSparse(p1), cmp.valueSparse(p2)); 
          p1++;p2++;
        } else if (firstI > secondI) {
          j = secondI;
          temp = difference(j, 0, cmp.valueSparse(p2));
          p2++;
        } else {
          j = firstI;
          temp = difference(j, inst.valueSparse(p1), 0);
          p1++;
        } 
       
        temp2 = temp_diffP_diffA_givNearest * temp; 
        // P of different prediction and different att value given
        // nearest instances
        if (m_weightByDistance) {
          temp2 *= (m_weightsByRank[i]/distNorm);
        }
        else {
          temp2 *= (1.0/(double)m_stored[0]); // equal influence
        }

        m_ndcda[j] += temp2;
       
        // P of different attribute val given nearest instances
        if (m_weightByDistance) {
          temp *= (m_weightsByRank[i]/distNorm);
        }
        else {
          temp *= (1.0/(double)m_stored[0]); // equal influence
        }

        m_nda[j] += temp;
      }
    }
  }


  /**
   * update attribute weights given an instance when the class is discrete
   *
   * @param instNum the index of the instance to use when updating weights
   */
  private void updateWeightsDiscreteClass (int instNum) {
    int i, j, k;
    int cl;
    double temp_diff, w_norm = 1.0;
    double[] tempDistClass;
    int[] tempSortedClass = null;
    double distNormClass = 1.0;
    double[] tempDistAtt;
    int[][] tempSortedAtt = null;
    double[] distNormAtt = null;
    int firstI, secondI;

    // store the indexes (sparse instances) of non-zero elements
    Instance inst = m_trainInstances.instance(instNum);

    // get the class of this instance
    cl = (int)m_trainInstances.instance(instNum).value(m_classIndex);

    // sort nearest neighbours and set up normalization variables
    if (m_weightByDistance) {
      // do class (hits) first
      // sort the distances
      tempDistClass = new double[m_stored[cl]];

      for (j = 0, distNormClass = 0; j < m_stored[cl]; j++) {
        // copy the distances
        tempDistClass[j] = m_karray[cl][j][0];
        // sum normalizer
        distNormClass += m_weightsByRank[j];
      }

      tempSortedClass = Utils.sort(tempDistClass);
      // do misses (other classes)
      tempSortedAtt = new int[m_numClasses][1];
      distNormAtt = new double[m_numClasses];

      for (k = 0; k < m_numClasses; k++) {
        if (k != cl) // already done cl
          {
            // sort the distances
            tempDistAtt = new double[m_stored[k]];

            for (j = 0, distNormAtt[k] = 0; j < m_stored[k]; j++) {
              // copy the distances
              tempDistAtt[j] = m_karray[k][j][0];
              // sum normalizer
              distNormAtt[k] += m_weightsByRank[j];
            }

            tempSortedAtt[k] = Utils.sort(tempDistAtt);
          }
      }
    }

    if (m_numClasses > 2) {
      // the amount of probability space left after removing the
      // probability of this instance's class value
      w_norm = (1.0 - m_classProbs[cl]);
    }
    
    // do the k nearest hits of the same class
    for (j = 0, temp_diff = 0.0; j < m_stored[cl]; j++) {
      Instance cmp;
      cmp = (m_weightByDistance) 
        ? m_trainInstances.
        instance((int)m_karray[cl][tempSortedClass[j]][1])
        : m_trainInstances.instance((int)m_karray[cl][j][1]);

      for (int p1 = 0, p2 = 0; 
           p1 < inst.numValues() || p2 < cmp.numValues();) {
        if (p1 >= inst.numValues()) {
          firstI = m_trainInstances.numAttributes();
        } else {
          firstI = inst.index(p1); 
        }
        if (p2 >= cmp.numValues()) {
          secondI = m_trainInstances.numAttributes();
        } else {
          secondI = cmp.index(p2);
        }
        if (firstI == m_trainInstances.classIndex()) {
          p1++; continue;
        } 
        if (secondI == m_trainInstances.classIndex()) {
          p2++; continue;
        } 
        if (firstI == secondI) {
          i = firstI;
          temp_diff = difference(i, inst.valueSparse(p1), 
                                 cmp.valueSparse(p2)); 
          p1++;p2++;
        } else if (firstI > secondI) {
          i = secondI;
          temp_diff = difference(i, 0, cmp.valueSparse(p2));
          p2++;
        } else {
          i = firstI;
          temp_diff = difference(i, inst.valueSparse(p1), 0);
          p1++;
        } 
        
        if (m_weightByDistance) {
          temp_diff *=
            (m_weightsByRank[j]/distNormClass);
        } else {
          if (m_stored[cl] > 0) {
            temp_diff /= (double)m_stored[cl];
          }
        }
        m_weights[i] -= temp_diff;

      }
    }
      

    // now do k nearest misses from each of the other classes
    temp_diff = 0.0;

    for (k = 0; k < m_numClasses; k++) {
      if (k != cl) // already done cl
        {
          for (j = 0; j < m_stored[k]; j++) {
            Instance cmp;
            cmp = (m_weightByDistance) 
              ? m_trainInstances.
              instance((int)m_karray[k][tempSortedAtt[k][j]][1])
              : m_trainInstances.instance((int)m_karray[k][j][1]);
        
            for (int p1 = 0, p2 = 0; 
                 p1 < inst.numValues() || p2 < cmp.numValues();) {
              if (p1 >= inst.numValues()) {
                firstI = m_trainInstances.numAttributes();
              } else {
                firstI = inst.index(p1); 
              }
              if (p2 >= cmp.numValues()) {
                secondI = m_trainInstances.numAttributes();
              } else {
                secondI = cmp.index(p2);
              }
              if (firstI == m_trainInstances.classIndex()) {
                p1++; continue;
              } 
              if (secondI == m_trainInstances.classIndex()) {
                p2++; continue;
              } 
              if (firstI == secondI) {
                i = firstI;
                temp_diff = difference(i, inst.valueSparse(p1), 
                                       cmp.valueSparse(p2)); 
                p1++;p2++;
              } else if (firstI > secondI) {
                i = secondI;
                temp_diff = difference(i, 0, cmp.valueSparse(p2));
                p2++;
              } else {
                i = firstI;
                temp_diff = difference(i, inst.valueSparse(p1), 0);
                p1++;
              } 

              if (m_weightByDistance) {
                temp_diff *=
                  (m_weightsByRank[j]/distNormAtt[k]);
              }
              else {
                if (m_stored[k] > 0) {
                  temp_diff /= (double)m_stored[k];
                }
              }
              if (m_numClasses > 2) {
                m_weights[i] += ((m_classProbs[k]/w_norm)*temp_diff);
              } else {
                m_weights[i] += temp_diff;
              }
            }
          }
        }
    }
  }


  /**
   * Find the K nearest instances to supplied instance if the class is numeric,
   * or the K nearest Hits (same class) and Misses (K from each of the other
   * classes) if the class is discrete.
   *
   * @param instNum the index of the instance to find nearest neighbours of
   */
  private void findKHitMiss (int instNum) {
    int i, j;
    int cl;
    double ww;
    double temp_diff = 0.0;
    Instance thisInst = m_trainInstances.instance(instNum);

    for (i = 0; i < m_numInstances; i++) {
      if (i != instNum) {
        Instance cmpInst = m_trainInstances.instance(i);
        temp_diff = distance(cmpInst, thisInst);

        // class of this training instance or 0 if numeric
        if (m_numericClass) {
          cl = 0;
        }
        else {
          cl = (int)m_trainInstances.instance(i).value(m_classIndex);
        }

        // add this diff to the list for the class of this instance
        if (m_stored[cl] < m_Knn) {
          m_karray[cl][m_stored[cl]][0] = temp_diff;
          m_karray[cl][m_stored[cl]][1] = i;
          m_stored[cl]++;

          // note the worst diff for this class
          for (j = 0, ww = -1.0; j < m_stored[cl]; j++) {
            if (m_karray[cl][j][0] > ww) {
              ww = m_karray[cl][j][0];
              m_index[cl] = j;
            }
          }

          m_worst[cl] = ww;
        }
        else 
          /* if we already have stored knn for this class then check to
             see if this instance is better than the worst */
          {
            if (temp_diff < m_karray[cl][m_index[cl]][0]) {
              m_karray[cl][m_index[cl]][0] = temp_diff;
              m_karray[cl][m_index[cl]][1] = i;

              for (j = 0, ww = -1.0; j < m_stored[cl]; j++) {
                if (m_karray[cl][j][0] > ww) {
                  ww = m_karray[cl][j][0];
                  m_index[cl] = j;
                }
              }

              m_worst[cl] = ww;
            }
          }
      }
    }
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5987 $");
  }

  // ============
  // Test method.
  // ============
  /**
   * Main method for testing this class.
   *
   * @param args the options
   */
  public static void main (String[] args) {
    runEvaluator(new ReliefFAttributeEval(), args);
  }
}
