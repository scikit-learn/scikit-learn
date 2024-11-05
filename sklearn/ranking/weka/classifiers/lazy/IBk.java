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
 *    IBk.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.lazy;

import weka.classifiers.Classifier;
import weka.classifiers.AbstractClassifier;
import weka.classifiers.UpdateableClassifier;
import weka.classifiers.rules.ZeroR;
import weka.core.Attribute;
import weka.core.Capabilities;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.labelranking.PreferenceAttribute;
import weka.core.neighboursearch.LinearNNSearch;
import weka.core.neighboursearch.NearestNeighbourSearch;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.SelectedTag;
import weka.core.Tag;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformationHandler;
import weka.core.Utils;
import weka.core.WeightedInstancesHandler;
import weka.core.Capabilities.Capability;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;
import weka.core.AdditionalMeasureProducer;

import java.util.Enumeration;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * K-nearest neighbours classifier. Can select appropriate value of K based on cross-validation. Can also do distance weighting.<br/>
 * <br/>
 * For more information, see<br/>
 * <br/>
 * D. Aha, D. Kibler (1991). Instance-based learning algorithms. Machine Learning. 6:37-66.
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;article{Aha1991,
 *    author = {D. Aha and D. Kibler},
 *    journal = {Machine Learning},
 *    pages = {37-66},
 *    title = {Instance-based learning algorithms},
 *    volume = {6},
 *    year = {1991}
 * }
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -I
 *  Weight neighbours by the inverse of their distance
 *  (use when k &gt; 1)</pre>
 * 
 * <pre> -F
 *  Weight neighbours by 1 - their distance
 *  (use when k &gt; 1)</pre>
 * 
 * <pre> -K &lt;number of neighbors&gt;
 *  Number of nearest neighbours (k) used in classification.
 *  (Default = 1)</pre>
 * 
 * <pre> -E
 *  Minimise mean squared error rather than mean absolute
 *  error when using -X option with numeric prediction.</pre>
 * 
 * <pre> -W &lt;window size&gt;
 *  Maximum number of training instances maintained.
 *  Training instances are dropped FIFO. (Default = no window)</pre>
 * 
 * <pre> -X
 *  Select the number of nearest neighbours between 1
 *  and the k value specified using hold-one-out evaluation
 *  on the training data (use when k &gt; 1)</pre>
 * 
 * <pre> -A
 *  The nearest neighbour search algorithm to use (default: weka.core.neighboursearch.LinearNNSearch).
 * </pre>
 * 
 <!-- options-end -->
 *
 * @author Stuart Inglis (singlis@cs.waikato.ac.nz)
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version $Revision: 6572 $
 */
public class IBk 
  extends AbstractClassifier 
  implements OptionHandler, UpdateableClassifier, WeightedInstancesHandler,
             TechnicalInformationHandler, AdditionalMeasureProducer {

  /** for serialization. */
  static final long serialVersionUID = -3080186098777067172L;

  /** The training instances used for classification. */
  protected Instances m_Train;

  /** The number of class values (or 1 if predicting numeric). */
  protected int m_NumClasses;

  /** The class attribute type. */
  protected int m_ClassType;

  /** The number of neighbours to use for classification (currently). */
  protected int m_kNN;

  /**
   * The value of kNN provided by the user. This may differ from
   * m_kNN if cross-validation is being used.
   */
  protected int m_kNNUpper;

  /**
   * Whether the value of k selected by cross validation has
   * been invalidated by a change in the training instances.
   */
  protected boolean m_kNNValid;

  /**
   * The maximum number of training instances allowed. When
   * this limit is reached, old training instances are removed,
   * so the training data is "windowed". Set to 0 for unlimited
   * numbers of instances.
   */
  protected int m_WindowSize;

  /** Whether the neighbours should be distance-weighted. */
  protected int m_DistanceWeighting;

  /** Whether to select k by cross validation. */
  protected boolean m_CrossValidate;

  /**
   * Whether to minimise mean squared error rather than mean absolute
   * error when cross-validating on numeric prediction tasks.
   */
  protected boolean m_MeanSquared;
  
  /** Default ZeroR model to use when there are no training instances */
  protected ZeroR m_defaultModel;

  /** no weighting. */
  public static final int WEIGHT_NONE = 1;
  /** weight by 1/distance. */
  public static final int WEIGHT_INVERSE = 2;
  /** weight by 1-distance. */
  public static final int WEIGHT_SIMILARITY = 4;
  /** possible instance weighting methods. */
  public static final Tag [] TAGS_WEIGHTING = {
    new Tag(WEIGHT_NONE, "No distance weighting"),
    new Tag(WEIGHT_INVERSE, "Weight by 1/distance"),
    new Tag(WEIGHT_SIMILARITY, "Weight by 1-distance")
  };
  
  /** for nearest-neighbor search. */
  protected NearestNeighbourSearch m_NNSearch = new LinearNNSearch();

  /** The number of attributes the contribute to a prediction. */
  protected double m_NumAttributesUsed;
  
  /**
   * IBk classifier. Simple instance-based learner that uses the class
   * of the nearest k training instances for the class of the test
   * instances.
   *
   * @param k the number of nearest neighbors to use for prediction
   */
  public IBk(int k) {

    init();
    setKNN(k);
  }  

  /**
   * IB1 classifer. Instance-based learner. Predicts the class of the
   * single nearest training instance for each test instance.
   */
  public IBk() {

    init();
  }
  
  /**
   * Returns a string describing classifier.
   * @return a description suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {

    return  "K-nearest neighbours classifier. Can "
      + "select appropriate value of K based on cross-validation. Can also do "
      + "distance weighting.\n\n"
      + "For more information, see\n\n"
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
    
    result = new TechnicalInformation(Type.ARTICLE);
    result.setValue(Field.AUTHOR, "D. Aha and D. Kibler");
    result.setValue(Field.YEAR, "1991");
    result.setValue(Field.TITLE, "Instance-based learning algorithms");
    result.setValue(Field.JOURNAL, "Machine Learning");
    result.setValue(Field.VOLUME, "6");
    result.setValue(Field.PAGES, "37-66");
    
    return result;
  }

  /**
   * Returns the tip text for this property.
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String KNNTipText() {
    return "The number of neighbours to use.";
  }
  
  /**
   * Set the number of neighbours the learner is to use.
   *
   * @param k the number of neighbours.
   */
  public void setKNN(int k) {
    m_kNN = k;
    m_kNNUpper = k;
    m_kNNValid = false;
  }

  /**
   * Gets the number of neighbours the learner will use.
   *
   * @return the number of neighbours.
   */
  public int getKNN() {

    return m_kNN;
  }

  /**
   * Returns the tip text for this property.
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String windowSizeTipText() {
    return "Gets the maximum number of instances allowed in the training " +
      "pool. The addition of new instances above this value will result " +
      "in old instances being removed. A value of 0 signifies no limit " +
      "to the number of training instances.";
  }
  
  /**
   * Gets the maximum number of instances allowed in the training
   * pool. The addition of new instances above this value will result
   * in old instances being removed. A value of 0 signifies no limit
   * to the number of training instances.
   *
   * @return Value of WindowSize.
   */
  public int getWindowSize() {
    
    return m_WindowSize;
  }
  
  /**
   * Sets the maximum number of instances allowed in the training
   * pool. The addition of new instances above this value will result
   * in old instances being removed. A value of 0 signifies no limit
   * to the number of training instances.
   *
   * @param newWindowSize Value to assign to WindowSize.
   */
  public void setWindowSize(int newWindowSize) {
    
    m_WindowSize = newWindowSize;
  }
  
  /**
   * Returns the tip text for this property.
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String distanceWeightingTipText() {

    return "Gets the distance weighting method used.";
  }
  
  /**
   * Gets the distance weighting method used. Will be one of
   * WEIGHT_NONE, WEIGHT_INVERSE, or WEIGHT_SIMILARITY
   *
   * @return the distance weighting method used.
   */
  public SelectedTag getDistanceWeighting() {

    return new SelectedTag(m_DistanceWeighting, TAGS_WEIGHTING);
  }
  
  /**
   * Sets the distance weighting method used. Values other than
   * WEIGHT_NONE, WEIGHT_INVERSE, or WEIGHT_SIMILARITY will be ignored.
   *
   * @param newMethod the distance weighting method to use
   */
  public void setDistanceWeighting(SelectedTag newMethod) {
    
    if (newMethod.getTags() == TAGS_WEIGHTING) {
      m_DistanceWeighting = newMethod.getSelectedTag().getID();
    }
  }
  
  /**
   * Returns the tip text for this property.
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String meanSquaredTipText() {

    return "Whether the mean squared error is used rather than mean "
      + "absolute error when doing cross-validation for regression problems.";
  }

  /**
   * Gets whether the mean squared error is used rather than mean
   * absolute error when doing cross-validation.
   *
   * @return true if so.
   */
  public boolean getMeanSquared() {
    
    return m_MeanSquared;
  }
  
  /**
   * Sets whether the mean squared error is used rather than mean
   * absolute error when doing cross-validation.
   *
   * @param newMeanSquared true if so.
   */
  public void setMeanSquared(boolean newMeanSquared) {
    
    m_MeanSquared = newMeanSquared;
  }
  
  /**
   * Returns the tip text for this property.
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String crossValidateTipText() {

    return "Whether hold-one-out cross-validation will be used " +
      "to select the best k value.";
  }
  
  /**
   * Gets whether hold-one-out cross-validation will be used
   * to select the best k value.
   *
   * @return true if cross-validation will be used.
   */
  public boolean getCrossValidate() {
    
    return m_CrossValidate;
  }
  
  /**
   * Sets whether hold-one-out cross-validation will be used
   * to select the best k value.
   *
   * @param newCrossValidate true if cross-validation should be used.
   */
  public void setCrossValidate(boolean newCrossValidate) {
    
    m_CrossValidate = newCrossValidate;
  }

  /**
   * Returns the tip text for this property.
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String nearestNeighbourSearchAlgorithmTipText() {
    return "The nearest neighbour search algorithm to use " +
    	   "(Default: weka.core.neighboursearch.LinearNNSearch).";
  }
  
  /**
   * Returns the current nearestNeighbourSearch algorithm in use.
   * @return the NearestNeighbourSearch algorithm currently in use.
   */
  public NearestNeighbourSearch getNearestNeighbourSearchAlgorithm() {
    return m_NNSearch;
  }
  
  /**
   * Sets the nearestNeighbourSearch algorithm to be used for finding nearest
   * neighbour(s).
   * @param nearestNeighbourSearchAlgorithm - The NearestNeighbourSearch class.
   */
  public void setNearestNeighbourSearchAlgorithm(NearestNeighbourSearch nearestNeighbourSearchAlgorithm) {
    m_NNSearch = nearestNeighbourSearchAlgorithm;
  }
   
  /**
   * Get the number of training instances the classifier is currently using.
   * 
   * @return the number of training instances the classifier is currently using
   */
  public int getNumTraining() {

    return m_Train.numInstances();
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
    result.enable(Capability.NOMINAL_CLASS);
    result.enable(Capability.NUMERIC_CLASS);
    result.enable(Capability.DATE_CLASS);
    result.enable(Capability.MISSING_CLASS_VALUES);

    // instances
    result.setMinimumNumberInstances(0);
    
    return result;
  }
  
  /**
   * Generates the classifier.
   *
   * @param instances set of instances serving as training data 
   * @throws Exception if the classifier has not been generated successfully
   */
  public void buildClassifier(Instances instances) throws Exception {
    
    // can classifier handle the data?
    getCapabilities().testWithFail(instances);

    // remove instances with missing class
    instances = new Instances(instances);
    instances.deleteWithMissingClass();
    
    m_NumClasses = instances.numClasses();
    m_ClassType = instances.classAttribute().type();
    m_Train = new Instances(instances, 0, instances.numInstances());

    // Throw away initial instances until within the specified window size
    if ((m_WindowSize > 0) && (instances.numInstances() > m_WindowSize)) {
      m_Train = new Instances(m_Train, 
			      m_Train.numInstances()-m_WindowSize, 
			      m_WindowSize);
    }

    m_NumAttributesUsed = 0.0;
    for (int i = 0; i < m_Train.numAttributes(); i++) {
      if ((i != m_Train.classIndex()) && 
	  (m_Train.attribute(i).isNominal() || m_Train.attribute(i).isRanking() ||
	   m_Train.attribute(i).isNumeric())) {
	m_NumAttributesUsed += 1.0;
      }
    }
    
    m_NNSearch.setInstances(m_Train);

    // Invalidate any currently cross-validation selected k
    m_kNNValid = false;
    
    m_defaultModel = new ZeroR();
    m_defaultModel.buildClassifier(instances);
  }

  /**
   * Adds the supplied instance to the training set.
   *
   * @param instance the instance to add
   * @throws Exception if instance could not be incorporated
   * successfully
   */
  public void updateClassifier(Instance instance) throws Exception {

    if (m_Train.equalHeaders(instance.dataset()) == false) {
      throw new Exception("Incompatible instance types\n" + m_Train.equalHeadersMsg(instance.dataset()));
    }
    if (instance.classIsMissing()) {
      return;
    }

    m_Train.add(instance);
    m_NNSearch.update(instance);
    m_kNNValid = false;
    if ((m_WindowSize > 0) && (m_Train.numInstances() > m_WindowSize)) {
      boolean deletedInstance=false;
      while (m_Train.numInstances() > m_WindowSize) {
	m_Train.delete(0);
        deletedInstance=true;
      }
      //rebuild datastructure KDTree currently can't delete
      if(deletedInstance==true)
        m_NNSearch.setInstances(m_Train);
    }
  }

  /**
   * Calculates the class membership probabilities for the given test instance.
   *
   * @param instance the instance to be classified
   * @return predicted class probability distribution
   * @throws Exception if an error occurred during the prediction
   */
  public double [] distributionForInstance(Instance instance) throws Exception {

    if (m_Train.numInstances() == 0) {
      //throw new Exception("No training instances!");
      return m_defaultModel.distributionForInstance(instance);
    }
    if ((m_WindowSize > 0) && (m_Train.numInstances() > m_WindowSize)) {
      m_kNNValid = false;
      boolean deletedInstance=false;
      while (m_Train.numInstances() > m_WindowSize) {
	m_Train.delete(0);
      }
      //rebuild datastructure KDTree currently can't delete
      if(deletedInstance==true)
        m_NNSearch.setInstances(m_Train);
    }

    // Select k by cross validation
    if (!m_kNNValid && (m_CrossValidate) && (m_kNNUpper >= 1)) {
      crossValidate();
    }

    m_NNSearch.addInstanceInfo(instance);

    Instances neighbours = m_NNSearch.kNearestNeighbours(instance, m_kNN);
    double [] distances = m_NNSearch.getDistances();
    double [] distribution = makeDistribution( neighbours, distances );

    return distribution;
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {

    Vector newVector = new Vector(8);

    newVector.addElement(new Option(
	      "\tWeight neighbours by the inverse of their distance\n"+
	      "\t(use when k > 1)",
	      "I", 0, "-I"));
    newVector.addElement(new Option(
	      "\tWeight neighbours by 1 - their distance\n"+
	      "\t(use when k > 1)",
	      "F", 0, "-F"));
    newVector.addElement(new Option(
	      "\tNumber of nearest neighbours (k) used in classification.\n"+
	      "\t(Default = 1)",
	      "K", 1,"-K <number of neighbors>"));
    newVector.addElement(new Option(
          "\tMinimise mean squared error rather than mean absolute\n"+
	      "\terror when using -X option with numeric prediction.",
	      "E", 0,"-E"));
    newVector.addElement(new Option(
          "\tMaximum number of training instances maintained.\n"+
	      "\tTraining instances are dropped FIFO. (Default = no window)",
	      "W", 1,"-W <window size>"));
    newVector.addElement(new Option(
	      "\tSelect the number of nearest neighbours between 1\n"+
	      "\tand the k value specified using hold-one-out evaluation\n"+
	      "\ton the training data (use when k > 1)",
	      "X", 0,"-X"));
    newVector.addElement(new Option(
	      "\tThe nearest neighbour search algorithm to use "+
          "(default: weka.core.neighboursearch.LinearNNSearch).\n",
	      "A", 0, "-A"));

    return newVector.elements();
  }

  /**
   * Parses a given list of options. <p/>
   *
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -I
   *  Weight neighbours by the inverse of their distance
   *  (use when k &gt; 1)</pre>
   * 
   * <pre> -F
   *  Weight neighbours by 1 - their distance
   *  (use when k &gt; 1)</pre>
   * 
   * <pre> -K &lt;number of neighbors&gt;
   *  Number of nearest neighbours (k) used in classification.
   *  (Default = 1)</pre>
   * 
   * <pre> -E
   *  Minimise mean squared error rather than mean absolute
   *  error when using -X option with numeric prediction.</pre>
   * 
   * <pre> -W &lt;window size&gt;
   *  Maximum number of training instances maintained.
   *  Training instances are dropped FIFO. (Default = no window)</pre>
   * 
   * <pre> -X
   *  Select the number of nearest neighbours between 1
   *  and the k value specified using hold-one-out evaluation
   *  on the training data (use when k &gt; 1)</pre>
   * 
   * <pre> -A
   *  The nearest neighbour search algorithm to use (default: weka.core.neighboursearch.LinearNNSearch).
   * </pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    
    String knnString = Utils.getOption('K', options);
    if (knnString.length() != 0) {
      setKNN(Integer.parseInt(knnString));
    } else {
      setKNN(1);
    }
    String windowString = Utils.getOption('W', options);
    if (windowString.length() != 0) {
      setWindowSize(Integer.parseInt(windowString));
    } else {
      setWindowSize(0);
    }
    if (Utils.getFlag('I', options)) {
      setDistanceWeighting(new SelectedTag(WEIGHT_INVERSE, TAGS_WEIGHTING));
    } else if (Utils.getFlag('F', options)) {
      setDistanceWeighting(new SelectedTag(WEIGHT_SIMILARITY, TAGS_WEIGHTING));
    } else {
      setDistanceWeighting(new SelectedTag(WEIGHT_NONE, TAGS_WEIGHTING));
    }
    setCrossValidate(Utils.getFlag('X', options));
    setMeanSquared(Utils.getFlag('E', options));

    String nnSearchClass = Utils.getOption('A', options);
    if(nnSearchClass.length() != 0) {
      String nnSearchClassSpec[] = Utils.splitOptions(nnSearchClass);
      if(nnSearchClassSpec.length == 0) { 
        throw new Exception("Invalid NearestNeighbourSearch algorithm " +
                            "specification string."); 
      }
      String className = nnSearchClassSpec[0];
      nnSearchClassSpec[0] = "";

      setNearestNeighbourSearchAlgorithm( (NearestNeighbourSearch)
                  Utils.forName( NearestNeighbourSearch.class, 
                                 className, 
                                 nnSearchClassSpec)
                                        );
    }
    else 
      this.setNearestNeighbourSearchAlgorithm(new LinearNNSearch());
    
    Utils.checkForRemainingOptions(options);
  }

  /**
   * Gets the current settings of IBk.
   *
   * @return an array of strings suitable for passing to setOptions()
   */
  public String [] getOptions() {

    String [] options = new String [11];
    int current = 0;
    options[current++] = "-K"; options[current++] = "" + getKNN();
    options[current++] = "-W"; options[current++] = "" + m_WindowSize;
    if (getCrossValidate()) {
      options[current++] = "-X";
    }
    if (getMeanSquared()) {
      options[current++] = "-E";
    }
    if (m_DistanceWeighting == WEIGHT_INVERSE) {
      options[current++] = "-I";
    } else if (m_DistanceWeighting == WEIGHT_SIMILARITY) {
      options[current++] = "-F";
    }

    options[current++] = "-A";
    options[current++] = m_NNSearch.getClass().getName()+" "+Utils.joinOptions(m_NNSearch.getOptions()); 
    
    while (current < options.length) {
      options[current++] = "";
    }
    
    return options;
  }

  /**
   * Returns an enumeration of the additional measure names 
   * produced by the neighbour search algorithm, plus the chosen K in case
   * cross-validation is enabled.
   * 
   * @return an enumeration of the measure names
   */
  public Enumeration enumerateMeasures() {
    if (m_CrossValidate) {
      Enumeration enm = m_NNSearch.enumerateMeasures();
      Vector measures = new Vector();
      while (enm.hasMoreElements())
	measures.add(enm.nextElement());
      measures.add("measureKNN");
      return measures.elements();
    }
    else {
      return m_NNSearch.enumerateMeasures();
    }
  }
  
  /**
   * Returns the value of the named measure from the 
   * neighbour search algorithm, plus the chosen K in case
   * cross-validation is enabled.
   * 
   * @param additionalMeasureName the name of the measure to query for its value
   * @return the value of the named measure
   * @throws IllegalArgumentException if the named measure is not supported
   */
  public double getMeasure(String additionalMeasureName) {
    if (additionalMeasureName.equals("measureKNN"))
      return m_kNN;
    else
      return m_NNSearch.getMeasure(additionalMeasureName);
  }
  
  
  /**
   * Returns a description of this classifier.
   *
   * @return a description of this classifier as a string.
   */
  public String toString() {

    if (m_Train == null) {
      return "IBk: No model built yet.";
    }
    
    if (m_Train.numInstances() == 0) {
      return "Warning: no training instances - ZeroR model used.";
    }

    if (!m_kNNValid && m_CrossValidate) {
      crossValidate();
    }
    
    String result = "IB1 instance-based classifier\n" +
      "using " + m_kNN;

    switch (m_DistanceWeighting) {
    case WEIGHT_INVERSE:
      result += " inverse-distance-weighted";
      break;
    case WEIGHT_SIMILARITY:
      result += " similarity-weighted";
      break;
    }
    result += " nearest neighbour(s) for classification\n";

    if (m_WindowSize != 0) {
      result += "using a maximum of " 
	+ m_WindowSize + " (windowed) training instances\n";
    }
    return result;
  }

  /**
   * Initialise scheme variables.
   */
  protected void init() {

    setKNN(1);
    m_WindowSize = 0;
    m_DistanceWeighting = WEIGHT_NONE;
    m_CrossValidate = false;
    m_MeanSquared = false;
  }
  
  /**
   * Turn the list of nearest neighbors into a probability distribution.
   *
   * @param neighbours the list of nearest neighboring instances
   * @param distances the distances of the neighbors
   * @return the probability distribution
   * @throws Exception if computation goes wrong or has no class attribute
   */
  protected double [] makeDistribution(Instances neighbours, double[] distances)
    throws Exception {

    double total = 0, weight;
    double [] distribution = new double [m_NumClasses];
    
    // Set up a correction to the estimator
    if (m_ClassType == Attribute.NOMINAL || m_ClassType == PreferenceAttribute.RANKING) {
      for(int i = 0; i < m_NumClasses; i++) {
	distribution[i] = 1.0 / Math.max(1,m_Train.numInstances());
      }
      total = (double)m_NumClasses / Math.max(1,m_Train.numInstances());
    }

    for(int i=0; i < neighbours.numInstances(); i++) {
      // Collect class counts
      Instance current = neighbours.instance(i);
      distances[i] = distances[i]*distances[i];
      distances[i] = Math.sqrt(distances[i]/m_NumAttributesUsed);
      switch (m_DistanceWeighting) {
        case WEIGHT_INVERSE:
          weight = 1.0 / (distances[i] + 0.001); // to avoid div by zero
          break;
        case WEIGHT_SIMILARITY:
          weight = 1.0 - distances[i];
          break;
        default:                                 // WEIGHT_NONE:
          weight = 1.0;
          break;
      }
      weight *= current.weight();
      try {
        switch (m_ClassType) {
          case Attribute.NOMINAL:
          case PreferenceAttribute.RANKING:
            distribution[(int)current.classValue()] += weight;
            break;
          case Attribute.NUMERIC:
            distribution[0] += current.classValue() * weight;
            break;
        }
      } catch (Exception ex) {
        throw new Error("Data has no class attribute!");
      }
      total += weight;      
    }

    // Normalise distribution
    if (total > 0) {
      Utils.normalize(distribution, total);
    }
    return distribution;
  }

  /**
   * Select the best value for k by hold-one-out cross-validation.
   * If the class attribute is nominal, classification error is
   * minimised. If the class attribute is numeric, mean absolute
   * error is minimised
   */
  protected void crossValidate() {

    try {
      if (m_NNSearch instanceof weka.core.neighboursearch.CoverTree)
	throw new Exception("CoverTree doesn't support hold-one-out "+
			    "cross-validation. Use some other NN " +
			    "method.");

      double [] performanceStats = new double [m_kNNUpper];
      double [] performanceStatsSq = new double [m_kNNUpper];

      for(int i = 0; i < m_kNNUpper; i++) {
	performanceStats[i] = 0;
	performanceStatsSq[i] = 0;
      }


      m_kNN = m_kNNUpper;
      Instance instance;
      Instances neighbours;
      double[] origDistances, convertedDistances;
      for(int i = 0; i < m_Train.numInstances(); i++) {
	if (m_Debug && (i % 50 == 0)) {
	  System.err.print("Cross validating "
			   + i + "/" + m_Train.numInstances() + "\r");
	}
	instance = m_Train.instance(i);
	neighbours = m_NNSearch.kNearestNeighbours(instance, m_kNN);
        origDistances = m_NNSearch.getDistances();
        
	for(int j = m_kNNUpper - 1; j >= 0; j--) {
	  // Update the performance stats
          convertedDistances = new double[origDistances.length];
          System.arraycopy(origDistances, 0, 
                           convertedDistances, 0, origDistances.length);
	  double [] distribution = makeDistribution(neighbours, 
                                                    convertedDistances);
          double thisPrediction = Utils.maxIndex(distribution);
	  if (m_Train.classAttribute().isNumeric()) {
	    thisPrediction = distribution[0];
	    double err = thisPrediction - instance.classValue();
	    performanceStatsSq[j] += err * err;   // Squared error
	    performanceStats[j] += Math.abs(err); // Absolute error
	  } else {
	    if (thisPrediction != instance.classValue()) {
	      performanceStats[j] ++;             // Classification error
	    }
	  }
	  if (j >= 1) {
	    neighbours = pruneToK(neighbours, convertedDistances, j);
	  }
	}
      }

      // Display the results of the cross-validation
      for(int i = 0; i < m_kNNUpper; i++) {
	if (m_Debug) {
	  System.err.print("Hold-one-out performance of " + (i + 1)
			   + " neighbors " );
	}
	if (m_Train.classAttribute().isNumeric()) {
	  if (m_Debug) {
	    if (m_MeanSquared) {
	      System.err.println("(RMSE) = "
				 + Math.sqrt(performanceStatsSq[i]
					     / m_Train.numInstances()));
	    } else {
	      System.err.println("(MAE) = "
				 + performanceStats[i]
				 / m_Train.numInstances());
	    }
	  }
	} else {
	  if (m_Debug) {
	    System.err.println("(%ERR) = "
			       + 100.0 * performanceStats[i]
			       / m_Train.numInstances());
	  }
	}
      }


      // Check through the performance stats and select the best
      // k value (or the lowest k if more than one best)
      double [] searchStats = performanceStats;
      if (m_Train.classAttribute().isNumeric() && m_MeanSquared) {
	searchStats = performanceStatsSq;
      }
      double bestPerformance = Double.NaN;
      int bestK = 1;
      for(int i = 0; i < m_kNNUpper; i++) {
	if (Double.isNaN(bestPerformance)
	    || (bestPerformance > searchStats[i])) {
	  bestPerformance = searchStats[i];
	  bestK = i + 1;
	}
      }
      m_kNN = bestK;
      if (m_Debug) {
	System.err.println("Selected k = " + bestK);
      }
      
      m_kNNValid = true;
    } catch (Exception ex) {
      throw new Error("Couldn't optimize by cross-validation: "
		      +ex.getMessage());
    }
  }
  
  /**
   * Prunes the list to contain the k nearest neighbors. If there are
   * multiple neighbors at the k'th distance, all will be kept.
   *
   * @param neighbours the neighbour instances.
   * @param distances the distances of the neighbours from target instance.
   * @param k the number of neighbors to keep.
   * @return the pruned neighbours.
   */
  public Instances pruneToK(Instances neighbours, double[] distances, int k) {
    
    if(neighbours==null || distances==null || neighbours.numInstances()==0) {
      return null;
    }
    if (k < 1) {
      k = 1;
    }
    
    int currentK = 0;
    double currentDist;
    for(int i=0; i < neighbours.numInstances(); i++) {
      currentK++;
      currentDist = distances[i];
      if(currentK>k && currentDist!=distances[i-1]) {
        currentK--;
        neighbours = new Instances(neighbours, 0, currentK);
        break;
      }
    }

    return neighbours;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6572 $");
  }
  
  /**
   * Main method for testing this class.
   *
   * @param argv should contain command line options (see setOptions)
   */
  public static void main(String [] argv) {
    runClassifier(new IBk(), argv);
  }
}
