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
 *    LWL.java
 *    Copyright (C) 1999, 2002, 2003 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.lazy;

import weka.classifiers.Classifier;
import weka.classifiers.AbstractClassifier;
import weka.classifiers.SingleClassifierEnhancer;
import weka.classifiers.UpdateableClassifier;
import weka.core.Capabilities;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.neighboursearch.LinearNNSearch;
import weka.core.neighboursearch.NearestNeighbourSearch;
import weka.core.Option;
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
 * Locally weighted learning. Uses an instance-based algorithm to assign instance weights which are then used by a specified WeightedInstancesHandler.<br/>
 * Can do classification (e.g. using naive Bayes) or regression (e.g. using linear regression).<br/>
 * <br/>
 * For more info, see<br/>
 * <br/>
 * Eibe Frank, Mark Hall, Bernhard Pfahringer: Locally Weighted Naive Bayes. In: 19th Conference in Uncertainty in Artificial Intelligence, 249-256, 2003.<br/>
 * <br/>
 * C. Atkeson, A. Moore, S. Schaal (1996). Locally weighted learning. AI Review..
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;inproceedings{Frank2003,
 *    author = {Eibe Frank and Mark Hall and Bernhard Pfahringer},
 *    booktitle = {19th Conference in Uncertainty in Artificial Intelligence},
 *    pages = {249-256},
 *    publisher = {Morgan Kaufmann},
 *    title = {Locally Weighted Naive Bayes},
 *    year = {2003}
 * }
 * 
 * &#64;article{Atkeson1996,
 *    author = {C. Atkeson and A. Moore and S. Schaal},
 *    journal = {AI Review},
 *    title = {Locally weighted learning},
 *    year = {1996}
 * }
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -A
 *  The nearest neighbour search algorithm to use (default: weka.core.neighboursearch.LinearNNSearch).
 * </pre>
 * 
 * <pre> -K &lt;number of neighbours&gt;
 *  Set the number of neighbours used to set the kernel bandwidth.
 *  (default all)</pre>
 * 
 * <pre> -U &lt;number of weighting method&gt;
 *  Set the weighting kernel shape to use. 0=Linear, 1=Epanechnikov,
 *  2=Tricube, 3=Inverse, 4=Gaussian.
 *  (default 0 = Linear)</pre>
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
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @author Ashraf M. Kibriya (amk14[at-the-rate]cs[dot]waikato[dot]ac[dot]nz)
 * @version $Revision: 6055 $ 
 */
public class LWL 
  extends SingleClassifierEnhancer
  implements UpdateableClassifier, WeightedInstancesHandler, 
             TechnicalInformationHandler {

  /** for serialization. */
  static final long serialVersionUID = 1979797405383665815L;

  /** The training instances used for classification. */
  protected Instances m_Train;
    
  /** The number of neighbours used to select the kernel bandwidth. */
  protected int m_kNN = -1;

  /** The weighting kernel method currently selected. */
  protected int m_WeightKernel = LINEAR;

  /** True if m_kNN should be set to all instances. */
  protected boolean m_UseAllK = true;
  
  /** The nearest neighbour search algorithm to use. 
   * (Default: weka.core.neighboursearch.LinearNNSearch) 
   */
  protected NearestNeighbourSearch m_NNSearch =  new LinearNNSearch();
  
  /** The available kernel weighting methods. */
  public static final int LINEAR       = 0;
  public static final int EPANECHNIKOV = 1;
  public static final int TRICUBE      = 2;  
  public static final int INVERSE      = 3;
  public static final int GAUSS        = 4;
  public static final int CONSTANT     = 5;

  /** a ZeroR model in case no model can be built from the data. */
  protected Classifier m_ZeroR;
    
  /**
   * Returns a string describing classifier.
   * @return a description suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "Locally weighted learning. Uses an instance-based algorithm to "
      + "assign instance weights which are then used by a specified "
      + "WeightedInstancesHandler.\n"
      + "Can do classification (e.g. using naive Bayes) or regression "
      + "(e.g. using linear regression).\n\n"
      + "For more info, see\n\n"
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
    TechnicalInformation 	additional;
    
    result = new TechnicalInformation(Type.INPROCEEDINGS);
    result.setValue(Field.AUTHOR, "Eibe Frank and Mark Hall and Bernhard Pfahringer");
    result.setValue(Field.YEAR, "2003");
    result.setValue(Field.TITLE, "Locally Weighted Naive Bayes");
    result.setValue(Field.BOOKTITLE, "19th Conference in Uncertainty in Artificial Intelligence");
    result.setValue(Field.PAGES, "249-256");
    result.setValue(Field.PUBLISHER, "Morgan Kaufmann");
    
    additional = result.add(Type.ARTICLE);
    additional.setValue(Field.AUTHOR, "C. Atkeson and A. Moore and S. Schaal");
    additional.setValue(Field.YEAR, "1996");
    additional.setValue(Field.TITLE, "Locally weighted learning");
    additional.setValue(Field.JOURNAL, "AI Review");
    
    return result;
  }
    
  /**
   * Constructor.
   */
  public LWL() {    
    m_Classifier = new weka.classifiers.trees.DecisionStump();
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
   * Returns an enumeration of the additional measure names 
   * produced by the neighbour search algorithm.
   * @return an enumeration of the measure names
   */
  public Enumeration enumerateMeasures() {
    return m_NNSearch.enumerateMeasures();
  }
  
  /**
   * Returns the value of the named measure from the 
   * neighbour search algorithm.
   * @param additionalMeasureName the name of the measure to query for its value
   * @return the value of the named measure
   * @throws IllegalArgumentException if the named measure is not supported
   */
  public double getMeasure(String additionalMeasureName) {
    return m_NNSearch.getMeasure(additionalMeasureName);
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    
    Vector newVector = new Vector(3);
    newVector.addElement(new Option("\tThe nearest neighbour search " +
                                    "algorithm to use " +
                                    "(default: weka.core.neighboursearch.LinearNNSearch).\n",
                                    "A", 0, "-A"));
    newVector.addElement(new Option("\tSet the number of neighbours used to set"
				    +" the kernel bandwidth.\n"
				    +"\t(default all)",
				    "K", 1, "-K <number of neighbours>"));
    newVector.addElement(new Option("\tSet the weighting kernel shape to use."
				    +" 0=Linear, 1=Epanechnikov,\n"
				    +"\t2=Tricube, 3=Inverse, 4=Gaussian.\n"
				    +"\t(default 0 = Linear)",
				    "U", 1,"-U <number of weighting method>"));
    
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
   * <pre> -A
   *  The nearest neighbour search algorithm to use (default: weka.core.neighboursearch.LinearNNSearch).
   * </pre>
   * 
   * <pre> -K &lt;number of neighbours&gt;
   *  Set the number of neighbours used to set the kernel bandwidth.
   *  (default all)</pre>
   * 
   * <pre> -U &lt;number of weighting method&gt;
   *  Set the weighting kernel shape to use. 0=Linear, 1=Epanechnikov,
   *  2=Tricube, 3=Inverse, 4=Gaussian.
   *  (default 0 = Linear)</pre>
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

    String knnString = Utils.getOption('K', options);
    if (knnString.length() != 0) {
      setKNN(Integer.parseInt(knnString));
    } else {
      setKNN(-1);
    }

    String weightString = Utils.getOption('U', options);
    if (weightString.length() != 0) {
      setWeightingKernel(Integer.parseInt(weightString));
    } else {
      setWeightingKernel(LINEAR);
    }
    
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

    super.setOptions(options);
  }

  /**
   * Gets the current settings of the classifier.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String [] getOptions() {

    String [] superOptions = super.getOptions();
    String [] options = new String [superOptions.length + 6];

    int current = 0;

    options[current++] = "-U"; options[current++] = "" + getWeightingKernel();
    if ( (getKNN() == 0) && m_UseAllK) {
      options[current++] = "-K"; options[current++] = "-1";
    }
    else {
      options[current++] = "-K"; options[current++] = "" + getKNN();
    }
    options[current++] = "-A";
    options[current++] = m_NNSearch.getClass().getName()+" "+Utils.joinOptions(m_NNSearch.getOptions()); 

    System.arraycopy(superOptions, 0, options, current,
                     superOptions.length);

    return options;
  }
  
  /**
   * Returns the tip text for this property.
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String KNNTipText() {
    return "How many neighbours are used to determine the width of the "
      + "weighting function (<= 0 means all neighbours).";
  }

  /**
   * Sets the number of neighbours used for kernel bandwidth setting.
   * The bandwidth is taken as the distance to the kth neighbour.
   *
   * @param knn the number of neighbours included inside the kernel
   * bandwidth, or 0 to specify using all neighbors.
   */
  public void setKNN(int knn) {

    m_kNN = knn;
    if (knn <= 0) {
      m_kNN = 0;
      m_UseAllK = true;
    } else {
      m_UseAllK = false;
    }
  }

  /**
   * Gets the number of neighbours used for kernel bandwidth setting.
   * The bandwidth is taken as the distance to the kth neighbour.
   *
   * @return the number of neighbours included inside the kernel
   * bandwidth, or 0 for all neighbours
   */
  public int getKNN() {

    return m_kNN;
  }

  /**
   * Returns the tip text for this property.
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String weightingKernelTipText() {
    return "Determines weighting function. [0 = Linear, 1 = Epnechnikov,"+
	   "2 = Tricube, 3 = Inverse, 4 = Gaussian and 5 = Constant. "+
	   "(default 0 = Linear)].";
  }

  /**
   * Sets the kernel weighting method to use. Must be one of LINEAR, 
   * EPANECHNIKOV,  TRICUBE, INVERSE, GAUSS or CONSTANT, other values
   * are ignored.
   *
   * @param kernel the new kernel method to use. Must be one of LINEAR,
   * EPANECHNIKOV,  TRICUBE, INVERSE, GAUSS or CONSTANT.
   */
  public void setWeightingKernel(int kernel) {

    if ((kernel != LINEAR)
	&& (kernel != EPANECHNIKOV)
	&& (kernel != TRICUBE)
	&& (kernel != INVERSE)
	&& (kernel != GAUSS)
	&& (kernel != CONSTANT)) {
      return;
    }
    m_WeightKernel = kernel;
  }

  /**
   * Gets the kernel weighting method to use.
   *
   * @return the new kernel method to use. Will be one of LINEAR,
   * EPANECHNIKOV,  TRICUBE, INVERSE, GAUSS or CONSTANT.
   */
  public int getWeightingKernel() {

    return m_WeightKernel;
  }

  /**
   * Returns the tip text for this property.
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String nearestNeighbourSearchAlgorithmTipText() {
    return "The nearest neighbour search algorithm to use (Default: LinearNN).";
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
   * Returns default capabilities of the classifier.
   *
   * @return      the capabilities of this classifier
   */
  public Capabilities getCapabilities() {
    Capabilities      result;
    
    
    
    if (m_Classifier != null) {
      result = m_Classifier.getCapabilities();
    } else {
      result = super.getCapabilities();
    }
    
    
    result.setMinimumNumberInstances(0);
    
    // set dependencies
    for (Capability cap: Capability.values())
      result.enableDependency(cap);
    
    result.enable(Capability.PREFERENCE_ATTRIBUTE);
    result.enable(Capability.RANKING);
    
    return result;
  }
  
  /**
   * Generates the classifier.
   *
   * @param instances set of instances serving as training data 
   * @throws Exception if the classifier has not been generated successfully
   */
  public void buildClassifier(Instances instances) throws Exception {

    if (!(m_Classifier instanceof WeightedInstancesHandler)) {
      throw new IllegalArgumentException("Classifier must be a "
					 + "WeightedInstancesHandler!");
    }

    // can classifier handle the data?
    getCapabilities().testWithFail(instances);

    // remove instances with missing class
    instances = new Instances(instances);
    instances.deleteWithMissingClass();
    
    // only class? -> build ZeroR model
    if (instances.numAttributes() == 1) {
      System.err.println(
	  "Cannot build model (only class attribute present in data!), "
	  + "using ZeroR model instead!");
      m_ZeroR = new weka.classifiers.rules.ZeroR();
      m_ZeroR.buildClassifier(instances);
      return;
    }
    else {
      m_ZeroR = null;
    }
    
    m_Train = new Instances(instances, 0, instances.numInstances());

    m_NNSearch.setInstances(m_Train);
  }

  /**
   * Adds the supplied instance to the training set.
   *
   * @param instance the instance to add
   * @throws Exception if instance could not be incorporated
   * successfully
   */
  public void updateClassifier(Instance instance) throws Exception {

    if (m_Train == null) {
      throw new Exception("No training instance structure set!");
    }
    else if (m_Train.equalHeaders(instance.dataset()) == false) {
      throw new Exception("Incompatible instance types\n" + m_Train.equalHeadersMsg(instance.dataset()));
    }
    if (!instance.classIsMissing()) {
      m_NNSearch.update(instance);
      m_Train.add(instance);
    }
  }
  
  /**
   * Calculates the class membership probabilities for the given test instance.
   *
   * @param instance the instance to be classified
   * @return preedicted class probability distribution
   * @throws Exception if distribution can't be computed successfully
   */
  public double[] distributionForInstance(Instance instance) throws Exception {
    
    // default model?
    if (m_ZeroR != null) {
      return m_ZeroR.distributionForInstance(instance);
    }
    
    if (m_Train.numInstances() == 0) {
      throw new Exception("No training instances!");
    }
    
    m_NNSearch.addInstanceInfo(instance);
    
    int k = m_Train.numInstances();
    if( (!m_UseAllK && (m_kNN < k)) /*&&
       !(m_WeightKernel==INVERSE ||
         m_WeightKernel==GAUSS)*/ ) {
      k = m_kNN;
    }
    
    Instances neighbours = m_NNSearch.kNearestNeighbours(instance, k);
    double distances[] = m_NNSearch.getDistances();

    if (m_Debug) {
      System.out.println("Test Instance: "+instance);
      System.out.println("For "+k+" kept " + neighbours.numInstances() + " out of " + 
                         m_Train.numInstances() + " instances.");
    }
    
    //IF LinearNN has skipped so much that <k neighbours are remaining.
    if(k>distances.length)
      k = distances.length;

    if (m_Debug) {
      System.out.println("Instance Distances");
      for (int i = 0; i < distances.length; i++) {
	System.out.println("" + distances[i]);
      }
    }

    // Determine the bandwidth
    double bandwidth = distances[k-1];

    // Check for bandwidth zero
    if (bandwidth <= 0) {
      //if the kth distance is zero than give all instances the same weight
      for(int i=0; i < distances.length; i++)
        distances[i] = 1;
    } else {
      // Rescale the distances by the bandwidth
      for (int i = 0; i < distances.length; i++)
        distances[i] = distances[i] / bandwidth;
    }
    
    // Pass the distances through a weighting kernel
    for (int i = 0; i < distances.length; i++) {
      switch (m_WeightKernel) {
        case LINEAR:
          distances[i] = 1.0001 - distances[i];
          break;
        case EPANECHNIKOV:
          distances[i] = 3/4D*(1.0001 - distances[i]*distances[i]);
          break;
        case TRICUBE:
          distances[i] = Math.pow( (1.0001 - Math.pow(distances[i], 3)), 3 );
          break;
        case CONSTANT:
          //System.err.println("using constant kernel");
          distances[i] = 1;
          break;
        case INVERSE:
          distances[i] = 1.0 / (1.0 + distances[i]);
          break;
        case GAUSS:
          distances[i] = Math.exp(-distances[i] * distances[i]);
          break;
      }
    }

    if (m_Debug) {
      System.out.println("Instance Weights");
      for (int i = 0; i < distances.length; i++) {
	System.out.println("" + distances[i]);
      }
    }
    
    // Set the weights on the training data
    double sumOfWeights = 0, newSumOfWeights = 0;
    for (int i = 0; i < distances.length; i++) {
      double weight = distances[i];
      Instance inst = (Instance) neighbours.instance(i);
      sumOfWeights += inst.weight();
      newSumOfWeights += inst.weight() * weight;
      inst.setWeight(inst.weight() * weight);
      //weightedTrain.add(newInst);
    }
    
    // Rescale weights
    for (int i = 0; i < neighbours.numInstances(); i++) {
      Instance inst = neighbours.instance(i);
      inst.setWeight(inst.weight() * sumOfWeights / newSumOfWeights);
    }

    // Create a weighted classifier
    m_Classifier.buildClassifier(neighbours);

    if (m_Debug) {
      System.out.println("Classifying test instance: " + instance);
      System.out.println("Built base classifier:\n" 
			 + m_Classifier.toString());
    }

    // Return the classifier's predictions
    return m_Classifier.distributionForInstance(instance);
  }
 
  /**
   * Returns a description of this classifier.
   *
   * @return a description of this classifier as a string.
   */
  public String toString() {

    // only ZeroR model?
    if (m_ZeroR != null) {
      StringBuffer buf = new StringBuffer();
      buf.append(this.getClass().getName().replaceAll(".*\\.", "") + "\n");
      buf.append(this.getClass().getName().replaceAll(".*\\.", "").replaceAll(".", "=") + "\n\n");
      buf.append("Warning: No model could be built, hence ZeroR model is used:\n\n");
      buf.append(m_ZeroR.toString());
      return buf.toString();
    }
    
    if (m_Train == null) {
      return "Locally weighted learning: No model built yet.";
    }
    String result = "Locally weighted learning\n"
      + "===========================\n";

    result += "Using classifier: " + m_Classifier.getClass().getName() + "\n";

    switch (m_WeightKernel) {
    case LINEAR:
      result += "Using linear weighting kernels\n";
      break;
    case EPANECHNIKOV:
      result += "Using epanechnikov weighting kernels\n";
      break;
    case TRICUBE:
      result += "Using tricube weighting kernels\n";
      break;
    case INVERSE:
      result += "Using inverse-distance weighting kernels\n";
      break;
    case GAUSS:
      result += "Using gaussian weighting kernels\n";
      break;
    case CONSTANT:
      result += "Using constant weighting kernels\n";
      break;
    }
    result += "Using " + (m_UseAllK ? "all" : "" + m_kNN) + " neighbours";
    return result;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6055 $");
  }
  
  /**
   * Main method for testing this class.
   *
   * @param argv the options
   */
  public static void main(String [] argv) {
    runClassifier(new LWL(), argv);
  }
}
