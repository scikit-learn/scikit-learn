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
 *    AttributeSelectedClassifier.java
 *    Copyright (C) 2000 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.meta;

import weka.attributeSelection.ASEvaluation;
import weka.attributeSelection.ASSearch;
import weka.attributeSelection.AttributeSelection;
import weka.classifiers.SingleClassifierEnhancer;
import weka.core.AdditionalMeasureProducer;
import weka.core.Capabilities;
import weka.core.Drawable;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.core.WeightedInstancesHandler;
import weka.core.Capabilities.Capability;

import java.util.Enumeration;
import java.util.Random;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * Dimensionality of training and test data is reduced by attribute selection before being passed on to a classifier.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -E &lt;attribute evaluator specification&gt;
 *  Full class name of attribute evaluator, followed
 *  by its options.
 *  eg: "weka.attributeSelection.CfsSubsetEval -L"
 *  (default weka.attributeSelection.CfsSubsetEval)</pre>
 * 
 * <pre> -S &lt;search method specification&gt;
 *  Full class name of search method, followed
 *  by its options.
 *  eg: "weka.attributeSelection.BestFirst -D 1"
 *  (default weka.attributeSelection.BestFirst)</pre>
 * 
 * <pre> -D
 *  If set, classifier is run in debug mode and
 *  may output additional info to the console</pre>
 * 
 * <pre> -W
 *  Full name of base classifier.
 *  (default: weka.classifiers.trees.J48)</pre>
 * 
 * <pre> 
 * Options specific to classifier weka.classifiers.trees.J48:
 * </pre>
 * 
 * <pre> -U
 *  Use unpruned tree.</pre>
 * 
 * <pre> -C &lt;pruning confidence&gt;
 *  Set confidence threshold for pruning.
 *  (default 0.25)</pre>
 * 
 * <pre> -M &lt;minimum number of instances&gt;
 *  Set minimum number of instances per leaf.
 *  (default 2)</pre>
 * 
 * <pre> -R
 *  Use reduced error pruning.</pre>
 * 
 * <pre> -N &lt;number of folds&gt;
 *  Set number of folds for reduced error
 *  pruning. One fold is used as pruning set.
 *  (default 3)</pre>
 * 
 * <pre> -B
 *  Use binary splits only.</pre>
 * 
 * <pre> -S
 *  Don't perform subtree raising.</pre>
 * 
 * <pre> -L
 *  Do not clean up after the tree has been built.</pre>
 * 
 * <pre> -A
 *  Laplace smoothing for predicted probabilities.</pre>
 * 
 * <pre> -Q &lt;seed&gt;
 *  Seed for random data shuffling (default 1).</pre>
 * 
 <!-- options-end -->
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @version $Revision: 1.26 $
 */
public class AttributeSelectedClassifier 
  extends SingleClassifierEnhancer
  implements OptionHandler, Drawable, AdditionalMeasureProducer,
             WeightedInstancesHandler {

  /** for serialization */
  static final long serialVersionUID = -5951805453487947577L;
  
  /** The attribute selection object */
  protected AttributeSelection m_AttributeSelection = null;

  /** The attribute evaluator to use */
  protected ASEvaluation m_Evaluator = 
    new weka.attributeSelection.CfsSubsetEval();

  /** The search method to use */
  protected ASSearch m_Search = new weka.attributeSelection.BestFirst();

  /** The header of the dimensionally reduced data */
  protected Instances m_ReducedHeader;

  /** The number of class vals in the training data (1 if class is numeric) */
  protected int m_numClasses;

  /** The number of attributes selected by the attribute selection phase */
  protected double m_numAttributesSelected;

  /** The time taken to select attributes in milliseconds */
  protected double m_selectionTime;

  /** The time taken to select attributes AND build the classifier */
  protected double m_totalTime;

  
  /**
   * String describing default classifier.
   * 
   * @return the default classifier classname
   */
  protected String defaultClassifierString() {
    
    return "weka.classifiers.trees.J48";
  }
  
  /**
   * Default constructor.
   */
  public AttributeSelectedClassifier() {
    m_Classifier = new weka.classifiers.trees.J48();
  }

  /**
   * Returns a string describing this search method
   * @return a description of the search method suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return "Dimensionality of training and test data is reduced by "
      +"attribute selection before being passed on to a classifier.";
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {
     Vector newVector = new Vector(3);
    
    newVector.addElement(new Option(
	      "\tFull class name of attribute evaluator, followed\n"
	      + "\tby its options.\n"
	      + "\teg: \"weka.attributeSelection.CfsSubsetEval -L\"\n"
	      + "\t(default weka.attributeSelection.CfsSubsetEval)",
	      "E", 1, "-E <attribute evaluator specification>"));

    newVector.addElement(new Option(
	      "\tFull class name of search method, followed\n"
	      + "\tby its options.\n"
	      + "\teg: \"weka.attributeSelection.BestFirst -D 1\"\n"
	      + "\t(default weka.attributeSelection.BestFirst)",
	      "S", 1, "-S <search method specification>"));
    
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
   * <pre> -E &lt;attribute evaluator specification&gt;
   *  Full class name of attribute evaluator, followed
   *  by its options.
   *  eg: "weka.attributeSelection.CfsSubsetEval -L"
   *  (default weka.attributeSelection.CfsSubsetEval)</pre>
   * 
   * <pre> -S &lt;search method specification&gt;
   *  Full class name of search method, followed
   *  by its options.
   *  eg: "weka.attributeSelection.BestFirst -D 1"
   *  (default weka.attributeSelection.BestFirst)</pre>
   * 
   * <pre> -D
   *  If set, classifier is run in debug mode and
   *  may output additional info to the console</pre>
   * 
   * <pre> -W
   *  Full name of base classifier.
   *  (default: weka.classifiers.trees.J48)</pre>
   * 
   * <pre> 
   * Options specific to classifier weka.classifiers.trees.J48:
   * </pre>
   * 
   * <pre> -U
   *  Use unpruned tree.</pre>
   * 
   * <pre> -C &lt;pruning confidence&gt;
   *  Set confidence threshold for pruning.
   *  (default 0.25)</pre>
   * 
   * <pre> -M &lt;minimum number of instances&gt;
   *  Set minimum number of instances per leaf.
   *  (default 2)</pre>
   * 
   * <pre> -R
   *  Use reduced error pruning.</pre>
   * 
   * <pre> -N &lt;number of folds&gt;
   *  Set number of folds for reduced error
   *  pruning. One fold is used as pruning set.
   *  (default 3)</pre>
   * 
   * <pre> -B
   *  Use binary splits only.</pre>
   * 
   * <pre> -S
   *  Don't perform subtree raising.</pre>
   * 
   * <pre> -L
   *  Do not clean up after the tree has been built.</pre>
   * 
   * <pre> -A
   *  Laplace smoothing for predicted probabilities.</pre>
   * 
   * <pre> -Q &lt;seed&gt;
   *  Seed for random data shuffling (default 1).</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {

    // same for attribute evaluator
    String evaluatorString = Utils.getOption('E', options);
    if (evaluatorString.length() == 0)
      evaluatorString = weka.attributeSelection.CfsSubsetEval.class.getName();
    String [] evaluatorSpec = Utils.splitOptions(evaluatorString);
    if (evaluatorSpec.length == 0) {
      throw new Exception("Invalid attribute evaluator specification string");
    }
    String evaluatorName = evaluatorSpec[0];
    evaluatorSpec[0] = "";
    setEvaluator(ASEvaluation.forName(evaluatorName, evaluatorSpec));

    // same for search method
    String searchString = Utils.getOption('S', options);
    if (searchString.length() == 0)
      searchString = weka.attributeSelection.BestFirst.class.getName();
    String [] searchSpec = Utils.splitOptions(searchString);
    if (searchSpec.length == 0) {
      throw new Exception("Invalid search specification string");
    }
    String searchName = searchSpec[0];
    searchSpec[0] = "";
    setSearch(ASSearch.forName(searchName, searchSpec));

    super.setOptions(options);
  }

  /**
   * Gets the current settings of the Classifier.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String [] getOptions() {

    String [] superOptions = super.getOptions();
    String [] options = new String [superOptions.length + 4];

    int current = 0;

    // same attribute evaluator
    options[current++] = "-E";
    options[current++] = "" +getEvaluatorSpec();
    
    // same for search
    options[current++] = "-S";
    options[current++] = "" + getSearchSpec();

    System.arraycopy(superOptions, 0, options, current, 
		     superOptions.length);
    
    return options;
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String evaluatorTipText() {
    return "Set the attribute evaluator to use. This evaluator is used "
      +"during the attribute selection phase before the classifier is "
      +"invoked.";
  }

  /**
   * Sets the attribute evaluator
   *
   * @param evaluator the evaluator with all options set.
   */
  public void setEvaluator(ASEvaluation evaluator) {
    m_Evaluator = evaluator;
  }

  /**
   * Gets the attribute evaluator used
   *
   * @return the attribute evaluator
   */
  public ASEvaluation getEvaluator() {
    return m_Evaluator;
  }

  /**
   * Gets the evaluator specification string, which contains the class name of
   * the attribute evaluator and any options to it
   *
   * @return the evaluator string.
   */
  protected String getEvaluatorSpec() {
    
    ASEvaluation e = getEvaluator();
    if (e instanceof OptionHandler) {
      return e.getClass().getName() + " "
	+ Utils.joinOptions(((OptionHandler)e).getOptions());
    }
    return e.getClass().getName();
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String searchTipText() {
    return "Set the search method. This search method is used "
      +"during the attribute selection phase before the classifier is "
      +"invoked.";
  }
  
  /**
   * Sets the search method
   *
   * @param search the search method with all options set.
   */
  public void setSearch(ASSearch search) {
    m_Search = search;
  }

  /**
   * Gets the search method used
   *
   * @return the search method
   */
  public ASSearch getSearch() {
    return m_Search;
  }

  /**
   * Gets the search specification string, which contains the class name of
   * the search method and any options to it
   *
   * @return the search string.
   */
  protected String getSearchSpec() {
    
    ASSearch s = getSearch();
    if (s instanceof OptionHandler) {
      return s.getClass().getName() + " "
	+ Utils.joinOptions(((OptionHandler)s).getOptions());
    }
    return s.getClass().getName();
  }

  /**
   * Returns default capabilities of the classifier.
   *
   * @return      the capabilities of this classifier
   */
  public Capabilities getCapabilities() {
    Capabilities	result;
    
    if (getEvaluator() == null)
      result = super.getCapabilities();
    else
      result = getEvaluator().getCapabilities();
    
    // set dependencies
    for (Capability cap: Capability.values())
      result.enableDependency(cap);
    
    //RANKING BEGIN
    result.disable(Capability.PREFERENCE_ATTRIBUTE);
    result.disable(Capability.RANKING);
    //RANKING END
    
    return result;
  }

  /**
   * Build the classifier on the dimensionally reduced data.
   *
   * @param data the training data
   * @throws Exception if the classifier could not be built successfully
   */
  public void buildClassifier(Instances data) throws Exception {
    if (m_Classifier == null) {
      throw new Exception("No base classifier has been set!");
    }

    if (m_Evaluator == null) {
      throw new Exception("No attribute evaluator has been set!");
    }

    if (m_Search == null) {
      throw new Exception("No search method has been set!");
    }
   
    // can classifier handle the data?
    getCapabilities().testWithFail(data);

    // remove instances with missing class
    Instances newData = new Instances(data);
    newData.deleteWithMissingClass();
    
    if (newData.numInstances() == 0) {
      m_Classifier.buildClassifier(newData);
      return;
    }
    if (newData.classAttribute().isNominal() || newData.classAttribute().isRanking()) {
      m_numClasses = newData.classAttribute().numValues();
    } else {
      m_numClasses = 1;
    }

    Instances resampledData = null;
    // check to see if training data has all equal weights
    double weight = newData.instance(0).weight();
    boolean ok = false;
    for (int i = 1; i < newData.numInstances(); i++) {
      if (newData.instance(i).weight() != weight) {
        ok = true;
        break;
      }
    }
    
    if (ok) {
      if (!(m_Evaluator instanceof WeightedInstancesHandler) || 
          !(m_Classifier instanceof WeightedInstancesHandler)) {
        Random r = new Random(1);
        for (int i = 0; i < 10; i++) {
          r.nextDouble();
        }
        resampledData = newData.resampleWithWeights(r);
      }
    } else {
      // all equal weights in the training data so just use as is
      resampledData = newData;
    }

    m_AttributeSelection = new AttributeSelection();
    m_AttributeSelection.setEvaluator(m_Evaluator);
    m_AttributeSelection.setSearch(m_Search);
    long start = System.currentTimeMillis();
    m_AttributeSelection.
      SelectAttributes((m_Evaluator instanceof WeightedInstancesHandler) 
                       ? newData
                       : resampledData);
    long end = System.currentTimeMillis();
    if (m_Classifier instanceof WeightedInstancesHandler) {
      newData = m_AttributeSelection.reduceDimensionality(newData);
      m_Classifier.buildClassifier(newData);
    } else {
      resampledData = m_AttributeSelection.reduceDimensionality(resampledData);
      m_Classifier.buildClassifier(resampledData);
    }

    long end2 = System.currentTimeMillis();
    m_numAttributesSelected = m_AttributeSelection.numberAttributesSelected();
    m_ReducedHeader = 
      new Instances((m_Classifier instanceof WeightedInstancesHandler) ?
                    newData
                    : resampledData, 0);
    m_selectionTime = (double)(end - start);
    m_totalTime = (double)(end2 - start);
  }

  /**
   * Classifies a given instance after attribute selection
   *
   * @param instance the instance to be classified
   * @return the class distribution
   * @throws Exception if instance could not be classified
   * successfully
   */
  public double [] distributionForInstance(Instance instance)
    throws Exception {

    Instance newInstance;
    if (m_AttributeSelection == null) {
      //      throw new Exception("AttributeSelectedClassifier: No model built yet!");
      newInstance = instance;
    } else {
      newInstance = m_AttributeSelection.reduceDimensionality(instance);
    }

    return m_Classifier.distributionForInstance(newInstance);
  }

  /**
   *  Returns the type of graph this classifier
   *  represents.
   *  
   *  @return the type of graph
   */   
  public int graphType() {
    
    if (m_Classifier instanceof Drawable)
      return ((Drawable)m_Classifier).graphType();
    else 
      return Drawable.NOT_DRAWABLE;
  }

  /**
   * Returns graph describing the classifier (if possible).
   *
   * @return the graph of the classifier in dotty format
   * @throws Exception if the classifier cannot be graphed
   */
  public String graph() throws Exception {
    
    if (m_Classifier instanceof Drawable)
      return ((Drawable)m_Classifier).graph();
    else throw new Exception("Classifier: " + getClassifierSpec()
			     + " cannot be graphed");
  }

  /**
   * Output a representation of this classifier
   * 
   * @return a representation of this classifier
   */
  public String toString() {
    if (m_AttributeSelection == null) {
      return "AttributeSelectedClassifier: No attribute selection possible.\n\n"
	+m_Classifier.toString();
    }

    StringBuffer result = new StringBuffer();
    result.append("AttributeSelectedClassifier:\n\n");
    result.append(m_AttributeSelection.toResultsString());
    result.append("\n\nHeader of reduced data:\n"+m_ReducedHeader.toString());
    result.append("\n\nClassifier Model\n"+m_Classifier.toString());

    return result.toString();
  }

  /**
   * Additional measure --- number of attributes selected
   * @return the number of attributes selected
   */
  public double measureNumAttributesSelected() {
    return m_numAttributesSelected;
  }

  /**
   * Additional measure --- time taken (milliseconds) to select the attributes
   * @return the time taken to select attributes
   */
  public double measureSelectionTime() {
    return m_selectionTime;
  }

  /**
   * Additional measure --- time taken (milliseconds) to select attributes
   * and build the classifier
   * @return the total time (select attributes + build classifier)
   */
  public double measureTime() {
    return m_totalTime;
  }

  /**
   * Returns an enumeration of the additional measure names
   * @return an enumeration of the measure names
   */
  public Enumeration enumerateMeasures() {
    Vector newVector = new Vector(3);
    newVector.addElement("measureNumAttributesSelected");
    newVector.addElement("measureSelectionTime");
    newVector.addElement("measureTime");
    if (m_Classifier instanceof AdditionalMeasureProducer) {
      Enumeration en = ((AdditionalMeasureProducer)m_Classifier).
	enumerateMeasures();
      while (en.hasMoreElements()) {
	String mname = (String)en.nextElement();
	newVector.addElement(mname);
      }
    }
    return newVector.elements();
  }
  
  /**
   * Returns the value of the named measure
   * @param additionalMeasureName the name of the measure to query for its value
   * @return the value of the named measure
   * @throws IllegalArgumentException if the named measure is not supported
   */
  public double getMeasure(String additionalMeasureName) {
    if (additionalMeasureName.compareToIgnoreCase("measureNumAttributesSelected") == 0) {
      return measureNumAttributesSelected();
    } else if (additionalMeasureName.compareToIgnoreCase("measureSelectionTime") == 0) {
      return measureSelectionTime();
    } else if (additionalMeasureName.compareToIgnoreCase("measureTime") == 0) {
      return measureTime();
    } else if (m_Classifier instanceof AdditionalMeasureProducer) {
      return ((AdditionalMeasureProducer)m_Classifier).
	getMeasure(additionalMeasureName);
    } else {
      throw new IllegalArgumentException(additionalMeasureName 
			  + " not supported (AttributeSelectedClassifier)");
    }
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 1.26 $");
  }

  /**
   * Main method for testing this class.
   *
   * @param argv should contain the following arguments:
   * -t training file [-T test file] [-c class index]
   */
  public static void main(String [] argv) {
    runClassifier(new AttributeSelectedClassifier(), argv);
  }
}
