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
 *    WrapperSubsetEval.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.attributeSelection;

import weka.classifiers.Classifier;
import weka.classifiers.AbstractClassifier;
import weka.classifiers.Evaluation;
import weka.classifiers.rules.ZeroR;
import weka.core.Capabilities;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.SelectedTag;
import weka.core.Tag;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformationHandler;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;
import weka.filters.Filter;
import weka.filters.unsupervised.attribute.Remove;

import java.util.BitSet;
import java.util.Enumeration;
import java.util.Random;
import java.util.Vector;

/** 
 <!-- globalinfo-start -->
 * WrapperSubsetEval:<br/>
 * <br/>
 * Evaluates attribute sets by using a learning scheme. Cross validation is used to estimate the accuracy of the learning scheme for a set of attributes.<br/>
 * <br/>
 * For more information see:<br/>
 * <br/>
 * Ron Kohavi, George H. John (1997). Wrappers for feature subset selection. Artificial Intelligence. 97(1-2):273-324.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;article{Kohavi1997,
 *    author = {Ron Kohavi and George H. John},
 *    journal = {Artificial Intelligence},
 *    note = {Special issue on relevance},
 *    number = {1-2},
 *    pages = {273-324},
 *    title = {Wrappers for feature subset selection},
 *    volume = {97},
 *    year = {1997},
 *    ISSN = {0004-3702}
 * }
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -B &lt;base learner&gt;
 *  class name of base learner to use for  accuracy estimation.
 *  Place any classifier options LAST on the command line
 *  following a "--". eg.:
 *   -B weka.classifiers.bayes.NaiveBayes ... -- -K
 *  (default: weka.classifiers.rules.ZeroR)</pre>
 * 
 * <pre> -F &lt;num&gt;
 *  number of cross validation folds to use for estimating accuracy.
 *  (default=5)</pre>
 * 
 * <pre> -R &lt;seed&gt;
 *  Seed for cross validation accuracy testimation.
 *  (default = 1)</pre>
 * 
 * <pre> -T &lt;num&gt;
 *  threshold by which to execute another cross validation
 *  (standard deviation---expressed as a percentage of the mean).
 *  (default: 0.01 (1%))</pre>
 * 
 * <pre> -E &lt;acc | rmse | mae | f-meas | auc&gt;
 *  Performance evaluation measure to use for selecting attributes.
 *  (Default = accuracy for discrete class and rmse for numeric class)</pre>
 * 
 * <pre> 
 * Options specific to scheme weka.classifiers.rules.ZeroR:
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
public class WrapperSubsetEval
  extends ASEvaluation
  implements SubsetEvaluator,
             OptionHandler, 
             TechnicalInformationHandler {
  
  /** for serialization */
  static final long serialVersionUID = -4573057658746728675L;

  /** training instances */
  private Instances m_trainInstances;
  /** class index */
  private int m_classIndex;
  /** number of attributes in the training data */
  private int m_numAttribs;
  /** number of instances in the training data */
  private int m_numInstances;
  /** holds an evaluation object */
  private Evaluation m_Evaluation;
  /** holds the base classifier object */
  private Classifier m_BaseClassifier;
  /** number of folds to use for cross validation */
  private int m_folds;
  /** random number seed */
  private int m_seed;
  /** 
   * the threshold by which to do further cross validations when
   * estimating the accuracy of a subset
   */
  private double m_threshold;
  
  public static final int EVAL_DEFAULT = 1;
  public static final int EVAL_ACCURACY = 2;
  public static final int EVAL_RMSE = 3;
  public static final int EVAL_MAE = 4;
  public static final int EVAL_FMEASURE = 5;
  public static final int EVAL_AUC = 6;
  
  public static final Tag[] TAGS_EVALUATION = {
    new Tag(EVAL_DEFAULT, "Default: accuracy (discrete class); RMSE (numeric class)"),
    new Tag(EVAL_ACCURACY, "Accuracy (discrete class only)"),
    new Tag(EVAL_RMSE, "RMSE (of the class probabilities for discrete class)"),
    new Tag(EVAL_MAE, "MAE (of the class probabilities for discrete class)"),
    new Tag(EVAL_FMEASURE, "F-measure (discrete class only)"),
    new Tag(EVAL_AUC, "AUC (area under the ROC curve - discrete class only)")
  };
  
  /** The evaluation measure to use */
  protected int m_evaluationMeasure = EVAL_DEFAULT;

  /**
   * Returns a string describing this attribute evaluator
   * @return a description of the evaluator suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return "WrapperSubsetEval:\n\n"
      +"Evaluates attribute sets by using a learning scheme. Cross "
      +"validation is used to estimate the accuracy of the learning "
      +"scheme for a set of attributes.\n\n"
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
    TechnicalInformation 	result;
    
    result = new TechnicalInformation(Type.ARTICLE);
    result.setValue(Field.AUTHOR, "Ron Kohavi and George H. John");
    result.setValue(Field.YEAR, "1997");
    result.setValue(Field.TITLE, "Wrappers for feature subset selection");
    result.setValue(Field.JOURNAL, "Artificial Intelligence");
    result.setValue(Field.VOLUME, "97");
    result.setValue(Field.NUMBER, "1-2");
    result.setValue(Field.PAGES, "273-324");
    result.setValue(Field.NOTE, "Special issue on relevance");
    result.setValue(Field.ISSN, "0004-3702");
    
    return result;
  }

  /**
   * Constructor. Calls restOptions to set default options
   **/
  public WrapperSubsetEval () {
    resetOptions();
  }


  /**
   * Returns an enumeration describing the available options.
   * @return an enumeration of all the available options.
   **/
  public Enumeration listOptions () {
    Vector newVector = new Vector(4);
    newVector.addElement(new Option(
	"\tclass name of base learner to use for \taccuracy estimation.\n"
	+ "\tPlace any classifier options LAST on the command line\n"
	+ "\tfollowing a \"--\". eg.:\n"
	+ "\t\t-B weka.classifiers.bayes.NaiveBayes ... -- -K\n"
	+ "\t(default: weka.classifiers.rules.ZeroR)", 
	"B", 1, "-B <base learner>"));
    
    newVector.addElement(new Option(
	"\tnumber of cross validation folds to use for estimating accuracy.\n" 
	+ "\t(default=5)", 
	"F", 1, "-F <num>"));
    
    newVector.addElement(new Option(
	"\tSeed for cross validation accuracy testimation.\n"
	+ "\t(default = 1)", 
	"R", 1,"-R <seed>"));
    
    newVector.addElement(new Option(
	"\tthreshold by which to execute another cross validation\n" 
	+ "\t(standard deviation---expressed as a percentage of the mean).\n"
	+ "\t(default: 0.01 (1%))", 
	"T", 1, "-T <num>"));
    
    newVector.addElement(new Option(
        "\tPerformance evaluation measure to use for selecting attributes.\n" +
        "\t(Default = accuracy for discrete class and rmse for numeric class)",
        "E", 1, "-E <acc | rmse | mae | f-meas | auc>"));

    if ((m_BaseClassifier != null) && 
	(m_BaseClassifier instanceof OptionHandler)) {
      newVector.addElement(new Option("", "", 0, "\nOptions specific to scheme " 
				      + m_BaseClassifier.getClass().getName() 
				      + ":"));
      Enumeration enu = ((OptionHandler)m_BaseClassifier).listOptions();

      while (enu.hasMoreElements()) {
        newVector.addElement(enu.nextElement());
      }
    }

    return  newVector.elements();
  }


  /**
   * Parses a given list of options. <p/>
   *
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -B &lt;base learner&gt;
   *  class name of base learner to use for  accuracy estimation.
   *  Place any classifier options LAST on the command line
   *  following a "--". eg.:
   *   -B weka.classifiers.bayes.NaiveBayes ... -- -K
   *  (default: weka.classifiers.rules.ZeroR)</pre>
   * 
   * <pre> -F &lt;num&gt;
   *  number of cross validation folds to use for estimating accuracy.
   *  (default=5)</pre>
   * 
   * <pre> -R &lt;seed&gt;
   *  Seed for cross validation accuracy testimation.
   *  (default = 1)</pre>
   * 
   * <pre> -T &lt;num&gt;
   *  threshold by which to execute another cross validation
   *  (standard deviation---expressed as a percentage of the mean).
   *  (default: 0.01 (1%))</pre>
   * 
   * <pre> -E &lt;acc | rmse | mae | f-meas | auc&gt;
   *  Performance evaluation measure to use for selecting attributes.
   *  (Default = accuracy for discrete class and rmse for numeric class)</pre>
   * 
   * <pre> 
   * Options specific to scheme weka.classifiers.rules.ZeroR:
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
  public void setOptions (String[] options)
    throws Exception {
    String optionString;
    resetOptions();
    optionString = Utils.getOption('B', options);

    if (optionString.length() == 0)
      optionString = ZeroR.class.getName();
    setClassifier(AbstractClassifier.forName(optionString, 
				     Utils.partitionOptions(options)));
    optionString = Utils.getOption('F', options);

    if (optionString.length() != 0) {
      setFolds(Integer.parseInt(optionString));
    }

    optionString = Utils.getOption('R', options);
    if (optionString.length() != 0) {
      setSeed(Integer.parseInt(optionString));
    }

    //       optionString = Utils.getOption('S',options);
    //       if (optionString.length() != 0)
    //         {
    //  	 seed = Integer.parseInt(optionString);
    //         }
    optionString = Utils.getOption('T', options);

    if (optionString.length() != 0) {
      Double temp;
      temp = Double.valueOf(optionString);
      setThreshold(temp.doubleValue());
    }
    
    optionString = Utils.getOption('E', options);
    if (optionString.length() != 0) {
      if (optionString.equals("acc")) {
        setEvaluationMeasure(new SelectedTag(EVAL_ACCURACY, TAGS_EVALUATION));
      } else if (optionString.equals("rmse")) {
        setEvaluationMeasure(new SelectedTag(EVAL_RMSE, TAGS_EVALUATION));
      } else if (optionString.equals("mae")) {
        setEvaluationMeasure(new SelectedTag(EVAL_MAE, TAGS_EVALUATION));
      } else if (optionString.equals("f-meas")) {
        setEvaluationMeasure(new SelectedTag(EVAL_FMEASURE, TAGS_EVALUATION));
      } else if (optionString.equals("auc")) {
        setEvaluationMeasure(new SelectedTag(EVAL_AUC, TAGS_EVALUATION));
      } else {
        throw new IllegalArgumentException("Invalid evaluation measure");
      }
    }
  }
  
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String evaluationMeasureTipText() {
    return "The measure used to evaluate the performance of attribute combinations.";
  }
  /**
   * Gets the currently set performance evaluation measure used for selecting
   * attributes for the decision table
   * 
   * @return the performance evaluation measure
   */
  public SelectedTag getEvaluationMeasure() {
    return new SelectedTag(m_evaluationMeasure, TAGS_EVALUATION);
  }

  /**
   * Sets the performance evaluation measure to use for selecting attributes
   * for the decision table
   * 
   * @param newMethod the new performance evaluation metric to use
   */
  public void setEvaluationMeasure(SelectedTag newMethod) {
    if (newMethod.getTags() == TAGS_EVALUATION) {
      m_evaluationMeasure = newMethod.getSelectedTag().getID();
    }
  }
  
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String thresholdTipText() {
    return "Repeat xval if stdev of mean exceeds this value.";
  }

  /**
   * Set the value of the threshold for repeating cross validation
   *
   * @param t the value of the threshold
   */
  public void setThreshold (double t) {
    m_threshold = t;
  }


  /**
   * Get the value of the threshold
   *
   * @return the threshold as a double
   */
  public double getThreshold () {
    return  m_threshold;
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String foldsTipText() {
    return "Number of xval folds to use when estimating subset accuracy.";
  }

  /**
   * Set the number of folds to use for accuracy estimation
   *
   * @param f the number of folds
   */
  public void setFolds (int f) {
    m_folds = f;
  }


  /**
   * Get the number of folds used for accuracy estimation
   *
   * @return the number of folds
   */
  public int getFolds () {
    return  m_folds;
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String seedTipText() {
    return "Seed to use for randomly generating xval splits.";
  }

  /**
   * Set the seed to use for cross validation
   *
   * @param s the seed
   */
  public void setSeed (int s) {
    m_seed = s;
  }


  /**
   * Get the random number seed used for cross validation
   *
   * @return the seed
   */
  public int getSeed () {
    return  m_seed;
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String classifierTipText() {
    return "Classifier to use for estimating the accuracy of subsets";
  }

  /**
   * Set the classifier to use for accuracy estimation
   *
   * @param newClassifier the Classifier to use.
   */
  public void setClassifier (Classifier newClassifier) {
    m_BaseClassifier = newClassifier;
  }


  /**
   * Get the classifier used as the base learner.
   *
   * @return the classifier used as the classifier
   */
  public Classifier getClassifier () {
    return  m_BaseClassifier;
  }


  /**
   * Gets the current settings of WrapperSubsetEval.
   *
   * @return an array of strings suitable for passing to setOptions()
   */
  public String[] getOptions () {
    String[] classifierOptions = new String[0];

    if ((m_BaseClassifier != null) && 
	(m_BaseClassifier instanceof OptionHandler)) {
      classifierOptions = ((OptionHandler)m_BaseClassifier).getOptions();
    }

    String[] options = new String[9 + classifierOptions.length];
    int current = 0;

    if (getClassifier() != null) {
      options[current++] = "-B";
      options[current++] = getClassifier().getClass().getName();
    }

    options[current++] = "-F";
    options[current++] = "" + getFolds();
    options[current++] = "-T";
    options[current++] = "" + getThreshold();
    options[current++] = "-R";
    options[current++] = "" + getSeed();
    options[current++] = "--";
    System.arraycopy(classifierOptions, 0, options, current, 
		     classifierOptions.length);
    current += classifierOptions.length;

    while (current < options.length) {
      options[current++] = "";
    }

    return  options;
  }


  protected void resetOptions () {
    m_trainInstances = null;
    m_Evaluation = null;
    m_BaseClassifier = new ZeroR();
    m_folds = 5;
    m_seed = 1;
    m_threshold = 0.01;
  }

  /**
   * Returns the capabilities of this evaluator.
   *
   * @return            the capabilities of this evaluator
   * @see               Capabilities
   */
  public Capabilities getCapabilities() {
    Capabilities	result;
    
    if (getClassifier() == null) {
      result = super.getCapabilities();
      result.disableAll();
    } else {
      result = getClassifier().getCapabilities();
    }
    
    // set dependencies
    for (Capability cap: Capability.values())
      result.enableDependency(cap);
    
    // adjustment for class based on selected evaluation metric
    result.disable(Capability.NUMERIC_CLASS);
    result.disable(Capability.DATE_CLASS);
    if (m_evaluationMeasure != EVAL_ACCURACY && m_evaluationMeasure != EVAL_FMEASURE &&
        m_evaluationMeasure != EVAL_AUC) {
      result.enable(Capability.NUMERIC_CLASS);
      result.enable(Capability.DATE_CLASS);
    }
    
    result.setMinimumNumberInstances(getFolds());
    
    return result;
  }

  /**
   * Generates a attribute evaluator. Has to initialize all fields of the 
   * evaluator that are not being set via options.
   *
   * @param data set of instances serving as training data 
   * @throws Exception if the evaluator has not been 
   * generated successfully
   */
  public void buildEvaluator (Instances data)
    throws Exception {

    // can evaluator handle data?
    getCapabilities().testWithFail(data);

    m_trainInstances = data;
    m_classIndex = m_trainInstances.classIndex();
    m_numAttribs = m_trainInstances.numAttributes();
    m_numInstances = m_trainInstances.numInstances();
  }


  /**
   * Evaluates a subset of attributes
   *
   * @param subset a bitset representing the attribute subset to be 
   * evaluated 
   * @return the error rate
   * @throws Exception if the subset could not be evaluated
   */
  public double evaluateSubset (BitSet subset)
    throws Exception {
    double evalMetric = 0;
    double[] repError = new double[5];
    int numAttributes = 0;
    int i, j;
    Random Rnd = new Random(m_seed);
    Remove delTransform = new Remove();
    delTransform.setInvertSelection(true);
    // copy the instances
    Instances trainCopy = new Instances(m_trainInstances);

    // count attributes set in the BitSet
    for (i = 0; i < m_numAttribs; i++) {
      if (subset.get(i)) {
        numAttributes++;
      }
    }

    // set up an array of attribute indexes for the filter (+1 for the class)
    int[] featArray = new int[numAttributes + 1];

    for (i = 0, j = 0; i < m_numAttribs; i++) {
      if (subset.get(i)) {
        featArray[j++] = i;
      }
    }

    featArray[j] = m_classIndex;
    delTransform.setAttributeIndicesArray(featArray);
    delTransform.setInputFormat(trainCopy);
    trainCopy = Filter.useFilter(trainCopy, delTransform);

    // max of 5 repetitions of cross validation
    for (i = 0; i < 5; i++) {
      m_Evaluation = new Evaluation(trainCopy);
      m_Evaluation.crossValidateModel(m_BaseClassifier, trainCopy, m_folds, Rnd);
      
      switch (m_evaluationMeasure) {
      case EVAL_DEFAULT:
        repError[i] = m_Evaluation.errorRate();
        break;
      case EVAL_ACCURACY:
        repError[i] = m_Evaluation.errorRate();
        break;
      case EVAL_RMSE:
        repError[i] = m_Evaluation.rootMeanSquaredError();
        break;
      case EVAL_MAE:
        repError[i] = m_Evaluation.meanAbsoluteError();
        break;
      case EVAL_FMEASURE:
        repError[i] = m_Evaluation.weightedFMeasure();
        break;
      case EVAL_AUC:
        repError[i] = m_Evaluation.weightedAreaUnderROC();
        break;
      }

      // check on the standard deviation
      if (!repeat(repError, i + 1)) {
        i++;
        break;
      }
    }

    for (j = 0; j < i; j++) {
      evalMetric += repError[j];
    }

    evalMetric /= (double)i;
    m_Evaluation = null;
    
    switch (m_evaluationMeasure) {
    case EVAL_DEFAULT:
    case EVAL_ACCURACY:
    case EVAL_RMSE:
    case EVAL_MAE:
      evalMetric = -evalMetric; // maximize
      break;
    }
    
    return evalMetric;
  }


  /**
   * Returns a string describing the wrapper
   *
   * @return the description as a string
   */
  public String toString () {
    StringBuffer text = new StringBuffer();

    if (m_trainInstances == null) {
      text.append("\tWrapper subset evaluator has not been built yet\n");
    }
    else {
      text.append("\tWrapper Subset Evaluator\n");
      text.append("\tLearning scheme: " 
		  + getClassifier().getClass().getName() + "\n");
      text.append("\tScheme options: ");
      String[] classifierOptions = new String[0];

      if (m_BaseClassifier instanceof OptionHandler) {
        classifierOptions = ((OptionHandler)m_BaseClassifier).getOptions();

        for (int i = 0; i < classifierOptions.length; i++) {
          text.append(classifierOptions[i] + " ");
        }
      }

      text.append("\n");
      switch (m_evaluationMeasure) {
      case EVAL_DEFAULT:
      case EVAL_ACCURACY:
        if (m_trainInstances.attribute(m_classIndex).isNumeric()) {
          text.append("\tSubset evaluation: RMSE\n");
        } else {
          text.append("\tSubset evaluation: classification error\n");
        }
        break;
      case EVAL_RMSE:
        if (m_trainInstances.attribute(m_classIndex).isNumeric()) {
          text.append("\tSubset evaluation: RMSE\n");
        } else {
          text.append("\tSubset evaluation: RMSE (probability estimates)\n");
        }
        break;
      case EVAL_MAE:
        if (m_trainInstances.attribute(m_classIndex).isNumeric()) {
          text.append("\tSubset evaluation: MAE\n");
        } else {
          text.append("\tSubset evaluation: MAE (probability estimates)\n");
        }
        break;
      case EVAL_FMEASURE:
        text.append("\tSubset evaluation: F-measure\n");
        break;
      case EVAL_AUC:
        text.append("\tSubset evaluation: area under the ROC curve\n");
        break;
      }
      
      text.append("\tNumber of folds for accuracy estimation: " 
          + m_folds 
          + "\n");
    }

    return  text.toString();
  }


  /**
   * decides whether to do another repeat of cross validation. If the
   * standard deviation of the cross validations
   * is greater than threshold% of the mean (default 1%) then another 
   * repeat is done. 
   *
   * @param repError an array of cross validation results
   * @param entries the number of cross validations done so far
   * @return true if another cv is to be done
   */
  private boolean repeat (double[] repError, int entries) {
    int i;
    double mean = 0;
    double variance = 0;

    if (entries == 1) {
      return  true;
    }

    for (i = 0; i < entries; i++) {
      mean += repError[i];
    }

    mean /= (double)entries;

    for (i = 0; i < entries; i++) {
      variance += ((repError[i] - mean)*(repError[i] - mean));
    }

    variance /= (double)entries;

    if (variance > 0) {
      variance = Math.sqrt(variance);
    }

    if ((variance/mean) > m_threshold) {
      return  true;
    }

    return  false;
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
   * @param args the options
   */
  public static void main (String[] args) {
    runEvaluator(new WrapperSubsetEval(), args);
  }
}

