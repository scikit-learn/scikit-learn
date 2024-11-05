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
 *    Evaluation.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers;

import weka.classifiers.evaluation.NominalPrediction;
import weka.classifiers.evaluation.NumericPrediction;
import weka.classifiers.evaluation.ThresholdCurve;
import weka.classifiers.evaluation.output.prediction.AbstractOutput;
import weka.classifiers.evaluation.output.prediction.PlainText;
import weka.classifiers.pmml.consumer.PMMLClassifier;
import weka.classifiers.xml.XMLClassifier;
import weka.core.Drawable;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.core.Summarizable;
import weka.core.Utils;
import weka.core.Version;
import weka.core.converters.ConverterUtils.DataSink;
import weka.core.converters.ConverterUtils.DataSource;
import weka.core.labelranking.PreferenceAttribute;
import weka.core.labelranking.PreferenceDenseInstance;
import weka.core.labelranking.RankUtilities;
import weka.core.pmml.PMMLFactory;
import weka.core.pmml.PMMLModel;
import weka.core.xml.KOML;
import weka.core.xml.XMLOptions;
import weka.core.xml.XMLSerialization;
import weka.estimators.UnivariateKernelEstimator;

import java.beans.BeanInfo;
import java.beans.Introspector;
import java.beans.MethodDescriptor;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStream;
import java.io.Reader;
import java.lang.reflect.Method;
import java.util.Date;
import java.util.Enumeration;
import java.util.Random;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 * Class for evaluating machine learning models. <p/>
 *
 * ------------------------------------------------------------------- <p/>
 *
 * General options when evaluating a learning scheme from the command-line: <p/>
 *
 * -t filename <br/>
 * Name of the file with the training data. (required) <p/>
 *
 * -T filename <br/>
 * Name of the file with the test data. If missing a cross-validation
 * is performed. <p/>
 *
 * -c index <br/>
 * Index of the class attribute (1, 2, ...; default: last). <p/>
 *
 * -x number <br/>
 * The number of folds for the cross-validation (default: 10). <p/>
 *
 * -no-cv <br/>
 * No cross validation.  If no test file is provided, no evaluation
 * is done. <p/>
 *
 * -split-percentage percentage <br/>
 * Sets the percentage for the train/test set split, e.g., 66. <p/>
 *
 * -preserve-order <br/>
 * Preserves the order in the percentage split instead of randomizing
 * the data first with the seed value ('-s'). <p/>
 *
 * -s seed <br/>
 * Random number seed for the cross-validation and percentage split
 * (default: 1). <p/>
 *
 * -m filename <br/>
 * The name of a file containing a cost matrix. <p/>
 *
 * -l filename <br/>
 * Loads classifier from the given file. In case the filename ends with ".xml",
 * a PMML file is loaded or, if that fails, options are loaded from XML. <p/>
 *
 * -d filename <br/>
 * Saves classifier built from the training data into the given file. In case
 * the filename ends with ".xml" the options are saved XML, not the model. <p/>
 *
 * -v <br/>
 * Outputs no statistics for the training data. <p/>
 *
 * -o <br/>
 * Outputs statistics only, not the classifier. <p/>
 *
 * -i <br/>
 * Outputs information-retrieval statistics per class. <p/>
 *
 * -k <br/>
 * Outputs information-theoretic statistics. <p/>
 *
 * -classifications "weka.classifiers.evaluation.output.prediction.AbstractOutput + options" <br/>
 * Uses the specified class for generating the classification output.
 * E.g.: weka.classifiers.evaluation.output.prediction.PlainText
 * or  : weka.classifiers.evaluation.output.prediction.CSV
 *
 * -p range <br/>
 * Outputs predictions for test instances (or the train instances if no test
 * instances provided and -no-cv is used), along with the attributes in the specified range
 * (and nothing else). Use '-p 0' if no attributes are desired. <p/>
 * Deprecated: use "-classifications ..." instead. <p/>
 *
 * -distribution <br/>
 * Outputs the distribution instead of only the prediction
 * in conjunction with the '-p' option (only nominal classes). <p/>
 * Deprecated: use "-classifications ..." instead. <p/>
 *
 * -r <br/>
 * Outputs cumulative margin distribution (and nothing else). <p/>
 *
 * -g <br/>
 * Only for classifiers that implement "Graphable." Outputs
 * the graph representation of the classifier (and nothing
 * else). <p/>
 *
 * -xml filename | xml-string <br/>
 * Retrieves the options from the XML-data instead of the command line. <p/>
 *
 * -threshold-file file <br/>
 * The file to save the threshold data to.
 * The format is determined by the extensions, e.g., '.arff' for ARFF
 * format or '.csv' for CSV. <p/>
 *
 * -threshold-label label <br/>
 * The class label to determine the threshold data for
 * (default is the first label) <p/>
 *
 * ------------------------------------------------------------------- <p/>
 *
 * Example usage as the main of a classifier (called FunkyClassifier):
 * <code> <pre>
 * public static void main(String [] args) {
 *   runClassifier(new FunkyClassifier(), args);
 * }
 * </pre> </code>
 * <p/>
 *
 * ------------------------------------------------------------------ <p/>
 *
 * Example usage from within an application:
 * <code> <pre>
 * Instances trainInstances = ... instances got from somewhere
 * Instances testInstances = ... instances got from somewhere
 * Classifier scheme = ... scheme got from somewhere
 *
 * Evaluation evaluation = new Evaluation(trainInstances);
 * evaluation.evaluateModel(scheme, testInstances);
 * System.out.println(evaluation.toSummaryString());
 * </pre> </code>
 *
 *
 * @author   Eibe Frank (eibe@cs.waikato.ac.nz)
 * @author   Len Trigg (trigg@cs.waikato.ac.nz)
 * @version  $Revision: 6804 $
 */
@SuppressWarnings("deprecation")
public class Evaluation


  implements Summarizable, RevisionHandler {

  //RANKING BEGIN
  /** the dataset used for ranking measures. */
  private Instances rnkData;
  private double tau;
  private double[] spearman = new double[2];
  int[][] rnkDataMatrix;
  public int train_testInst=0;
  //RANKING END
  /** The number of classes. */
  protected int m_NumClasses;

  /** The number of folds for a cross-validation. */
  protected int m_NumFolds;

  /** The weight of all incorrectly classified instances. */
  protected double m_Incorrect;

  /** The weight of all correctly classified instances. */
  protected double m_Correct;

  /** The weight of all unclassified instances. */
  protected double m_Unclassified;

  /*** The weight of all instances that had no class assigned to them. */
  protected double m_MissingClass;

  /** The weight of all instances that had a class assigned to them. */
  protected double m_WithClass;

  /** Array for storing the confusion matrix. */
  protected double [][] m_ConfusionMatrix;

  /** The names of the classes. */
  protected String [] m_ClassNames;

  /** Is the class nominal or numeric? */
  protected boolean m_ClassIsNominal;
  
  //RANKING BEGIN
  /** Is the class a Ranking?*/
  protected boolean m_ClassIsRanking;
  //RANKING END
  /** The prior probabilities of the classes. */
  protected double [] m_ClassPriors;

  /** The sum of counts for priors. */
  protected double m_ClassPriorsSum;

  /** The cost matrix (if given). */
  protected CostMatrix m_CostMatrix;

  /** The total cost of predictions (includes instance weights). */
  protected double m_TotalCost;

  /** Sum of errors. */
  protected double m_SumErr;

  /** Sum of absolute errors. */
  protected double m_SumAbsErr;

  /** Sum of squared errors. */
  protected double m_SumSqrErr;

  /** Sum of class values. */
  protected double m_SumClass;

  /** Sum of squared class values. */
  protected double m_SumSqrClass;

  /*** Sum of predicted values. */
  protected double m_SumPredicted;

  /** Sum of squared predicted values. */
  protected double m_SumSqrPredicted;

  /** Sum of predicted * class values. */
  protected double m_SumClassPredicted;

  /** Sum of absolute errors of the prior. */
  protected double m_SumPriorAbsErr;

  /** Sum of absolute errors of the prior. */
  protected double m_SumPriorSqrErr;

  /** Total Kononenko & Bratko Information. */
  protected double m_SumKBInfo;

  /*** Resolution of the margin histogram. */
  protected static int k_MarginResolution = 500;

  /** Cumulative margin distribution. */
  protected double m_MarginCounts [];

  /** Number of non-missing class training instances seen. */
  protected int m_NumTrainClassVals;

  /** Array containing all numeric training class values seen. */
  protected double [] m_TrainClassVals;

  /** Array containing all numeric training class weights. */
  protected double [] m_TrainClassWeights;

  /** Numeric class estimator for prior. */
  protected UnivariateKernelEstimator m_PriorEstimator;

  /** Whether complexity statistics are available. */
  protected boolean m_ComplexityStatisticsAvailable = true;

  /**
   * The minimum probablility accepted from an estimator to avoid
   * taking log(0) in Sf calculations.
   */
  protected static final double MIN_SF_PROB = Double.MIN_VALUE;

  /** Total entropy of prior predictions. */
  protected double m_SumPriorEntropy;

  /** Total entropy of scheme predictions. */
  protected double m_SumSchemeEntropy;

  /** Whether coverage statistics are available. */
  protected boolean m_CoverageStatisticsAvailable = true;

  /**  The confidence level used for coverage statistics. */
  protected double m_ConfLevel = 0.95;

  /** Total size of predicted regions at the given confidence level. */
  protected double m_TotalSizeOfRegions;

  /** Total coverage of test cases at the given confidence level. */
  protected double m_TotalCoverage;

  /** Minimum target value. */
  protected double m_MinTarget;

  /** Maximum target value. */
  protected double m_MaxTarget;

  /** The list of predictions that have been generated (for computing AUC). */
  protected FastVector m_Predictions;

  /** enables/disables the use of priors, e.g., if no training set is
   * present in case of de-serialized schemes. */
  protected boolean m_NoPriors = false;

  /** The header of the training set. */
  protected Instances m_Header;
  
  /**
   * Initializes all the counters for the evaluation.
   * Use <code>useNoPriors()</code> if the dataset is the test set and you
   * can't initialize with the priors from the training set via
   * <code>setPriors(Instances)</code>.
   *
   * @param data 	set of training instances, to get some header
   * 			information and prior class distribution information
   * @throws Exception 	if the class is not defined
   * @see 		#useNoPriors()
   * @see 		#setPriors(Instances)
   */
  public Evaluation(Instances data) throws Exception {

    this(data, null);
  }

  /**
   * Initializes all the counters for the evaluation and also takes a
   * cost matrix as parameter.
   * Use <code>useNoPriors()</code> if the dataset is the test set and you
   * can't initialize with the priors from the training set via
   * <code>setPriors(Instances)</code>.
   *
   * @param data 	set of training instances, to get some header
   * 			information and prior class distribution information
   * @param costMatrix 	the cost matrix---if null, default costs will be used
   * @throws Exception 	if cost matrix is not compatible with
   * 			data, the class is not defined or the class is numeric
   * @see 		#useNoPriors()
   * @see 		#setPriors(Instances)
   */
  public Evaluation(Instances data, CostMatrix costMatrix)
  throws Exception {
    //RANKING BEGIN
	  try{
	this.rnkData = data;
	if(rnkData.getLabels().size()!=0)
		rnkDataMatrix = new int[rnkData.getLabels().size()][2*rnkData.getLabels().size()];
	else{
		rnkDataMatrix = new int[RankUtilities.labels.size()][2*RankUtilities.labels.size()];
		rnkData.setLabels(RankUtilities.labels);
		data.setLabels(RankUtilities.labels);
	}
	m_ClassIsRanking = data.classAttribute() instanceof PreferenceAttribute || data.instance(0) instanceof PreferenceDenseInstance;
	
    m_Header = new Instances(data, 0);
    m_NumClasses = data.numClasses();
    m_NumFolds = 1;
    m_ClassIsNominal = data.classAttribute().isNominal();

    if(data.classAttribute().isRanking()) m_ClassIsNominal = true;

    

    if (m_ClassIsNominal) {
      m_ConfusionMatrix = new double [m_NumClasses][m_NumClasses];
      
      m_ClassNames = new String [m_NumClasses];
      
      
      for(int i = 0; i < m_NumClasses; i++) {
    	  m_ClassNames[i] = data.classAttribute().value(i);
      }
    }
    m_CostMatrix = costMatrix;
    if (m_CostMatrix != null) {
      if (!m_ClassIsNominal) {
        throw new Exception("Class has to be nominal if cost matrix given!");
      }
      if (m_CostMatrix.size() != m_NumClasses) {
        throw new Exception("Cost matrix not compatible with data!");
      }
    }
    m_ClassPriors = new double [m_NumClasses];
    setPriors(data);
    m_MarginCounts = new double [k_MarginResolution + 1];
	  }
	  catch(Exception e){//RANKING END

		    m_Header = new Instances(data, 0);
		    m_NumClasses = data.numClasses();
		    m_NumFolds = 1;
		    m_ClassIsNominal = data.classAttribute().isNominal();

		    if (m_ClassIsNominal) {
		      m_ConfusionMatrix = new double [m_NumClasses][m_NumClasses];
		      m_ClassNames = new String [m_NumClasses];
		      for(int i = 0; i < m_NumClasses; i++) {
		        m_ClassNames[i] = data.classAttribute().value(i);
		      }
		    }
		    m_CostMatrix = costMatrix;
		    if (m_CostMatrix != null) {
		      if (!m_ClassIsNominal) {
		        throw new Exception("Class has to be nominal if cost matrix given!");
		      }
		      if (m_CostMatrix.size() != m_NumClasses) {
		        throw new Exception("Cost matrix not compatible with data!");
		      }
		    }
		    m_ClassPriors = new double [m_NumClasses];
		    setPriors(data);
		    m_MarginCounts = new double [k_MarginResolution + 1];
	  }
  }

  /**
   * Returns the header of the underlying dataset.
   *
   * @return		the header information
   */
  public Instances getHeader() {
    return m_Header;
  }

  /**
   * Returns the area under ROC for those predictions that have been collected
   * in the evaluateClassifier(Classifier, Instances) method. Returns
   * Utils.missingValue() if the area is not available.
   *
   * @param classIndex the index of the class to consider as "positive"
   * @return the area under the ROC curve or not a number
   */
  public double areaUnderROC(int classIndex) {

    // Check if any predictions have been collected
    if (m_Predictions == null) {
      return Utils.missingValue();
    } else {
      ThresholdCurve tc = new ThresholdCurve();
      Instances result = tc.getCurve(m_Predictions, classIndex);
      return ThresholdCurve.getROCArea(result);
    }
  }

  /**
   * Calculates the weighted (by class size) AUC.
   *
   * @return the weighted AUC.
   */
  public double weightedAreaUnderROC() {
    double[] classCounts = new double[m_NumClasses];
    double classCountSum = 0;

    for (int i = 0; i < m_NumClasses; i++) {
      for (int j = 0; j < m_NumClasses; j++) {
        classCounts[i] += m_ConfusionMatrix[i][j];
      }
      classCountSum += classCounts[i];
    }

    double aucTotal = 0;
    for(int i = 0; i < m_NumClasses; i++) {
      double temp = areaUnderROC(i);
      if (!Utils.isMissingValue(temp)) {
        aucTotal += (temp * classCounts[i]);
      }
    }

    return aucTotal / classCountSum;
  }

  /**
   * Returns a copy of the confusion matrix.
   *
   * @return a copy of the confusion matrix as a two-dimensional array
   */
  public double[][] confusionMatrix() {

    double[][] newMatrix = new double[m_ConfusionMatrix.length][0];

    for (int i = 0; i < m_ConfusionMatrix.length; i++) {
      newMatrix[i] = new double[m_ConfusionMatrix[i].length];
      System.arraycopy(m_ConfusionMatrix[i], 0, newMatrix[i], 0,
          m_ConfusionMatrix[i].length);
    }
    return newMatrix;
  }

  /**
   * Performs a (stratified if class is nominal) cross-validation
   * for a classifier on a set of instances. Now performs
   * a deep copy of the classifier before each call to
   * buildClassifier() (just in case the classifier is not
   * initialized properly).
   *
   * @param classifier the classifier with any options set.
   * @param data the data on which the cross-validation is to be
   * performed
   * @param numFolds the number of folds for the cross-validation
   * @param random random number generator for randomization
   * @param forPredictionsPrinting varargs parameter that, if supplied, is
   * expected to hold a weka.classifiers.evaluation.output.prediction.AbstractOutput
   * object
   * @throws Exception if a classifier could not be generated
   * successfully or the class is not defined
   */
  public void crossValidateModel(Classifier classifier,
                                 Instances data, int numFolds, Random random,
                                 Object... forPredictionsPrinting)
  throws Exception {

    // Make a copy of the data we can reorder
    data = new Instances(data);
    data.randomize(random);
    if (data.classAttribute().isNominal() || data.classAttribute().isRanking()) {
      data.stratify(numFolds);
    }

    // We assume that the first element is a
    // weka.classifiers.evaluation.output.prediction.AbstractOutput object
    AbstractOutput classificationOutput = null;
    if (forPredictionsPrinting.length > 0) {
      // print the header first
      classificationOutput = (AbstractOutput) forPredictionsPrinting[0];
      classificationOutput.setHeader(data);
      classificationOutput.printHeader();
    }

    // Do the folds
    for (int i = 0; i < numFolds; i++) {
      Instances train = data.trainCV(numFolds, i, random);
      setPriors(train);
      Classifier copiedClassifier = AbstractClassifier.makeCopy(classifier);
      copiedClassifier.buildClassifier(train);
      Instances test = data.testCV(numFolds, i);
      evaluateModel(copiedClassifier, test, forPredictionsPrinting);
    }
    m_NumFolds = numFolds;
   
    if (classificationOutput != null)
      classificationOutput.printFooter();
    
  }

  /**
   * Performs a (stratified if class is nominal) cross-validation
   * for a classifier on a set of instances.
   *
   * @param classifierString a string naming the class of the classifier
   * @param data the data on which the cross-validation is to be
   * performed
   * @param numFolds the number of folds for the cross-validation
   * @param options the options to the classifier. Any options
   * @param random the random number generator for randomizing the data
   * accepted by the classifier will be removed from this array.
   * @throws Exception if a classifier could not be generated
   * successfully or the class is not defined
   */
  public void crossValidateModel(String classifierString,
      Instances data, int numFolds,
      String[] options, Random random)
    throws Exception {

    crossValidateModel(AbstractClassifier.forName(classifierString, options),
        data, numFolds, random);
  }

  /**
   * Evaluates a classifier with the options given in an array of
   * strings. <p/>
   *
   * Valid options are: <p/>
   *
   * -t filename <br/>
   * Name of the file with the training data. (required) <p/>
   *
   * -T filename <br/>
   * Name of the file with the test data. If missing a cross-validation
   * is performed. <p/>
   *
   * -c index <br/>
   * Index of the class attribute (1, 2, ...; default: last). <p/>
   *
   * -x number <br/>
   * The number of folds for the cross-validation (default: 10). <p/>
   *
   * -no-cv <br/>
   * No cross validation.  If no test file is provided, no evaluation
   * is done. <p/>
   *
   * -split-percentage percentage <br/>
   * Sets the percentage for the train/test set split, e.g., 66. <p/>
   *
   * -preserve-order <br/>
   * Preserves the order in the percentage split instead of randomizing
   * the data first with the seed value ('-s'). <p/>
   *
   * -s seed <br/>
   * Random number seed for the cross-validation and percentage split
   * (default: 1). <p/>
   *
   * -m filename <br/>
   * The name of a file containing a cost matrix. <p/>
   *
   * -l filename <br/>
   * Loads classifier from the given file. In case the filename ends with
   * ".xml",a PMML file is loaded or, if that fails, options are loaded from XML. <p/>
   *
   * -d filename <br/>
   * Saves classifier built from the training data into the given file. In case
   * the filename ends with ".xml" the options are saved XML, not the model. <p/>
   *
   * -v <br/>
   * Outputs no statistics for the training data. <p/>
   *
   * -o <br/>
   * Outputs statistics only, not the classifier. <p/>
   *
   * -i <br/>
   * Outputs detailed information-retrieval statistics per class. <p/>
   *
   * -k <br/>
   * Outputs information-theoretic statistics. <p/>
   *
   * -classifications "weka.classifiers.evaluation.output.prediction.AbstractOutput + options" <br/>
   * Uses the specified class for generating the classification output.
   * E.g.: weka.classifiers.evaluation.output.prediction.PlainText
   * or  : weka.classifiers.evaluation.output.prediction.CSV
   *
   * -p range <br/>
   * Outputs predictions for test instances (or the train instances if no test
   * instances provided and -no-cv is used), along with the attributes in the specified range
   * (and nothing else). Use '-p 0' if no attributes are desired. <p/>
   * Deprecated: use "-classifications ..." instead. <p/>
   *
   * -distribution <br/>
   * Outputs the distribution instead of only the prediction
   * in conjunction with the '-p' option (only nominal classes). <p/>
   * Deprecated: use "-classifications ..." instead. <p/>
   *
   * -r <br/>
   * Outputs cumulative margin distribution (and nothing else). <p/>
   *
   * -g <br/>
   * Only for classifiers that implement "Graphable." Outputs
   * the graph representation of the classifier (and nothing
   * else). <p/>
   *
   * -xml filename | xml-string <br/>
   * Retrieves the options from the XML-data instead of the command line. <p/>
   *
   * -threshold-file file <br/>
   * The file to save the threshold data to.
   * The format is determined by the extensions, e.g., '.arff' for ARFF
   * format or '.csv' for CSV. <p/>
   *
   * -threshold-label label <br/>
   * The class label to determine the threshold data for
   * (default is the first label) <p/>
   *
   * @param classifierString class of machine learning classifier as a string
   * @param options the array of string containing the options
   * @throws Exception if model could not be evaluated successfully
   * @return a string describing the results
   */
  public static String evaluateModel(String classifierString,
      String [] options) throws Exception {

    Classifier classifier;

    // Create classifier
    try {
      classifier =
        //  (Classifier)Class.forName(classifierString).newInstance();
        AbstractClassifier.forName(classifierString, null);
    } catch (Exception e) {
      throw new Exception("Can't find class with name "
          + classifierString + '.');
    }
    return evaluateModel(classifier, options);
  }

  /**
   * A test method for this class. Just extracts the first command line
   * argument as a classifier class name and calls evaluateModel.
   * @param args an array of command line arguments, the first of which
   * must be the class name of a classifier.
   */
  public static void main(String [] args) {

    try {
      if (args.length == 0) {
        throw new Exception("The first argument must be the class name"
            + " of a classifier");
      }
      String classifier = args[0];
      args[0] = "";
      System.out.println(evaluateModel(classifier, args));
    } catch (Exception ex) {
      ex.printStackTrace();
      System.err.println(ex.getMessage());
    }
  }

  /**
   * Evaluates a classifier with the options given in an array of
   * strings. <p/>
   *
   * Valid options are: <p/>
   *
   * -t name of training file <br/>
   * Name of the file with the training data. (required) <p/>
   *
   * -T name of test file <br/>
   * Name of the file with the test data. If missing a cross-validation
   * is performed. <p/>
   *
   * -c class index <br/>
   * Index of the class attribute (1, 2, ...; default: last). <p/>
   *
   * -x number of folds <br/>
   * The number of folds for the cross-validation (default: 10). <p/>
   *
   * -no-cv <br/>
   * No cross validation.  If no test file is provided, no evaluation
   * is done. <p/>
   *
   * -split-percentage percentage <br/>
   * Sets the percentage for the train/test set split, e.g., 66. <p/>
   *
   * -preserve-order <br/>
   * Preserves the order in the percentage split instead of randomizing
   * the data first with the seed value ('-s'). <p/>
   *
   * -s seed <br/>
   * Random number seed for the cross-validation and percentage split
   * (default: 1). <p/>
   *
   * -m file with cost matrix <br/>
   * The name of a file containing a cost matrix. <p/>
   *
   * -l filename <br/>
   * Loads classifier from the given file. In case the filename ends with
   * ".xml",a PMML file is loaded or, if that fails, options are loaded from XML. <p/>
   *
   * -d filename <br/>
   * Saves classifier built from the training data into the given file. In case
   * the filename ends with ".xml" the options are saved XML, not the model. <p/>
   *
   * -v <br/>
   * Outputs no statistics for the training data. <p/>
   *
   * -o <br/>
   * Outputs statistics only, not the classifier. <p/>
   *
   * -i <br/>
   * Outputs detailed information-retrieval statistics per class. <p/>
   *
   * -k <br/>
   * Outputs information-theoretic statistics. <p/>
   *
   * -classifications "weka.classifiers.evaluation.output.prediction.AbstractOutput + options" <br/>
   * Uses the specified class for generating the classification output.
   * E.g.: weka.classifiers.evaluation.output.prediction.PlainText
   * or  : weka.classifiers.evaluation.output.prediction.CSV
   *
   * -p range <br/>
   * Outputs predictions for test instances (or the train instances if no test
   * instances provided and -no-cv is used), along with the attributes in the specified range
   * (and nothing else). Use '-p 0' if no attributes are desired. <p/>
   * Deprecated: use "-classifications ..." instead. <p/>
   *
   * -distribution <br/>
   * Outputs the distribution instead of only the prediction
   * in conjunction with the '-p' option (only nominal classes). <p/>
   * Deprecated: use "-classifications ..." instead. <p/>
   *
   * -r <br/>
   * Outputs cumulative margin distribution (and nothing else). <p/>
   *
   * -g <br/>
   * Only for classifiers that implement "Graphable." Outputs
   * the graph representation of the classifier (and nothing
   * else). <p/>
   *
   * -xml filename | xml-string <br/>
   * Retrieves the options from the XML-data instead of the command line. <p/>
   *
   * @param classifier machine learning classifier
   * @param options the array of string containing the options
   * @throws Exception if model could not be evaluated successfully
   * @return a string describing the results
   */
  public static String evaluateModel(Classifier classifier,
      String [] options) throws Exception {

    Instances train = null, tempTrain, test = null, template = null;
    int seed = 1, folds = 10, classIndex = -1;
    boolean noCrossValidation = false;
    String trainFileName, testFileName, sourceClass,
    classIndexString, seedString, foldsString, objectInputFileName,
    objectOutputFileName;
    boolean noOutput = false,
    trainStatistics = true,
    printMargins = false, printComplexityStatistics = false,
    printGraph = false, classStatistics = false, printSource = false;
    StringBuffer text = new StringBuffer();
    DataSource trainSource = null, testSource = null;
    ObjectInputStream objectInputStream = null;
    BufferedInputStream xmlInputStream = null;
    CostMatrix costMatrix = null;
    StringBuffer schemeOptionsText = null;
    long trainTimeStart = 0, trainTimeElapsed = 0,
    testTimeStart = 0, testTimeElapsed = 0;
    String xml = "";
    String[] optionsTmp = null;
    Classifier classifierBackup;
    Classifier classifierClassifications = null;
    int actualClassIndex = -1;  // 0-based class index
    String splitPercentageString = "";
    double splitPercentage = -1;
    boolean preserveOrder = false;
    boolean trainSetPresent = false;
    boolean testSetPresent = false;
    String thresholdFile;
    String thresholdLabel;
    StringBuffer predsBuff = null; // predictions from cross-validation
    AbstractOutput classificationOutput = null;

    // help requested?
    if (Utils.getFlag("h", options) || Utils.getFlag("help", options)) {

      // global info requested as well?
      boolean globalInfo = Utils.getFlag("synopsis", options) ||
        Utils.getFlag("info", options);

      throw new Exception("\nHelp requested."
          + makeOptionString(classifier, globalInfo));
    }

    try {
      // do we get the input from XML instead of normal parameters?
      xml = Utils.getOption("xml", options);
      if (!xml.equals(""))
        options = new XMLOptions(xml).toArray();

      // is the input model only the XML-Options, i.e. w/o built model?
      optionsTmp = new String[options.length];
      for (int i = 0; i < options.length; i++)
        optionsTmp[i] = options[i];

      String tmpO = Utils.getOption('l', optionsTmp);
      //if (Utils.getOption('l', optionsTmp).toLowerCase().endsWith(".xml")) {
      if (tmpO.endsWith(".xml")) {
        // try to load file as PMML first
        boolean success = false;
        try {
          PMMLModel pmmlModel = PMMLFactory.getPMMLModel(tmpO);
          if (pmmlModel instanceof PMMLClassifier) {
            classifier = ((PMMLClassifier)pmmlModel);
            success = true;
          }
        } catch (IllegalArgumentException ex) {
          success = false;
        }
        if (!success) {
          // load options from serialized data  ('-l' is automatically erased!)
          XMLClassifier xmlserial = new XMLClassifier();
          OptionHandler cl = (OptionHandler) xmlserial.read(Utils.getOption('l', options));

          // merge options
          optionsTmp = new String[options.length + cl.getOptions().length];
          System.arraycopy(cl.getOptions(), 0, optionsTmp, 0, cl.getOptions().length);
          System.arraycopy(options, 0, optionsTmp, cl.getOptions().length, options.length);
          options = optionsTmp;
        }
      }

      noCrossValidation = Utils.getFlag("no-cv", options);
      // Get basic options (options the same for all schemes)
      classIndexString = Utils.getOption('c', options);
      if (classIndexString.length() != 0) {
        if (classIndexString.equals("first"))
          classIndex = 1;
        else if (classIndexString.equals("last"))
          classIndex = -1;
        else
          classIndex = Integer.parseInt(classIndexString);
      }
      trainFileName = Utils.getOption('t', options);
      objectInputFileName = Utils.getOption('l', options);
      objectOutputFileName = Utils.getOption('d', options);
      testFileName = Utils.getOption('T', options);
      foldsString = Utils.getOption('x', options);
      if (foldsString.length() != 0) {
        folds = Integer.parseInt(foldsString);
      }
      seedString = Utils.getOption('s', options);
      if (seedString.length() != 0) {
        seed = Integer.parseInt(seedString);
      }
      if (trainFileName.length() == 0) {
        if (objectInputFileName.length() == 0) {
          throw new Exception("No training file and no object input file given.");
        }
        if (testFileName.length() == 0) {
          throw new Exception("No training file and no test file given.");
        }
      } else if ((objectInputFileName.length() != 0) &&
          ((!(classifier instanceof UpdateableClassifier)) ||
           (testFileName.length() == 0))) {
        throw new Exception("Classifier not incremental, or no " +
            "test file provided: can't "+
            "use both train and model file.");
      }
      try {
        if (trainFileName.length() != 0) {
          trainSetPresent = true;
          trainSource = new DataSource(trainFileName);
        }
        if (testFileName.length() != 0) {
          testSetPresent = true;
          testSource = new DataSource(testFileName);
        }
        if (objectInputFileName.length() != 0) {
          if (objectInputFileName.endsWith(".xml")) {
            // if this is the case then it means that a PMML classifier was
            // successfully loaded earlier in the code
            objectInputStream = null;
            xmlInputStream = null;
          } else {
            InputStream is = new FileInputStream(objectInputFileName);
            if (objectInputFileName.endsWith(".gz")) {
              is = new GZIPInputStream(is);
            }
            // load from KOML?
            if (!(objectInputFileName.endsWith(".koml") && KOML.isPresent()) ) {
              objectInputStream = new ObjectInputStream(is);
              xmlInputStream    = null;
            }
            else {
              objectInputStream = null;
              xmlInputStream    = new BufferedInputStream(is);
            }
          }
        }
      } catch (Exception e) {
        throw new Exception("Can't open file " + e.getMessage() + '.');
      }
      if (testSetPresent) {
        template = test = testSource.getStructure();
        if (classIndex != -1) {
          test.setClassIndex(classIndex - 1);
        } else {
          if ( (test.classIndex() == -1) || (classIndexString.length() != 0) )
            test.setClassIndex(test.numAttributes() - 1);
        }
        actualClassIndex = test.classIndex();
      }
      else {
        // percentage split
        splitPercentageString = Utils.getOption("split-percentage", options);
        if (splitPercentageString.length() != 0) {
          if (foldsString.length() != 0)
            throw new Exception(
                "Percentage split cannot be used in conjunction with "
                + "cross-validation ('-x').");
          splitPercentage = Double.parseDouble(splitPercentageString);
          if ((splitPercentage <= 0) || (splitPercentage >= 100))
            throw new Exception("Percentage split value needs be >0 and <100.");
        }
        else {
          splitPercentage = -1;
        }
        preserveOrder = Utils.getFlag("preserve-order", options);
        if (preserveOrder) {
          if (splitPercentage == -1)
            throw new Exception("Percentage split ('-percentage-split') is missing.");
        }
        // create new train/test sources
        if (splitPercentage > 0) {
          testSetPresent = true;
          Instances tmpInst = trainSource.getDataSet(actualClassIndex);
          if (!preserveOrder)
            tmpInst.randomize(new Random(seed));
          int trainSize = 
            (int) Math.round(tmpInst.numInstances() * splitPercentage / 100);
          int testSize  = tmpInst.numInstances() - trainSize;
          Instances trainInst = new Instances(tmpInst, 0, trainSize);
          Instances testInst  = new Instances(tmpInst, trainSize, testSize);
          trainSource = new DataSource(trainInst);
          testSource  = new DataSource(testInst);
          template = test = testSource.getStructure();
          if (classIndex != -1) {
            test.setClassIndex(classIndex - 1);
          } else {
            if ( (test.classIndex() == -1) || (classIndexString.length() != 0) )
              test.setClassIndex(test.numAttributes() - 1);
          }
          actualClassIndex = test.classIndex();
        }
      }
      if (trainSetPresent) {
        template = train = trainSource.getStructure();
        if (classIndex != -1) {
          train.setClassIndex(classIndex - 1);
        } else {
          if ( (train.classIndex() == -1) || (classIndexString.length() != 0) )
            train.setClassIndex(train.numAttributes() - 1);
        }
        actualClassIndex = train.classIndex();
        if (!(classifier instanceof weka.classifiers.misc.InputMappedClassifier)) {
          if ((testSetPresent) && !test.equalHeaders(train)) {
            throw new IllegalArgumentException("Train and test file not compatible!\n" + test.equalHeadersMsg(train));
          }
        }
      }
      if (template == null) {
        throw new Exception("No actual dataset provided to use as template");
      }
      costMatrix = handleCostOption(
          Utils.getOption('m', options), template.numClasses());

      classStatistics = Utils.getFlag('i', options);
      noOutput = Utils.getFlag('o', options);
      trainStatistics = !Utils.getFlag('v', options);
      printComplexityStatistics = Utils.getFlag('k', options);
      printMargins = Utils.getFlag('r', options);
      printGraph = Utils.getFlag('g', options);
      sourceClass = Utils.getOption('z', options);
      printSource = (sourceClass.length() != 0);
      thresholdFile = Utils.getOption("threshold-file", options);
      thresholdLabel = Utils.getOption("threshold-label", options);

      String classifications = Utils.getOption("classifications", options);
      String classificationsOld = Utils.getOption("p", options);
      if (classifications.length() > 0) {
        noOutput = true;
        classificationOutput = AbstractOutput.fromCommandline(classifications);
        classificationOutput.setHeader(template);
      }
      // backwards compatible with old "-p range" and "-distribution" options
      else if (classificationsOld.length() > 0) {
        noOutput = true;
        classificationOutput = new PlainText();
        classificationOutput.setHeader(template);
        if (!classificationsOld.equals("0"))
          classificationOutput.setAttributes(classificationsOld);
        classificationOutput.setOutputDistribution(Utils.getFlag("distribution", options));
      }
      // -distribution flag needs -p option
      else {
        if (Utils.getFlag("distribution", options))
          throw new Exception("Cannot print distribution without '-p' option!");
      }

      // if no training file given, we don't have any priors
      if ( (!trainSetPresent) && (printComplexityStatistics) )
        throw new Exception("Cannot print complexity statistics ('-k') without training file ('-t')!");

      // If a model file is given, we can't process
      // scheme-specific options
      if (objectInputFileName.length() != 0) {
        Utils.checkForRemainingOptions(options);
      } else {

        // Set options for classifier
        if (classifier instanceof OptionHandler) {
          for (int i = 0; i < options.length; i++) {
            if (options[i].length() != 0) {
              if (schemeOptionsText == null) {
                schemeOptionsText = new StringBuffer();
              }
              if (options[i].indexOf(' ') != -1) {
                schemeOptionsText.append('"' + options[i] + "\" ");
              } else {
                schemeOptionsText.append(options[i] + " ");
              }
            }
          }
          ((OptionHandler)classifier).setOptions(options);
        }
      }

      Utils.checkForRemainingOptions(options);
    } catch (Exception e) {
      throw new Exception("\nWeka exception: " + e.getMessage()
          + makeOptionString(classifier, false));
    }

    if (objectInputFileName.length() != 0) {
      // Load classifier from file
      if (objectInputStream != null) {
        classifier = (Classifier) objectInputStream.readObject();
        // try and read a header (if present)
        Instances savedStructure = null;
        try {
          savedStructure = (Instances) objectInputStream.readObject();
        } catch (Exception ex) {
          // don't make a fuss
        }
        if (savedStructure != null) {
          // test for compatibility with template
          if (!template.equalHeaders(savedStructure)) {
            throw new Exception("training and test set are not compatible\n" + template.equalHeadersMsg(savedStructure));
          }
        }
        objectInputStream.close();
      }
      else if (xmlInputStream != null) {
        // whether KOML is available has already been checked (objectInputStream would null otherwise)!
        classifier = (Classifier) KOML.read(xmlInputStream);
        xmlInputStream.close();
      }
    }
    
    // Setup up evaluation objects
    Evaluation trainingEvaluation = new Evaluation(new Instances(template, 0), costMatrix);
    Evaluation testingEvaluation = new Evaluation(new Instances(template, 0), costMatrix);
    if (classifier instanceof weka.classifiers.misc.InputMappedClassifier) {
      Instances mappedClassifierHeader = 
        ((weka.classifiers.misc.InputMappedClassifier)classifier).
          getModelHeader(new Instances(template, 0));
            
      trainingEvaluation = new Evaluation(new Instances(mappedClassifierHeader, 0), costMatrix);
      testingEvaluation = new Evaluation(new Instances(mappedClassifierHeader, 0), costMatrix);
    }

    // disable use of priors if no training file given
    if (!trainSetPresent)
      testingEvaluation.useNoPriors();

    // backup of fully setup classifier for cross-validation
    classifierBackup = AbstractClassifier.makeCopy(classifier);

    // Build the classifier if no object file provided
    if ((classifier instanceof UpdateableClassifier) &&
        (testSetPresent || noCrossValidation) &&
        (costMatrix == null) &&
        (trainSetPresent)) {
      // Build classifier incrementally
      trainingEvaluation.setPriors(train);
      testingEvaluation.setPriors(train);
      trainTimeStart = System.currentTimeMillis();
      if (objectInputFileName.length() == 0) {
        classifier.buildClassifier(train);
      }
      Instance trainInst;
      while (trainSource.hasMoreElements(train)) {
        trainInst = trainSource.nextElement(train);
        trainingEvaluation.updatePriors(trainInst);
        testingEvaluation.updatePriors(trainInst);
        ((UpdateableClassifier)classifier).updateClassifier(trainInst);
      }
      trainTimeElapsed = System.currentTimeMillis() - trainTimeStart;
    } else if (objectInputFileName.length() == 0) {
      // Build classifier in one go
      tempTrain = trainSource.getDataSet(actualClassIndex);
      
      if (classifier instanceof weka.classifiers.misc.InputMappedClassifier &&
          !trainingEvaluation.getHeader().equalHeaders(tempTrain)) {
        // we need to make a new dataset that maps the training instances to
        // the structure expected by the mapped classifier - this is only
        // to ensure that the structure and priors computed by the *testing*
        // evaluation object is correct with respect to the mapped classifier
        Instances mappedClassifierDataset = 
          ((weka.classifiers.misc.InputMappedClassifier)classifier).
            getModelHeader(new Instances(template, 0));
        for (int zz = 0; zz < tempTrain.numInstances(); zz++) {
          Instance mapped = ((weka.classifiers.misc.InputMappedClassifier)classifier).
            constructMappedInstance(tempTrain.instance(zz));
          mappedClassifierDataset.add(mapped);
        }
        tempTrain = mappedClassifierDataset;
      }
      
      trainingEvaluation.setPriors(tempTrain);
      testingEvaluation.setPriors(tempTrain);
      trainTimeStart = System.currentTimeMillis();
      classifier.buildClassifier(tempTrain);
      trainTimeElapsed = System.currentTimeMillis() - trainTimeStart;
    }

    // backup of fully trained classifier for printing the classifications
    if (classificationOutput != null) {
      classifierClassifications = AbstractClassifier.makeCopy(classifier);
      if (classifier instanceof weka.classifiers.misc.InputMappedClassifier) {
        classificationOutput.setHeader(trainingEvaluation.getHeader());
      }
    }

    // Save the classifier if an object output file is provided
    if (objectOutputFileName.length() != 0) {
      OutputStream os = new FileOutputStream(objectOutputFileName);
      // binary
      if (!(objectOutputFileName.endsWith(".xml") || (objectOutputFileName.endsWith(".koml") && KOML.isPresent()))) {
        if (objectOutputFileName.endsWith(".gz")) {
          os = new GZIPOutputStream(os);
        }
        ObjectOutputStream objectOutputStream = new ObjectOutputStream(os);
        objectOutputStream.writeObject(classifier);
        if (template != null) {
          objectOutputStream.writeObject(template);
        }
        objectOutputStream.flush();
        objectOutputStream.close();
      }
      // KOML/XML
      else {
        BufferedOutputStream xmlOutputStream = new BufferedOutputStream(os);
        if (objectOutputFileName.endsWith(".xml")) {
          XMLSerialization xmlSerial = new XMLClassifier();
          xmlSerial.write(xmlOutputStream, classifier);
        }
        else
          // whether KOML is present has already been checked
          // if not present -> ".koml" is interpreted as binary - see above
          if (objectOutputFileName.endsWith(".koml")) {
            KOML.write(xmlOutputStream, classifier);
          }
        xmlOutputStream.close();
      }
    }

    // If classifier is drawable output string describing graph
    if ((classifier instanceof Drawable) && (printGraph)){
      return ((Drawable)classifier).graph();
    }

    // Output the classifier as equivalent source
    if ((classifier instanceof Sourcable) && (printSource)){
      return wekaStaticWrapper((Sourcable) classifier, sourceClass);
    }

    // Output model
    if (!(noOutput || printMargins)) {
      if (classifier instanceof OptionHandler) {
        if (schemeOptionsText != null) {
          text.append("\nOptions: "+schemeOptionsText);
          text.append("\n");
        }
      }
      text.append("\n" + classifier.toString() + "\n");
    }

    if (!printMargins && (costMatrix != null)) {
      text.append("\n=== Evaluation Cost Matrix ===\n\n");
      text.append(costMatrix.toString());
    }

    // Output test instance predictions only
    if (classificationOutput != null) {
      DataSource source = testSource;
      predsBuff = new StringBuffer();
      classificationOutput.setBuffer(predsBuff);
      // no test set -> use train set
      if (source == null && noCrossValidation) {
        source = trainSource;
        predsBuff.append("\n=== Predictions on training data ===\n\n");
      } else {
        predsBuff.append("\n=== Predictions on test data ===\n\n");
      }
      if (source != null)
        classificationOutput.print(classifierClassifications, source);
    }

    // Compute error estimate from training data
    if ((trainStatistics) && (trainSetPresent)) {

      if ((classifier instanceof UpdateableClassifier) &&
          (testSetPresent) &&
          (costMatrix == null)) {

        // Classifier was trained incrementally, so we have to
        // reset the source.
        trainSource.reset();

        // Incremental testing
        train = trainSource.getStructure(actualClassIndex);
        testTimeStart = System.currentTimeMillis();
        Instance trainInst;
        while (trainSource.hasMoreElements(train)) {
          trainInst = trainSource.nextElement(train);
          trainingEvaluation.evaluateModelOnce((Classifier)classifier, trainInst);
        }
        testTimeElapsed = System.currentTimeMillis() - testTimeStart;
      } else {
        testTimeStart = System.currentTimeMillis();
        trainingEvaluation.evaluateModel(
            classifier, trainSource.getDataSet(actualClassIndex));
        testTimeElapsed = System.currentTimeMillis() - testTimeStart;
      }

      // Print the results of the training evaluation
      if (printMargins) {
        return trainingEvaluation.toCumulativeMarginDistributionString();
      } else {
        if (classificationOutput == null) {
          text.append("\nTime taken to build model: "
              + Utils.doubleToString(trainTimeElapsed / 1000.0,2)
              + " seconds");

          if (splitPercentage > 0)
            text.append("\nTime taken to test model on training split: ");
          else
            text.append("\nTime taken to test model on training data: ");
          text.append(Utils.doubleToString(testTimeElapsed / 1000.0,2) + " seconds");

          if (splitPercentage > 0)
            text.append(trainingEvaluation.toSummaryString("\n\n=== Error on training"
                  + " split ===\n", printComplexityStatistics));
          else
            text.append(trainingEvaluation.toSummaryString("\n\n=== Error on training"
                  + " data ===\n", printComplexityStatistics));

          //RANKING BEGIN
          if (template.classAttribute().isNominal() || template.classAttribute().isRanking()) {
          //RANKING END
            if (classStatistics) {
              text.append("\n\n" + trainingEvaluation.toClassDetailsString());
            }
            if (!noCrossValidation)
              text.append("\n\n" + trainingEvaluation.toMatrixString());
          }
        }
      }
    }

    // Compute proper error estimates
    if (testSource != null) {
      // Testing is on the supplied test data
      testSource.reset();
      test = testSource.getStructure(test.classIndex());
      Instance testInst;
      while (testSource.hasMoreElements(test)) {        
        testInst = testSource.nextElement(test);
        testingEvaluation.evaluateModelOnceAndRecordPrediction(
            (Classifier)classifier, testInst);
      }

      if (splitPercentage > 0) {
        if (classificationOutput == null) {
          text.append("\n\n" + testingEvaluation.
              toSummaryString("=== Error on test split ===\n",
                  printComplexityStatistics));
        }
      } else {
        if (classificationOutput == null) {
          text.append("\n\n" + testingEvaluation.
              toSummaryString("=== Error on test data ===\n",
                  printComplexityStatistics));
        }
      }

    } else if (trainSource != null) {
      if (!noCrossValidation) {
        // Testing is via cross-validation on training data
        Random random = new Random(seed);
        // use untrained (!) classifier for cross-validation
        classifier = AbstractClassifier.makeCopy(classifierBackup);
        if (classificationOutput == null) {
          testingEvaluation.crossValidateModel(classifier,
                                               trainSource.getDataSet(actualClassIndex),
                                               folds, random);
          if (template.classAttribute().isNumeric()) {
            text.append("\n\n\n" + testingEvaluation.
                        toSummaryString("=== Cross-validation ===\n",
                                        printComplexityStatistics));
          } else {
            text.append("\n\n\n" + testingEvaluation.
                        toSummaryString("=== Stratified " +
                                        "cross-validation ===\n",
                                        printComplexityStatistics));
          }
        } else {
          predsBuff = new StringBuffer();
          classificationOutput.setBuffer(predsBuff);
          predsBuff.append("\n=== Predictions under cross-validation ===\n\n");
          testingEvaluation.crossValidateModel(classifier,
                                               trainSource.getDataSet(actualClassIndex),
                                               folds, random, classificationOutput);
        }
      }
    }
    //RANKING BEGIN
    if (template.classAttribute().isNominal() || template.classAttribute().isRanking()) {
    //RANKING END
      if (classStatistics && !noCrossValidation && (classificationOutput == null)) {
        text.append("\n\n" + testingEvaluation.toClassDetailsString());
      }
      if (!noCrossValidation && (classificationOutput == null))
        text.append("\n\n" + testingEvaluation.toMatrixString());

    }

    // predictions from cross-validation?
    if (predsBuff != null) {
      text.append("\n" + predsBuff);
    }

    //RANKING BEGIN
    if ((thresholdFile.length() != 0) && (template.classAttribute().isNominal() || template.classAttribute().isRanking())) {
    //RANKING END
      int labelIndex = 0;
      if (thresholdLabel.length() != 0)
        labelIndex = template.classAttribute().indexOfValue(thresholdLabel);
      if (labelIndex == -1)
        throw new IllegalArgumentException(
            "Class label '" + thresholdLabel + "' is unknown!");
      ThresholdCurve tc = new ThresholdCurve();
      Instances result = tc.getCurve(testingEvaluation.predictions(), labelIndex);
      DataSink.write(thresholdFile, result);
    }
    
    return text.toString();
   // return "";
  }

  /**
   * Attempts to load a cost matrix.
   *
   * @param costFileName the filename of the cost matrix
   * @param numClasses the number of classes that should be in the cost matrix
   * (only used if the cost file is in old format).
   * @return a <code>CostMatrix</code> value, or null if costFileName is empty
   * @throws Exception if an error occurs.
   */
  protected static CostMatrix handleCostOption(String costFileName,
      int numClasses)
    throws Exception {

    if ((costFileName != null) && (costFileName.length() != 0)) {
      System.out.println(
          "NOTE: The behaviour of the -m option has changed between WEKA 3.0"
          +" and WEKA 3.1. -m now carries out cost-sensitive *evaluation*"
          +" only. For cost-sensitive *prediction*, use one of the"
          +" cost-sensitive metaschemes such as"
          +" weka.classifiers.meta.CostSensitiveClassifier or"
          +" weka.classifiers.meta.MetaCost");

      Reader costReader = null;
      try {
        costReader = new BufferedReader(new FileReader(costFileName));
      } catch (Exception e) {
        throw new Exception("Can't open file " + e.getMessage() + '.');
      }
      try {
        // First try as a proper cost matrix format
        return new CostMatrix(costReader);
      } catch (Exception ex) {
        try {
          // Now try as the poxy old format :-)
          //System.err.println("Attempting to read old format cost file");
          try {
            costReader.close(); // Close the old one
            costReader = new BufferedReader(new FileReader(costFileName));
          } catch (Exception e) {
            throw new Exception("Can't open file " + e.getMessage() + '.');
          }
          CostMatrix costMatrix = new CostMatrix(numClasses);
          //System.err.println("Created default cost matrix");
          costMatrix.readOldFormat(costReader);
          return costMatrix;
          //System.err.println("Read old format");
        } catch (Exception e2) {
          // re-throw the original exception
          //System.err.println("Re-throwing original exception");
          throw ex;
        }
      }
    } else {
      return null;
    }
  }

  /**
   * Evaluates the classifier on a given set of instances. Note that
   * the data must have exactly the same format (e.g. order of
   * attributes) as the data used to train the classifier! Otherwise
   * the results will generally be meaningless.
   *
   * @param classifier machine learning classifier
   * @param data set of test instances for evaluation
   * @param forPredictionsPrinting varargs parameter that, if supplied, is
   * expected to hold a weka.classifiers.evaluation.output.prediction.AbstractOutput
   * object
   * @return the predictions
   * @throws Exception if model could not be evaluated
   * successfully
   */
  public double[] evaluateModel(Classifier classifier,
                                Instances data,
                                Object... forPredictionsPrinting) throws Exception {
	  
    // for predictions printing
    AbstractOutput classificationOutput = null;

    double predictions[] = new double[data.numInstances()];

    if (forPredictionsPrinting.length > 0) {
      classificationOutput = (AbstractOutput) forPredictionsPrinting[0];
    }

    // Need to be able to collect predictions if appropriate (for AUC)

    for (int i = 0; i < data.numInstances(); i++) {
      predictions[i] = evaluateModelOnceAndRecordPrediction((Classifier)classifier,
          data.instance(i));
      if (classificationOutput != null)
        classificationOutput.printClassification(classifier, data.instance(i), i);
    }

    return predictions;
  }

  /**
   * Evaluates the supplied distribution on a single instance.
   *
   * @param dist the supplied distribution
   * @param instance the test instance to be classified
   * @param storePredictions whether to store predictions for nominal classifier
   * @return the prediction
   * @throws Exception if model could not be evaluated successfully
   */
  public double evaluationForSingleInstance(double[] dist, Instance instance,
                                            boolean storePredictions) throws Exception {

    double pred;

    
    
    //RANKING BEGIN
    //In case of a Ranking, the distribution array contains a ranking, because
    //Label Ranking algorithms calculate Rankings in the DistributionForInstance method and
    //returns them...
    //Here, it would be the other way round, the minIndex would be stored in the prediction,
    //because the highest ordered label in ranking has the lowest value (0).
    
    if(m_ClassIsRanking){    	
   // 	 PreferenceDenseInstance pdi = (PreferenceDenseInstance)instance;
    	
    	pred = Utils.minIndex(dist);
        if (dist[(int)pred] < 0 || dist[(int)pred]>=Double.POSITIVE_INFINITY) {
          pred = Utils.missingValue();
        }
        updateStatsForClassifier(dist, instance);
        if (storePredictions) {
          if (m_Predictions == null)
            m_Predictions = new FastVector();
          m_Predictions.addElement(new NominalPrediction(instance.classValue(), dist,
                                                         instance.weight()));
        }

    }   
    //RANKING END
    if (m_ClassIsNominal && !m_ClassIsRanking) {
      pred = Utils.maxIndex(dist);
      if (dist[(int)pred] <= 0) {
        pred = Utils.missingValue();
      }
      updateStatsForClassifier(dist, instance);
      if (storePredictions) {
        if (m_Predictions == null)
          m_Predictions = new FastVector();
        m_Predictions.addElement(new NominalPrediction(instance.classValue(), dist,
                                                       instance.weight()));
      }
    } else {
      pred = dist[0];
      updateStatsForPredictor(pred, instance);
      if (storePredictions) {
        if (m_Predictions == null)
          m_Predictions = new FastVector();
        m_Predictions.addElement(new NumericPrediction(instance.classValue(), pred,
                                                       instance.weight()));
      }
    }
    return pred;
  }

  /**
   * Evaluates the classifier on a single instance and records the
   * prediction.
   *
   * @param classifier machine learning classifier
   * @param instance the test instance to be classified
   * @param storePredictions whether to store predictions for nominal classifier
   * @return the prediction made by the clasifier
   * @throws Exception if model could not be evaluated
   * successfully or the data contains string attributes
   */
  protected double evaluationForSingleInstance(Classifier classifier,
                                               Instance instance,
                                               boolean storePredictions) throws Exception {
	  //RANKING BEGIN
	  int key=0;
	  //Calculating kendall's Tau, Spearman footrule and 
	  //the relative error matrix for pairwise preferences.
	  train_testInst++;
	  if(instance instanceof PreferenceDenseInstance){
		  PreferenceDenseInstance pdi = (PreferenceDenseInstance) instance;
		  for(int i=0; ;i++){
			  if(pdi.getHashMap().get(i)!=null){
				  RankUtilities.noPartialRanking(pdi.getHashMap().get(i));
				  key=i;
				  break;
			  }
		  }
		  
		  double[] pred = classifier.distributionForInstance(instance);
		  
		  int[] predInt = new int[pred.length];
		  for(int z=0; z<predInt.length; z++){
			  predInt[z]=(int)pred[z];
		  }
		  
		  pdi.insertPrediction(key, predInt);
		  double[] tmp = spearman(RankUtilities.triangleLabels, pred, instance.weight());
		  
		  int[][] out = relativeErrorLabelPairs(pdi.getHashMap().get(key), predInt);
		  
		  for(int a=0; a<rnkDataMatrix.length; a++){
			  for(int b=0; b<rnkDataMatrix[0].length; b++){
				  rnkDataMatrix[a][b]+=out[a][b];
			  }
		  }
		  
		  tau+=kendallsTau(RankUtilities.triangleLabels, pred,instance.weight());		  
		  spearman[0]+= tmp[0];
		  spearman[1]+= tmp[1];
	  }
	//RANKING END 
    
    Instance classMissing = (Instance)instance.copy();
    classMissing.setDataset(instance.dataset());
    
    if (classifier instanceof weka.classifiers.misc.InputMappedClassifier) {
      instance = (Instance)instance.copy();
      instance = 
        ((weka.classifiers.misc.InputMappedClassifier)classifier).
          constructMappedInstance(instance);
//      System.out.println("Mapped instance " + instance);
      int mappedClass = 
        ((weka.classifiers.misc.InputMappedClassifier)classifier).getMappedClassIndex();
      classMissing.setMissing(mappedClass);
    } else {
      classMissing.setClassMissing();
    }
    
//    System.out.println("instance (to predict)" + classMissing);
    double pred = evaluationForSingleInstance(classifier.distributionForInstance(classMissing),
                                              instance, storePredictions);

    // We don't need to do the following if the class is nominal because in that case
    // entropy and coverage statistics are always computed.
    if (!m_ClassIsNominal) {
      if (!instance.classIsMissing() && !Utils.isMissingValue(pred)) {
        if (classifier instanceof IntervalEstimator) {
          updateStatsForIntervalEstimator((IntervalEstimator)classifier, classMissing,
                                          instance.classValue());
        } else {
          m_CoverageStatisticsAvailable = false;
        }
        if (classifier instanceof ConditionalDensityEstimator) {
          updateStatsForConditionalDensityEstimator((ConditionalDensityEstimator)classifier,
                                                    classMissing, instance.classValue());
        } else {
          m_ComplexityStatisticsAvailable = false;
        }
      }
    }
    return pred;
  }

  /**
   * Evaluates the classifier on a single instance and records the
   * prediction.
   *
   * @param classifier machine learning classifier
   * @param instance the test instance to be classified
   * @return the prediction made by the clasifier
   * @throws Exception if model could not be evaluated
   * successfully or the data contains string attributes
   */
  public double evaluateModelOnceAndRecordPrediction(Classifier classifier,
      Instance instance) throws Exception {

    return evaluationForSingleInstance(classifier, instance, true);
  }

  /**
   * Evaluates the classifier on a single instance.
   *
   * @param classifier machine learning classifier
   * @param instance the test instance to be classified
   * @return the prediction made by the clasifier
   * @throws Exception if model could not be evaluated
   * successfully or the data contains string attributes
   */
  public double evaluateModelOnce(Classifier classifier, Instance instance) throws Exception {

    return evaluationForSingleInstance(classifier, instance, false);
  }

  /**
   * Evaluates the supplied distribution on a single instance.
   *
   * @param dist the supplied distribution
   * @param instance the test instance to be classified
   * @return the prediction
   * @throws Exception if model could not be evaluated
   * successfully
   */
  public double evaluateModelOnce(double [] dist, Instance instance) throws Exception {

    return evaluationForSingleInstance(dist, instance, false);
  }

  /**
   * Evaluates the supplied distribution on a single instance.
   *
   * @param dist the supplied distribution
   * @param instance the test instance to be classified
   * @return the prediction
   * @throws Exception if model could not be evaluated
   * successfully
   */
  public double evaluateModelOnceAndRecordPrediction(double [] dist,
      Instance instance) throws Exception {

    return evaluationForSingleInstance(dist, instance, true);
  }

  /**
   * Evaluates the supplied prediction on a single instance.
   *
   * @param prediction the supplied prediction
   * @param instance the test instance to be classified
   * @throws Exception if model could not be evaluated
   * successfully
   */
  public void evaluateModelOnce(double prediction,
      Instance instance) throws Exception {

    evaluateModelOnce(makeDistribution(prediction), instance);
  }

  /**
   * Returns the predictions that have been collected.
   *
   * @return a reference to the FastVector containing the predictions
   * that have been collected. This should be null if no predictions
   * have been collected.
   */
  public FastVector predictions() {
    return m_Predictions;
  }

  /**
   * Wraps a static classifier in enough source to test using the weka
   * class libraries.
   *
   * @param classifier a Sourcable Classifier
   * @param className the name to give to the source code class
   * @return the source for a static classifier that can be tested with
   * weka libraries.
   * @throws Exception if code-generation fails
   */
  public static String wekaStaticWrapper(Sourcable classifier, String className)
    throws Exception {

    StringBuffer result = new StringBuffer();
    String staticClassifier = classifier.toSource(className);

    result.append("// Generated with Weka " + Version.VERSION + "\n");
    result.append("//\n");
    result.append("// This code is public domain and comes with no warranty.\n");
    result.append("//\n");
    result.append("// Timestamp: " + new Date() + "\n");
    result.append("\n");
    result.append("package weka.classifiers;\n");
    result.append("\n");
    result.append("import weka.core.Attribute;\n");
    result.append("import weka.core.Capabilities;\n");
    result.append("import weka.core.Capabilities.Capability;\n");
    result.append("import weka.core.Instance;\n");
    result.append("import weka.core.Instances;\n");
    result.append("import weka.core.RevisionUtils;\n");
    result.append("import weka.classifiers.Classifier;\nimport weka.classifiers.AbstractClassifier;\n");
    result.append("\n");
    result.append("public class WekaWrapper\n");
    result.append("  extends AbstractClassifier {\n");

    // globalInfo
    result.append("\n");
    result.append("  /**\n");
    result.append("   * Returns only the toString() method.\n");
    result.append("   *\n");
    result.append("   * @return a string describing the classifier\n");
    result.append("   */\n");
    result.append("  public String globalInfo() {\n");
    result.append("    return toString();\n");
    result.append("  }\n");

    // getCapabilities
    result.append("\n");
    result.append("  /**\n");
    result.append("   * Returns the capabilities of this classifier.\n");
    result.append("   *\n");
    result.append("   * @return the capabilities\n");
    result.append("   */\n");
    result.append("  public Capabilities getCapabilities() {\n");
    result.append(((Classifier) classifier).getCapabilities().toSource("result", 4));
    result.append("    return result;\n");
    result.append("  }\n");

    // buildClassifier
    result.append("\n");
    result.append("  /**\n");
    result.append("   * only checks the data against its capabilities.\n");
    result.append("   *\n");
    result.append("   * @param i the training data\n");
    result.append("   */\n");
    result.append("  public void buildClassifier(Instances i) throws Exception {\n");
    result.append("    // can classifier handle the data?\n");
    result.append("    getCapabilities().testWithFail(i);\n");
    result.append("  }\n");

    // classifyInstance
    result.append("\n");
    result.append("  /**\n");
    result.append("   * Classifies the given instance.\n");
    result.append("   *\n");
    result.append("   * @param i the instance to classify\n");
    result.append("   * @return the classification result\n");
    result.append("   */\n");
    result.append("  public double classifyInstance(Instance i) throws Exception {\n");
    result.append("    Object[] s = new Object[i.numAttributes()];\n");
    result.append("    \n");
    result.append("    for (int j = 0; j < s.length; j++) {\n");
    result.append("      if (!i.isMissing(j)) {\n");
    result.append("        if (i.attribute(j).isNominal() || i.attribute(j).isRanking())\n");
    result.append("          s[j] = new String(i.stringValue(j));\n");
    result.append("        else if (i.attribute(j).isNumeric())\n");
    result.append("          s[j] = new Double(i.value(j));\n");
    result.append("      }\n");
    result.append("    }\n");
    result.append("    \n");
    result.append("    // set class value to missing\n");
    result.append("    s[i.classIndex()] = null;\n");
    result.append("    \n");
    result.append("    return " + className + ".classify(s);\n");
    result.append("  }\n");

    // getRevision
    result.append("\n");
    result.append("  /**\n");
    result.append("   * Returns the revision string.\n");
    result.append("   * \n");
    result.append("   * @return        the revision\n");
    result.append("   */\n");
    result.append("  public String getRevision() {\n");
    result.append("    return RevisionUtils.extract(\"1.0\");\n");
    result.append("  }\n");

    // toString
    result.append("\n");
    result.append("  /**\n");
    result.append("   * Returns only the classnames and what classifier it is based on.\n");
    result.append("   *\n");
    result.append("   * @return a short description\n");
    result.append("   */\n");
    result.append("  public String toString() {\n");
    result.append("    return \"Auto-generated classifier wrapper, based on "
        + classifier.getClass().getName() + " (generated with Weka " + Version.VERSION + ").\\n"
        + "\" + this.getClass().getName() + \"/" + className + "\";\n");
    result.append("  }\n");

    // main
    result.append("\n");
    result.append("  /**\n");
    result.append("   * Runs the classfier from commandline.\n");
    result.append("   *\n");
    result.append("   * @param args the commandline arguments\n");
    result.append("   */\n");
    result.append("  public static void main(String args[]) {\n");
    result.append("    runClassifier(new WekaWrapper(), args);\n");
    result.append("  }\n");
    result.append("}\n");

    // actual classifier code
    result.append("\n");
    result.append(staticClassifier);

    return result.toString();
  }

  /**
   * Gets the number of test instances that had a known class value
   * (actually the sum of the weights of test instances with known
   * class value).
   *
   * @return the number of test instances with known class
   */
  public final double numInstances() {

    return m_WithClass;
  }

  /**
   * Gets the coverage of the test cases by the predicted regions at
   * the confidence level specified when evaluation was performed.
   *
   * @return the coverage of the test cases by the predicted regions
   */
  public final double coverageOfTestCasesByPredictedRegions() {

    if (!m_CoverageStatisticsAvailable)
      return Double.NaN;

    return 100 * m_TotalCoverage / m_WithClass;
  }

  /**
   * Gets the average size of the predicted regions, relative to the
   * range of the target in the training data, at the confidence level
   * specified when evaluation was performed.
   *
   * @return the average size of the predicted regions
   */
  public final double sizeOfPredictedRegions() {

    if (m_NoPriors || !m_CoverageStatisticsAvailable)
      return Double.NaN;

    return 100 * m_TotalSizeOfRegions / m_WithClass;
  }

  /**
   * Gets the number of instances incorrectly classified (that is, for
   * which an incorrect prediction was made). (Actually the sum of the
   * weights of these instances)
   *
   * @return the number of incorrectly classified instances
   */
  public final double incorrect() {

    return m_Incorrect;
  }

  /**
   * Gets the percentage of instances incorrectly classified (that is,
   * for which an incorrect prediction was made).
   *
   * @return the percent of incorrectly classified instances
   * (between 0 and 100)
   */
  public final double pctIncorrect() {

    return 100 * m_Incorrect / m_WithClass;
  }

  /**
   * Gets the total cost, that is, the cost of each prediction times the
   * weight of the instance, summed over all instances.
   *
   * @return the total cost
   */
  public final double totalCost() {

    return m_TotalCost;
  }

  /**
   * Gets the average cost, that is, total cost of misclassifications
   * (incorrect plus unclassified) over the total number of instances.
   *
   * @return the average cost.
   */
  public final double avgCost() {

    return m_TotalCost / m_WithClass;
  }

  /**
   * Gets the number of instances correctly classified (that is, for
   * which a correct prediction was made). (Actually the sum of the weights
   * of these instances)
   *
   * @return the number of correctly classified instances
   */
  public final double correct() {

    return m_Correct;
  }

  /**
   * Gets the percentage of instances correctly classified (that is, for
   * which a correct prediction was made).
   *
   * @return the percent of correctly classified instances (between 0 and 100)
   */
  public final double pctCorrect() {

    return 100 * m_Correct / m_WithClass;
  }

  /**
   * Gets the number of instances not classified (that is, for
   * which no prediction was made by the classifier). (Actually the sum
   * of the weights of these instances)
   *
   * @return the number of unclassified instances
   */
  public final double unclassified() {

    return m_Unclassified;
  }

  /**
   * Gets the percentage of instances not classified (that is, for
   * which no prediction was made by the classifier).
   *
   * @return the percent of unclassified instances (between 0 and 100)
   */
  public final double pctUnclassified() {

    return 100 * m_Unclassified / m_WithClass;
  }

  /**
   * Returns the estimated error rate or the root mean squared error
   * (if the class is numeric). If a cost matrix was given this
   * error rate gives the average cost.
   *
   * @return the estimated error rate (between 0 and 1, or between 0 and
   * maximum cost)
   */
  public final double errorRate() {

    if (!m_ClassIsNominal) {
      return Math.sqrt(m_SumSqrErr / (m_WithClass - m_Unclassified));
    }
    if (m_CostMatrix == null) {
      return m_Incorrect / m_WithClass;
    } else {
      return avgCost();
    }
  }

  /**
   * Returns value of kappa statistic if class is nominal.
   *
   * @return the value of the kappa statistic
   */
  public final double kappa() {


    double[] sumRows = new double[m_ConfusionMatrix.length];
    double[] sumColumns = new double[m_ConfusionMatrix.length];
    double sumOfWeights = 0;
    for (int i = 0; i < m_ConfusionMatrix.length; i++) {
      for (int j = 0; j < m_ConfusionMatrix.length; j++) {
        sumRows[i] += m_ConfusionMatrix[i][j];
        sumColumns[j] += m_ConfusionMatrix[i][j];
        sumOfWeights += m_ConfusionMatrix[i][j];
      }
    }
    double correct = 0, chanceAgreement = 0;
    for (int i = 0; i < m_ConfusionMatrix.length; i++) {
      chanceAgreement += (sumRows[i] * sumColumns[i]);
      correct += m_ConfusionMatrix[i][i];
    }
    chanceAgreement /= (sumOfWeights * sumOfWeights);
    correct /= sumOfWeights;

    if (chanceAgreement < 1) {
      return (correct - chanceAgreement) / (1 - chanceAgreement);
    } else {
      return 1;
    }
  }

  /**
   * Returns the correlation coefficient if the class is numeric.
   *
   * @return the correlation coefficient
   * @throws Exception if class is not numeric
   */
  public final double correlationCoefficient() throws Exception {

    if (m_ClassIsNominal) {
      throw
      new Exception("Can't compute correlation coefficient: " +
      "class is nominal!");
    }

    double correlation = 0;
    double varActual =
      m_SumSqrClass - m_SumClass * m_SumClass /
      (m_WithClass - m_Unclassified);
    double varPredicted =
      m_SumSqrPredicted - m_SumPredicted * m_SumPredicted /
      (m_WithClass - m_Unclassified);
    double varProd =
      m_SumClassPredicted - m_SumClass * m_SumPredicted /
      (m_WithClass - m_Unclassified);

    if (varActual * varPredicted <= 0) {
      correlation = 0.0;
    } else {
      correlation = varProd / Math.sqrt(varActual * varPredicted);
    }

    return correlation;
  }

  /**
   * Returns the mean absolute error. Refers to the error of the
   * predicted values for numeric classes, and the error of the
   * predicted probability distribution for nominal classes.
   *
   * @return the mean absolute error
   */
  public final double meanAbsoluteError() {

    return m_SumAbsErr / (m_WithClass - m_Unclassified);
  }

  /**
   * Returns the mean absolute error of the prior.
   *
   * @return the mean absolute error
   */
  public final double meanPriorAbsoluteError() {

    if (m_NoPriors)
      return Double.NaN;

    return m_SumPriorAbsErr / m_WithClass;
  }

  /**
   * Returns the relative absolute error.
   *
   * @return the relative absolute error
   * @throws Exception if it can't be computed
   */
  public final double relativeAbsoluteError() throws Exception {

    if (m_NoPriors)
      return Double.NaN;

    return 100 * meanAbsoluteError() / meanPriorAbsoluteError();
  }

  /**
   * Returns the root mean squared error.
   *
   * @return the root mean squared error
   */
  public final double rootMeanSquaredError() {

    return Math.sqrt(m_SumSqrErr / (m_WithClass - m_Unclassified));
  }

  /**
   * Returns the root mean prior squared error.
   *
   * @return the root mean prior squared error
   */
  public final double rootMeanPriorSquaredError() {

    if (m_NoPriors)
      return Double.NaN;

    return Math.sqrt(m_SumPriorSqrErr / m_WithClass);
  }

  /**
   * Returns the root relative squared error if the class is numeric.
   *
   * @return the root relative squared error
   */
  public final double rootRelativeSquaredError() {

    if (m_NoPriors)
      return Double.NaN;

    return 100.0 * rootMeanSquaredError() / rootMeanPriorSquaredError();
  }

  /**
   * Calculate the entropy of the prior distribution.
   *
   * @return the entropy of the prior distribution
   * @throws Exception if the class is not nominal
   */
  public final double priorEntropy() throws Exception {

    if (!m_ClassIsNominal) {
      throw
      new Exception("Can't compute entropy of class prior: " +
      "class numeric!");
    }

    if (m_NoPriors)
      return Double.NaN;

    double entropy = 0;
    for(int i = 0; i < m_NumClasses; i++) {
      entropy -= m_ClassPriors[i] / m_ClassPriorsSum *
        Utils.log2(m_ClassPriors[i] / m_ClassPriorsSum);
    }
    return entropy;
  }

  /**
   * Return the total Kononenko & Bratko Information score in bits.
   *
   * @return the K&B information score
   * @throws Exception if the class is not nominal
   */
  public final double KBInformation() throws Exception {

    if (!m_ClassIsNominal) {
      throw
      new Exception("Can't compute K&B Info score: " +
      "class numeric!");
    }

    if (m_NoPriors)
      return Double.NaN;

    return m_SumKBInfo;
  }

  /**
   * Return the Kononenko & Bratko Information score in bits per
   * instance.
   *
   * @return the K&B information score
   * @throws Exception if the class is not nominal
   */
  public final double KBMeanInformation() throws Exception {

    if (!m_ClassIsNominal) {
      throw
      new Exception("Can't compute K&B Info score: class numeric!");
    }

    if (m_NoPriors)
      return Double.NaN;

    return m_SumKBInfo / (m_WithClass - m_Unclassified);
  }

  /**
   * Return the Kononenko & Bratko Relative Information score.
   *
   * @return the K&B relative information score
   * @throws Exception if the class is not nominal
   */
  public final double KBRelativeInformation() throws Exception {

    if (!m_ClassIsNominal) {
      throw
      new Exception("Can't compute K&B Info score: " +
      "class numeric!");
    }

    if (m_NoPriors)
      return Double.NaN;

    return 100.0 * KBInformation() / priorEntropy();
  }

  /**
   * Returns the total entropy for the null model.
   *
   * @return the total null model entropy
   */
  public final double SFPriorEntropy() {

    if (m_NoPriors || !m_ComplexityStatisticsAvailable)
      return Double.NaN;

    return m_SumPriorEntropy;
  }

  /**
   * Returns the entropy per instance for the null model.
   *
   * @return the null model entropy per instance
   */
  public final double SFMeanPriorEntropy() {

    if (m_NoPriors || !m_ComplexityStatisticsAvailable)
      return Double.NaN;

    return m_SumPriorEntropy / m_WithClass;
  }

  /**
   * Returns the total entropy for the scheme.
   *
   * @return the total scheme entropy
   */
  public final double SFSchemeEntropy() {

    if (!m_ComplexityStatisticsAvailable)
      return Double.NaN;

    return m_SumSchemeEntropy;
  }

  /**
   * Returns the entropy per instance for the scheme.
   *
   * @return the scheme entropy per instance
   */
  public final double SFMeanSchemeEntropy() {

    if (!m_ComplexityStatisticsAvailable)
      return Double.NaN;

    return m_SumSchemeEntropy / (m_WithClass - m_Unclassified);
  }

  /**
   * Returns the total SF, which is the null model entropy minus
   * the scheme entropy.
   *
   * @return the total SF
   */
  public final double SFEntropyGain() {

    if (m_NoPriors || !m_ComplexityStatisticsAvailable)
      return Double.NaN;

    return m_SumPriorEntropy - m_SumSchemeEntropy;
  }

  /**
   * Returns the SF per instance, which is the null model entropy
   * minus the scheme entropy, per instance.
   *
   * @return the SF per instance
   */
  public final double SFMeanEntropyGain() {

    if (m_NoPriors || !m_ComplexityStatisticsAvailable)
      return Double.NaN;

    return (m_SumPriorEntropy - m_SumSchemeEntropy) /
      (m_WithClass - m_Unclassified);
  }

  /**
   * Output the cumulative margin distribution as a string suitable
   * for input for gnuplot or similar package.
   *
   * @return the cumulative margin distribution
   * @throws Exception if the class attribute is nominal
   */
  public String toCumulativeMarginDistributionString() throws Exception {

    if (!m_ClassIsNominal) {
      throw new Exception("Class must be nominal for margin distributions");
    }
    String result = "";
    double cumulativeCount = 0;
    double margin;
    for(int i = 0; i <= k_MarginResolution; i++) {
      if (m_MarginCounts[i] != 0) {
        cumulativeCount += m_MarginCounts[i];
        margin = (double)i * 2.0 / k_MarginResolution - 1.0;
        result = result + Utils.doubleToString(margin, 7, 3) + ' '
          + Utils.doubleToString(cumulativeCount * 100
              / m_WithClass, 7, 3) + '\n';
      } else if (i == 0) {
        result = Utils.doubleToString(-1.0, 7, 3) + ' '
          + Utils.doubleToString(0, 7, 3) + '\n';
      }
    }
    return result;
  }

  /**
   * Calls toSummaryString() with no title and no complexity stats.
   *
   * @return a summary description of the classifier evaluation
   */
  public String toSummaryString() {

    return toSummaryString("", false);
  }

  /**
   * Calls toSummaryString() with a default title.
   *
   * @param printComplexityStatistics if true, complexity statistics are
   * returned as well
   * @return the summary string
   */
  public String toSummaryString(boolean printComplexityStatistics) {

    return toSummaryString("=== Summary ===\n", printComplexityStatistics);
  }

  /**
   * Outputs the performance statistics in summary form. Lists
   * number (and percentage) of instances classified correctly,
   * incorrectly and unclassified. Outputs the total number of
   * instances classified, and the number of instances (if any)
   * that had no class value provided.
   *
   * @param title the title for the statistics
   * @param printComplexityStatistics if true, complexity statistics are
   * returned as well
   * @return the summary as a String
   */
  public String toSummaryString(String title,
      boolean printComplexityStatistics) {

    StringBuffer text = new StringBuffer();

    if (printComplexityStatistics && m_NoPriors) {
      printComplexityStatistics = false;
      System.err.println("Priors disabled, cannot print complexity statistics!");
    }

    text.append(title + "\n");
    try {
      if (m_WithClass > 0) {
        if (m_ClassIsNominal) {

          text.append("Correctly Classified Instances     ");
          text.append(Utils.doubleToString(correct(), 12, 4) + "     " +
              Utils.doubleToString(pctCorrect(),
                12, 4) + " %\n");
          text.append("Incorrectly Classified Instances   ");
          text.append(Utils.doubleToString(incorrect(), 12, 4) + "     " +
              Utils.doubleToString(pctIncorrect(),
                12, 4) + " %\n");
          text.append("Kappa statistic                    ");
          text.append(Utils.doubleToString(kappa(), 12, 4) + "\n");

          if (m_CostMatrix != null) {
            text.append("Total Cost                         ");
            text.append(Utils.doubleToString(totalCost(), 12, 4) + "\n");
            text.append("Average Cost                       ");
            text.append(Utils.doubleToString(avgCost(), 12, 4) + "\n");
          }
          if (printComplexityStatistics) {
            text.append("K&B Relative Info Score            ");
            text.append(Utils.doubleToString(KBRelativeInformation(), 12, 4)
                + " %\n");
            text.append("K&B Information Score              ");
            text.append(Utils.doubleToString(KBInformation(), 12, 4)
                + " bits");
            text.append(Utils.doubleToString(KBMeanInformation(), 12, 4)
                + " bits/instance\n");
          }
        } else {
          text.append("Correlation coefficient            ");
          text.append(Utils.doubleToString(correlationCoefficient(), 12 , 4) +
              "\n");
        }
        if (printComplexityStatistics && m_ComplexityStatisticsAvailable) {
          text.append("Class complexity | order 0         ");
          text.append(Utils.doubleToString(SFPriorEntropy(), 12, 4)
              + " bits");
          text.append(Utils.doubleToString(SFMeanPriorEntropy(), 12, 4)
              + " bits/instance\n");
          text.append("Class complexity | scheme          ");
          text.append(Utils.doubleToString(SFSchemeEntropy(), 12, 4)
              + " bits");
          text.append(Utils.doubleToString(SFMeanSchemeEntropy(), 12, 4)
              + " bits/instance\n");
          text.append("Complexity improvement     (Sf)    ");
          text.append(Utils.doubleToString(SFEntropyGain(), 12, 4) + " bits");
          text.append(Utils.doubleToString(SFMeanEntropyGain(), 12, 4)
              + " bits/instance\n");
        }

        text.append("Mean absolute error                ");
        text.append(Utils.doubleToString(meanAbsoluteError(), 12, 4)
            + "\n");
        text.append("Root mean squared error            ");
        text.append(Utils.
            doubleToString(rootMeanSquaredError(), 12, 4)
            + "\n");
        if (!m_NoPriors) {
          text.append("Relative absolute error            ");
          text.append(Utils.doubleToString(relativeAbsoluteError(),
                12, 4) + " %\n");
          text.append("Root relative squared error        ");
          text.append(Utils.doubleToString(rootRelativeSquaredError(),
                12, 4) + " %\n");
        }
        if (m_CoverageStatisticsAvailable) {
          text.append("Coverage of cases (" + Utils.doubleToString(m_ConfLevel, 4, 2) + " level)     ");
          text.append(Utils.doubleToString(coverageOfTestCasesByPredictedRegions(),
                12, 4) + " %\n");
          if (!m_NoPriors) {
            text.append("Mean rel. region size (" + Utils.doubleToString(m_ConfLevel, 4, 2) + " level) ");
            text.append(Utils.doubleToString(sizeOfPredictedRegions(), 12, 4) + " %\n");
          }
        }
      }
      if (Utils.gr(unclassified(), 0)) {
        text.append("UnClassified Instances             ");
        text.append(Utils.doubleToString(unclassified(), 12,4) +  "     " +
            Utils.doubleToString(pctUnclassified(),
              12, 4) + " %\n");
      }
      text.append("Total Number of Instances          ");
    	  text.append(Utils.doubleToString(m_WithClass, 12, 4) + "\n");
      if (m_MissingClass > 0) {
        text.append("Ignored Class Unknown Instances            ");
        text.append(Utils.doubleToString(m_MissingClass, 12, 4) + "\n");
      }
    } catch (Exception ex) {
      // Should never occur since the class is known to be nominal
      // here
      System.err.println("Arggh - Must be a bug in Evaluation class");
    }

    return text.toString();
  }

  /**
   * Calls toMatrixString() with a default title.
   *
   * @return the confusion matrix as a string
   * @throws Exception if the class is numeric
   */
  public String toMatrixString() throws Exception {

    return toMatrixString("=== Confusion Matrix ===\n");
  }

  /**
   * Outputs the performance statistics as a classification confusion
   * matrix. For each class value, shows the distribution of
   * predicted class values.
   *
   * @param title the title for the confusion matrix
   * @return the confusion matrix as a String
   * @throws Exception if the class is numeric
   */
  public String toMatrixString(String title) throws Exception {

    StringBuffer text = new StringBuffer();
    char [] IDChars = {'a','b','c','d','e','f','g','h','i','j',
      'k','l','m','n','o','p','q','r','s','t',
      'u','v','w','x','y','z'};
    int IDWidth;
    boolean fractional = false;

    if (!m_ClassIsNominal) {
      throw new Exception("Evaluation: No confusion matrix possible!");
    }

    // Find the maximum value in the matrix
    // and check for fractional display requirement
    double maxval = 0;
//    for(int i=0; i < LabelList.rankInstances.size(); i++){
//    	for(int j=0; j < LabelList.rankInstances.size(); j++){   
//    for(int i=0; i < RankUtilities.indexedInstances.size(); i++){
//    	for(int j=0; j < RankUtilities.indexedInstances.size(); j++){   
    for(int i = 0; i < m_NumClasses; i++) {
      for(int j = 0; j < m_NumClasses; j++) {
        double current = m_ConfusionMatrix[i][j];
        if (current < 0) {
          current *= -10;
        }
        if (current > maxval) {
          maxval = current;
        }
        double fract = current - Math.rint(current);
        if (!fractional && ((Math.log(fract) / Math.log(10)) >= -2)) {
          fractional = true;
        }
      }
    }

    IDWidth = 1 + Math.max((int)(Math.log(maxval) / Math.log(10)
          + (fractional ? 3 : 0)),
        (int)(Math.log(m_NumClasses) /
          Math.log(IDChars.length)));
    text.append(title).append("\n");
//    for(int i=0; i < RankUtilities.indexedInstances.size(); i++){
    for(int i = 0; i < m_NumClasses; i++) {
      if (fractional) {
        text.append(" ").append(num2ShortID(i,IDChars,IDWidth - 3))
          .append("   ");
      } else {
        text.append(" ").append(num2ShortID(i,IDChars,IDWidth));
      }
    }
    text.append("   <-- classified as\n");
//  for(int i=0; i < LabelList.rankInstances.size(); i++){
//	for(int j=0; j < LabelList.rankInstances.size(); j++){   
//    for(int i=0; i < RankUtilities.indexedInstances.size(); i++){
//    	for(int j=0; j < RankUtilities.indexedInstances.size(); j++){   
    for(int i = 0; i< m_NumClasses; i++) {
      for(int j = 0; j < m_NumClasses; j++) {

        text.append(" ").append(
            Utils.doubleToString(m_ConfusionMatrix[i][j],
              IDWidth,
              (fractional ? 2 : 0)));
      }
      text.append(" | ").append(num2ShortID(i,IDChars,IDWidth))
        .append(" = ").append(m_ClassNames[i]).append("\n");
    }
    
    
    
    return text.toString();
  }

  /**
   * Generates a breakdown of the accuracy for each class (with default title),
   * incorporating various information-retrieval statistics, such as
   * true/false positive rate, precision/recall/F-Measure.  Should be
   * useful for ROC curves, recall/precision curves.
   *
   * @return the statistics presented as a string
   * @throws Exception if class is not nominal
   */
  public String toClassDetailsString() throws Exception {

    return toClassDetailsString("=== Detailed Accuracy By Class ===\n");
  }

  /**
   * Generates a breakdown of the accuracy for each class,
   * incorporating various information-retrieval statistics, such as
   * true/false positive rate, precision/recall/F-Measure.  Should be
   * useful for ROC curves, recall/precision curves.
   *
   * @param title the title to prepend the stats string with
   * @return the statistics presented as a string
   * @throws Exception if class is not nominal
   */
  public String toClassDetailsString(String title) throws Exception {

    if (!m_ClassIsNominal) {
      throw new Exception("Evaluation: No per class statistics possible!");
    }

    StringBuffer text = new StringBuffer(title
        + "\n               TP Rate   FP Rate"
        + "   Precision   Recall"
        + "  F-Measure   ROC Area  Class\n");
    for(int i = 0; i < m_NumClasses; i++) {
      text.append("               " + Utils.doubleToString(truePositiveRate(i), 7, 3))
        .append("   ");
      text.append(Utils.doubleToString(falsePositiveRate(i), 7, 3))
        .append("    ");
      text.append(Utils.doubleToString(precision(i), 7, 3))
        .append("   ");
      text.append(Utils.doubleToString(recall(i), 7, 3))
        .append("   ");
      text.append(Utils.doubleToString(fMeasure(i), 7, 3))
        .append("    ");

      double rocVal = areaUnderROC(i);
      if (Utils.isMissingValue(rocVal)) {
        text.append("  ?    ")
          .append("    ");
      } else {
        text.append(Utils.doubleToString(rocVal, 7, 3))
          .append("    ");
      }
      text.append(m_ClassNames[i]).append('\n');
    }

    text.append("Weighted Avg.  " + Utils.doubleToString(weightedTruePositiveRate(), 7, 3));
    text.append("   " + Utils.doubleToString(weightedFalsePositiveRate(), 7 ,3));
    text.append("    " + Utils.doubleToString(weightedPrecision(), 7 ,3));
    text.append("   " + Utils.doubleToString(weightedRecall(), 7 ,3));
    text.append("   " + Utils.doubleToString(weightedFMeasure(), 7 ,3));
    text.append("    " + Utils.doubleToString(weightedAreaUnderROC(), 7 ,3));
    text.append("\n");

    return text.toString();
  }

  /**
   * Calculate the number of true positives with respect to a particular class.
   * This is defined as<p/>
   * <pre>
   * correctly classified positives
   * </pre>
   *
   * @param classIndex the index of the class to consider as "positive"
   * @return the true positive rate
   */
  public double numTruePositives(int classIndex) {

    double correct = 0;
    for (int j = 0; j < m_NumClasses; j++) {
      if (j == classIndex) {
        correct += m_ConfusionMatrix[classIndex][j];
      }
    }
    return correct;
  }

  /**
   * Calculate the true positive rate with respect to a particular class.
   * This is defined as<p/>
   * <pre>
   * correctly classified positives
   * ------------------------------
   *       total positives
   * </pre>
   *
   * @param classIndex the index of the class to consider as "positive"
   * @return the true positive rate
   */
  public double truePositiveRate(int classIndex) {

    double correct = 0, total = 0;
    for (int j = 0; j < m_NumClasses; j++) {
      if (j == classIndex) {
        correct += m_ConfusionMatrix[classIndex][j];
      }
      total += m_ConfusionMatrix[classIndex][j];
    }
    if (total == 0) {
      return 0;
    }
    return correct / total;
  }

  /**
   * Calculates the weighted (by class size) true positive rate.
   *
   * @return the weighted true positive rate.
   */
  public double weightedTruePositiveRate() {
    double[] classCounts = new double[m_NumClasses];
    double classCountSum = 0;

    for (int i = 0; i < m_NumClasses; i++) {
      for (int j = 0; j < m_NumClasses; j++) {
        classCounts[i] += m_ConfusionMatrix[i][j];
      }
      classCountSum += classCounts[i];
    }

    double truePosTotal = 0;
    for(int i = 0; i < m_NumClasses; i++) {
      double temp = truePositiveRate(i);
      truePosTotal += (temp * classCounts[i]);
    }

    return truePosTotal / classCountSum;
  }

  /**
   * Calculate the number of true negatives with respect to a particular class.
   * This is defined as<p/>
   * <pre>
   * correctly classified negatives
   * </pre>
   *
   * @param classIndex the index of the class to consider as "positive"
   * @return the true positive rate
   */
  public double numTrueNegatives(int classIndex) {

    double correct = 0;
    for (int i = 0; i < m_NumClasses; i++) {
      if (i != classIndex) {
        for (int j = 0; j < m_NumClasses; j++) {
          if (j != classIndex) {
            correct += m_ConfusionMatrix[i][j];
          }
        }
      }
    }
    return correct;
  }

  /**
   * Calculate the true negative rate with respect to a particular class.
   * This is defined as<p/>
   * <pre>
   * correctly classified negatives
   * ------------------------------
   *       total negatives
   * </pre>
   *
   * @param classIndex the index of the class to consider as "positive"
   * @return the true positive rate
   */
  public double trueNegativeRate(int classIndex) {

    double correct = 0, total = 0;
    for (int i = 0; i < m_NumClasses; i++) {
      if (i != classIndex) {
        for (int j = 0; j < m_NumClasses; j++) {
          if (j != classIndex) {
            correct += m_ConfusionMatrix[i][j];
          }
          total += m_ConfusionMatrix[i][j];
        }
      }
    }
    if (total == 0) {
      return 0;
    }
    return correct / total;
  }

  /**
   * Calculates the weighted (by class size) true negative rate.
   *
   * @return the weighted true negative rate.
   */
  public double weightedTrueNegativeRate() {
    double[] classCounts = new double[m_NumClasses];
    double classCountSum = 0;

    for (int i = 0; i < m_NumClasses; i++) {
      for (int j = 0; j < m_NumClasses; j++) {
        classCounts[i] += m_ConfusionMatrix[i][j];
      }
      classCountSum += classCounts[i];
    }

    double trueNegTotal = 0;
    for(int i = 0; i < m_NumClasses; i++) {
      double temp = trueNegativeRate(i);
      trueNegTotal += (temp * classCounts[i]);
    }

    return trueNegTotal / classCountSum;
  }

  /**
   * Calculate number of false positives with respect to a particular class.
   * This is defined as<p/>
   * <pre>
   * incorrectly classified negatives
   * </pre>
   *
   * @param classIndex the index of the class to consider as "positive"
   * @return the false positive rate
   */
  public double numFalsePositives(int classIndex) {

    double incorrect = 0;
    for (int i = 0; i < m_NumClasses; i++) {
      if (i != classIndex) {
        for (int j = 0; j < m_NumClasses; j++) {
          if (j == classIndex) {
            incorrect += m_ConfusionMatrix[i][j];
          }
        }
      }
    }
    return incorrect;
  }

  /**
   * Calculate the false positive rate with respect to a particular class.
   * This is defined as<p/>
   * <pre>
   * incorrectly classified negatives
   * --------------------------------
   *        total negatives
   * </pre>
   *
   * @param classIndex the index of the class to consider as "positive"
   * @return the false positive rate
   */
  public double falsePositiveRate(int classIndex) {

    double incorrect = 0, total = 0;
    for (int i = 0; i < m_NumClasses; i++) {
      if (i != classIndex) {
        for (int j = 0; j < m_NumClasses; j++) {
          if (j == classIndex) {
            incorrect += m_ConfusionMatrix[i][j];
          }
          total += m_ConfusionMatrix[i][j];
        }
      }
    }
    if (total == 0) {
      return 0;
    }
    return incorrect / total;
  }

  /**
   * Calculates the weighted (by class size) false positive rate.
   *
   * @return the weighted false positive rate.
   */
  public double weightedFalsePositiveRate() {
    double[] classCounts = new double[m_NumClasses];
    double classCountSum = 0;

    for (int i = 0; i < m_NumClasses; i++) {
      for (int j = 0; j < m_NumClasses; j++) {
        classCounts[i] += m_ConfusionMatrix[i][j];
      }
      classCountSum += classCounts[i];
    }

    double falsePosTotal = 0;
    for(int i = 0; i < m_NumClasses; i++) {
      double temp = falsePositiveRate(i);
      falsePosTotal += (temp * classCounts[i]);
    }

    return falsePosTotal / classCountSum;
  }



  /**
   * Calculate number of false negatives with respect to a particular class.
   * This is defined as<p/>
   * <pre>
   * incorrectly classified positives
   * </pre>
   *
   * @param classIndex the index of the class to consider as "positive"
   * @return the false positive rate
   */
  public double numFalseNegatives(int classIndex) {

    double incorrect = 0;
    for (int i = 0; i < m_NumClasses; i++) {
      if (i == classIndex) {
        for (int j = 0; j < m_NumClasses; j++) {
          if (j != classIndex) {
            incorrect += m_ConfusionMatrix[i][j];
          }
        }
      }
    }
    return incorrect;
  }

  /**
   * Calculate the false negative rate with respect to a particular class.
   * This is defined as<p/>
   * <pre>
   * incorrectly classified positives
   * --------------------------------
   *        total positives
   * </pre>
   *
   * @param classIndex the index of the class to consider as "positive"
   * @return the false positive rate
   */
  public double falseNegativeRate(int classIndex) {

    double incorrect = 0, total = 0;
    for (int i = 0; i < m_NumClasses; i++) {
      if (i == classIndex) {
        for (int j = 0; j < m_NumClasses; j++) {
          if (j != classIndex) {
            incorrect += m_ConfusionMatrix[i][j];
          }
          total += m_ConfusionMatrix[i][j];
        }
      }
    }
    if (total == 0) {
      return 0;
    }
    return incorrect / total;
  }

  /**
   * Calculates the weighted (by class size) false negative rate.
   *
   * @return the weighted false negative rate.
   */
  public double weightedFalseNegativeRate() {
    double[] classCounts = new double[m_NumClasses];
    double classCountSum = 0;

    for (int i = 0; i < m_NumClasses; i++) {
      for (int j = 0; j < m_NumClasses; j++) {
        classCounts[i] += m_ConfusionMatrix[i][j];
      }
      classCountSum += classCounts[i];
    }

    double falseNegTotal = 0;
    for(int i = 0; i < m_NumClasses; i++) {
      double temp = falseNegativeRate(i);
      falseNegTotal += (temp * classCounts[i]);
    }

    return falseNegTotal / classCountSum;
  }

  /**
   * Calculate the recall with respect to a particular class.
   * This is defined as<p/>
   * <pre>
   * correctly classified positives
   * ------------------------------
   *       total positives
   * </pre><p/>
   * (Which is also the same as the truePositiveRate.)
   *
   * @param classIndex the index of the class to consider as "positive"
   * @return the recall
   */
  public double recall(int classIndex) {

    return truePositiveRate(classIndex);
  }

  /**
   * Calculates the weighted (by class size) recall.
   *
   * @return the weighted recall.
   */
  public double weightedRecall() {
    return weightedTruePositiveRate();
  }

  /**
   * Calculate the precision with respect to a particular class.
   * This is defined as<p/>
   * <pre>
   * correctly classified positives
   * ------------------------------
   *  total predicted as positive
   * </pre>
   *
   * @param classIndex the index of the class to consider as "positive"
   * @return the precision
   */
  public double precision(int classIndex) {

    double correct = 0, total = 0;
    for (int i = 0; i < m_NumClasses; i++) {
      if (i == classIndex) {
        correct += m_ConfusionMatrix[i][classIndex];
      }
      total += m_ConfusionMatrix[i][classIndex];
    }
    if (total == 0) {
      return 0;
    }
    return correct / total;
  }

  /**
   * Calculates the weighted (by class size) false precision.
   *
   * @return the weighted precision.
   */
  public double weightedPrecision() {
    double[] classCounts = new double[m_NumClasses];
    double classCountSum = 0;

    for (int i = 0; i < m_NumClasses; i++) {
      for (int j = 0; j < m_NumClasses; j++) {
        classCounts[i] += m_ConfusionMatrix[i][j];
      }
      classCountSum += classCounts[i];
    }

    double precisionTotal = 0;
    for(int i = 0; i < m_NumClasses; i++) {
      double temp = precision(i);
      precisionTotal += (temp * classCounts[i]);
    }

    return precisionTotal / classCountSum;
  }

  /**
   * Calculate the F-Measure with respect to a particular class.
   * This is defined as<p/>
   * <pre>
   * 2 * recall * precision
   * ----------------------
   *   recall + precision
   * </pre>
   *
   * @param classIndex the index of the class to consider as "positive"
   * @return the F-Measure
   */
  public double fMeasure(int classIndex) {

    double precision = precision(classIndex);
    double recall = recall(classIndex);
    if ((precision + recall) == 0) {
      return 0;
    }
    return 2 * precision * recall / (precision + recall);
  }

  /**
   * Calculates the macro weighted (by class size) average
   * F-Measure.
   *
   * @return the weighted F-Measure.
   */
  public double weightedFMeasure() {
    double[] classCounts = new double[m_NumClasses];
    double classCountSum = 0;

    for (int i = 0; i < m_NumClasses; i++) {
      for (int j = 0; j < m_NumClasses; j++) {
        classCounts[i] += m_ConfusionMatrix[i][j];
      }
      classCountSum += classCounts[i];
    }

    double fMeasureTotal = 0;
    for(int i = 0; i < m_NumClasses; i++) {
      double temp = fMeasure(i);
      fMeasureTotal += (temp * classCounts[i]);
    }

    return fMeasureTotal / classCountSum;
  }

  /**
   * Unweighted macro-averaged F-measure. If some classes not present in the
   * test set, they're just skipped (since recall is undefined there anyway) .
   *
   * @return unweighted macro-averaged F-measure.
   * */
  public double unweightedMacroFmeasure() {
    weka.experiment.Stats rr = new weka.experiment.Stats();
    for (int c = 0; c < m_NumClasses; c++) {
      // skip if no testing positive cases of this class
      if (numTruePositives(c)+numFalseNegatives(c) > 0) {
        rr.add(fMeasure(c));
      }
    }
    rr.calculateDerived();
    return rr.mean;
  }

  /**
   * Unweighted micro-averaged F-measure. If some classes not present in the
   * test set, they have no effect.
   *
   * Note: if the test set is *single-label*, then this is the same as accuracy.
   *
   * @return unweighted micro-averaged F-measure.
   */
  public double unweightedMicroFmeasure() {
    double tp = 0;
    double fn = 0;
    double fp = 0;
    for (int c = 0; c < m_NumClasses; c++) {
      tp += numTruePositives(c);
      fn += numFalseNegatives(c);
      fp += numFalsePositives(c);
    }
    return 2*tp / (2*tp + fn + fp);
  }

  /**
   * Sets the class prior probabilities.
   *
   * @param train the training instances used to determine the prior probabilities
   * @throws Exception if the class attribute of the instances is not set
   */

  public void setPriors(Instances train) throws Exception {

    m_NoPriors = false;

    if (!m_ClassIsNominal) {

      m_NumTrainClassVals = 0;
      m_TrainClassVals = null;
      m_TrainClassWeights = null;
      m_PriorEstimator = null;

      m_MinTarget = Double.MAX_VALUE;
      m_MaxTarget = -Double.MAX_VALUE;

      for (int i = 0; i < train.numInstances(); i++) {
        Instance currentInst = train.instance(i);
        if (!currentInst.classIsMissing()) {
          addNumericTrainClass(currentInst.classValue(), currentInst.weight());
        }
      }

      m_ClassPriors[0] = m_ClassPriorsSum = 0;
      for (int i = 0; i < train.numInstances(); i++) {
        if (!train.instance(i).classIsMissing()) {
          m_ClassPriors[0] += train.instance(i).classValue() * train.instance(i).weight();
          m_ClassPriorsSum += train.instance(i).weight();
        }
      }

    } else {
      for (int i = 0; i < m_NumClasses; i++) {
        m_ClassPriors[i] = 1;
      }
      m_ClassPriorsSum = m_NumClasses;
      for (int i = 0; i < train.numInstances(); i++) {
        if (!train.instance(i).classIsMissing()) {
 
        }
      }
      m_MaxTarget = m_NumClasses;
      m_MinTarget = 0;
    }
  }

  /**
   * Get the current weighted class counts.
   *
   * @return the weighted class counts
   */
  public double [] getClassPriors() {
    return m_ClassPriors;
  }

  /**
   * Updates the class prior probabilities or the mean respectively (when incrementally
   * training).
   *
   * @param instance the new training instance seen
   * @throws Exception if the class of the instance is not set
   */
  public void updatePriors(Instance instance) throws Exception {
    if (!instance.classIsMissing()) {
      if (!m_ClassIsNominal) {
        addNumericTrainClass(instance.classValue(), instance.weight());
        m_ClassPriors[0] += instance.classValue() * instance.weight();
        m_ClassPriorsSum += instance.weight();
      } else {
        m_ClassPriors[(int)instance.classValue()] += instance.weight();
        m_ClassPriorsSum += instance.weight();
      }
    }
  }

  /**
   * disables the use of priors, e.g., in case of de-serialized schemes
   * that have no access to the original training set, but are evaluated
   * on a set set.
   */
  public void useNoPriors() {
    m_NoPriors = true;
  }

  /**
   * Tests whether the current evaluation object is equal to another
   * evaluation object.
   *
   * @param obj the object to compare against
   * @return true if the two objects are equal
   */
  public boolean equals(Object obj) {

    if ((obj == null) || !(obj.getClass().equals(this.getClass()))) {
      return false;
    }
    Evaluation cmp = (Evaluation) obj;
    if (m_ClassIsNominal != cmp.m_ClassIsNominal) return false;
    if (m_NumClasses != cmp.m_NumClasses) return false;

    if (m_Incorrect != cmp.m_Incorrect) return false;
    if (m_Correct != cmp.m_Correct) return false;
    if (m_Unclassified != cmp.m_Unclassified) return false;
    if (m_MissingClass != cmp.m_MissingClass) return false;
    if (m_WithClass != cmp.m_WithClass) return false;

    if (m_SumErr != cmp.m_SumErr) return false;
    if (m_SumAbsErr != cmp.m_SumAbsErr) return false;
    if (m_SumSqrErr != cmp.m_SumSqrErr) return false;
    if (m_SumClass != cmp.m_SumClass) return false;
    if (m_SumSqrClass != cmp.m_SumSqrClass) return false;
    if (m_SumPredicted != cmp.m_SumPredicted) return false;
    if (m_SumSqrPredicted != cmp.m_SumSqrPredicted) return false;
    if (m_SumClassPredicted != cmp.m_SumClassPredicted) return false;

    if (m_ClassIsNominal) {
      for (int i = 0; i < m_NumClasses; i++) {
        for (int j = 0; j < m_NumClasses; j++) {
          if (m_ConfusionMatrix[i][j] != cmp.m_ConfusionMatrix[i][j]) {
            return false;
          }
        }
      }
    }

    return true;
  }

  /**
   * Make up the help string giving all the command line options.
   *
   * @param classifier the classifier to include options for
   * @param globalInfo include the global information string
   * for the classifier (if available).
   * @return a string detailing the valid command line options
   */
  protected static String makeOptionString(Classifier classifier,
                                           boolean globalInfo) {

    StringBuffer optionsText = new StringBuffer("");

    // General options
    optionsText.append("\n\nGeneral options:\n\n");
    optionsText.append("-h or -help\n");
    optionsText.append("\tOutput help information.\n");
    optionsText.append("-synopsis or -info\n");
    optionsText.append("\tOutput synopsis for classifier (use in conjunction "
        + " with -h)\n");
    optionsText.append("-t <name of training file>\n");
    optionsText.append("\tSets training file.\n");
    optionsText.append("-T <name of test file>\n");
    optionsText.append("\tSets test file. If missing, a cross-validation will be performed\n");
    optionsText.append("\ton the training data.\n");
    optionsText.append("-c <class index>\n");
    optionsText.append("\tSets index of class attribute (default: last).\n");
    optionsText.append("-x <number of folds>\n");
    optionsText.append("\tSets number of folds for cross-validation (default: 10).\n");
    optionsText.append("-no-cv\n");
    optionsText.append("\tDo not perform any cross validation.\n");
    optionsText.append("-split-percentage <percentage>\n");
    optionsText.append("\tSets the percentage for the train/test set split, e.g., 66.\n");
    optionsText.append("-preserve-order\n");
    optionsText.append("\tPreserves the order in the percentage split.\n");
    optionsText.append("-s <random number seed>\n");
    optionsText.append("\tSets random number seed for cross-validation or percentage split\n");
    optionsText.append("\t(default: 1).\n");
    optionsText.append("-m <name of file with cost matrix>\n");
    optionsText.append("\tSets file with cost matrix.\n");
    optionsText.append("-l <name of input file>\n");
    optionsText.append("\tSets model input file. In case the filename ends with '.xml',\n");
    optionsText.append("\ta PMML file is loaded or, if that fails, options are loaded\n");
    optionsText.append("\tfrom the XML file.\n");
    optionsText.append("-d <name of output file>\n");
    optionsText.append("\tSets model output file. In case the filename ends with '.xml',\n");
    optionsText.append("\tonly the options are saved to the XML file, not the model.\n");
    optionsText.append("-v\n");
    optionsText.append("\tOutputs no statistics for training data.\n");
    optionsText.append("-o\n");
    optionsText.append("\tOutputs statistics only, not the classifier.\n");
    optionsText.append("-i\n");
    optionsText.append("\tOutputs detailed information-retrieval");
    optionsText.append(" statistics for each class.\n");
    optionsText.append("-k\n");
    optionsText.append("\tOutputs information-theoretic statistics.\n");
    optionsText.append("-classifications \"weka.classifiers.evaluation.output.prediction.AbstractOutput + options\"\n");
    optionsText.append("\tUses the specified class for generating the classification output.\n");
    optionsText.append("\tE.g.: " + PlainText.class.getName() + "\n");
    optionsText.append("-p range\n");
    optionsText.append("\tOutputs predictions for test instances (or the train instances if\n");
    optionsText.append("\tno test instances provided and -no-cv is used), along with the \n");
    optionsText.append("\tattributes in the specified range (and nothing else). \n");
    optionsText.append("\tUse '-p 0' if no attributes are desired.\n");
    optionsText.append("\tDeprecated: use \"-classifications ...\" instead.\n");
    optionsText.append("-distribution\n");
    optionsText.append("\tOutputs the distribution instead of only the prediction\n");
    optionsText.append("\tin conjunction with the '-p' option (only nominal classes).\n");
    optionsText.append("\tDeprecated: use \"-classifications ...\" instead.\n");
    optionsText.append("-r\n");
    optionsText.append("\tOnly outputs cumulative margin distribution.\n");
    if (classifier instanceof Sourcable) {
      optionsText.append("-z <class name>\n");
      optionsText.append("\tOnly outputs the source representation"
          + " of the classifier,\n\tgiving it the supplied"
          + " name.\n");
    }
    if (classifier instanceof Drawable) {
      optionsText.append("-g\n");
      optionsText.append("\tOnly outputs the graph representation"
          + " of the classifier.\n");
    }
    optionsText.append("-xml filename | xml-string\n");
    optionsText.append("\tRetrieves the options from the XML-data instead of the "
        + "command line.\n");
    optionsText.append("-threshold-file <file>\n");
    optionsText.append("\tThe file to save the threshold data to.\n"
        + "\tThe format is determined by the extensions, e.g., '.arff' for ARFF \n"
        + "\tformat or '.csv' for CSV.\n");
    optionsText.append("-threshold-label <label>\n");
    optionsText.append("\tThe class label to determine the threshold data for\n"
        + "\t(default is the first label)\n");

    // Get scheme-specific options
    if (classifier instanceof OptionHandler) {
      optionsText.append("\nOptions specific to "
          + classifier.getClass().getName()
          + ":\n\n");
      Enumeration enu = ((OptionHandler)classifier).listOptions();
      while (enu.hasMoreElements()) {
        Option option = (Option) enu.nextElement();
        optionsText.append(option.synopsis() + '\n');
        optionsText.append(option.description() + "\n");
      }
    }

    // Get global information (if available)
    if (globalInfo) {
      try {
        String gi = getGlobalInfo(classifier);
        optionsText.append(gi);
      } catch (Exception ex) {
        // quietly ignore
      }
    }
    return optionsText.toString();
  }

  /**
   * Return the global info (if it exists) for the supplied classifier.
   *
   * @param classifier the classifier to get the global info for
   * @return the global info (synopsis) for the classifier
   * @throws Exception if there is a problem reflecting on the classifier
   */
  protected static String getGlobalInfo(Classifier classifier) throws Exception {
    BeanInfo bi = Introspector.getBeanInfo(classifier.getClass());
    MethodDescriptor[] methods;
    methods = bi.getMethodDescriptors();
    Object[] args = {};
    String result = "\nSynopsis for " + classifier.getClass().getName()
      + ":\n\n";

    for (int i = 0; i < methods.length; i++) {
      String name = methods[i].getDisplayName();
      Method meth = methods[i].getMethod();
      if (name.equals("globalInfo")) {
        String globalInfo = (String)(meth.invoke(classifier, args));
        result += globalInfo;
        break;
      }
    }

    return result;
  }

  /**
   * Method for generating indices for the confusion matrix.
   *
   * @param num 	integer to format
   * @param IDChars	the characters to use
   * @param IDWidth	the width of the entry
   * @return 		the formatted integer as a string
   */
  protected String num2ShortID(int num, char[] IDChars, int IDWidth) {

    char ID [] = new char [IDWidth];
    int i;

    for(i = IDWidth - 1; i >=0; i--) {
      ID[i] = IDChars[num % IDChars.length];
      num = num / IDChars.length - 1;
      if (num < 0) {
        break;
      }
    }
    for(i--; i >= 0; i--) {
      ID[i] = ' ';
    }

    return new String(ID);
  }

  /**
   * Convert a single prediction into a probability distribution
   * with all zero probabilities except the predicted value which
   * has probability 1.0.
   *
   * @param predictedClass the index of the predicted class
   * @return the probability distribution
   */
  protected double [] makeDistribution(double predictedClass) {

    double [] result = new double [m_NumClasses];
    if (Utils.isMissingValue(predictedClass)) {
      return result;
    }
    if (m_ClassIsNominal) {
      result[(int)predictedClass] = 1.0;
    } else {
      result[0] = predictedClass;
    }
    return result;
  }

  /**
   * Updates all the statistics about a classifiers performance for
   * the current test instance.
   *
   * @param predictedDistribution the probabilities assigned to
   * each class
   * @param instance the instance to be classified
   * @throws Exception if the class of the instance is not
   * set
   */
  protected void updateStatsForClassifier(double [] predictedDistribution,
      Instance instance)
  throws Exception {

    int actualClass = (int)instance.classValue();
    
    if (!instance.classIsMissing()) {
      updateMargins(predictedDistribution, actualClass, instance.weight());

      // Determine the predicted class (doesn't detect multiple
      // classifications)
      int predictedClass = -1;
      double bestProb = 0.0;
      for(int i = 0; i < m_NumClasses; i++) {
        if (predictedDistribution[i] > bestProb) {
          predictedClass = i;
          bestProb = predictedDistribution[i];
        }
      }

      m_WithClass += instance.weight();

      // Determine misclassification cost
      if (m_CostMatrix != null) {
        if (predictedClass < 0) {
          // For missing predictions, we assume the worst possible cost.
          // This is pretty harsh.
          // Perhaps we could take the negative of the cost of a correct
          // prediction (-m_CostMatrix.getElement(actualClass,actualClass)),
          // although often this will be zero
          m_TotalCost += instance.weight() * m_CostMatrix.getMaxCost(actualClass, instance);
        } else {
          m_TotalCost += instance.weight() * m_CostMatrix.getElement(actualClass, predictedClass,
              instance);
        }
      }

      // Update counts when no class was predicted
      if (predictedClass < 0) {
        m_Unclassified += instance.weight();
        return;
      }
      double predictedProb = Math.max(MIN_SF_PROB, predictedDistribution[actualClass]);
      double priorProb = Math.max(MIN_SF_PROB, m_ClassPriors[actualClass] / m_ClassPriorsSum);
      if (predictedProb >= priorProb) {
        m_SumKBInfo += (Utils.log2(predictedProb) - Utils.log2(priorProb)) * instance.weight();
      } else {
        m_SumKBInfo -= (Utils.log2(1.0-predictedProb) - Utils.log2(1.0-priorProb))
          * instance.weight();
      }

      m_SumSchemeEntropy -= Utils.log2(predictedProb) * instance.weight();
      m_SumPriorEntropy -= Utils.log2(priorProb) * instance.weight();

      updateNumericScores(predictedDistribution,
          makeDistribution(instance.classValue()),
          instance.weight());

      // Update coverage stats
      int[] indices = Utils.sort(predictedDistribution);
      double sum = 0, sizeOfRegions = 0;
      for (int i = predictedDistribution.length - 1; i >= 0; i--) {
        if (sum >= m_ConfLevel) {
          break;
        }
        sum += predictedDistribution[indices[i]];
        sizeOfRegions++;
        if (actualClass == indices[i]) {
          m_TotalCoverage += instance.weight();
        }
      }
      m_TotalSizeOfRegions += sizeOfRegions / (m_MaxTarget - m_MinTarget);

      // Update other stats
      m_ConfusionMatrix[actualClass][predictedClass] += instance.weight();
      if (predictedClass != actualClass) {
        m_Incorrect += instance.weight();
      } else {
        m_Correct += instance.weight();
      }
    } else {
      m_MissingClass += instance.weight();
    }
  }

  /**
   * Updates stats for interval estimator based on current test instance.
   *
   * @param classifier the interval estimator
   * @param classMissing the instance for which the intervals are computed, without a class value
   * @param classValue the class value of this instance
   * @throws Exception if intervals could not be computed successfully
   */
  protected void updateStatsForIntervalEstimator(IntervalEstimator  classifier, Instance classMissing,
                                                 double classValue) throws Exception {

    double[][] preds = classifier.predictIntervals(classMissing, m_ConfLevel);
    if (m_Predictions != null)
      ((NumericPrediction) m_Predictions.lastElement()).setPredictionIntervals(preds);
    for (int i = 0; i < preds.length; i++) {
      m_TotalSizeOfRegions += (preds[i][1] - preds[i][0]) / (m_MaxTarget - m_MinTarget);
    }
    for (int i = 0; i < preds.length; i++) {
      if ((preds[i][1] >= classValue) && (preds[i][0] <= classValue)) {
        m_TotalCoverage += classMissing.weight();
        break;
      }
    }
  }

  /**
   * Updates stats for conditional density estimator based on current test instance.
   *
   * @param classifier the conditional density estimator
   * @param classMissing the instance for which density is to be computed, without a class value
   * @param classValue the class value of this instance
   * @throws Exception if density could not be computed successfully
   */
  protected void updateStatsForConditionalDensityEstimator(ConditionalDensityEstimator classifier,
                                                           Instance classMissing,
                                                           double classValue) throws Exception {

    if (m_PriorEstimator == null) {
      setNumericPriorsFromBuffer();
    }
    m_SumSchemeEntropy -= classifier.logDensity(classMissing, classValue) * classMissing.weight() /
      Utils.log2;
    m_SumPriorEntropy -= m_PriorEstimator.logDensity(classValue) * classMissing.weight() /
      Utils.log2;
  }

  /**
   * Updates all the statistics about a predictors performance for
   * the current test instance.
   *
   * @param predictedValue the numeric value the classifier predicts
   * @param instance the instance to be classified
   * @throws Exception if the class of the instance is not set
   */
  protected void updateStatsForPredictor(double predictedValue, Instance instance)
    throws Exception {

    if (!instance.classIsMissing()){

      // Update stats
      m_WithClass += instance.weight();
      if (Utils.isMissingValue(predictedValue)) {
        m_Unclassified += instance.weight();
        return;
      }
      m_SumClass += instance.weight() * instance.classValue();
      m_SumSqrClass += instance.weight() * instance.classValue() * instance.classValue();
      m_SumClassPredicted += instance.weight() * instance.classValue() * predictedValue;
      m_SumPredicted += instance.weight() * predictedValue;
      m_SumSqrPredicted += instance.weight() * predictedValue * predictedValue;

      updateNumericScores(makeDistribution(predictedValue),
          makeDistribution(instance.classValue()),
          instance.weight());

    } else
      m_MissingClass += instance.weight();
  }

  /**
   * Update the cumulative record of classification margins.
   *
   * @param predictedDistribution the probability distribution predicted for
   * the current instance
   * @param actualClass the index of the actual instance class
   * @param weight the weight assigned to the instance
   */
  protected void updateMargins(double [] predictedDistribution,int actualClass, double weight){
	  //RANKING BEGIN
		  boolean ranking=false;
		  for(int i=0; i<predictedDistribution.length; i++){
			  if(predictedDistribution[i]>1){
				  ranking=true;
				  break;
			  }
		  }
		  if(!ranking){			  
	  			double probActual = predictedDistribution[actualClass];
	  			double probNext = 0;

	  			for(int i = 0; i < m_NumClasses; i++)
	  				if ((i != actualClass) &&
	  						(predictedDistribution[i] > probNext))
	  					probNext = predictedDistribution[i];
	 
	    
	  			double margin = probActual - probNext;
	  			int bin = (int)((margin + 1.0) / 2.0 * k_MarginResolution);
	  			m_MarginCounts[bin] += weight;
  		}
		//RANKING END
  }

  /**
   * Update the numeric accuracy measures. For numeric classes, the
   * accuracy is between the actual and predicted class values. For
   * nominal classes, the accuracy is between the actual and
   * predicted class probabilities.
   *
   * @param predicted the predicted values
   * @param actual the actual value
   * @param weight the weight associated with this prediction
   */
  protected void updateNumericScores(double [] predicted,
      double [] actual, double weight) {

    double diff;
    double sumErr = 0, sumAbsErr = 0, sumSqrErr = 0;
    double sumPriorAbsErr = 0, sumPriorSqrErr = 0;
    for(int i = 0; i < m_NumClasses; i++) {
      diff = predicted[i] - actual[i];
      sumErr += diff;
      sumAbsErr += Math.abs(diff);
      sumSqrErr += diff * diff;
      diff = (m_ClassPriors[i] / m_ClassPriorsSum) - actual[i];
      sumPriorAbsErr += Math.abs(diff);
      sumPriorSqrErr += diff * diff;
    }
    m_SumErr += weight * sumErr / m_NumClasses;
    m_SumAbsErr += weight * sumAbsErr / m_NumClasses;
    m_SumSqrErr += weight * sumSqrErr / m_NumClasses;
    m_SumPriorAbsErr += weight * sumPriorAbsErr / m_NumClasses;
    m_SumPriorSqrErr += weight * sumPriorSqrErr / m_NumClasses;
  }

  /**
   * Adds a numeric (non-missing) training class value and weight to
   * the buffer of stored values. Also updates minimum and maximum target value.
   *
   * @param classValue the class value
   * @param weight the instance weight
   */
  protected void addNumericTrainClass(double classValue, double weight) {

    // Update minimum and maximum target value
    if (classValue > m_MaxTarget) {
      m_MaxTarget = classValue;
    }
    if (classValue < m_MinTarget) {
      m_MinTarget = classValue;
    }

    // Update buffer
    if (m_TrainClassVals == null) {
      m_TrainClassVals = new double [100];
      m_TrainClassWeights = new double [100];
    }
    if (m_NumTrainClassVals == m_TrainClassVals.length) {
      double [] temp = new double [m_TrainClassVals.length * 2];
      System.arraycopy(m_TrainClassVals, 0,
          temp, 0, m_TrainClassVals.length);
      m_TrainClassVals = temp;

      temp = new double [m_TrainClassWeights.length * 2];
      System.arraycopy(m_TrainClassWeights, 0,
          temp, 0, m_TrainClassWeights.length);
      m_TrainClassWeights = temp;
    }
    m_TrainClassVals[m_NumTrainClassVals] = classValue;
    m_TrainClassWeights[m_NumTrainClassVals] = weight;
    m_NumTrainClassVals++;
  }

  /**
   * Sets up the priors for numeric class attributes from the
   * training class values that have been seen so far.
   */
  protected void setNumericPriorsFromBuffer() {

    m_PriorEstimator = new UnivariateKernelEstimator();
    for (int i = 0; i < m_NumTrainClassVals; i++) {
      m_PriorEstimator.addValue(m_TrainClassVals[i], m_TrainClassWeights[i]);
    }
  }

  /**
   * Returns the revision string.
   *
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6804 $");
  }
  
  //RANKING BEGIN
   
  /**
   * The following method computes the Spearman measures for label rankings.
   * Spearman footrule calculates the distance between two rankings by adding all differences between label positions.
   * Spearman correlation calculates a value between -1 and 1 giving a measure in how far two rankings correlate.
   * @param ranking1 (out of dataset)
   * @param ranking2 (prediction)
   * @return array position 0: Spearman footrule, array position 1: Spearman correlation.
   * @throws Exception 
   */
  public double[] spearman(double[] rnking1, double[] rnking2, double weight) throws Exception{	  
	  double output[] = new double[2];
	  double n=0;
	  double footrule=0;		//spearman Footrule
	  double d_iSquare=0;		//(d_i)^2 for spearman correlation
	  double infinity=0;
	  double[] ranking1 = new double[rnking1.length];
	  double[] ranking2 = new double[rnking2.length];
	  	  	  
	  //Counting missing labels.
	  for(int i=0; i<ranking1.length; i++){
		  if(rnking1[i]==Double.POSITIVE_INFINITY || rnking1[i]==Double.NaN)infinity++;
	  }
	  
	  
	  //By definition, spearman rank values start with 1, not with 0. So we have to add a one to all values in our rankings.
	  //We use new arrays in order to keep the input arrays as they are.

	  for(int a=0; a<ranking1.length; a++){
		  ranking1[a] = rnking1[a];
		  if(ranking1[a]!=Double.POSITIVE_INFINITY)
			  ranking1[a]++;
	  }
	  for(int a=0; a<ranking2.length; a++){
		  ranking2[a] = rnking2[a];
			  ranking2[a]++;
	  }
	  	  	  	  
	  //calculating footrule and spearman correlation.
	  for(int k=0; k<ranking1.length; k++){
		  //With this if-Block we leave out missing values, so they won't tear our values down.
		  if(ranking1[k]!=Double.POSITIVE_INFINITY && ranking1[k]!=0 && ranking2[k]!=0){
		      footrule+=Math.abs(ranking1[k]-ranking2[k]);
		      d_iSquare+=Math.pow(ranking1[k]-ranking2[k],2);  
		      n++;
		  }
	  }
	  	
	  
	  output[0]=footrule;
	  output[1]=( 1-(6*d_iSquare/(n*(Math.pow(n, 2)-1))));
	  if(infinity!=0){
		  //weighting of incomplete rankings.
		  output[0]*=(ranking1.length-infinity)/ranking1.length*weight;
		  output[1]*=(ranking1.length-infinity)/ranking1.length*weight;
	  }
	  
	  
	  
	  
	  //another version of spearman correlation.
	  if(rnking1.length!=rnking2.length){
			throw new Exception("The sizes of given rankings are not same.");
		}
  /*
		else{
			
			double dif_squ_sum = 0;
			for(int i=0; i<ranking1.length; i++)
				dif_squ_sum = dif_squ_sum + (ranking1[i] - ranking2[i]) * (ranking1[i] - ranking2[i]);
			//Careful! a.size()*a.size()*a.size() is huge.
			output[1]= 1-((6*dif_squ_sum/ranking1.length)/(ranking1.length+1))/(ranking1.length-1);//rho
		}*/
	  
	  
	  return output;
  }
   
  /**
   * Position error calculates the distance between the label with highest order
   * to its position inside of the prediction.
   * Example:
   *    given ranking:			prediction:
   *     L1>L2>L3				   L2>L3>L1			positionError = 2
   * @param label
   * @throws Exception 
   * @return the summed position error for one specified label.
   */
  /*public double positionError() throws Exception{
	  
	  int[] result = new int[rnkData.numInstances()];
	  double label = 0;
	  double predictedPos = 0;
	  int key = 0;
	  	  	  
	  PreferenceDenseInstance givenInst;

	  for(int i=0; i<rnkData.size(); i++){
		  givenInst = (PreferenceDenseInstance)rnkData.get(i);
		  
		  givenInst.g
		  
		  for(int j=0; ; j++){
			//Find position of label in rankings given from dataset.
			if(givenInst.getHashMap().get(j)!=null){
				double[] rank = givenInst.getHashMap().get(j));
				key = j;
		      		      
				label = rank[0];
		      
		      break;
			}
		  }
	      
		  //Find position of label in predicted rankings.
		  for(int l=0; l<givenInst.getPrediction().size(); l++){
			  int[] rnk = givenInst.getPrediction().get(key);
			  for(int m=0; m<rnk.length; m++){
				  if((int)label == rnk[m]){
					  predictedPos = m;
					  break;
				  }
			  }
		  }
		  
		  result[(int)label]+=predictedPos;  
	  }	  
	  
	  double output=0;
	  for(int i=0; i<result.length; i++){
		  output+=result[i];
	  }
	  return output/data.numInstances();
  }*/
  
  public Instances getRnkData(){
	  return rnkData;
  }
  
  /**
   * Calculates the relative error matrix for pairwise rankings and displays
   * them on the classifierPanel.
   * @return
 * @throws Exception 
   */
  public int[][] relativeErrorLabelPairs(int[][] pairRanks, int[] pred) throws Exception{

	  int[][] output;
	  if(rnkData.getLabels().size()!=0)
		  output = new int[rnkData.getLabels().size()][2*rnkData.getLabels().size()];
	  else
		  output = new int[RankUtilities.labels.size()][2*RankUtilities.labels.size()];

			  boolean empty = true;

		    	for(int k=0; k<pairRanks.length; k++){
		    		for(int m=0; m<pairRanks.length; m++){
		    			if(pairRanks[k][m]==1){
		    				empty=false;

		    				if(pred[k]<pred[m]){
		    					output[k][2*m]+=1;
				    		}

		    				else if(pred[m]<pred[k]){
		    					output[k][2*m+1]+=1;
				    		}
		    				else if(pred[m]==pred[k]&&!empty){
		    					output[k][2*m+1]+=1;
				    		}	
		    			}		
		    		}		    		
		    	}

	  return output;
  }
  
  /**
	 * This function is the same with kendall(ArrayList<Double> a, ArrayList<Double> b), 
	 * only takes different type of arguments. Note: the index of the ranking is started from 0. 
	 * @param a
	 * @param b
	 * @return Kendall tau rank correlation coefficient.
	 * @throws Exception 
	 * @author Weiwei Cheng
	 */
  public static double kendallsTau(double[] a, double[] b, double weight) throws Exception {

	  int infty=0;
	  for(int i=0; i<a.length; i++){
		  if(a[i]==Double.POSITIVE_INFINITY || a[i]==Double.NaN){
			  infty++;
		  }
	  }
	  
      double tau = -1;
      
      if(a.length!=b.length){
          throw new Exception("The sizes of given rankings are not same.");
      }
      
      else{
          /**
           * 1. Sort array a and get the new vision of b.
           */         
          double[] b_adjust = new double[a.length];
          for(int i=0; i<a.length; i++){
              int[] sort = Utils.sort(a);
              b_adjust[i] = b[sort[i]];
          }
          
          /**
           * 2. Compute the times of pairwise adjacent transpositions.
           */
          int transposition = 0;
          for(int i=0; i<a.length-1; i++)
              for(int j=i+1; j<a.length; j++){
                  if(b_adjust[i]>b_adjust[j])
                      transposition++;
          }
          
          /**
           * 3. Compute tau.
           */
          tau = 1-(4D*transposition/(a.length)/(a.length-1));
          if(infty==0)
        	  return tau*weight;
          else
        	  return tau*weight*((a.length-infty)/a.length);
  }
      
  }
  

 /**
  * Returns the kendall's Tau measure. 
  * @return
  */
  public double tau(){
	  return tau/train_testInst;
  }
  
  /**
   * Returning Spearman footrule[0] and correlation[1].
   * @return
   */
  public double[] spearman(){
	  double[] spear = new double[2];
	  spear[0] = spearman[0]/train_testInst;
	  spear[1] = spearman[1]/train_testInst;
	  return spear;
  }
  
  /**
   * Returning the matrix of pairwise label errors.
   * @return
   */
  public int[][] getRnkDataMatrix(){
	  return rnkDataMatrix;
  }
  
  //RANKING END
}
