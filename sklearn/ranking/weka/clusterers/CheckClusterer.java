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
 * CheckClusterer.java
 * Copyright (C) 2006 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.clusterers;

import weka.core.CheckScheme;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.MultiInstanceCapabilitiesHandler;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.SerializationHelper;
import weka.core.TestInstances;
import weka.core.Utils;
import weka.core.WeightedInstancesHandler;

import java.util.Enumeration;
import java.util.Random;
import java.util.Vector;

/**
 * Class for examining the capabilities and finding problems with 
 * clusterers. If you implement a clusterer using the WEKA.libraries,
 * you should run the checks on it to ensure robustness and correct
 * operation. Passing all the tests of this object does not mean
 * bugs in the clusterer don't exist, but this will help find some
 * common ones. <p/>
 * 
 * Typical usage: <p/>
 * <code>java weka.clusterers.CheckClusterer -W clusterer_name 
 * -- clusterer_options </code><p/>
 * 
 * CheckClusterer reports on the following:
 * <ul>
 *    <li> Clusterer abilities 
 *      <ul>
 *         <li> Possible command line options to the clusterer </li>
 *         <li> Whether the clusterer can predict nominal, numeric, string, 
 *              date or relational class attributes.</li>
 *         <li> Whether the clusterer can handle numeric predictor attributes </li>
 *         <li> Whether the clusterer can handle nominal predictor attributes </li>
 *         <li> Whether the clusterer can handle string predictor attributes </li>
 *         <li> Whether the clusterer can handle date predictor attributes </li>
 *         <li> Whether the clusterer can handle relational predictor attributes </li>
 *         <li> Whether the clusterer can handle multi-instance data </li>
 *         <li> Whether the clusterer can handle missing predictor values </li>
 *         <li> Whether the clusterer can handle instance weights </li>
 *      </ul>
 *    </li>
 *    <li> Correct functioning 
 *      <ul>
 *         <li> Correct initialisation during buildClusterer (i.e. no result
 *              changes when buildClusterer called repeatedly) </li>
 *         <li> Whether the clusterer alters the data pased to it 
 *              (number of instances, instance order, instance weights, etc) </li>
 *      </ul>
 *    </li>
 *    <li> Degenerate cases 
 *      <ul>
 *         <li> building clusterer with zero training instances </li>
 *         <li> all but one predictor attribute values missing </li>
 *         <li> all predictor attribute values missing </li>
 *         <li> all but one class values missing </li>
 *         <li> all class values missing </li>
 *      </ul>
 *    </li>
 * </ul>
 * Running CheckClusterer with the debug option set will output the 
 * training dataset for any failed tests.<p/>
 *
 * The <code>weka.clusterers.AbstractClustererTest</code> uses this
 * class to test all the clusterers. Any changes here, have to be 
 * checked in that abstract test class, too. <p/>
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -D
 *  Turn on debugging output.</pre>
 * 
 * <pre> -S
 *  Silent mode - prints nothing to stdout.</pre>
 * 
 * <pre> -N &lt;num&gt;
 *  The number of instances in the datasets (default 20).</pre>
 * 
 * <pre> -nominal &lt;num&gt;
 *  The number of nominal attributes (default 2).</pre>
 * 
 * <pre> -nominal-values &lt;num&gt;
 *  The number of values for nominal attributes (default 1).</pre>
 * 
 * <pre> -numeric &lt;num&gt;
 *  The number of numeric attributes (default 1).</pre>
 * 
 * <pre> -string &lt;num&gt;
 *  The number of string attributes (default 1).</pre>
 * 
 * <pre> -date &lt;num&gt;
 *  The number of date attributes (default 1).</pre>
 * 
 * <pre> -relational &lt;num&gt;
 *  The number of relational attributes (default 1).</pre>
 * 
 * <pre> -num-instances-relational &lt;num&gt;
 *  The number of instances in relational/bag attributes (default 10).</pre>
 * 
 * <pre> -words &lt;comma-separated-list&gt;
 *  The words to use in string attributes.</pre>
 * 
 * <pre> -word-separators &lt;chars&gt;
 *  The word separators to use in string attributes.</pre>
 * 
 * <pre> -W
 *  Full name of the clusterer analyzed.
 *  eg: weka.clusterers.SimpleKMeans
 *  (default weka.clusterers.SimpleKMeans)</pre>
 * 
 * <pre> 
 * Options specific to clusterer weka.clusterers.SimpleKMeans:
 * </pre>
 * 
 * <pre> -N &lt;num&gt;
 *  number of clusters.
 *  (default 2).</pre>
 * 
 * <pre> -V
 *  Display std. deviations for centroids.
 * </pre>
 * 
 * <pre> -M
 *  Replace missing values with mean/mode.
 * </pre>
 * 
 * <pre> -S &lt;num&gt;
 *  Random number seed.
 *  (default 10)</pre>
 * 
 <!-- options-end -->
 *
 * Options after -- are passed to the designated clusterer.<p/>
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 1.11 $
 * @see TestInstances
 */
public class CheckClusterer 
  extends CheckScheme {

  /*
   * Note about test methods:
   * - methods return array of booleans
   * - first index: success or not
   * - second index: acceptable or not (e.g., Exception is OK)
   *
   * FracPete (fracpete at waikato dot ac dot nz)
   */
  
  /*** The clusterer to be examined */
  protected Clusterer m_Clusterer = new SimpleKMeans();
  
  /**
   * default constructor
   */
  public CheckClusterer() {
    super();
    
    setNumInstances(40);
  }
  
  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector result = new Vector();
    
    Enumeration en = super.listOptions();
    while (en.hasMoreElements())
      result.addElement(en.nextElement());
    
    result.addElement(new Option(
        "\tFull name of the clusterer analyzed.\n"
        +"\teg: weka.clusterers.SimpleKMeans\n"
        + "\t(default weka.clusterers.SimpleKMeans)",
        "W", 1, "-W"));
    
    if ((m_Clusterer != null) 
        && (m_Clusterer instanceof OptionHandler)) {
      result.addElement(new Option("", "", 0, 
          "\nOptions specific to clusterer "
          + m_Clusterer.getClass().getName()
          + ":"));
      Enumeration enu = ((OptionHandler)m_Clusterer).listOptions();
      while (enu.hasMoreElements())
        result.addElement(enu.nextElement());
    }
    
    return result.elements();
  }
  
  /**
   * Parses a given list of options. <p/>
   *
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -D
   *  Turn on debugging output.</pre>
   * 
   * <pre> -S
   *  Silent mode - prints nothing to stdout.</pre>
   * 
   * <pre> -N &lt;num&gt;
   *  The number of instances in the datasets (default 20).</pre>
   * 
   * <pre> -nominal &lt;num&gt;
   *  The number of nominal attributes (default 2).</pre>
   * 
   * <pre> -nominal-values &lt;num&gt;
   *  The number of values for nominal attributes (default 1).</pre>
   * 
   * <pre> -numeric &lt;num&gt;
   *  The number of numeric attributes (default 1).</pre>
   * 
   * <pre> -string &lt;num&gt;
   *  The number of string attributes (default 1).</pre>
   * 
   * <pre> -date &lt;num&gt;
   *  The number of date attributes (default 1).</pre>
   * 
   * <pre> -relational &lt;num&gt;
   *  The number of relational attributes (default 1).</pre>
   * 
   * <pre> -num-instances-relational &lt;num&gt;
   *  The number of instances in relational/bag attributes (default 10).</pre>
   * 
   * <pre> -words &lt;comma-separated-list&gt;
   *  The words to use in string attributes.</pre>
   * 
   * <pre> -word-separators &lt;chars&gt;
   *  The word separators to use in string attributes.</pre>
   * 
   * <pre> -W
   *  Full name of the clusterer analyzed.
   *  eg: weka.clusterers.SimpleKMeans
   *  (default weka.clusterers.SimpleKMeans)</pre>
   * 
   * <pre> 
   * Options specific to clusterer weka.clusterers.SimpleKMeans:
   * </pre>
   * 
   * <pre> -N &lt;num&gt;
   *  number of clusters.
   *  (default 2).</pre>
   * 
   * <pre> -V
   *  Display std. deviations for centroids.
   * </pre>
   * 
   * <pre> -M
   *  Replace missing values with mean/mode.
   * </pre>
   * 
   * <pre> -S &lt;num&gt;
   *  Random number seed.
   *  (default 10)</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String      tmpStr;
    
    tmpStr = Utils.getOption('N', options);
    
    super.setOptions(options);
    
    if (tmpStr.length() != 0)
      setNumInstances(Integer.parseInt(tmpStr));
    else
      setNumInstances(40);

    tmpStr = Utils.getOption('W', options);
    if (tmpStr.length() == 0)
      tmpStr = weka.clusterers.SimpleKMeans.class.getName();
    setClusterer(
	(Clusterer) forName(
	    "weka.clusterers", 
	    Clusterer.class, 
	    tmpStr, 
	    Utils.partitionOptions(options)));
  }
  
  /**
   * Gets the current settings of the CheckClusterer.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    Vector        result;
    String[]      options;
    int           i;
    
    result = new Vector();
    
    options = super.getOptions();
    for (i = 0; i < options.length; i++)
      result.add(options[i]);
    
    if (getClusterer() != null) {
      result.add("-W");
      result.add(getClusterer().getClass().getName());
    }
    
    if ((m_Clusterer != null) && (m_Clusterer instanceof OptionHandler))
      options = ((OptionHandler) m_Clusterer).getOptions();
    else
      options = new String[0];
    
    if (options.length > 0) {
      result.add("--");
      for (i = 0; i < options.length; i++)
        result.add(options[i]);
    }
    
    return (String[]) result.toArray(new String[result.size()]);
  }
  
  /**
   * Begin the tests, reporting results to System.out
   */
  public void doTests() {
    
    if (getClusterer() == null) {
      println("\n=== No clusterer set ===");
      return;
    }
    println("\n=== Check on Clusterer: "
        + getClusterer().getClass().getName()
        + " ===\n");
    
    // Start tests
    println("--> Checking for interfaces");
    canTakeOptions();
    boolean updateable = updateableClusterer()[0];
    boolean weightedInstancesHandler = weightedInstancesHandler()[0];
    boolean multiInstanceHandler = multiInstanceHandler()[0];
    println("--> Clusterer tests");
    declaresSerialVersionUID();
    runTests(weightedInstancesHandler, multiInstanceHandler, updateable);
  }
  
  /**
   * Set the clusterer for testing. 
   *
   * @param newClusterer the Clusterer to use.
   */
  public void setClusterer(Clusterer newClusterer) {
    m_Clusterer = newClusterer;
  }
  
  /**
   * Get the clusterer used as the clusterer
   *
   * @return the clusterer used as the clusterer
   */
  public Clusterer getClusterer() {
    return m_Clusterer;
  }
  
  /**
   * Run a battery of tests
   *
   * @param weighted true if the clusterer says it handles weights
   * @param multiInstance true if the clusterer is a multi-instance clusterer
   * @param updateable true if the classifier is updateable
   */
  protected void runTests(boolean weighted, boolean multiInstance, boolean updateable) {
    
    boolean PNom = canPredict(true,  false, false, false, false, multiInstance)[0];
    boolean PNum = canPredict(false, true,  false, false, false, multiInstance)[0];
    boolean PStr = canPredict(false, false, true,  false, false, multiInstance)[0];
    boolean PDat = canPredict(false, false, false, true,  false, multiInstance)[0];
    boolean PRel;
    if (!multiInstance)
      PRel = canPredict(false, false, false, false,  true, multiInstance)[0];
    else
      PRel = false;

    if (PNom || PNum || PStr || PDat || PRel) {
      if (weighted)
        instanceWeights(PNom, PNum, PStr, PDat, PRel, multiInstance);
      
      canHandleZeroTraining(PNom, PNum, PStr, PDat, PRel, multiInstance);
      boolean handleMissingPredictors = canHandleMissing(PNom, PNum, PStr, PDat, PRel, 
          multiInstance, true, 20)[0];
      if (handleMissingPredictors)
        canHandleMissing(PNom, PNum, PStr, PDat, PRel, multiInstance, true, 100);
      
      correctBuildInitialisation(PNom, PNum, PStr, PDat, PRel, multiInstance);
      datasetIntegrity(PNom, PNum, PStr, PDat, PRel, multiInstance, handleMissingPredictors);
      if (updateable)
        updatingEquality(PNom, PNum, PStr, PDat, PRel, multiInstance);
    }
  }
  
  /**
   * Checks whether the scheme can take command line options.
   *
   * @return index 0 is true if the clusterer can take options
   */
  protected boolean[] canTakeOptions() {
    
    boolean[] result = new boolean[2];
    
    print("options...");
    if (m_Clusterer instanceof OptionHandler) {
      println("yes");
      if (m_Debug) {
        println("\n=== Full report ===");
        Enumeration enu = ((OptionHandler)m_Clusterer).listOptions();
        while (enu.hasMoreElements()) {
          Option option = (Option) enu.nextElement();
          print(option.synopsis() + "\n" 
              + option.description() + "\n");
        }
        println("\n");
      }
      result[0] = true;
    }
    else {
      println("no");
      result[0] = false;
    }
    
    return result;
  }
  
  /**
   * Checks whether the scheme can build models incrementally.
   *
   * @return index 0 is true if the clusterer can train incrementally
   */
  protected boolean[] updateableClusterer() {
    
    boolean[] result = new boolean[2];
    
    print("updateable clusterer...");
    if (m_Clusterer instanceof UpdateableClusterer) {
      println("yes");
      result[0] = true;
    }
    else {
      println("no");
      result[0] = false;
    }
    
    return result;
  }
  
  /**
   * Checks whether the scheme says it can handle instance weights.
   *
   * @return true if the clusterer handles instance weights
   */
  protected boolean[] weightedInstancesHandler() {
    
    boolean[] result = new boolean[2];
    
    print("weighted instances clusterer...");
    if (m_Clusterer instanceof WeightedInstancesHandler) {
      println("yes");
      result[0] = true;
    }
    else {
      println("no");
      result[0] = false;
    }
    
    return result;
  }
  
  /**
   * Checks whether the scheme handles multi-instance data.
   * 
   * @return true if the clusterer handles multi-instance data
   */
  protected boolean[] multiInstanceHandler() {
    boolean[] result = new boolean[2];
    
    print("multi-instance clusterer...");
    if (m_Clusterer instanceof MultiInstanceCapabilitiesHandler) {
      println("yes");
      result[0] = true;
    }
    else {
      println("no");
      result[0] = false;
    }
    
    return result;
  }
  
  /**
   * tests for a serialVersionUID. Fails in case the scheme doesn't declare
   * a UID.
   *
   * @return index 0 is true if the scheme declares a UID
   */
  protected boolean[] declaresSerialVersionUID() {
    boolean[] result = new boolean[2];
    
    print("serialVersionUID...");
    
    result[0] = !SerializationHelper.needsUID(m_Clusterer.getClass());
    
    if (result[0])
      println("yes");
    else
      println("no");
    
    return result;
  }
  
  /**
   * Checks basic prediction of the scheme, for simple non-troublesome
   * datasets.
   *
   * @param nominalPredictor if true use nominal predictor attributes
   * @param numericPredictor if true use numeric predictor attributes
   * @param stringPredictor if true use string predictor attributes
   * @param datePredictor if true use date predictor attributes
   * @param relationalPredictor if true use relational predictor attributes
   * @param multiInstance whether multi-instance is needed
   * @return index 0 is true if the test was passed, index 1 is true if test 
   *         was acceptable
   */
  protected boolean[] canPredict(
      boolean nominalPredictor,
      boolean numericPredictor, 
      boolean stringPredictor, 
      boolean datePredictor,
      boolean relationalPredictor,
      boolean multiInstance) {
    
    print("basic predict");
    printAttributeSummary(
        nominalPredictor, numericPredictor, stringPredictor, datePredictor, relationalPredictor, multiInstance);
    print("...");
    FastVector accepts = new FastVector();
    accepts.addElement("unary");
    accepts.addElement("binary");
    accepts.addElement("nominal");
    accepts.addElement("numeric");
    accepts.addElement("string");
    accepts.addElement("date");
    accepts.addElement("relational");
    accepts.addElement("multi-instance");
    accepts.addElement("not in classpath");
    int numTrain = getNumInstances(), missingLevel = 0;
    boolean predictorMissing = false;
    
    return runBasicTest(nominalPredictor, numericPredictor, stringPredictor, 
        datePredictor, relationalPredictor, 
        multiInstance,
        missingLevel, predictorMissing, 
        numTrain, 
        accepts);
  }
  
  /**
   * Checks whether the scheme can handle zero training instances.
   *
   * @param nominalPredictor if true use nominal predictor attributes
   * @param numericPredictor if true use numeric predictor attributes
   * @param stringPredictor if true use string predictor attributes
   * @param datePredictor if true use date predictor attributes
   * @param relationalPredictor if true use relational predictor attributes
   * @param multiInstance whether multi-instance is needed
   * @return index 0 is true if the test was passed, index 1 is true if test 
   *         was acceptable
   */
  protected boolean[] canHandleZeroTraining(
      boolean nominalPredictor,
      boolean numericPredictor, 
      boolean stringPredictor, 
      boolean datePredictor,
      boolean relationalPredictor,
      boolean multiInstance) {
    
    print("handle zero training instances");
    printAttributeSummary(
        nominalPredictor, numericPredictor, stringPredictor, datePredictor, relationalPredictor, multiInstance);
    print("...");
    FastVector accepts = new FastVector();
    accepts.addElement("train");
    accepts.addElement("value");
    int numTrain = 0, missingLevel = 0;
    boolean predictorMissing = false;
    
    return runBasicTest(
              nominalPredictor, numericPredictor, stringPredictor, 
              datePredictor, relationalPredictor, 
              multiInstance,
              missingLevel, predictorMissing,
              numTrain, 
              accepts);
  }
  
  /**
   * Checks whether the scheme correctly initialises models when 
   * buildClusterer is called. This test calls buildClusterer with
   * one training dataset. buildClusterer is then called on a training set 
   * with different structure, and then again with the original training set. 
   * If the equals method of the ClusterEvaluation class returns 
   * false, this is noted as incorrect build initialisation.
   * 
   * @param nominalPredictor if true use nominal predictor attributes
   * @param numericPredictor if true use numeric predictor attributes
   * @param stringPredictor if true use string predictor attributes
   * @param datePredictor if true use date predictor attributes
   * @param relationalPredictor if true use relational predictor attributes
   * @param multiInstance whether multi-instance is needed
   * @return index 0 is true if the test was passed
   */
  protected boolean[] correctBuildInitialisation(
      boolean nominalPredictor,
      boolean numericPredictor, 
      boolean stringPredictor, 
      boolean datePredictor,
      boolean relationalPredictor,
      boolean multiInstance) {

    boolean[] result = new boolean[2];
    
    print("correct initialisation during buildClusterer");
    printAttributeSummary(
        nominalPredictor, numericPredictor, stringPredictor, datePredictor, relationalPredictor, multiInstance);
    print("...");
    int numTrain = getNumInstances(), missingLevel = 0;
    boolean predictorMissing = false;
    
    Instances train1 = null;
    Instances train2 = null;
    Clusterer clusterer = null;
    ClusterEvaluation evaluation1A = null;
    ClusterEvaluation evaluation1B = null;
    ClusterEvaluation evaluation2 = null;
    boolean built = false;
    int stage = 0;
    try {
      
      // Make two train sets with different numbers of attributes
      train1 = makeTestDataset(42, numTrain, 
                               nominalPredictor    ? getNumNominal()    : 0,
                               numericPredictor    ? getNumNumeric()    : 0, 
                               stringPredictor     ? getNumString()     : 0, 
                               datePredictor       ? getNumDate()       : 0, 
                               relationalPredictor ? getNumRelational() : 0, 
                               multiInstance);
      train2 = makeTestDataset(84, numTrain, 
                               nominalPredictor    ? getNumNominal() + 1 : 0,
                               numericPredictor    ? getNumNumeric() + 1 : 0, 
                               stringPredictor     ? getNumString()      : 0, 
                               datePredictor       ? getNumDate()        : 0, 
                               relationalPredictor ? getNumRelational()  : 0, 
                               multiInstance);
      if (nominalPredictor && !multiInstance) {
        train1.deleteAttributeAt(0);
        train2.deleteAttributeAt(0);
      }
      if (missingLevel > 0) {
        addMissing(train1, missingLevel, predictorMissing);
        addMissing(train2, missingLevel, predictorMissing);
      }
      
      clusterer = AbstractClusterer.makeCopies(getClusterer(), 1)[0];
      evaluation1A = new ClusterEvaluation();
      evaluation1B = new ClusterEvaluation();
      evaluation2 = new ClusterEvaluation();
    } catch (Exception ex) {
      throw new Error("Error setting up for tests: " + ex.getMessage());
    }
    try {
      stage = 0;
      clusterer.buildClusterer(train1);
      built = true;
      evaluation1A.setClusterer(clusterer);
      evaluation1A.evaluateClusterer(train1);
      
      stage = 1;
      built = false;
      clusterer.buildClusterer(train2);
      built = true;
      evaluation2.setClusterer(clusterer);
      evaluation2.evaluateClusterer(train2);
      
      stage = 2;
      built = false;
      clusterer.buildClusterer(train1);
      built = true;
      evaluation1B.setClusterer(clusterer);
      evaluation1B.evaluateClusterer(train1);
      
      stage = 3;
      if (!evaluation1A.equals(evaluation1B)) {
        if (m_Debug) {
          println("\n=== Full report ===\n");
          println("First buildClusterer()");
          println(evaluation1A.clusterResultsToString() + "\n\n");
          println("Second buildClusterer()");
          println(evaluation1B.clusterResultsToString() + "\n\n");
        }
        throw new Exception("Results differ between buildClusterer calls");
      }
      println("yes");
      result[0] = true;
      
      if (false && m_Debug) {
        println("\n=== Full report ===\n");
        println("First buildClusterer()");
        println(evaluation1A.clusterResultsToString() + "\n\n");
        println("Second buildClusterer()");
        println(evaluation1B.clusterResultsToString() + "\n\n");
      }
    } 
    catch (Exception ex) {
      println("no");
      result[0] = false;
      if (m_Debug) {
        println("\n=== Full Report ===");
        print("Problem during");
        if (built) {
          print(" testing");
        } else {
          print(" training");
        }
        switch (stage) {
          case 0:
            print(" of dataset 1");
            break;
          case 1:
            print(" of dataset 2");
            break;
          case 2:
            print(" of dataset 1 (2nd build)");
            break;
          case 3:
            print(", comparing results from builds of dataset 1");
            break;	  
        }
        println(": " + ex.getMessage() + "\n");
        println("here are the datasets:\n");
        println("=== Train1 Dataset ===\n"
            + train1.toString() + "\n");
        println("=== Train2 Dataset ===\n"
            + train2.toString() + "\n");
      }
    }
    
    return result;
  }
  
  /**
   * Checks basic missing value handling of the scheme. If the missing
   * values cause an exception to be thrown by the scheme, this will be
   * recorded.
   *
   * @param nominalPredictor if true use nominal predictor attributes
   * @param numericPredictor if true use numeric predictor attributes
   * @param stringPredictor if true use string predictor attributes
   * @param datePredictor if true use date predictor attributes
   * @param relationalPredictor if true use relational predictor attributes
   * @param multiInstance whether multi-instance is needed
   * @param predictorMissing true if the missing values may be in 
   * the predictors
   * @param missingLevel the percentage of missing values
   * @return index 0 is true if the test was passed, index 1 is true if test 
   *         was acceptable
   */
  protected boolean[] canHandleMissing(
      boolean nominalPredictor,
      boolean numericPredictor, 
      boolean stringPredictor, 
      boolean datePredictor,
      boolean relationalPredictor,
      boolean multiInstance,
      boolean predictorMissing,
      int missingLevel) {
    
    if (missingLevel == 100)
      print("100% ");
    print("missing");
    if (predictorMissing) {
      print(" predictor");
    }
    print(" values");
    printAttributeSummary(
        nominalPredictor, numericPredictor, stringPredictor, datePredictor, relationalPredictor, multiInstance);
    print("...");
    FastVector accepts = new FastVector();
    accepts.addElement("missing");
    accepts.addElement("value");
    accepts.addElement("train");
    int numTrain = getNumInstances();
    
    return runBasicTest(nominalPredictor, numericPredictor, stringPredictor, 
        datePredictor, relationalPredictor, 
        multiInstance,
        missingLevel, predictorMissing,
        numTrain, 
        accepts);
  }
  
  /**
   * Checks whether the clusterer can handle instance weights.
   * This test compares the clusterer performance on two datasets
   * that are identical except for the training weights. If the 
   * results change, then the clusterer must be using the weights. It
   * may be possible to get a false positive from this test if the 
   * weight changes aren't significant enough to induce a change
   * in clusterer performance (but the weights are chosen to minimize
   * the likelihood of this).
   *
   * @param nominalPredictor if true use nominal predictor attributes
   * @param numericPredictor if true use numeric predictor attributes
   * @param stringPredictor if true use string predictor attributes
   * @param datePredictor if true use date predictor attributes
   * @param relationalPredictor if true use relational predictor attributes
   * @param multiInstance whether multi-instance is needed
   * @return index 0 true if the test was passed
   */
  protected boolean[] instanceWeights(
      boolean nominalPredictor,
      boolean numericPredictor, 
      boolean stringPredictor, 
      boolean datePredictor,
      boolean relationalPredictor,
      boolean multiInstance) {
    
    print("clusterer uses instance weights");
    printAttributeSummary(
        nominalPredictor, numericPredictor, stringPredictor, datePredictor, relationalPredictor, multiInstance);
    print("...");
    int numTrain = 2*getNumInstances(), missingLevel = 0;
    boolean predictorMissing = false;
    
    boolean[] result = new boolean[2];
    Instances train = null;
    Clusterer [] clusterers = null;
    ClusterEvaluation evaluationB = null;
    ClusterEvaluation evaluationI = null;
    boolean built = false;
    boolean evalFail = false;
    try {
      train = makeTestDataset(42, numTrain, 
                              nominalPredictor    ? getNumNominal() + 1 : 0,
                              numericPredictor    ? getNumNumeric() + 1 : 0, 
                              stringPredictor     ? getNumString()      : 0, 
                              datePredictor       ? getNumDate()        : 0, 
                              relationalPredictor ? getNumRelational()  : 0, 
                              multiInstance);
      if (nominalPredictor && !multiInstance)
        train.deleteAttributeAt(0);
      if (missingLevel > 0)
        addMissing(train, missingLevel, predictorMissing);
      clusterers = AbstractClusterer.makeCopies(getClusterer(), 2);
      evaluationB = new ClusterEvaluation();
      evaluationI = new ClusterEvaluation();
      clusterers[0].buildClusterer(train);
      evaluationB.setClusterer(clusterers[0]);
    } catch (Exception ex) {
      throw new Error("Error setting up for tests: " + ex.getMessage());
    }
    try {
      
      // Now modify instance weights and re-built/test
      for (int i = 0; i < train.numInstances(); i++) {
        train.instance(i).setWeight(0);
      }
      Random random = new Random(1);
      for (int i = 0; i < train.numInstances() / 2; i++) {
        int inst = Math.abs(random.nextInt()) % train.numInstances();
        int weight = Math.abs(random.nextInt()) % 10 + 1;
        train.instance(inst).setWeight(weight);
      }
      clusterers[1].buildClusterer(train);
      built = true;
      evaluationI.setClusterer(clusterers[1]);
      if (evaluationB.equals(evaluationI)) {
        //	println("no");
        evalFail = true;
        throw new Exception("evalFail");
      }
      
      println("yes");
      result[0] = true;
    } catch (Exception ex) {
      println("no");
      result[0] = false;
      
      if (m_Debug) {
        println("\n=== Full Report ===");
        
        if (evalFail) {
          println("Results don't differ between non-weighted and "
              + "weighted instance models.");
          println("Here are the results:\n");
          println("\nboth methods\n");
          println(evaluationB.clusterResultsToString());
        } else {
          print("Problem during");
          if (built) {
            print(" testing");
          } else {
            print(" training");
          }
          println(": " + ex.getMessage() + "\n");
        }
        println("Here is the dataset:\n");
        println("=== Train Dataset ===\n"
            + train.toString() + "\n");
        println("=== Train Weights ===\n");
        for (int i = 0; i < train.numInstances(); i++) {
          println(" " + (i + 1) 
              + "    " + train.instance(i).weight());
        }
      }
    }
    
    return result;
  }
  
  /**
   * Checks whether the scheme alters the training dataset during
   * training. If the scheme needs to modify the training
   * data it should take a copy of the training data. Currently checks
   * for changes to header structure, number of instances, order of
   * instances, instance weights.
   *
   * @param nominalPredictor if true use nominal predictor attributes
   * @param numericPredictor if true use numeric predictor attributes
   * @param stringPredictor if true use string predictor attributes
   * @param datePredictor if true use date predictor attributes
   * @param relationalPredictor if true use relational predictor attributes
   * @param multiInstance whether multi-instance is needed
   * @param predictorMissing true if we know the clusterer can handle
   * (at least) moderate missing predictor values
   * @return index 0 is true if the test was passed
   */
  protected boolean[] datasetIntegrity(
      boolean nominalPredictor,
      boolean numericPredictor, 
      boolean stringPredictor, 
      boolean datePredictor,
      boolean relationalPredictor,
      boolean multiInstance,
      boolean predictorMissing) {
    
    print("clusterer doesn't alter original datasets");
    printAttributeSummary(
        nominalPredictor, numericPredictor, stringPredictor, datePredictor, relationalPredictor, multiInstance);
    print("...");
    int numTrain = getNumInstances(), missingLevel = 20;
    
    boolean[] result = new boolean[2];
    Instances train = null;
    Clusterer clusterer = null;
    try {
      train = makeTestDataset(42, numTrain, 
                              nominalPredictor    ? getNumNominal()    : 0,
                              numericPredictor    ? getNumNumeric()    : 0, 
                              stringPredictor     ? getNumString()     : 0, 
                              datePredictor       ? getNumDate()       : 0, 
                              relationalPredictor ? getNumRelational() : 0, 
                              multiInstance);
      if (nominalPredictor && !multiInstance)
        train.deleteAttributeAt(0);
      if (missingLevel > 0)
        addMissing(train, missingLevel, predictorMissing);
      clusterer = AbstractClusterer.makeCopies(getClusterer(), 1)[0];
    } catch (Exception ex) {
      throw new Error("Error setting up for tests: " + ex.getMessage());
    }
    try {
      Instances trainCopy = new Instances(train);
      clusterer.buildClusterer(trainCopy);
      compareDatasets(train, trainCopy);
      
      println("yes");
      result[0] = true;
    } catch (Exception ex) {
      println("no");
      result[0] = false;
      
      if (m_Debug) {
        println("\n=== Full Report ===");
        print("Problem during training");
        println(": " + ex.getMessage() + "\n");
        println("Here is the dataset:\n");
        println("=== Train Dataset ===\n"
            + train.toString() + "\n");
      }
    }
    
    return result;
  }
  
  /**
   * Checks whether an updateable scheme produces the same model when
   * trained incrementally as when batch trained. The model itself
   * cannot be compared, so we compare the evaluation on test data
   * for both models. It is possible to get a false positive on this
   * test (likelihood depends on the classifier).
   *
   * @param nominalPredictor if true use nominal predictor attributes
   * @param numericPredictor if true use numeric predictor attributes
   * @param stringPredictor if true use string predictor attributes
   * @param datePredictor if true use date predictor attributes
   * @param relationalPredictor if true use relational predictor attributes
   * @param multiInstance whether multi-instance is needed
   * @return index 0 is true if the test was passed
   */
  protected boolean[] updatingEquality(
      boolean nominalPredictor,
      boolean numericPredictor, 
      boolean stringPredictor, 
      boolean datePredictor,
      boolean relationalPredictor,
      boolean multiInstance) {
    
    print("incremental training produces the same results"
        + " as batch training");
    printAttributeSummary(
        nominalPredictor, numericPredictor, stringPredictor, datePredictor, relationalPredictor, multiInstance);
    print("...");
    int numTrain = getNumInstances(), missingLevel = 0;
    boolean predictorMissing = false, classMissing = false;
    
    boolean[] result = new boolean[2];
    Instances train = null;
    Clusterer[] clusterers = null;
    ClusterEvaluation evaluationB = null;
    ClusterEvaluation evaluationI = null;
    boolean built = false;
    try {
      train = makeTestDataset(42, numTrain, 
                              nominalPredictor    ? getNumNominal()    : 0,
                              numericPredictor    ? getNumNumeric()    : 0, 
                              stringPredictor     ? getNumString()     : 0, 
                              datePredictor       ? getNumDate()       : 0, 
                              relationalPredictor ? getNumRelational() : 0, 
                              multiInstance);
      if (missingLevel > 0)
        addMissing(train, missingLevel, predictorMissing, classMissing);
      clusterers = AbstractClusterer.makeCopies(getClusterer(), 2);
      evaluationB = new ClusterEvaluation();
      evaluationI = new ClusterEvaluation();
      clusterers[0].buildClusterer(train);
      evaluationB.setClusterer(clusterers[0]);
    } catch (Exception ex) {
      throw new Error("Error setting up for tests: " + ex.getMessage());
    }
    try {
      clusterers[1].buildClusterer(new Instances(train, 0));
      for (int i = 0; i < train.numInstances(); i++) {
        ((UpdateableClusterer)clusterers[1]).updateClusterer(
            train.instance(i));
      }
      built = true;
      evaluationI.setClusterer(clusterers[1]);
      if (!evaluationB.equals(evaluationI)) {
        println("no");
        result[0] = false;
        
        if (m_Debug) {
          println("\n=== Full Report ===");
          println("Results differ between batch and "
              + "incrementally built models.\n"
              + "Depending on the classifier, this may be OK");
          println("Here are the results:\n");
          println("\nbatch built results\n" + evaluationB.clusterResultsToString());
          println("\nincrementally built results\n" + evaluationI.clusterResultsToString());
          println("Here are the datasets:\n");
          println("=== Train Dataset ===\n"
              + train.toString() + "\n");
        }
      }
      else {
        println("yes");
        result[0] = true;
      }
    } catch (Exception ex) {
      result[0] = false;
      
      print("Problem during");
      if (built)
        print(" testing");
      else
        print(" training");
      println(": " + ex.getMessage() + "\n");
    }
    
    return result;
  }
  
  /**
   * Runs a text on the datasets with the given characteristics.
   * 
   * @param nominalPredictor if true use nominal predictor attributes
   * @param numericPredictor if true use numeric predictor attributes
   * @param stringPredictor if true use string predictor attributes
   * @param datePredictor if true use date predictor attributes
   * @param relationalPredictor if true use relational predictor attributes
   * @param multiInstance whether multi-instance is needed
   * @param missingLevel the percentage of missing values
   * @param predictorMissing true if the missing values may be in 
   * the predictors
   * @param numTrain the number of instances in the training set
   * @param accepts the acceptable string in an exception
   * @return index 0 is true if the test was passed, index 1 is true if test 
   *         was acceptable
   */
  protected boolean[] runBasicTest(boolean nominalPredictor,
      boolean numericPredictor, 
      boolean stringPredictor,
      boolean datePredictor,
      boolean relationalPredictor,
      boolean multiInstance,
      int missingLevel,
      boolean predictorMissing,
      int numTrain,
      FastVector accepts) {
    
    boolean[] result = new boolean[2];
    Instances train = null;
    Clusterer clusterer = null;
    try {
      train = makeTestDataset(42, numTrain, 
                              nominalPredictor    ? getNumNominal()    : 0,
                              numericPredictor    ? getNumNumeric()    : 0, 
                              stringPredictor     ? getNumString()     : 0,
                              datePredictor       ? getNumDate()       : 0,
                              relationalPredictor ? getNumRelational() : 0,
                              multiInstance);
      if (nominalPredictor && !multiInstance)
        train.deleteAttributeAt(0);
      if (missingLevel > 0)
        addMissing(train, missingLevel, predictorMissing);
      clusterer = AbstractClusterer.makeCopies(getClusterer(), 1)[0];
    } catch (Exception ex) {
      ex.printStackTrace();
      throw new Error("Error setting up for tests: " + ex.getMessage());
    }
    try {
      clusterer.buildClusterer(train);
      println("yes");
      result[0] = true;
    } 
    catch (Exception ex) {
      boolean acceptable = false;
      String msg = ex.getMessage().toLowerCase();
      for (int i = 0; i < accepts.size(); i++) {
        if (msg.indexOf((String)accepts.elementAt(i)) >= 0) {
          acceptable = true;
        }
      }
      
      println("no" + (acceptable ? " (OK error message)" : ""));
      result[1] = acceptable;
      
      if (m_Debug) {
        println("\n=== Full Report ===");
        print("Problem during training");
        println(": " + ex.getMessage() + "\n");
        if (!acceptable) {
          if (accepts.size() > 0) {
            print("Error message doesn't mention ");
            for (int i = 0; i < accepts.size(); i++) {
              if (i != 0) {
                print(" or ");
              }
              print('"' + (String)accepts.elementAt(i) + '"');
            }
          }
          println("here is the dataset:\n");
          println("=== Train Dataset ===\n"
              + train.toString() + "\n");
        }
      }
    }
    
    return result;
  }
  
  /**
   * Add missing values to a dataset.
   *
   * @param data the instances to add missing values to
   * @param level the level of missing values to add (if positive, this
   * is the probability that a value will be set to missing, if negative
   * all but one value will be set to missing (not yet implemented))
   * @param predictorMissing if true, predictor attributes will be modified
   */
  protected void addMissing(Instances data, int level, boolean predictorMissing) {
    
    Random random = new Random(1);
    for (int i = 0; i < data.numInstances(); i++) {
      Instance current = data.instance(i);
      for (int j = 0; j < data.numAttributes(); j++) {
        if (predictorMissing) {
          if (Math.abs(random.nextInt()) % 100 < level)
            current.setMissing(j);
        }
      }
    }
  }
  
  /**
   * Make a simple set of instances with variable position of the class 
   * attribute, which can later be modified for use in specific tests.
   *
   * @param seed the random number seed
   * @param numInstances the number of instances to generate
   * @param numNominal the number of nominal attributes
   * @param numNumeric the number of numeric attributes
   * @param numString the number of string attributes
   * @param numDate the number of date attributes
   * @param numRelational the number of relational attributes
   * @param multiInstance whether the dataset should a multi-instance dataset
   * @return the test dataset
   * @throws Exception if the dataset couldn't be generated
   * @see TestInstances#CLASS_IS_LAST
   */
  protected Instances makeTestDataset(int seed, int numInstances, 
                                      int numNominal, int numNumeric, 
                                      int numString, int numDate,
                                      int numRelational,
                                      boolean multiInstance)
  throws Exception {
    
    TestInstances dataset = new TestInstances();
    
    dataset.setSeed(seed);
    dataset.setNumInstances(numInstances);
    dataset.setNumNominal(numNominal);
    dataset.setNumNumeric(numNumeric);
    dataset.setNumString(numString);
    dataset.setNumDate(numDate);
    dataset.setNumRelational(numRelational);
    dataset.setClassIndex(TestInstances.NO_CLASS);
    dataset.setMultiInstance(multiInstance);
    
    return dataset.generate();
  }
  
  /**
   * Print out a short summary string for the dataset characteristics
   *
   * @param nominalPredictor true if nominal predictor attributes are present
   * @param numericPredictor true if numeric predictor attributes are present
   * @param stringPredictor true if string predictor attributes are present
   * @param datePredictor true if date predictor attributes are present
   * @param relationalPredictor true if relational predictor attributes are present
   * @param multiInstance whether multi-instance is needed
   */
  protected void printAttributeSummary(boolean nominalPredictor, 
                                       boolean numericPredictor, 
                                       boolean stringPredictor, 
                                       boolean datePredictor, 
                                       boolean relationalPredictor, 
                                       boolean multiInstance) {
    
    String str = "";

    if (numericPredictor)
      str += "numeric";
    
    if (nominalPredictor) {
      if (str.length() > 0)
        str += " & ";
      str += "nominal";
    }
    
    if (stringPredictor) {
      if (str.length() > 0)
        str += " & ";
      str += "string";
    }
    
    if (datePredictor) {
      if (str.length() > 0)
        str += " & ";
      str += "date";
    }
    
    if (relationalPredictor) {
      if (str.length() > 0)
        str += " & ";
      str += "relational";
    }
    
    str = " (" + str + " predictors)";
    
    print(str);
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 1.11 $");
  }
  
  /**
   * Test method for this class
   * 
   * @param args the commandline options
   */
  public static void main(String [] args) {
    runCheck(new CheckClusterer(), args);
  }
}

