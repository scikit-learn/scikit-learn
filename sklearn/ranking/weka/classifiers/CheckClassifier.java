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
 *    CheckClassifier.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers;

import weka.core.Attribute;
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
import weka.core.labelranking.PreferenceAttribute;

import java.util.Enumeration;
import java.util.Random;
import java.util.Vector;

/**
 * Class for examining the capabilities and finding problems with
 * classifiers. If you implement a classifier using the WEKA.libraries,
 * you should run the checks on it to ensure robustness and correct
 * operation. Passing all the tests of this object does not mean
 * bugs in the classifier don't exist, but this will help find some
 * common ones. <p/>
 *
 * Typical usage: <p/>
 * <code>java weka.classifiers.CheckClassifier -W classifier_name
 * classifier_options </code><p/>
 *
 * CheckClassifier reports on the following:
 * <ul>
 *    <li> Classifier abilities
 *      <ul>
 *         <li> Possible command line options to the classifier </li>
 *         <li> Whether the classifier can predict nominal, numeric, string,
 *              date or relational class attributes. Warnings will be displayed if
 *              performance is worse than ZeroR </li>
 *         <li> Whether the classifier can be trained incrementally </li>
 *         <li> Whether the classifier can handle numeric predictor attributes </li>
 *         <li> Whether the classifier can handle nominal predictor attributes </li>
 *         <li> Whether the classifier can handle string predictor attributes </li>
 *         <li> Whether the classifier can handle date predictor attributes </li>
 *         <li> Whether the classifier can handle relational predictor attributes </li>
 *         <li> Whether the classifier can handle multi-instance data </li>
 *         <li> Whether the classifier can handle missing predictor values </li>
 *         <li> Whether the classifier can handle missing class values </li>
 *         <li> Whether a nominal classifier only handles 2 class problems </li>
 *         <li> Whether the classifier can handle instance weights </li>
 *      </ul>
 *    </li>
 *    <li> Correct functioning
 *      <ul>
 *         <li> Correct initialisation during buildClassifier (i.e. no result
 *              changes when buildClassifier called repeatedly) </li>
 *         <li> Whether incremental training produces the same results
 *              as during non-incremental training (which may or may not
 *              be OK) </li>
 *         <li> Whether the classifier alters the data pased to it
 *              (number of instances, instance order, instance weights, etc) </li>
 *         <li> Whether the toString() method works correctly before the
 *              classifier has been built. </li>
 *      </ul>
 *    </li>
 *    <li> Degenerate cases
 *      <ul>
 *         <li> building classifier with zero training instances </li>
 *         <li> all but one predictor attribute values missing </li>
 *         <li> all predictor attribute values missing </li>
 *         <li> all but one class values missing </li>
 *         <li> all class values missing </li>
 *      </ul>
 *    </li>
 * </ul>
 * Running CheckClassifier with the debug option set will output the
 * training and test datasets for any failed tests.<p/>
 *
 * The <code>weka.classifiers.AbstractClassifierTest</code> uses this
 * class to test all the classifiers. Any changes here, have to be
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
 *  Full name of the classifier analysed.
 *  eg: weka.classifiers.bayes.NaiveBayes
 *  (default weka.classifiers.rules.ZeroR)</pre>
 *
 * <pre>
 * Options specific to classifier weka.classifiers.rules.ZeroR:
 * </pre>
 *
 * <pre> -D
 *  If set, classifier is run in debug mode and
 *  may output additional info to the console</pre>
 *
 <!-- options-end -->
 *
 * Options after -- are passed to the designated classifier.<p/>
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 6041 $
 * @see TestInstances
 */
public class CheckClassifier
  extends CheckScheme {

  /*
   * Note about test methods:
   * - methods return array of booleans
   * - first index: success or not
   * - second index: acceptable or not (e.g., Exception is OK)
   * - in case the performance is worse than that of ZeroR both indices are true
   *
   * FracPete (fracpete at waikato dot ac dot nz)
   */

  /*** The classifier to be examined */
  protected Classifier m_Classifier = new weka.classifiers.rules.ZeroR();

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
        "\tFull name of the classifier analysed.\n"
        +"\teg: weka.classifiers.bayes.NaiveBayes\n"
        + "\t(default weka.classifiers.rules.ZeroR)",
        "W", 1, "-W"));

    if ((m_Classifier != null)
        && (m_Classifier instanceof OptionHandler)) {
      result.addElement(new Option("", "", 0,
          "\nOptions specific to classifier "
          + m_Classifier.getClass().getName()
          + ":"));
      Enumeration enu = ((OptionHandler)m_Classifier).listOptions();
      while (enu.hasMoreElements())
        result.addElement(enu.nextElement());
    }

    return result.elements();
  }

  /**
   * Parses a given list of options.
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
   *  Full name of the classifier analysed.
   *  eg: weka.classifiers.bayes.NaiveBayes
   *  (default weka.classifiers.rules.ZeroR)</pre>
   *
   * <pre>
   * Options specific to classifier weka.classifiers.rules.ZeroR:
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
    String      tmpStr;

    super.setOptions(options);

    tmpStr = Utils.getOption('W', options);
    if (tmpStr.length() == 0)
      tmpStr = weka.classifiers.rules.ZeroR.class.getName();
    setClassifier(
        (Classifier) forName(
          "weka.classifiers",
          Classifier.class,
          tmpStr,
          Utils.partitionOptions(options)));
  }

  /**
   * Gets the current settings of the CheckClassifier.
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

    if (getClassifier() != null) {
      result.add("-W");
      result.add(getClassifier().getClass().getName());
    }

    if ((m_Classifier != null) && (m_Classifier instanceof OptionHandler))
      options = ((OptionHandler) m_Classifier).getOptions();
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

    if (getClassifier() == null) {
      println("\n=== No classifier set ===");
      return;
    }
    println("\n=== Check on Classifier: "
        + getClassifier().getClass().getName()
        + " ===\n");

    // Start tests
    m_ClasspathProblems = false;
    println("--> Checking for interfaces");
    canTakeOptions();
    boolean updateableClassifier = updateableClassifier()[0];
    boolean weightedInstancesHandler = weightedInstancesHandler()[0];
    boolean multiInstanceHandler = multiInstanceHandler()[0];
    println("--> Classifier tests");
    declaresSerialVersionUID();
    testToString();
    testsPerClassType(PreferenceAttribute.RANKING,    updateableClassifier, weightedInstancesHandler, multiInstanceHandler);
    testsPerClassType(Attribute.NOMINAL,    updateableClassifier, weightedInstancesHandler, multiInstanceHandler);
    testsPerClassType(Attribute.NUMERIC,    updateableClassifier, weightedInstancesHandler, multiInstanceHandler);
    testsPerClassType(Attribute.DATE,       updateableClassifier, weightedInstancesHandler, multiInstanceHandler);
    testsPerClassType(Attribute.STRING,     updateableClassifier, weightedInstancesHandler, multiInstanceHandler);
    testsPerClassType(Attribute.RELATIONAL, updateableClassifier, weightedInstancesHandler, multiInstanceHandler);
  }

  /**
   * Set the classifier for boosting.
   *
   * @param newClassifier the Classifier to use.
   */
  public void setClassifier(Classifier newClassifier) {
    m_Classifier = newClassifier;
  }

  /**
   * Get the classifier used as the classifier
   *
   * @return the classifier used as the classifier
   */
  public Classifier getClassifier() {
    return m_Classifier;
  }

  /**
   * Run a battery of tests for a given class attribute type
   *
   * @param classType true if the class attribute should be numeric
   * @param updateable true if the classifier is updateable
   * @param weighted true if the classifier says it handles weights
   * @param multiInstance true if the classifier is a multi-instance classifier
   */
  protected void testsPerClassType(int classType,
                                   boolean updateable,
                                   boolean weighted,
                                   boolean multiInstance) {

    boolean PNom = canPredict(true,  false, false, false, false, multiInstance, classType)[0];
    boolean PNum = canPredict(false, true,  false, false, false, multiInstance, classType)[0];
    boolean PStr = canPredict(false, false, true,  false, false, multiInstance, classType)[0];
    boolean PDat = canPredict(false, false, false, true,  false, multiInstance, classType)[0];
    boolean PRel;
    if (!multiInstance)
      PRel = canPredict(false, false, false, false,  true, multiInstance, classType)[0];
    else
      PRel = false;

    if (PNom || PNum || PStr || PDat || PRel) {
      if (weighted)
        instanceWeights(PNom, PNum, PStr, PDat, PRel, multiInstance, classType);

      canHandleOnlyClass(PNom, PNum, PStr, PDat, PRel, classType);

      if (classType == Attribute.NOMINAL || classType == PreferenceAttribute.RANKING)
        canHandleNClasses(PNom, PNum, PStr, PDat, PRel, multiInstance, 4);

      if (!multiInstance) {
        canHandleClassAsNthAttribute(PNom, PNum, PStr, PDat, PRel, multiInstance, classType, 0);
        canHandleClassAsNthAttribute(PNom, PNum, PStr, PDat, PRel, multiInstance, classType, 1);
      }

      canHandleZeroTraining(PNom, PNum, PStr, PDat, PRel, multiInstance, classType);
      boolean handleMissingPredictors = canHandleMissing(PNom, PNum, PStr, PDat, PRel,
          multiInstance, classType,
          true, false, 20)[0];
      if (handleMissingPredictors)
        canHandleMissing(PNom, PNum, PStr, PDat, PRel, multiInstance, classType, true, false, 100);

      boolean handleMissingClass = canHandleMissing(PNom, PNum, PStr, PDat, PRel,
          multiInstance, classType,
          false, true, 20)[0];
      if (handleMissingClass)
        canHandleMissing(PNom, PNum, PStr, PDat, PRel, multiInstance, classType, false, true, 100);

      correctBuildInitialisation(PNom, PNum, PStr, PDat, PRel, multiInstance, classType);
      datasetIntegrity(PNom, PNum, PStr, PDat, PRel, multiInstance, classType,
          handleMissingPredictors, handleMissingClass);
      doesntUseTestClassVal(PNom, PNum, PStr, PDat, PRel, multiInstance, classType);
      if (updateable)
        updatingEquality(PNom, PNum, PStr, PDat, PRel, multiInstance, classType);
    }
  }

  /**
   * Checks whether the scheme's toString() method works even though the
   * classifies hasn't been built yet.
   *
   * @return index 0 is true if the toString() method works fine
   */
  protected boolean[] testToString() {
    boolean[] result = new boolean[2];

    print("toString...");

    try {
      Classifier copy = (Classifier) m_Classifier.getClass().newInstance();
      copy.toString();
      result[0] = true;
      println("yes");
    }
    catch (Exception e) {
      result[0] = false;
      println("no");
      if (m_Debug) {
        println("\n=== Full report ===");
        e.printStackTrace();
        println("\n");
      }
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

    result[0] = !SerializationHelper.needsUID(m_Classifier.getClass());

    if (result[0])
      println("yes");
    else
      println("no");

    return result;
  }

  /**
   * Checks whether the scheme can take command line options.
   *
   * @return index 0 is true if the classifier can take options
   */
  protected boolean[] canTakeOptions() {

    boolean[] result = new boolean[2];

    print("options...");
    if (m_Classifier instanceof OptionHandler) {
      println("yes");
      if (m_Debug) {
        println("\n=== Full report ===");
        Enumeration enu = ((OptionHandler)m_Classifier).listOptions();
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
   * @return index 0 is true if the classifier can train incrementally
   */
  protected boolean[] updateableClassifier() {

    boolean[] result = new boolean[2];

    print("updateable classifier...");
    if (m_Classifier instanceof UpdateableClassifier) {
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
   * @return true if the classifier handles instance weights
   */
  protected boolean[] weightedInstancesHandler() {

    boolean[] result = new boolean[2];

    print("weighted instances classifier...");
    if (m_Classifier instanceof WeightedInstancesHandler) {
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
   * @return true if the classifier handles multi-instance data
   */
  protected boolean[] multiInstanceHandler() {
    boolean[] result = new boolean[2];

    print("multi-instance classifier...");
    if (m_Classifier instanceof MultiInstanceCapabilitiesHandler) {
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
   * Checks basic prediction of the scheme, for simple non-troublesome
   * datasets.
   *
   * @param nominalPredictor if true use nominal predictor attributes
   * @param numericPredictor if true use numeric predictor attributes
   * @param stringPredictor if true use string predictor attributes
   * @param datePredictor if true use date predictor attributes
   * @param relationalPredictor if true use relational predictor attributes
   * @param multiInstance whether multi-instance is needed
   * @param classType the class type (NOMINAL, NUMERIC, etc.)
   * @return index 0 is true if the test was passed, index 1 is true if test
   *         was acceptable
   */
  protected boolean[] canPredict(
      boolean nominalPredictor,
      boolean numericPredictor,
      boolean stringPredictor,
      boolean datePredictor,
      boolean relationalPredictor,
      boolean multiInstance,
      int classType) {

    print("basic predict");
    printAttributeSummary(
        nominalPredictor, numericPredictor, stringPredictor, datePredictor, relationalPredictor, multiInstance, classType);
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
    int numTrain = getNumInstances(), numTest = getNumInstances(),
    numClasses = 2, missingLevel = 0;
    boolean predictorMissing = false, classMissing = false;

    return runBasicTest(nominalPredictor, numericPredictor, stringPredictor,
        datePredictor, relationalPredictor,
        multiInstance,
        classType,
        missingLevel, predictorMissing, classMissing,
        numTrain, numTest, numClasses,
        accepts);
  }

  /**
   * Checks whether the scheme can handle data that contains only the class
   * attribute. If a scheme cannot build a proper model with that data, it
   * should default back to a ZeroR model.
   *
   * @param nominalPredictor if true use nominal predictor attributes
   * @param numericPredictor if true use numeric predictor attributes
   * @param stringPredictor if true use string predictor attributes
   * @param datePredictor if true use date predictor attributes
   * @param relationalPredictor if true use relational predictor attributes
   * @param classType the class type (NOMINAL, NUMERIC, etc.)
   * @return index 0 is true if the test was passed, index 1 is true if test
   *         was acceptable
   */
  protected boolean[] canHandleOnlyClass(
      boolean nominalPredictor,
      boolean numericPredictor,
      boolean stringPredictor,
      boolean datePredictor,
      boolean relationalPredictor,
      int classType) {

    print("only class in data");
    printAttributeSummary(
        nominalPredictor, numericPredictor, stringPredictor, datePredictor, relationalPredictor, false, classType);
    print("...");
    FastVector accepts = new FastVector();
    accepts.addElement("class");
    accepts.addElement("zeror");
    int numTrain = getNumInstances(), numTest = getNumInstances(),
    missingLevel = 0;
    boolean predictorMissing = false, classMissing = false;

    return runBasicTest(false, false, false, false, false,
                        false,
                        classType,
                        missingLevel, predictorMissing, classMissing,
                        numTrain, numTest, 2,
                        accepts);
  }

  /**
   * Checks whether nominal schemes can handle more than two classes.
   * If a scheme is only designed for two-class problems it should
   * throw an appropriate exception for multi-class problems.
   *
   * @param nominalPredictor if true use nominal predictor attributes
   * @param numericPredictor if true use numeric predictor attributes
   * @param stringPredictor if true use string predictor attributes
   * @param datePredictor if true use date predictor attributes
   * @param relationalPredictor if true use relational predictor attributes
   * @param multiInstance whether multi-instance is needed
   * @param numClasses the number of classes to test
   * @return index 0 is true if the test was passed, index 1 is true if test
   *         was acceptable
   */
  protected boolean[] canHandleNClasses(
      boolean nominalPredictor,
      boolean numericPredictor,
      boolean stringPredictor,
      boolean datePredictor,
      boolean relationalPredictor,
      boolean multiInstance,
      int numClasses) {

    print("more than two class problems");
    printAttributeSummary(
        nominalPredictor, numericPredictor, stringPredictor, datePredictor, relationalPredictor, multiInstance, Attribute.NOMINAL);
    print("...");
    FastVector accepts = new FastVector();
    accepts.addElement("number");
    accepts.addElement("class");
    int numTrain = getNumInstances(), numTest = getNumInstances(),
    missingLevel = 0;
    boolean predictorMissing = false, classMissing = false;

    return runBasicTest(nominalPredictor, numericPredictor, stringPredictor,
                        datePredictor, relationalPredictor,
                        multiInstance,
                        Attribute.NOMINAL,
                        missingLevel, predictorMissing, classMissing,
                        numTrain, numTest, numClasses,
                        accepts);
  }

  /**
   * Checks whether the scheme can handle class attributes as Nth attribute.
   *
   * @param nominalPredictor if true use nominal predictor attributes
   * @param numericPredictor if true use numeric predictor attributes
   * @param stringPredictor if true use string predictor attributes
   * @param datePredictor if true use date predictor attributes
   * @param relationalPredictor if true use relational predictor attributes
   * @param multiInstance whether multi-instance is needed
   * @param classType the class type (NUMERIC, NOMINAL, etc.)
   * @param classIndex the index of the class attribute (0-based, -1 means last attribute)
   * @return index 0 is true if the test was passed, index 1 is true if test
   *         was acceptable
   * @see TestInstances#CLASS_IS_LAST
   */
  protected boolean[] canHandleClassAsNthAttribute(
      boolean nominalPredictor,
      boolean numericPredictor,
      boolean stringPredictor,
      boolean datePredictor,
      boolean relationalPredictor,
      boolean multiInstance,
      int classType,
      int classIndex) {

    if (classIndex == TestInstances.CLASS_IS_LAST)
      print("class attribute as last attribute");
    else
      print("class attribute as " + (classIndex + 1) + ". attribute");
    printAttributeSummary(
        nominalPredictor, numericPredictor, stringPredictor, datePredictor, relationalPredictor, multiInstance, classType);
    print("...");
    FastVector accepts = new FastVector();
    int numTrain = getNumInstances(), numTest = getNumInstances(), numClasses = 2,
    missingLevel = 0;
    boolean predictorMissing = false, classMissing = false;

    return runBasicTest(nominalPredictor, numericPredictor, stringPredictor,
                        datePredictor, relationalPredictor,
                        multiInstance,
                        classType,
                        classIndex,
                        missingLevel, predictorMissing, classMissing,
                        numTrain, numTest, numClasses,
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
   * @param classType the class type (NUMERIC, NOMINAL, etc.)
   * @return index 0 is true if the test was passed, index 1 is true if test
   *         was acceptable
   */
  protected boolean[] canHandleZeroTraining(
      boolean nominalPredictor,
      boolean numericPredictor,
      boolean stringPredictor,
      boolean datePredictor,
      boolean relationalPredictor,
      boolean multiInstance,
      int classType) {

    print("handle zero training instances");
    printAttributeSummary(
        nominalPredictor, numericPredictor, stringPredictor, datePredictor, relationalPredictor, multiInstance, classType);
    print("...");
    FastVector accepts = new FastVector();
    accepts.addElement("train");
    accepts.addElement("value");
    int numTrain = 0, numTest = getNumInstances(), numClasses = 2,
    missingLevel = 0;
    boolean predictorMissing = false, classMissing = false;

    return runBasicTest(
              nominalPredictor, numericPredictor, stringPredictor,
              datePredictor, relationalPredictor,
              multiInstance,
              classType,
              missingLevel, predictorMissing, classMissing,
              numTrain, numTest, numClasses,
              accepts);
  }

  /**
   * Checks whether the scheme correctly initialises models when
   * buildClassifier is called. This test calls buildClassifier with
   * one training dataset and records performance on a test set.
   * buildClassifier is then called on a training set with different
   * structure, and then again with the original training set. The
   * performance on the test set is compared with the original results
   * and any performance difference noted as incorrect build initialisation.
   *
   * @param nominalPredictor if true use nominal predictor attributes
   * @param numericPredictor if true use numeric predictor attributes
   * @param stringPredictor if true use string predictor attributes
   * @param datePredictor if true use date predictor attributes
   * @param relationalPredictor if true use relational predictor attributes
   * @param multiInstance whether multi-instance is needed
   * @param classType the class type (NUMERIC, NOMINAL, etc.)
   * @return index 0 is true if the test was passed, index 1 is true if the
   *         scheme performs worse than ZeroR, but without error (index 0 is
   *         false)
   */
  protected boolean[] correctBuildInitialisation(
      boolean nominalPredictor,
      boolean numericPredictor,
      boolean stringPredictor,
      boolean datePredictor,
      boolean relationalPredictor,
      boolean multiInstance,
      int classType) {

    boolean[] result = new boolean[2];

    print("correct initialisation during buildClassifier");
    printAttributeSummary(
        nominalPredictor, numericPredictor, stringPredictor, datePredictor, relationalPredictor, multiInstance, classType);
    print("...");
    int numTrain = getNumInstances(), numTest = getNumInstances(),
    numClasses = 2, missingLevel = 0;
    boolean predictorMissing = false, classMissing = false;

    Instances train1 = null;
    Instances test1 = null;
    Instances train2 = null;
    Instances test2 = null;
    Classifier classifier = null;
    Evaluation evaluation1A = null;
    Evaluation evaluation1B = null;
    Evaluation evaluation2 = null;
    boolean built = false;
    int stage = 0;
    try {

      // Make two sets of train/test splits with different
      // numbers of attributes
      train1 = makeTestDataset(42, numTrain,
                               nominalPredictor    ? getNumNominal()    : 0,
                               numericPredictor    ? getNumNumeric()    : 0,
                               stringPredictor     ? getNumString()     : 0,
                               datePredictor       ? getNumDate()       : 0,
                               relationalPredictor ? getNumRelational() : 0,
                               numClasses,
                               classType,
                               multiInstance);
      train2 = makeTestDataset(84, numTrain,
                               nominalPredictor    ? getNumNominal() + 1 : 0,
                               numericPredictor    ? getNumNumeric() + 1 : 0,
                               stringPredictor     ? getNumString()      : 0,
                               datePredictor       ? getNumDate()        : 0,
                               relationalPredictor ? getNumRelational()  : 0,
                               numClasses,
                               classType,
                               multiInstance);
      test1 = makeTestDataset(24, numTest,
                              nominalPredictor     ? getNumNominal()    : 0,
                              numericPredictor     ? getNumNumeric()    : 0,
                              stringPredictor      ? getNumString()     : 0,
                              datePredictor        ? getNumDate()       : 0,
                              relationalPredictor  ? getNumRelational() : 0,
                              numClasses,
                              classType,
                              multiInstance);
      test2 = makeTestDataset(48, numTest,
                              nominalPredictor    ? getNumNominal() + 1 : 0,
                              numericPredictor    ? getNumNumeric() + 1 : 0,
                              stringPredictor     ? getNumString()      : 0,
                              datePredictor       ? getNumDate()        : 0,
                              relationalPredictor ? getNumRelational()  : 0,
                              numClasses,
                              classType,
                              multiInstance);
      if (missingLevel > 0) {
        addMissing(train1, missingLevel, predictorMissing, classMissing);
        addMissing(test1, Math.min(missingLevel,50), predictorMissing,
            classMissing);
        addMissing(train2, missingLevel, predictorMissing, classMissing);
        addMissing(test2, Math.min(missingLevel,50), predictorMissing,
            classMissing);
      }

      classifier = AbstractClassifier.makeCopies(getClassifier(), 1)[0];
      evaluation1A = new Evaluation(train1);
      evaluation1B = new Evaluation(train1);
      evaluation2 = new Evaluation(train2);
    } catch (Exception ex) {
      throw new Error("Error setting up for tests: " + ex.getMessage());
    }
    try {
      stage = 0;
      classifier.buildClassifier(train1);
      built = true;
      if (!testWRTZeroR(classifier, evaluation1A, train1, test1)[0]) {
        throw new Exception("Scheme performs worse than ZeroR");
      }

      stage = 1;
      built = false;
      classifier.buildClassifier(train2);
      built = true;
      if (!testWRTZeroR(classifier, evaluation2, train2, test2)[0]) {
        throw new Exception("Scheme performs worse than ZeroR");
      }

      stage = 2;
      built = false;
      classifier.buildClassifier(train1);
      built = true;
      if (!testWRTZeroR(classifier, evaluation1B, train1, test1)[0]) {
        throw new Exception("Scheme performs worse than ZeroR");
      }

      stage = 3;
      if (!evaluation1A.equals(evaluation1B)) {
        if (m_Debug) {
          println("\n=== Full report ===\n"
              + evaluation1A.toSummaryString("\nFirst buildClassifier()",
                  true)
                  + "\n\n");
          println(
              evaluation1B.toSummaryString("\nSecond buildClassifier()",
                  true)
                  + "\n\n");
        }
        throw new Exception("Results differ between buildClassifier calls");
      }
      println("yes");
      result[0] = true;

      if (false && m_Debug) {
        println("\n=== Full report ===\n"
            + evaluation1A.toSummaryString("\nFirst buildClassifier()",
                true)
                + "\n\n");
        println(
            evaluation1B.toSummaryString("\nSecond buildClassifier()",
                true)
                + "\n\n");
      }
    }
    catch (Exception ex) {
      String msg = ex.getMessage().toLowerCase();
      if (msg.indexOf("worse than zeror") >= 0) {
        println("warning: performs worse than ZeroR");
        result[0] = (stage < 1);
        result[1] = (stage < 1);
      } else {
        println("no");
        result[0] = false;
      }
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
        println("=== Test1 Dataset ===\n"
            + test1.toString() + "\n\n");
        println("=== Train2 Dataset ===\n"
            + train2.toString() + "\n");
        println("=== Test2 Dataset ===\n"
            + test2.toString() + "\n\n");
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
   * @param classType the class type (NUMERIC, NOMINAL, etc.)
   * @param predictorMissing true if the missing values may be in
   * the predictors
   * @param classMissing true if the missing values may be in the class
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
      int classType,
      boolean predictorMissing,
      boolean classMissing,
      int missingLevel) {

    if (missingLevel == 100)
      print("100% ");
    print("missing");
    if (predictorMissing) {
      print(" predictor");
      if (classMissing)
        print(" and");
    }
    if (classMissing)
      print(" class");
    print(" values");
    printAttributeSummary(
        nominalPredictor, numericPredictor, stringPredictor, datePredictor, relationalPredictor, multiInstance, classType);
    print("...");
    FastVector accepts = new FastVector();
    accepts.addElement("missing");
    accepts.addElement("value");
    accepts.addElement("train");
    int numTrain = getNumInstances(), numTest = getNumInstances(),
    numClasses = 2;

    return runBasicTest(nominalPredictor, numericPredictor, stringPredictor,
        datePredictor, relationalPredictor,
        multiInstance,
        classType,
        missingLevel, predictorMissing, classMissing,
        numTrain, numTest, numClasses,
        accepts);
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
   * @param classType the class type (NUMERIC, NOMINAL, etc.)
   * @return index 0 is true if the test was passed
   */
  protected boolean[] updatingEquality(
      boolean nominalPredictor,
      boolean numericPredictor,
      boolean stringPredictor,
      boolean datePredictor,
      boolean relationalPredictor,
      boolean multiInstance,
      int classType) {

    print("incremental training produces the same results"
        + " as batch training");
    printAttributeSummary(
        nominalPredictor, numericPredictor, stringPredictor, datePredictor, relationalPredictor, multiInstance, classType);
    print("...");
    int numTrain = getNumInstances(), numTest = getNumInstances(),
    numClasses = 2, missingLevel = 0;
    boolean predictorMissing = false, classMissing = false;

    boolean[] result = new boolean[2];
    Instances train = null;
    Instances test = null;
    Classifier [] classifiers = null;
    Evaluation evaluationB = null;
    Evaluation evaluationI = null;
    boolean built = false;
    try {
      train = makeTestDataset(42, numTrain,
                              nominalPredictor    ? getNumNominal()    : 0,
                              numericPredictor    ? getNumNumeric()    : 0,
                              stringPredictor     ? getNumString()     : 0,
                              datePredictor       ? getNumDate()       : 0,
                              relationalPredictor ? getNumRelational() : 0,
                              numClasses,
                              classType,
                              multiInstance);
      test = makeTestDataset(24, numTest,
                             nominalPredictor    ? getNumNominal()    : 0,
                             numericPredictor    ? getNumNumeric()    : 0,
                             stringPredictor     ? getNumString()     : 0,
                             datePredictor       ? getNumDate()       : 0,
                             relationalPredictor ? getNumRelational() : 0,
                             numClasses,
                             classType,
                             multiInstance);
      if (missingLevel > 0) {
        addMissing(train, missingLevel, predictorMissing, classMissing);
        addMissing(test, Math.min(missingLevel, 50), predictorMissing,
            classMissing);
      }
      classifiers = AbstractClassifier.makeCopies(getClassifier(), 2);
      evaluationB = new Evaluation(train);
      evaluationI = new Evaluation(train);
      classifiers[0].buildClassifier(train);
      testWRTZeroR(classifiers[0], evaluationB, train, test);
    } catch (Exception ex) {
      throw new Error("Error setting up for tests: " + ex.getMessage());
    }
    try {
      classifiers[1].buildClassifier(new Instances(train, 0));
      for (int i = 0; i < train.numInstances(); i++) {
        ((UpdateableClassifier)classifiers[1]).updateClassifier(
            train.instance(i));
      }
      built = true;
      testWRTZeroR(classifiers[1], evaluationI, train, test);
      if (!evaluationB.equals(evaluationI)) {
        println("no");
        result[0] = false;

        if (m_Debug) {
          println("\n=== Full Report ===");
          println("Results differ between batch and "
              + "incrementally built models.\n"
              + "Depending on the classifier, this may be OK");
          println("Here are the results:\n");
          println(evaluationB.toSummaryString(
              "\nbatch built results\n", true));
          println(evaluationI.toSummaryString(
              "\nincrementally built results\n", true));
          println("Here are the datasets:\n");
          println("=== Train Dataset ===\n"
              + train.toString() + "\n");
          println("=== Test Dataset ===\n"
              + test.toString() + "\n\n");
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
   * Checks whether the classifier erroneously uses the class
   * value of test instances (if provided). Runs the classifier with
   * test instance class values set to missing and compares with results
   * when test instance class values are left intact.
   *
   * @param nominalPredictor if true use nominal predictor attributes
   * @param numericPredictor if true use numeric predictor attributes
   * @param stringPredictor if true use string predictor attributes
   * @param datePredictor if true use date predictor attributes
   * @param relationalPredictor if true use relational predictor attributes
   * @param multiInstance whether multi-instance is needed
   * @param classType the class type (NUMERIC, NOMINAL, etc.)
   * @return index 0 is true if the test was passed
   */
  protected boolean[] doesntUseTestClassVal(
      boolean nominalPredictor,
      boolean numericPredictor,
      boolean stringPredictor,
      boolean datePredictor,
      boolean relationalPredictor,
      boolean multiInstance,
      int classType) {

    print("classifier ignores test instance class vals");
    printAttributeSummary(
        nominalPredictor, numericPredictor, stringPredictor, datePredictor, relationalPredictor, multiInstance, classType);
    print("...");
    int numTrain = 2*getNumInstances(), numTest = getNumInstances(),
    numClasses = 2, missingLevel = 0;
    boolean predictorMissing = false, classMissing = false;

    boolean[] result = new boolean[2];
    Instances train = null;
    Instances test = null;
    Classifier [] classifiers = null;
    boolean evalFail = false;
    try {
      train = makeTestDataset(42, numTrain,
                              nominalPredictor    ? getNumNominal() + 1 : 0,
                              numericPredictor    ? getNumNumeric() + 1 : 0,
                              stringPredictor     ? getNumString()      : 0,
                              datePredictor       ? getNumDate()        : 0,
                              relationalPredictor ? getNumRelational()  : 0,
                              numClasses,
                              classType,
                              multiInstance);
      test = makeTestDataset(24, numTest,
                             nominalPredictor    ? getNumNominal() + 1 : 0,
                             numericPredictor    ? getNumNumeric() + 1 : 0,
                             stringPredictor     ? getNumString()      : 0,
                             datePredictor       ? getNumDate()        : 0,
                             relationalPredictor ? getNumRelational()  : 0,
                             numClasses,
                             classType,
                             multiInstance);
      if (missingLevel > 0) {
        addMissing(train, missingLevel, predictorMissing, classMissing);
        addMissing(test, Math.min(missingLevel, 50), predictorMissing,
            classMissing);
      }
      classifiers = AbstractClassifier.makeCopies(getClassifier(), 2);
      classifiers[0].buildClassifier(train);
      classifiers[1].buildClassifier(train);
    } catch (Exception ex) {
      throw new Error("Error setting up for tests: " + ex.getMessage());
    }
    try {

      // Now set test values to missing when predicting
      for (int i = 0; i < test.numInstances(); i++) {
        Instance testInst = test.instance(i);
        Instance classMissingInst = (Instance)testInst.copy();
        classMissingInst.setDataset(test);
        classMissingInst.setClassMissing();
        double [] dist0 = classifiers[0].distributionForInstance(testInst);
        double [] dist1 = classifiers[1].distributionForInstance(classMissingInst);
        for (int j = 0; j < dist0.length; j++) {
          // ignore, if both are NaNs
          if (Double.isNaN(dist0[j]) && Double.isNaN(dist1[j])) {
            if (getDebug())
              System.out.println("Both predictions are NaN!");
            continue;
          }
          // distribution different?
          if (dist0[j] != dist1[j]) {
            throw new Exception("Prediction different for instance " + (i + 1));
          }
        }
      }

      println("yes");
      result[0] = true;
    } catch (Exception ex) {
      println("no");
      result[0] = false;

      if (m_Debug) {
        println("\n=== Full Report ===");

        if (evalFail) {
          println("Results differ between non-missing and "
              + "missing test class values.");
        } else {
          print("Problem during testing");
          println(": " + ex.getMessage() + "\n");
        }
        println("Here are the datasets:\n");
        println("=== Train Dataset ===\n"
            + train.toString() + "\n");
        println("=== Train Weights ===\n");
        for (int i = 0; i < train.numInstances(); i++) {
          println(" " + (i + 1)
              + "    " + train.instance(i).weight());
        }
        println("=== Test Dataset ===\n"
            + test.toString() + "\n\n");
        println("(test weights all 1.0\n");
      }
    }

    return result;
  }

  /**
   * Checks whether the classifier can handle instance weights.
   * This test compares the classifier performance on two datasets
   * that are identical except for the training weights. If the
   * results change, then the classifier must be using the weights. It
   * may be possible to get a false positive from this test if the
   * weight changes aren't significant enough to induce a change
   * in classifier performance (but the weights are chosen to minimize
   * the likelihood of this).
   *
   * @param nominalPredictor if true use nominal predictor attributes
   * @param numericPredictor if true use numeric predictor attributes
   * @param stringPredictor if true use string predictor attributes
   * @param datePredictor if true use date predictor attributes
   * @param relationalPredictor if true use relational predictor attributes
   * @param multiInstance whether multi-instance is needed
   * @param classType the class type (NUMERIC, NOMINAL, etc.)
   * @return index 0 true if the test was passed
   */
  protected boolean[] instanceWeights(
      boolean nominalPredictor,
      boolean numericPredictor,
      boolean stringPredictor,
      boolean datePredictor,
      boolean relationalPredictor,
      boolean multiInstance,
      int classType) {

    print("classifier uses instance weights");
    printAttributeSummary(
        nominalPredictor, numericPredictor, stringPredictor, datePredictor, relationalPredictor, multiInstance, classType);
    print("...");
    int numTrain = 2*getNumInstances(), numTest = getNumInstances(),
    numClasses = 2, missingLevel = 0;
    boolean predictorMissing = false, classMissing = false;

    boolean[] result = new boolean[2];
    Instances train = null;
    Instances test = null;
    Classifier [] classifiers = null;
    Evaluation evaluationB = null;
    Evaluation evaluationI = null;
    boolean built = false;
    boolean evalFail = false;
    try {
      train = makeTestDataset(42, numTrain,
                              nominalPredictor    ? getNumNominal() + 1 : 0,
                              numericPredictor    ? getNumNumeric() + 1 : 0,
                              stringPredictor     ? getNumString()      : 0,
                              datePredictor       ? getNumDate()        : 0,
                              relationalPredictor ? getNumRelational()  : 0,
                              numClasses,
                              classType,
                              multiInstance);
      test = makeTestDataset(24, numTest,
                             nominalPredictor    ? getNumNominal() + 1 : 0,
                             numericPredictor    ? getNumNumeric() + 1 : 0,
                             stringPredictor     ? getNumString()      : 0,
                             datePredictor       ? getNumDate()        : 0,
                             relationalPredictor ? getNumRelational()  : 0,
                             numClasses,
                             classType,
                             multiInstance);
      if (missingLevel > 0) {
        addMissing(train, missingLevel, predictorMissing, classMissing);
        addMissing(test, Math.min(missingLevel, 50), predictorMissing,
            classMissing);
      }
      classifiers = AbstractClassifier.makeCopies(getClassifier(), 2);
      evaluationB = new Evaluation(train);
      evaluationI = new Evaluation(train);
      classifiers[0].buildClassifier(train);
      testWRTZeroR(classifiers[0], evaluationB, train, test);
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
      classifiers[1].buildClassifier(train);
      built = true;
      testWRTZeroR(classifiers[1], evaluationI, train, test);
      if (evaluationB.equals(evaluationI)) {
        // println("no");
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
          println(evaluationB.toSummaryString("\nboth methods\n",
              true));
        } else {
          print("Problem during");
          if (built) {
            print(" testing");
          } else {
            print(" training");
          }
          println(": " + ex.getMessage() + "\n");
        }
        println("Here are the datasets:\n");
        println("=== Train Dataset ===\n"
            + train.toString() + "\n");
        println("=== Train Weights ===\n");
        for (int i = 0; i < train.numInstances(); i++) {
          println(" " + (i + 1)
              + "    " + train.instance(i).weight());
        }
        println("=== Test Dataset ===\n"
            + test.toString() + "\n\n");
        println("(test weights all 1.0\n");
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
   * @param classType the class type (NUMERIC, NOMINAL, etc.)
   * @param predictorMissing true if we know the classifier can handle
   * (at least) moderate missing predictor values
   * @param classMissing true if we know the classifier can handle
   * (at least) moderate missing class values
   * @return index 0 is true if the test was passed
   */
  protected boolean[] datasetIntegrity(
      boolean nominalPredictor,
      boolean numericPredictor,
      boolean stringPredictor,
      boolean datePredictor,
      boolean relationalPredictor,
      boolean multiInstance,
      int classType,
      boolean predictorMissing,
      boolean classMissing) {

    print("classifier doesn't alter original datasets");
    printAttributeSummary(
        nominalPredictor, numericPredictor, stringPredictor, datePredictor, relationalPredictor, multiInstance, classType);
    print("...");
    int numTrain = getNumInstances(), numTest = getNumInstances(),
    numClasses = 2, missingLevel = 20;

    boolean[] result = new boolean[2];
    Instances train = null;
    Instances test = null;
    Classifier classifier = null;
    Evaluation evaluation = null;
    boolean built = false;
    try {
      train = makeTestDataset(42, numTrain,
                              nominalPredictor    ? getNumNominal()    : 0,
                              numericPredictor    ? getNumNumeric()    : 0,
                              stringPredictor     ? getNumString()     : 0,
                              datePredictor       ? getNumDate()       : 0,
                              relationalPredictor ? getNumRelational() : 0,
                              numClasses,
                              classType,
                              multiInstance);
      test = makeTestDataset(24, numTest,
                             nominalPredictor     ? getNumNominal()    : 0,
                             numericPredictor     ? getNumNumeric()    : 0,
                             stringPredictor      ? getNumString()     : 0,
                             datePredictor        ? getNumDate()       : 0,
                             relationalPredictor  ? getNumRelational() : 0,
                             numClasses,
                             classType,
                             multiInstance);
      if (missingLevel > 0) {
        addMissing(train, missingLevel, predictorMissing, classMissing);
        addMissing(test, Math.min(missingLevel, 50), predictorMissing,
            classMissing);
      }
      classifier = AbstractClassifier.makeCopies(getClassifier(), 1)[0];
      evaluation = new Evaluation(train);
    } catch (Exception ex) {
      throw new Error("Error setting up for tests: " + ex.getMessage());
    }
    try {
      Instances trainCopy = new Instances(train);
      Instances testCopy = new Instances(test);
      classifier.buildClassifier(trainCopy);
      compareDatasets(train, trainCopy);
      built = true;
      testWRTZeroR(classifier, evaluation, trainCopy, testCopy);
      compareDatasets(test, testCopy);

      println("yes");
      result[0] = true;
    } catch (Exception ex) {
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
        println(": " + ex.getMessage() + "\n");
        println("Here are the datasets:\n");
        println("=== Train Dataset ===\n"
            + train.toString() + "\n");
        println("=== Test Dataset ===\n"
            + test.toString() + "\n\n");
      }
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
   * @param classType the class type (NUMERIC, NOMINAL, etc.)
   * @param missingLevel the percentage of missing values
   * @param predictorMissing true if the missing values may be in
   * the predictors
   * @param classMissing true if the missing values may be in the class
   * @param numTrain the number of instances in the training set
   * @param numTest the number of instaces in the test set
   * @param numClasses the number of classes
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
      int classType,
      int missingLevel,
      boolean predictorMissing,
      boolean classMissing,
      int numTrain,
      int numTest,
      int numClasses,
      FastVector accepts) {

    return runBasicTest(
        nominalPredictor,
        numericPredictor,
        stringPredictor,
        datePredictor,
        relationalPredictor,
        multiInstance,
        classType,
        TestInstances.CLASS_IS_LAST,
        missingLevel,
        predictorMissing,
        classMissing,
        numTrain,
        numTest,
        numClasses,
        accepts);
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
   * @param classType the class type (NUMERIC, NOMINAL, etc.)
   * @param classIndex the attribute index of the class
   * @param missingLevel the percentage of missing values
   * @param predictorMissing true if the missing values may be in
   * the predictors
   * @param classMissing true if the missing values may be in the class
   * @param numTrain the number of instances in the training set
   * @param numTest the number of instaces in the test set
   * @param numClasses the number of classes
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
      int classType,
      int classIndex,
      int missingLevel,
      boolean predictorMissing,
      boolean classMissing,
      int numTrain,
      int numTest,
      int numClasses,
      FastVector accepts) {

    boolean[] result = new boolean[2];
    Instances train = null;
    Instances test = null;
    Classifier classifier = null;
    Evaluation evaluation = null;
    boolean built = false;
    try {
      train = makeTestDataset(42, numTrain,
                              nominalPredictor     ? getNumNominal()    : 0,
                              numericPredictor     ? getNumNumeric()    : 0,
                              stringPredictor      ? getNumString()     : 0,
                              datePredictor        ? getNumDate()       : 0,
                              relationalPredictor  ? getNumRelational() : 0,
                              numClasses,
                              classType,
                              classIndex,
                              multiInstance);
      test = makeTestDataset(24, numTest,
                             nominalPredictor     ? getNumNominal()    : 0,
                             numericPredictor     ? getNumNumeric()    : 0,
                             stringPredictor      ? getNumString()     : 0,
                             datePredictor        ? getNumDate()       : 0,
                             relationalPredictor  ? getNumRelational() : 0,
                             numClasses,
                             classType,
                             classIndex,
                             multiInstance);
      if (missingLevel > 0) {
        addMissing(train, missingLevel, predictorMissing, classMissing);
        addMissing(test, Math.min(missingLevel, 50), predictorMissing,
            classMissing);
      }
      classifier = AbstractClassifier.makeCopies(getClassifier(), 1)[0];
      evaluation = new Evaluation(train);
    } catch (Exception ex) {
      ex.printStackTrace();
      throw new Error("Error setting up for tests: " + ex.getMessage());
    }
    try {
      classifier.buildClassifier(train);
      built = true;
      if (!testWRTZeroR(classifier, evaluation, train, test)[0]) {
        result[0] = true;
        result[1] = true;
        throw new Exception("Scheme performs worse than ZeroR");
      }

      println("yes");
      result[0] = true;
    }
    catch (Exception ex) {
      boolean acceptable = false;
      String msg;
      if (ex.getMessage() == null)
        msg = "";
      else
        msg = ex.getMessage().toLowerCase();
      if (msg.indexOf("not in classpath") > -1)
        m_ClasspathProblems = true;
      if (msg.indexOf("worse than zeror") >= 0) {
        println("warning: performs worse than ZeroR");
        result[0] = true;
        result[1] = true;
      } else {
        for (int i = 0; i < accepts.size(); i++) {
          if (msg.indexOf((String)accepts.elementAt(i)) >= 0) {
            acceptable = true;
          }
        }

        println("no" + (acceptable ? " (OK error message)" : ""));
        result[1] = acceptable;
      }

      if (m_Debug) {
        println("\n=== Full Report ===");
        print("Problem during");
        if (built) {
          print(" testing");
        } else {
          print(" training");
        }
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
          println("here are the datasets:\n");
          println("=== Train Dataset ===\n"
              + train.toString() + "\n");
          println("=== Test Dataset ===\n"
              + test.toString() + "\n\n");
        }
      }
    }

    return result;
  }

  /**
   * Determine whether the scheme performs worse than ZeroR during testing
   *
   * @param classifier the pre-trained classifier
   * @param evaluation the classifier evaluation object
   * @param train the training data
   * @param test the test data
   * @return index 0 is true if the scheme performs better than ZeroR
   * @throws Exception if there was a problem during the scheme's testing
   */
  protected boolean[] testWRTZeroR(Classifier classifier,
                                   Evaluation evaluation,
                                   Instances train, Instances test)
  throws Exception {

    boolean[] result = new boolean[2];

    evaluation.evaluateModel(classifier, test);
    try {

      // Tested OK, compare with ZeroR
      Classifier zeroR = new weka.classifiers.rules.ZeroR();
      zeroR.buildClassifier(train);
      Evaluation zeroREval = new Evaluation(train);
      zeroREval.evaluateModel(zeroR, test);
      result[0] = Utils.grOrEq(zeroREval.errorRate(), evaluation.errorRate());
    }
    catch (Exception ex) {
      throw new Error("Problem determining ZeroR performance: "
          + ex.getMessage());
    }

    return result;
  }

  /**
   * Make a simple set of instances, which can later be modified
   * for use in specific tests.
   *
   * @param seed the random number seed
   * @param numInstances the number of instances to generate
   * @param numNominal the number of nominal attributes
   * @param numNumeric the number of numeric attributes
   * @param numString the number of string attributes
   * @param numDate the number of date attributes
   * @param numRelational the number of relational attributes
   * @param numClasses the number of classes (if nominal class)
   * @param classType the class type (NUMERIC, NOMINAL, etc.)
   * @param multiInstance whether the dataset should a multi-instance dataset
   * @return the test dataset
   * @throws Exception if the dataset couldn't be generated
   * @see #process(Instances)
   */
  protected Instances makeTestDataset(int seed, int numInstances,
                                      int numNominal, int numNumeric,
                                      int numString, int numDate,
                                      int numRelational,
                                      int numClasses, int classType,
                                      boolean multiInstance)
    throws Exception {

    return makeTestDataset(
        seed,
        numInstances,
        numNominal,
        numNumeric,
        numString,
        numDate,
        numRelational,
        numClasses,
        classType,
        TestInstances.CLASS_IS_LAST,
        multiInstance);
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
   * @param numClasses the number of classes (if nominal class)
   * @param classType the class type (NUMERIC, NOMINAL, etc.)
   * @param classIndex the index of the class (0-based, -1 as last)
   * @param multiInstance whether the dataset should a multi-instance dataset
   * @return the test dataset
   * @throws Exception if the dataset couldn't be generated
   * @see TestInstances#CLASS_IS_LAST
   * @see #process(Instances)
   */
  protected Instances makeTestDataset(int seed, int numInstances,
                                      int numNominal, int numNumeric,
                                      int numString, int numDate,
                                      int numRelational,
                                      int numClasses, int classType,
                                      int classIndex,
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
    dataset.setNumClasses(numClasses);
    dataset.setClassType(classType);
    dataset.setClassIndex(classIndex);
    dataset.setNumClasses(numClasses);
    dataset.setMultiInstance(multiInstance);
    dataset.setWords(getWords());
    dataset.setWordSeparators(getWordSeparators());

    return process(dataset.generate());
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
   * @param classType the class type (NUMERIC, NOMINAL, etc.)
   */
  protected void printAttributeSummary(boolean nominalPredictor,
                                       boolean numericPredictor,
                                       boolean stringPredictor,
                                       boolean datePredictor,
                                       boolean relationalPredictor,
                                       boolean multiInstance,
                                       int classType) {

    String str = "";

    if (numericPredictor)
      str += " numeric";

    if (nominalPredictor) {
      if (str.length() > 0)
        str += " &";
      str += " nominal";
    }

    if (stringPredictor) {
      if (str.length() > 0)
        str += " &";
      str += " string";
    }

    if (datePredictor) {
      if (str.length() > 0)
        str += " &";
      str += " date";
    }

    if (relationalPredictor) {
      if (str.length() > 0)
        str += " &";
      str += " relational";
    }

    str += " predictors)";

    switch (classType) {
      case Attribute.NUMERIC:
        str = " (numeric class," + str;
        break;
      //RANKING BEGIN
      case PreferenceAttribute.RANKING:
    	str = " (ranking," + str;
    	break;
      //RANKING END
      case Attribute.NOMINAL:
        str = " (nominal class," + str;
        break;
      case Attribute.STRING:
        str = " (string class," + str;
        break;
      case Attribute.DATE:
        str = " (date class," + str;
        break;
      case Attribute.RELATIONAL:
        str = " (relational class," + str;
        break;
    }

    print(str);
  }

  /**
   * Returns the revision string.
   *
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6041 $");
  }

  /**
   * Test method for this class
   *
   * @param args the commandline parameters
   */
  public static void main(String [] args) {
    runCheck(new CheckClassifier(), args);
  }
}
