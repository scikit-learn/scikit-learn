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
 * CheckAttributeSelection.java
 * Copyright (C) 2006 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.attributeSelection;

import weka.core.Attribute;
import weka.core.CheckScheme;
import weka.core.FastVector;
import weka.core.Instances;
import weka.core.MultiInstanceCapabilitiesHandler;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.SerializationHelper;
import weka.core.SerializedObject;
import weka.core.TestInstances;
import weka.core.Utils;
import weka.core.WeightedInstancesHandler;
import weka.core.labelranking.PreferenceAttribute;

import java.util.Enumeration;
import java.util.Random;
import java.util.Vector;

/**
 * Class for examining the capabilities and finding problems with 
 * attribute selection schemes. If you implement an attribute selection using 
 * the WEKA.libraries, you should run the checks on it to ensure robustness 
 * and correct operation. Passing all the tests of this object does not mean
 * bugs in the attribute selection don't exist, but this will help find some
 * common ones. <p/>
 * 
 * Typical usage: <p/>
 * <code>java weka.attributeSelection.CheckAttributeSelection -W ASscheme_name 
 * -- ASscheme_options </code><p/>
 * 
 * CheckAttributeSelection reports on the following:
 * <ul>
 *    <li> Scheme abilities 
 *      <ul>
 *         <li> Possible command line options to the scheme </li>
 *         <li> Whether the scheme can predict nominal, numeric, string, 
 *              date or relational class attributes. </li>
 *         <li> Whether the scheme can handle numeric predictor attributes </li>
 *         <li> Whether the scheme can handle nominal predictor attributes </li>
 *         <li> Whether the scheme can handle string predictor attributes </li>
 *         <li> Whether the scheme can handle date predictor attributes </li>
 *         <li> Whether the scheme can handle relational predictor attributes </li>
 *         <li> Whether the scheme can handle multi-instance data </li>
 *         <li> Whether the scheme can handle missing predictor values </li>
 *         <li> Whether the scheme can handle missing class values </li>
 *         <li> Whether a nominal scheme only handles 2 class problems </li>
 *         <li> Whether the scheme can handle instance weights </li>
 *      </ul>
 *    </li>
 *    <li> Correct functioning 
 *      <ul>
 *         <li> Correct initialisation during search (i.e. no result
 *              changes when search is performed repeatedly) </li>
 *         <li> Whether the scheme alters the data pased to it 
 *              (number of instances, instance order, instance weights, etc) </li>
 *      </ul>
 *    </li>
 *    <li> Degenerate cases 
 *      <ul>
 *         <li> building scheme with zero instances </li>
 *         <li> all but one predictor attribute values missing </li>
 *         <li> all predictor attribute values missing </li>
 *         <li> all but one class values missing </li>
 *         <li> all class values missing </li>
 *      </ul>
 *    </li>
 * </ul>
 * Running CheckAttributeSelection with the debug option set will output the 
 * training dataset for any failed tests.<p/>
 *
 * The <code>weka.attributeSelection.AbstractAttributeSelectionTest</code> 
 * uses this class to test all the schemes. Any changes here, have to be 
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
 * <pre> -eval name [options]
 *  Full name and options of the evaluator analyzed.
 *  eg: weka.attributeSelection.CfsSubsetEval</pre>
 * 
 * <pre> -search name [options]
 *  Full name and options of the search method analyzed.
 *  eg: weka.attributeSelection.Ranker</pre>
 * 
 * <pre> -test &lt;eval|search&gt;
 *  The scheme to test, either the evaluator or the search method.
 *  (Default: eval)</pre>
 * 
 * <pre> 
 * Options specific to evaluator weka.attributeSelection.CfsSubsetEval:
 * </pre>
 * 
 * <pre> -M
 *  Treat missing values as a seperate value.</pre>
 * 
 * <pre> -L
 *  Don't include locally predictive attributes.</pre>
 * 
 * <pre> 
 * Options specific to search method weka.attributeSelection.Ranker:
 * </pre>
 * 
 * <pre> -P &lt;start set&gt;
 *  Specify a starting set of attributes.
 *  Eg. 1,3,5-7.
 *  Any starting attributes specified are
 *  ignored during the ranking.</pre>
 * 
 * <pre> -T &lt;threshold&gt;
 *  Specify a theshold by which attributes
 *  may be discarded from the ranking.</pre>
 * 
 * <pre> -N &lt;num to select&gt;
 *  Specify number of attributes to select</pre>
 * 
 <!-- options-end -->
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 4783 $
 * @see TestInstances
 */
public class CheckAttributeSelection 
  extends CheckScheme {

  /*
   * Note about test methods:
   * - methods return array of booleans
   * - first index: success or not
   * - second index: acceptable or not (e.g., Exception is OK)
   *
   * FracPete (fracpete at waikato dot ac dot nz)
   */
  
  /*** The evaluator to be examined */
  protected ASEvaluation m_Evaluator = new CfsSubsetEval();
  
  /*** The search method to be used */
  protected ASSearch m_Search = new Ranker();
  
  /** whether to test the evaluator (default) or the search method */
  protected boolean m_TestEvaluator = true;
  
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
        "\tFull name and options of the evaluator analyzed.\n"
        +"\teg: weka.attributeSelection.CfsSubsetEval",
        "eval", 1, "-eval name [options]"));
    
    result.addElement(new Option(
        "\tFull name and options of the search method analyzed.\n"
        +"\teg: weka.attributeSelection.Ranker",
        "search", 1, "-search name [options]"));
    
    result.addElement(new Option(
        "\tThe scheme to test, either the evaluator or the search method.\n"
        +"\t(Default: eval)",
        "test", 1, "-test <eval|search>"));
    
    if ((m_Evaluator != null) && (m_Evaluator instanceof OptionHandler)) {
      result.addElement(new Option("", "", 0, 
          "\nOptions specific to evaluator "
          + m_Evaluator.getClass().getName()
          + ":"));
      Enumeration enm = ((OptionHandler) m_Evaluator).listOptions();
      while (enm.hasMoreElements())
        result.addElement(enm.nextElement());
    }
    
    if ((m_Search != null) && (m_Search instanceof OptionHandler)) {
      result.addElement(new Option("", "", 0, 
          "\nOptions specific to search method "
          + m_Search.getClass().getName()
          + ":"));
      Enumeration enm = ((OptionHandler) m_Search).listOptions();
      while (enm.hasMoreElements())
        result.addElement(enm.nextElement());
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
   * <pre> -eval name [options]
   *  Full name and options of the evaluator analyzed.
   *  eg: weka.attributeSelection.CfsSubsetEval</pre>
   * 
   * <pre> -search name [options]
   *  Full name and options of the search method analyzed.
   *  eg: weka.attributeSelection.Ranker</pre>
   * 
   * <pre> -test &lt;eval|search&gt;
   *  The scheme to test, either the evaluator or the search method.
   *  (Default: eval)</pre>
   * 
   * <pre> 
   * Options specific to evaluator weka.attributeSelection.CfsSubsetEval:
   * </pre>
   * 
   * <pre> -M
   *  Treat missing values as a seperate value.</pre>
   * 
   * <pre> -L
   *  Don't include locally predictive attributes.</pre>
   * 
   * <pre> 
   * Options specific to search method weka.attributeSelection.Ranker:
   * </pre>
   * 
   * <pre> -P &lt;start set&gt;
   *  Specify a starting set of attributes.
   *  Eg. 1,3,5-7.
   *  Any starting attributes specified are
   *  ignored during the ranking.</pre>
   * 
   * <pre> -T &lt;threshold&gt;
   *  Specify a theshold by which attributes
   *  may be discarded from the ranking.</pre>
   * 
   * <pre> -N &lt;num to select&gt;
   *  Specify number of attributes to select</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String      tmpStr;
    String[]	tmpOptions;
    
    super.setOptions(options);
    
    tmpStr     = Utils.getOption("eval", options);
    tmpOptions = Utils.splitOptions(tmpStr);
    if (tmpOptions.length != 0) {
      tmpStr        = tmpOptions[0];
      tmpOptions[0] = "";
      setEvaluator(
	  (ASEvaluation) forName(
	      "weka.attributeSelection", 
	      ASEvaluation.class, 
	      tmpStr, 
	      tmpOptions));
    }
    
    tmpStr     = Utils.getOption("search", options);
    tmpOptions = Utils.splitOptions(tmpStr);
    if (tmpOptions.length != 0) {
      tmpStr        = tmpOptions[0];
      tmpOptions[0] = "";
      setSearch(
	  (ASSearch) forName(
	      "weka.attributeSelection", 
	      ASSearch.class, 
	      tmpStr, 
	      tmpOptions));
    }

    tmpStr = Utils.getOption("test", options);
    setTestEvaluator(!tmpStr.equalsIgnoreCase("search"));
  }
  
  /**
   * Gets the current settings of the CheckAttributeSelection.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    Vector	result;
    String[]	options;
    int		i;
    
    result = new Vector();
    
    options = super.getOptions();
    for (i = 0; i < options.length; i++)
      result.add(options[i]);
    
    result.add("-eval");
    if (getEvaluator() instanceof OptionHandler)
      result.add(
	  getEvaluator().getClass().getName() 
	  + " " 
	  + Utils.joinOptions(((OptionHandler) getEvaluator()).getOptions()));
    else
      result.add(
	  getEvaluator().getClass().getName());

    result.add("-search");
    if (getSearch() instanceof OptionHandler)
      result.add(
	  getSearch().getClass().getName() 
	  + " " 
	  + Utils.joinOptions(((OptionHandler) getSearch()).getOptions()));
    else
      result.add(
	  getSearch().getClass().getName());
    
    result.add("-test");
    if (getTestEvaluator())
      result.add("eval");
    else
      result.add("search");
    
    return (String[]) result.toArray(new String[result.size()]);
  }
  
  /**
   * Begin the tests, reporting results to System.out
   */
  public void doTests() {
    
    if (getTestObject() == null) {
      println("\n=== No scheme set ===");
      return;
    }
    println("\n=== Check on scheme: "
        + getTestObject().getClass().getName()
        + " ===\n");
    
    // Start tests
    m_ClasspathProblems = false;
    println("--> Checking for interfaces");
    canTakeOptions();
    boolean weightedInstancesHandler = weightedInstancesHandler()[0];
    boolean multiInstanceHandler = multiInstanceHandler()[0];
    println("--> Scheme tests");
    declaresSerialVersionUID();
    testsPerClassType(PreferenceAttribute.RANKING,    weightedInstancesHandler, multiInstanceHandler);
    testsPerClassType(Attribute.NOMINAL,    weightedInstancesHandler, multiInstanceHandler);
    testsPerClassType(Attribute.NUMERIC,    weightedInstancesHandler, multiInstanceHandler);
    testsPerClassType(Attribute.DATE,       weightedInstancesHandler, multiInstanceHandler);
    testsPerClassType(Attribute.STRING,     weightedInstancesHandler, multiInstanceHandler);
    testsPerClassType(Attribute.RELATIONAL, weightedInstancesHandler, multiInstanceHandler);
  }
  
  /**
   * Set the evaluator to test. 
   *
   * @param value	the evaluator to use.
   */
  public void setEvaluator(ASEvaluation value) {
    m_Evaluator = value;
  }
  
  /**
   * Get the current evaluator
   *
   * @return 		the current evaluator
   */
  public ASEvaluation getEvaluator() {
    return m_Evaluator;
  }
  
  /**
   * Set the search method to test. 
   *
   * @param value	the search method to use.
   */
  public void setSearch(ASSearch value) {
    m_Search = value;
  }
  
  /**
   * Get the current search method
   *
   * @return 		the current search method
   */
  public ASSearch getSearch() {
    return m_Search;
  }

  /**
   * Sets whether the evaluator or the search method is being tested.
   * 
   * @param value	if true then the evaluator will be tested
   */
  public void setTestEvaluator(boolean value) {
    m_TestEvaluator = value;
  }
  
  /**
   * Gets whether the evaluator is being tested or the search method.
   * 
   * @return		true if the evaluator is being tested
   */
  public boolean getTestEvaluator() {
    return m_TestEvaluator;
  }
  
  /**
   * returns either the evaluator or the search method.
   * 
   * @return		the object to be tested
   * @see		#m_TestEvaluator
   */
  protected Object getTestObject() {
    if (getTestEvaluator())
      return getEvaluator();
    else
      return getSearch();
  }
  
  /**
   * returns deep copies of the given object
   * 
   * @param obj		the object to copy
   * @param num		the number of copies
   * @return		the deep copies
   * @throws Exception	if copying fails
   */
  protected Object[] makeCopies(Object obj, int num) throws Exception {
    if (obj == null)
      throw new Exception("No object set");

    Object[] objs = new Object[num];
    SerializedObject so = new SerializedObject(obj);
    for(int i = 0; i < objs.length; i++) {
      objs[i] = so.getObject();
    }
    
    return objs;
  }
  
  /**
   * Performs a attribute selection with the given search and evaluation scheme 
   * on the provided data. The generated AttributeSelection object is returned.
   * 
   * @param search	the search scheme to use
   * @param eval	the evaluator to use
   * @param data	the data to work on
   * @return		the used attribute selection object
   * @throws Exception	if the attribute selection fails
   */
  protected AttributeSelection search(ASSearch search, ASEvaluation eval, 
      Instances data) throws Exception {
    
    AttributeSelection	result;
    
    result = new AttributeSelection();
    result.setSeed(42);
    result.setSearch(search);
    result.setEvaluator(eval);
    result.SelectAttributes(data);
    
    return result;
  }
  
  /**
   * Run a battery of tests for a given class attribute type
   *
   * @param classType true if the class attribute should be numeric
   * @param weighted true if the scheme says it handles weights
   * @param multiInstance true if the scheme handles multi-instance data
   */
  protected void testsPerClassType(int classType, 
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
      
      correctSearchInitialisation(PNom, PNum, PStr, PDat, PRel, multiInstance, classType);
      datasetIntegrity(PNom, PNum, PStr, PDat, PRel, multiInstance, classType,
          handleMissingPredictors, handleMissingClass);
    }
  }
  
  /**
   * Checks whether the scheme can take command line options.
   *
   * @return index 0 is true if the scheme can take options
   */
  protected boolean[] canTakeOptions() {
    
    boolean[] result = new boolean[2];
    
    print("options...");
    if (getTestObject() instanceof OptionHandler) {
      println("yes");
      if (m_Debug) {
        println("\n=== Full report ===");
        Enumeration enu = ((OptionHandler) getTestObject()).listOptions();
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
   * Checks whether the scheme says it can handle instance weights.
   *
   * @return true if the scheme handles instance weights
   */
  protected boolean[] weightedInstancesHandler() {
    
    boolean[] result = new boolean[2];
    
    print("weighted instances scheme...");
    if (getTestObject() instanceof WeightedInstancesHandler) {
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
   * @return true if the scheme handles multi-instance data
   */
  protected boolean[] multiInstanceHandler() {
    boolean[] result = new boolean[2];
    
    print("multi-instance scheme...");
    if (getTestObject() instanceof MultiInstanceCapabilitiesHandler) {
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
   * tests for a serialVersionUID. Fails in case the schemes don't declare
   * a UID (both must!).
   *
   * @return index 0 is true if the scheme declares a UID
   */
  protected boolean[] declaresSerialVersionUID() {
    boolean[] result = new boolean[2];
    boolean eval;
    boolean search;
    
    print("serialVersionUID...");
    
    eval   = !SerializationHelper.needsUID(m_Evaluator.getClass());
    search = !SerializationHelper.needsUID(m_Search.getClass());
    
    result[0] = eval && search;
    
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
    int numTrain = getNumInstances(), numClasses = 2, missingLevel = 0;
    boolean predictorMissing = false, classMissing = false;
    
    return runBasicTest(nominalPredictor, numericPredictor, stringPredictor, 
        datePredictor, relationalPredictor, 
        multiInstance,
        classType, 
        missingLevel, predictorMissing, classMissing,
        numTrain, numClasses, 
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
    int numTrain = getNumInstances(), missingLevel = 0;
    boolean predictorMissing = false, classMissing = false;
    
    return runBasicTest(nominalPredictor, numericPredictor, stringPredictor, 
                        datePredictor, relationalPredictor, 
                        multiInstance,
                        Attribute.NOMINAL,
                        missingLevel, predictorMissing, classMissing,
                        numTrain, numClasses, 
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
    int numTrain = getNumInstances(), numClasses = 2, missingLevel = 0;
    boolean predictorMissing = false, classMissing = false;
    
    return runBasicTest(nominalPredictor, numericPredictor, stringPredictor, 
                        datePredictor, relationalPredictor, 
                        multiInstance,
                        classType,
                        classIndex,
                        missingLevel, predictorMissing, classMissing,
                        numTrain, numClasses, 
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
    int numTrain = 0, numClasses = 2, missingLevel = 0;
    boolean predictorMissing = false, classMissing = false;
    
    return runBasicTest(
              nominalPredictor, numericPredictor, stringPredictor, 
              datePredictor, relationalPredictor, 
              multiInstance,
              classType, 
              missingLevel, predictorMissing, classMissing,
              numTrain, numClasses, 
              accepts);
  }
  
  /**
   * Checks whether the scheme correctly initialises models when 
   * ASSearch.search is called. This test calls search with
   * one training dataset. ASSearch is then called on a training set with 
   * different structure, and then again with the original training set. 
   * If the equals method of the ASEvaluation class returns false, this is 
   * noted as incorrect search initialisation.
   *
   * @param nominalPredictor if true use nominal predictor attributes
   * @param numericPredictor if true use numeric predictor attributes
   * @param stringPredictor if true use string predictor attributes
   * @param datePredictor if true use date predictor attributes
   * @param relationalPredictor if true use relational predictor attributes
   * @param multiInstance whether multi-instance is needed
   * @param classType the class type (NUMERIC, NOMINAL, etc.)
   * @return index 0 is true if the test was passed, index 1 is always false
   */
  protected boolean[] correctSearchInitialisation(
      boolean nominalPredictor,
      boolean numericPredictor, 
      boolean stringPredictor, 
      boolean datePredictor,
      boolean relationalPredictor,
      boolean multiInstance,
      int classType) {

    boolean[] result = new boolean[2];
    print("correct initialisation during search");
    printAttributeSummary(
        nominalPredictor, numericPredictor, stringPredictor, datePredictor, relationalPredictor, multiInstance, classType);
    print("...");
    int numTrain = getNumInstances(), 
    numClasses = 2, missingLevel = 0;
    boolean predictorMissing = false, classMissing = false;
    
    Instances train1 = null;
    Instances train2 = null;
    ASSearch search = null;
    ASEvaluation evaluation1A = null;
    ASEvaluation evaluation1B = null;
    ASEvaluation evaluation2 = null;
    AttributeSelection attsel1A = null;
    AttributeSelection attsel1B = null;
    int stage = 0;
    try {
      
      // Make two train sets with different numbers of attributes
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
      if (missingLevel > 0) {
        addMissing(train1, missingLevel, predictorMissing, classMissing);
        addMissing(train2, missingLevel, predictorMissing, classMissing);
      }
      
      search = ASSearch.makeCopies(getSearch(), 1)[0];
      evaluation1A = ASEvaluation.makeCopies(getEvaluator(), 1)[0];
      evaluation1B = ASEvaluation.makeCopies(getEvaluator(), 1)[0];
      evaluation2 = ASEvaluation.makeCopies(getEvaluator(), 1)[0];
    } catch (Exception ex) {
      throw new Error("Error setting up for tests: " + ex.getMessage());
    }
    try {
      stage = 0;
      attsel1A = search(search, evaluation1A, train1);
      
      stage = 1;
      search(search, evaluation2, train2);
      
      stage = 2;
      attsel1B = search(search, evaluation1B, train1);
      
      stage = 3;
      if (!attsel1A.toResultsString().equals(attsel1B.toResultsString())) {
        if (m_Debug) {
          println(
              "\n=== Full report ===\n"
              + "\nFirst search\n"
              + attsel1A.toResultsString()
              + "\n\n");
          println(
              "\nSecond search\n"
              + attsel1B.toResultsString()
              + "\n\n");
        }
        throw new Exception("Results differ between search calls");
      }
      println("yes");
      result[0] = true;
      
      if (false && m_Debug) {
        println(
            "\n=== Full report ===\n"
            + "\nFirst search\n"
            + evaluation1A.toString()
            + "\n\n");
        println(
            "\nSecond search\n"
            + evaluation1B.toString()
            + "\n\n");
      }
    } 
    catch (Exception ex) {
      println("no");
      result[0] = false;
      if (m_Debug) {
        println("\n=== Full Report ===");
        print("Problem during  training");
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
    accepts.addElement("no attributes");
    int numTrain = getNumInstances(), numClasses = 2;
    
    return runBasicTest(nominalPredictor, numericPredictor, stringPredictor, 
        datePredictor, relationalPredictor, 
        multiInstance,
        classType, 
        missingLevel, predictorMissing, classMissing,
        numTrain, numClasses, 
        accepts);
  }
  
  /**
   * Checks whether the scheme can handle instance weights.
   * This test compares the scheme performance on two datasets
   * that are identical except for the training weights. If the 
   * results change, then the scheme must be using the weights. It
   * may be possible to get a false positive from this test if the 
   * weight changes aren't significant enough to induce a change
   * in scheme performance (but the weights are chosen to minimize
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
    
    print("scheme uses instance weights");
    printAttributeSummary(
        nominalPredictor, numericPredictor, stringPredictor, datePredictor, relationalPredictor, multiInstance, classType);
    print("...");
    int numTrain = 2*getNumInstances(), 
    numClasses = 2, missingLevel = 0;
    boolean predictorMissing = false, classMissing = false;
    
    boolean[] result = new boolean[2];
    Instances train = null;
    ASSearch[] search = null;
    ASEvaluation evaluationB = null;
    ASEvaluation evaluationI = null;
    AttributeSelection attselB = null;
    AttributeSelection attselI = null;
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
      if (missingLevel > 0)
        addMissing(train, missingLevel, predictorMissing, classMissing);
      search = ASSearch.makeCopies(getSearch(), 2);
      evaluationB = ASEvaluation.makeCopies(getEvaluator(), 1)[0];
      evaluationI = ASEvaluation.makeCopies(getEvaluator(), 1)[0];
      attselB = search(search[0], evaluationB, train);
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
      attselI = search(search[1], evaluationI, train);
      if (attselB.toResultsString().equals(attselI.toResultsString())) {
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
          println(evaluationB.toString());
        } else {
          print("Problem during training");
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
   * @param classType the class type (NUMERIC, NOMINAL, etc.)
   * @param predictorMissing true if we know the scheme can handle
   * (at least) moderate missing predictor values
   * @param classMissing true if we know the scheme can handle
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
    
    print("scheme doesn't alter original datasets");
    printAttributeSummary(
        nominalPredictor, numericPredictor, stringPredictor, datePredictor, relationalPredictor, multiInstance, classType);
    print("...");
    int numTrain = getNumInstances(), 
    numClasses = 2, missingLevel = 20;
    
    boolean[] result = new boolean[2];
    Instances train = null;
    Instances trainCopy = null;
    ASSearch search = null;
    ASEvaluation evaluation = null;
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
      if (missingLevel > 0)
        addMissing(train, missingLevel, predictorMissing, classMissing);
      search = ASSearch.makeCopies(getSearch(), 1)[0];
      evaluation = ASEvaluation.makeCopies(getEvaluator(), 1)[0];
      trainCopy = new Instances(train);
    } catch (Exception ex) {
      throw new Error("Error setting up for tests: " + ex.getMessage());
    }
    try {
      search(search, evaluation, trainCopy);
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
        println("Here are the datasets:\n");
        println("=== Train Dataset (original) ===\n"
            + trainCopy.toString() + "\n");
        println("=== Train Dataset ===\n"
            + train.toString() + "\n");
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
      int numClasses,
      FastVector accepts) {
    
    boolean[] result = new boolean[2];
    Instances train = null;
    ASSearch search = null;
    ASEvaluation evaluation = null;
    try {
      train = makeTestDataset(42, numTrain, 
                              nominalPredictor    ? getNumNominal()    : 0,
                              numericPredictor    ? getNumNumeric()    : 0, 
                              stringPredictor     ? getNumString()     : 0,
                              datePredictor       ? getNumDate()       : 0,
                              relationalPredictor ? getNumRelational() : 0,
                              numClasses, 
                              classType,
                              classIndex,
                              multiInstance);
      if (missingLevel > 0)
        addMissing(train, missingLevel, predictorMissing, classMissing);
      search = ASSearch.makeCopies(getSearch(), 1)[0];
      evaluation = ASEvaluation.makeCopies(getEvaluator(), 1)[0];
    } catch (Exception ex) {
      ex.printStackTrace();
      throw new Error("Error setting up for tests: " + ex.getMessage());
    }
    try {
      search(search, evaluation, train);
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
    	str = " (ranking," +str;
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
    return RevisionUtils.extract("$Revision: 4783 $");
  }
  
  /**
   * Test method for this class
   * 
   * @param args the commandline parameters
   */
  public static void main(String [] args) {
    runCheck(new CheckAttributeSelection(), args);
  }
}

