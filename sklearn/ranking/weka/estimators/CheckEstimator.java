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
 *    CheckEstimator.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.estimators;

import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.core.TestInstances;
import weka.core.Utils;
import weka.core.WeightedInstancesHandler;
import weka.core.labelranking.PreferenceAttribute;

import java.util.Enumeration;
import java.util.Random;
import java.util.Vector;

/**
 * Class for examining the capabilities and finding problems with 
 * estimators. If you implement a estimator using the WEKA.libraries,
 * you should run the checks on it to ensure robustness and correct
 * operation. Passing all the tests of this object does not mean
 * bugs in the estimator don't exist, but this will help find some
 * common ones. <p/>
 * 
 * Typical usage: <p/>
 * <code>java weka.estimators.CheckEstimator -W estimator_name 
 * estimator_options </code><p/>
 * 
 * This class uses code from the CheckEstimatorClass
 * ATTENTION! Current estimators can only 
 * 1. split on a nominal class attribute
 * 2. build estimators for nominal and numeric attributes
 * 3. build estimators independendly of the class type
 * The functionality to test on other class and attribute types
 * is left in big parts in the code. 
 * 
 * CheckEstimator reports on the following:
 * <ul>
 *    <li> Estimator abilities 
 *      <ul>
 *         <li> Possible command line options to the estimator </li>
 *         <li> Whether the estimator can predict nominal, numeric, string, 
 *              date or relational class attributes. Warnings will be displayed if 
 *              performance is worse than ZeroR </li>
 *         <li> Whether the estimator can be trained incrementally </li>
 *         <li> Whether the estimator can build estimates for numeric attributes </li>
 *         <li> Whether the estimator can handle nominal attributes </li>
 *         <li> Whether the estimator can handle string attributes </li>
 *         <li> Whether the estimator can handle date attributes </li>
 *         <li> Whether the estimator can handle relational  attributes </li>
 *         <li> Whether the estimator build estimates for multi-instance data </li>
 *         <li> Whether the estimator can handle missing attribute values </li>
 *         <li> Whether the estimator can handle missing class values </li>
 *         <li> Whether a nominal estimator only handles 2 class problems </li>
 *         <li> Whether the estimator can handle instance weights </li>
 *      </ul>
 *    </li>
 *    <li> Correct functioning 
 *      <ul>
 *         <li> Correct initialisation during addvalues (i.e. no result
 *              changes when addValues called repeatedly) </li>
 *         <li> Whether incremental training produces the same results
 *              as during non-incremental training (which may or may not 
 *              be OK) </li>
 *         <li> Whether the estimator alters the data pased to it 
 *              (number of instances, instance order, instance weights, etc) </li>
 *      </ul>
 *    </li>
 *    <li> Degenerate cases 
 *      <ul>
 *         <li> building estimator with zero training instances </li>
 *         <li> all but one attribute attribute values missing </li>
 *         <li> all attribute attribute values missing </li>
 *         <li> all but one class values missing </li>
 *         <li> all class values missing </li>
 *      </ul>
 *    </li>
 * </ul>
 * Running CheckEstimator with the debug option set will output the 
 * training and test datasets for any failed tests.<p/>
 *
 * The <code>weka.estimators.AbstractEstimatorTest</code> uses this
 * class to test all the estimators. Any changes here, have to be 
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
 *  The number of instances in the datasets (default 100).</pre>
 * 
 * <pre> -W
 *  Full name of the estimator analysed.
 *  eg: weka.estimators.NormalEstimator</pre>
 * 
 * <pre> 
 * Options specific to estimator weka.estimators.NormalEstimator:
 * </pre>
 * 
 * <pre> -D
 *  If set, estimator is run in debug mode and
 *  may output additional info to the console</pre>
 * 
 <!-- options-end -->
 *
 * Options after -- are passed to the designated estimator.<p/>
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 4997 $
 * @see TestInstances
 */
public class CheckEstimator implements OptionHandler, RevisionHandler {

  /*
   * Note about test methods:
   * - methods return array of booleans
   * - first index: success or not
   * - second index: acceptable or not (e.g., Exception is OK)
   * - in case the performance is worse than that of ZeroR both indices are true
   *
   * FracPete (fracpete at waikato dot ac dot nz)
   */
  
  /** a class for postprocessing the test-data 
   */
  public class PostProcessor
    implements RevisionHandler {
    /**
     * Provides a hook for derived classes to further modify the data. Currently,
     * the data is just passed through.
     * 
     * @param data	the data to process
     * @return		the processed data
     */
    protected Instances process(Instances data) {
      return data;
    }
    
    /**
     * Returns the revision string.
     * 
     * @return		the revision
     */
    public String getRevision() {
      return RevisionUtils.extract("$Revision: 4997 $");
    }
  }
  
  /*** The estimator to be examined */
  protected Estimator m_Estimator = (Estimator) new weka.estimators.NormalEstimator(0.000001);
  
  /** The options to be passed to the base estimator. */
  protected String[] m_EstimatorOptions;
  
  /** The results of the analysis as a string */
  protected String m_AnalysisResults;
  
  /** Debugging mode, gives extra output if true */
  protected boolean m_Debug = false;
  
  /** Silent mode, for no output at all to stdout */
  protected boolean m_Silent = false;
  
  /** The number of instances in the datasets */
  protected int m_NumInstances = 100;
  
  /** for post-processing the data even further */
  protected PostProcessor m_PostProcessor = null;
  
  /** whether classpath problems occurred */
  protected boolean m_ClasspathProblems = false;
  
  /**
   * class that contains info about the attribute types the estimator can estimate
   * estimator work on one attribute only
   */
  public static class AttrTypes
    implements RevisionHandler {
    
    boolean nominal = false;
    boolean numeric = false; 
    boolean string = false;
    boolean date = false;
    boolean relational = false;
    boolean ranking = false;
	
    AttrTypes() {
    }

    AttrTypes (AttrTypes newTypes) {
      nominal = newTypes.nominal;
      numeric = newTypes.numeric;
      string = newTypes.string;
      date = newTypes.date;
      relational = newTypes.relational;
      ranking = newTypes.ranking;
    }
			
    AttrTypes (int type) {
      if (type == Attribute.NOMINAL) nominal = true;
      if (type == Attribute.NUMERIC) numeric = true;
      if (type == Attribute.STRING) string = true;
      if (type == Attribute.DATE) date = true;
      if (type == Attribute.RELATIONAL) relational = true;
      if (type == PreferenceAttribute.RANKING) ranking = true;
    }

    int getSetType() throws Exception {			
      int sum = 0;
      int type = -1;
      if (nominal) { sum ++; type = Attribute.NOMINAL; }
      if (ranking) { sum ++; type = PreferenceAttribute.RANKING;}
      if (numeric) { sum ++; type = Attribute.NUMERIC; }
      if (string) { sum ++; type = Attribute.STRING; }
      if (date) { sum ++; type = Attribute.DATE; }
      if (relational) { sum ++; type = Attribute.RELATIONAL; }
      if (sum > 1)
	throw new Exception("Expected to have only one type set used wrongly.");
      if (type < 0)
	throw new Exception("No type set.");
      return type;
    }

    boolean oneIsSet() {
      return (ranking || nominal || numeric || string || date || relational);
    }

    public Vector getVectorOfAttrTypes() {
      Vector attrs = new Vector();
      if (ranking) attrs.add(new Integer(PreferenceAttribute.RANKING));
      if (nominal) attrs.add(new Integer(Attribute.NOMINAL));
      if (numeric) attrs.add(new Integer(Attribute.NUMERIC));
      if (string) attrs.add(new Integer(Attribute.STRING));
      if (date) attrs.add(new Integer(Attribute.DATE));
      if (relational) attrs.add(new Integer(Attribute.RELATIONAL));
      return attrs;
    }   
    
    /**
     * Returns the revision string.
     * 
     * @return		the revision
     */
    public String getRevision() {
      return RevisionUtils.extract("$Revision: 4997 $");
    }
  }

  /**
   * public class that contains info about the chosen attribute type
   * estimator work on one attribute only
   */
  public static class EstTypes
    implements RevisionHandler {
    
    boolean incremental = false;
    boolean weighted = false;
    boolean supervised = false;

    /**
     * Constructor
     */
    public EstTypes () {
    }

    /**
     * Constructor
     */
    public EstTypes (boolean i, boolean w, boolean s) {
      incremental = i;
      weighted    = w;
      supervised  = s;
    }
    
    /**
     * Returns the revision string.
     * 
     * @return		the revision
     */
    public String getRevision() {
      return RevisionUtils.extract("$Revision: 4997 $");
    }
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    
    Vector newVector = new Vector(2);
    
    newVector.addElement(new Option(
        "\tTurn on debugging output.",
        "D", 0, "-D"));
    
    newVector.addElement(new Option(
        "\tSilent mode - prints nothing to stdout.",
        "S", 0, "-S"));
    
    newVector.addElement(new Option(
        "\tThe number of instances in the datasets (default 100).",
        "N", 1, "-N <num>"));
    
    newVector.addElement(new Option(
        "\tFull name of the estimator analysed.\n"
        +"\teg: weka.estimators.NormalEstimator",
        "W", 1, "-W"));
    
    if ((m_Estimator != null) 
        && (m_Estimator instanceof OptionHandler)) {
      newVector.addElement(new Option("", "", 0, 
          "\nOptions specific to estimator "
          + m_Estimator.getClass().getName()
          + ":"));
      Enumeration enu = ((OptionHandler)m_Estimator).listOptions();
      while (enu.hasMoreElements())
        newVector.addElement(enu.nextElement());
    }
    
    return newVector.elements();
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
   *  The number of instances in the datasets (default 100).</pre>
   * 
   * <pre> -W
   *  Full name of the estimator analysed.
   *  eg: weka.estimators.NormalEstimator</pre>
   * 
   * <pre> 
   * Options specific to estimator weka.estimators.NormalEstimator:
   * </pre>
   * 
   * <pre> -D
   *  If set, estimator is run in debug mode and
   *  may output additional info to the console</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String      tmpStr;
    
    setDebug(Utils.getFlag('D', options));
    
    setSilent(Utils.getFlag('S', options));
    
    tmpStr = Utils.getOption('N', options);
    if (tmpStr.length() != 0)
      setNumInstances(Integer.parseInt(tmpStr));
    else
      setNumInstances(100);
    
    tmpStr = Utils.getOption('W', options);
    if (tmpStr.length() == 0)
      throw new Exception("A estimator must be specified with the -W option.");
    setEstimator(Estimator.forName(tmpStr, Utils.partitionOptions(options)));
  }
  
  /**
   * Gets the current settings of the CheckEstimator.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    Vector        result;
    String[]      options;
    int           i;
    
    result = new Vector();
    
    if (getDebug())
      result.add("-D");
    
    if (getSilent())
      result.add("-S");
    
    result.add("-N");
    result.add("" + getNumInstances());
    
    if (getEstimator() != null) {
      result.add("-W");
      result.add(getEstimator().getClass().getName());
    }
    
    if ((m_Estimator != null) && (m_Estimator instanceof OptionHandler))
      options = ((OptionHandler) m_Estimator).getOptions();
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
   * sets the PostProcessor to use
   * 
   * @param value	the new PostProcessor
   * @see #m_PostProcessor
   */
  public void setPostProcessor(PostProcessor value) {
    m_PostProcessor = value;
  }
  
  /**
   * returns the current PostProcessor, can be null
   * 
   * @return		the current PostProcessor
   */
  public PostProcessor getPostProcessor() {
    return m_PostProcessor;
  }
  
  /**
   * returns TRUE if the estimator returned a "not in classpath" Exception
   * 
   * @return	true if CLASSPATH problems occurred
   */
  public boolean hasClasspathProblems() {
    return m_ClasspathProblems;
  }
  
  /**
   * Begin the tests, reporting results to System.out
   */
  public void doTests() {
    
    if (getEstimator() == null) {
      println("\n=== No estimator set ===");
      return;
    }
    println("\n=== Check on Estimator: "
        + getEstimator().getClass().getName()
        + " ===\n");
    
    m_ClasspathProblems = false;

    // Start tests with test for options
    canTakeOptions();

    // test what type of estimator it is 
    EstTypes estTypes = new EstTypes();
    estTypes.incremental = incrementalEstimator()[0];
    estTypes.weighted = weightedInstancesHandler()[0];
    estTypes.supervised = supervisedEstimator()[0];
   
    // in none of the estimators yet the functionality is depending on the class type
    // since this could change the basic structure taken from checkclassifiers is kept here
    int classType = Attribute.NOMINAL;
    AttrTypes attrTypes = testsPerClassType(classType, estTypes);
    
 
    // only nominal class can be split up so far
    canSplitUpClass(attrTypes, classType);
 }
  
  
  /**
   * Set debugging mode
   *
   * @param debug true if debug output should be printed
   */
  public void setDebug(boolean debug) {
    m_Debug = debug;

    // disable silent mode, if necessary
    if (getDebug())
      setSilent(false);
  }
  
  /**
   * Get whether debugging is turned on
   *
   * @return true if debugging output is on
   */
  public boolean getDebug() {
    return m_Debug;
  }
  
  /**
   * Set slient mode, i.e., no output at all to stdout
   *
   * @param value whether silent mode is active or not
   */
  public void setSilent(boolean value) {
    m_Silent = value;
  }
  
  /**
   * Get whether silent mode is turned on
   *
   * @return true if silent mode is on
   */
  public boolean getSilent() {
    return m_Silent;
  }
  
  /**
   * Sets the number of instances to use in the datasets (some estimators
   * might require more instances).
   *
   * @param value the number of instances to use
   */
  public void setNumInstances(int value) {
    m_NumInstances = value;
  }
  
  /**
   * Gets the current number of instances to use for the datasets.
   *
   * @return the number of instances
   */
  public int getNumInstances() {
    return m_NumInstances;
  }
  
  /**
   * Set the estimator for boosting. 
   *
   * @param newEstimator the Estimator to use.
   */
  public void setEstimator(Estimator newEstimator) {
    m_Estimator = newEstimator;
  }
  
  /**
   * Get the estimator used as the estimator
   *
   * @return the estimator used as the estimator
   */
  public Estimator getEstimator() {
    return m_Estimator;
  }
  
  /**
   * prints the given message to stdout, if not silent mode
   * 
   * @param msg         the text to print to stdout
   */
  protected void print(Object msg) {
    if (!getSilent())
      System.out.print(msg);
  }
  
  /**
   * prints the given message (+ LF) to stdout, if not silent mode
   * 
   * @param msg         the message to println to stdout
   */
  protected void println(Object msg) {
    print(msg + "\n");
  }
  
  /**
   * prints a LF to stdout, if not silent mode
   */
  protected void println() {
    print("\n");
  }
  
  /**
   * Run a battery of tests for a given class attribute type
   *
   * @param classType true if the class attribute should be numeric
   * @param estTypes types the estimator is, like incremental, weighted, supervised etc
   * @return attribute types estimator can work with
   */
  protected AttrTypes testsPerClassType(int classType, EstTypes estTypes) {
    
    // in none of the estimators yet is the estimation depending on the class type
    // since this could change the basic structure taken from checkclassifiers is kept here
    
    // test A: simple test - if can estimate
    AttrTypes attrTypes = new AttrTypes();
    AttrTypes at = new AttrTypes(Attribute.NOMINAL);
    attrTypes.nominal = canEstimate(at, estTypes.supervised, classType)[0];
    at = new AttrTypes(Attribute.NUMERIC);
    attrTypes.numeric = canEstimate(at, estTypes.supervised, classType)[0];
    attrTypes.string = false;
    attrTypes.date = false;
    attrTypes.relational = false;
    
//  if (!multiInstance)
//  PRel = canEstimate(false, false, false, false,  true, classType)[0];
//  else
//  PRel = false;
    
//  one of the attribute types succeeded
    
    if (attrTypes.oneIsSet()) {
      Vector attributesSet = attrTypes.getVectorOfAttrTypes();
      
      // make tests for each attribute
      for (int i = 0; i < attributesSet.size(); i++) {
        AttrTypes workAttrTypes = new AttrTypes(((Integer) attributesSet.elementAt(i)).intValue());
        
        // test B: weights change estimate or not
        if (estTypes.weighted)
          instanceWeights(workAttrTypes, classType);
        
        if (classType == Attribute.NOMINAL || classType == PreferenceAttribute.RANKING) {
          int numClasses = 4;
          canHandleNClasses(workAttrTypes, numClasses);
        }
        
        // tests with class not the last attribute and the attribute not the first
        
        //   if (!multiInstance) {
        int numAtt = 4; 
        
        canHandleClassAsNthAttribute(workAttrTypes, numAtt, 0, classType, 1);
        
        //TODOTODOcanHandleAttrAsNthAttribute(workAttrTypes, numAtt, 2, classType);
        //}
        
        canHandleZeroTraining(workAttrTypes, classType);
        boolean handleMissingAttributes = canHandleMissing(workAttrTypes, 
            classType, true, false, 20)[0];
        if (handleMissingAttributes)
          canHandleMissing(workAttrTypes, classType, true, false, 100);
        
        boolean handleMissingClass = canHandleMissing(workAttrTypes, 
            classType, 
            false, true, 20)[0];
        if (handleMissingClass)
          canHandleMissing(workAttrTypes, classType, false, true, 100);
        
        correctBuildInitialisation(workAttrTypes, classType);
        datasetIntegrity(workAttrTypes, classType,
            handleMissingAttributes, handleMissingClass);
        
        if (estTypes.incremental)
          incrementingEquality(workAttrTypes, classType);
      }
    }
    return attrTypes;
  }
  
  /**
   * Checks whether the scheme can take command line options.
   *
   * @return index 0 is true if the estimator can take options
   */
  protected boolean[] canTakeOptions() {
    
    boolean[] result = new boolean[2];
    
    print("options...");
    if (m_Estimator instanceof OptionHandler) {
      println("yes");
      if (m_Debug) {
        println("\n=== Full report ===");
        Enumeration enu = ((OptionHandler)m_Estimator).listOptions();
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
   * @return index 0 is true if the estimator can train incrementally
   */
  protected boolean[] incrementalEstimator() {
    
    boolean[] result = new boolean[2];
    
    print("incremental estimator...");
    if (m_Estimator instanceof IncrementalEstimator) {
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
   * @return true if the estimator handles instance weights
   */
  protected boolean[] weightedInstancesHandler() {
    
    boolean[] result = new boolean[2];
    
    print("weighted instances estimator...");
    if (m_Estimator instanceof WeightedInstancesHandler) {
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
   * Checks whether the estimator is supervised.
   *
   * @return true if the estimator handles instance weights
   */
  protected boolean[] supervisedEstimator() {
    boolean[] result = new boolean[2];
    result[0] = false;
    return result;
  }

  /**
   * Checks basic estimation of one attribute of the scheme, for simple non-troublesome
   * datasets.
   *
   * @param attrTypes the types the estimator can work with
   * @param classType the class type (NOMINAL, NUMERIC, etc.)
   * @return index 0 is true if the test was passed, index 1 is true if test 
   *         was acceptable
   */
  protected boolean[] canEstimate(AttrTypes attrTypes, boolean supervised, int classType) {
    
  // supervised is ignored, no supervised estimators used yet
    
    print("basic estimation");
    printAttributeSummary(attrTypes, classType);
    print("...");
    FastVector accepts = new FastVector();
    accepts.addElement("nominal");
    accepts.addElement("numeric");
    accepts.addElement("string");
    accepts.addElement("date");
    accepts.addElement("relational");
    accepts.addElement("not in classpath");
    int numTrain = getNumInstances(), numTest = getNumInstances(), 
    numClasses = 2, missingLevel = 0;
    boolean attributeMissing = false, classMissing = false;
    int numAtts = 1, attrIndex = 0;

    return runBasicTest(attrTypes, numAtts, attrIndex,
			classType, 
			missingLevel, attributeMissing, classMissing,
			numTrain, numTest, numClasses, 
			accepts);
  }
  
  /**
   * Checks basic estimation of one attribute of the scheme, for simple non-troublesome
   * datasets.
   *
   * @param attrTypes the types the estimator can work with
   * @param classType the class type (NOMINAL, NUMERIC, etc.)
    */
  protected void canSplitUpClass(AttrTypes attrTypes, int classType) {
    
    if (attrTypes.nominal)
      canSplitUpClass(Attribute.NOMINAL, classType);
    if (attrTypes.ranking)
        canSplitUpClass(PreferenceAttribute.RANKING, classType);
    if (attrTypes.numeric)
      canSplitUpClass(Attribute.NUMERIC, classType);
  }
  
  /**
   * Checks basic estimation of one attribute of the scheme, for simple non-troublesome
   * datasets.
   *
   * @param attrType the type of the estimator
   * @param classType the class type (NOMINAL, NUMERIC, etc.)
   * @return index 0 is true if the test was passed, index 1 is true if test 
   *         was acceptable
   */
  protected boolean[] canSplitUpClass(int attrType, int classType) {
    
    boolean[] result = new boolean[2];

    FastVector accepts = new FastVector();
    accepts.addElement("not in classpath");

    // supervised is ignored, no supervised estimators used yet
    print("split per class type ");
    printAttributeSummary(attrType, Attribute.NOMINAL);
    print("...");
      
    int numTrain = getNumInstances(), numTest = getNumInstances(), 
    numClasses = 2;
    boolean attributeMissing = false, classMissing = false;
    int numAtts = 3, attrIndex = 0, classIndex = 1;
    Instances train = null;
    Vector test;
    Estimator estimator = null;
    boolean built = false;
    
    try {
      AttrTypes at = new AttrTypes(attrType);
      train = makeTestDataset(42, numTrain, numAtts, at,
          numClasses, classType, classIndex);
      
       // prepare training data set and test value list
      test = makeTestValueList(24, numTest, train, attrIndex,
          attrType);
      
       estimator = Estimator.makeCopies(getEstimator(), 1)[0];
    } catch (Exception ex) {
      ex.printStackTrace();
      throw new Error("Error setting up for tests: " + ex.getMessage());
    }
    try {
      estimator.addValues(train, attrIndex, classType, classIndex);
      built = true;
      
      testWithTestValues(estimator, test);
      
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
   * Checks whether nominal schemes can handle more than two classes.
   * If a scheme is only designed for two-class problems it should
   * throw an appropriate exception for multi-class problems.
   *
   * @param attrTypes attribute types the estimator excepts 
   * @param numClasses the number of classes to test
   * @return index 0 is true if the test was passed, index 1 is true if test 
   *         was acceptable
   */
  protected boolean[] canHandleNClasses(AttrTypes attrTypes, int numClasses) {
    
    print("more than two class problems");
    printAttributeSummary(attrTypes, Attribute.NOMINAL);
    print("...");

    FastVector accepts = new FastVector();
    accepts.addElement("number");
    accepts.addElement("class");

    int numTrain = getNumInstances(), numTest = getNumInstances(), 
      missingLevel = 0;
    boolean attributeMissing = false, classMissing = false;
    int numAttr = 1, attrIndex = 0;

    return runBasicTest(attrTypes,
                        numAttr, attrIndex,
                        Attribute.NOMINAL,
                        missingLevel, attributeMissing, classMissing,
                        numTrain, numTest, numClasses, 
                        accepts);
  }
  
  /**
   * Checks whether the scheme can handle class attributes as Nth attribute.
   *
   * @param attrTypes the attribute types the estimator accepts
   * @param numAtts of attributes
   * @param attrIndex the index of the attribute
   * @param classType the class type (NUMERIC, NOMINAL, etc.)
   * @param classIndex the index of the class attribute (0-based, -1 means last attribute)
   * @return index 0 is true if the test was passed, index 1 is true if test 
   *         was acceptable
   * @see TestInstances#CLASS_IS_LAST
   */
  protected boolean[] canHandleClassAsNthAttribute(AttrTypes attrTypes,
						   int numAtts,
						   int attrIndex,
						   int classType,
						   int classIndex) {
    
    if (classIndex == TestInstances.CLASS_IS_LAST)
      print("class attribute as last attribute");
    else
      print("class attribute as " + (classIndex + 1) + ". attribute");
    printAttributeSummary(attrTypes, classType);
    print("...");
    FastVector accepts = new FastVector();
    int numTrain = getNumInstances(), numTest = getNumInstances(), numClasses = 2, 
    missingLevel = 0;
    boolean attributeMissing = false, classMissing = false;
    
    return runBasicTest(attrTypes,
			numAtts, attrIndex,
                        classType, classIndex,
                        missingLevel, attributeMissing, classMissing,
                        numTrain, numTest, numClasses, 
                        accepts);
  }
  
  /**
   * Checks whether the scheme can handle zero training instances.
   *
   * @param attrTypes attribute types that can be estimated
   * @param classType the class type (NUMERIC, NOMINAL, etc.)
   * @return index 0 is true if the test was passed, index 1 is true if test 
   *         was acceptable
   */
  protected boolean[] canHandleZeroTraining(AttrTypes attrTypes, int classType) {
    
    print("handle zero training instances");
    printAttributeSummary(attrTypes, classType);

    print("...");
    FastVector accepts = new FastVector();
    accepts.addElement("train");
    accepts.addElement("value");
    int numTrain = 0, numTest = getNumInstances(), numClasses = 2, 
    missingLevel = 0;
    boolean attributeMissing = false, classMissing = false;
    int numAtts = 1;
    int attrIndex = 0;
    return runBasicTest(
              attrTypes, numAtts, attrIndex,
              classType, 
              missingLevel, attributeMissing, classMissing,
              numTrain, numTest, numClasses, 
              accepts);
  }
  
  /**
   * Checks whether the scheme correctly initialises models when 
   * buildEstimator is called. This test calls buildEstimator with
   * one training dataset and records performance on a test set. 
   * buildEstimator is then called on a training set with different
   * structure, and then again with the original training set. The
   * performance on the test set is compared with the original results
   * and any performance difference noted as incorrect build initialisation.
   *
   * @param attrTypes attribute types that can be estimated
   * @param classType the class type (NUMERIC, NOMINAL, etc.)
   * @return index 0 is true if the test was passed, index 1 is true if the
   *         scheme performs worse than ZeroR, but without error (index 0 is
   *         false)
   */
  protected boolean[] correctBuildInitialisation(AttrTypes attrTypes,
						 int classType) {

    boolean[] result = new boolean[2];
    
    print("correct initialisation during buildEstimator");
    printAttributeSummary(attrTypes, classType);

    print("...");
    int numTrain = getNumInstances(), numTest = getNumInstances(), 
    numClasses = 2, missingLevel = 0;
    boolean attributeMissing = false, classMissing = false;
    
    Instances train1 = null;
    Instances test1 = null;
    Instances train2 = null;
    Instances test2 = null;
    Estimator estimator = null;
    Estimator estimator1 = null;
    
    boolean built = false;
    int stage = 0;
    int attrIndex1 = 1;
    int attrIndex2 = 2;

    try {
      
      // Make two sets of train/test splits with different 
      // numbers of attributes
      train1 = makeTestDataset(42, numTrain, 2, attrTypes,
                               numClasses, 
                               classType);
      train2 = makeTestDataset(84, numTrain, 3, attrTypes,
                               numClasses, 
                               classType);
      if (missingLevel > 0) {
        addMissing(train1, missingLevel, attributeMissing, classMissing, attrIndex1);
        addMissing(train2, missingLevel, attributeMissing, classMissing, attrIndex2);
      }
      
      estimator = Estimator.makeCopies(getEstimator(), 1)[0];
    } catch (Exception ex) {
      throw new Error("Error setting up for tests: " + ex.getMessage());
    }
    try {
      //TESTING??
      stage = 0;
      estimator.addValues(train1, attrIndex1);
      built = true;

      estimator1 = estimator.makeCopies(getEstimator(), 1)[0];
      
      stage = 1;
      built = false;
      estimator.addValues(train2, attrIndex2);
      built = true;
       
      stage = 2;
      built = false;
      estimator.addValues(train1, attrIndex1);
      built = true;
      
      stage = 3;
      if (!estimator.equals(estimator1)) {
        if (m_Debug) {
          println("\n=== Full report ===\n"
		  + "\nFirst build estimator\n"+
                  estimator.toString() + "\n\n");
          println("\nSecond build estimator\n"+
		  estimator.toString() + "\n\n");
	}
        throw new Exception("Results differ between buildEstimator calls");
      }
      println("yes");
      result[0] = true;
      
      if (false && m_Debug) {
        println("\n=== Full report ===\n"
		+ "\nFirst buildEstimator()"
                + "\n\n");
        println("\nSecond buildEstimator()" 
		+ "\n\n");
      }
    }
    catch (Exception ex) {
      String msg = ex.getMessage().toLowerCase();
      if (msg.indexOf("worse than zeror") >= 0) {
        println("warning: performs worse than ZeroR");
        result[0] = true;
        result[1] = true;
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
   * @param attrTypes attribute types that can be estimated
   * @param classType the class type (NUMERIC, NOMINAL, etc.)
   * @param attributeMissing true if the missing values may be in 
   * the attributes
   * @param classMissing true if the missing values may be in the class
   * @param missingLevel the percentage of missing values
   * @return index 0 is true if the test was passed, index 1 is true if test 
   *         was acceptable
   */
  protected boolean[] canHandleMissing(AttrTypes attrTypes,
				       int classType,
				       boolean attributeMissing,
				       boolean classMissing,
				       int missingLevel) {
    
    if (missingLevel == 100)
      print("100% ");
    print("missing");
    if (attributeMissing) {
      print(" attribute");
      if (classMissing)
        print(" and");
    }
    if (classMissing)
      print(" class");
    print(" values");
    printAttributeSummary(attrTypes, classType);

    print("...");
    FastVector accepts = new FastVector();
    accepts.addElement("missing");
    accepts.addElement("value");
    accepts.addElement("train");
    int numTrain = getNumInstances(), numTest = getNumInstances(), 
    numClasses = 2;
    
    int numAtts = 1, attrIndex = 0;
    return runBasicTest(attrTypes,
			numAtts, attrIndex,
			classType, 
			missingLevel, attributeMissing, classMissing,
			numTrain, numTest, numClasses, 
			accepts);
  }
  
  /**
   * Checks whether an incremental scheme produces the same model when
   * trained incrementally as when batch trained. The model itself
   * cannot be compared, so we compare the evaluation on test data
   * for both models. It is possible to get a false positive on this
   * test (likelihood depends on the estimator).
   *
   * @param attrTypes attribute types that can be estimated
   * @param classType the class type (NUMERIC, NOMINAL, etc.)
   * @return index 0 is true if the test was passed
   */
  protected boolean[] incrementingEquality(AttrTypes attrTypes,
					   int classType) {
    
    print("incremental training produces the same results"
        + " as batch training");
    printAttributeSummary(attrTypes, classType);

    print("...");
    int numTrain = getNumInstances(), numTest = getNumInstances(), 
    numClasses = 2, missingLevel = 0;
    boolean attributeMissing = false, classMissing = false;
    
    boolean[] result = new boolean[2];
    Instances train = null;
    Estimator [] estimators = null;
    boolean built = false;
    int attrIndex = 0;
    Vector test;
    try {
      train = makeTestDataset(42, numTrain, 1, attrTypes,
                              numClasses, 
                              classType
                              );

        // prepare training data set and test value list
      test = makeTestValueList(24, numTest, train, attrIndex,
			       attrTypes.getSetType());

      if (missingLevel > 0) {
        addMissing(train, missingLevel, attributeMissing, classMissing, attrIndex);
      }
      estimators = Estimator.makeCopies(getEstimator(), 2);
      estimators[0].addValues(train, attrIndex);
    } catch (Exception ex) {
      throw new Error("Error setting up for tests: " + ex.getMessage());
    }
    try {
      for (int i = 0; i < train.numInstances(); i++) {
        ((IncrementalEstimator)estimators[1]).addValue(train.instance(i).value(attrIndex), 1.0);
      }
      built = true;
      if (!estimators[0].equals(estimators[1])) {
        println("no");
        result[0] = false;
       
        if (m_Debug) {
          println("\n=== Full Report ===");
          println("Results differ between batch and "
              + "incrementally built models.\n"
              + "Depending on the estimator, this may be OK");
          println("Here are the results:\n");
          println("batch built results\n" + estimators[0].toString());
          println("incrementally built results\n" + estimators[1].toString());
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
   * Checks whether the estimator can handle instance weights.
   * This test compares the estimator performance on two datasets
   * that are identical except for the training weights. If the 
   * results change, then the estimator must be using the weights. It
   * may be possible to get a false positive from this test if the 
   * weight changes aren't significant enough to induce a change
   * in estimator performance (but the weights are chosen to minimize
   * the likelihood of this).
   *
   * @param attrTypes attribute types that can be estimated
   * @param classType the class type (NUMERIC, NOMINAL, etc.)
   * @return index 0 true if the test was passed
   */
  protected boolean[] instanceWeights(AttrTypes attrTypes,
				      int classType) {
    
    print("estimator uses instance weights");
    printAttributeSummary(attrTypes, classType);

    print("...");

    int numTrain = 2 * getNumInstances(), numTest = getNumInstances(), 
      numClasses = 2, missingLevel = 0;
    boolean attributeMissing = false, classMissing = false;
    
    boolean[] result = new boolean[2];
    Instances train = null;
    Vector test = null;
    Estimator [] estimators = null;
    
    Vector resultProbsO = null;
    Vector resultProbsW = null;
    boolean built = false;
    boolean evalFail = false;
    int attrIndex = 0;
    try {
      train = makeTestDataset(42, numTrain, 1, 
                              attrTypes, numClasses, 
                              classType);
  
      // prepare training data set and test value list
      test = makeTestValueList(24, numTest, train, attrIndex,
			       attrTypes.getSetType());

      if (missingLevel > 0) {
        addMissing(train, missingLevel, attributeMissing, classMissing, attrIndex);
      }

      estimators = Estimator.makeCopies(getEstimator(), 2);

      estimators[0].addValues(train, attrIndex);
      resultProbsO = testWithTestValues(estimators[0], test);

    } catch (Exception ex) {
      throw new Error("Error setting up for tests: " + ex.getMessage());
    }
    try {
            
      // Now modify instance weights and re-built
      for (int i = 0; i < train.numInstances(); i++) {
        train.instance(i).setWeight(0);
      }
      Random random = new Random(1);
      for (int i = 0; i < train.numInstances() / 2; i++) {
        int inst = Math.abs(random.nextInt()) % train.numInstances();
        int weight = Math.abs(random.nextInt()) % 10 + 1;
        train.instance(inst).setWeight(weight);
      }
      estimators[1].addValues(train, attrIndex);
      resultProbsW = testWithTestValues(estimators[1], test);

      built = true;
      if (resultProbsO.equals(resultProbsW)) {
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
          println(probsToString(resultProbsO));
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
   * @param attrTypes attribute types that can be estimated
   * @param classType the class type (NUMERIC, NOMINAL, etc.)
   * @param attributeMissing true if we know the estimator can handle
   * (at least) moderate missing attribute values
   * @param classMissing true if we know the estimator can handle
   * (at least) moderate missing class values
   * @return index 0 is true if the test was passed
   */
  protected boolean[] datasetIntegrity(AttrTypes attrTypes,
				       int classType,
				       boolean attributeMissing,
				       boolean classMissing) {
    
    Estimator estimator = null;
    print("estimator doesn't alter original datasets");
    printAttributeSummary(attrTypes, classType);
    print("...");
    int numTrain = getNumInstances(), numTest = getNumInstances(), 
    numClasses = 2, missingLevel = 100;
    
    boolean[] result = new boolean[2];
    Instances train = null;
     boolean built = false;
    try {
      train = makeTestDataset(42, numTrain, 1, attrTypes,
                              numClasses, 
                              classType);
      int attrIndex = 0;
 
      if (missingLevel > 0) {
        addMissing(train, missingLevel, attributeMissing, classMissing, attrIndex);
      }
      estimator = Estimator.makeCopies(getEstimator(), 1)[0];
    } catch (Exception ex) {
      throw new Error("Error setting up for tests: " + ex.getMessage());
    }
    try {
      Instances trainCopy = new Instances(train);
      int attrIndex = 0;
      estimator.addValues(trainCopy, attrIndex);
      compareDatasets(train, trainCopy);
      built = true;
      
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
      }
    }
    
    return result;
  }
  
  /**
   * Runs a text on the datasets with the given characteristics.
   * 
   * @param attrTypes attribute types that can be estimated
   * @param numAtts number of attributes
   * @param attrIndex attribute index 
   * @param classType the class type (NUMERIC, NOMINAL, etc.)
   * @param missingLevel the percentage of missing values
   * @param attributeMissing true if the missing values may be in 
   * the attributes
   * @param classMissing true if the missing values may be in the class
   * @param numTrain the number of instances in the training set
   * @param numTest the number of instaces in the test set
   * @param numClasses the number of classes
   * @param accepts the acceptable string in an exception
   * @return index 0 is true if the test was passed, index 1 is true if test 
   *         was acceptable
   */
  protected boolean[] runBasicTest(AttrTypes attrTypes,
				   int numAtts,
				   int attrIndex,
				   int classType,
				   int missingLevel,
				   boolean attributeMissing,
				   boolean classMissing,
				   int numTrain,
				   int numTest,
				   int numClasses,
				   FastVector accepts) {
    
    return runBasicTest(attrTypes,
			numAtts,
			attrIndex,
			classType, 
			TestInstances.CLASS_IS_LAST,
			missingLevel,
			attributeMissing,
			classMissing,
			numTrain,
			numTest,
			numClasses,
		accepts);
  }
  
  /**
   * Runs a text on the datasets with the given characteristics.
   * 
   * @param attrTypes attribute types that can be estimated
   * @param numAtts number of attributes
   * @param classType the class type (NUMERIC, NOMINAL, etc.)
   * @param classIndex the attribute index of the class
   * @param missingLevel the percentage of missing values
   * @param attributeMissing true if the missing values may be in 
   * the attributes
   * @param classMissing true if the missing values may be in the class
   * @param numTrain the number of instances in the training set
   * @param numTest the number of instaces in the test set
   * @param numClasses the number of classes
   * @param accepts the acceptable string in an exception
   * @return index 0 is true if the test was passed, index 1 is true if test 
   *         was acceptable
   */
  protected boolean[] runBasicTest(AttrTypes attrTypes,
				   int numAtts,
				   int attrIndex,
				   int classType,
				   int classIndex,
				   int missingLevel,
				   boolean attributeMissing,
				   boolean classMissing,
				   int numTrain,
				   int numTest,
				   int numClasses,
				   FastVector accepts) {
    
    boolean[] result = new boolean[2];
    Instances train = null;
    Vector test = null;
    Estimator estimator = null;
    boolean built = false;
   
    try {
      train = makeTestDataset(42, numTrain, numAtts, attrTypes,
          numClasses, 
          classType,
          classIndex);
            
      // prepare training data set and test value list
      if (numTrain > 0) {
        test = makeTestValueList(24, numTest, train, attrIndex,
            attrTypes.getSetType());
     
      } else {
        double min = -10.0;
        double max = 8.0;
        test = makeTestValueList(24, numTest, min, max,
            attrTypes.getSetType());
     }
      
      if (missingLevel > 0) {
        addMissing(train, missingLevel, attributeMissing, classMissing, attrIndex);
      }
      estimator = Estimator.makeCopies(getEstimator(), 1)[0];
    } catch (Exception ex) {
      ex.printStackTrace();
      throw new Error("Error setting up for tests: " + ex.getMessage());
    }
    try {
      estimator.addValues(train, attrIndex);
      built = true;
      
      testWithTestValues(estimator, test);
      
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
   * Compare two datasets to see if they differ.
   *
   * @param data1 one set of instances
   * @param data2 the other set of instances
   * @throws Exception if the datasets differ
   */
  protected void compareDatasets(Instances data1, Instances data2)
  throws Exception {
    if (!data2.equalHeaders(data1)) {
      throw new Exception("header has been modified\n" + data2.equalHeadersMsg(data1));
    }
    if (!(data2.numInstances() == data1.numInstances())) {
      throw new Exception("number of instances has changed");
    }
    for (int i = 0; i < data2.numInstances(); i++) {
      Instance orig = data1.instance(i);
      Instance copy = data2.instance(i);
      for (int j = 0; j < orig.numAttributes(); j++) {
        if (orig.isMissing(j)) {
          if (!copy.isMissing(j)) {
            throw new Exception("instances have changed");
          }
        } else if (orig.value(j) != copy.value(j)) {
          throw new Exception("instances have changed");
        }
        if (orig.weight() != copy.weight()) {
          throw new Exception("instance weights have changed");
        }	  
      }
    }
  }
  
  /**
   * Add missing values to a dataset.
   *
   * @param data the instances to add missing values to
   * @param level the level of missing values to add (if positive, this
   * is the probability that a value will be set to missing, if negative
   * all but one value will be set to missing (not yet implemented))
   * @param attributeMissing if true, attributes will be modified
   * @param classMissing if true, the class attribute will be modified
   * @param attrIndex index of the attribute
   */
  protected void addMissing(Instances data, int level,
			    boolean attributeMissing, boolean classMissing,
			    int attrIndex) {
    
    int classIndex = data.classIndex();
    Random random = new Random(1);
    for (int i = 0; i < data.numInstances(); i++) {
      Instance current = data.instance(i);

      for (int j = 0; j < data.numAttributes(); j++) {
        if (((j == classIndex) && classMissing) ||
            ((j == attrIndex) && attributeMissing)) {
          if (Math.abs(random.nextInt()) % 100 < level)
            current.setMissing(j);
        }
      }
    }
  }
  
  /**
   * Make a simple set of instances, which can later be modified
   * for use in specific tests.
   *
   * @param seed the random number seed
   * @param numInstances the number of instances to generate
   * @param numAttr the number of attributes
   * @param attrTypes the attribute types
   * @param numClasses the number of classes (if nominal class)
   * @param classType the class type (NUMERIC, NOMINAL, etc.)
   * @return the test dataset
   * @throws Exception if the dataset couldn't be generated
   * @see #process(Instances)
   */
  protected Instances makeTestDataset(int seed, 
				      int numInstances, 
				      int numAttr,
				      AttrTypes attrTypes,
				      int numClasses, 
				      int classType)
    throws Exception {
    
    return makeTestDataset(
			   seed,
			   numInstances,
			   numAttr,
			   attrTypes,
			   numClasses, 
			   classType,
			   TestInstances.CLASS_IS_LAST);
  }


  /**
   * Make a simple set of instances with variable position of the class 
   * attribute, which can later be modified for use in specific tests.
   *
   * @param seed the random number seed
   * @param numInstances the number of instances to generate
   * @param numAttr the number of attributes to generate
   * @param attrTypes the type of attrbute that is excepted
   * @param numClasses the number of classes (if nominal class)
   * @param classType the class type (NUMERIC, NOMINAL, etc.)
   * @param classIndex the index of the class (0-based, -1 as last)
   * @return the test dataset
   * @throws Exception if the dataset couldn't be generated
   * @see TestInstances#CLASS_IS_LAST
   * @see #process(Instances)
   */
  protected Instances makeTestDataset(int seed, int numInstances, 
				      int numAttr, AttrTypes attrTypes,
				      int numClasses, int classType,
				      int classIndex)
    throws Exception {
    
    TestInstances dataset = new TestInstances();
    
    dataset.setSeed(seed);
    dataset.setNumInstances(numInstances);
    dataset.setNumNominal   (attrTypes.nominal     ? numAttr : 0);
    dataset.setNumNumeric   (attrTypes.numeric     ? numAttr : 0);
    dataset.setNumString    (attrTypes.string      ? numAttr : 0);
    dataset.setNumDate      (attrTypes.date        ? numAttr : 0);
    dataset.setNumRelational(attrTypes.relational  ? numAttr : 0);
    dataset.setNumClasses(numClasses);
    dataset.setClassType(classType);
    dataset.setClassIndex(classIndex);
    
    return process(dataset.generate());
  }

  /**
   * Make a simple set of values. Only one of the num'type' parameters should be larger 0.
   * (just to make parameter similar to the makeTestDataset parameters)
   *
   * @param seed the random number seed
   * @param numValues the number of values to generate
   * @param data the dataset to make test examples for
   * @param attrIndex index of the attribute
   * @param attrType the class type (NUMERIC, NOMINAL, etc.)
   * @throws Exception if the dataset couldn't be generated
   * @see #process(Instances)
   */
  protected Vector makeTestValueList(int seed, int numValues, 
      Instances data, int attrIndex, int attrType)
  throws Exception {
    
    // get min max
    double []minMax = getMinimumMaximum(data, attrIndex);
    double minValue = minMax[0];
    double maxValue = minMax[1];
    
    // make value list and put into a VECTOR
    double range = maxValue - minValue; 
    Vector values = new Vector(numValues); 
    Random random = new Random(seed);
    
    if (attrType == Attribute.NOMINAL || attrType == PreferenceAttribute.RANKING) {
      for (int i = 0; i < numValues; i++) {
        Double v = new Double((Math.abs(random.nextInt()) % (int)range)+ (int)minValue);
        values.add(v);
      }
    }
    if (attrType == Attribute.NUMERIC) {
      for (int i = 0; i < numValues; i++) {
        Double v = new Double(random.nextDouble() * range + minValue);
        values.add(v);
      }
    }
    return values;
  }

  /**
   * Make a simple set of values. Only one of the num'type' parameters should be larger 0.
   * (just to make parameter similar to the makeTestDataset parameters)
   *
   * @param seed the random number seed
   * @param numValues the number of values to generate
   * @param minValue the minimal data value
   * @param maxValue the maximal data value
   * @param attrType the class type (NUMERIC, NOMINAL, etc.)
   * @throws Exception if the dataset couldn't be generated
   * @see #process(Instances)
   */
  protected Vector makeTestValueList(int seed, int numValues, 
      double minValue, double maxValue, int attrType)
  throws Exception {
    
      
    // make value list and put into a VECTOR
    double range = maxValue - minValue; 
    Vector values = new Vector(numValues); 
    Random random = new Random(seed);
    
    if (attrType == Attribute.NOMINAL || attrType == PreferenceAttribute.RANKING) {
      for (int i = 0; i < numValues; i++) {
        Double v = new Double((Math.abs(random.nextInt()) % (int)range)+ (int)minValue);
        values.add(v);
      }
    }
    if (attrType == Attribute.NUMERIC) {
      for (int i = 0; i < numValues; i++) {
        Double v = new Double(random.nextDouble() * range + minValue);
        values.add(v);
      }
    }
    return values;
  }

  /**
   * Test with test values.
   *
   * @param est estimator to be tested
   * @param test vector with test values
   *
   **/
  protected Vector testWithTestValues(Estimator est, Vector test) {
    
    Vector results = new Vector();
    for (int i = 0; i < test.size(); i++) {
      double testValue = ((Double)(test.elementAt(i))).doubleValue();
      double prob = est.getProbability(testValue);
      Double p = new Double(prob);
      results.add(p);
    }
    return results;
  }

  /**
   * Gets the minimum and maximum of the values a the first attribute
   * of the given data set
   *
   * @param inst the instance
   * @param attrIndex the index of the attribut to find min and max
   * @return the array with the minimum value on index 0 and the max on index 1
   */
  
  protected double[] getMinimumMaximum(Instances inst, int attrIndex) {
    double []minMax = new double[2];
    
    try {
      int num = getMinMax(inst, attrIndex, minMax);
    } catch (Exception ex) {
      ex.printStackTrace();
      System.out.println(ex.getMessage());
    }
    return minMax;
    //      double minValue = minMax[0];
    //      double maxValue = minMax[1];
  }
  
  /** 
   * Find the minimum and the maximum of the attribute and return it in 
   * the last parameter..
   * @param inst instances used to build the estimator
   * @param attrIndex index of the attribute
   * @param minMax the array to return minimum and maximum in
   * @return number of not missing values
   * @exception Exception if parameter minMax wasn't initialized properly
   */
  public static int getMinMax(Instances inst, int attrIndex, double [] minMax) 
    throws Exception {
    double min = Double.NaN;
    double max = Double.NaN;
    Instance instance = null;
    int numNotMissing = 0;
    if ((minMax == null) || (minMax.length < 2)) {
      throw new Exception("Error in Program, privat method getMinMax");
    }
    
    Enumeration enumInst = inst.enumerateInstances();
    if (enumInst.hasMoreElements()) {
      do {
	instance = (Instance) enumInst.nextElement();
      } while (instance.isMissing(attrIndex) && (enumInst.hasMoreElements()));
      
      // add values if not  missing
      if (!instance.isMissing(attrIndex)) {
	numNotMissing++;
	min = instance.value(attrIndex);
	max = instance.value(attrIndex);
      }
      while (enumInst.hasMoreElements()) {
	instance = (Instance) enumInst.nextElement();
	if (!instance.isMissing(attrIndex)) {
	  numNotMissing++;
	  if (instance.value(attrIndex) < min) {
	    min = (instance.value(attrIndex));
	  } else {
	    if (instance.value(attrIndex) > max) {	      
	      max = (instance.value(attrIndex));
	    }
	  }
	}
      }
    }
    minMax[0] = min;
    minMax[1] = max;
    return numNotMissing;
  }

  /**
   * Print the probabilities after testing
   * @param probs vector with probability values
   * @return string with probability values printed
   */ 
  private String probsToString(Vector probs) {
    StringBuffer txt = new StringBuffer (" ");
    for (int i = 0; i < probs.size(); i++) {
      txt.append("" + ((Double)(probs.elementAt(i))).doubleValue() + " ");
    }
    return txt.toString();
  }
  
  /**
   * Provides a hook for derived classes to further modify the data. 
   * 
   * @param data	the data to process
   * @return		the processed data
   * @see #m_PostProcessor
   */
  protected Instances process(Instances data) {
    if (getPostProcessor() == null)
      return data;
    else
      return getPostProcessor().process(data);
  }
  
  /**
   * Print out a short summary string for the dataset characteristics
   *
   * @param attrTypes the attribute types used (NUMERIC, NOMINAL, etc.)
   * @param classType the class type (NUMERIC, NOMINAL, etc.)
   */
  protected void printAttributeSummary(AttrTypes attrTypes, int classType) {
    
    String str = "";
    
    if (attrTypes.numeric)
      str += " numeric";
    
    if (attrTypes.nominal) {
      if (str.length() > 0)
        str += " &";
      str += " nominal";
    }
    
    if (attrTypes.string) {
      if (str.length() > 0)
        str += " &";
      str += " string";
    }
    
    if (attrTypes.date) {
      if (str.length() > 0)
        str += " &";
      str += " date";
    }
    
    if (attrTypes.relational) {
      if (str.length() > 0)
        str += " &";
      str += " relational";
    }
    
    str += " attributes)";
    
    switch (classType) {
      case Attribute.NUMERIC:
        str = " (numeric class," + str;
        break;
      case Attribute.NOMINAL:
        str = " (nominal class," + str;
        break;
      case PreferenceAttribute.RANKING:
          str = " (ranking," + str;
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
   * Print out a short summary string for the dataset characteristics
   *
   * @param attrType the attribute type (NUMERIC, NOMINAL, etc.)
   * @param classType the class type (NUMERIC, NOMINAL, etc.)
   */
  protected void printAttributeSummary(int attrType, int classType) {
    
    String str = "";
    
    switch (attrType) {
    case Attribute.NUMERIC:
      str = " numeric" + str;
      break;
    case Attribute.NOMINAL:
      str = " nominal" + str;
      break;
    case PreferenceAttribute.RANKING:
      str = " ranking," + str;
      break;
    case Attribute.STRING:
      str = " string" + str;
      break;
    case Attribute.DATE:
      str = " date" + str;
      break;
    case Attribute.RELATIONAL:
      str = " relational" + str;
      break;
    }
    str += " attribute(s))";
    
    switch (classType) {
    case Attribute.NUMERIC:
      str = " (numeric class," + str;
      break;
    case Attribute.NOMINAL:
      str = " (nominal class," + str;
      break;
    case PreferenceAttribute.RANKING:
        str = " (rannking," + str;
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
    return RevisionUtils.extract("$Revision: 4997 $");
  }

  /**
   * Test method for this class
   * 
   * @param args the commandline parameters
   */
  public static void main(String [] args) {
    try {
      CheckEstimator check = new CheckEstimator();
      
      try {
        check.setOptions(args);
        Utils.checkForRemainingOptions(args);
      } catch (Exception ex) {
        String result = ex.getMessage() + "\n\n" + check.getClass().getName().replaceAll(".*\\.", "") + " Options:\n\n";
        Enumeration enu = check.listOptions();
        while (enu.hasMoreElements()) {
          Option option = (Option) enu.nextElement();
          result += option.synopsis() + "\n" + option.description() + "\n";
        }
        throw new Exception(result);
      }
      
      check.doTests();
    } catch (Exception ex) {
      System.err.println(ex.getMessage());
    }
  }
}

