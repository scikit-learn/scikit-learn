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
 *    AveragingResultProducer.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.experiment;

import weka.core.AdditionalMeasureProducer;
import weka.core.FastVector;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.core.Utils;

import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * Takes the results from a ResultProducer and submits the average to the result listener. Normally used with a CrossValidationResultProducer to perform n x m fold cross validation. For non-numeric result fields, the first value is used.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -F &lt;field name&gt;
 *  The name of the field to average over.
 *  (default "Fold")</pre>
 * 
 * <pre> -X &lt;num results&gt;
 *  The number of results expected per average.
 *  (default 10)</pre>
 * 
 * <pre> -S
 *  Calculate standard deviations.
 *  (default only averages)</pre>
 * 
 * <pre> -W &lt;class name&gt;
 *  The full class name of a ResultProducer.
 *  eg: weka.experiment.CrossValidationResultProducer</pre>
 * 
 * <pre> 
 * Options specific to result producer weka.experiment.CrossValidationResultProducer:
 * </pre>
 * 
 * <pre> -X &lt;number of folds&gt;
 *  The number of folds to use for the cross-validation.
 *  (default 10)</pre>
 * 
 * <pre> -D
 * Save raw split evaluator output.</pre>
 * 
 * <pre> -O &lt;file/directory name/path&gt;
 *  The filename where raw output will be stored.
 *  If a directory name is specified then then individual
 *  outputs will be gzipped, otherwise all output will be
 *  zipped to the named file. Use in conjuction with -D. (default splitEvalutorOut.zip)</pre>
 * 
 * <pre> -W &lt;class name&gt;
 *  The full class name of a SplitEvaluator.
 *  eg: weka.experiment.ClassifierSplitEvaluator</pre>
 * 
 * <pre> 
 * Options specific to split evaluator weka.experiment.ClassifierSplitEvaluator:
 * </pre>
 * 
 * <pre> -W &lt;class name&gt;
 *  The full class name of the classifier.
 *  eg: weka.classifiers.bayes.NaiveBayes</pre>
 * 
 * <pre> -C &lt;index&gt;
 *  The index of the class for which IR statistics
 *  are to be output. (default 1)</pre>
 * 
 * <pre> -I &lt;index&gt;
 *  The index of an attribute to output in the
 *  results. This attribute should identify an
 *  instance in order to know which instances are
 *  in the test set of a cross validation. if 0
 *  no output (default 0).</pre>
 * 
 * <pre> -P
 *  Add target and prediction columns to the result
 *  for each fold.</pre>
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
 * All options after -- will be passed to the result producer.
 * 
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 6420 $
 */
public class AveragingResultProducer 
  implements ResultListener, ResultProducer, OptionHandler,
	     AdditionalMeasureProducer, RevisionHandler {

  /** for serialization */
  static final long serialVersionUID = 2551284958501991352L;
  
  /** The dataset of interest */
  protected Instances m_Instances;

  /** The ResultListener to send results to */
  protected ResultListener m_ResultListener = new CSVResultListener();

  /** The ResultProducer used to generate results */
  protected ResultProducer m_ResultProducer
    = new CrossValidationResultProducer();

  /** The names of any additional measures to look for in SplitEvaluators */
  protected String [] m_AdditionalMeasures = null;
  
  /** The number of results expected to average over for each run */
  protected int m_ExpectedResultsPerAverage = 10;

  /** True if standard deviation fields should be produced */
  protected boolean m_CalculateStdDevs;
    
  /**
   * The name of the field that will contain the number of results
   * averaged over.
   */
  protected String m_CountFieldName = "Num_" + CrossValidationResultProducer
    .FOLD_FIELD_NAME;

  /** The name of the key field to average over */
  protected String m_KeyFieldName = CrossValidationResultProducer
    .FOLD_FIELD_NAME;

  /** The index of the field to average over in the resultproducers key */
  protected int m_KeyIndex = -1;

  /** Collects the keys from a single run */
  protected FastVector m_Keys = new FastVector();
  
  /** Collects the results from a single run */
  protected FastVector m_Results = new FastVector();

  /**
   * Returns a string describing this result producer
   * @return a description of the result producer suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return "Takes the results from a ResultProducer "
      +"and submits the average to the result listener. Normally used with "
      +"a CrossValidationResultProducer to perform n x m fold cross "
      +"validation. For non-numeric result fields, the first value is used.";
  }

  /**
   * Scans through the key field names of the result producer to find
   * the index of the key field to average over. Sets the value of
   * m_KeyIndex to the index, or -1 if no matching key field was found.
   *
   * @return the index of the key field to average over
   */
  protected int findKeyIndex() {

    m_KeyIndex = -1;
    try {
      if (m_ResultProducer != null) {
	String [] keyNames = m_ResultProducer.getKeyNames();
	for (int i = 0; i < keyNames.length; i++) {
	  if (keyNames[i].equals(m_KeyFieldName)) {
	    m_KeyIndex = i;
	    break;
	  }
	}
      }
    } catch (Exception ex) {
    }
    return m_KeyIndex;
  }

  /**
   * Determines if there are any constraints (imposed by the
   * destination) on the result columns to be produced by
   * resultProducers. Null should be returned if there are NO
   * constraints, otherwise a list of column names should be
   * returned as an array of Strings.
   * @param rp the ResultProducer to which the constraints will apply
   * @return an array of column names to which resutltProducer's
   * results will be restricted.
   * @throws Exception if constraints can't be determined
   */
  public String [] determineColumnConstraints(ResultProducer rp) 
    throws Exception {
    return null;
  }

  /**
   * Simulates a run to collect the keys the sub-resultproducer could
   * generate. Does some checking on the keys and determines the 
   * template key.
   *
   * @param run the run number
   * @return a template key (null for the field being averaged)
   * @throws Exception if an error occurs
   */
  protected Object [] determineTemplate(int run) throws Exception {

    if (m_Instances == null) {
      throw new Exception("No Instances set");
    }
    m_ResultProducer.setInstances(m_Instances);

    // Clear the collected results
    m_Keys.removeAllElements();
    m_Results.removeAllElements();
    
    m_ResultProducer.doRunKeys(run);
    checkForMultipleDifferences();

    Object [] template = (Object [])((Object [])m_Keys.elementAt(0)).clone();
    template[m_KeyIndex] = null;
    // Check for duplicate keys
    checkForDuplicateKeys(template);

    return template;
  }

  /**
   * Gets the keys for a specified run number. Different run
   * numbers correspond to different randomizations of the data. Keys
   * produced should be sent to the current ResultListener
   *
   * @param run the run number to get keys for.
   * @throws Exception if a problem occurs while getting the keys
   */
  public void doRunKeys(int run) throws Exception {

    // Generate the template
    Object [] template = determineTemplate(run);
    String [] newKey = new String [template.length - 1];
    System.arraycopy(template, 0, newKey, 0, m_KeyIndex);
    System.arraycopy(template, m_KeyIndex + 1,
		     newKey, m_KeyIndex,
		     template.length - m_KeyIndex - 1);
    m_ResultListener.acceptResult(this, newKey, null);      
  }

  /**
   * Gets the results for a specified run number. Different run
   * numbers correspond to different randomizations of the data. Results
   * produced should be sent to the current ResultListener
   *
   * @param run the run number to get results for.
   * @throws Exception if a problem occurs while getting the results
   */
  public void doRun(int run) throws Exception {

    // Generate the key and ask whether the result is required
    Object [] template = determineTemplate(run);
    String [] newKey = new String [template.length - 1];
    System.arraycopy(template, 0, newKey, 0, m_KeyIndex);
    System.arraycopy(template, m_KeyIndex + 1,
		     newKey, m_KeyIndex,
		     template.length - m_KeyIndex - 1);

    if (m_ResultListener.isResultRequired(this, newKey)) {
      // Clear the collected keys
      m_Keys.removeAllElements();
      m_Results.removeAllElements();
      
      m_ResultProducer.doRun(run);
      
      // Average the results collected
      //System.err.println("Number of results collected: " + m_Keys.size());
      
      // Check that the keys only differ on the selected key field
      checkForMultipleDifferences();
      
      template = (Object [])((Object [])m_Keys.elementAt(0)).clone();
      template[m_KeyIndex] = null;
      // Check for duplicate keys
      checkForDuplicateKeys(template);
      // Calculate the average and submit it if necessary
      doAverageResult(template);
    }
  }

  
  /**
   * Compares a key to a template to see whether they match. Null
   * fields in the template are ignored in the matching.
   *
   * @param template the template to match against
   * @param test the key to test
   * @return true if the test key matches the template on all non-null template
   * fields
   */
  protected boolean matchesTemplate(Object [] template, Object [] test) {
    
    if (template.length != test.length) {
      return false;
    }
    for (int i = 0; i < test.length; i++) {
      if ((template[i] != null) && (!template[i].equals(test[i]))) {
	return false;
      }
    }
    return true;
  }
  
  /**
   * Asks the resultlistener whether an average result is required, and
   * if so, calculates it.
   *
   * @param template the template to match keys against when calculating the
   * average
   * @throws Exception if an error occurs
   */
  protected void doAverageResult(Object [] template) throws Exception {

    // Generate the key and ask whether the result is required
    String [] newKey = new String [template.length - 1];
    System.arraycopy(template, 0, newKey, 0, m_KeyIndex);
    System.arraycopy(template, m_KeyIndex + 1,
		     newKey, m_KeyIndex,
		     template.length - m_KeyIndex - 1);
    if (m_ResultListener.isResultRequired(this, newKey)) {
      Object [] resultTypes = m_ResultProducer.getResultTypes();
      Stats [] stats = new Stats [resultTypes.length];
      for (int i = 0; i < stats.length; i++) {
	stats[i] = new Stats();
      }
      Object [] result = getResultTypes();
      int numMatches = 0;
      for (int i = 0; i < m_Keys.size(); i++) {
	Object [] currentKey = (Object [])m_Keys.elementAt(i);
	// Skip non-matching keys
	if (!matchesTemplate(template, currentKey)) {
	  continue;
	}
	// Add the results to the stats accumulator
	Object [] currentResult = (Object [])m_Results.elementAt(i);
	numMatches++;
	for (int j = 0; j < resultTypes.length; j++) {
	  if (resultTypes[j] instanceof Double) {
	    if (currentResult[j] == null) {

	      // set the stats object for this result to null---
	      // more than likely this is an additional measure field
	      // not supported by the low level split evaluator
	      if (stats[j] != null) {
		stats[j] = null;
	      }
	      
	      /* throw new Exception("Null numeric result field found:\n"
		 + DatabaseUtils.arrayToString(currentKey)
		 + " -- "
		 + DatabaseUtils
		 .arrayToString(currentResult)); */
	    }
	    if (stats[j] != null) {
	      double currentVal = ((Double)currentResult[j]).doubleValue();
	      stats[j].add(currentVal);
	    }
	  }
	}
      }
      if (numMatches != m_ExpectedResultsPerAverage) {
	throw new Exception("Expected " + m_ExpectedResultsPerAverage
			    + " results matching key \""
			    + DatabaseUtils.arrayToString(template)
			    + "\" but got "
			    + numMatches);
      }
      result[0] = new Double(numMatches);
      Object [] currentResult = (Object [])m_Results.elementAt(0);
      int k = 1;
      for (int j = 0; j < resultTypes.length; j++) {
	if (resultTypes[j] instanceof Double) {
	  if (stats[j] != null) {
	    stats[j].calculateDerived();
	    result[k++] = new Double(stats[j].mean);
	  } else {
	    result[k++] = null;
	  }
	  if (getCalculateStdDevs()) {
	    if (stats[j] != null) {
	      result[k++] = new Double(stats[j].stdDev);
	    } else {
	      result[k++] = null;
	    }
	  }
	} else {
	  result[k++] = currentResult[j];
	}
      }
      m_ResultListener.acceptResult(this, newKey, result);      
    }
  }
  
  /**
   * Checks whether any duplicate results (with respect to a key template)
   * were received.
   *
   * @param template the template key.
   * @throws Exception if duplicate results are detected
   */
  protected void checkForDuplicateKeys(Object [] template) throws Exception {

    Hashtable hash = new Hashtable();
    int numMatches = 0;
    for (int i = 0; i < m_Keys.size(); i++) {
      Object [] current = (Object [])m_Keys.elementAt(i);
      // Skip non-matching keys
      if (!matchesTemplate(template, current)) {
	continue;
      }
      if (hash.containsKey(current[m_KeyIndex])) {
	throw new Exception("Duplicate result received:"
			    + DatabaseUtils.arrayToString(current));
      }
      numMatches++;
      hash.put(current[m_KeyIndex], current[m_KeyIndex]);
    }
    if (numMatches != m_ExpectedResultsPerAverage) {
      throw new Exception("Expected " + m_ExpectedResultsPerAverage
			  + " results matching key \""
			  + DatabaseUtils.arrayToString(template)
			  + "\" but got "
			  + numMatches);
    }
  }
  
  /**
   * Checks that the keys for a run only differ in one key field. If they
   * differ in more than one field, a more sophisticated averager will submit
   * multiple results - for now an exception is thrown. Currently assumes that
   * the most differences will be shown between the first and last
   * result received.
   *
   * @throws Exception if the keys differ on fields other than the
   * key averaging field
   */
  protected void checkForMultipleDifferences() throws Exception {
    
    Object [] firstKey = (Object [])m_Keys.elementAt(0);
    Object [] lastKey = (Object [])m_Keys.elementAt(m_Keys.size() - 1);
    /*
    System.err.println("First key:" +  DatabaseUtils.arrayToString(firstKey));
    System.err.println("Last key :" + DatabaseUtils.arrayToString(lastKey));
    */
    for (int i = 0; i < firstKey.length; i++) {
      if ((i != m_KeyIndex) && !firstKey[i].equals(lastKey[i])) {
	throw new Exception("Keys differ on fields other than \""
			    + m_KeyFieldName
			    + "\" -- time to implement multiple averaging");
      }
    }
  }
  
  /**
   * Prepare for the results to be received.
   *
   * @param rp the ResultProducer that will generate the results
   * @throws Exception if an error occurs during preprocessing.
   */
  public void preProcess(ResultProducer rp) throws Exception {

    if (m_ResultListener == null) {
      throw new Exception("No ResultListener set");
    }
    m_ResultListener.preProcess(this);
  }

  /**
   * Prepare to generate results. The ResultProducer should call
   * preProcess(this) on the ResultListener it is to send results to.
   *
   * @throws Exception if an error occurs during preprocessing.
   */
  public void preProcess() throws Exception {
    
    if (m_ResultProducer == null) {
      throw new Exception("No ResultProducer set");
    }
    // Tell the resultproducer to send results to us
    m_ResultProducer.setResultListener(this);
    findKeyIndex();
    if (m_KeyIndex == -1) {
      throw new Exception("No key field called " + m_KeyFieldName
			  + " produced by "
			  + m_ResultProducer.getClass().getName());
    }
    m_ResultProducer.preProcess();
  }
  
  /**
   * When this method is called, it indicates that no more results
   * will be sent that need to be grouped together in any way.
   *
   * @param rp the ResultProducer that generated the results
   * @throws Exception if an error occurs
   */
  public void postProcess(ResultProducer rp) throws Exception {

    m_ResultListener.postProcess(this);
  }

  /**
   * When this method is called, it indicates that no more requests to
   * generate results for the current experiment will be sent. The
   * ResultProducer should call preProcess(this) on the
   * ResultListener it is to send results to.
   *
   * @throws Exception if an error occurs
   */
  public void postProcess() throws Exception {

    m_ResultProducer.postProcess();
  }
  
  /**
   * Accepts results from a ResultProducer.
   *
   * @param rp the ResultProducer that generated the results
   * @param key an array of Objects (Strings or Doubles) that uniquely
   * identify a result for a given ResultProducer with given compatibilityState
   * @param result the results stored in an array. The objects stored in
   * the array may be Strings, Doubles, or null (for the missing value).
   * @throws Exception if the result could not be accepted.
   */
  public void acceptResult(ResultProducer rp, Object [] key, Object [] result)
    throws Exception {

    if (m_ResultProducer != rp) {
      throw new Error("Unrecognized ResultProducer sending results!!");
    }
    m_Keys.addElement(key);
    m_Results.addElement(result);
  }

  /**
   * Determines whether the results for a specified key must be
   * generated.
   *
   * @param rp the ResultProducer wanting to generate the results
   * @param key an array of Objects (Strings or Doubles) that uniquely
   * identify a result for a given ResultProducer with given compatibilityState
   * @return true if the result should be generated
   * @throws Exception if it could not be determined if the result 
   * is needed.
   */
  public boolean isResultRequired(ResultProducer rp, Object [] key) 
    throws Exception {

    if (m_ResultProducer != rp) {
      throw new Error("Unrecognized ResultProducer sending results!!");
    }
    return true;
  }

  /**
   * Gets the names of each of the columns produced for a single run.
   *
   * @return an array containing the name of each column
   * @throws Exception if key names cannot be generated
   */
  public String [] getKeyNames() throws Exception {

    if (m_KeyIndex == -1) {
      throw new Exception("No key field called " + m_KeyFieldName
			  + " produced by "
			  + m_ResultProducer.getClass().getName());
    }
    String [] keyNames = m_ResultProducer.getKeyNames();
    String [] newKeyNames = new String [keyNames.length - 1];
    System.arraycopy(keyNames, 0, newKeyNames, 0, m_KeyIndex);
    System.arraycopy(keyNames, m_KeyIndex + 1,
		     newKeyNames, m_KeyIndex,
		     keyNames.length - m_KeyIndex - 1);
    return newKeyNames;
  }

  /**
   * Gets the data types of each of the columns produced for a single run.
   * This method should really be static.
   *
   * @return an array containing objects of the type of each column. The 
   * objects should be Strings, or Doubles.
   * @throws Exception if the key types could not be determined (perhaps
   * because of a problem from a nested sub-resultproducer)
   */
  public Object [] getKeyTypes() throws Exception {

    if (m_KeyIndex == -1) {
      throw new Exception("No key field called " + m_KeyFieldName
			  + " produced by "
			  + m_ResultProducer.getClass().getName());
    }
    Object [] keyTypes = m_ResultProducer.getKeyTypes();
    // Find and remove the key field that is being averaged over
    Object [] newKeyTypes = new String [keyTypes.length - 1];
    System.arraycopy(keyTypes, 0, newKeyTypes, 0, m_KeyIndex);
    System.arraycopy(keyTypes, m_KeyIndex + 1,
		     newKeyTypes, m_KeyIndex,
		     keyTypes.length - m_KeyIndex - 1);
    return newKeyTypes;
  }

  /**
   * Gets the names of each of the columns produced for a single run.
   * A new result field is added for the number of results used to
   * produce each average.
   * If only averages are being produced the names are not altered, if
   * standard deviations are produced then "Dev_" and "Avg_" are prepended
   * to each result deviation and average field respectively.
   *
   * @return an array containing the name of each column
   * @throws Exception if the result names could not be determined (perhaps
   * because of a problem from a nested sub-resultproducer)
   */
  public String [] getResultNames() throws Exception {

    String [] resultNames = m_ResultProducer.getResultNames();
    // Add in the names of our extra Result fields
    if (getCalculateStdDevs()) {
      Object [] resultTypes = m_ResultProducer.getResultTypes();
      int numNumeric = 0;
      for (int i = 0; i < resultTypes.length; i++) {
	if (resultTypes[i] instanceof Double) {
	  numNumeric++;
	}
      }
      String [] newResultNames = new String [resultNames.length +
					    1 + numNumeric];
      newResultNames[0] = m_CountFieldName;
      int j = 1;
      for (int i = 0; i < resultNames.length; i++) {
	newResultNames[j++] = "Avg_" + resultNames[i];
	if (resultTypes[i] instanceof Double) {
	  newResultNames[j++] = "Dev_" + resultNames[i];
	}
      }
      return newResultNames;
    } else {
      String [] newResultNames = new String [resultNames.length + 1];
      newResultNames[0] = m_CountFieldName;
      System.arraycopy(resultNames, 0, newResultNames, 1, resultNames.length);
      return newResultNames;
    }
  }

  /**
   * Gets the data types of each of the columns produced for a single run.
   *
   * @return an array containing objects of the type of each column. The 
   * objects should be Strings, or Doubles.
   * @throws Exception if the result types could not be determined (perhaps
   * because of a problem from a nested sub-resultproducer)
   */
  public Object [] getResultTypes() throws Exception {

    Object [] resultTypes = m_ResultProducer.getResultTypes();
    // Add in the types of our extra Result fields
    if (getCalculateStdDevs()) {
      int numNumeric = 0;
      for (int i = 0; i < resultTypes.length; i++) {
	if (resultTypes[i] instanceof Double) {
	  numNumeric++;
	}
      }
      Object [] newResultTypes = new Object [resultTypes.length +
					    1 + numNumeric];
      newResultTypes[0] = new Double(0);
      int j = 1;
      for (int i = 0; i < resultTypes.length; i++) {
	newResultTypes[j++] = resultTypes[i];
	if (resultTypes[i] instanceof Double) {
	  newResultTypes[j++] = new Double(0);
	}
      }
      return newResultTypes;
    } else {
      Object [] newResultTypes = new Object [resultTypes.length + 1];
      newResultTypes[0] = new Double(0);
      System.arraycopy(resultTypes, 0, newResultTypes, 1, resultTypes.length);
      return newResultTypes;
    }
  }

  /**
   * Gets a description of the internal settings of the result
   * producer, sufficient for distinguishing a ResultProducer
   * instance from another with different settings (ignoring
   * those settings set through this interface). For example,
   * a cross-validation ResultProducer may have a setting for the
   * number of folds. For a given state, the results produced should
   * be compatible. Typically if a ResultProducer is an OptionHandler,
   * this string will represent the command line arguments required
   * to set the ResultProducer to that state.
   *
   * @return the description of the ResultProducer state, or null
   * if no state is defined
   */
  public String getCompatibilityState() {

    String result = // "-F " + Utils.quote(getKeyFieldName())
      " -X " + getExpectedResultsPerAverage() + " ";
    if (getCalculateStdDevs()) {
      result += "-S ";
    }
    if (m_ResultProducer == null) {
      result += "<null ResultProducer>";
    } else {
      result += "-W " + m_ResultProducer.getClass().getName();
      result  += " -- " + m_ResultProducer.getCompatibilityState();
    }
    
    return result.trim();
  }


  /**
   * Returns an enumeration describing the available options..
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {

    Vector newVector = new Vector(2);

    newVector.addElement(new Option(
	     "\tThe name of the field to average over.\n"
	      +"\t(default \"Fold\")", 
	     "F", 1, 
	     "-F <field name>"));
    newVector.addElement(new Option(
	     "\tThe number of results expected per average.\n"
	      +"\t(default 10)", 
	     "X", 1, 
	     "-X <num results>"));
    newVector.addElement(new Option(
	     "\tCalculate standard deviations.\n"
	      +"\t(default only averages)", 
	     "S", 0, 
	     "-S"));
    newVector.addElement(new Option(
	     "\tThe full class name of a ResultProducer.\n"
	      +"\teg: weka.experiment.CrossValidationResultProducer", 
	     "W", 1, 
	     "-W <class name>"));

    if ((m_ResultProducer != null) &&
	(m_ResultProducer instanceof OptionHandler)) {
      newVector.addElement(new Option(
	     "",
	     "", 0, "\nOptions specific to result producer "
	     + m_ResultProducer.getClass().getName() + ":"));
      Enumeration enu = ((OptionHandler)m_ResultProducer).listOptions();
      while (enu.hasMoreElements()) {
	newVector.addElement(enu.nextElement());
      }
    }
    return newVector.elements();
  }

  /**
   * Parses a given list of options. <p/>
   *
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -F &lt;field name&gt;
   *  The name of the field to average over.
   *  (default "Fold")</pre>
   * 
   * <pre> -X &lt;num results&gt;
   *  The number of results expected per average.
   *  (default 10)</pre>
   * 
   * <pre> -S
   *  Calculate standard deviations.
   *  (default only averages)</pre>
   * 
   * <pre> -W &lt;class name&gt;
   *  The full class name of a ResultProducer.
   *  eg: weka.experiment.CrossValidationResultProducer</pre>
   * 
   * <pre> 
   * Options specific to result producer weka.experiment.CrossValidationResultProducer:
   * </pre>
   * 
   * <pre> -X &lt;number of folds&gt;
   *  The number of folds to use for the cross-validation.
   *  (default 10)</pre>
   * 
   * <pre> -D
   * Save raw split evaluator output.</pre>
   * 
   * <pre> -O &lt;file/directory name/path&gt;
   *  The filename where raw output will be stored.
   *  If a directory name is specified then then individual
   *  outputs will be gzipped, otherwise all output will be
   *  zipped to the named file. Use in conjuction with -D. (default splitEvalutorOut.zip)</pre>
   * 
   * <pre> -W &lt;class name&gt;
   *  The full class name of a SplitEvaluator.
   *  eg: weka.experiment.ClassifierSplitEvaluator</pre>
   * 
   * <pre> 
   * Options specific to split evaluator weka.experiment.ClassifierSplitEvaluator:
   * </pre>
   * 
   * <pre> -W &lt;class name&gt;
   *  The full class name of the classifier.
   *  eg: weka.classifiers.bayes.NaiveBayes</pre>
   * 
   * <pre> -C &lt;index&gt;
   *  The index of the class for which IR statistics
   *  are to be output. (default 1)</pre>
   * 
   * <pre> -I &lt;index&gt;
   *  The index of an attribute to output in the
   *  results. This attribute should identify an
   *  instance in order to know which instances are
   *  in the test set of a cross validation. if 0
   *  no output (default 0).</pre>
   * 
   * <pre> -P
   *  Add target and prediction columns to the result
   *  for each fold.</pre>
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
   * All options after -- will be passed to the result producer.
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    
    String keyFieldName = Utils.getOption('F', options);
    if (keyFieldName.length() != 0) {
      setKeyFieldName(keyFieldName);
    } else {
      setKeyFieldName(CrossValidationResultProducer.FOLD_FIELD_NAME);
    }

    String numResults = Utils.getOption('X', options);
    if (numResults.length() != 0) {
      setExpectedResultsPerAverage(Integer.parseInt(numResults));
    } else {
      setExpectedResultsPerAverage(10);
    }

    setCalculateStdDevs(Utils.getFlag('S', options));
    
    String rpName = Utils.getOption('W', options);
    if (rpName.length() == 0) {
      throw new Exception("A ResultProducer must be specified with"
			  + " the -W option.");
    }
    // Do it first without options, so if an exception is thrown during
    // the option setting, listOptions will contain options for the actual
    // RP.
    setResultProducer((ResultProducer)Utils.forName(
		      ResultProducer.class,
		      rpName,
		      null));
    if (getResultProducer() instanceof OptionHandler) {
      ((OptionHandler) getResultProducer())
	.setOptions(Utils.partitionOptions(options));
    }
  }

  /**
   * Gets the current settings of the result producer.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String [] getOptions() {

    String [] seOptions = new String [0];
    if ((m_ResultProducer != null) && 
	(m_ResultProducer instanceof OptionHandler)) {
      seOptions = ((OptionHandler)m_ResultProducer).getOptions();
    }
    
    String [] options = new String [seOptions.length + 8];
    int current = 0;

    options[current++] = "-F";
    options[current++] = "" + getKeyFieldName();
    options[current++] = "-X";
    options[current++] = "" + getExpectedResultsPerAverage();
    if (getCalculateStdDevs()) {
      options[current++] = "-S";
    }
    if (getResultProducer() != null) {
      options[current++] = "-W";
      options[current++] = getResultProducer().getClass().getName();
    }
    options[current++] = "--";

    System.arraycopy(seOptions, 0, options, current, 
		     seOptions.length);
    current += seOptions.length;
    while (current < options.length) {
      options[current++] = "";
    }
    return options;
  }

  /**
   * Set a list of method names for additional measures to look for
   * in SplitEvaluators. This could contain many measures (of which only a
   * subset may be produceable by the current resultProducer) if an experiment
   * is the type that iterates over a set of properties.
   * @param additionalMeasures an array of measure names, null if none
   */
  public void setAdditionalMeasures(String [] additionalMeasures) {
    m_AdditionalMeasures = additionalMeasures;

    if (m_ResultProducer != null) {
      System.err.println("AveragingResultProducer: setting additional "
			 +"measures for "
			 +"ResultProducer");
      m_ResultProducer.setAdditionalMeasures(m_AdditionalMeasures);
    }
  }

  /**
   * Returns an enumeration of any additional measure names that might be
   * in the result producer
   * @return an enumeration of the measure names
   */
  public Enumeration enumerateMeasures() {
    Vector newVector = new Vector();
    if (m_ResultProducer instanceof AdditionalMeasureProducer) {
      Enumeration en = ((AdditionalMeasureProducer)m_ResultProducer).
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
    if (m_ResultProducer instanceof AdditionalMeasureProducer) {
      return ((AdditionalMeasureProducer)m_ResultProducer).
	getMeasure(additionalMeasureName);
    } else {
      throw new IllegalArgumentException("AveragingResultProducer: "
			  +"Can't return value for : "+additionalMeasureName
			  +". "+m_ResultProducer.getClass().getName()+" "
			  +"is not an AdditionalMeasureProducer");
    }
  }

  /**
   * Sets the dataset that results will be obtained for.
   *
   * @param instances a value of type 'Instances'.
   */
  public void setInstances(Instances instances) {
    
    m_Instances = instances;
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String calculateStdDevsTipText() {
    return "Record standard deviations for each run.";
  }

  /**
   * Get the value of CalculateStdDevs.
   *
   * @return Value of CalculateStdDevs.
   */
  public boolean getCalculateStdDevs() {
    
    return m_CalculateStdDevs;
  }
  
  /**
   * Set the value of CalculateStdDevs.
   *
   * @param newCalculateStdDevs Value to assign to CalculateStdDevs.
   */
  public void setCalculateStdDevs(boolean newCalculateStdDevs) {
    
    m_CalculateStdDevs = newCalculateStdDevs;
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String expectedResultsPerAverageTipText() {
    return "Set the expected number of results to average per run. "
      +"For example if a CrossValidationResultProducer is being used "
      +"(with the number of folds set to 10), then the expected number "
      +"of results per run is 10.";
  }

  /**
   * Get the value of ExpectedResultsPerAverage.
   *
   * @return Value of ExpectedResultsPerAverage.
   */
  public int getExpectedResultsPerAverage() {
    
    return m_ExpectedResultsPerAverage;
  }
  
  /**
   * Set the value of ExpectedResultsPerAverage.
   *
   * @param newExpectedResultsPerAverage Value to assign to
   * ExpectedResultsPerAverage.
   */
  public void setExpectedResultsPerAverage(int newExpectedResultsPerAverage) {
    
    m_ExpectedResultsPerAverage = newExpectedResultsPerAverage;
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String keyFieldNameTipText() {
    return "Set the field name that will be unique for a run.";
  }

  /**
   * Get the value of KeyFieldName.
   *
   * @return Value of KeyFieldName.
   */
  public String getKeyFieldName() {
    
    return m_KeyFieldName;
  }
  
  /**
   * Set the value of KeyFieldName.
   *
   * @param newKeyFieldName Value to assign to KeyFieldName.
   */
  public void setKeyFieldName(String newKeyFieldName) {
    
    m_KeyFieldName = newKeyFieldName;
    m_CountFieldName = "Num_" + m_KeyFieldName;
    findKeyIndex();
  }
  
  /**
   * Sets the object to send results of each run to.
   *
   * @param listener a value of type 'ResultListener'
   */
  public void setResultListener(ResultListener listener) {

    m_ResultListener = listener;
  }
  
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String resultProducerTipText() {
    return "Set the resultProducer for which results are to be averaged.";
  }

  /**
   * Get the ResultProducer.
   *
   * @return the ResultProducer.
   */
  public ResultProducer getResultProducer() {
    
    return m_ResultProducer;
  }
  
  /**
   * Set the ResultProducer.
   *
   * @param newResultProducer new ResultProducer to use.
   */
  public void setResultProducer(ResultProducer newResultProducer) {

    m_ResultProducer = newResultProducer;
    m_ResultProducer.setResultListener(this);
    findKeyIndex();
  }

  /**
   * Gets a text descrption of the result producer.
   *
   * @return a text description of the result producer.
   */
  public String toString() {

    String result = "AveragingResultProducer: ";
    result += getCompatibilityState();
    if (m_Instances == null) {
      result += ": <null Instances>";
    } else {
      result += ": " + Utils.backQuoteChars(m_Instances.relationName());
    }
    return result;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6420 $");
  }
} // AveragingResultProducer
