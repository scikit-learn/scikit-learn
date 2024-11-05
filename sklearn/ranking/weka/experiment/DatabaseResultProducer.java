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
 *    DatabaseResultProducer.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.experiment;

import weka.core.AdditionalMeasureProducer;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.Utils;

import java.util.Enumeration;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * Examines a database and extracts out the results produced by the specified ResultProducer and submits them to the specified ResultListener. If a result needs to be generated, the ResultProducer is used to obtain the result.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -F &lt;field name&gt;
 *  The name of the database field to cache over.
 *  eg: "Fold" (default none)</pre>
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
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 6423 $
 */
public class DatabaseResultProducer 
  extends DatabaseResultListener
  implements ResultProducer, OptionHandler, AdditionalMeasureProducer {

  /** for serialization */
  static final long serialVersionUID = -5620660780203158666L;
  
  /** The dataset of interest */
  protected Instances m_Instances;

  /** The ResultListener to send results to */
  protected ResultListener m_ResultListener = new CSVResultListener();

  /** The ResultProducer used to generate results */
  protected ResultProducer m_ResultProducer
    = new CrossValidationResultProducer();

  /** The names of any additional measures to look for in SplitEvaluators */
  protected String [] m_AdditionalMeasures = null;

  /**
   * Returns a string describing this result producer
   * @return a description of the result producer suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return "Examines a database and extracts out "
      +"the results produced by the specified ResultProducer "
      +"and submits them to the specified ResultListener. If a result needs "
      +"to be generated, the ResultProducer is used to obtain the result.";
  }

  /**
   * Creates the DatabaseResultProducer, letting the parent constructor do
   * it's thing.
   *
   * @throws Exception if an error occurs
   */
  public DatabaseResultProducer() throws Exception {

    super();
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

    if (m_ResultProducer == null) {
      throw new Exception("No ResultProducer set");
    }
    if (m_ResultListener == null) {
      throw new Exception("No ResultListener set");
    }
    if (m_Instances == null) {
      throw new Exception("No Instances set");
    }

    // Tell the resultproducer to send results to us
    m_ResultProducer.setResultListener(this);
    m_ResultProducer.setInstances(m_Instances);
    m_ResultProducer.doRunKeys(run);
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

    if (m_ResultProducer == null) {
      throw new Exception("No ResultProducer set");
    }
    if (m_ResultListener == null) {
      throw new Exception("No ResultListener set");
    }
    if (m_Instances == null) {
      throw new Exception("No Instances set");
    }

    // Tell the resultproducer to send results to us
    m_ResultProducer.setResultListener(this);
    m_ResultProducer.setInstances(m_Instances);
    m_ResultProducer.doRun(run);

  }
  
  /**
   * Prepare for the results to be received.
   *
   * @param rp the ResultProducer that will generate the results
   * @throws Exception if an error occurs during preprocessing.
   */
  public void preProcess(ResultProducer rp) throws Exception {

    super.preProcess(rp);
    if (m_ResultListener == null) {
      throw new Exception("No ResultListener set");
    }
    m_ResultListener.preProcess(this);
  }

  /**
   * When this method is called, it indicates that no more results
   * will be sent that need to be grouped together in any way.
   *
   * @param rp the ResultProducer that generated the results
   * @throws Exception if an error occurs
   */
  public void postProcess(ResultProducer rp) throws Exception {

    super.postProcess(rp);
    m_ResultListener.postProcess(this);
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
    m_ResultProducer.setResultListener(this);
    m_ResultProducer.preProcess();
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
    //    System.err.println("DBRP::acceptResult");

    // Is the result needed by the listener?
    boolean isRequiredByListener = m_ResultListener.isResultRequired(this,
								     key);
    // Is the result already in the database?
    boolean isRequiredByDatabase = super.isResultRequired(rp, key);

    // Insert it into the database here
    if (isRequiredByDatabase) {
      // We could alternatively throw an exception if we only want values
      // that are already in the database
      if (result != null) {

	// null result could occur from a chain of doRunKeys calls
	super.acceptResult(rp, key, result);
      }
    }

    // Pass it on
    if (isRequiredByListener) {
      m_ResultListener.acceptResult(this, key, result);
    }
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
    //    System.err.println("DBRP::isResultRequired");

    // Is the result needed by the listener?
    boolean isRequiredByListener = m_ResultListener.isResultRequired(this,
								     key);
    // Is the result already in the database?
    boolean isRequiredByDatabase = super.isResultRequired(rp, key);

    if (!isRequiredByDatabase && isRequiredByListener) {
      // Pass the result through to the listener
      Object [] result = getResultFromTable(m_ResultsTableName,
					    rp, key);
      System.err.println("Got result from database: "
			 + DatabaseUtils.arrayToString(result));
      m_ResultListener.acceptResult(this, key, result);
      return false;
    }

    return (isRequiredByListener || isRequiredByDatabase);
  }

  /**
   * Gets the names of each of the columns produced for a single run.
   *
   * @return an array containing the name of each column
   * @throws Exception if something goes wrong.
   */
  public String [] getKeyNames() throws Exception {

    return m_ResultProducer.getKeyNames();
  }

  /**
   * Gets the data types of each of the columns produced for a single run.
   * This method should really be static.
   *
   * @return an array containing objects of the type of each column. The 
   * objects should be Strings, or Doubles.
   * @throws Exception if something goes wrong.
   */
  public Object [] getKeyTypes() throws Exception {

    return m_ResultProducer.getKeyTypes();
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
   * @throws Exception if something goes wrong.
   */
  public String [] getResultNames() throws Exception {

    return m_ResultProducer.getResultNames();
  }

  /**
   * Gets the data types of each of the columns produced for a single run.
   *
   * @return an array containing objects of the type of each column. The 
   * objects should be Strings, or Doubles.
   * @throws Exception if something goes wrong.
   */
  public Object [] getResultTypes() throws Exception {

    return m_ResultProducer.getResultTypes();
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

    String result = "";
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
	     "\tThe name of the database field to cache over.\n"
	      +"\teg: \"Fold\" (default none)", 
	     "F", 1, 
	     "-F <field name>"));
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
   *  The name of the database field to cache over.
   *  eg: "Fold" (default none)</pre>
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
   * All option after -- will be passed to the result producer.
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    
    setCacheKeyName(Utils.getOption('F', options));
    
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

    if (!getCacheKeyName().equals("")) {
      options[current++] = "-F";
      options[current++] = getCacheKeyName();
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
      System.err.println("DatabaseResultProducer: setting additional "
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
      throw new IllegalArgumentException("DatabaseResultProducer: "
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
    return "Set the result producer to use. If some results are not found "
      +"in the source database then this result producer is used to generate "
      +"them.";
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
  }

  /**
   * Gets a text descrption of the result producer.
   *
   * @return a text description of the result producer.
   */
  public String toString() {

    String result = "DatabaseResultProducer: ";
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
    return RevisionUtils.extract("$Revision: 6423 $");
  }
} // DatabaseResultProducer
