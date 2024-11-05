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
 *    CrossValidationResultProducer.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.experiment;

import weka.core.AdditionalMeasureProducer;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.core.Utils;

import java.io.File;
import java.util.Calendar;
import java.util.Enumeration;
import java.util.Random;
import java.util.TimeZone;
import java.util.Vector;


/**
 <!-- globalinfo-start -->
 * Generates for each run, carries out an n-fold cross-validation, using the set SplitEvaluator to generate some results. If the class attribute is nominal, the dataset is stratified. Results for each fold are generated, so you may wish to use this in addition with an AveragingResultProducer to obtain averages for each run.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
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
 * All options after -- will be passed to the split evaluator.
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 1.17 $
 */
public class CrossValidationResultProducer 
  implements ResultProducer, OptionHandler, AdditionalMeasureProducer, 
             RevisionHandler {
  
  /** for serialization */
  static final long serialVersionUID = -1580053925080091917L;
  
  /** The dataset of interest */
  protected Instances m_Instances;

  /** The ResultListener to send results to */
  protected ResultListener m_ResultListener = new CSVResultListener();

  /** The number of folds in the cross-validation */
  protected int m_NumFolds = 10;

  /** Save raw output of split evaluators --- for debugging purposes */
  protected boolean m_debugOutput = false;

  /** The output zipper to use for saving raw splitEvaluator output */
  protected OutputZipper m_ZipDest = null;

  /** The destination output file/directory for raw output */
  protected File m_OutputFile = new File(
				new File(System.getProperty("user.dir")), 
				"splitEvalutorOut.zip");

  /** The SplitEvaluator used to generate results */
  protected SplitEvaluator m_SplitEvaluator = new ClassifierSplitEvaluator();

  /** The names of any additional measures to look for in SplitEvaluators */
  protected String [] m_AdditionalMeasures = null;

  /** The name of the key field containing the dataset name */
  public static String DATASET_FIELD_NAME = "Dataset";

  /** The name of the key field containing the run number */
  public static String RUN_FIELD_NAME = "Run";

  /** The name of the key field containing the fold number */
  public static String FOLD_FIELD_NAME = "Fold";

  /** The name of the result field containing the timestamp */
  public static String TIMESTAMP_FIELD_NAME = "Date_time";

  /**
   * Returns a string describing this result producer
   * @return a description of the result producer suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "Generates for each run, carries out an n-fold cross-validation, "
      + "using the set SplitEvaluator to generate some results. If the class "
      + "attribute is nominal, the dataset is stratified. Results for each fold "
      + "are generated, so you may wish to use this in addition with an "
      + "AveragingResultProducer to obtain averages for each run.";
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
   * Set a list of method names for additional measures to look for
   * in SplitEvaluators. This could contain many measures (of which only a
   * subset may be produceable by the current SplitEvaluator) if an experiment
   * is the type that iterates over a set of properties.
   * @param additionalMeasures an array of measure names, null if none
   */
  public void setAdditionalMeasures(String [] additionalMeasures) {
    m_AdditionalMeasures = additionalMeasures;

    if (m_SplitEvaluator != null) {
      System.err.println("CrossValidationResultProducer: setting additional "
			 +"measures for "
			 +"split evaluator");
      m_SplitEvaluator.setAdditionalMeasures(m_AdditionalMeasures);
    }
  }

  /**
   * Returns an enumeration of any additional measure names that might be
   * in the SplitEvaluator
   * @return an enumeration of the measure names
   */
  public Enumeration enumerateMeasures() {
    Vector newVector = new Vector();
    if (m_SplitEvaluator instanceof AdditionalMeasureProducer) {
      Enumeration en = ((AdditionalMeasureProducer)m_SplitEvaluator).
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
    if (m_SplitEvaluator instanceof AdditionalMeasureProducer) {
      return ((AdditionalMeasureProducer)m_SplitEvaluator).
	getMeasure(additionalMeasureName);
    } else {
      throw new IllegalArgumentException("CrossValidationResultProducer: "
			  +"Can't return value for : "+additionalMeasureName
			  +". "+m_SplitEvaluator.getClass().getName()+" "
			  +"is not an AdditionalMeasureProducer");
    }
  }
  
  /**
   * Gets a Double representing the current date and time.
   * eg: 1:46pm on 20/5/1999 -> 19990520.1346
   *
   * @return a value of type Double
   */
  public static Double getTimestamp() {

    Calendar now = Calendar.getInstance(TimeZone.getTimeZone("UTC"));
    double timestamp = now.get(Calendar.YEAR) * 10000
      + (now.get(Calendar.MONTH) + 1) * 100
      + now.get(Calendar.DAY_OF_MONTH)
      + now.get(Calendar.HOUR_OF_DAY) / 100.0
      + now.get(Calendar.MINUTE) / 10000.0;
    return new Double(timestamp);
  }
  
  /**
   * Prepare to generate results.
   *
   * @throws Exception if an error occurs during preprocessing.
   */
  public void preProcess() throws Exception {

    if (m_SplitEvaluator == null) {
      throw new Exception("No SplitEvalutor set");
    }
    if (m_ResultListener == null) {
      throw new Exception("No ResultListener set");
    }
    m_ResultListener.preProcess(this);
  }
  
  /**
   * Perform any postprocessing. When this method is called, it indicates
   * that no more requests to generate results for the current experiment
   * will be sent.
   *
   * @throws Exception if an error occurs
   */
  public void postProcess() throws Exception {

    m_ResultListener.postProcess(this);

    if (m_debugOutput) {
      if (m_ZipDest != null) {
	m_ZipDest.finished();
	m_ZipDest = null;
      }
    }
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
    if (m_Instances == null) {
      throw new Exception("No Instances set");
    }
    /*    // Randomize on a copy of the original dataset
    Instances runInstances = new Instances(m_Instances);
    runInstances.randomize(new Random(run));
    if (runInstances.classAttribute().isNominal()) {
      runInstances.stratify(m_NumFolds);
      } */
    for (int fold = 0; fold < m_NumFolds; fold++) {
      // Add in some fields to the key like run and fold number, dataset name
      Object [] seKey = m_SplitEvaluator.getKey();
      Object [] key = new Object [seKey.length + 3];
      key[0] = Utils.backQuoteChars(m_Instances.relationName());
      key[1] = "" + run;
      key[2] = "" + (fold + 1);
      System.arraycopy(seKey, 0, key, 3, seKey.length);
      if (m_ResultListener.isResultRequired(this, key)) {
	try {
	  m_ResultListener.acceptResult(this, key, null);
	} catch (Exception ex) {
	  // Save the train and test datasets for debugging purposes?
	  throw ex;
	}
      }
    }
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

    if (getRawOutput()) {
      if (m_ZipDest == null) {
	m_ZipDest = new OutputZipper(m_OutputFile);
      }
    }

    if (m_Instances == null) {
      throw new Exception("No Instances set");
    }
    // Randomize on a copy of the original dataset
    Instances runInstances = new Instances(m_Instances);
    Random random = new Random(run);
    runInstances.randomize(random);
    if (runInstances.classAttribute().isNominal() || runInstances.classAttribute().isRanking()) {
      runInstances.stratify(m_NumFolds);
    }
    for (int fold = 0; fold < m_NumFolds; fold++) {
      // Add in some fields to the key like run and fold number, dataset name
      Object [] seKey = m_SplitEvaluator.getKey();
      Object [] key = new Object [seKey.length + 3];
      key[0] =  Utils.backQuoteChars(m_Instances.relationName());
      key[1] = "" + run;
      key[2] = "" + (fold + 1);
      System.arraycopy(seKey, 0, key, 3, seKey.length);
      if (m_ResultListener.isResultRequired(this, key)) {
	Instances train = runInstances.trainCV(m_NumFolds, fold, random);
	Instances test = runInstances.testCV(m_NumFolds, fold);
	try {
	  Object [] seResults = m_SplitEvaluator.getResult(train, test);
	  Object [] results = new Object [seResults.length + 1];
	  results[0] = getTimestamp();
	  System.arraycopy(seResults, 0, results, 1,
			   seResults.length);
	  if (m_debugOutput) {
	    String resultName = (""+run+"."+(fold+1)+"."
	      + Utils.backQuoteChars(runInstances.relationName())
	      +"."
	      +m_SplitEvaluator.toString()).replace(' ','_');
	    resultName = Utils.removeSubstring(resultName, 
					       "weka.classifiers.");
	    resultName = Utils.removeSubstring(resultName, 
					       "weka.filters.");
	    resultName = Utils.removeSubstring(resultName, 
					       "weka.attributeSelection.");
	    m_ZipDest.zipit(m_SplitEvaluator.getRawResultOutput(), resultName);
	  }
	  m_ResultListener.acceptResult(this, key, results);
	} catch (Exception ex) {
	  // Save the train and test datasets for debugging purposes?
	  throw ex;
	}
      }
    }
  }

  /**
   * Gets the names of each of the columns produced for a single run.
   * This method should really be static.
   *
   * @return an array containing the name of each column
   */
  public String [] getKeyNames() {

    String [] keyNames = m_SplitEvaluator.getKeyNames();
    // Add in the names of our extra key fields
    String [] newKeyNames = new String [keyNames.length + 3];
    newKeyNames[0] = DATASET_FIELD_NAME;
    newKeyNames[1] = RUN_FIELD_NAME;
    newKeyNames[2] = FOLD_FIELD_NAME;
    System.arraycopy(keyNames, 0, newKeyNames, 3, keyNames.length);
    return newKeyNames;
  }

  /**
   * Gets the data types of each of the columns produced for a single run.
   * This method should really be static.
   *
   * @return an array containing objects of the type of each column. The 
   * objects should be Strings, or Doubles.
   */
  public Object [] getKeyTypes() {

    Object [] keyTypes = m_SplitEvaluator.getKeyTypes();
    // Add in the types of our extra fields
    Object [] newKeyTypes = new String [keyTypes.length + 3];
    newKeyTypes[0] = new String();
    newKeyTypes[1] = new String();
    newKeyTypes[2] = new String();
    System.arraycopy(keyTypes, 0, newKeyTypes, 3, keyTypes.length);
    return newKeyTypes;
  }

  /**
   * Gets the names of each of the columns produced for a single run.
   * This method should really be static.
   *
   * @return an array containing the name of each column
   */
  public String [] getResultNames() {

    String [] resultNames = m_SplitEvaluator.getResultNames();
    // Add in the names of our extra Result fields
    String [] newResultNames = new String [resultNames.length + 1];
    newResultNames[0] = TIMESTAMP_FIELD_NAME;
    System.arraycopy(resultNames, 0, newResultNames, 1, resultNames.length);
    return newResultNames;
  }

  /**
   * Gets the data types of each of the columns produced for a single run.
   * This method should really be static.
   *
   * @return an array containing objects of the type of each column. The 
   * objects should be Strings, or Doubles.
   */
  public Object [] getResultTypes() {

    Object [] resultTypes = m_SplitEvaluator.getResultTypes();
    // Add in the types of our extra Result fields
    Object [] newResultTypes = new Object [resultTypes.length + 1];
    newResultTypes[0] = new Double(0);
    System.arraycopy(resultTypes, 0, newResultTypes, 1, resultTypes.length);
    return newResultTypes;
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

    String result = "-X " + m_NumFolds + " " ;
    if (m_SplitEvaluator == null) {
      result += "<null SplitEvaluator>";
    } else {
      result += "-W " + m_SplitEvaluator.getClass().getName();
    }
    return result + " --";
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String outputFileTipText() {
    return "Set the destination for saving raw output. If the rawOutput "
      +"option is selected, then output from the splitEvaluator for "
      +"individual folds is saved. If the destination is a directory, "
      +"then each output is saved to an individual gzip file; if the "
      +"destination is a file, then each output is saved as an entry "
      +"in a zip file.";
  }

  /**
   * Get the value of OutputFile.
   *
   * @return Value of OutputFile.
   */
  public File getOutputFile() {
    
    return m_OutputFile;
  }
  
  /**
   * Set the value of OutputFile.
   *
   * @param newOutputFile Value to assign to OutputFile.
   */
  public void setOutputFile(File newOutputFile) {
    
    m_OutputFile = newOutputFile;
  }  

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String numFoldsTipText() {
    return "Number of folds to use in cross validation.";
  }

  /**
   * Get the value of NumFolds.
   *
   * @return Value of NumFolds.
   */
  public int getNumFolds() {
    
    return m_NumFolds;
  }
  
  /**
   * Set the value of NumFolds.
   *
   * @param newNumFolds Value to assign to NumFolds.
   */
  public void setNumFolds(int newNumFolds) {
    
    m_NumFolds = newNumFolds;
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String rawOutputTipText() {
    return "Save raw output (useful for debugging). If set, then output is "
      +"sent to the destination specified by outputFile";
  }

  /**
   * Get if raw split evaluator output is to be saved
   * @return true if raw split evalutor output is to be saved
   */
  public boolean getRawOutput() {
    return m_debugOutput;
  }
  
  /**
   * Set to true if raw split evaluator output is to be saved
   * @param d true if output is to be saved
   */
  public void setRawOutput(boolean d) {
    m_debugOutput = d;
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String splitEvaluatorTipText() {
    return "The evaluator to apply to the cross validation folds. "
      +"This may be a classifier, regression scheme etc.";
  }
 
  /**
   * Get the SplitEvaluator.
   *
   * @return the SplitEvaluator.
   */
  public SplitEvaluator getSplitEvaluator() {
    
    return m_SplitEvaluator;
  }
  
  /**
   * Set the SplitEvaluator.
   *
   * @param newSplitEvaluator new SplitEvaluator to use.
   */
  public void setSplitEvaluator(SplitEvaluator newSplitEvaluator) {
    
    m_SplitEvaluator = newSplitEvaluator;
    m_SplitEvaluator.setAdditionalMeasures(m_AdditionalMeasures);
  }

  /**
   * Returns an enumeration describing the available options..
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {

    Vector newVector = new Vector(4);

    newVector.addElement(new Option(
	     "\tThe number of folds to use for the cross-validation.\n"
	      +"\t(default 10)", 
	     "X", 1, 
	     "-X <number of folds>"));

    newVector.addElement(new Option(
	     "Save raw split evaluator output.",
	     "D",0,"-D"));

    newVector.addElement(new Option(
	     "\tThe filename where raw output will be stored.\n"
	     +"\tIf a directory name is specified then then individual\n"
	     +"\toutputs will be gzipped, otherwise all output will be\n"
	     +"\tzipped to the named file. Use in conjuction with -D."
	     +"\t(default splitEvalutorOut.zip)", 
	     "O", 1, 
	     "-O <file/directory name/path>"));

    newVector.addElement(new Option(
	     "\tThe full class name of a SplitEvaluator.\n"
	      +"\teg: weka.experiment.ClassifierSplitEvaluator", 
	     "W", 1, 
	     "-W <class name>"));

    if ((m_SplitEvaluator != null) &&
	(m_SplitEvaluator instanceof OptionHandler)) {
      newVector.addElement(new Option(
	     "",
	     "", 0, "\nOptions specific to split evaluator "
	     + m_SplitEvaluator.getClass().getName() + ":"));
      Enumeration enu = ((OptionHandler)m_SplitEvaluator).listOptions();
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
   * All options after -- will be passed to the split evaluator.
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    
    setRawOutput(Utils.getFlag('D', options));

    String fName = Utils.getOption('O', options);
    if (fName.length() != 0) {
      setOutputFile(new File(fName));
    }

    String numFolds = Utils.getOption('X', options);
    if (numFolds.length() != 0) {
      setNumFolds(Integer.parseInt(numFolds));
    } else {
      setNumFolds(10);
    }

    String seName = Utils.getOption('W', options);
    if (seName.length() == 0) {
      throw new Exception("A SplitEvaluator must be specified with"
			  + " the -W option.");
    }
    // Do it first without options, so if an exception is thrown during
    // the option setting, listOptions will contain options for the actual
    // SE.
    setSplitEvaluator((SplitEvaluator)Utils.forName(
		      SplitEvaluator.class,
		      seName,
		      null));
    if (getSplitEvaluator() instanceof OptionHandler) {
      ((OptionHandler) getSplitEvaluator())
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
    if ((m_SplitEvaluator != null) && 
	(m_SplitEvaluator instanceof OptionHandler)) {
      seOptions = ((OptionHandler)m_SplitEvaluator).getOptions();
    }
    
    String [] options = new String [seOptions.length + 8];
    int current = 0;

    options[current++] = "-X"; options[current++] = "" + getNumFolds();

    if (getRawOutput()) {
      options[current++] = "-D";
    }

    options[current++] = "-O"; 
    options[current++] = getOutputFile().getName();
    
    if (getSplitEvaluator() != null) {
      options[current++] = "-W";
      options[current++] = getSplitEvaluator().getClass().getName();
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
   * Gets a text descrption of the result producer.
   *
   * @return a text description of the result producer.
   */
  public String toString() {

    String result = "CrossValidationResultProducer: ";
    result += getCompatibilityState();
    if (m_Instances == null) {
      result += ": <null Instances>";
    } else {
      result += ": " +  Utils.backQuoteChars(m_Instances.relationName());
    }
    return result;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 1.17 $");
  }
    
  /** 
   * Quick test of timestamp
   * 
   * @param args	the commandline options
   */
  public static void main(String [] args) {
    
    System.err.println(Utils.doubleToString(getTimestamp().doubleValue(), 4));
  }
} // CrossValidationResultProducer
