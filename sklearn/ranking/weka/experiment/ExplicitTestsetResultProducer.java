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
 *    ExplicitTestsetResultProducer.java
 *    Copyright (C) 2009 University of Waikato, Hamilton, New Zealand
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
import weka.core.WekaException;
import weka.core.converters.ConverterUtils.DataSource;

import java.io.File;
import java.util.Calendar;
import java.util.Enumeration;
import java.util.Random;
import java.util.TimeZone;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * Loads the external test set and calls the appropriate SplitEvaluator to generate some results.<br/>
 * The filename of the test set is constructed as follows:<br/>
 *    &lt;dir&gt; + / + &lt;prefix&gt; + &lt;relation-name&gt; + &lt;suffix&gt;<br/>
 * The relation-name can be modified by using the regular expression to replace the matching sub-string with a specified replacement string. In order to get rid of the string that the Weka filters add to the end of the relation name, just use '.*-weka' as the regular expression to find.<br/>
 * The suffix determines the type of file to load, i.e., one is not restricted to ARFF files. As long as Weka recognizes the extension specified in the suffix, the data will be loaded with one of Weka's converters.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -D
 * Save raw split evaluator output.</pre>
 * 
 * <pre> -O &lt;file/directory name/path&gt;
 *  The filename where raw output will be stored.
 *  If a directory name is specified then then individual
 *  outputs will be gzipped, otherwise all output will be
 *  zipped to the named file. Use in conjuction with -D.
 *  (default: splitEvalutorOut.zip)</pre>
 * 
 * <pre> -W &lt;class name&gt;
 *  The full class name of a SplitEvaluator.
 *  eg: weka.experiment.ClassifierSplitEvaluator</pre>
 * 
 * <pre> -R
 *  Set when data is to be randomized.</pre>
 * 
 * <pre> -dir &lt;directory&gt;
 *  The directory containing the test sets.
 *  (default: current directory)</pre>
 * 
 * <pre> -prefix &lt;string&gt;
 *  An optional prefix for the test sets (before the relation name).
 * (default: empty string)</pre>
 * 
 * <pre> -suffix &lt;string&gt;
 *  The suffix to append to the test set.
 *  (default: _test.arff)</pre>
 * 
 * <pre> -find &lt;regular expression&gt;
 *  The regular expression to search the relation name with.
 *  Not used if an empty string.
 *  (default: empty string)</pre>
 * 
 * <pre> -replace &lt;string&gt;
 *  The replacement string for the all the matches of '-find'.
 *  (default: empty string)</pre>
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
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5353 $
 */
public class ExplicitTestsetResultProducer 
  implements ResultProducer, OptionHandler, AdditionalMeasureProducer, 
             RevisionHandler {
  
  /** for serialization. */
  private static final long serialVersionUID = 2613585409333652530L;

  /** the default suffix. */
  public final static String DEFAULT_SUFFIX = "_test.arff";
  
  /** The dataset of interest. */
  protected Instances m_Instances;

  /** The ResultListener to send results to. */
  protected ResultListener m_ResultListener = new CSVResultListener();

  /** The directory containing all the test sets. */
  protected File m_TestsetDir = new File(System.getProperty("user.dir"));

  /** The prefix for all the test sets. */
  protected String m_TestsetPrefix = "";

  /** The suffix for all the test sets. */
  protected String m_TestsetSuffix = DEFAULT_SUFFIX;

  /** The regular expression to search for in the relation name. */
  protected String m_RelationFind = "";

  /** The string to use to replace the matches of the regular expression. */
  protected String m_RelationReplace = "";

  /** Whether dataset is to be randomized. */
  protected boolean m_randomize = false;

  /** The SplitEvaluator used to generate results. */
  protected SplitEvaluator m_SplitEvaluator = new ClassifierSplitEvaluator();

  /** The names of any additional measures to look for in SplitEvaluators. */
  protected String[] m_AdditionalMeasures = null;

  /** Save raw output of split evaluators --- for debugging purposes. */
  protected boolean m_debugOutput = false;

  /** The output zipper to use for saving raw splitEvaluator output. */
  protected OutputZipper m_ZipDest = null;

  /** The destination output file/directory for raw output. */
  protected File m_OutputFile = new File(
			        new File(System.getProperty("user.dir")), 
				"splitEvalutorOut.zip");

  /** The name of the key field containing the dataset name. */
  public static String DATASET_FIELD_NAME = "Dataset";

  /** The name of the key field containing the run number. */
  public static String RUN_FIELD_NAME = "Run";

  /** The name of the result field containing the timestamp. */
  public static String TIMESTAMP_FIELD_NAME = "Date_time";

  /**
   * Returns a string describing this result producer.
   * 
   * @return 		a description of the result producer suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return
        "Loads the external test set and calls the appropriate "
      + "SplitEvaluator to generate some results.\n"
      + "The filename of the test set is constructed as follows:\n"
      + "   <dir> + / + <prefix> + <relation-name> + <suffix>\n"
      + "The relation-name can be modified by using the regular expression "
      + "to replace the matching sub-string with a specified replacement "
      + "string. In order to get rid of the string that the Weka filters "
      + "add to the end of the relation name, just use '.*-weka' as the "
      + "regular expression to find.\n"
      + "The suffix determines the type of file to load, i.e., one is "
      + "not restricted to ARFF files. As long as Weka recognizes the "
      + "extension specified in the suffix, the data will be loaded with "
      + "one of Weka's converters.";
  }

  /**
   * Returns an enumeration describing the available options..
   *
   * @return 		an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector result = new Vector();

    result.addElement(new Option(
	"Save raw split evaluator output.",
	"D", 0, "-D"));

    result.addElement(new Option(
	"\tThe filename where raw output will be stored.\n"
	+"\tIf a directory name is specified then then individual\n"
	+"\toutputs will be gzipped, otherwise all output will be\n"
	+"\tzipped to the named file. Use in conjuction with -D.\n"
	+"\t(default: splitEvalutorOut.zip)", 
	"O", 1, "-O <file/directory name/path>"));

    result.addElement(new Option(
	"\tThe full class name of a SplitEvaluator.\n"
	+"\teg: weka.experiment.ClassifierSplitEvaluator", 
	"W", 1, "-W <class name>"));

    result.addElement(new Option(
	"\tSet when data is to be randomized.",
	"R", 0 ,"-R"));

    result.addElement(new Option(
	"\tThe directory containing the test sets.\n"
	+ "\t(default: current directory)", 
	"dir", 1, "-dir <directory>"));

    result.addElement(new Option(
	"\tAn optional prefix for the test sets (before the relation name).\n"
	+ "(default: empty string)", 
	"prefix", 1, "-prefix <string>"));

    result.addElement(new Option(
	"\tThe suffix to append to the test set.\n"
	+ "\t(default: " + DEFAULT_SUFFIX + ")", 
	"suffix", 1, "-suffix <string>"));

    result.addElement(new Option(
	"\tThe regular expression to search the relation name with.\n"
	+ "\tNot used if an empty string.\n"
	+ "\t(default: empty string)", 
	"find", 1, "-find <regular expression>"));

    result.addElement(new Option(
	"\tThe replacement string for the all the matches of '-find'.\n"
	+ "\t(default: empty string)", 
	"replace", 1, "-replace <string>"));
    
    if ((m_SplitEvaluator != null) && (m_SplitEvaluator instanceof OptionHandler)) {
      result.addElement(new Option(
	  "",
	  "", 0, "\nOptions specific to split evaluator "
	  + m_SplitEvaluator.getClass().getName() + ":"));
      Enumeration enu = ((OptionHandler)m_SplitEvaluator).listOptions();
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
   * Save raw split evaluator output.</pre>
   * 
   * <pre> -O &lt;file/directory name/path&gt;
   *  The filename where raw output will be stored.
   *  If a directory name is specified then then individual
   *  outputs will be gzipped, otherwise all output will be
   *  zipped to the named file. Use in conjuction with -D.
   *  (default: splitEvalutorOut.zip)</pre>
   * 
   * <pre> -W &lt;class name&gt;
   *  The full class name of a SplitEvaluator.
   *  eg: weka.experiment.ClassifierSplitEvaluator</pre>
   * 
   * <pre> -R
   *  Set when data is to be randomized.</pre>
   * 
   * <pre> -dir &lt;directory&gt;
   *  The directory containing the test sets.
   *  (default: current directory)</pre>
   * 
   * <pre> -prefix &lt;string&gt;
   *  An optional prefix for the test sets (before the relation name).
   * (default: empty string)</pre>
   * 
   * <pre> -suffix &lt;string&gt;
   *  The suffix to append to the test set.
   *  (default: _test.arff)</pre>
   * 
   * <pre> -find &lt;regular expression&gt;
   *  The regular expression to search the relation name with.
   *  Not used if an empty string.
   *  (default: empty string)</pre>
   * 
   * <pre> -replace &lt;string&gt;
   *  The replacement string for the all the matches of '-find'.
   *  (default: empty string)</pre>
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
   * @param options 	the list of options as an array of strings
   * @throws Exception 	if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String	tmpStr;
    
    setRawOutput(Utils.getFlag('D', options));
    setRandomizeData(!Utils.getFlag('R', options));

    tmpStr = Utils.getOption('O', options);
    if (tmpStr.length() != 0)
      setOutputFile(new File(tmpStr));

    tmpStr = Utils.getOption("dir", options);
    if (tmpStr.length() > 0)
      setTestsetDir(new File(tmpStr));
    else
      setTestsetDir(new File(System.getProperty("user.dir")));

    tmpStr = Utils.getOption("prefix", options);
    if (tmpStr.length() > 0)
      setTestsetPrefix(tmpStr);
    else
      setTestsetPrefix("");

    tmpStr = Utils.getOption("suffix", options);
    if (tmpStr.length() > 0)
      setTestsetSuffix(tmpStr);
    else
      setTestsetSuffix(DEFAULT_SUFFIX);
    
    tmpStr = Utils.getOption("find", options);
    if (tmpStr.length() > 0)
      setRelationFind(tmpStr);
    else
      setRelationFind("");
    
    tmpStr = Utils.getOption("replace", options);
    if ((tmpStr.length() > 0) && (getRelationFind().length() > 0))
      setRelationReplace(tmpStr);
    else
      setRelationReplace("");
    
    tmpStr = Utils.getOption('W', options);
    if (tmpStr.length() == 0)
      throw new Exception("A SplitEvaluator must be specified with the -W option.");
    
    // Do it first without options, so if an exception is thrown during
    // the option setting, listOptions will contain options for the actual
    // SE.
    setSplitEvaluator((SplitEvaluator)Utils.forName(SplitEvaluator.class, tmpStr, null));
    if (getSplitEvaluator() instanceof OptionHandler)
      ((OptionHandler) getSplitEvaluator()).setOptions(Utils.partitionOptions(options));
  }

  /**
   * Gets the current settings of the result producer.
   *
   * @return 		an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    Vector<String>	result;
    String[] 		seOptions;
    int			i;
    
    result = new Vector<String>();
    
    seOptions = new String [0];
    if ((m_SplitEvaluator != null) && (m_SplitEvaluator instanceof OptionHandler))
      seOptions = ((OptionHandler)m_SplitEvaluator).getOptions();

    if (getRawOutput())
      result.add("-D");
    
    if (!getRandomizeData())
      result.add("-R");

    result.add("-O"); 
    result.add(getOutputFile().getName());

    result.add("-dir");
    result.add(getTestsetDir().getPath());
    
    if (getTestsetPrefix().length() > 0) {
      result.add("-prefix");
      result.add(getTestsetPrefix());
    }

    result.add("-suffix");
    result.add(getTestsetSuffix());
    
    if (getRelationFind().length() > 0) {
      result.add("-find");
      result.add(getRelationFind());
      
      if (getRelationReplace().length() > 0) {
	result.add("-replace");
	result.add(getRelationReplace());
      }
    }

    if (getSplitEvaluator() != null) {
      result.add("-W");
      result.add(getSplitEvaluator().getClass().getName());
    }
    
    if (seOptions.length > 0) {
      result.add("--");
      for (i = 0; i < seOptions.length; i++)
	result.add(seOptions[i]);
    }

    return result.toArray(new String[result.size()]);
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
   * Set a list of method names for additional measures to look for
   * in SplitEvaluators. This could contain many measures (of which only a
   * subset may be produceable by the current SplitEvaluator) if an experiment
   * is the type that iterates over a set of properties.
   * 
   * @param additionalMeasures 	an array of measure names, null if none
   */
  public void setAdditionalMeasures(String[] additionalMeasures) {
    m_AdditionalMeasures = additionalMeasures;

    if (m_SplitEvaluator != null) {
      System.err.println(
	  "ExplicitTestsetResultProducer: setting additional "
	  + "measures for split evaluator");
      m_SplitEvaluator.setAdditionalMeasures(m_AdditionalMeasures);
    }
  }
  
  /**
   * Returns an enumeration of any additional measure names that might be
   * in the SplitEvaluator.
   * 
   * @return 		an enumeration of the measure names
   */
  public Enumeration enumerateMeasures() {
    Vector result = new Vector();
    if (m_SplitEvaluator instanceof AdditionalMeasureProducer) {
      Enumeration en = ((AdditionalMeasureProducer)m_SplitEvaluator).enumerateMeasures();
      while (en.hasMoreElements()) {
	String mname = (String) en.nextElement();
	result.addElement(mname);
      }
    }
    return result.elements();
  }
  
  /**
   * Returns the value of the named measure.
   * 
   * @param additionalMeasureName 	the name of the measure to query for its value
   * @return 				the value of the named measure
   * @throws IllegalArgumentException 	if the named measure is not supported
   */
  public double getMeasure(String additionalMeasureName) {
    if (m_SplitEvaluator instanceof AdditionalMeasureProducer)
      return ((AdditionalMeasureProducer)m_SplitEvaluator).getMeasure(additionalMeasureName);
    else
      throw new IllegalArgumentException(
	  "ExplicitTestsetResultProducer: "
	  + "Can't return value for : " + additionalMeasureName
	  + ". " + m_SplitEvaluator.getClass().getName() + " "
	  + "is not an AdditionalMeasureProducer");
  }
  
  /**
   * Sets the object to send results of each run to.
   *
   * @param listener 	a value of type 'ResultListener'
   */
  public void setResultListener(ResultListener listener) {
    m_ResultListener = listener;
  }

  /**
   * Gets a Double representing the current date and time.
   * eg: 1:46pm on 20/5/1999 -> 19990520.1346
   *
   * @return 		a value of type Double
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
   * @throws Exception 	if an error occurs during preprocessing.
   */
  public void preProcess() throws Exception {
    if (m_SplitEvaluator == null)
      throw new Exception("No SplitEvalutor set");

    if (m_ResultListener == null)
      throw new Exception("No ResultListener set");

    m_ResultListener.preProcess(this);
  }
  
  /**
   * Perform any postprocessing. When this method is called, it indicates
   * that no more requests to generate results for the current experiment
   * will be sent.
   *
   * @throws Exception 	if an error occurs
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
   * @param run 	the run number to get keys for.
   * @throws Exception 	if a problem occurs while getting the keys
   */
  public void doRunKeys(int run) throws Exception {
    if (m_Instances == null)
      throw new Exception("No Instances set");

    // Add in some fields to the key like run number, dataset name
    Object[] seKey = m_SplitEvaluator.getKey();
    Object[] key = new Object [seKey.length + 2];
    key[0] = Utils.backQuoteChars(m_Instances.relationName());
    key[1] = "" + run;
    System.arraycopy(seKey, 0, key, 2, seKey.length);
    if (m_ResultListener.isResultRequired(this, key)) {
      try {
	m_ResultListener.acceptResult(this, key, null);
      }
      catch (Exception ex) {
	// Save the train and test datasets for debugging purposes?
	throw ex;
      }
    }
  }

  /**
   * Generates a new filename for the given relation based on the current 
   * setup.
   * 
   * @param inst	the instances to create the filename for
   * @return		the generated filename
   */
  protected String createFilename(Instances inst) {
    String	result;
    String	name;

    name = inst.relationName();
    if (getRelationFind().length() > 0)
      name = name.replaceAll(getRelationFind(), getRelationReplace());
    
    result  = getTestsetDir().getPath() + File.separator;
    result += getTestsetPrefix() + name + getTestsetSuffix();
    
    return result;
  }
  
  /**
   * Gets the results for a specified run number. Different run
   * numbers correspond to different randomizations of the data. Results
   * produced should be sent to the current ResultListener
   *
   * @param run 	the run number to get results for.
   * @throws Exception 	if a problem occurs while getting the results
   */
  public void doRun(int run) throws Exception {
    if (getRawOutput()) {
      if (m_ZipDest == null)
	m_ZipDest = new OutputZipper(m_OutputFile);
    }

    if (m_Instances == null)
      throw new Exception("No Instances set");
    
    // Add in some fields to the key like run number, dataset name
    Object[] seKey = m_SplitEvaluator.getKey();
    Object[] key = new Object [seKey.length + 2];
    key[0] = Utils.backQuoteChars(m_Instances.relationName());
    key[1] = "" + run;
    System.arraycopy(seKey, 0, key, 2, seKey.length);
    if (m_ResultListener.isResultRequired(this, key)) {
      // training set
      Instances train = new Instances(m_Instances);
      if (m_randomize) {
	Random rand = new Random(run);
	train.randomize(rand);
      }

      // test set
      String filename = createFilename(train);
      File file = new File(filename);
      if (!file.exists())
	throw new WekaException("Test set '" + filename + "' not found!");
      Instances test = DataSource.read(filename);
      // can we set the class attribute safely?
      if (train.numAttributes() == test.numAttributes())
	test.setClassIndex(train.classIndex());
      else
	throw new WekaException(
	    "Train and test set (= " + filename + ") "
	    + "differ in number of attributes: "
	    + train.numAttributes() + " != " + test.numAttributes());
      // test headers
      if (!train.equalHeaders(test))
	throw new WekaException(
	    "Train and test set (= " + filename + ") "
	    + "are not compatible:\n"
	    + train.equalHeadersMsg(test));
      
      try {
	Object[] seResults = m_SplitEvaluator.getResult(train, test);
	Object[] results = new Object [seResults.length + 1];
	results[0] = getTimestamp();
	System.arraycopy(seResults, 0, results, 1,
			 seResults.length);
	if (m_debugOutput) {
	  String resultName = 
	    (""+run+"."+
	     Utils.backQuoteChars(train.relationName())
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
      }
      catch (Exception e) {
	// Save the train and test datasets for debugging purposes?
	throw e;
      }
    }
  }

  /**
   * Gets the names of each of the columns produced for a single run.
   * This method should really be static.
   *
   * @return 		an array containing the name of each column
   */
  public String[] getKeyNames() {
    String[] keyNames = m_SplitEvaluator.getKeyNames();
    // Add in the names of our extra key fields
    String[] newKeyNames = new String [keyNames.length + 2];
    newKeyNames[0] = DATASET_FIELD_NAME;
    newKeyNames[1] = RUN_FIELD_NAME;
    System.arraycopy(keyNames, 0, newKeyNames, 2, keyNames.length);
    return newKeyNames;
  }

  /**
   * Gets the data types of each of the columns produced for a single run.
   * This method should really be static.
   *
   * @return 		an array containing objects of the type of each column. 
   * 			The objects should be Strings, or Doubles.
   */
  public Object[] getKeyTypes() {
    Object[] keyTypes = m_SplitEvaluator.getKeyTypes();
    // Add in the types of our extra fields
    Object[] newKeyTypes = new String [keyTypes.length + 2];
    newKeyTypes[0] = new String();
    newKeyTypes[1] = new String();
    System.arraycopy(keyTypes, 0, newKeyTypes, 2, keyTypes.length);
    return newKeyTypes;
  }

  /**
   * Gets the names of each of the columns produced for a single run.
   * This method should really be static.
   *
   * @return 		an array containing the name of each column
   */
  public String[] getResultNames() {
    String[] resultNames = m_SplitEvaluator.getResultNames();
    // Add in the names of our extra Result fields
    String[] newResultNames = new String [resultNames.length + 1];
    newResultNames[0] = TIMESTAMP_FIELD_NAME;
    System.arraycopy(resultNames, 0, newResultNames, 1, resultNames.length);
    return newResultNames;
  }

  /**
   * Gets the data types of each of the columns produced for a single run.
   * This method should really be static.
   *
   * @return 		an array containing objects of the type of each column. 
   * 			The objects should be Strings, or Doubles.
   */
  public Object[] getResultTypes() {
    Object[] resultTypes = m_SplitEvaluator.getResultTypes();
    // Add in the types of our extra Result fields
    Object[] newResultTypes = new Object [resultTypes.length + 1];
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
   * @return 		the description of the ResultProducer state, or null
   * 			if no state is defined
   */
  public String getCompatibilityState() {
    String 	result;
    
    result = "";
    if (getRandomizeData())
      result += " -R";

    result += " -dir " + getTestsetDir();
    
    if (getTestsetPrefix().length() > 0)
      result += " -prefix " + getTestsetPrefix();
    
    result += " -suffix " + getTestsetSuffix();
    
    if (getRelationFind().length() > 0) {
      result += " -find " + getRelationFind();
      
      if (getRelationReplace().length() > 0)
        result += " -replace " + getRelationReplace();
    }
    
    if (m_SplitEvaluator == null)
      result += " <null SplitEvaluator>";
    else
      result += " -W " + m_SplitEvaluator.getClass().getName();

    return result + " --";
  }

  /**
   * Returns the tip text for this property.
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String outputFileTipText() {
    return "Set the destination for saving raw output. If the rawOutput "
      +"option is selected, then output from the splitEvaluator for "
      +"individual train-test splits is saved. If the destination is a "
      +"directory, "
      +"then each output is saved to an individual gzip file; if the "
      +"destination is a file, then each output is saved as an entry "
      +"in a zip file.";
  }

  /**
   * Get the value of OutputFile.
   *
   * @return 		Value of OutputFile.
   */
  public File getOutputFile() {
    return m_OutputFile;
  }
  
  /**
   * Set the value of OutputFile.
   *
   * @param value 	Value to assign to OutputFile.
   */
  public void setOutputFile(File value) {
    m_OutputFile = value;
  }  

  /**
   * Returns the tip text for this property.
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String randomizeDataTipText() {
    return "Do not randomize dataset and do not perform probabilistic rounding " +
      "if true";
  }

  /**
   * Get if dataset is to be randomized.
   * 
   * @return 		true if dataset is to be randomized
   */
  public boolean getRandomizeData() {
    return m_randomize;
  }
  
  /**
   * Set to true if dataset is to be randomized.
   * 
   * @param value 		true if dataset is to be randomized
   */
  public void setRandomizeData(boolean value) {
    m_randomize = value;
  }

  /**
   * Returns the tip text for this property.
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String rawOutputTipText() {
    return "Save raw output (useful for debugging). If set, then output is "
      +"sent to the destination specified by outputFile";
  }

  /**
   * Get if raw split evaluator output is to be saved.
   * 
   * @return 		true if raw split evalutor output is to be saved
   */
  public boolean getRawOutput() {
    return m_debugOutput;
  }
  
  /**
   * Set to true if raw split evaluator output is to be saved.
   * 
   * @param value 		true if output is to be saved
   */
  public void setRawOutput(boolean value) {
    m_debugOutput = value;
  }

  /**
   * Returns the tip text for this property.
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String splitEvaluatorTipText() {
    return "The evaluator to apply to the test data. "
      +"This may be a classifier, regression scheme etc.";
  }

  /**
   * Get the SplitEvaluator.
   *
   * @return 		the SplitEvaluator.
   */
  public SplitEvaluator getSplitEvaluator() {
    return m_SplitEvaluator;
  }
  
  /**
   * Set the SplitEvaluator.
   *
   * @param value 	new SplitEvaluator to use.
   */
  public void setSplitEvaluator(SplitEvaluator value) {
    m_SplitEvaluator = value;
    m_SplitEvaluator.setAdditionalMeasures(m_AdditionalMeasures);
  }

  /**
   * Returns the tip text for this property.
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String testsetDirTipText() {
    return "The directory containing the test sets.";
  }

  /**
   * Returns the currently set directory for the test sets.
   *
   * @return 		the directory
   */
  public File getTestsetDir() {
    return m_TestsetDir;
  }
  
  /**
   * Sets the directory to use for the test sets.
   *
   * @param value 	the directory to use
   */
  public void setTestsetDir(File value) {
    m_TestsetDir = value;
  }  

  /**
   * Returns the tip text for this property.
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String testsetPrefixTipText() {
    return "The prefix to use for the filename of the test sets.";
  }

  /**
   * Returns the currently set prefix.
   *
   * @return 		the prefix
   */
  public String getTestsetPrefix() {
    return m_TestsetPrefix;
  }
  
  /**
   * Sets the prefix to use for the test sets.
   *
   * @param value 	the prefix
   */
  public void setTestsetPrefix(String value) {
    m_TestsetPrefix = value;
  }  

  /**
   * Returns the tip text for this property.
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String testsetSuffixTipText() {
    return 
        "The suffix to use for the filename of the test sets - must contain "
      + "the file extension.";
  }

  /**
   * Returns the currently set suffix.
   *
   * @return 		the suffix
   */
  public String getTestsetSuffix() {
    return m_TestsetSuffix;
  }
  
  /**
   * Sets the suffix to use for the test sets.
   *
   * @param value 	the suffix
   */
  public void setTestsetSuffix(String value) {
    if ((value == null) || (value.length() == 0))
      value = DEFAULT_SUFFIX;
    m_TestsetSuffix = value;
  }  

  /**
   * Returns the tip text for this property.
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String relationFindTipText() {
    return 
        "The regular expression to use for removing parts of the relation "
      + "name, ignored if empty.";
  }

  /**
   * Returns the currently set regular expression to use on the relation name.
   *
   * @return 		the regular expression
   */
  public String getRelationFind() {
    return m_RelationFind;
  }
  
  /**
   * Sets the regular expression to use on the relation name.
   *
   * @param value 	the regular expression
   */
  public void setRelationFind(String value) {
    m_RelationFind = value;
  }  

  /**
   * Returns the tip text for this property.
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String relationReplaceTipText() {
    return "The string to replace all matches of the regular expression with.";
  }

  /**
   * Returns the currently set replacement string to use on the relation name.
   *
   * @return 		the replacement string
   */
  public String getRelationReplace() {
    return m_RelationReplace;
  }
  
  /**
   * Sets the replacement string to use on the relation name.
   *
   * @param value 	the regular expression
   */
  public void setRelationReplace(String value) {
    m_RelationReplace = value;
  }  

  /**
   * Gets a text descrption of the result producer.
   *
   * @return 		a text description of the result producer.
   */
  public String toString() {
    String result = "ExplicitTestsetResultProducer: ";
    result += getCompatibilityState();
    if (m_Instances == null)
      result += ": <null Instances>";
    else
      result += ": " + Utils.backQuoteChars(m_Instances.relationName());
    return result;
  }

  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5353 $");
  }
}
