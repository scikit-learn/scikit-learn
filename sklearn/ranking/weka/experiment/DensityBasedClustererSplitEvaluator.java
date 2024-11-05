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
 *    DensityBasedClustererSplitEvaluator.java
 *    Copyright (C) 2008 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.experiment;

import weka.clusterers.ClusterEvaluation;
import weka.clusterers.Clusterer;
import weka.clusterers.AbstractClusterer;
import weka.clusterers.AbstractDensityBasedClusterer;
import weka.clusterers.DensityBasedClusterer;
import weka.clusterers.EM;
import weka.core.AdditionalMeasureProducer;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.filters.Filter;
import weka.filters.unsupervised.attribute.Remove;

import java.io.ObjectStreamClass;
import java.io.Serializable;
import java.util.Enumeration;
import java.util.Vector;

/**
 * A SplitEvaluator that produces results for a density based clusterer.
 *
 * -W classname <br>
 * Specify the full class name of the clusterer to evaluate. <p>
 *
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}org
 * @version $Revision: 5563 $
 */

public class DensityBasedClustererSplitEvaluator 
  implements SplitEvaluator,
	     OptionHandler,
	     AdditionalMeasureProducer,
	     RevisionHandler {

  /** Remove the class column (if set) from the data */
  protected boolean m_removeClassColumn = true;

  /** The clusterer used for evaluation */
  protected DensityBasedClusterer m_clusterer = new EM();

  /** The names of any additional measures to look for in SplitEvaluators */
  protected String [] m_additionalMeasures = null;

  /** Array of booleans corresponding to the measures in m_AdditionalMeasures
      indicating which of the AdditionalMeasures the current clusterer
      can produce */
  protected boolean [] m_doesProduce = null;

  /** The number of additional measures that need to be filled in
      after taking into account column constraints imposed by the final
      destination for results */
  protected int m_numberAdditionalMeasures = 0;

  /** Holds the statistics for the most recent application of the clusterer */
  protected String m_result = null;

  /** The clusterer options (if any) */
  protected String m_clustererOptions = "";

  /** The clusterer version */
  protected String m_clustererVersion = "";

  /** The length of a key */
  private static final int KEY_SIZE = 3;

  /** The length of a result */
  private static final int RESULT_SIZE = 6;

  
  public DensityBasedClustererSplitEvaluator() {
    updateOptions();
  }

  /**
   * Returns a string describing this split evaluator
   * @return a description of the split evaluator suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return " A SplitEvaluator that produces results for a density based clusterer. ";
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {

    Vector newVector = new Vector(1);

    newVector.addElement(new Option(
				    "\tThe full class name of the density based clusterer.\n"
				    +"\teg: weka.clusterers.EM", 
				    "W", 1, 
				    "-W <class name>"));

    if ((m_clusterer != null) &&
	(m_clusterer instanceof OptionHandler)) {
      newVector.addElement(new Option(
				      "",
				      "", 0, "\nOptions specific to clusterer "
				      + m_clusterer.getClass().getName() + ":"));
      Enumeration enu = ((OptionHandler)m_clusterer).listOptions();
      while (enu.hasMoreElements()) {
	newVector.addElement(enu.nextElement());
      }
    }
    return newVector.elements();
  }

  /**
   * Parses a given list of options. Valid options are:<p>
   *
   * -W classname <br>
   * Specify the full class name of the clusterer to evaluate. <p>
   *
   * All option after -- will be passed to the classifier.
   *
   * @param options the list of options as an array of strings
   * @exception Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    
    String cName = Utils.getOption('W', options);
    if (cName.length() == 0) {
      throw new Exception("A clusterer must be specified with"
			  + " the -W option.");
    }
    // Do it first without options, so if an exception is thrown during
    // the option setting, listOptions will contain options for the actual
    // Classifier.
    setClusterer((DensityBasedClusterer)AbstractClusterer.forName(cName, null));
    if (getClusterer() instanceof OptionHandler) {
      ((OptionHandler) getClusterer())
	.setOptions(Utils.partitionOptions(options));
      updateOptions();
    }
  }

  /**
   * Gets the current settings of the Classifier.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String [] getOptions() {

    String [] clustererOptions = new String [0];
    if ((m_clusterer != null) && 
	(m_clusterer instanceof OptionHandler)) {
      clustererOptions = ((OptionHandler)m_clusterer).getOptions();
    }
    
    String [] options = new String [clustererOptions.length + 3];
    int current = 0;

    if (getClusterer() != null) {
      options[current++] = "-W";
      options[current++] = getClusterer().getClass().getName();
    }

    options[current++] = "--";

    System.arraycopy(clustererOptions, 0, options, current, 
		     clustererOptions.length);
    current += clustererOptions.length;
    while (current < options.length) {
      options[current++] = "";
    }
    return options;
  }

  /**
   * Set a list of method names for additional measures to look for
   * in Classifiers. This could contain many measures (of which only a
   * subset may be produceable by the current Classifier) if an experiment
   * is the type that iterates over a set of properties.
   * @param additionalMeasures a list of method names
   */
  public void setAdditionalMeasures(String [] additionalMeasures) {
    // System.err.println("ClassifierSplitEvaluator: setting additional measures");
    m_additionalMeasures = additionalMeasures;
    
    // determine which (if any) of the additional measures this clusterer
    // can produce
    if (m_additionalMeasures != null && m_additionalMeasures.length > 0) {
      m_doesProduce = new boolean [m_additionalMeasures.length];

      if (m_clusterer instanceof AdditionalMeasureProducer) {
	Enumeration en = ((AdditionalMeasureProducer)m_clusterer).
	  enumerateMeasures();
	while (en.hasMoreElements()) {
	  String mname = (String)en.nextElement();
	  for (int j=0;j<m_additionalMeasures.length;j++) {
	    if (mname.compareToIgnoreCase(m_additionalMeasures[j]) == 0) {
	      m_doesProduce[j] = true;
	    }
	  }
	}
      }
    } else {
      m_doesProduce = null;
    }
  }

  /**
   * Returns an enumeration of any additional measure names that might be
   * in the classifier
   * @return an enumeration of the measure names
   */
  public Enumeration enumerateMeasures() {
    Vector newVector = new Vector();
    if (m_clusterer instanceof AdditionalMeasureProducer) {
      Enumeration en = ((AdditionalMeasureProducer)m_clusterer).
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
   * @exception IllegalArgumentException if the named measure is not supported
   */
  public double getMeasure(String additionalMeasureName) {
    if (m_clusterer instanceof AdditionalMeasureProducer) {
      return ((AdditionalMeasureProducer)m_clusterer).
	getMeasure(additionalMeasureName);
    } else {
      throw new IllegalArgumentException("DensityBasedClustererSplitEvaluator: "
					 +"Can't return value for : "+additionalMeasureName
					 +". "+m_clusterer.getClass().getName()+" "
					 +"is not an AdditionalMeasureProducer");
    }
  }

  /**
   * Gets the data types of each of the key columns produced for a single run.
   * The number of key fields must be constant
   * for a given SplitEvaluator.
   *
   * @return an array containing objects of the type of each key column. The 
   * objects should be Strings, or Doubles.
   */
  public Object [] getKeyTypes() {

    Object [] keyTypes = new Object[KEY_SIZE];
    keyTypes[0] = "";
    keyTypes[1] = "";
    keyTypes[2] = "";
    return keyTypes;
  }

  /**
   * Gets the names of each of the key columns produced for a single run.
   * The number of key fields must be constant
   * for a given SplitEvaluator.
   *
   * @return an array containing the name of each key column
   */
  public String [] getKeyNames() {

    String [] keyNames = new String[KEY_SIZE];
    keyNames[0] = "Scheme";
    keyNames[1] = "Scheme_options";
    keyNames[2] = "Scheme_version_ID";
    return keyNames;
  }

  /**
   * Gets the key describing the current SplitEvaluator. For example
   * This may contain the name of the classifier used for classifier
   * predictive evaluation. The number of key fields must be constant
   * for a given SplitEvaluator.
   *
   * @return an array of objects containing the key.
   */
  public Object [] getKey(){

    Object [] key = new Object[KEY_SIZE];
    key[0] = m_clusterer.getClass().getName();
    key[1] = m_clustererOptions;
    key[2] = m_clustererVersion;
    return key;
  }

  /**
   * Gets the data types of each of the result columns produced for a 
   * single run. The number of result fields must be constant
   * for a given SplitEvaluator.
   *
   * @return an array containing objects of the type of each result column. 
   * The objects should be Strings, or Doubles.
   */
  public Object [] getResultTypes() {
    int addm = (m_additionalMeasures != null) 
      ? m_additionalMeasures.length 
      : 0;
    int overall_length = RESULT_SIZE+addm;

    Object [] resultTypes = new Object[overall_length];
    Double doub = new Double(0);
    int current = 0;
    
    // number of training and testing instances
    resultTypes[current++] = doub;
    resultTypes[current++] = doub;
    
    // log liklihood
    resultTypes[current++] = doub;
    // number of clusters
    resultTypes[current++] = doub;

    // timing stats
    resultTypes[current++] = doub;
    resultTypes[current++] = doub;


    //    resultTypes[current++] = "";

    // add any additional measures
    for (int i=0;i<addm;i++) {
      resultTypes[current++] = doub;
    }
    if (current != overall_length) {
      throw new Error("ResultTypes didn't fit RESULT_SIZE");
    }
    return resultTypes;
  }

  /**
   * Gets the names of each of the result columns produced for a single run.
   * The number of result fields must be constant
   * for a given SplitEvaluator.
   *
   * @return an array containing the name of each result column
   */
  public String [] getResultNames() {
    int addm = (m_additionalMeasures != null) 
      ? m_additionalMeasures.length 
      : 0;
    int overall_length = RESULT_SIZE+addm;
   
    String [] resultNames = new String[overall_length];
    int current = 0;
    resultNames[current++] = "Number_of_training_instances";
    resultNames[current++] = "Number_of_testing_instances";

    // Basic performance stats
    resultNames[current++] = "Log_likelihood";
    resultNames[current++] = "Number_of_clusters";

    // Timing stats
    resultNames[current++] = "Time_training";
    resultNames[current++] = "Time_testing";

    // Classifier defined extras
    //    resultNames[current++] = "Summary";
    // add any additional measures
    for (int i=0;i<addm;i++) {
      resultNames[current++] = m_additionalMeasures[i];
    }
    if (current != overall_length) {
      throw new Error("ResultNames didn't fit RESULT_SIZE");
    }
    return resultNames;
  }

  /**
   * Gets the results for the supplied train and test datasets.
   *
   * @param train the training Instances.
   * @param test the testing Instances.
   * @return the results stored in an array. The objects stored in
   * the array may be Strings, Doubles, or null (for the missing value).
   * @exception Exception if a problem occurs while getting the results
   */
  public Object [] getResult(Instances train, Instances test) 
    throws Exception {
    
    if (m_clusterer == null) {
      throw new Exception("No clusterer has been specified");
    }
    int addm = (m_additionalMeasures != null) 
      ? m_additionalMeasures.length 
      : 0;
    int overall_length = RESULT_SIZE+addm;

    if (m_removeClassColumn && train.classIndex() != -1) {
      // remove the class column from the training and testing data
      Remove r = new Remove();
      r.setAttributeIndicesArray(new int [] {train.classIndex()});
      r.setInvertSelection(false);
      r.setInputFormat(train);
      train = Filter.useFilter(train, r);
      
      test = Filter.useFilter(test, r);
    }
    train.setClassIndex(-1);
    test.setClassIndex(-1);
      

    ClusterEvaluation eval = new ClusterEvaluation();

    Object [] result = new Object[overall_length];
    long trainTimeStart = System.currentTimeMillis();
    m_clusterer.buildClusterer(train);
    double numClusters = m_clusterer.numberOfClusters();
    eval.setClusterer(m_clusterer);
    long trainTimeElapsed = System.currentTimeMillis() - trainTimeStart;
    long testTimeStart = System.currentTimeMillis();
    eval.evaluateClusterer(test);
    long testTimeElapsed = System.currentTimeMillis() - testTimeStart;
    //    m_result = eval.toSummaryString();

    // The results stored are all per instance -- can be multiplied by the
    // number of instances to get absolute numbers
    int current = 0;
    result[current++] = new Double(train.numInstances());
    result[current++] = new Double(test.numInstances());

    result[current++] = new Double(eval.getLogLikelihood());
    result[current++] = new Double(numClusters);
    
    // Timing stats
    result[current++] = new Double(trainTimeElapsed / 1000.0);
    result[current++] = new Double(testTimeElapsed / 1000.0);
    
    for (int i=0;i<addm;i++) {
      if (m_doesProduce[i]) {
	try {
	  double dv = ((AdditionalMeasureProducer)m_clusterer).
	    getMeasure(m_additionalMeasures[i]);
	  Double value = new Double(dv);
	  
	  result[current++] = value;
	} catch (Exception ex) {
	  System.err.println(ex);
	}
      } else {
	result[current++] = null;
      }
    }
    
    if (current != overall_length) {
      throw new Error("Results didn't fit RESULT_SIZE");
    }
    return result;
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String removeClassColumnTipText() {
    return "Remove the class column (if set) from the data.";
  }

  /**
   * Set whether the class column should be removed from the data.
   *
   * @param r true if the class column is to be removed.
   */
  public void setRemoveClassColumn(boolean r) {
    m_removeClassColumn = r;
  }

  /**
   * Get whether the class column is to be removed.
   *
   * @return true if the class column is to be removed.
   */
  public boolean getRemoveClassColumn() {
    return m_removeClassColumn;
  }
  
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String clustererTipText() {
    return "The density based clusterer to use.";
  }

  /**
   * Get the value of clusterer
   *
   * @return Value of clusterer.
   */
  public DensityBasedClusterer getClusterer() {
    
    return m_clusterer;
  }
  
  /**
   * Sets the clusterer.
   *
   * @param newClusterer the new clusterer to use.
   */
  public void setClusterer(DensityBasedClusterer newClusterer) {
    
    m_clusterer = newClusterer;
    updateOptions();
  }


  protected void updateOptions() {
    
    if (m_clusterer instanceof OptionHandler) {
      m_clustererOptions = Utils.joinOptions(((OptionHandler)m_clusterer)
					     .getOptions());
    } else {
      m_clustererOptions = "";
    }
    if (m_clusterer instanceof Serializable) {
      ObjectStreamClass obs = ObjectStreamClass.lookup(m_clusterer
						       .getClass());
      m_clustererVersion = "" + obs.getSerialVersionUID();
    } else {
      m_clustererVersion = "";
    }
  }

  /**
   * Set the Clusterer to use, given it's class name. A new clusterer will be
   * instantiated.
   *
   * @param newClustererName the clusterer class name.
   * @exception Exception if the class name is invalid.
   */
  public void setClustererName(String newClustererName) throws Exception {

    try {
      setClusterer((DensityBasedClusterer)Class.forName(newClustererName)
		    .newInstance());
    } catch (Exception ex) {
      throw new Exception("Can't find Clusterer with class name: "
			  + newClustererName);
    }
  }

  /**
   * Gets the raw output from the classifier
   * @return the raw output from the classifier
   */
  public String getRawResultOutput() {
    StringBuffer result = new StringBuffer();

    if (m_clusterer == null) {
      return "<null> clusterer";
    }
    result.append(toString());
    result.append("Clustering model: \n"+m_clusterer.toString()+'\n');

    // append the performance statistics
    if (m_result != null) {
      //      result.append(m_result);
      
      if (m_doesProduce != null) {
	for (int i=0;i<m_doesProduce.length;i++) {
	  if (m_doesProduce[i]) {
	    try {
	      double dv = ((AdditionalMeasureProducer)m_clusterer).
		getMeasure(m_additionalMeasures[i]);
	      Double value = new Double(dv);
	      
	      result.append(m_additionalMeasures[i]+" : "+value+'\n');
	    } catch (Exception ex) {
	      System.err.println(ex);
	    }
	  } 
	}
      }
    }
    return result.toString();
  }

  /**
   * Returns a text description of the split evaluator.
   *
   * @return a text description of the split evaluator.
   */
  public String toString() {

    String result = "DensityBasedClustererSplitEvaluator: ";
    if (m_clusterer == null) {
      return result + "<null> clusterer";
    }
    return result + m_clusterer.getClass().getName() + " " 
      + m_clustererOptions + "(version " + m_clustererVersion + ")";
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5563 $");
  }
}
