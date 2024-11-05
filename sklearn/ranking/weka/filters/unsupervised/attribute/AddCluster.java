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
 *    AddCluster.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.filters.unsupervised.attribute;

import weka.clusterers.AbstractClusterer;
import weka.clusterers.Clusterer;
import weka.core.Attribute;
import weka.core.Capabilities;
import weka.core.FastVector;
import weka.core.Instance; 
import weka.core.DenseInstance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.Range;
import weka.core.RevisionUtils;
import weka.core.SparseInstance;
import weka.core.Utils;
import weka.core.WekaException;
import weka.core.labelranking.PreferenceDenseInstance;
import weka.filters.Filter;
import weka.filters.UnsupervisedFilter;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.ObjectInputStream;
import java.util.Enumeration;
import java.util.Vector;

/** 
 <!-- globalinfo-start -->
 * A filter that adds a new nominal attribute representing the cluster assigned to each instance by the specified clustering algorithm.<br/>
 * Either the clustering algorithm gets built with the first batch of data or one specifies are serialized clusterer model file to use instead.
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -W &lt;clusterer specification&gt;
 *  Full class name of clusterer to use, followed
 *  by scheme options. eg:
 *   "weka.clusterers.SimpleKMeans -N 3"
 *  (default: weka.clusterers.SimpleKMeans)</pre>
 * 
 * <pre> -serialized &lt;file&gt;
 *  Instead of building a clusterer on the data, one can also provide
 *  a serialized model and use that for adding the clusters.</pre>
 * 
 * <pre> -I &lt;att1,att2-att4,...&gt;
 *  The range of attributes the clusterer should ignore.
 * </pre>
 * 
 <!-- options-end -->
 *
 * @author Richard Kirkby (rkirkby@cs.waikato.ac.nz)
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5987 $
 */
public class AddCluster 
  extends Filter 
  implements UnsupervisedFilter, OptionHandler {
  
  /** for serialization. */
  static final long serialVersionUID = 7414280611943807337L;

  /** The clusterer used to do the cleansing. */
  protected Clusterer m_Clusterer = new weka.clusterers.SimpleKMeans();

  /** The file from which to load a serialized clusterer. */
  protected File m_SerializedClustererFile = new File(System.getProperty("user.dir"));
  
  /** The actual clusterer used to do the clustering. */
  protected Clusterer m_ActualClusterer = null;

  /** Range of attributes to ignore. */
  protected Range m_IgnoreAttributesRange = null;

  /** Filter for removing attributes. */
  protected Filter m_removeAttributes = new Remove();

  /** 
   * Returns the Capabilities of this filter, makes sure that the class is
   * never set (for the clusterer).
   *
   * @param data	the data to use for customization
   * @return            the capabilities of this object, based on the data
   * @see               #getCapabilities()
   */
  public Capabilities getCapabilities(Instances data) {
    Instances	newData;
    
    newData = new Instances(data, 0);
    newData.setClassIndex(-1);
    
    return super.getCapabilities(newData);
  }

  /** 
   * Returns the Capabilities of this filter.
   *
   * @return            the capabilities of this object
   * @see               Capabilities
   */
  public Capabilities getCapabilities() {
    Capabilities result = m_Clusterer.getCapabilities();
    
    result.setMinimumNumberInstances(0);
    
    return result;
  }
  
  /**
   * tests the data whether the filter can actually handle it.
   * 
   * @param instanceInfo	the data to test
   * @throws Exception		if the test fails
   */
  protected void testInputFormat(Instances instanceInfo) throws Exception {
    getCapabilities(instanceInfo).testWithFail(removeIgnored(instanceInfo));
  }

  /**
   * Sets the format of the input instances.
   *
   * @param instanceInfo an Instances object containing the input instance
   * structure (any instances contained in the object are ignored - only the
   * structure is required).
   * @return true if the outputFormat may be collected immediately
   * @throws Exception if the inputFormat can't be set successfully 
   */ 
  public boolean setInputFormat(Instances instanceInfo) throws Exception {
    super.setInputFormat(instanceInfo);

    m_removeAttributes = null;

    return false;
  }

  /**
   * filters all attributes that should be ignored.
   * 
   * @param data	the data to filter
   * @return		the filtered data
   * @throws Exception	if filtering fails
   */
  protected Instances removeIgnored(Instances data) throws Exception {
    Instances result = data;
    
    if (m_IgnoreAttributesRange != null || data.classIndex() >= 0) {
      m_removeAttributes = new Remove();
      String rangeString = "";
      if (m_IgnoreAttributesRange != null) {
	rangeString += m_IgnoreAttributesRange.getRanges();
      }
      if (data.classIndex() >= 0) {
	if (rangeString.length() > 0) {
	  rangeString += "," + (data.classIndex() + 1);
	} else {
	  rangeString = "" + (data.classIndex() + 1);
	}
      }
      ((Remove) m_removeAttributes).setAttributeIndices(rangeString);
      ((Remove) m_removeAttributes).setInvertSelection(false);
      m_removeAttributes.setInputFormat(data);
      result = Filter.useFilter(data, m_removeAttributes);
    }
    
    return result;
  }
  
  /**
   * Signify that this batch of input to the filter is finished.
   *
   * @return true if there are instances pending output
   * @throws IllegalStateException if no input structure has been defined 
   */  
  public boolean batchFinished() throws Exception {
    if (getInputFormat() == null)
      throw new IllegalStateException("No input instance format defined");

    Instances toFilter = getInputFormat();
    
    if (!isFirstBatchDone()) {
      // filter out attributes if necessary
      Instances toFilterIgnoringAttributes = removeIgnored(toFilter);

      // serialized model or build clusterer from scratch?
      File file = getSerializedClustererFile();
      if (!file.isDirectory()) {
	ObjectInputStream ois = new ObjectInputStream(new FileInputStream(file));
	m_ActualClusterer = (Clusterer) ois.readObject();
	Instances header = null;
	// let's see whether there's an Instances header stored as well
	try {
	  header = (Instances) ois.readObject();
	}
	catch (Exception e) {
	  // ignored
	}
	ois.close();
	// same dataset format?
	if ((header != null) && (!header.equalHeaders(toFilterIgnoringAttributes)))
	  throw new WekaException(
	      "Training header of clusterer and filter dataset don't match:\n"
	      + header.equalHeadersMsg(toFilterIgnoringAttributes));
      }
      else {
	m_ActualClusterer = AbstractClusterer.makeCopy(m_Clusterer);
	m_ActualClusterer.buildClusterer(toFilterIgnoringAttributes);
      }

      // create output dataset with new attribute
      Instances filtered = new Instances(toFilter, 0); 
      FastVector nominal_values = new FastVector(m_ActualClusterer.numberOfClusters());
      for (int i = 0; i < m_ActualClusterer.numberOfClusters(); i++) {
	nominal_values.addElement("cluster" + (i+1)); 
      }
      filtered.insertAttributeAt(new Attribute("cluster", nominal_values),
	  filtered.numAttributes());

      setOutputFormat(filtered);
    }

    // build new dataset
    for (int i=0; i<toFilter.numInstances(); i++) {
      convertInstance(toFilter.instance(i));
    }
    
    flushInput();
    m_NewBatch = true;
    m_FirstBatchDone = true;

    return (numPendingOutput() != 0);
  }

  /**
   * Input an instance for filtering. Ordinarily the instance is processed
   * and made available for output immediately. Some filters require all
   * instances be read before producing output.
   *
   * @param instance the input instance
   * @return true if the filtered instance may now be
   * collected with output().
   * @throws IllegalStateException if no input format has been defined.
   */
  public boolean input(Instance instance) throws Exception {
    if (getInputFormat() == null)
      throw new IllegalStateException("No input instance format defined");

    if (m_NewBatch) {
      resetQueue();
      m_NewBatch = false;
    }
    
    if (outputFormatPeek() != null) {
      convertInstance(instance);
      return true;
    }

    bufferInput(instance);
    return false;
  }

  /**
   * Convert a single instance over. The converted instance is added to 
   * the end of the output queue.
   *
   * @param instance the instance to convert
   * @throws Exception if something goes wrong
   */
  protected void convertInstance(Instance instance) throws Exception {
    Instance original, processed;
    original = instance;

    // copy values
    double[] instanceVals = new double[instance.numAttributes()+1];
    for(int j = 0; j < instance.numAttributes(); j++) {
      instanceVals[j] = original.value(j);
    }
    Instance filteredI = null;
    if (m_removeAttributes != null) {
      m_removeAttributes.input(instance);
      filteredI = m_removeAttributes.output();
    } else {
      filteredI = instance;
    }

    // add cluster to end
    try {
      instanceVals[instance.numAttributes()] = m_ActualClusterer.clusterInstance(filteredI);
    }
    catch (Exception e) {
      // clusterer couldn't cluster instance -> missing
      instanceVals[instance.numAttributes()] = Utils.missingValue();
    }

    // create new instance
    if (original instanceof SparseInstance) {
      processed = new SparseInstance(original.weight(), instanceVals);
    } 
    //RANKING BEGIN
    if(original instanceof PreferenceDenseInstance){
    	PreferenceDenseInstance pdi = (PreferenceDenseInstance) original;
    	processed = new PreferenceDenseInstance(original.weight(), instanceVals, pdi.getHashMap());
    }
    //RANKING END
    else {
      processed = new DenseInstance(original.weight(), instanceVals);
    }

    processed.setDataset(instance.dataset());
    copyValues(processed, false, instance.dataset(), getOutputFormat());
    processed.setDataset(getOutputFormat());
      
    push(processed);
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector result = new Vector();
    
    result.addElement(new Option(
	"\tFull class name of clusterer to use, followed\n"
	+ "\tby scheme options. eg:\n"
	+ "\t\t\"weka.clusterers.SimpleKMeans -N 3\"\n"
	+ "\t(default: weka.clusterers.SimpleKMeans)",
	"W", 1, "-W <clusterer specification>"));

    result.addElement(new Option(
	"\tInstead of building a clusterer on the data, one can also provide\n"
	+ "\ta serialized model and use that for adding the clusters.",
	"serialized", 1, "-serialized <file>"));
    
    result.addElement(new Option(
	"\tThe range of attributes the clusterer should ignore.\n",
	"I", 1,"-I <att1,att2-att4,...>"));

    return result.elements();
  }


  /**
   * Parses a given list of options. <p/>
   * 
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -W &lt;clusterer specification&gt;
   *  Full class name of clusterer to use, followed
   *  by scheme options. eg:
   *   "weka.clusterers.SimpleKMeans -N 3"
   *  (default: weka.clusterers.SimpleKMeans)</pre>
   * 
   * <pre> -serialized &lt;file&gt;
   *  Instead of building a clusterer on the data, one can also provide
   *  a serialized model and use that for adding the clusters.</pre>
   * 
   * <pre> -I &lt;att1,att2-att4,...&gt;
   *  The range of attributes the clusterer should ignore.
   * </pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String	tmpStr;
    String[] 	tmpOptions;
    File	file;
    boolean 	serializedModel;
    
    serializedModel = false;
    tmpStr = Utils.getOption("serialized", options);
    if (tmpStr.length() != 0) {
      file = new File(tmpStr);
      if (!file.exists())
	throw new FileNotFoundException(
	    "File '" + file.getAbsolutePath() + "' not found!");
      if (file.isDirectory())
	throw new FileNotFoundException(
	    "'" + file.getAbsolutePath() + "' points to a directory not a file!");
      setSerializedClustererFile(file);
      serializedModel = true;
    }
    else {
      setSerializedClustererFile(null);
    }

    if (!serializedModel) {
      tmpStr = Utils.getOption('W', options);
      if (tmpStr.length() == 0)
	tmpStr = weka.clusterers.SimpleKMeans.class.getName();
      tmpOptions = Utils.splitOptions(tmpStr);
      if (tmpOptions.length == 0) {
	throw new Exception("Invalid clusterer specification string");
      }
      tmpStr = tmpOptions[0];
      tmpOptions[0] = "";
      setClusterer(AbstractClusterer.forName(tmpStr, tmpOptions));
    }
        
    setIgnoredAttributeIndices(Utils.getOption('I', options));

    Utils.checkForRemainingOptions(options);
  }

  /**
   * Gets the current settings of the filter.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    Vector<String>	result;
    File		file;

    result = new Vector<String>();

    file = getSerializedClustererFile();
    if ((file != null) && (!file.isDirectory())) {
      result.add("-serialized");
      result.add(file.getAbsolutePath());
    }
    else {
      result.add("-W");
      result.add(getClustererSpec());
    }
    
    if (!getIgnoredAttributeIndices().equals("")) {
      result.add("-I");
      result.add(getIgnoredAttributeIndices());
    }
    
    return result.toArray(new String[result.size()]);
  }

  /**
   * Returns a string describing this filter.
   *
   * @return a description of the filter suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "A filter that adds a new nominal attribute representing the cluster "
      + "assigned to each instance by the specified clustering algorithm.\n"
      + "Either the clustering algorithm gets built with the first batch of "
      + "data or one specifies are serialized clusterer model file to use "
      + "instead.";
  }

  /**
   * Returns the tip text for this property.
   *
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String clustererTipText() {
    return "The clusterer to assign clusters with.";
  }

  /**
   * Sets the clusterer to assign clusters with.
   *
   * @param clusterer The clusterer to be used (with its options set).
   */
  public void setClusterer(Clusterer clusterer) {
    m_Clusterer = clusterer;
  }
  
  /**
   * Gets the clusterer used by the filter.
   *
   * @return The clusterer being used.
   */
  public Clusterer getClusterer() {
    return m_Clusterer;
  }

  /**
   * Gets the clusterer specification string, which contains the class name of
   * the clusterer and any options to the clusterer.
   *
   * @return the clusterer string.
   */
  protected String getClustererSpec() {
    Clusterer c = getClusterer();
    if (c instanceof OptionHandler) {
      return c.getClass().getName() + " "
	+ Utils.joinOptions(((OptionHandler)c).getOptions());
    }
    return c.getClass().getName();
  }

  /**
   * Returns the tip text for this property.
   *
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String ignoredAttributeIndicesTipText() {
    return "The range of attributes to be ignored by the clusterer. eg: first-3,5,9-last";
  }

  /**
   * Gets ranges of attributes to be ignored.
   *
   * @return a string containing a comma-separated list of ranges
   */
  public String getIgnoredAttributeIndices() {
    if (m_IgnoreAttributesRange == null)
      return "";
    else
      return m_IgnoreAttributesRange.getRanges();
  }

  /**
   * Sets the ranges of attributes to be ignored. If provided string
   * is null, no attributes will be ignored.
   *
   * @param rangeList a string representing the list of attributes. 
   * eg: first-3,5,6-last
   * @throws IllegalArgumentException if an invalid range list is supplied 
   */
  public void setIgnoredAttributeIndices(String rangeList) {
    if ((rangeList == null) || (rangeList.length() == 0)) {
      m_IgnoreAttributesRange = null;
    } else {
      m_IgnoreAttributesRange = new Range();
      m_IgnoreAttributesRange.setRanges(rangeList);
    }
  }

  /**
   * Gets the file pointing to a serialized, built clusterer. If it is
   * null or pointing to a directory it will not be used.
   * 
   * @return		the file the serialized, built clusterer is located in
   */
  public File getSerializedClustererFile() {
    return m_SerializedClustererFile;
  }

  /**
   * Sets the file pointing to a serialized, built clusterer. If the
   * argument is null, doesn't exist or pointing to a directory, then the 
   * value is ignored.
   * 
   * @param value	the file pointing to the serialized, built clusterer
   */
  public void setSerializedClustererFile(File value) {
    if ((value == null) || (!value.exists()))
      value = new File(System.getProperty("user.dir"));

    m_SerializedClustererFile = value;
  }
  
  /**
   * Returns the tip text for this property.
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String serializedClustererFileTipText() {
    return "A file containing the serialized model of a built clusterer.";
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5987 $");
  }

  /**
   * Main method for testing this class.
   *
   * @param argv should contain arguments to the filter: use -h for help
   */
  public static void main(String[] argv) {
    runFilter(new AddCluster(), argv);
  }
}
