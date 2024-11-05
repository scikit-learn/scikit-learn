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
 *    ClusterMembership.java
 *    Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.filters.unsupervised.attribute;

import weka.clusterers.DensityBasedClusterer;
import weka.clusterers.AbstractDensityBasedClusterer;
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
import weka.core.Utils;
import weka.core.labelranking.PreferenceDenseInstance;
import weka.filters.Filter;
import weka.filters.UnsupervisedFilter;

import java.util.Enumeration;
import java.util.Vector;

/** 
 <!-- globalinfo-start -->
 * A filter that uses a density-based clusterer to generate cluster membership values; filtered instances are composed of these values plus the class attribute (if set in the input data). If a (nominal) class attribute is set, the clusterer is run separately for each class. The class attribute (if set) and any user-specified attributes are ignored during the clustering operation
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -W &lt;clusterer name&gt;
 *  Full name of clusterer to use. eg:
 *   weka.clusterers.EM
 *  Additional options after the '--'.
 *  (default: weka.clusterers.EM)</pre>
 * 
 * <pre> -I &lt;att1,att2-att4,...&gt;
 *  The range of attributes the clusterer should ignore.
 *  (the class attribute is automatically ignored)</pre>
 * 
 <!-- options-end -->
 *
 * Options after the -- are passed on to the clusterer.
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @author Eibe Frank
 * @version $Revision: 5987 $
 */
public class ClusterMembership 
  extends Filter 
  implements UnsupervisedFilter, OptionHandler {
  
  /** for serialization */
  static final long serialVersionUID = 6675702504667714026L;

  /** The clusterer */
  protected DensityBasedClusterer m_clusterer = new weka.clusterers.EM();

  /** Array for storing the clusterers */
  protected DensityBasedClusterer[] m_clusterers;

  /** Range of attributes to ignore */
  protected Range m_ignoreAttributesRange;

  /** Filter for removing attributes */
  protected Filter m_removeAttributes;

  /** The prior probability for each class */
  protected double[] m_priors;

  /** 
   * Returns the Capabilities of this filter.
   *
   * @return            the capabilities of this object
   * @see               Capabilities
   */
  public Capabilities getCapabilities() {
    Capabilities result = m_clusterer.getCapabilities();
    
    result.setMinimumNumberInstances(0);
    
    return result;
  }

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
   * tests the data whether the filter can actually handle it
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
    m_priors = null;

    return false;
  }

  /**
   * filters all attributes that should be ignored
   * 
   * @param data	the data to filter
   * @return		the filtered data
   * @throws Exception	if filtering fails
   */
  protected Instances removeIgnored(Instances data) throws Exception {
    Instances result = data;
    
    if (m_ignoreAttributesRange != null || data.classIndex() >= 0) {
      result = new Instances(data);
      m_removeAttributes = new Remove();
      String rangeString = "";
      if (m_ignoreAttributesRange != null) {
	rangeString += m_ignoreAttributesRange.getRanges();
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

    if (getInputFormat() == null) {
      throw new IllegalStateException("No input instance format defined");
    }

    if (outputFormatPeek() == null) {
      Instances toFilter = getInputFormat();
      Instances[] toFilterIgnoringAttributes;

      // Make subsets if class is nominal
      if ((toFilter.classIndex() >= 0) && (toFilter.classAttribute().isNominal() || toFilter.classAttribute().isRanking())) {
	toFilterIgnoringAttributes = new Instances[toFilter.numClasses()];
	for (int i = 0; i < toFilter.numClasses(); i++) {
	  toFilterIgnoringAttributes[i] = new Instances(toFilter, toFilter.numInstances());
	}
	for (int i = 0; i < toFilter.numInstances(); i++) {
	  toFilterIgnoringAttributes[(int)toFilter.instance(i).classValue()].add(toFilter.instance(i));
	}
	m_priors = new double[toFilter.numClasses()];
	for (int i = 0; i < toFilter.numClasses(); i++) {
	  toFilterIgnoringAttributes[i].compactify();
	  m_priors[i] = toFilterIgnoringAttributes[i].sumOfWeights();
	}
	Utils.normalize(m_priors);
      } else {
	toFilterIgnoringAttributes = new Instances[1];
	toFilterIgnoringAttributes[0] = toFilter;
	m_priors = new double[1];
	m_priors[0] = 1;
      }

      // filter out attributes if necessary
      for (int i = 0; i < toFilterIgnoringAttributes.length; i++)
	toFilterIgnoringAttributes[i] = removeIgnored(toFilterIgnoringAttributes[i]);

      // build the clusterers
      if ((toFilter.classIndex() <= 0) || !toFilter.classAttribute().isNominal() && !toFilter.classAttribute().isRanking()) {
	m_clusterers = AbstractDensityBasedClusterer.makeCopies(m_clusterer, 1);
	m_clusterers[0].buildClusterer(toFilterIgnoringAttributes[0]);
      } else {
	m_clusterers = AbstractDensityBasedClusterer.makeCopies(m_clusterer, toFilter.numClasses());
	for (int i = 0; i < m_clusterers.length; i++) {
	  if (toFilterIgnoringAttributes[i].numInstances() == 0) {
	    m_clusterers[i] = null;
	  } else {
	    m_clusterers[i].buildClusterer(toFilterIgnoringAttributes[i]);
	  }
	}
      }
      
      // create output dataset
      FastVector attInfo = new FastVector();
      for (int j = 0; j < m_clusterers.length; j++) {
	if (m_clusterers[j] != null) {
	  for (int i = 0; i < m_clusterers[j].numberOfClusters(); i++) {
	    attInfo.addElement(new Attribute("pCluster_" + j + "_" + i));
	  }
	}
      }
      if (toFilter.classIndex() >= 0) {
	attInfo.addElement(toFilter.classAttribute().copy());
      }
      attInfo.trimToSize();
      Instances filtered = new Instances(toFilter.relationName()+"_clusterMembership",
					 attInfo, 0);
      if (toFilter.classIndex() >= 0) {
	filtered.setClassIndex(filtered.numAttributes() - 1);
      }
      setOutputFormat(filtered);

      // build new dataset
      for (int i = 0; i < toFilter.numInstances(); i++) {
	convertInstance(toFilter.instance(i));
      }
    }
    flushInput();

    m_NewBatch = true;
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

    if (getInputFormat() == null) {
      throw new IllegalStateException("No input instance format defined");
    }
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
   * Converts logs back to density values.
   * 
   * @param j the index of the clusterer
   * @param in the instance to convert the logs back
   * @return the densities
   * @throws Exception if something goes wrong
   */
  protected double[] logs2densities(int j, Instance in) throws Exception {

    double[] logs = m_clusterers[j].logJointDensitiesForInstance(in);

    for (int i = 0; i < logs.length; i++) {
      logs[i] += Math.log(m_priors[j]);
    }
    return logs;
  }

  /**
   * Convert a single instance over. The converted instance is added to 
   * the end of the output queue.
   *
   * @param instance the instance to convert
   * @throws Exception if something goes wrong
   */
  protected void convertInstance(Instance instance) throws Exception {
    
    // set up values
    double [] instanceVals = new double[outputFormatPeek().numAttributes()];
    double [] tempvals;
    if (instance.classIndex() >= 0) {
      tempvals = new double[outputFormatPeek().numAttributes() - 1];
    } else {
      tempvals = new double[outputFormatPeek().numAttributes()];
    }
    int pos = 0;
    for (int j = 0; j < m_clusterers.length; j++) {
      if (m_clusterers[j] != null) {
	double [] probs;
	if (m_removeAttributes != null) {
	  m_removeAttributes.input(instance);
	  probs = logs2densities(j, m_removeAttributes.output());
	} else {
	  probs = logs2densities(j, instance);
	}
	System.arraycopy(probs, 0, tempvals, pos, probs.length);
	pos += probs.length;
      }
    }
    tempvals = Utils.logs2probs(tempvals);
    System.arraycopy(tempvals, 0, instanceVals, 0, tempvals.length);
    if (instance.classIndex() >= 0) {
      instanceVals[instanceVals.length - 1] = instance.classValue();
    }
    //RANKING BEGIN
    if(instance instanceof PreferenceDenseInstance){
    	PreferenceDenseInstance pdi = (PreferenceDenseInstance) instance;
    	push(new PreferenceDenseInstance(instance.weight(), instanceVals, pdi.getHashMap()));
    }
    else
    	push(new DenseInstance(instance.weight(), instanceVals));
    //RANKING END
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    
    Vector newVector = new Vector(2);
    
    newVector.
      addElement(new Option("\tFull name of clusterer to use. eg:\n"
	                    + "\t\tweka.clusterers.EM\n"
			    + "\tAdditional options after the '--'.\n"
			    + "\t(default: weka.clusterers.EM)",
			    "W", 1, "-W <clusterer name>"));

    newVector.
      addElement(new Option("\tThe range of attributes the clusterer should ignore."
			    +"\n\t(the class attribute is automatically ignored)",
			    "I", 1,"-I <att1,att2-att4,...>"));

    return newVector.elements();
  }

  /**
   * Parses a given list of options. <p/>
   * 
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -W &lt;clusterer name&gt;
   *  Full name of clusterer to use. eg:
   *   weka.clusterers.EM
   *  Additional options after the '--'.
   *  (default: weka.clusterers.EM)</pre>
   * 
   * <pre> -I &lt;att1,att2-att4,...&gt;
   *  The range of attributes the clusterer should ignore.
   *  (the class attribute is automatically ignored)</pre>
   * 
   <!-- options-end -->
   *
   * Options after the -- are passed on to the clusterer.
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {

    String clustererString = Utils.getOption('W', options);
    if (clustererString.length() == 0)
      clustererString = weka.clusterers.EM.class.getName();
    setDensityBasedClusterer((DensityBasedClusterer)Utils.
			     forName(DensityBasedClusterer.class, clustererString,
				     Utils.partitionOptions(options)));

    setIgnoredAttributeIndices(Utils.getOption('I', options));
    Utils.checkForRemainingOptions(options);
  }

  /**
   * Gets the current settings of the filter.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String [] getOptions() {

    String [] clustererOptions = new String [0];
    if ((m_clusterer != null) &&
	(m_clusterer instanceof OptionHandler)) {
      clustererOptions = ((OptionHandler)m_clusterer).getOptions();
    }
    String [] options = new String [clustererOptions.length + 5];
    int current = 0;

    if (!getIgnoredAttributeIndices().equals("")) {
      options[current++] = "-I";
      options[current++] = getIgnoredAttributeIndices();
    }
    
    if (m_clusterer != null) {
      options[current++] = "-W"; 
      options[current++] = getDensityBasedClusterer().getClass().getName();
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
   * Returns a string describing this filter
   *
   * @return a description of the filter suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {

    return "A filter that uses a density-based clusterer to generate cluster "
      + "membership values; filtered instances are composed of these values "
      + "plus the class attribute (if set in the input data). If a (nominal) "
      + "class attribute is set, the clusterer is run separately for each "
      + "class. The class attribute (if set) and any user-specified "
      + "attributes are ignored during the clustering operation";
  }
  
  /**
   * Returns a description of this option suitable for display
   * as a tip text in the gui.
   *
   * @return description of this option
   */
  public String densityBasedClustererTipText() {
    return "The clusterer that will generate membership values for the instances.";
  }

  /**
   * Set the clusterer for use in filtering
   *
   * @param newClusterer the clusterer to use
   */
  public void setDensityBasedClusterer(DensityBasedClusterer newClusterer) {
    m_clusterer = newClusterer;
  }

  /**
   * Get the clusterer used by this filter
   *
   * @return the clusterer used
   */
  public DensityBasedClusterer getDensityBasedClusterer() {
    return m_clusterer;
  }

  /**
   * Returns the tip text for this property
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

    if (m_ignoreAttributesRange == null) {
      return "";
    } else {
      return m_ignoreAttributesRange.getRanges();
    }
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
      m_ignoreAttributesRange = null;
    } else {
      m_ignoreAttributesRange = new Range();
      m_ignoreAttributesRange.setRanges(rangeList);
    }
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
  public static void main(String [] argv) {
    runFilter(new ClusterMembership(), argv);
  }
}
