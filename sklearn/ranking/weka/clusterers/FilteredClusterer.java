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
 * FilteredClusterer.java
 * Copyright (C) 2006 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.clusterers;

import weka.core.Capabilities;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.filters.Filter;
import weka.filters.SupervisedFilter;

import java.util.Enumeration;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * Class for running an arbitrary clusterer on data that has been passed through an arbitrary filter. Like the clusterer, the structure of the filter is based exclusively on the training data and test instances will be processed by the filter without changing their structure.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -F &lt;filter specification&gt;
 *  Full class name of filter to use, followed
 *  by filter options.
 *  eg: "weka.filters.unsupervised.attribute.Remove -V -R 1,2"
 * (default: weka.filters.AllFilter)</pre>
 * 
 * <pre> -W
 *  Full name of base clusterer.
 *  (default: weka.clusterers.SimpleKMeans)</pre>
 * 
 * <pre> 
 * Options specific to clusterer weka.clusterers.SimpleKMeans:
 * </pre>
 * 
 * <pre> -N &lt;num&gt;
 *  number of clusters.
 *  (default 2).</pre>
 * 
 * <pre> -V
 *  Display std. deviations for centroids.
 * </pre>
 * 
 * <pre> -M
 *  Replace missing values with mean/mode.
 * </pre>
 * 
 * <pre> -S &lt;num&gt;
 *  Random number seed.
 *  (default 10)</pre>
 * 
 <!-- options-end -->
 *
 * Based on code from the FilteredClassifier by Len Trigg.
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5488 $
 * @see weka.classifiers.meta.FilteredClassifier
 */
public class FilteredClusterer
  extends SingleClustererEnhancer {

  /** for serialization. */
  private static final long serialVersionUID = 1420005943163412943L;
  
  /** The filter. */
  protected Filter m_Filter;

  /** The instance structure of the filtered instances. */
  protected Instances m_FilteredInstances;

  /**
   * Default constructor.
   */
  public FilteredClusterer() {
    m_Clusterer = new SimpleKMeans();
    m_Filter    = new weka.filters.AllFilter();
  }

  /**
   * Returns a string describing this clusterer.
   * 
   * @return 		a description of the clusterer suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return
        "Class for running an arbitrary clusterer on data that has been passed "
      + "through an arbitrary filter. Like the clusterer, the structure of the filter "
      + "is based exclusively on the training data and test instances will be processed "
      + "by the filter without changing their structure.";
  }

  /**
   * String describing default filter.
   * 
   * @return 		the default filter classname
   */
  protected String defaultFilterString() {
    return weka.filters.AllFilter.class.getName();
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return 		an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector result = new Vector();

    result.addElement(new Option(
	"\tFull class name of filter to use, followed\n"
	+ "\tby filter options.\n"
	+ "\teg: \"weka.filters.unsupervised.attribute.Remove -V -R 1,2\"\n"
	+ "(default: " + defaultFilterString() + ")",
	"F", 1, "-F <filter specification>"));

    Enumeration enm = super.listOptions();
    while (enm.hasMoreElements())
      result.addElement(enm.nextElement());

    return result.elements();
  }

  /**
   * Parses a given list of options. <p/>
   *
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -F &lt;filter specification&gt;
   *  Full class name of filter to use, followed
   *  by filter options.
   *  eg: "weka.filters.unsupervised.attribute.Remove -V -R 1,2"
   * (default: weka.filters.AllFilter)</pre>
   * 
   * <pre> -W
   *  Full name of base clusterer.
   *  (default: weka.clusterers.SimpleKMeans)</pre>
   * 
   * <pre> 
   * Options specific to clusterer weka.clusterers.SimpleKMeans:
   * </pre>
   * 
   * <pre> -N &lt;num&gt;
   *  number of clusters.
   *  (default 2).</pre>
   * 
   * <pre> -V
   *  Display std. deviations for centroids.
   * </pre>
   * 
   * <pre> -M
   *  Replace missing values with mean/mode.
   * </pre>
   * 
   * <pre> -S &lt;num&gt;
   *  Random number seed.
   *  (default 10)</pre>
   * 
   <!-- options-end -->
   *
   * @param options 	the list of options as an array of strings
   * @throws Exception 	if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String	tmpStr;
    String[]	tmpOptions;
    
    tmpStr = Utils.getOption('F', options);
    if (tmpStr.length() > 0) {
      tmpOptions = Utils.splitOptions(tmpStr);
      if (tmpOptions.length == 0)
	throw new IllegalArgumentException("Invalid filter specification string");
      tmpStr = tmpOptions[0];
      tmpOptions[0] = "";
      setFilter((Filter) Utils.forName(Filter.class, tmpStr, tmpOptions));
    } 
    else {
      setFilter(new weka.filters.AllFilter());
    }
    
    super.setOptions(options);
  }

  /**
   * Gets the current settings of the clusterer.
   *
   * @return 		an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    Vector	result;
    String[]	options;
    int		i;
    
    result = new Vector();
    
    result.add("-F");
    result.add(getFilterSpec());
    
    options = super.getOptions();
    for (i = 0; i < options.length; i++)
      result.add(options[i]);

    return (String[]) result.toArray(new String[result.size()]);
  }
  
  /**
   * Returns the tip text for this property.
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String filterTipText() {
    return "The filter to be used.";
  }

  /**
   * Sets the filter.
   *
   * @param filter 	the filter with all options set.
   */
  public void setFilter(Filter filter) {
    m_Filter = filter;
    
    if (m_Filter instanceof SupervisedFilter)
      System.out.println(
	  "WARNING: you are using a supervised filter, which will leak "
	  + "information about the class attribute!");
  }

  /**
   * Gets the filter used.
   *
   * @return 		the filter
   */
  public Filter getFilter() {
    return m_Filter;
  }
  
  /**
   * Gets the filter specification string, which contains the class name of
   * the filter and any options to the filter.
   *
   * @return 		the filter string.
   */
  protected String getFilterSpec() {
    String	result;
    Filter 	filter;
    
    filter = getFilter();
    result = filter.getClass().getName();
    
    if (filter instanceof OptionHandler)
      result += " " + Utils.joinOptions(((OptionHandler) filter).getOptions());
    
    return result;
  }

  /**
   * Returns default capabilities of the clusterer.
   *
   * @return		the capabilities of this clusterer
   */
  public Capabilities getCapabilities() {
    Capabilities	result;
    
    if (getFilter() == null) {
      result = super.getCapabilities();
      result.disableAll();
      result.enable(Capability.NO_CLASS);
    } else {
      result = getFilter().getCapabilities();
    }
    
    // set dependencies
    for (Capability cap: Capability.values())
      result.enableDependency(cap);
    
    return result;
  }

  /**
   * Build the clusterer on the filtered data.
   *
   * @param data 	the training data
   * @throws Exception 	if the clusterer could not be built successfully
   */
  public void buildClusterer(Instances data) throws Exception {
    if (m_Clusterer == null)
      throw new Exception("No base clusterer has been set!");

    // remove instances with missing class
    if (data.classIndex() > -1) {
      data = new Instances(data);
      data.deleteWithMissingClass();
    }
    
    m_Filter.setInputFormat(data);  // filter capabilities are checked here
    data = Filter.useFilter(data, m_Filter);

    // can clusterer handle the data?
    getClusterer().getCapabilities().testWithFail(data);

    m_FilteredInstances = data.stringFreeStructure();
    m_Clusterer.buildClusterer(data);
  }

  /**
   * Classifies a given instance after filtering.
   *
   * @param instance 	the instance to be classified
   * @return 		the class distribution for the given instance
   * @throws Exception 	if instance could not be classified
   * 			successfully
   */
  public double[] distributionForInstance(Instance instance)
    throws Exception {

    if (m_Filter.numPendingOutput() > 0)
      throw new Exception("Filter output queue not empty!");
    
    if (!m_Filter.input(instance))
      throw new Exception(
	  "Filter didn't make the test instance immediately available!");
    
    m_Filter.batchFinished();
    Instance newInstance = m_Filter.output();

    return m_Clusterer.distributionForInstance(newInstance);
  }

  /**
   * Output a representation of this clusterer.
   * 
   * @return 		a representation of this clusterer
   */
  public String toString() {
    String 	result;
    
    if (m_FilteredInstances == null)
      result = "FilteredClusterer: No model built yet.";
    else
      result = "FilteredClusterer using "
	+ getClustererSpec()
	+ " on data filtered through "
	+ getFilterSpec()
	+ "\n\nFiltered Header\n"
	+ m_FilteredInstances.toString()
	+ "\n\nClusterer Model\n"
	+ m_Clusterer.toString();
    
    return result;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5488 $");
  }

  /**
   * Main method for testing this class.
   *
   * @param args 	the commandline options, use "-h" for help
   */
  public static void main(String [] args) {
    runClusterer(new FilteredClusterer(), args);
  }
}

