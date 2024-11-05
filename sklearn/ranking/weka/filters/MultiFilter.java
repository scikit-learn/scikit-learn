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
 * MultiFilter.java
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.filters;

import weka.core.Capabilities;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.Utils;

import java.util.Enumeration;
import java.util.Vector;

/** 
 <!-- globalinfo-start -->
 * Applies several filters successively. In case all supplied filters are StreamableFilters, it will act as a streamable one, too.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -D
 *  Turns on output of debugging information.</pre>
 * 
 * <pre> -F &lt;classname [options]&gt;
 *  A filter to apply (can be specified multiple times).</pre>
 * 
 <!-- options-end -->
 *
 * @author  FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5501 $
 * @see     weka.filters.StreamableFilter
 */
public class MultiFilter
  extends SimpleStreamFilter {

  /** for serialization */
  private static final long serialVersionUID = -6293720886005713120L;

  /** The filters */
  protected Filter m_Filters[] = {new AllFilter()};

  /** caches the streamable state */
  protected boolean m_Streamable = false;

  /** whether we already checked the streamable state */
  protected boolean m_StreamableChecked = false;
  
  /**
   * Returns a string describing this filter
   * @return 		a description of the filter suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "Applies several filters successively. In case all supplied filters " 
      + "are StreamableFilters, it will act as a streamable one, too.";
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return 		an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector result = new Vector();
    Enumeration enm = super.listOptions();
    while (enm.hasMoreElements())
      result.add(enm.nextElement());
      
    result.addElement(new Option(
              "\tA filter to apply (can be specified multiple times).",
              "F", 1, "-F <classname [options]>"));

    return result.elements();
  }

  /**
   * Parses a list of options for this object. <p/>
   *
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -D
   *  Turns on output of debugging information.</pre>
   * 
   * <pre> -F &lt;classname [options]&gt;
   *  A filter to apply (can be specified multiple times).</pre>
   * 
   <!-- options-end -->
   *
   * @param options 	the list of options as an array of strings
   * @throws Exception 	if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String        tmpStr;
    String        filter;
    String[]      options2;
    Vector        filters;

    super.setOptions(options);
    
    filters = new Vector();
    while ((tmpStr = Utils.getOption("F", options)).length() != 0) {
      options2    = Utils.splitOptions(tmpStr);
      filter      = options2[0];
      options2[0] = "";
      filters.add(Utils.forName(Filter.class, filter, options2));
    }

    // at least one filter
    if (filters.size() == 0)
      filters.add(new AllFilter());

    setFilters((Filter[]) filters.toArray(new Filter[filters.size()]));
  }

  /**
   * Gets the current settings of the filter.
   *
   * @return 		an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    Vector        result;
    String[]      options;
    int           i;

    result = new Vector();

    options = super.getOptions();
    for (i = 0; i < options.length; i++)
      result.add(options[i]);
    
    for (i = 0; i < getFilters().length; i++) {
      result.add("-F");
      result.add(getFilterSpec(getFilter(i)));
    }

    return (String[]) result.toArray(new String[result.size()]);
  }

  /** 
   * Returns the Capabilities of this filter.
   *
   * @return            the capabilities of this object
   * @see               Capabilities
   */
  public Capabilities getCapabilities() {
    if (getFilters().length == 0) {
      Capabilities result = super.getCapabilities();
      result.disableAll();
      
      return result;
    } else {
      return getFilters()[0].getCapabilities();
    }
  }

  /**
   * resets the filter, i.e., m_NewBatch to true and m_FirstBatchDone to
   * false.
   *
   * @see #m_NewBatch
   * @see #m_FirstBatchDone
   */
  protected void reset() {
    super.reset();
    m_StreamableChecked = false;
  }

  /**
   * Sets the list of possible filters to choose from.
   * Also resets the state of the filter (this reset doesn't affect the 
   * options).
   *
   * @param filters	an array of filters with all options set.
   * @see #reset()
   */
  public void setFilters(Filter[] filters) {
    m_Filters           = filters;
    reset();
  }

  /**
   * Gets the list of possible filters to choose from.
   *
   * @return 	the array of Filters
   */
  public Filter[] getFilters() {
    return m_Filters;
  }
  
  /**
   * Returns the tip text for this property
   * @return    tip text for this property suitable for
   *            displaying in the explorer/experimenter gui
   */
  public String filtersTipText() {
    return "The base filters to be used.";
  }
  
  /**
   * Gets a single filter from the set of available filters.
   *
   * @param 	index the index of the filter wanted
   * @return 	the Filter
   */
  public Filter getFilter(int index) {
    return m_Filters[index];
  }

  /**
   * returns the filter classname and the options as one string
   * 
   * @param filter	the filter to get the specs for
   * @return		the classname plus options
   */
  protected String getFilterSpec(Filter filter) {
    String        result;

    if (filter == null) {
      result = "";
    }
    else {
      result  = filter.getClass().getName();
      if (filter instanceof OptionHandler)
        result += " " 
          + Utils.joinOptions(((OptionHandler) filter).getOptions());
    }

    return result;
  }

  /**
   * tests whether all the enclosed filters are streamable
   *
   * @return	true if all the enclosed filters are streamable
   */
  public boolean isStreamableFilter() {
    int           i;

    if (!m_StreamableChecked) {
      m_Streamable        = true;
      m_StreamableChecked = true;

      for (i = 0; i < getFilters().length; i++) {
        if (getFilter(i) instanceof MultiFilter)
          m_Streamable = ((MultiFilter) getFilter(i)).isStreamableFilter();
        else if (getFilter(i) instanceof StreamableFilter)
          m_Streamable = true;
        else
          m_Streamable = false;

        if (!m_Streamable)
          break;
      }

      if (getDebug())
        System.out.println("Streamable: " + m_Streamable);
    }

    return m_Streamable;
  }

  /**
   * Returns true if the output format is immediately available after the
   * input format has been set and not only after all the data has been
   * seen (see batchFinished()). This method should normally return true
   * for a stream filter, since the data will be processed in a batch
   * manner instead (or at least for the second batch of files, see
   * m_FirstBatchDone).
   *
   * @return      true if the output format is immediately available
   * @see         #batchFinished()
   * @see         #setInputFormat(Instances)
   * @see         #m_FirstBatchDone
   */
  protected boolean hasImmediateOutputFormat() {
    return isStreamableFilter();
  }
  
  /**
   * Determines the output format based on the input format and returns 
   * this. In case the output format cannot be returned immediately, i.e.,
   * hasImmediateOutputFormat() returns false, then this method will called
   * from batchFinished() after the call of preprocess(Instances), in which,
   * e.g., statistics for the actual processing step can be gathered.
   *
   * @param inputFormat     the input format to base the output format on
   * @return                the output format
   * @throws Exception      in case the determination goes wrong
   * @see                   #hasImmediateOutputFormat()
   * @see                   #batchFinished()
   * @see                   #preprocess(Instances)
   */
  protected Instances determineOutputFormat(Instances inputFormat) throws Exception {
    Instances   result;
    int         i;
    
    result = getInputFormat();
    
    for (i = 0; i < getFilters().length; i++) {
      if (!isFirstBatchDone())
	getFilter(i).setInputFormat(result);
      result = getFilter(i).getOutputFormat();
    }

    return result;
  }

  /**
   * processes the given instance (may change the provided instance) and
   * returns the modified version.
   *
   * @param instance    the instance to process
   * @return            the modified data
   * @throws Exception  in case the processing goes wrong
   */
  protected Instance process(Instance instance) throws Exception {
    Instance    result;
    int         i;
    
    result = (Instance) instance.copy();

    for (i = 0; i < getFilters().length; i++) {
      getFilter(i).input(result);
      result = getFilter(i).output();
    }

    return result;
  }

  /**
   * Processes the given data (may change the provided dataset) and returns
   * the modified version. This method is called in batchFinished().
   * This implementation only calls process(Instance) for each instance
   * in the given dataset.
   *
   * @param instances   the data to process
   * @return            the modified data
   * @throws Exception  in case the processing goes wrong
   * @see               #batchFinished()
   * @see               #process(Instance)
   */
  protected Instances process(Instances instances) throws Exception {
    Instances     result;
    int           i;

    result = instances;
    
    for (i = 0; i < getFilters().length; i++) {
      if (!isFirstBatchDone())
	getFilter(i).setInputFormat(result);
      result = Filter.useFilter(result, getFilter(i));
    }
    
    return result;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5501 $");
  }

  /**
   * Main method for executing this class.
   *
   * @param args should contain arguments for the filter: use -h for help
   */
  public static void main(String[] args) {
    runFilter(new MultiFilter(), args);
  }
}
