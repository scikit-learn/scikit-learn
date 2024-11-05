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
 * CheckSource.java
 * Copyright (C) 2007 University of Waikato, Hamilton, New Zealand
 */

package weka.filters;

import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.core.converters.ConverterUtils.DataSource;

import java.io.File;
import java.util.Enumeration;
import java.util.Vector;

/**
 * A simple class for checking the source generated from Filters
 * implementing the <code>weka.filters.Sourcable</code> interface.
 * It takes a filter, the classname of the generated source
 * and the dataset the source was generated with as parameters and tests
 * the output of the built filter against the output of the generated
 * source. Use option '-h' to display all available commandline options.
 * 
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -W &lt;classname and options&gt;
 *  The filter (incl. options) that was used to generate
 *  the source code.</pre>
 * 
 * <pre> -S &lt;classname&gt;
 *  The classname of the generated source code.</pre>
 * 
 * <pre> -t &lt;file&gt;
 *  The training set with which the source code was generated.</pre>
 * 
 * <pre> -c &lt;index&gt;
 *  The class index of the training set. 'first' and 'last' are
 *  valid indices.
 *  (default: none)</pre>
 * 
 <!-- options-end -->
 *
 * Options after -- are passed to the designated filter.<p>
 *
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 1.3 $
 * @see     weka.filters.Sourcable
 */
public class CheckSource
  implements OptionHandler, RevisionHandler {

  /** the classifier used for generating the source code */
  protected Filter m_Filter = null;
  
  /** the generated source code */
  protected Filter m_SourceCode = null;
  
  /** the dataset to use for testing */
  protected File m_Dataset = null;
  
  /** the class index */
  protected int m_ClassIndex = -1;
  
  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector result = new Vector();
    
    result.addElement(new Option(
        "\tThe filter (incl. options) that was used to generate\n"
        + "\tthe source code.",
        "W", 1, "-W <classname and options>"));
    
    result.addElement(new Option(
        "\tThe classname of the generated source code.",
        "S", 1, "-S <classname>"));
    
    result.addElement(new Option(
        "\tThe training set with which the source code was generated.",
        "t", 1, "-t <file>"));
    
    result.addElement(new Option(
        "\tThe class index of the training set. 'first' and 'last' are\n"
        + "\tvalid indices.\n"
        + "\t(default: none)",
        "c", 1, "-c <index>"));
    
    return result.elements();
  }

  /**
   * Parses a given list of options. <p/>
   *
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -W &lt;classname and options&gt;
   *  The filter (incl. options) that was used to generate
   *  the source code.</pre>
   * 
   * <pre> -S &lt;classname&gt;
   *  The classname of the generated source code.</pre>
   * 
   * <pre> -t &lt;file&gt;
   *  The training set with which the source code was generated.</pre>
   * 
   * <pre> -c &lt;index&gt;
   *  The class index of the training set. 'first' and 'last' are
   *  valid indices.
   *  (default: none)</pre>
   * 
   <!-- options-end -->
   *
   * Options after -- are passed to the designated filter.<p>
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String      tmpStr;
    String[]    spec;
    String      classname;

    tmpStr = Utils.getOption('W', options);
    if (tmpStr.length() > 0) {
      spec = Utils.splitOptions(tmpStr);
      if (spec.length == 0)
        throw new IllegalArgumentException("Invalid filter specification string");
      classname = spec[0];
      spec[0]   = "";
      setFilter((Filter) Utils.forName(Filter.class, classname, spec));
    }
    else {
      throw new Exception("No filter (classname + options) provided!");
    }
    
    tmpStr = Utils.getOption('S', options);
    if (tmpStr.length() > 0) {
      spec = Utils.splitOptions(tmpStr);
      if (spec.length != 1)
        throw new IllegalArgumentException("Invalid source code specification string");
      classname = spec[0];
      spec[0]   = "";
      setSourceCode((Filter) Utils.forName(Filter.class, classname, spec));
    }
    else {
      throw new Exception("No source code (classname) provided!");
    }

    tmpStr = Utils.getOption('t', options);
    if (tmpStr.length() != 0)
      setDataset(new File(tmpStr));
    else
      throw new Exception("No dataset provided!");

    tmpStr = Utils.getOption('c', options);
    if (tmpStr.length() != 0) {
      if (tmpStr.equals("first"))
        setClassIndex(0);
      else if (tmpStr.equals("last"))
        setClassIndex(-2);
      else 
        setClassIndex(Integer.parseInt(tmpStr) - 1);
    }
    else {
      setClassIndex(-1);
    }
  }

  /**
   * Gets the current settings of the filter.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    Vector<String>      result;
    
    result  = new Vector<String>();

    if (getFilter() != null) {
      result.add("-W");
      result.add(getFilter().getClass().getName() + " " 
          + Utils.joinOptions(((OptionHandler) getFilter()).getOptions()));
    }

    if (getSourceCode() != null) {
      result.add("-S");
      result.add(getSourceCode().getClass().getName());
    }

    if (getDataset() != null) {
      result.add("-t");
      result.add(m_Dataset.getAbsolutePath());
    }

    if (getClassIndex() != -1) {
      result.add("-c");
      if (getClassIndex() == -2)
	result.add("last");
      else if (getClassIndex() == 0)
	result.add("first");
      else 
	result.add("" + (getClassIndex() + 1));
    }
    
    return result.toArray(new String[result.size()]);
  }

  /**
   * Sets the filter to use for the comparison.
   * 
   * @param value       the filter to use
   */
  public void setFilter(Filter value) {
    m_Filter = value;
  }
  
  /**
   * Gets the filter being used for the tests, can be null.
   * 
   * @return            the currently set filter
   */
  public Filter getFilter() {
    return m_Filter;
  }
  
  /**
   * Sets the class to test.
   * 
   * @param value       the class to test
   */
  public void setSourceCode(Filter value) {
    m_SourceCode = value;
  }
  
  /**
   * Gets the class to test.
   * 
   * @return            the currently set class, can be null.
   */
  public Filter getSourceCode() {
    return m_SourceCode;
  }
  
  /**
   * Sets the dataset to use for testing.
   * 
   * @param value       the dataset to use.
   */
  public void setDataset(File value) {
    if (!value.exists())
      throw new IllegalArgumentException(
          "Dataset '" + value.getAbsolutePath() + "' does not exist!");
    else
      m_Dataset = value;
  }
  
  /**
   * Gets the dataset to use for testing, can be null.
   * 
   * @return            the dataset to use.
   */
  public File getDataset() {
    return m_Dataset;
  }

  /**
   * Sets the class index of the dataset.
   * 
   * @param value       the class index of the dataset.
   */
  public void setClassIndex(int value) {
    m_ClassIndex = value;
  }
  
  /**
   * Gets the class index of the dataset.
   * 
   * @return            the current class index.
   */
  public int getClassIndex() {
    return m_ClassIndex;
  }

  /**
   * compares two Instance
   * 
   * @param inst1	the first Instance object to compare
   * @param inst2	the second Instance object to compare
   * @return		true if both are the same
   */
  protected boolean compare(Instance inst1, Instance inst2) {
    boolean	result;
    int		i;
    
    // check dimension
    result = (inst1.numAttributes() == inst2.numAttributes());
    
    // check content
    if (result) {
      for (i = 0; i < inst1.numAttributes(); i++) {
	if (Double.isNaN(inst1.value(i)) && (Double.isNaN(inst2.value(i))))
	  continue;
	
	if (inst1.value(i) != inst2.value(i)) {
	  result = false;
	  System.out.println(
	      "Values at position " + (i+1) + " differ (Filter/Source code): " 
	      + inst1.value(i) + " != " + inst2.value(i));
	  break;
	}
      }
    }
    
    return result;
  }

  /**
   * compares the two Instances objects
   * 
   * @param inst1	the first Instances object to compare
   * @param inst2	the second Instances object to compare
   * @return		true if both are the same
   */
  protected boolean compare(Instances inst1, Instances inst2) {
    boolean	result;
    int		i;
    
    // check dimensions
    result = (inst1.numInstances() == inst2.numInstances());

    // check content
    if (result) {
      for (i = 0; i < inst1.numInstances(); i++) {
	result = compare(inst1.instance(i), inst2.instance(i));
	if (!result) {
	  System.out.println(
	      "Values in line " + (i+1) + " differ!");
	  break;
	}
      }
    }
    
    return result;
  }
  
  /**
   * performs the comparison test
   * 
   * @return            true if tests were successful
   * @throws Exception  if tests fail
   */
  public boolean execute() throws Exception {
    boolean     result;
    Instances   data;
    Instance	filteredInstance;
    Instances	filteredInstances;
    Instance	filteredInstanceSource;
    Instances	filteredInstancesSource;
    DataSource  source;
    Filter	filter;
    Filter	filterSource;
    int		i;
    
    result = true;
    
    // a few checks
    if (getFilter() == null)
      throw new Exception("No filter set!");
    if (getSourceCode() == null)
      throw new Exception("No source code set!");
    if (getDataset() == null)
      throw new Exception("No dataset set!");
    if (!getDataset().exists())
      throw new Exception(
          "Dataset '" + getDataset().getAbsolutePath() + "' does not exist!");
    
    // load data
    source = new DataSource(getDataset().getAbsolutePath());
    data   = source.getDataSet();
    if (getClassIndex() == -2)
      data.setClassIndex(data.numAttributes() - 1);
    else
      data.setClassIndex(getClassIndex());
    
    // compare output
    // 1. batch filtering
    filter = Filter.makeCopy(getFilter());
    filter.setInputFormat(data);
    filteredInstances = Filter.useFilter(data, filter);

    filterSource = Filter.makeCopy(getSourceCode());
    filterSource.setInputFormat(data);
    filteredInstancesSource = Filter.useFilter(data, filterSource);

    result = compare(filteredInstances, filteredInstancesSource);
    
    // 2. instance by instance
    if (result) {
      filter = Filter.makeCopy(getFilter());
      filter.setInputFormat(data);
      Filter.useFilter(data, filter);
      
      filterSource = Filter.makeCopy(getSourceCode());
      filterSource.setInputFormat(data);

      for (i = 0; i < data.numInstances(); i++) {
	filter.input(data.instance(i));
	filter.batchFinished();
	filteredInstance = filter.output();
	
	filterSource.input(data.instance(i));
	filterSource.batchFinished();
	filteredInstanceSource = filterSource.output();
	
	if (!compare(filteredInstance, filteredInstanceSource))
	  System.out.println(
	      (i+1) + ". instance (Filter/Source code): " 
	      + filteredInstance + " != " + filteredInstanceSource);
      }
    }
    
    return result;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 1.3 $");
  }
  
  /**
   * Executes the tests, use "-h" to list the commandline options.
   * 
   * @param args        the commandline parameters
   * @throws Exception  if something goes wrong
   */
  public static void main(String[] args) throws Exception{
    CheckSource         check;
    StringBuffer        text;
    Enumeration         enm;
    
    check = new CheckSource();
    if (Utils.getFlag('h', args)) {
      text = new StringBuffer();   
      text.append("\nHelp requested:\n\n");
      enm = check.listOptions();
      while (enm.hasMoreElements()) {
        Option option = (Option) enm.nextElement();
        text.append(option.synopsis() + "\n");
        text.append(option.description() + "\n");
      }
      System.out.println("\n" + text + "\n");
    }
    else {
      check.setOptions(args);
      if (check.execute())
        System.out.println("Tests OK!");
      else
        System.out.println("Tests failed!");
    }
  }
}
