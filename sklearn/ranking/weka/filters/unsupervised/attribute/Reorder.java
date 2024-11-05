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
 * Reorder.java
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.filters.unsupervised.attribute;

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
import weka.core.Capabilities.Capability;
import weka.core.labelranking.PreferenceDenseInstance;
import weka.filters.Filter;
import weka.filters.StreamableFilter;
import weka.filters.UnsupervisedFilter;

import java.util.Enumeration;
import java.util.StringTokenizer;
import java.util.Vector;

/** 
 <!-- globalinfo-start -->
 * A filter that generates output with a new order of the attributes. Useful if one wants to move an attribute to the end to use it as class attribute (e.g. with using "-R 2-last,1").<br/>
 * But it's not only possible to change the order of all the attributes, but also to leave out attributes. E.g. if you have 10 attributes, you can generate the following output order: 1,3,5,7,9,10 or 10,1-5.<br/>
 * You can also duplicate attributes, e.g. for further processing later on: e.g. 1,1,1,4,4,4,2,2,2 where the second and the third column of each attribute are processed differently and the first one, i.e. the original one is kept.<br/>
 * One can simply inverse the order of the attributes via 'last-first'.<br/>
 * After appyling the filter, the index of the class attribute is the last attribute.
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -R &lt;index1,index2-index4,...&gt;
 *  Specify list of columns to copy. First and last are valid
 *  indexes. (default first-last)</pre>
 * 
 <!-- options-end -->
 *
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 6250 $
 */
public class Reorder 
  extends Filter 
  implements UnsupervisedFilter, StreamableFilter, OptionHandler {
  
  /** for serialization */
  static final long serialVersionUID = -1135571321097202292L;

  /** Stores which columns to reorder */
  protected String m_NewOrderCols = "first-last";

  /**
   * Stores the indexes of the selected attributes in order, once the
   * dataset is seen
   */
  protected int[] m_SelectedAttributes;

  /** 
   * Contains an index of string attributes in the input format
   * that survive the filtering process -- some entries may be duplicated 
   */
  protected int[] m_InputStringIndex;

  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector newVector = new Vector();

    newVector.addElement(new Option(
              "\tSpecify list of columns to copy. First and last are valid\n"
	      +"\tindexes. (default first-last)",
              "R", 1, "-R <index1,index2-index4,...>"));

    return newVector.elements();
  }

  /**
   * Parses a given list of options. <p/>
   * 
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -R &lt;index1,index2-index4,...&gt;
   *  Specify list of columns to copy. First and last are valid
   *  indexes. (default first-last)</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String orderList = Utils.getOption('R', options);
    if (orderList.length() != 0) {
      setAttributeIndices(orderList);
    }
    
    if (getInputFormat() != null) {
      setInputFormat(getInputFormat());
    }
  }

  /**
   * Gets the current settings of the filter.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String [] getOptions() {
    String [] options = new String [2];
    int current = 0;

    if (!getAttributeIndices().equals("")) {
      options[current++] = "-R"; 
      options[current++] = getAttributeIndices();
    }

    while (current < options.length) {
      options[current++] = "";
    }
    return options;
  }
  
  /**
   * parses the index string and returns the corresponding int index
   * 
   * @param s			the index string to parse
   * @param numAttributes	necessary for "last" and OutOfBounds checks
   * @return			the int index determined form the index string
   * @throws Exception		if index is not valid
   */
  protected int determineIndex(String s, int numAttributes) throws Exception {
    int		result;
    
    if (s.equals("first"))
      result = 0;
    else if (s.equals("last"))
      result = numAttributes - 1;
    else
      result = Integer.parseInt(s) - 1;
    
    // out of bounds?
    if ( (result < 0) || (result > numAttributes - 1) )
      throw new IllegalArgumentException(
	  "'" + s + "' is not a valid index for the range '1-" + numAttributes + "'!");
    
    return result;
  }
  
  /**
   * parses the range string and returns an array with the indices
   * 
   * @param numAttributes	necessary for "last" and OutOfBounds checks
   * @return			the indices determined form the range string
   * @see			#m_NewOrderCols
   * @throws Exception		if range is not valid
   */
  protected int[] determineIndices(int numAttributes) throws Exception {
    int[]		result;
    Vector<Integer>	list;
    int			i;
    StringTokenizer	tok;
    String		token;
    String[]		range;
    int			from;
    int			to;
    
    list = new Vector<Integer>();
    
    // parse range
    tok = new StringTokenizer(m_NewOrderCols, ",");
    while (tok.hasMoreTokens()) {
      token = tok.nextToken();
      if (token.indexOf("-") > -1) {
	range = token.split("-");
	if (range.length != 2)
	  throw new IllegalArgumentException("'" + token + "' is not a valid range!");
	from = determineIndex(range[0], numAttributes);
	to   = determineIndex(range[1], numAttributes);

	if (from <= to) {
	  for (i = from; i <= to; i++)
	    list.add(i);
	}
	else {
	  for (i = from; i >= to; i--)
	    list.add(i);
	}
      }
      else {
	list.add(determineIndex(token, numAttributes));
      }
    }
    
    // turn vector into int array
    result = new int[list.size()];
    for (i = 0; i < list.size(); i++)
      result[i] = list.get(i);
    
    return result;
  }

  /** 
   * Returns the Capabilities of this filter.
   *
   * @return            the capabilities of this object
   * @see               Capabilities
   */
  public Capabilities getCapabilities() {
    Capabilities result = super.getCapabilities();
    result.disableAll();

    // attribute
    result.enableAllAttributes();
    result.enable(Capability.MISSING_VALUES);
    
    // class
    result.enableAllClasses();
    result.enable(Capability.NO_CLASS);
    result.enable(Capability.MISSING_CLASS_VALUES);
    
    return result;
  }

  /**
   * Sets the format of the input instances.
   *
   * @param instanceInfo an Instances object containing the input instance
   * structure (any instances contained in the object are ignored - only the
   * structure is required).
   * @return true if the outputFormat may be collected immediately
   * @throws Exception if a problem occurs setting the input format
   */
  public boolean setInputFormat(Instances instanceInfo) throws Exception {
    super.setInputFormat(instanceInfo);
    
    FastVector attributes = new FastVector();
    int outputClass = -1;
    m_SelectedAttributes = determineIndices(instanceInfo.numAttributes());
    for (int i = 0; i < m_SelectedAttributes.length; i++) {
      int current = m_SelectedAttributes[i];
      if (instanceInfo.classIndex() == current) {
	outputClass = attributes.size();
      }
      Attribute keep = (Attribute)instanceInfo.attribute(current).copy();
      attributes.addElement(keep);
    }
    
    initInputLocators(instanceInfo, m_SelectedAttributes);

    Instances outputFormat = new Instances(instanceInfo.relationName(),
					   attributes, 0); 
    outputFormat.setClassIndex(outputClass);
    setOutputFormat(outputFormat);
    
    return true;
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
  public boolean input(Instance instance) {
    if (getInputFormat() == null) {
      throw new IllegalStateException("No input instance format defined");
    }
    if (m_NewBatch) {
      resetQueue();
      m_NewBatch = false;
    }

    double[] vals = new double[outputFormatPeek().numAttributes()];
    for (int i = 0; i < m_SelectedAttributes.length; i++) {
      int current = m_SelectedAttributes[i];
      vals[i] = instance.value(current);
    }
    Instance inst = null;
    if (instance instanceof SparseInstance)
      inst = new SparseInstance(instance.weight(), vals);
    //RANKING BEGIN
    if(instance instanceof PreferenceDenseInstance){
    	PreferenceDenseInstance pdi = (PreferenceDenseInstance) instance;
    	inst = new PreferenceDenseInstance(instance.weight(), vals, pdi.getHashMap());
    }
    //RANKING END
    else
      inst = new DenseInstance(instance.weight(), vals);

    inst.setDataset(getOutputFormat());
    copyValues(inst, false, instance.dataset(), getOutputFormat());
    inst.setDataset(getOutputFormat());
    
    push(inst);
    
    return true;
  }

  /**
   * Returns a string describing this filter
   *
   * @return a description of the filter suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "A filter that generates output with a new order of the "
      + "attributes. Useful if one wants to move an attribute to the end to "
      + "use it as class attribute (e.g. with using \"-R 2-last,1\").\n"
      + "But it's not only possible to change the order of all the attributes, "
      + "but also to leave out attributes. E.g. if you have 10 attributes, you "
      + "can generate the following output order: 1,3,5,7,9,10 or 10,1-5.\n"
      + "You can also duplicate attributes, e.g. for further processing later "
      + "on: e.g. 1,1,1,4,4,4,2,2,2 where the second and the third column of "
      + "each attribute are processed differently and the first one, i.e. the "
      + "original one is kept.\n"
      + "One can simply inverse the order of the attributes via 'last-first'.\n"
      + "After appyling the filter, the index of the class attribute is the "
      + "last attribute.";
  }

  /**
   * Get the current range selection
   *
   * @return a string containing a comma separated list of ranges
   */
  public String getAttributeIndices() {
    return m_NewOrderCols;
  }

  /**
   * Returns the tip text for this property
   *
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String attributeIndicesTipText() {
    return "Specify range of attributes to act on."
      + " This is a comma separated list of attribute indices, with"
      + " \"first\" and \"last\" valid values. Specify an inclusive"
      + " range with \"-\". E.g: \"first-3,5,6-10,last\".";
  }

  /**
   * Set which attributes are to be copied (or kept if invert is true)
   *
   * @param rangeList a string representing the list of attributes.  Since
   * the string will typically come from a user, attributes are indexed from
   * 1. <br>
   * eg: first-3,5,6-last<br>
   * Note: use this method before you call 
   * <code>setInputFormat(Instances)</code>, since the output format is
   * determined in that method.
   * @throws Exception if an invalid range list is supplied
   */
  public void setAttributeIndices(String rangeList) throws Exception {
    // simple test
    if (rangeList.replaceAll("[afilrst0-9\\-,]*", "").length() != 0)
      throw new IllegalArgumentException("Not a valid range string!");
    
    m_NewOrderCols = rangeList;
  }

  /**
   * Set which attributes are to be copied (or kept if invert is true)
   *
   * @param attributes an array containing indexes of attributes to select.
   * Since the array will typically come from a program, attributes are indexed
   * from 0.<br>
   * Note: use this method before you call 
   * <code>setInputFormat(Instances)</code>, since the output format is
   * determined in that method.
   * @throws Exception if an invalid set of ranges is supplied
   */
  public void setAttributeIndicesArray(int [] attributes) throws Exception {
    setAttributeIndices(Range.indicesToRangeList(attributes));
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6250 $");
  }

  /**
   * Main method for testing this class.
   *
   * @param argv should contain arguments to the filter: use -h for help
   */
  public static void main(String [] argv) {
    runFilter(new Reorder(), argv);
  }
}
