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
 * RandomSubset.java
 * Copyright (C) 2007 University of Waikato, Hamilton, New Zealand
 */

package weka.filters.unsupervised.attribute;

import weka.core.Capabilities;
import weka.core.FastVector;
import weka.core.Instance; 
import weka.core.DenseInstance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.filters.SimpleStreamFilter;

import java.util.Collections;
import java.util.Enumeration;
import java.util.Random;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * Chooses a random subset of attributes, either an absolute number or a percentage. The class is always included in the output (as the last attribute).
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -D
 *  Turns on output of debugging information.</pre>
 * 
 * <pre> -N &lt;double&gt;
 *  The number of attributes to randomly select.
 *  If &lt; 1 then percentage, &gt;= 1 absolute number.
 *  (default: 0.5)</pre>
 * 
 * <pre> -S &lt;int&gt;
 *  The seed value.
 *  (default: 1)</pre>
 * 
 <!-- options-end -->
 *
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5987 $
 */
public class RandomSubset
  extends SimpleStreamFilter {

  /** for serialization. */
  private static final long serialVersionUID = 2911221724251628050L;

  /** The number of attributes to randomly choose (&gt;= 1 absolute number of
   * attributes, &lt; 1 percentage). */
  protected double m_NumAttributes = 0.5;
  
  /** The seed value. */
  protected int m_Seed = 1;
  
  /** The indices of the attributes that got selected. */
  protected int[] m_Indices = null;
  
  /**
   * Returns a string describing this filter.
   *
   * @return 		a description of the filter suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "Chooses a random subset of attributes, either an absolute number "
      + "or a percentage. The class is always included in the output ("
      + "as the last attribute).";
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector        result;
    Enumeration   enm;

    result = new Vector();

    enm = super.listOptions();
    while (enm.hasMoreElements())
      result.addElement(enm.nextElement());

    result.addElement(new Option(
	"\tThe number of attributes to randomly select.\n"
	+ "\tIf < 1 then percentage, >= 1 absolute number.\n"
	+ "\t(default: 0.5)",
	"N", 1, "-N <double>"));
    
    result.addElement(new Option(
	"\tThe seed value.\n"
	+ "\t(default: 1)",
	"S", 1, "-S <int>"));

    return result.elements();
  }	  

  /**
   * Gets the current settings of the filter.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    int			i;
    Vector<String>	result;
    String[]		options;

    result  = new Vector<String>();
    options = super.getOptions();
    for (i = 0; i < options.length; i++)
      result.add(options[i]);

    result.add("-N"); 
    result.add("" + m_NumAttributes);

    result.add("-S"); 
    result.add("" + m_Seed);

    return result.toArray(new String[result.size()]);	  
  }	  

  /**
   * Parses a given list of options. <p/>
   *
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -D
   *  Turns on output of debugging information.</pre>
   * 
   * <pre> -N &lt;double&gt;
   *  The number of attributes to randomly select.
   *  If &lt; 1 then percentage, &gt;= 1 absolute number.
   *  (default: 0.5)</pre>
   * 
   * <pre> -S &lt;int&gt;
   *  The seed value.
   *  (default: 1)</pre>
   * 
   <!-- options-end -->
   * 
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported 
   */
  public void setOptions(String[] options) throws Exception {
    String	tmpStr;

    tmpStr = Utils.getOption("N", options);
    if (tmpStr.length() != 0)
      setNumAttributes(Double.parseDouble(tmpStr));
    else
      setNumAttributes(0.5);
    
    tmpStr = Utils.getOption("S", options);
    if (tmpStr.length() != 0)
      setSeed(Integer.parseInt(tmpStr));
    else
      setSeed(1);
    
    super.setOptions(options);
  }	  

  /**
   * Returns the tip text for this property.
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String numAttributesTipText() {
    return "The number of attributes to choose: < 1 percentage, >= 1 absolute number.";
  }

  /**
   * Get the number of attributes (&lt; 1 percentage, &gt;= 1 absolute number).
   *
   * @return 		the number of attributes.
   */
  public double getNumAttributes() {
    return m_NumAttributes;
  }

  /**
   * Set the number of attributes. 
   *
   * @param value	the number of attributes to use.
   */
  public void setNumAttributes(double value) {
    m_NumAttributes = value;
  }

  /**
   * Returns the tip text for this property.
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String seedTipText() {
    return "The seed value for the random number generator.";
  }

  /**
   * Get the seed value for the random number generator.
   *
   * @return 		the seed value.
   */
  public int getSeed() {
    return m_Seed;
  }

  /**
   * Set the seed value for the random number generator. 
   *
   * @param value	the seed value.
   */
  public void setSeed(int value) {
    m_Seed = value;
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

    // attributes
    result.enableAllAttributes();
    result.enable(Capability.MISSING_VALUES);
    
    // class
    result.enableAllClasses();
    result.enable(Capability.MISSING_CLASS_VALUES);
    result.enable(Capability.NO_CLASS);
    
    //RANKING BEGIN
    result.disable(Capability.RANKING);
    result.disable(Capability.PREFERENCE_ATTRIBUTE);
    //RANKING END
    
    return result;
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
   */
  protected Instances determineOutputFormat(Instances inputFormat) throws Exception {
    Instances		result;
    FastVector		atts;
    int			i;
    int			numAtts;
    Vector<Integer>	indices;
    Vector<Integer>	subset;
    Random		rand;
    int			index;
 
    // determine the number of attributes
    numAtts = inputFormat.numAttributes();
    if (inputFormat.classIndex() > -1)
      numAtts--;
    
    if (m_NumAttributes < 1) {
      numAtts = (int) Math.round((double) numAtts * m_NumAttributes);
    }
    else {
      if (m_NumAttributes < numAtts)
	numAtts = (int) m_NumAttributes;
    }
    if (getDebug())
      System.out.println("# of atts: " + numAtts);
    
    // determine random indices
    indices = new Vector<Integer>();
    for (i = 0; i < inputFormat.numAttributes(); i++) {
      if (i == inputFormat.classIndex())
	continue;
      indices.add(i);
    }
    
    subset = new Vector<Integer>();
    rand   = new Random(m_Seed);
    for (i = 0; i < numAtts; i++) {
      index = rand.nextInt(indices.size());
      subset.add(indices.get(index));
      indices.remove(index);
    }
    Collections.sort(subset);
    if (inputFormat.classIndex() > -1)
      subset.add(inputFormat.classIndex());
    if (getDebug())
      System.out.println("indices: " + subset);
    
    // generate output format
    atts      = new FastVector();
    m_Indices = new int[subset.size()];
    for (i = 0; i < subset.size(); i++) {
      atts.addElement(inputFormat.attribute(subset.get(i)));
      m_Indices[i] = subset.get(i);
    }
    result = new Instances(inputFormat.relationName(), atts, 0);
    if (inputFormat.classIndex() > -1)
      result.setClassIndex(result.numAttributes() - 1);

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
    Instance	result;
    double[]	values;
    int		i;
    
    values = new double[m_Indices.length];
    for (i = 0; i < m_Indices.length; i++)
      values[i] = instance.value(m_Indices[i]);

    result = new DenseInstance(instance.weight(), values);
    result.setDataset(getOutputFormat());
    
    copyValues(result, false, instance.dataset(), getOutputFormat());
    result.setDataset(getOutputFormat());
    
    return result;
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
   * Runs the filter with the given parameters. Use -h to list options.
   * 
   * @param args	the commandline options
   */
  public static void main(String[] args) {
    runFilter(new RandomSubset(), args);
  }
}
