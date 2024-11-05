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
 *    PolyKernel.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.functions.supportVector;

import weka.core.Capabilities;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.core.Capabilities.Capability;

import java.util.Enumeration;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * The polynomial kernel : K(x, y) = &lt;x, y&gt;^p or K(x, y) = (&lt;x, y&gt;+1)^p
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -D
 *  Enables debugging output (if available) to be printed.
 *  (default: off)</pre>
 * 
 * <pre> -no-checks
 *  Turns off all checks - use with caution!
 *  (default: checks on)</pre>
 * 
 * <pre> -C &lt;num&gt;
 *  The size of the cache (a prime number), 0 for full cache and 
 *  -1 to turn it off.
 *  (default: 250007)</pre>
 * 
 * <pre> -E &lt;num&gt;
 *  The Exponent to use.
 *  (default: 1.0)</pre>
 * 
 * <pre> -L
 *  Use lower-order terms.
 *  (default: no)</pre>
 * 
 <!-- options-end -->
 *
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @author Shane Legg (shane@intelligenesis.net) (sparse vector code)
 * @author Stuart Inglis (stuart@reeltwo.com) (sparse vector code)
 * @version $Revision: 5807 $
 */
public class PolyKernel 
  extends CachedKernel {

  /** for serialization */
  static final long serialVersionUID = -321831645846363201L;
  
  /** Use lower-order terms? */
  protected boolean m_lowerOrder = false;

  /** The exponent for the polynomial kernel. */
  protected double m_exponent = 1.0;

  /**
   * default constructor - does nothing.
   */
  public PolyKernel() {
    super();
  }

  /**
   * Frees the cache used by the kernel.
   */
  public void clean() {
    if (getExponent() == 1.0) {
      m_data = null;
    }    
    super.clean();
  }
  
  /**
   * Creates a new <code>PolyKernel</code> instance.
   * 
   * @param data	the training dataset used.
   * @param cacheSize	the size of the cache (a prime number)
   * @param exponent	the exponent to use
   * @param lowerOrder	whether to use lower-order terms
   * @throws Exception	if something goes wrong
   */
  public PolyKernel(Instances data, int cacheSize, double exponent,
		    boolean lowerOrder) throws Exception {
		
    super();
    
    setCacheSize(cacheSize);
    setExponent(exponent);
    setUseLowerOrder(lowerOrder);

    buildKernel(data);
  }
  
  /**
   * Returns a string describing the kernel
   * 
   * @return a description suitable for displaying in the
   *         explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "The polynomial kernel : K(x, y) = <x, y>^p or K(x, y) = (<x, y>+1)^p";
  }
  
  /**
   * Returns an enumeration describing the available options.
   *
   * @return 		an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector		result;
    Enumeration		en;
    
    result = new Vector();

    en = super.listOptions();
    while (en.hasMoreElements())
      result.addElement(en.nextElement());

    result.addElement(new Option(
	"\tThe Exponent to use.\n"
	+ "\t(default: 1.0)",
	"E", 1, "-E <num>"));

    result.addElement(new Option(
	"\tUse lower-order terms.\n"
	+ "\t(default: no)",
	"L", 0, "-L"));

    return result.elements();
  }

  /**
   * Parses a given list of options. <p/>
   * 
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -D
   *  Enables debugging output (if available) to be printed.
   *  (default: off)</pre>
   * 
   * <pre> -no-checks
   *  Turns off all checks - use with caution!
   *  (default: checks on)</pre>
   * 
   * <pre> -C &lt;num&gt;
   *  The size of the cache (a prime number), 0 for full cache and 
   *  -1 to turn it off.
   *  (default: 250007)</pre>
   * 
   * <pre> -E &lt;num&gt;
   *  The Exponent to use.
   *  (default: 1.0)</pre>
   * 
   * <pre> -L
   *  Use lower-order terms.
   *  (default: no)</pre>
   * 
   <!-- options-end -->
   * 
   * @param options 	the list of options as an array of strings
   * @throws Exception 	if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String	tmpStr;
    
    tmpStr = Utils.getOption('E', options);
    if (tmpStr.length() != 0)
      setExponent(Double.parseDouble(tmpStr));
    else
      setExponent(1.0);

    setUseLowerOrder(Utils.getFlag('L', options));
    
    super.setOptions(options);
  }

  /**
   * Gets the current settings of the Kernel.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    int       i;
    Vector    result;
    String[]  options;

    result = new Vector();
    options = super.getOptions();
    for (i = 0; i < options.length; i++)
      result.add(options[i]);

    result.add("-E");
    result.add("" + getExponent());

    if (getUseLowerOrder())
      result.add("-L");

    return (String[]) result.toArray(new String[result.size()]);	  
  }

  /**
   * 
   * @param id1   	the index of instance 1
   * @param id2		the index of instance 2
   * @param inst1	the instance 1 object
   * @return 		the dot product
   * @throws Exception 	if something goes wrong
   */
  protected double evaluate(int id1, int id2, Instance inst1)
    throws Exception {
		
    double result;
    if (id1 == id2) {
      result = dotProd(inst1, inst1);
    } else {
      result = dotProd(inst1, m_data.instance(id2));
    }
    // Use lower order terms?
    if (m_lowerOrder) {
      result += 1.0;
    }
    if (m_exponent != 1.0) {
      result = Math.pow(result, m_exponent);
    }
    return result;
  }

  /** 
   * Returns the Capabilities of this kernel.
   *
   * @return            the capabilities of this object
   * @see               Capabilities
   */
  public Capabilities getCapabilities() {
    Capabilities result = super.getCapabilities();
    result.disableAll();
    
    result.enable(Capability.NUMERIC_ATTRIBUTES);
    result.enableAllClasses();
    result.enable(Capability.MISSING_CLASS_VALUES);
    
    return result;
  }
  
  /**
   * Sets the exponent value.
   * 
   * @param value	the exponent value
   */
  public void setExponent(double value) {
    m_exponent = value;
  }
  
  /**
   * Gets the exponent value.
   * 
   * @return		the exponent value
   */
  public double getExponent() {
    return m_exponent;
  }

  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String exponentTipText() {
    return "The exponent value.";
  }
  
  /**
   * Sets whether to use lower-order terms.
   * 
   * @param value	true if lower-order terms will be used
   */
  public void setUseLowerOrder(boolean value) {
    m_lowerOrder = value;
  }
  
  /**
   * Gets whether lower-order terms are used.
   * 
   * @return		true if lower-order terms are used
   */
  public boolean getUseLowerOrder() {
    return m_lowerOrder;
  }

  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String useLowerOrderTipText() {
    return "Whether to use lower-order terms.";
  }
  
  /**
   * returns a string representation for the Kernel
   * 
   * @return 		a string representaiton of the kernel
   */
  public String toString() {
    String	result;
    
    if (getExponent() == 1.0) {
      if (getUseLowerOrder())
        result = "Linear Kernel with lower order: K(x,y) = <x,y> + 1";
      else
        result = "Linear Kernel: K(x,y) = <x,y>";
    }
    else {
      if (getUseLowerOrder())
	result = "Poly Kernel with lower order: K(x,y) = (<x,y> + 1)^" + getExponent();
      else
	result = "Poly Kernel: K(x,y) = <x,y>^" + getExponent();
    }
    
    return result;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5807 $");
  }
}

