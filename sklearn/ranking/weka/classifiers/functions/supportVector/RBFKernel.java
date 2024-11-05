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
 *    RBFKernel.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *    Copyright (C) 2005 J. Lindgren
 *
 */

package weka.classifiers.functions.supportVector;

import weka.core.Capabilities;
import weka.core.Capabilities.Capability;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.RevisionUtils;
import weka.core.Utils;

import java.util.Enumeration;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * The RBF kernel. K(x, y) = e^-(gamma * &lt;x-y, x-y&gt;^2)
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
 * <pre> -G &lt;num&gt;
 *  The Gamma parameter.
 *  (default: 0.01)</pre>
 * 
 <!-- options-end -->
 *
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @author Shane Legg (shane@intelligenesis.net) (sparse vector code)
 * @author Stuart Inglis (stuart@reeltwo.com) (sparse vector code)
 * @author J. Lindgren (jtlindgr{at}cs.helsinki.fi) (RBF kernel)
 * @version $Revision: 5450 $
 */
public class RBFKernel 
  extends CachedKernel {
  
  /** for serialization */
  static final long serialVersionUID = 5247117544316387852L;

  /** The precalculated dotproducts of &lt;inst_i,inst_i&gt; */
  protected double m_kernelPrecalc[];

  /** Gamma for the RBF kernel. */
  protected double m_gamma = 0.01;

  /**
   * default constructor - does nothing.
   */
  public RBFKernel() {
    super();
  }
  
  /**
   * Constructor. Initializes m_kernelPrecalc[].
   * 
   * @param data	the data to use
   * @param cacheSize	the size of the cache
   * @param gamma	the bandwidth
   * @throws Exception	if something goes wrong
   */
  public RBFKernel(Instances data, int cacheSize, double gamma)
    throws Exception {

    super();
    
    setCacheSize(cacheSize);
    setGamma(gamma);
    
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
        "The RBF kernel. K(x, y) = e^-(gamma * <x-y, x-y>^2)";
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
	"\tThe Gamma parameter.\n"
	+ "\t(default: 0.01)",
	"G", 1, "-G <num>"));

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
   * <pre> -G &lt;num&gt;
   *  The Gamma parameter.
   *  (default: 0.01)</pre>
   * 
   <!-- options-end -->
   * 
   * @param options 	the list of options as an array of strings
   * @throws Exception 	if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String	tmpStr;
    
    tmpStr = Utils.getOption('G', options);
    if (tmpStr.length() != 0)
      setGamma(Double.parseDouble(tmpStr));
    else
      setGamma(0.01);
    
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

    result.add("-G");
    result.add("" + getGamma());

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

    if (id1 == id2) {
      return 1.0;
    } else {
      double precalc1;
      if (id1 == -1)
        precalc1 = dotProd(inst1, inst1);
      else
        precalc1 = m_kernelPrecalc[id1];
      Instance inst2 = m_data.instance(id2);
      double result = Math.exp(m_gamma
  			     * (2. * dotProd(inst1, inst2) - precalc1 - m_kernelPrecalc[id2]));
      
      return result;
    }
  }
    
  /**
   * Sets the gamma value.
   * 
   * @param value	the gamma value
   */
  public void setGamma(double value) {
    m_gamma = value;
  }
  
  /**
   * Gets the gamma value.
   * 
   * @return		the gamma value
   */
  public double getGamma() {
    return m_gamma;
  }

  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String gammaTipText() {
    return "The Gamma value.";
  }

  /**
   * initializes variables etc.
   * 
   * @param data	the data to use
   */
  protected void initVars(Instances data) {
    super.initVars(data);
    
    m_kernelPrecalc = new double[data.numInstances()];
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
   * builds the kernel with the given data. Initializes the kernel cache. 
   * The actual size of the cache in bytes is (64 * cacheSize).
   * 
   * @param data	the data to base the kernel on
   * @throws Exception	if something goes wrong
   */
  public void buildKernel(Instances data) throws Exception {
    // does kernel handle the data?
    if (!getChecksTurnedOff())
      getCapabilities().testWithFail(data);
    
    initVars(data);
    
    for (int i = 0; i < data.numInstances(); i++)
      m_kernelPrecalc[i] = dotProd(data.instance(i), data.instance(i));
  }
  
  /**
   * returns a string representation for the Kernel
   * 
   * @return 		a string representaiton of the kernel
   */
  public String toString() {
    return "RBF kernel: K(x,y) = e^-(" + getGamma() + "* <x-y,x-y>^2)";
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5450 $");
  }
}

