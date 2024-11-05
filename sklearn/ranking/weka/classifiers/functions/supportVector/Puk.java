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
 *    Puk.java
 *    Copyright (C) 2007 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.functions.supportVector;

import weka.core.Capabilities;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.RevisionUtils;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformationHandler;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;

import java.util.Enumeration;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * The Pearson VII function-based universal kernel.<br/>
 * <br/>
 * For more information see:<br/>
 * <br/>
 * B. Uestuen, W.J. Melssen, L.M.C. Buydens (2006). Facilitating the application of Support Vector Regression by using a universal Pearson VII function based kernel. Chemometrics and Intelligent Laboratory Systems. 81:29-40.
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
 * <pre> -O &lt;num&gt;
 *  The Omega parameter.
 *  (default: 1.0)</pre>
 * 
 * <pre> -S &lt;num&gt;
 *  The Sigma parameter.
 *  (default: 1.0)</pre>
 * 
 <!-- options-end -->
 *
 * @author Bernhard Pfahringer (bernhard@cs.waikato.ac.nz)
 * @version $Revision: 5450 $
 */
public class Puk 
  extends CachedKernel
  implements TechnicalInformationHandler {
  
  /** for serialization */
  private static final long serialVersionUID = 1682161522559978851L;

  /** The precalculated dotproducts of &lt;inst_i,inst_i&gt; */
  protected double m_kernelPrecalc[];

  /** Omega for the Puk kernel. */
  protected double m_omega = 1.0;

  /** Sigma for the Puk kernel. */
  protected double m_sigma = 1.0;

  /** Cached factor for the Puk kernel. */
  protected double m_factor = 1.0;

  /**
   * default constructor - does nothing.
   */
  public Puk() {
    super();
  }
  
  /**
   * Constructor. Initializes m_kernelPrecalc[].
   * 
   * @param data	the data to use
   * @param cacheSize	the size of the cache
   * @param omega	the exponent
   * @param sigma	the bandwidth
   * @throws Exception	if something goes wrong
   */
  public Puk(Instances data, int cacheSize, double omega, double sigma)
    throws Exception {

    super();
    
    setCacheSize(cacheSize);
    setOmega(omega);
    setSigma(sigma);
    
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
        "The Pearson VII function-based universal kernel.\n\n"
      + "For more information see:\n\n"
      + getTechnicalInformation().toString();
  }

  /**
   * Returns an instance of a TechnicalInformation object, containing 
   * detailed information about the technical background of this class,
   * e.g., paper reference or book this class is based on.
   * 
   * @return the technical information about this class
   */
  public TechnicalInformation getTechnicalInformation() {
    TechnicalInformation 	result;
    
    result = new TechnicalInformation(Type.ARTICLE);
    result.setValue(Field.AUTHOR, "B. Uestuen and W.J. Melssen and L.M.C. Buydens");
    result.setValue(Field.YEAR, "2006");
    result.setValue(Field.TITLE, "Facilitating the application of Support Vector Regression by using a universal Pearson VII function based kernel");
    result.setValue(Field.JOURNAL, "Chemometrics and Intelligent Laboratory Systems");
    result.setValue(Field.VOLUME, "81");
    result.setValue(Field.PAGES, "29-40");
    result.setValue(Field.PDF, "http://www.cac.science.ru.nl/research/publications/PDFs/ustun2006.pdf");

    return result;
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
	"\tThe Omega parameter.\n"
	+ "\t(default: 1.0)",
	"O", 1, "-O <num>"));

    result.addElement(new Option(
	"\tThe Sigma parameter.\n"
	+ "\t(default: 1.0)",
	"S", 1, "-S <num>"));

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
   * <pre> -O &lt;num&gt;
   *  The Omega parameter.
   *  (default: 1.0)</pre>
   * 
   * <pre> -S &lt;num&gt;
   *  The Sigma parameter.
   *  (default: 1.0)</pre>
   * 
   <!-- options-end -->
   * 
   * @param options 	the list of options as an array of strings
   * @throws Exception 	if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String	tmpStr;
    
    tmpStr = Utils.getOption('O', options);
    if (tmpStr.length() != 0)
      setOmega(Double.parseDouble(tmpStr));
    else
      setOmega(1.0);
    
    tmpStr = Utils.getOption('S', options);
    if (tmpStr.length() != 0)
      setSigma(Double.parseDouble(tmpStr));
    else
      setSigma(1.0);
    
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

    result.add("-O");
    result.add("" + getOmega());

    result.add("-S");
    result.add("" + getSigma());

    return (String[]) result.toArray(new String[result.size()]);	  
  }

  /**
   * returns the dot product
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
      double squaredDifference = -2.0 * dotProd(inst1, inst2) + precalc1 + m_kernelPrecalc[id2];
      double intermediate = m_factor * Math.sqrt(squaredDifference);
      double result = 1.0 / Math.pow(1.0 + intermediate * intermediate, getOmega());
      return result;
    }
  }
    
  /**
   * Sets the omega value.
   * 
   * @param value	the omega value
   */
  public void setOmega(double value) {
    m_omega  = value;
    m_factor = computeFactor(m_omega, m_sigma);
  }
  
  /**
   * Gets the omega value.
   * 
   * @return		the omega value
   */
  public double getOmega() {
    return m_omega;
  }

  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String omegaTipText() {
    return "The Omega value.";
  }

  /**
   * Sets the sigma value.
   * 
   * @param value	the sigma value
   */
  public void setSigma(double value) {
    m_sigma  = value;
    m_factor = computeFactor(m_omega, m_sigma);
  }
  
  /**
   * Gets the sigma value.
   * 
   * @return		the sigma value
   */
  public double getSigma() {
    return m_sigma;
  }

  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String sigmaTipText() {
    return "The Sigma value.";
  }

  /**
   * computes the factor for curve-fitting (see equation (13) in paper)
   * 
   * @param omega	the omega to use
   * @param sigma	the sigma to use
   * @return		the factor for curve-fitting
   */
  protected double computeFactor(double omega, double sigma) {
    double root = Math.sqrt(Math.pow(2.0, 1.0 / omega) - 1);
    return 2.0 * root / sigma;
  }

  /**
   * initializes variables etc.
   * 
   * @param data	the data to use
   */
  protected void initVars(Instances data) {
    super.initVars(data);
    
    m_factor        = computeFactor(m_omega, m_sigma);
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
    return "Puk kernel";
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
