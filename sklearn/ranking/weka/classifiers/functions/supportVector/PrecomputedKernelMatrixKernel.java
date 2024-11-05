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
 *    PrecomputedKernelMatrixKernel.java
 *    Copyright (C) 2008 University of Waikato, Hamilton, New Zealand
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
import weka.core.matrix.Matrix;

import java.io.File;
import java.io.FileReader;
import java.util.Enumeration;
import java.util.Vector;

/**
 *
 <!-- globalinfo-start -->
 * This kernel is based on a static kernel matrix that is read from a file. Instances must have a single nominal attribute (excluding the class). This attribute must be the first attribute in the file and its values are used to reference rows/columns in the kernel matrix. The second attribute must be the class attribute.
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
 * <pre> -M &lt;file name&gt;
 *  The file name of the file that holds the kernel matrix.
 *  (default: kernelMatrix.matrix)</pre>
 * 
 <!-- options-end -->
 *
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version $Revision: 5450 $
 */
public class PrecomputedKernelMatrixKernel extends Kernel {

  /** for serialization */
  static final long serialVersionUID = -321831645846363333L;

  /** The file holding the kernel matrix. */
  protected File m_KernelMatrixFile = new File("kernelMatrix.matrix");
	  
  /** The kernel matrix. */
  protected Matrix m_KernelMatrix;
	  
  /** A classifier counter. */
  protected int m_Counter;
 
  /**
   * Returns a string describing the kernel
   * 
   * @return a description suitable for displaying in the
   *         explorer/experimenter gui
   */
  public String globalInfo() {
    return 
      "This kernel is based on a static kernel matrix that is read from a file. " +
      "Instances must have a single nominal attribute (excluding the class). " +
      "This attribute must be the first attribute in the file and its values are " +
      "used to reference rows/columns in the kernel matrix. The second attribute " +
      "must be the class attribute.";
  }
	  
  /**
   * Returns an enumeration describing the available options.
   *
   * @return 		an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector		result;
	    
    result = new Vector();

    Enumeration en = super.listOptions();
    while (en.hasMoreElements())
      result.addElement(en.nextElement());

    result.addElement(new Option(
                                 "\tThe file name of the file that holds the kernel matrix.\n"
                                 + "\t(default: kernelMatrix.matrix)",
                                 "M", 1, "-M <file name>"));

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
   * <pre> -M &lt;file name&gt;
   *  The file name of the file that holds the kernel matrix.
   *  (default: kernelMatrix.matrix)</pre>
   * 
   <!-- options-end -->
   * 
   * @param options 	the list of options as an array of strings
   * @throws Exception 	if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String	tmpStr;
	    
    tmpStr = Utils.getOption('M', options);
    if (tmpStr.length() != 0)
      setKernelMatrixFile(new File(tmpStr));
    else
      setKernelMatrixFile(new File("kernelMatrix.matrix"));
	    
    super.setOptions(options);
  }

  /**
   * Gets the current settings of the Kernel.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    int       i;
    Vector<String>    result;
    String[]  options;

    result = new Vector<String>();
	    
    options = super.getOptions();
    for (i = 0; i < options.length; i++)
      result.add(options[i]);

    result.add("-M");
    result.add("" + getKernelMatrixFile());

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
  public double eval(int id1, int id2, Instance inst1)
    throws Exception {
			
    if (m_KernelMatrix == null) {
      throw new IllegalArgumentException("Kernel matrix has not been loaded successfully.");
    }
    int index1 = -1;
    if (id1 > -1) { 
      index1 = (int)m_data.instance(id1).value(0);
    } else {
      index1 = (int)inst1.value(0);
    }
    int index2 = (int)m_data.instance(id2).value(0);
    return m_KernelMatrix.get(index1, index2);
  }
	  
  /**
   * initializes variables etc.
   * 
   * @param data	the data to use
   */
  protected void initVars(Instances data) {
    super.initVars(data);

    try {
      if (m_KernelMatrix == null) {
        m_KernelMatrix = new Matrix(new FileReader(m_KernelMatrixFile));
        //        System.err.println("Read kernel matrix.");
      }
    } catch (Exception e) {
      System.err.println("Problem reading matrix from " + m_KernelMatrixFile);
    }
    m_Counter++;
    //    System.err.print("Building classifier: " + m_Counter + "\r");
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
	    
    result.enable(Capability.NOMINAL_ATTRIBUTES);
    result.enableAllClasses();
    result.enable(Capability.MISSING_CLASS_VALUES);
	    
    return result;
  }
	  
  /**
   * Sets the file holding the kernel matrix
   * 
   * @param f	the file holding the matrix
   */
  public void setKernelMatrixFile(File f) {
    m_KernelMatrixFile = f;
  }
	  
  /**
   * Gets the file containing the kernel matrix.
   * 
   * @return		the exponent value
   */
  public File getKernelMatrixFile() {
    return m_KernelMatrixFile;
  }

  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String kernelMatrixFileTipText() {
    return "The file holding the kernel matrix.";
  }

  /**
   * Set the kernel matrix. This method is used by the
   * unit test for this class, as it loads at test matrix
   * as a system resource.
   *
   * @param km the kernel matrix to use
   */
  protected void setKernelMatrix(Matrix km) {
    m_KernelMatrix = km;
  }
		  
  /**
   * returns a string representation for the Kernel
   * 
   * @return 		a string representaiton of the kernel
   */
  public String toString() {
    return "Using kernel matrix from file with name: " + getKernelMatrixFile();
  }
	  
  /**
   * Frees the memory used by the kernel.
   * (Useful with kernels which use cache.)
   * This function is called when the training is done.
   * i.e. after that, eval will be called with id1 == -1.
   */
  public void clean() {
    // do nothing
  }

  /**
   * Returns the number of kernel evaluation performed.
   *
   * @return the number of kernel evaluation performed.
   */
  public int numEvals() {
    return 0;
  }

  /**
   * Returns the number of dot product cache hits.
   *
   * @return the number of dot product cache hits, or -1 if not supported by this kernel.
   */
  public int numCacheHits() {
    return 0;
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

