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
 *    Normalize.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.filters.unsupervised.attribute;

import weka.core.Capabilities;
import weka.core.Instance; 
import weka.core.DenseInstance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.SparseInstance;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.core.labelranking.PreferenceDenseInstance;
import weka.filters.Sourcable;
import weka.filters.UnsupervisedFilter;

import java.util.Enumeration;
import java.util.Vector;

/** 
 <!-- globalinfo-start -->
 * Normalizes all numeric values in the given dataset (apart from the class attribute, if set). The resulting values are by default in [0,1] for the data used to compute the normalization intervals. But with the scale and translation parameters one can change that, e.g., with scale = 2.0 and translation = -1.0 you get values in the range [-1,+1].
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -unset-class-temporarily
 *  Unsets the class index temporarily before the filter is
 *  applied to the data.
 *  (default: no)</pre>
 * 
 * <pre> -S &lt;num&gt;
 *  The scaling factor for the output range.
 *  (default: 1.0)</pre>
 * 
 * <pre> -T &lt;num&gt;
 *  The translation of the output range.
 *  (default: 0.0)</pre>
 * 
 <!-- options-end -->
 * 
 * @author Eibe Frank (eibe@cs.waikato.ac.nz) 
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5987 $
 */
public class Normalize 
  extends PotentialClassIgnorer 
  implements UnsupervisedFilter, Sourcable, OptionHandler {
  
  /** for serialization. */
  static final long serialVersionUID = -8158531150984362898L;

  /** The minimum values for numeric attributes. */
  protected double[] m_MinArray;
  
  /** The maximum values for numeric attributes. */
  protected double[] m_MaxArray;

  /** The translation of the output range. */
  protected double m_Translation = 0;
  
  /** The scaling factor of the output range. */
  protected double m_Scale = 1.0;
  
  /**
   * Returns a string describing this filter.
   *
   * @return 		a description of the filter suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "Normalizes all numeric values in the given dataset (apart from the "
      + "class attribute, if set). The resulting values are by default "
      + "in [0,1] for the data used to compute the normalization intervals. "
      + "But with the scale and translation parameters one can change that, "
      + "e.g., with scale = 2.0 and translation = -1.0 you get values in the "
      + "range [-1,+1].";
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return 		an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector result = new Vector();

    Enumeration en = super.listOptions();
    while (en.hasMoreElements())
      result.addElement(en.nextElement());

    result.addElement(new Option(
	"\tThe scaling factor for the output range.\n"
	+ "\t(default: 1.0)",
	"S", 1, "-S <num>"));

    result.addElement(new Option(
	"\tThe translation of the output range.\n"
	+"\t(default: 0.0)",
	"T", 1,"-T <num>"));

    return result.elements();
  }


  /**
   * Parses a given list of options. <p/>
   * 
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -unset-class-temporarily
   *  Unsets the class index temporarily before the filter is
   *  applied to the data.
   *  (default: no)</pre>
   * 
   * <pre> -S &lt;num&gt;
   *  The scaling factor for the output range.
   *  (default: 1.0)</pre>
   * 
   * <pre> -T &lt;num&gt;
   *  The translation of the output range.
   *  (default: 0.0)</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String      tmpStr;

    tmpStr = Utils.getOption('S', options);
    if (tmpStr.length() != 0)
      setScale(Double.parseDouble(tmpStr));
    else
      setScale(1.0);
    
    tmpStr = Utils.getOption('T', options);
    if (tmpStr.length() != 0)
      setTranslation(Double.parseDouble(tmpStr));
    else
      setTranslation(0.0);

    if (getInputFormat() != null)
      setInputFormat(getInputFormat());
  }

  /**
   * Gets the current settings of the filter.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    Vector<String>	result;
    
    result = new Vector<String>();

    result.add("-S");
    result.add("" + getScale());

    result.add("-T");
    result.add("" + getTranslation());
    
    return result.toArray(new String[result.size()]);
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
    
    return result;
  }

  /**
   * Sets the format of the input instances.
   *
   * @param instanceInfo 	an Instances object containing the input 
   * 				instance structure (any instances contained in 
   * 				the object are ignored - only the structure is 
   * 				required).
   * @return 			true if the outputFormat may be collected 
   * 				immediately
   * @throws Exception 		if the input format can't be set successfully
   */
  public boolean setInputFormat(Instances instanceInfo) 
       throws Exception {

    super.setInputFormat(instanceInfo);
    setOutputFormat(instanceInfo);
    m_MinArray = m_MaxArray = null;
    return true;
  }

  /**
   * Input an instance for filtering. Filter requires all
   * training instances be read before producing output.
   *
   * @param instance 	the input instance
   * @return 		true if the filtered instance may now be
   * 			collected with output().
   * @throws Exception 	if an error occurs
   * @throws IllegalStateException 	if no input format has been set.
   */
  public boolean input(Instance instance) throws Exception {
    if (getInputFormat() == null)
      throw new IllegalStateException("No input instance format defined");
    
    if (m_NewBatch) {
      resetQueue();
      m_NewBatch = false;
    }
    if (m_MinArray == null) {
      bufferInput(instance);
      return false;
    }
    else {
      convertInstance(instance);
      return true;
    }
  }

  /**
   * Signify that this batch of input to the filter is finished. 
   * If the filter requires all instances prior to filtering,
   * output() may now be called to retrieve the filtered instances.
   *
   * @return 		true if there are instances pending output
   * @throws Exception 	if an error occurs
   * @throws IllegalStateException 	if no input structure has been defined
   */
  public boolean batchFinished() throws Exception {
    if (getInputFormat() == null)
      throw new IllegalStateException("No input instance format defined");

    if (m_MinArray == null) {
      Instances input = getInputFormat();
      // Compute minimums and maximums
      m_MinArray = new double[input.numAttributes()];
      m_MaxArray = new double[input.numAttributes()];
      for (int i = 0; i < input.numAttributes(); i++)
	m_MinArray[i] = Double.NaN;

      for (int j = 0; j < input.numInstances(); j++) {
	double[] value = input.instance(j).toDoubleArray();
	for (int i = 0; i < input.numAttributes(); i++) {
	  if (input.attribute(i).isNumeric() &&
	      (input.classIndex() != i)) {
	    if (!Utils.isMissingValue(value[i])) {
	      if (Double.isNaN(m_MinArray[i])) {
		m_MinArray[i] = m_MaxArray[i] = value[i];
	      }
	      else {
		if (value[i] < m_MinArray[i])
		  m_MinArray[i] = value[i];
		if (value[i] > m_MaxArray[i])
		  m_MaxArray[i] = value[i];
	      }
	    }
	  }
	} 
      }

      // Convert pending input instances
      for (int i = 0; i < input.numInstances(); i++)
	convertInstance(input.instance(i));
    } 
    // Free memory
    flushInput();

    m_NewBatch = true;
    return (numPendingOutput() != 0);
  }

  /**
   * Convert a single instance over. The converted instance is 
   * added to the end of the output queue.
   *
   * @param instance 	the instance to convert
   * @throws Exception 	if conversion fails
   */
  protected void convertInstance(Instance instance) throws Exception {
    Instance inst = null;
    if (instance instanceof SparseInstance) {
      double[] newVals = new double[instance.numAttributes()];
      int[] newIndices = new int[instance.numAttributes()];
      double[] vals = instance.toDoubleArray();
      int ind = 0;
      for (int j = 0; j < instance.numAttributes(); j++) {
	double value;
	if (instance.attribute(j).isNumeric() &&
	    (!Utils.isMissingValue(vals[j])) &&
	    (getInputFormat().classIndex() != j)) {
	  if (Double.isNaN(m_MinArray[j]) ||
	      (m_MaxArray[j] == m_MinArray[j])) {
	    value = 0;
	  }
	  else {
	    value = (vals[j] - m_MinArray[j]) / 
	      (m_MaxArray[j] - m_MinArray[j]) * m_Scale + m_Translation;
            if (Double.isNaN(value)) {
              throw new Exception("A NaN value was generated "
                                  + "while normalizing " 
                                  + instance.attribute(j).name());
            }
	  }
	  if (value != 0.0) {
	    newVals[ind] = value;
	    newIndices[ind] = j;
	    ind++;
	  }
	}
	else {
	  value = vals[j];
	  if (value != 0.0) {
	    newVals[ind] = value;
	    newIndices[ind] = j;
	    ind++;
	  }
	}
      }	
      double[] tempVals = new double[ind];
      int[] tempInd = new int[ind];
      System.arraycopy(newVals, 0, tempVals, 0, ind);
      System.arraycopy(newIndices, 0, tempInd, 0, ind);
      inst = new SparseInstance(instance.weight(), tempVals, tempInd,
                                instance.numAttributes());
    }
    else {
      double[] vals = instance.toDoubleArray();
      for (int j = 0; j < getInputFormat().numAttributes(); j++) {
	if (instance.attribute(j).isNumeric() &&
	    (!Utils.isMissingValue(vals[j])) &&
	    (getInputFormat().classIndex() != j)) {
	  if (Double.isNaN(m_MinArray[j]) ||
	      (m_MaxArray[j] == m_MinArray[j])) {
	    vals[j] = 0;
	  }
	  else {
	    vals[j] = (vals[j] - m_MinArray[j]) / 
	      (m_MaxArray[j] - m_MinArray[j]) * m_Scale + m_Translation;
            if (Double.isNaN(vals[j])) {
              throw new Exception("A NaN value was generated "
                                  + "while normalizing " 
                                  + instance.attribute(j).name());
            }
	  }
	}
      }	
      //RANKING BEGIN
      if(instance instanceof PreferenceDenseInstance){
    	  PreferenceDenseInstance pdi = (PreferenceDenseInstance) instance;
    	  inst = new PreferenceDenseInstance(instance.weight(), vals, pdi.getHashMap());
      }
      else
    	  inst = new DenseInstance(instance.weight(), vals);
      //RANKING END
    }
    inst.setDataset(instance.dataset());
    push(inst);
  }
  
  /**
   * Returns a string that describes the filter as source. The
   * filter will be contained in a class with the given name (there may
   * be auxiliary classes),
   * and will contain two methods with these signatures:
   * <pre><code>
   * // converts one row
   * public static Object[] filter(Object[] i);
   * // converts a full dataset (first dimension is row index)
   * public static Object[][] filter(Object[][] i);
   * </code></pre>
   * where the array <code>i</code> contains elements that are either
   * Double, String, with missing values represented as null. The generated
   * code is public domain and comes with no warranty.
   *
   * @param className   the name that should be given to the source class.
   * @param data	the dataset used for initializing the filter
   * @return            the object source described by a string
   * @throws Exception  if the source can't be computed
   */
  public String toSource(String className, Instances data) throws Exception {
    StringBuffer        result;
    boolean[]		process;
    int			i;
    
    result = new StringBuffer();
    
    // determine what attributes were processed
    process = new boolean[data.numAttributes()];
    for (i = 0; i < data.numAttributes(); i++) 
      process[i] = (data.attribute(i).isNumeric() && (i != data.classIndex()));
  
    result.append("class " + className + " {\n");
    result.append("\n");
    result.append("  /** lists which attributes will be processed */\n");
    result.append("  protected final static boolean[] PROCESS = new boolean[]{" + Utils.arrayToString(process) + "};\n");
    result.append("\n");
    result.append("  /** the minimum values for numeric values */\n");
    result.append("  protected final static double[] MIN = new double[]{" + Utils.arrayToString(m_MinArray).replaceAll("NaN", "Double.NaN") + "};\n");
    result.append("\n");
    result.append("  /** the maximum values for numeric values */\n");
    result.append("  protected final static double[] MAX = new double[]{" + Utils.arrayToString(m_MaxArray) + "};\n");
    result.append("\n");
    result.append("  /** the scale factor */\n");
    result.append("  protected final static double SCALE = " + m_Scale + ";\n");
    result.append("\n");
    result.append("  /** the translation */\n");
    result.append("  protected final static double TRANSLATION = " + m_Translation + ";\n");
    result.append("\n");
    result.append("  /**\n");
    result.append("   * filters a single row\n");
    result.append("   * \n");
    result.append("   * @param i the row to process\n");
    result.append("   * @return the processed row\n");
    result.append("   */\n");
    result.append("  public static Object[] filter(Object[] i) {\n");
    result.append("    Object[] result;\n");
    result.append("\n");
    result.append("    result = new Object[i.length];\n");
    result.append("    for (int n = 0; n < i.length; n++) {\n");
    result.append("      if (PROCESS[n] && (i[n] != null)) {\n");
    result.append("        if (Double.isNaN(MIN[n]) || (MIN[n] == MAX[n]))\n");
    result.append("          result[n] = 0;\n");
    result.append("        else\n");
    result.append("          result[n] = (((Double) i[n]) - MIN[n]) / (MAX[n] - MIN[n]) * SCALE + TRANSLATION;\n");
    result.append("      }\n");
    result.append("      else {\n");
    result.append("        result[n] = i[n];\n");
    result.append("      }\n");
    result.append("    }\n");
    result.append("\n");
    result.append("    return result;\n");
    result.append("  }\n");
    result.append("\n");
    result.append("  /**\n");
    result.append("   * filters multiple rows\n");
    result.append("   * \n");
    result.append("   * @param i the rows to process\n");
    result.append("   * @return the processed rows\n");
    result.append("   */\n");
    result.append("  public static Object[][] filter(Object[][] i) {\n");
    result.append("    Object[][] result;\n");
    result.append("\n");
    result.append("    result = new Object[i.length][];\n");
    result.append("    for (int n = 0; n < i.length; n++) {\n");
    result.append("      result[n] = filter(i[n]);\n");
    result.append("    }\n");
    result.append("\n");
    result.append("    return result;\n");
    result.append("  }\n");
    result.append("}\n");
    
    return result.toString();
  }

  /**
   * Returns the calculated minimum values for the attributes in the data.
   * 
   * @return		the array with the minimum values
   */
  public double[] getMinArray() {
    return m_MinArray;
  }

  /**
   * Returns the calculated maximum values for the attributes in the data.
   * 
   * @return		the array with the maximum values
   */
  public double[] getMaxArray() {
    return m_MaxArray;
  }

  /**
   * Returns the tip text for this property.
   *
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String scaleTipText() {
    return "The factor for scaling the output range (default: 1).";
  }

  /**
   * Get the scaling factor.
   *
   * @return 		the factor
   */
  public double getScale() {
    return m_Scale;
  }

  /**
   * Sets the scaling factor.
   *
   * @param value 	the scaling factor
   */
  public void setScale(double value) {
    m_Scale = value;
  }

  /**
   * Returns the tip text for this property.
   *
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String translationTipText() {
    return "The translation of the output range (default: 0).";
  }

  /**
   * Get the translation.
   *
   * @return 		the translation
   */
  public double getTranslation() {
    return m_Translation;
  }

  /**
   * Sets the translation.
   *
   * @param value 	the translation
   */
  public void setTranslation(double value) {
    m_Translation = value;
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
   * Main method for running this filter.
   *
   * @param args 	should contain arguments to the filter, use -h for help
   */
  public static void main(String[] args) {
    runFilter(new Normalize(), args);
  }
}
