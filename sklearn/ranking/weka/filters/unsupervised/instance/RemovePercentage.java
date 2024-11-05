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
 *    RemovePercentage.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.filters.unsupervised.instance;

import weka.core.Capabilities;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.filters.Filter;
import weka.filters.UnsupervisedFilter;

import java.util.Enumeration;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * A filter that removes a given percentage of a dataset.
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -P &lt;percentage&gt;
 *  Specifies percentage of instances to select. (default 50)
 * </pre>
 * 
 * <pre> -V
 *  Specifies if inverse of selection is to be output.
 * </pre>
 * 
 <!-- options-end -->
 *
 * @author Richard Kirkby (eibe@cs.waikato.ac.nz)
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version $Revision: 5499 $ 
*/
public class RemovePercentage 
  extends Filter
  implements UnsupervisedFilter, OptionHandler {

  /** for serialization */
  static final long serialVersionUID = 2150341191158533133L;
  
  /** Percentage of instances to select. */
  private double m_Percentage = 50;

  /** Indicates if inverse of selection is to be output. */
  private boolean m_Inverse = false;

  /**
   * Gets an enumeration describing the available options..
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {

    Vector newVector = new Vector(2);

    newVector.addElement(new Option(
              "\tSpecifies percentage of instances to select. (default 50)\n",
              "P", 1, "-P <percentage>"));

    newVector.addElement(new Option(
	      "\tSpecifies if inverse of selection is to be output.\n",
	      "V", 0, "-V"));

    return newVector.elements();
  }

  /**
   * Parses a given list of options. <p/>
   * 
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -P &lt;percentage&gt;
   *  Specifies percentage of instances to select. (default 50)
   * </pre>
   * 
   * <pre> -V
   *  Specifies if inverse of selection is to be output.
   * </pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {

    String percent = Utils.getOption('P', options);
    if (percent.length() != 0) {
      setPercentage(Double.parseDouble(percent));
    } else {
      setPercentage(50.0);
    }
    setInvertSelection(Utils.getFlag('V', options));

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

    String [] options = new String [5];
    int current = 0;

    options[current++] = "-P"; options[current++] = "" + getPercentage();
    if (getInvertSelection()) {
      options[current++] = "-V";
    }

    while (current < options.length) {
      options[current++] = "";
    }
    return options;
  }

  /**
   * Returns a string describing this filter
   *
   * @return a description of the filter suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {

    return "A filter that removes a given percentage of a dataset.";
  }

  /**
   * Returns the tip text for this property
   *
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String percentageTipText() {

    return "The percentage of the data to select.";
  }

  /**
   * Gets the percentage of instances to select.
   * 
   * @return the percentage.
   */
  public double getPercentage() {

    return m_Percentage;
  }

  /**
   * Sets the percentage of intances to select.
   *
   * @param percent the percentage
   * @throws IllegalArgumentException if percentage out of range
   */
  public void setPercentage(double percent) {

    if (percent < 0 || percent > 100) {
      throw new IllegalArgumentException("Percentage must be between 0 and 100.");
    }
    m_Percentage = percent;
  }

  /**
   * Returns the tip text for this property
   *
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String invertSelectionTipText() {

    return "Whether to invert the selection.";
  }

  /**
   * Gets if selection is to be inverted.
   *
   * @return true if the selection is to be inverted
   */
  public boolean getInvertSelection() {

    return m_Inverse;
  }

  /**
   * Sets if selection is to be inverted.
   *
   * @param inverse true if inversion is to be performed
   */
  public void setInvertSelection(boolean inverse) {
    
    m_Inverse = inverse;
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
   * @param instanceInfo an Instances object containing the input instance
   * structure (any instances contained in the object are ignored - only the
   * structure is required).
   * @return true because outputFormat can be collected immediately
   * @throws Exception if the input format can't be set successfully
   */  
  public boolean setInputFormat(Instances instanceInfo) throws Exception {

    super.setInputFormat(instanceInfo);
    setOutputFormat(instanceInfo);
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
   * @throws IllegalStateException if no input format has been set.
   */
  public boolean input(Instance instance) {
     if (getInputFormat() == null) {
        throw new IllegalStateException("No input instance format defined");
     }
     
     if (m_NewBatch) {
        resetQueue();
        m_NewBatch = false;
     }

     if (isFirstBatchDone()) {
       push(instance);
       return true;
     } 
     else {
       bufferInput(instance);
       return false;
     }
  }

  /**
   * Signify that this batch of input to the filter is
   * finished. Output() may now be called to retrieve the filtered
   * instances.
   *
   * @return true if there are instances pending output
   * @throws IllegalStateException if no input structure has been defined 
   */
  public boolean batchFinished() {

    if (getInputFormat() == null) {
      throw new IllegalStateException("No input instance format defined");
    }

    // Push instances for output into output queue
    Instances toFilter = getInputFormat();
    int cutOff = (int) Math.round(toFilter.numInstances() * m_Percentage / 100);
    
    if (m_Inverse) {
      for (int i = 0; i < cutOff; i++) {
	push(toFilter.instance(i));
      }
    } else {
      for (int i = cutOff; i < toFilter.numInstances(); i++) {
	push(toFilter.instance(i));
      }
    }
    flushInput();
    
    m_NewBatch = true;
    m_FirstBatchDone = true;
    
    return (numPendingOutput() != 0);
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5499 $");
  }

  /**
   * Main method for testing this class.
   *
   * @param argv should contain arguments to the filter: use -h for help
   */
  public static void main(String [] argv) {
    runFilter(new RemovePercentage(), argv);
  }
}
