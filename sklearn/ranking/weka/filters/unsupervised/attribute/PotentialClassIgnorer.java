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
 *    PotentialClassIgnorer.java
 *    Copyright (C) 2003 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.filters.unsupervised.attribute;

import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.Utils;
import weka.filters.Filter;

import java.util.Enumeration;
import java.util.Vector;

/**
 * This filter should be extended by other unsupervised attribute
 * filters to allow processing of the class attribute if that's
 * required. It the class is to be ignored it is essential that the
 * extending filter does not change the position (i.e. index) of the
 * attribute that is originally the class attribute !
 *
 * @author Eibe Frank (eibe@cs.waikato.ac.nz), Mark Hall (mhall@cs.waikato.ac.nz)
 * @version $Revision: 1.7 $ 
 */
public abstract class PotentialClassIgnorer
  extends Filter
  implements OptionHandler {

  /** for serialization */
  private static final long serialVersionUID = 8625371119276845454L;

  /** True if the class is to be unset */
  protected boolean m_IgnoreClass = false;

  /** Storing the class index */
  protected int m_ClassIndex = -1;

  /**
   * Returns an enumeration describing the available options.
   *
   * @return 		an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector result = new Vector();
      
    result.addElement(new Option(
	"\tUnsets the class index temporarily before the filter is\n"
	+ "\tapplied to the data.\n"
	+ "\t(default: no)",
	"unset-class-temporarily", 1, "-unset-class-temporarily"));

    return result.elements();
  }

  /**
   * Parses a list of options for this object.
   *
   * @param options 	the list of options as an array of strings
   * @throws Exception 	if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    setIgnoreClass(Utils.getFlag("unset-class-temporarily", options));
  }

  /**
   * Gets the current settings of the filter.
   *
   * @return 		an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    Vector        result;

    result = new Vector();

    if (getIgnoreClass())
      result.add("-unset-class-temporarily");

    return (String[]) result.toArray(new String[result.size()]);
  }

  /**
   * Sets the format of the input instances. If the filter is able to
   * determine the output format before seeing any input instances, it
   * does so here. This default implementation clears the output format
   * and output queue, and the new batch flag is set. Overriders should
   * call <code>super.setInputFormat(Instances)</code>
   *
   * @param instanceInfo an Instances object containing the input instance
   * structure (any instances contained in the object are ignored - only the
   * structure is required).
   * @return true if the outputFormat may be collected immediately
   * @throws Exception if the inputFormat can't be set successfully 
   */
  public boolean setInputFormat(Instances instanceInfo) throws Exception {

    boolean result = super.setInputFormat(instanceInfo);
    if (m_IgnoreClass) {
      m_ClassIndex = inputFormatPeek().classIndex();
      inputFormatPeek().setClassIndex(-1);
    }      
    return result;
  }

  /**
   * Gets the format of the output instances. This should only be called
   * after input() or batchFinished() has returned true. The relation
   * name of the output instances should be changed to reflect the
   * action of the filter (eg: add the filter name and options).
   *
   * @return an Instances object containing the output instance
   * structure only.
   * @throws NullPointerException if no input structure has been
   * defined (or the output format hasn't been determined yet) 
   */
  public Instances getOutputFormat() {

    if (m_IgnoreClass) {
      outputFormatPeek().setClassIndex(m_ClassIndex);
    }
    return super.getOutputFormat();
  }

  /**
   * Returns the tip text for this property
   *
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String ignoreClassTipText() {
    return "The class index will be unset temporarily before the filter is applied.";
  }

  /**
   * Set the IgnoreClass value. Set this to true if the
   * class index is to be unset before the filter is applied.
   * 
   * @param newIgnoreClass The new IgnoreClass value.
   */
  public void setIgnoreClass(boolean newIgnoreClass) {
    m_IgnoreClass = newIgnoreClass;
  }
  
  /**
   * Gets the IgnoreClass value. If this to true then the
   * class index is to unset before the filter is applied.
   * 
   * @return the current IgnoreClass value.
   */
  public boolean getIgnoreClass() {
    return m_IgnoreClass;
  }
}
