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
 *    SingleAssociatorEnhancer.java
 *    Copyright (C) 2007 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.associations;

import weka.core.Capabilities;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.Utils;
import weka.core.Capabilities.Capability;

import java.util.Enumeration;
import java.util.Vector;

/**
 * Abstract utility class for handling settings common to meta
 * associators that use a single base associator.  
 *
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 1.2 $
 */
public abstract class SingleAssociatorEnhancer
  extends AbstractAssociator
  implements OptionHandler {

  /** for serialization */
  private static final long serialVersionUID = -3665885256363525164L;

  /** The base associator to use */
  protected Associator m_Associator = new Apriori();

  /**
   * String describing default Associator.
   * 
   * @return		default classname
   */
  protected String defaultAssociatorString() {
    return Apriori.class.getName();
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return 		an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector result = new Vector();

    result.addElement(new Option(
	      "\tFull name of base associator.\n"
	      + "\t(default: " + defaultAssociatorString() +")",
	      "W", 1, "-W"));

    if (m_Associator instanceof OptionHandler) {
      result.addElement(new Option(
	  "",
	  "", 0, "\nOptions specific to associator "
	  + m_Associator.getClass().getName() + ":"));

      Enumeration enm = ((OptionHandler) m_Associator).listOptions();
      while (enm.hasMoreElements())
	result.addElement(enm.nextElement());
    }

    return result.elements();
  }

  /**
   * Parses a given list of options. Valid options are:<p>
   *
   * -W classname <br>
   * Specify the full class name of the base associator.<p>
   *
   * Options after -- are passed to the designated associator.<p>
   *
   * @param options 	the list of options as an array of strings
   * @throws Exception 	if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String	tmpStr;
    
    tmpStr = Utils.getOption('W', options);
    if (tmpStr.length() > 0) { 
      // This is just to set the associator in case the option 
      // parsing fails.
      setAssociator(AbstractAssociator.forName(tmpStr, null));
      setAssociator(AbstractAssociator.forName(tmpStr, Utils.partitionOptions(options)));
    }
    else {
      // This is just to set the associator in case the option 
      // parsing fails.
      setAssociator(AbstractAssociator.forName(defaultAssociatorString(), null));
      setAssociator(AbstractAssociator.forName(defaultAssociatorString(), Utils.partitionOptions(options)));
    }
  }

  /**
   * Gets the current settings of the associator.
   *
   * @return 		an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    int       		i;
    Vector<String>    	result;
    String[]		options;
    
    result = new Vector<String>();

    result.add("-W");
    result.add(getAssociator().getClass().getName());

    if (getAssociator() instanceof OptionHandler) {
      options = ((OptionHandler) getAssociator()).getOptions();
      result.add("--");
      for (i = 0; i < options.length; i++)
	result.add(options[i]);
    }

    return result.toArray(new String[result.size()]);
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String associatorTipText() {
    return "The base associator to be used.";
  }

  /**
   * Set the base associator.
   *
   * @param value 	the associator to use.
   */
  public void setAssociator(Associator value) {
    m_Associator = value;
  }

  /**
   * Get the associator used as the base associator.
   *
   * @return 		the currently used associator
   */
  public Associator getAssociator() {
    return m_Associator;
  }
  
  /**
   * Gets the associator specification string, which contains the class name of
   * the associator and any options to the associator
   *
   * @return the associator string
   */
  protected String getAssociatorSpec() {
    Associator c = getAssociator();
    return c.getClass().getName() + " "
      + Utils.joinOptions(((OptionHandler)c).getOptions());
  }

  /**
   * Returns default capabilities of the base associator.
   *
   * @return      the capabilities of the base associator
   */
  public Capabilities getCapabilities() {
    Capabilities        result;

    if (getAssociator() != null)
      result = getAssociator().getCapabilities();
    else
      result = new Capabilities(this);
    
    // set dependencies
    for (Capability cap: Capability.values())
      result.enableDependency(cap);
    
    result.setOwner(this);
    
    return result;
  }
}
