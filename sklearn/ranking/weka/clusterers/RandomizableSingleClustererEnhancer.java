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
 * RandomizableSingleClustererEnhancer.java
 * Copyright (C) 2006 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.clusterers;

import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.Randomizable;
import weka.core.Utils;

import java.util.Enumeration;
import java.util.Vector;

/**
 * Abstract utility class for handling settings common to randomizable
 * clusterers.
 *
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 1.3 $
 */
public abstract class RandomizableSingleClustererEnhancer
  extends AbstractClusterer
  implements OptionHandler, Randomizable {

  /** for serialization */
  private static final long serialVersionUID = -644847037106316249L;
  
  /** the default seed value */
  protected int m_SeedDefault = 1;
  
  /** The random number seed. */
  protected int m_Seed = m_SeedDefault;

  /**
   * Returns an enumeration describing the available options.
   *
   * @return 		an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector result = new Vector();

    result.addElement(new Option(
	"\tRandom number seed.\n"
	+ "\t(default " + m_SeedDefault + ")",
	"S", 1, "-S <num>"));

    return result.elements();
  }

  /**
   * Parses a given list of options. Valid options are:<p>
   *
   * @param options 	the list of options as an array of strings
   * @throws Exception 	if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String	tmpStr;
    
    tmpStr = Utils.getOption('S', options);
    if (tmpStr.length() != 0)
      setSeed(Integer.parseInt(tmpStr));
    else
      setSeed(m_SeedDefault);
  }

  /**
   * Gets the current settings of the classifier.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    Vector	result;
    
    result = new Vector();

    result.add("-S");
    result.add("" + getSeed());

    return (String[]) result.toArray(new String[result.size()]);
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String seedTipText() {
    return "The random number seed to be used.";
  }

  /**
   * Set the seed for random number generation.
   *
   * @param value 	the seed to use
   */
  public void setSeed(int value) {
    m_Seed = value;
  }

  /**
   * Gets the seed for the random number generations
   *
   * @return 		the seed for the random number generation
   */
  public int getSeed() {
    return m_Seed;
  }
}
