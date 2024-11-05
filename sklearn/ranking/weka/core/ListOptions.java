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
 * ListOptions.java
 * Copyright (C) 2006 University of Waikato, Hamilton, New Zealand
 */

package weka.core;

import java.util.Enumeration;
import java.util.Vector;

/**
 * Lists the options of an OptionHandler
 * 
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5953 $
 */
public class ListOptions
  implements OptionHandler, RevisionHandler {
  
  /** the classname */
  protected String m_Classname = ListOptions.class.getName();
  
  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector<Option> result = new Vector<Option>();

    result.addElement(new Option(
        "\tThe class to load.",
        "W", 1, "-W <classname>"));
    
    return result.elements();
  }
  
  /**
   * Parses a given list of options. 
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String      		tmpStr;
    
    tmpStr = Utils.getOption('W', options);
    if (tmpStr.length() > 0)
      setClassname(tmpStr);
    else
      setClassname(this.getClass().getName());
  }
  
  /**
   * Gets the current settings of this object.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    Vector<String> 	result;

    result = new Vector<String>();
    
    result.add("-W");
    result.add(getClassname());
    
    return (String[]) result.toArray(new String[result.size()]);
  }
  
  /**
   * sets the classname of the class to generate the Javadoc for
   * 
   * @param value	the new classname
   */
  public void setClassname(String value) {
    m_Classname = value;
  }
  
  /**
   * returns the current classname
   * 
   * @return	the current classname
   */
  public String getClassname() {
    return m_Classname;
  }
  
  /**
   * generates a string to print as help on the console
   * 
   * @return 	the generated help
   */
  public String generateHelp() {
    String 		result;
    Enumeration 	enm;
    Option 		option;
    
    result = getClass().getName().replaceAll(".*\\.", "") + " Options:\n\n";
    enm = listOptions();
    while (enm.hasMoreElements()) {
      option = (Option) enm.nextElement();
      result += option.synopsis() + "\n" + option.description() + "\n";
    }
    
    return result;
  }
  
  /**
   * generates the options string.
   * 
   * @return 		the options string
   * @throws Exception 	in case the generation fails
   */
  public String generate() throws Exception {
    StringBuffer 	result;
    OptionHandler	handler;
    Enumeration 	enm;
    Option 		option;
    
    result = new StringBuffer();
    
    handler = (OptionHandler) Class.forName(getClassname()).newInstance();
    
    enm = ((OptionHandler) handler).listOptions();
    while (enm.hasMoreElements()) {
      option = (Option) enm.nextElement();
      result.append(option.synopsis() + '\n');
      result.append(option.description() + "\n");
    }
    
    return result.toString();
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5953 $");
  }
  
  /**
   * runs the javadoc producer with the given commandline options
   * 
   * @param options	the commandline options
   */
  public static void main(String[] options) {
    ListOptions list = new ListOptions();
    
    try {
      try {
	if (Utils.getFlag('h', options))
	  throw new Exception("Help requested");

	list.setOptions(options);
        Utils.checkForRemainingOptions(options);
      } 
      catch (Exception ex) {
        String result = "\n" + ex.getMessage() + "\n\n" + list.generateHelp();
        throw new Exception(result);
      }

      System.out.println("\n" + list.generate());
    } 
    catch (Exception ex) {
      System.err.println(ex.getMessage());
    }
  }
}
