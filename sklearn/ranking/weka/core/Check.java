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
 * CheckScheme.java
 * Copyright (C) 2006 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core;

import java.util.Enumeration;
import java.util.Vector;

/**
 * Abstract general class for testing in Weka.
 *
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5953 $
 */
public abstract class Check
  implements OptionHandler, RevisionHandler {
  
  /** Debugging mode, gives extra output if true */
  protected boolean m_Debug = false;
  
  /** Silent mode, for no output at all to stdout */
  protected boolean m_Silent = false;
  
  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector<Option> result = new Vector<Option>();
    
    result.addElement(new Option(
        "\tTurn on debugging output.",
        "D", 0, "-D"));
    
    result.addElement(new Option(
        "\tSilent mode - prints nothing to stdout.",
        "S", 0, "-S"));
    
    return result.elements();
  }
  
  /**
   * Parses a given list of options. 
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    setDebug(Utils.getFlag('D', options));

    setSilent(Utils.getFlag('S', options));
  }
  
  /**
   * Gets the current settings of the CheckClassifier.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    Vector<String>        result;
    
    result = new Vector<String>();
    
    if (getDebug())
      result.add("-D");
    
    if (getSilent())
      result.add("-S");
    
    return (String[]) result.toArray(new String[result.size()]);
  }
  
  /**
   * Tries to instantiate a new instance of the given class and checks whether
   * it is an instance of the specified class. For convenience one can also
   * specify a classname prefix (e.g., "weka.classifiers") to avoid long 
   * classnames and then instantiate it with the shortened classname (e.g.,
   * "trees.J48").
   * 
   * @param prefix	the classname prefix (without trailing dot)
   * @param cls		the class to check whether the generated object is an
   * 			instance of
   * @param classname	the classname to instantiate
   * @param options	optional options for the object
   * @return		the configured object
   * @throws Exception	if instantiation fails
   */
  protected Object forName(String prefix, Class cls, String classname, 
      String[] options) throws Exception {

    Object	result;
    
    result = null;

    try {
      result = Utils.forName(cls, classname, options);
    }
    catch (Exception e) {
      // shall we try with prefix?
      if (e.getMessage().toLowerCase().indexOf("can't find") > -1) {
	try {
	  result = Utils.forName(cls, prefix + "." + classname, options);
	}
	catch (Exception ex) {
	  if (e.getMessage().toLowerCase().indexOf("can't find") > -1) {
	    throw new Exception(
		"Can't find class called '" + classname 
		+ "' or '" + prefix + "." + classname + "'!");
	  }
	  else {
	    throw new Exception(ex);
	  }
	}
      }
      else {
	throw new Exception(e);
      }
    }
    
    return result;
  }
  
  /**
   * Begin the tests, reporting results to System.out
   */
  public abstract void doTests();
  
  /**
   * Set debugging mode
   *
   * @param debug true if debug output should be printed
   */
  public void setDebug(boolean debug) {
    m_Debug = debug;
    // disable silent mode, if necessary
    if (getDebug())
      setSilent(false);
  }
  
  /**
   * Get whether debugging is turned on
   *
   * @return true if debugging output is on
   */
  public boolean getDebug() {
    return m_Debug;
  }
  
  /**
   * Set slient mode, i.e., no output at all to stdout
   *
   * @param value whether silent mode is active or not
   */
  public void setSilent(boolean value) {
    m_Silent = value;
  }
  
  /**
   * Get whether silent mode is turned on
   *
   * @return true if silent mode is on
   */
  public boolean getSilent() {
    return m_Silent;
  }
  
  /**
   * prints the given message to stdout, if not silent mode
   * 
   * @param msg         the text to print to stdout
   */
  protected void print(Object msg) {
    if (!getSilent())
      System.out.print(msg);
  }
  
  /**
   * prints the given message (+ LF) to stdout, if not silent mode
   * 
   * @param msg         the message to println to stdout
   */
  protected void println(Object msg) {
    print(msg + "\n");
  }
  
  /**
   * prints a LF to stdout, if not silent mode
   */
  protected void println() {
    print("\n");
  }
  
  /**
   * runs the CheckScheme with the given options
   * 
   * @param check	the checkscheme to setup and run
   * @param options	the commandline parameters to use
   */
  protected static void runCheck(Check check, String[] options) {
    try {
      try {
        check.setOptions(options);
        Utils.checkForRemainingOptions(options);
      }
      catch (Exception ex) {
        String result = ex.getMessage() + "\n\n" + check.getClass().getName().replaceAll(".*\\.", "") + " Options:\n\n";
        Enumeration enm = check.listOptions();
        while (enm.hasMoreElements()) {
          Option option = (Option) enm.nextElement();
          result += option.synopsis() + "\n" + option.description() + "\n";
        }
        throw new Exception(result);
      }
      
      check.doTests();
    }
    catch (Exception ex) {
      System.err.println(ex.getMessage());
    }
  }
}
