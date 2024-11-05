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
 *    CheckGOE.java
 *    Copyright (C) 2007 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core;

import java.beans.BeanInfo;
import java.beans.Introspector;
import java.beans.PropertyDescriptor;
import java.util.Collections;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Vector;

/**
 * Simple command line checking of classes that are editable in the GOE.<p/>
 *
 * Usage: <p/>
 * <code>
 *     CheckGOE -W classname -- test options
 * </code> <p/>
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -D
 *  Turn on debugging output.</pre>
 * 
 * <pre> -S
 *  Silent mode - prints nothing to stdout.</pre>
 * 
 * <pre> -ignored &lt;comma-separated list of properties&gt;
 *  Skipped properties.
 *  (default: capabilities,options)</pre>
 * 
 * <pre> -W
 *  Full name of the class analysed.
 *  eg: weka.classifiers.rules.ZeroR
 *  (default weka.classifiers.rules.ZeroR)</pre>
 * 
 <!-- options-end -->
 *
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5953 $
 */
public class CheckGOE
  extends Check {

  /** the object to test */
  protected Object m_Object = new weka.classifiers.rules.ZeroR();
  
  /** whether the tests were successful */
  protected boolean m_Success;
  
  /** properties that are skipped in the checkToolTips method 
   * @see #checkToolTips() */
  protected HashSet<String> m_IgnoredProperties = new HashSet<String>();
  
  /**
   * default constructor
   */
  public CheckGOE() {
    super();
    
    // set default options
    try {
      setOptions(new String[0]);
    }
    catch (Exception e) {
      e.printStackTrace();
    }
  }
  
  /**
   * Returns an enumeration describing the available options.
   *
   * @return 		an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector<Option> result = new Vector<Option>();
    
    Enumeration en = super.listOptions();
    while (en.hasMoreElements())
      result.addElement((Option)en.nextElement());
    
    result.addElement(new Option(
        "\tSkipped properties.\n"
	+ "\t(default: capabilities,options)",
        "ignored", 1, "-ignored <comma-separated list of properties>"));
    
    result.addElement(new Option(
        "\tFull name of the class analysed.\n"
        +"\teg: weka.classifiers.rules.ZeroR\n"
        + "\t(default weka.classifiers.rules.ZeroR)",
        "W", 1, "-W"));
    
    return result.elements();
  }
  
  /**
   * Parses a given list of options. <p/>
   *
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -D
   *  Turn on debugging output.</pre>
   * 
   * <pre> -S
   *  Silent mode - prints nothing to stdout.</pre>
   * 
   * <pre> -ignored &lt;comma-separated list of properties&gt;
   *  Skipped properties.
   *  (default: capabilities,options)</pre>
   * 
   * <pre> -W
   *  Full name of the class analysed.
   *  eg: weka.classifiers.rules.ZeroR
   *  (default weka.classifiers.rules.ZeroR)</pre>
   * 
   <!-- options-end -->
   *
   * @param options 	the list of options as an array of strings
   * @throws Exception 	if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String      tmpStr;
    
    super.setOptions(options);
    
    tmpStr = Utils.getOption('W', options);
    if (tmpStr.length() == 0)
      tmpStr = weka.classifiers.rules.ZeroR.class.getName();
    setObject(Utils.forName(Object.class, tmpStr, null));
    
    tmpStr = Utils.getOption("ignored", options);
    if (tmpStr.length() == 0)
      tmpStr = "capabilities,options";
    setIgnoredProperties(tmpStr);
  }
  
  /**
   * Gets the current settings of the object.
   *
   * @return 		an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    Vector<String>        result;
    String[]      options;
    int           i;
    
    result = new Vector<String>();
    
    options = super.getOptions();
    for (i = 0; i < options.length; i++)
      result.add(options[i]);

    result.add("-ignored");
    result.add(getIgnoredProperties());
    
    if (getObject() != null) {
      result.add("-W");
      result.add(getObject().getClass().getName());
    }
    
    return (String[]) result.toArray(new String[result.size()]);
  }
  
  /**
   * Set the object to work on.. 
   *
   * @param value 	the object to use.
   */
  public void setObject(Object value) {
    m_Object = value;
  }
  
  /**
   * Get the object used in the tests.
   *
   * @return 		the object used in the tests.
   */
  public Object getObject() {
    return m_Object;
  }

  /**
   * Sets the properties to ignore in checkToolTips(). Comma-separated list.
   * 
   * @param value	the list of properties
   * @see 		#checkToolTips()
   */
  public void setIgnoredProperties(String value) {
    String[] 	props;
    int		i;
    
    m_IgnoredProperties.clear();
    props = value.split(",");
    for (i = 0; i < props.length; i++)
      m_IgnoredProperties.add(props[i]);
  }

  /**
   * Get the ignored properties used in checkToolTips() as comma-separated 
   * list (sorted).
   * 
   * @return		the ignored properties
   * @see		#checkToolTips()
   */
  public String getIgnoredProperties() {
    String		result;
    Vector<String>	list;
    Iterator		iter;
    int			i;
    
    list = new Vector<String>();
    iter = m_IgnoredProperties.iterator();
    while (iter.hasNext())
      list.add((String) iter.next());
    
    // sort
    if (list.size() > 1)
      Collections.sort(list);
    
    result = "";
    for (i = 0; i < list.size(); i++) {
      if (i > 0)
	result += ",";
      result += list.get(i);
    }
     
    return result;
  }
  
  /**
   * returns the success of the tests
   * 
   * @return		true if the tests were successful
   */
  public boolean getSuccess() {
    return m_Success;
  }

  /**
   * checks whether the object declares a globalInfo method.
   * 
   * @return 		true if the test was passed
   */
  public boolean checkGlobalInfo() {
    boolean 	result;
    Class<?>	cls;
    
    print("Global info...");
    
    result = true;
    cls    = getObject().getClass();
    
    // test for globalInfo method
    try {
      cls.getMethod("globalInfo", (Class[]) null);
    }
    catch (Exception e) {
      result = false;
    }

    if (result)
      println("yes");
    else
      println("no");

    return result;
  }

  /**
   * checks whether the object declares tip text method for all its
   * properties.
   * 
   * @return 		true if the test was passed
   */
  public boolean checkToolTips() {
    boolean 			result;
    Class<?>			cls;
    BeanInfo			info;
    PropertyDescriptor[]	desc;
    int				i;
    Vector<String>		missing;
    String			suffix;
    
    print("Tool tips...");
    
    result = true;
    suffix = "TipText";
    cls    = getObject().getClass();
    
    // get properties
    try {
      info = Introspector.getBeanInfo(cls, Object.class);
      desc = info.getPropertyDescriptors();
    }
    catch (Exception e) {
      e.printStackTrace();
      desc = null;
    }

    // test for TipText methods
    if (desc != null) {
      missing = new Vector<String>();
      
      for (i = 0; i < desc.length; i++) {
	// skip property?
	if (m_IgnoredProperties.contains(desc[i].getName()))
	  continue;
	if ((desc[i].getReadMethod() == null) || (desc[i].getWriteMethod() == null))
	  continue;
	  
	try {
	  cls.getMethod(desc[i].getName() + suffix, (Class[]) null);
	}
	catch (Exception e) {
	  result = false;
	  missing.add(desc[i].getName() + suffix);
	}
      }
      
      if (result)
	println("yes");
      else
	println("no (missing: " + missing + ")");

    }
    else {
      println("maybe");
    }
    
    return result;
  }
  
  /**
   * Runs some diagnostic tests on the object. Output is
   * printed to System.out (if not silent).
   */
  public void doTests() {
    println("Object: " + m_Object.getClass().getName() + "\n");

    println("--> Tests");

    m_Success = checkGlobalInfo();

    if (m_Success)
      m_Success = checkToolTips();
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
   * Main method for using the CheckGOE.
   *
   * @param args 	the options to the CheckGOE
   */
  public static void main(String[] args) {
    CheckGOE check = new CheckGOE();
    runCheck(check, args);
    if (check.getSuccess())
      System.exit(0);
    else
      System.exit(1);
  }
}
