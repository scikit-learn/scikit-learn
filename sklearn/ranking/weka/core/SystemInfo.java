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
 *    SystemInfo.java
 *    Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */
package weka.core;

import weka.gui.LookAndFeel;

import java.util.Collections;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Properties;
import java.util.Vector;

/**
 * This class prints some information about the system setup, like Java
 * version, JVM settings etc. Useful for Bug-Reports.
 *
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5953 $
 */
public class SystemInfo
  implements RevisionHandler {
  
  /** for storing the information */
  private Hashtable<String,String> m_Info = null;
  
  /**
   * initializes the object and reads the system information
   */
  public SystemInfo() {
    m_Info = new Hashtable<String,String>();
    readProperties();
  }

  /**
   * reads all the properties and stores them in the hashtable
   */
  private void readProperties() {
    Properties          props;
    Enumeration         enm;
    String              name;
    String[]            laf;
    String              tmpStr;
    int                 i;
    Memory              mem;
    
    m_Info.clear();

    // System information
    props = System.getProperties();
    enm   = props.propertyNames();
    while (enm.hasMoreElements()) {
      name = (String)enm.nextElement();
      m_Info.put(name, (String)props.get(name));
    }

    // additional WEKA info
    m_Info.put("weka.version", Version.VERSION);

    // look and feel info
    laf    = LookAndFeel.getInstalledLookAndFeels();
    tmpStr = "";
    for (i = 0; i < laf.length; i++) {
      if (i > 0)
        tmpStr += ",";
      tmpStr += laf[i];
    }
    m_Info.put("ui.installedLookAndFeels", tmpStr);
    m_Info.put("ui.currentLookAndFeel", LookAndFeel.getSystemLookAndFeel());

    // memory info
    mem = new Memory();
    m_Info.put(
        "memory.initial", 
        "" + Utils.doubleToString(Memory.toMegaByte(mem.getInitial()), 1) + "MB" 
        + " (" + mem.getInitial() + ")");
    m_Info.put(
        "memory.max", 
        "" + Utils.doubleToString(Memory.toMegaByte(mem.getMax()), 1) + "MB"
        + " (" + mem.getMax() + ")");
  }

  /**
   * returns a copy of the system info. the key is the name of the property
   * and the associated object is the value of the property (a string).
   */
  public Hashtable getSystemInfo() {
    return (Hashtable) m_Info.clone();
  }

  /**
   * returns a string representation of all the system properties
   */
  public String toString() {
    Enumeration<String>     enm;
    String          result;
    String          key;
    Vector<String>          keys;
    int             i;
    String          value;

    result = "";
    keys   = new Vector<String>();
    
    // get names and sort them
    enm = m_Info.keys();
    while (enm.hasMoreElements())
      keys.add(enm.nextElement());
    Collections.sort(keys);
    
    // generate result
    for (i = 0; i < keys.size(); i++) {
      key   = keys.get(i).toString();
      value = m_Info.get(key).toString();
      if (key.equals("line.separator"))
        value = Utils.backQuoteChars(value);
      result += key + ": " + value + "\n";
    }

    return result;
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
   * for printing the system info to stdout.
   */
  public static void main(String[] args) {
    System.out.println(new SystemInfo());
  }
}
