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
 * Copyright.java
 * Copyright (C) 2007 University of Waikato, Hamilton, New Zealand
 */

package weka.core;

import java.util.Calendar;
import java.util.Properties;

/**
 * A class for providing centralized Copyright information.
 * 
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5953 $
 */
public class Copyright {
  
  /** the copyright file */
  public final static String PROPERTY_FILE = "weka/core/Copyright.props";
  
  /** Contains the properties */
  protected static Properties PROPERTIES;

  static {
    PROPERTIES = new Properties();
    
    try {
      //      PROPERTIES.load(ClassLoader.getSystemResourceAsStream(PROPERTY_FILE));
      PROPERTIES.
        load((new Copyright()).getClass().getClassLoader().getResourceAsStream(PROPERTY_FILE));
    }
    catch (Exception e) {
      System.err.println(
	  "Could not read configuration file for the copyright "
	  + "information - using default.");
    }
  }
  
  /**
   * returns the start year of the copyright
   * 
   * @return		the start year
   */
  public static String getFromYear() {
    return PROPERTIES.getProperty("FromYear", "1999");
  }
  
  /**
   * returns the end year of the copyright (i.e., current year)
   * 
   * @return		the end/current year
   */
  public static String getToYear() {
    return PROPERTIES.getProperty("ToYear", "" + Calendar.getInstance().get(Calendar.YEAR));
  }
  
  /**
   * returns the entity owning the copyright
   * 
   * @return		the owner
   */
  public static String getOwner() {
    return PROPERTIES.getProperty("Owner", "The University of Waikato");
  }
  
  /**
   * returns the address of the owner
   * 
   * @return		the address
   */
  public static String getAddress() {
    return PROPERTIES.getProperty("Address", "Hamilton, New Zealand");
  }
  
  /**
   * returns the URL of the owner
   * 
   * @return		the URL
   */
  public static String getURL() {
    return PROPERTIES.getProperty("URL", "http://www.cs.waikato.ac.nz/~ml/");
  }
  
  /**
   * Only for testing
   * 
   * @param args	ignored
   */
  public static void main(String[] args) {
    System.out.println(PROPERTIES);
  }
}
