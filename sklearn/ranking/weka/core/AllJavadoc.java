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
 * AllJavadoc.java
 * Copyright (C) 2006 University of Waikato, Hamilton, New Zealand
 */

package weka.core;

import java.util.HashSet;
import java.util.Vector;

/**
 * Applies all known Javadoc-derived classes to a source file.
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -W &lt;classname&gt;
 *  The class to load.</pre>
 * 
 * <pre> -nostars
 *  Suppresses the '*' in the Javadoc.</pre>
 * 
 * <pre> -dir &lt;dir&gt;
 *  The directory above the package hierarchy of the class.</pre>
 * 
 * <pre> -silent
 *  Suppresses printing in the console.</pre>
 * 
 <!-- options-end -->
 * 
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 6111 $
 */
public class AllJavadoc
  extends Javadoc {

  /** contains all the  */
  protected static Vector<Javadoc> m_Javadocs;
  
  /** determine all classes derived from Javadoc and instantiate them */
  static {
    // get all classnames, besides this one
    HashSet<String> set = new HashSet<String>(ClassDiscovery.find(Javadoc.class, Javadoc.class.getPackage().getName()));
    if (set.contains(AllJavadoc.class.getName()))
      set.remove(AllJavadoc.class.getName());
    
    // instantiate them
    m_Javadocs = new Vector<Javadoc>();
    for (String classname: set) {
      try {
	Class cls = Class.forName(classname);
	m_Javadocs.add((Javadoc)cls.newInstance());
      }
      catch (Exception e) {
	e.printStackTrace();
      }
    }
  }
  
  /**
   * sets the classname of the class to generate the Javadoc for
   * 
   * @param value	the new classname
   */
  public void setClassname(String value) {
    super.setClassname(value);
    for (int i = 0; i < m_Javadocs.size(); i++)
      ((Javadoc) m_Javadocs.get(i)).setClassname(value);
  }
  
  /**
   * sets whether to prefix the Javadoc with "*"
   * 
   * @param value	true if stars are to be used
   */
  public void setUseStars(boolean value) {
    super.setUseStars(value);
    for (int i = 0; i < m_Javadocs.size(); i++)
      ((Javadoc) m_Javadocs.get(i)).setUseStars(value);
  }
  
  /**
   * sets whether to suppress output in the console
   * 
   * @param value	true if output is to be suppressed
   */
  public void setSilent(boolean value) {
    super.setSilent(value);
    for (int i = 0; i < m_Javadocs.size(); i++)
      ((Javadoc) m_Javadocs.get(i)).setSilent(value);
  }

  /**
   * generates and returns the Javadoc for the specified start/end tag pair.
   * 
   * @param index	the index in the start/end tag array
   * @return		the generated Javadoc
   * @throws Exception 	in case the generation fails
   */
  protected String generateJavadoc(int index) throws Exception {
    throw new Exception("Not used!");
  }
  
  /**
   * updates the Javadoc in the given source code, using all the found
   * Javadoc updaters.
   * 
   * @param content	the source code
   * @return		the updated source code
   * @throws Exception 	in case the generation fails
   */
  protected String updateJavadoc(String content) throws Exception {
    String	result;
    int		i;
    
    result = content;
    
    for (i = 0; i < m_Javadocs.size(); i++) {
      result = ((Javadoc) m_Javadocs.get(i)).updateJavadoc(result);
    }
    
    return result;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6111 $");
  }

  /**
   * Parses the given commandline parameters and generates the Javadoc.
   * 
   * @param args	the commandline parameters for the object
   */
  public static void main(String[] args) {
    runJavadoc(new AllJavadoc(), args);
  }
}
