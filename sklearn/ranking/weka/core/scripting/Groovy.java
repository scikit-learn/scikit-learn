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
 * Groovy.java
 * Copyright (C) 2009 University of Waikato, Hamilton, New Zealand
 */

package weka.core.scripting;

import weka.core.RevisionHandler;
import weka.core.RevisionUtils;

import java.io.File;
import java.io.Serializable;
import java.lang.reflect.Constructor;
import java.lang.reflect.Method;

/**
 * A helper class for <a href="http://groovy.codehaus.org/" target="_blank">Groovy</a>.
 * <p/>
 * In order to use Groovy, the jar containing all the classes must be present
 * in the CLASSPATH. This jar is normally found in the <i>embeddable</i>
 * sub-directory of the Groovy installation.
 * <p/>
 * Tested with Groovy 1.5.7.
 * 
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5953 $
 */
public class Groovy
  implements Serializable, RevisionHandler {

  /** for serialization. */
  private static final long serialVersionUID = -2628766602043134673L;

  /** the classname of the Groovy classloader. */
  public final static String CLASS_GROOVYCLASSLOADER = "groovy.lang.GroovyClassLoader";
  
  /** whether the Groovy classes are in the Classpath. */
  protected static boolean m_Present = false;
  static {
    try {
      Class.forName(CLASS_GROOVYCLASSLOADER);
      m_Present = true;
    }
    catch (Exception e) {
      m_Present = false;
    }
  }
  
  /** the classloader. */
  protected Object m_ClassLoader;
  
  /**
   * default constructor, tries to instantiate a Groovy classloader.
   */
  public Groovy() {
    m_ClassLoader = newClassLoader();
  }
  
  /**
   * returns the currently used Groovy classloader.
   * 
   * @return		the classloader, can be null
   */
  public Object getClassLoader() {
    return m_ClassLoader;
  }
  
  /**
   * executes the specified method on the current interpreter and returns the 
   * result, if any.
   * 
   * @param methodName		the name of the method
   * @param paramClasses	the classes of the parameters
   * @param paramValues		the values of the parameters
   * @return			the return value of the method, if any (in that case null)
   */
  public Object invoke(String methodName, Class[] paramClasses, Object[] paramValues) {
    Object	result;
    
    result = null;
    if (getClassLoader() != null)
      result = invoke(getClassLoader(), methodName, paramClasses, paramValues);
    
    return result;
  }
  
  /**
   * returns whether the Groovy classes are present or not, i.e. whether the 
   * classes are in the classpath or not
   *
   * @return 			whether the Groovy classes are available
   */
  public static boolean isPresent() {
    return m_Present;
  }

  /**
   * initializes and returns a Groovy Interpreter.
   * 
   * @return			the interpreter or null if Groovy classes not present
   */
  public static Object newClassLoader() {
    Object	result;
    Class<?>	cls;
    Constructor	constr;
    
    result = null;
    
    if (isPresent()) {
      try {
	cls    = Class.forName(CLASS_GROOVYCLASSLOADER);
	constr = cls.getConstructor(new Class[]{ClassLoader.class});
	result = constr.newInstance(Groovy.class.getClassLoader());
      }
      catch (Exception e) {
	e.printStackTrace();
	result = null;
      }
    }

    return result;
  }

  /**
   * loads the module and returns a new instance of it as instance of the
   * provided Java class template.
   * 
   * @param file		the Groovy module file
   * @param template		the template for the returned Java object
   * @return			the Groovy object
   */
  public static Object newInstance(File file, Class template) {
    Object 	result;
    Object	interpreter;
    Class	cls;

    result = null;

    if (!isPresent())
      return result;
    
    interpreter = newClassLoader();
    if (interpreter == null)
      return result;

    try {
      cls    = (Class) invoke(interpreter, "parseClass", new Class[]{File.class}, new Object[]{file});
      result = cls.newInstance();
    }
    catch (Exception e) {
      e.printStackTrace();
    }

    return result;
  }
  
  /**
   * executes the specified method and returns the result, if any.
   * 
   * @param o			the object the method should be called from,
   * 				e.g., a Groovy Interpreter
   * @param methodName		the name of the method
   * @param paramClasses	the classes of the parameters
   * @param paramValues		the values of the parameters
   * @return			the return value of the method, if any (in that case null)
   */
  public static Object invoke(Object o, String methodName, Class[] paramClasses, Object[] paramValues) {
    Method      m;
    Object      result;
    
    result = null;
    
    try {
      m      = o.getClass().getMethod(methodName, paramClasses);
      result = m.invoke(o, paramValues);
    }
    catch (Exception e) {
      e.printStackTrace();
      result = null;
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
   * If no arguments are given, it just prints the presence of the Groovy
   * classes, otherwise it expects a Groovy filename to execute.
   * 
   * @param args		commandline arguments
   */
  public static void main(String[] args) {
    if (args.length == 0) {
      System.out.println("Groovy present: " + isPresent());
    }
    else {
      Groovy groovy = new Groovy();
      if (groovy.getClassLoader() == null) {
	System.err.println("Cannot instantiate Groovy ClassLoader!");
      }
      else {
	Object groovyObject = Groovy.newInstance(new File(args[0]), Object.class);
	Groovy.invoke(groovyObject, "run", new Class[]{}, new Object[]{});
      }
    }
  }
}
