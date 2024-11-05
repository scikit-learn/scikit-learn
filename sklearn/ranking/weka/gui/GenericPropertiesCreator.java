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
 * GenericPropertiesCreator.java
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */
package weka.gui;

import weka.core.ClassDiscovery;
import weka.core.Utils;
import weka.core.ClassDiscovery.StringCompare;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.util.Collections;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Properties;
import java.util.StringTokenizer;
import java.util.Vector;

/**
 * This class can generate the properties object that is normally loaded from
 * the <code>GenericObjectEditor.props</code> file (= PROPERTY_FILE). It takes
 * the <code>GenericPropertiesCreator.props</code> file as a template to
 * determine all the derived classes by checking the classes in the given
 * packages (a file with the same name in your home directory overrides the
 * the one in the weka/gui directory/package). <br>
 * E.g. if we want to have all the subclasses of the <code>Classifier</code>
 * class then we specify the superclass ("weka.classifiers.Classifier") and the
 * packages where to look for ("weka.classifiers.bayes" etc.):
 * 
 * <pre>
 * 
 *   weka.classifiers.Classifier=\
 *     weka.classifiers.bayes,\
 *     weka.classifiers.functions,\
 *     weka.classifiers.lazy,\
 *     weka.classifiers.meta,\
 *     weka.classifiers.trees,\
 *     weka.classifiers.rules
 *  
 * </pre>
 * 
 * This creates the same list as stored in the
 * <code>GenericObjectEditor.props</code> file, but it will also add
 * additional classes, that are not listed in the static list (e.g. a newly
 * developed Classifier), but still in the classpath. <br>
 * <br>
 * For discovering the subclasses the whole classpath is inspected, which means
 * that you can have several parallel directories with the same package
 * structure (e.g. a release directory and a developer directory with additional
 * classes). <br>
 * <br>
 * The dynamic discovery can be turned off via the <code>UseDyanmic</code>
 * property in the props file (values: true|false).
 * 
 * @see #CREATOR_FILE
 * @see #PROPERTY_FILE
 * @see #USE_DYNAMIC
 * @see GenericObjectEditor
 * @see weka.core.ClassDiscovery
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5968 $
 */
public class GenericPropertiesCreator {
  
  /** whether to output some debug information */
  public final static boolean VERBOSE = false;

  /** name of property whether to use the dynamic approach or the old 
   * GenericObjectEditor.props file */
  public final static String USE_DYNAMIC = "UseDynamic";
  
  /** The name of the properties file to use as a template. Contains the 
   * packages in which to look for derived classes. It has the same structure
   * as the <code>PROPERTY_FILE</code>
   * @see #PROPERTY_FILE
   */
  protected static String CREATOR_FILE = "weka/gui/GenericPropertiesCreator.props";
  
  /** The name of the properties file that lists classes/interfaces/superclasses
   * to exclude from being shown in the GUI. See the file for more information.
   */
  protected static String EXCLUDE_FILE = "weka/gui/GenericPropertiesCreator.excludes";

  /** the prefix for an interface exclusion */
  protected static String EXCLUDE_INTERFACE = "I";

  /** the prefix for an (exact) class exclusion */
  protected static String EXCLUDE_CLASS = "C";

  /** the prefix for a superclass exclusion */
  protected static String EXCLUDE_SUPERCLASS = "S";
  
  /** The name of the properties file for the static GenericObjectEditor 
   * (<code>USE_DYNAMIC</code> = <code>false</code>) 
   * @see GenericObjectEditor 
   * @see #USE_DYNAMIC 
   */
  protected static String PROPERTY_FILE = "weka/gui/GenericObjectEditor.props";
  
  /** the input file with the packages */
  protected String m_InputFilename;
  
  /** the output props file for the GenericObjectEditor */
  protected String m_OutputFilename;
  
  /** the "template" properties file with the layout and the packages */
  protected Properties m_InputProperties;
  
  /** the output properties file with the filled in classes */
  protected Properties m_OutputProperties;
  
  /** whether an explicit input file was given - if false, the Utils class
   * is used to locate the props-file */
  protected boolean m_ExplicitPropsFile;
  
  /** the hashtable that stores the excludes: 
   * key -&gt; Hashtable(prefix -&gt; Vector of classnames) */
  protected Hashtable m_Excludes;
  
  /**
   * initializes the creator, locates the props file with the Utils class.
   * 
   * @throws Exception if loading of CREATOR_FILE fails
   * @see #CREATOR_FILE
   * @see Utils#readProperties(String)
   * @see #loadInputProperties()
   */
  public GenericPropertiesCreator() throws Exception {
    this(CREATOR_FILE);
    m_ExplicitPropsFile = false;
  }

  /**
   * initializes the creator, the given file overrides the props-file search
   * of the Utils class
   * 
   * @param filename the file containing the packages to create a props file from
   * @throws Exception if loading of the file fails
   * @see #CREATOR_FILE
   * @see Utils#readProperties(String)
   * @see #loadInputProperties()
   */
  public GenericPropertiesCreator(String filename) throws Exception {
    super();
    m_InputFilename     = filename;
    m_OutputFilename    = PROPERTY_FILE;
    m_InputProperties   = null;
    m_OutputProperties  = null;
    m_ExplicitPropsFile = true;
    m_Excludes          = new Hashtable();
  }

  /**
   * if FALSE, the locating of a props-file of the Utils-class is used, 
   * otherwise it's tried to load the specified file
   * 
   * @param value 	if true the specified file will be loaded not via the 
   * 			Utils-class
   * @see Utils#readProperties(String)
   * @see #loadInputProperties()
   */
  public void setExplicitPropsFile(boolean value) {
    m_ExplicitPropsFile = value;
  }
  
  /**
   * returns TRUE, if a file is loaded and not the Utils class used for 
   * locating the props file.
   * 
   * @return	true if the specified file is used and not the one found by
   * 		the Utils class
   * @see Utils#readProperties(String)
   * @see #loadInputProperties()
   */
  public boolean getExplicitPropsFile() {
    return m_ExplicitPropsFile;
  }
  
  /**
   * returns the name of the output file
   * 
   * @return the name of the output file 
   */
  public String getOutputFilename() {
    return m_OutputFilename;
  }
  
  /**
   * sets the file to output the properties for the GEO to
   * 
   * @param filename the filename for the output 
   */
  public void setOutputFilename(String filename) {
    m_OutputFilename = filename;
  }

  /**
   * returns the name of the input file
   * 
   * @return the name of the input file 
   */
  public String getInputFilename() {
    return m_InputFilename;
  }
  
  /**
   * sets the file to get the information about the packages from. 
   * automatically sets explicitPropsFile to TRUE.
   * 
   * @param filename the filename for the input 
   * @see #setExplicitPropsFile(boolean)
   */
  public void setInputFilename(String filename) {
    m_InputFilename = filename;
    setExplicitPropsFile(true);
  }
  
  /**
   * returns the input properties object (template containing the packages)
   * 
   * @return		the input properties (the template)
   */
  public Properties getInputProperties() {
    return m_InputProperties;
  }
  
  /**
   * returns the output properties object (structure like the template, but
   * filled with classes instead of packages)
   * 
   * @return		the output properties (filled with classes)
   */
  public Properties getOutputProperties() {
    return m_OutputProperties;
  }
  
  /**
   * loads the property file containing the layout and the packages of 
   * the output-property-file. The exlcude property file is also read here.
   * 
   * @see #m_InputProperties
   * @see #m_InputFilename
   */
  protected void loadInputProperties()  {
    if (VERBOSE)
      System.out.println("Loading '" + getInputFilename() + "'...");
    m_InputProperties = new Properties();
    try {
      File f = new File(getInputFilename());
      if (getExplicitPropsFile() && f.exists())
        m_InputProperties.load(new FileInputStream(getInputFilename()));
      else
        m_InputProperties = Utils.readProperties(getInputFilename());
      
      // excludes
      m_Excludes.clear();
      Properties p = Utils.readProperties(EXCLUDE_FILE);
      Enumeration enm = p.propertyNames();
      while (enm.hasMoreElements()) {
	String name = enm.nextElement().toString();
	// new Hashtable for key
	Hashtable t = new Hashtable();
	m_Excludes.put(name, t);
	t.put(EXCLUDE_INTERFACE, new Vector());
	t.put(EXCLUDE_CLASS, new Vector());
	t.put(EXCLUDE_SUPERCLASS, new Vector());
	
	// process entries
	StringTokenizer tok = new StringTokenizer(p.getProperty(name), ",");
	while (tok.hasMoreTokens()) {
	  String item = tok.nextToken();
	  // get list
	  Vector list = new Vector();
	  if (item.startsWith(EXCLUDE_INTERFACE + ":"))
	    list = (Vector) t.get(EXCLUDE_INTERFACE);
	  else if (item.startsWith(EXCLUDE_CLASS + ":"))
	    list = (Vector) t.get(EXCLUDE_CLASS);
	  else if (item.startsWith(EXCLUDE_SUPERCLASS))
	    list = (Vector) t.get(EXCLUDE_SUPERCLASS);
	  // add to list
	  list.add(item.substring(item.indexOf(":") + 1));
	}
      }
    }
    catch (Exception e) {
      e.printStackTrace();
    }
  }

  /**
   * gets whether the dynamic approach should be used or not
   * 
   * @return		true if the dynamic approach is to be used
   */
  public boolean useDynamic() {
    if (getInputProperties() == null)
      loadInputProperties();
    
    // check our classloader against the system one - if different then
    // return false (as dynamic classloading only works for classes discoverable
    // in the system classpath
    /*if (!ClassLoader.getSystemClassLoader().equals(this.getClass().getClassLoader())) {
      if (Boolean.parseBoolean(getInputProperties().getProperty(USE_DYNAMIC, "true")) == true) {
        System.out.println("[GenericPropertiesCreator] classloader in use is not the system "
            + "classloader: using static entries in weka/gui/GenericObjectEditor.props rather "
            + "than dynamic class discovery.");
      }
      return false;
    }*/
    
    return Boolean.parseBoolean(
	getInputProperties().getProperty(USE_DYNAMIC, "true"));
  }
  
  /**
   * checks whether the classname is a valid one, i.e., from a public class
   *
   * @param classname   the classname to check
   * @return            whether the classname is a valid one
   */
  protected boolean isValidClassname(String classname) {
    return (classname.indexOf("$") == -1);
  }

  /**
   * Checks whether the classname is a valid one for the given key. This is
   * based on the settings in the Exclude file.
   *
   * @param key         the property key
   * @param classname   the classname to check
   * @return            whether the classname is a valid one
   * @see		#EXCLUDE_FILE
   */
  protected boolean isValidClassname(String key, String classname) {
    boolean	result;
    Class	cls;
    Class	clsCurrent;
    Vector	list;
    int		i;

    result = true;

    // are there excludes for this key?
    if (m_Excludes.containsKey(key)) {
      try {
	clsCurrent = Class.forName(classname);
      }
      catch (Exception e) {
	// we ignore this Exception
	clsCurrent = null;
      }
      
      // interface
      if ((clsCurrent != null) && result) {
	list = (Vector) ((Hashtable) m_Excludes.get(key)).get(EXCLUDE_INTERFACE);
	for (i = 0; i < list.size(); i++) {
	  try {
	    cls = Class.forName(list.get(i).toString());
	    if (ClassDiscovery.hasInterface(cls, clsCurrent)) {
	      result = false;
	      break;
	    }
	  }
	  catch (Exception e) {
	    // we ignore this Exception
	  }
	}
      }
      
      // superclass
      if ((clsCurrent != null) && result) {
	list = (Vector) ((Hashtable) m_Excludes.get(key)).get(EXCLUDE_SUPERCLASS);
	for (i = 0; i < list.size(); i++) {
	  try {
	    cls = Class.forName(list.get(i).toString());
	    if (ClassDiscovery.isSubclass(cls, clsCurrent)) {
	      result = false;
	      break;
	    }
	  }
	  catch (Exception e) {
	    // we ignore this Exception
	  }
	}
      }
      
      // class
      if ((clsCurrent != null) && result) {
	list = (Vector) ((Hashtable) m_Excludes.get(key)).get(EXCLUDE_CLASS);
	for (i = 0; i < list.size(); i++) {
	  try {
	    cls = Class.forName(list.get(i).toString());
	    if (cls.getName().equals(clsCurrent.getName()))
	      result = false;
	  }
	  catch (Exception e) {
	    // we ignore this Exception
	  }
	}
      }
    }

    return result;
  }
  
  /**
   * fills in all the classes (based on the packages in the input properties 
   * file) into the output properties file
   *     
   * @throws Exception 	if something goes wrong
   * @see #m_OutputProperties
   */
  protected void generateOutputProperties() throws Exception {
    Enumeration       keys;
    String            key;
    String            value;
    String            pkg;
    StringTokenizer   tok;
    Vector            classes;
    HashSet           names;
    int               i;
    
    m_OutputProperties = new Properties();
    keys               = m_InputProperties.propertyNames();
    while (keys.hasMoreElements()) {
      key   = keys.nextElement().toString();
      if (key.equals(USE_DYNAMIC))
	continue;
      tok   = new StringTokenizer(m_InputProperties.getProperty(key), ",");
      names = new HashSet();
      
      // get classes for all packages
      while (tok.hasMoreTokens()) {
        pkg = tok.nextToken().trim();
        
        try {
          classes = ClassDiscovery.find(Class.forName(key), pkg);
        }
        catch (Exception e) {
          System.out.println("Problem with '" + key + "': " + e);
          classes = new Vector();
        }
        
        for (i = 0; i < classes.size(); i++) {
          // skip non-public classes
          if (!isValidClassname(classes.get(i).toString()))
            continue;
          // some classes should not be listed for some keys
          if (!isValidClassname(key, classes.get(i).toString()))
            continue;
          names.add(classes.get(i));
        }
      }
      
      // generate list
      value   = "";
      classes = new Vector();
      classes.addAll(names);
      Collections.sort(classes, new StringCompare());
      for (i = 0; i < classes.size(); i++) {
        if (!value.equals(""))
          value += ",";
        value += classes.get(i).toString();
      }
      if (VERBOSE)
        System.out.println(pkg + " -> " + value);
      
      // set value
      m_OutputProperties.setProperty(key, value);
    }
  }
  
  /**
   * stores the generated output properties file
   * 
   * @throws Exception	if the saving fails
   * @see #m_OutputProperties
   * @see #m_OutputFilename 
   */
  protected void storeOutputProperties() throws Exception {
    if (VERBOSE)
      System.out.println("Saving '" + getOutputFilename() + "'...");
    m_OutputProperties.store(
        new FileOutputStream(getOutputFilename()), 
        " Customises the list of options given by the GenericObjectEditor\n# for various superclasses.");
  }
  
  /**
   * generates the props-file for the GenericObjectEditor and stores it
   * 
   * @throws Exception	if something goes wrong
   * @see #execute(boolean)
   */
  public void execute() throws Exception {
    execute(true);
  }
  
  /**
   * generates the props-file for the GenericObjectEditor and stores it only
   * if the the param <code>store</code> is TRUE. If it is FALSE then the
   * generated properties file can be retrieved via the <code>getOutputProperties</code>
   * method. 
   * 
   * @param store     	if TRUE then the properties file is stored to the stored 
   *                  	filename
   * @throws Exception	if something goes wrong
   * @see #getOutputFilename()
   * @see #setOutputFilename(String)
   * @see #getOutputProperties()
   */
  public void execute(boolean store) throws Exception {
    // read properties file
    loadInputProperties();
    
    // generate the props file
    generateOutputProperties();
    
    // write properties file
    if (store)
      storeOutputProperties();
  }
  
  /**
   * for generating props file:
   * <ul>
   *   <li>
   *     no parameter: 
   *     see default constructor
   *   </li>
   *   <li>
   *     1 parameter (i.e., filename): 
   *     see default constructor + setOutputFilename(String)
   *   </li>
   *   <li>
   *     2 parameters (i.e, filenames): 
   *     see constructor with String argument + setOutputFilename(String)
   *   </li>
   * </ul>
   *
   * @param args	the commandline arguments
   * @throws Exception 	if something goes wrong
   * @see #GenericPropertiesCreator()
   * @see #GenericPropertiesCreator(String)
   * @see #setOutputFilename(String)
   */
  public static void main(String[] args) throws Exception {
    GenericPropertiesCreator   c = null;
    
    if (args.length == 0) {
      c = new GenericPropertiesCreator();
    }
    else if (args.length == 1) {
      c = new GenericPropertiesCreator();
      c.setOutputFilename(args[0]);
    }
    else if (args.length == 2) {
      c = new GenericPropertiesCreator(args[0]);
      c.setOutputFilename(args[1]);
    }
    else {
      System.out.println("usage: " + GenericPropertiesCreator.class.getName() + " [<input.props>] [<output.props>]");
      System.exit(1);
    }
    
    c.execute(true);
  }
}
