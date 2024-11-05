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
 * Capabilities.java
 * Copyright (C) 2006 University of Waikato, Hamilton, New Zealand
 */

package weka.core;

import weka.core.converters.ConverterUtils.DataSource;
import weka.core.labelranking.PreferenceAttribute;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Properties;
import java.util.Vector;

/**
 * A class that describes the capabilites (e.g., handling certain types of
 * attributes, missing values, types of classes, etc.) of a specific
 * classifier. By default, the classifier is capable of nothing. This
 * ensures that new features have to be enabled explicitly. <p/>
 * 
 * A common code fragment for making use of the capabilities in a classifier 
 * would be this:
 * <pre>
 * public void <b>buildClassifier</b>(Instances instances) throws Exception {
 *   // can the classifier handle the data?
 *   getCapabilities().<b>testWithFail(instances)</b>;
 *   ...
 *   // possible deletion of instances with missing class labels, etc.
 * </pre>
 * For only testing a single attribute, use this:
 * <pre>
 *   ...
 *   Attribute att = instances.attribute(0);
 *   getCapabilities().<b>testWithFail(att)</b>;
 *   ...
 * </pre>
 * Or for testing the class attribute (uses the capabilities that are 
 * especially for the class):
 * <pre>
 *   ...
 *   Attribute att = instances.classAttribute();
 *   getCapabilities().<b>testWithFail(att, <i>true</i>)</b>;
 *   ...
 * </pre>
 * 
 * @author  FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 6447 $
 */
public class Capabilities 
  implements Cloneable, Serializable, RevisionHandler {
  
  /** serialversion UID */
  static final long serialVersionUID = -5478590032325567849L;  

  /** the properties file for managing the tests */
  public final static String PROPERTIES_FILE = "weka/core/Capabilities.props";

  /** the actual properties */
  protected static Properties PROPERTIES;
  
  /** defines an attribute type */
  private final static int ATTRIBUTE = 1;
  
  /** defines a class type */
  private final static int CLASS = 2;
  
  /** defines an attribute capability */
  private final static int ATTRIBUTE_CAPABILITY = 4;
  
  /** defines a class capability */
  private final static int CLASS_CAPABILITY = 8;
  
  /** defines a other capability */
  private final static int OTHER_CAPABILITY = 16;

  /** enumeration of all capabilities */
  public enum Capability {
    // attributes
    /** can handle nominal attributes */
    NOMINAL_ATTRIBUTES(ATTRIBUTE + ATTRIBUTE_CAPABILITY, "Nominal attributes"),
    /** can handle binary attributes */
    BINARY_ATTRIBUTES(ATTRIBUTE + ATTRIBUTE_CAPABILITY, "Binary attributes"),
    /** can handle unary attributes */
    UNARY_ATTRIBUTES(ATTRIBUTE + ATTRIBUTE_CAPABILITY, "Unary attributes"),
    /** can handle empty nominal attributes */
    EMPTY_NOMINAL_ATTRIBUTES(ATTRIBUTE + ATTRIBUTE_CAPABILITY, "Empty nominal attributes"),
    /** can handle numeric attributes */
    NUMERIC_ATTRIBUTES(ATTRIBUTE + ATTRIBUTE_CAPABILITY, "Numeric attributes"),
    /** can handle date attributes */
    DATE_ATTRIBUTES(ATTRIBUTE + ATTRIBUTE_CAPABILITY, "Date attributes"),
    /** can handle string attributes */
    STRING_ATTRIBUTES(ATTRIBUTE + ATTRIBUTE_CAPABILITY, "String attributes"),
    /** can handle relational attributes */
    RELATIONAL_ATTRIBUTES(ATTRIBUTE + ATTRIBUTE_CAPABILITY, "Relational attributes"),
    /** can handle missing values in attributes */
    MISSING_VALUES(ATTRIBUTE_CAPABILITY, "Missing values"),
    //RANKING BEGIN
    /** can handle preference attributes.*/
    PREFERENCE_ATTRIBUTE(ATTRIBUTE + ATTRIBUTE_CAPABILITY, "Preference attributes"),
    //RANKING END
    // class
    /** can handle data without class attribute, eg clusterers */
    NO_CLASS(CLASS_CAPABILITY, "No class"),
    /** can handle nominal classes */
    NOMINAL_CLASS(CLASS + CLASS_CAPABILITY, "Nominal class"),   
    /** can handle binary classes */
    BINARY_CLASS(CLASS + CLASS_CAPABILITY, "Binary class"),
    /** can handle unary classes */
    UNARY_CLASS(CLASS + CLASS_CAPABILITY, "Unary class"),
    /** can handle empty nominal classes */
    EMPTY_NOMINAL_CLASS(CLASS + CLASS_CAPABILITY, "Empty nominal class"),
    /** can handle numeric classes */
    NUMERIC_CLASS(CLASS + CLASS_CAPABILITY, "Numeric class"),
    /** can handle date classes */
    DATE_CLASS(CLASS + CLASS_CAPABILITY, "Date class"),
    /** can handle string classes */
    STRING_CLASS(CLASS + CLASS_CAPABILITY, "String class"),
    /** can handle relational classes */
    RELATIONAL_CLASS(CLASS + CLASS_CAPABILITY, "Relational class"),
    /** can handle missing values in class attribute */
    MISSING_CLASS_VALUES(CLASS_CAPABILITY, "Missing class values"),
    //RANKING BEGIN
    /** can handle rankings */
    RANKING(CLASS + CLASS_CAPABILITY,"Rankings"),
    //RANKING END
    // other
    /** can handle multi-instance data */
    ONLY_MULTIINSTANCE(OTHER_CAPABILITY, "Only multi-Instance data");

    /** the flags for the capabilities */
    private int m_Flags = 0;
    
    /** the display string */
    private String m_Display;
    
    /**
     * initializes the capability with the given flags
     * 
     * @param flags	"meta-data" for the capability
     * @param display	the display string (must be unique!)
     */
    private Capability(int flags, String display) {
      m_Flags   = flags;
      m_Display = display;
    }
    
    /**
     * returns true if the capability is an attribute
     * 
     * @return true if the capability is an attribute
     */
    public boolean isAttribute() {
      return ((m_Flags & ATTRIBUTE) == ATTRIBUTE);
    }
    
    /**
     * returns true if the capability is a class
     * 
     * @return true if the capability is a class
     */
    public boolean isClass() {
      return ((m_Flags & CLASS) == CLASS);
    }
    
    /**
     * returns true if the capability is an attribute capability
     * 
     * @return true if the capability is an attribute capability
     */
    public boolean isAttributeCapability() {
      return ((m_Flags & ATTRIBUTE_CAPABILITY) == ATTRIBUTE_CAPABILITY);
    }
    
    /**
     * returns true if the capability is a class capability
     * 
     * @return true if the capability is a class capability
     */
    public boolean isOtherCapability() {
      return ((m_Flags & OTHER_CAPABILITY) == OTHER_CAPABILITY);
    }
    
    /**
     * returns true if the capability is a other capability
     * 
     * @return true if the capability is a other capability
     */
    public boolean isClassCapability() {
      return ((m_Flags & CLASS_CAPABILITY) == CLASS_CAPABILITY);
    }
    
    /**
     * returns the display string of the capability
     * 
     * @return the display string
     */
    public String toString() {
      return m_Display;
    }
  };

  /** the object that owns this capabilities instance */
  protected CapabilitiesHandler m_Owner;
  
  /** the hashset for storing the active capabilities */
  protected HashSet<Capability> m_Capabilities;
  
  /** the hashset for storing dependent capabilities, eg for meta-classifiers */
  protected HashSet<Capability> m_Dependencies;
  
  /** the reason why the test failed, used to throw an exception */
  protected Exception m_FailReason = null;

  /** the minimum number of instances in a dataset */
  protected int m_MinimumNumberInstances = 1;

  /** whether to perform any tests at all */
  protected boolean m_Test;

  /** whether to perform data based tests */
  protected boolean m_InstancesTest;

  /** whether to perform attribute based tests */
  protected boolean m_AttributeTest;

  /** whether to test for missing values */
  protected boolean m_MissingValuesTest;

  /** whether to test for missing class values */
  protected boolean m_MissingClassValuesTest;

  /** whether to test for minimum number of instances */
  protected boolean m_MinimumNumberInstancesTest;
  
  /**
   * initializes the capabilities for the given owner
   * 
   * @param owner       the object that produced this Capabilities instance
   */
  public Capabilities(CapabilitiesHandler owner) {
    super();

    setOwner(owner);
    m_Capabilities = new HashSet<Capability>();
    m_Dependencies = new HashSet<Capability>();

    // load properties
    if (PROPERTIES == null) {
      try {
        PROPERTIES = Utils.readProperties(PROPERTIES_FILE);
      }
      catch (Exception e) {
	e.printStackTrace();
	PROPERTIES = new Properties();
      }
    }
    
    m_Test                       = Boolean.parseBoolean(PROPERTIES.getProperty("Test", "true"));
    m_InstancesTest              = Boolean.parseBoolean(PROPERTIES.getProperty("InstancesTest", "true")) && m_Test;
    m_AttributeTest              = Boolean.parseBoolean(PROPERTIES.getProperty("AttributeTest", "true")) && m_Test;
    m_MissingValuesTest          = Boolean.parseBoolean(PROPERTIES.getProperty("MissingValuesTest", "true")) && m_Test;
    m_MissingClassValuesTest     = Boolean.parseBoolean(PROPERTIES.getProperty("MissingClassValuesTest", "true")) && m_Test;
    m_MinimumNumberInstancesTest = Boolean.parseBoolean(PROPERTIES.getProperty("MinimumNumberInstancesTest", "true")) && m_Test;
  }
  
  /**
   * Creates and returns a copy of this object.
   * 
   * @return	a clone of this object
   */
  public Object clone() {
    Capabilities    result;

    result = new Capabilities(m_Owner);
    result.assign(this);

    return result;
  }
  
  /**
   * retrieves the data from the given Capabilities object
   * 
   * @param c	  the capabilities object to initialize with
   */
  public void assign(Capabilities c) {
    for (Capability cap: Capability.values()) {
      // capability
      if (c.handles(cap))
        enable(cap);
      else
	disable(cap);
      // dependency
      if (c.hasDependency(cap))
        enableDependency(cap);
      else
	disableDependency(cap);
    }

    setMinimumNumberInstances(c.getMinimumNumberInstances());
  }

  /**
   * performs an AND conjunction with the capabilities of the given 
   * Capabilities object and updates itself
   *
   * @param c     the capabilities to AND with
   */
  public void and(Capabilities c) {
    for (Capability cap: Capability.values()) {
      // capability
      if (handles(cap) && c.handles(cap))
        m_Capabilities.add(cap);
      else
        m_Capabilities.remove(cap);
      // dependency
      if (hasDependency(cap) && c.hasDependency(cap))
        m_Dependencies.add(cap);
      else
        m_Dependencies.remove(cap);
    }
    
    // minimum number of instances that both handlers need at least to work
    if (c.getMinimumNumberInstances() > getMinimumNumberInstances())
      setMinimumNumberInstances(c.getMinimumNumberInstances());
  }

  /**
   * performs an OR conjunction with the capabilities of the given 
   * Capabilities object and updates itself
   *
   * @param c     the capabilities to OR with
   */
  public void or(Capabilities c) {
    for (Capability cap: Capability.values()) {
      // capability
      if (handles(cap) || c.handles(cap))
        m_Capabilities.add(cap);
      else
        m_Capabilities.remove(cap);
      // dependency
      if (hasDependency(cap) || c.hasDependency(cap))
        m_Dependencies.add(cap);
      else
        m_Dependencies.remove(cap);
    }
    
    if (c.getMinimumNumberInstances() < getMinimumNumberInstances())
      setMinimumNumberInstances(c.getMinimumNumberInstances());
  }
  
  /**
   * Returns true if the currently set capabilities support at least all of
   * the capabiliites of the given Capabilities object (checks only the enum!)
   * 
   * @param c	the capabilities to support at least
   * @return	true if all the requested capabilities are supported
   */
  public boolean supports(Capabilities c) {
    boolean	result;
    
    result = true;
    
    for (Capability cap: Capability.values()) {
      if (c.handles(cap) && !handles(cap)) {
	result = false;
	break;
      }
    }

    return result;
  }
  
  /**
   * Returns true if the currently set capabilities support (or have a 
   * dependency) at least all of the capabilities of the given Capabilities 
   * object (checks only the enum!)
   * 
   * @param c	the capabilities (or dependencies) to support at least
   * @return	true if all the requested capabilities are supported (or at 
   * 		least have a dependency)
   */
  public boolean supportsMaybe(Capabilities c) {
    boolean	result;
    
    result = true;
    
    for (Capability cap: Capability.values()) {
      if (c.handles(cap) && !(handles(cap) || hasDependency(cap))) {
	result = false;
	break;
      }
    }

    return result;
  }

  /**
   * sets the owner of this capabilities object
   * 
   * @param value       the new owner
   */
  public void setOwner(CapabilitiesHandler value) {
    m_Owner = value;
  }
  
  /**
   * returns the owner of this capabilities object
   * 
   * @return            the current owner of this capabilites object
   */
  public CapabilitiesHandler getOwner() {
    return m_Owner;
  }

  /**
   * sets the minimum number of instances that have to be in the dataset
   * 
   * @param value       the minimum number of instances
   */
  public void setMinimumNumberInstances(int value) {
    if (value >= 0)
      m_MinimumNumberInstances = value;
  }
  
  /**
   * returns the minimum number of instances that have to be in the dataset
   * 
   * @return            the minimum number of instances
   */
  public int getMinimumNumberInstances() {
    return m_MinimumNumberInstances;
  }
  
  /**
   * Returns an Iterator over the stored capabilities
   * 
   * @return iterator over the current capabilities
   */
  public Iterator capabilities() {
    return m_Capabilities.iterator();
  }
  
  /**
   * Returns an Iterator over the stored dependencies
   * 
   * @return iterator over the current dependencies
   */
  public Iterator dependencies() {
    return m_Dependencies.iterator();
  }
  
  /**
   * enables the given capability. 
   * Enabling NOMINAL_ATTRIBUTES also enables BINARY_ATTRIBUTES, 
   * UNARY_ATTRIBUTES and EMPTY_NOMINAL_ATTRIBUTES. 
   * Enabling BINARY_ATTRIBUTES also enables UNARY_ATTRIBUTES and 
   * EMPTY_NOMINAL_ATTRIBUTES. 
   * Enabling UNARY_ATTRIBUTES also enables EMPTY_NOMINAL_ATTRIBUTES.
   * But NOMINAL_CLASS only enables BINARY_CLASS, since normal schemes in Weka
   * don't work with datasets that have only 1 class label (or none).
   *
   * @param c     the capability to enable
   */
  public void enable(Capability c) {
    // attributes
    if (c == Capability.NOMINAL_ATTRIBUTES) {
      enable(Capability.BINARY_ATTRIBUTES);
    }
    else if (c == Capability.BINARY_ATTRIBUTES) {
      enable(Capability.UNARY_ATTRIBUTES);
    }
    else if (c == Capability.UNARY_ATTRIBUTES) {
      enable(Capability.EMPTY_NOMINAL_ATTRIBUTES);
    }
    //RANKING BEGIN
    else if (c == Capability.PREFERENCE_ATTRIBUTE){
    	enable(Capability.NOMINAL_ATTRIBUTES);
    }
    //RANKING END
    // class
    else if (c == Capability.NOMINAL_CLASS) {
      enable(Capability.BINARY_CLASS);
    }
    //RANKING BEGIN
    else if (c == Capability.RANKING){
    	//enable(Capability.NOMINAL_CLASS);
    	enable(Capability.PREFERENCE_ATTRIBUTE);
    }
    //RANKING END

    m_Capabilities.add(c);
  }
  
  /**
   * enables the dependency flag for the given capability
   * Enabling NOMINAL_ATTRIBUTES also enables BINARY_ATTRIBUTES, 
   * UNARY_ATTRIBUTES and EMPTY_NOMINAL_ATTRIBUTES. 
   * Enabling BINARY_ATTRIBUTES also enables UNARY_ATTRIBUTES and 
   * EMPTY_NOMINAL_ATTRIBUTES. 
   * Enabling UNARY_ATTRIBUTES also enables EMPTY_NOMINAL_ATTRIBUTES.
   * But NOMINAL_CLASS only enables BINARY_CLASS, since normal schemes in Weka
   * don't work with datasets that have only 1 class label (or none).
   *
   * @param c     the capability to enable the dependency flag for
   */
  public void enableDependency(Capability c) {
    // attributes
    if (c == Capability.NOMINAL_ATTRIBUTES) {
      enableDependency(Capability.BINARY_ATTRIBUTES);
    }
    else if (c == Capability.BINARY_ATTRIBUTES) {
      enableDependency(Capability.UNARY_ATTRIBUTES);
    }
    else if (c == Capability.UNARY_ATTRIBUTES) {
      enableDependency(Capability.EMPTY_NOMINAL_ATTRIBUTES);
    }
    //RANKING BEGIN
    //else if (c == Capability.PREFERENCE_ATTRIBUTES){
    //	enable(Capability.NOMINAL_ATTRIBUTES);
    //}
    //RANKING END
    // class
    else if (c == Capability.NOMINAL_CLASS) {
      enableDependency(Capability.BINARY_CLASS);
    }
    //RANKING BEGIN
   // else if (c == Capability.RANKING){
   // 	enable(Capability.NOMINAL_CLASS);
   // }
    //RANKING END

    m_Dependencies.add(c);
  }
  
  /**
   * enables all class types
   * 
   * @see #disableAllClasses()
   * @see #getClassCapabilities()
   */
  public void enableAllClasses() {
    for (Capability cap: Capability.values()) {
      if (cap.isClass())
	enable(cap);
    }
  }
  
  /**
   * enables all class type dependencies
   * 
   * @see #disableAllClassDependencies()
   * @see #getClassCapabilities()
   */
  public void enableAllClassDependencies() {
    for (Capability cap: Capability.values()) {
      if (cap.isClass())
	enableDependency(cap);
    }
  }
  
  /**
   * enables all attribute types
   * 
   * @see #disableAllAttributes()
   * @see #getAttributeCapabilities()
   */
  public void enableAllAttributes() {
    for (Capability cap: Capability.values()) {
      if (cap.isAttribute())
	enable(cap);
    }
  }
  
  /**
   * enables all attribute type dependencies
   * 
   * @see #disableAllAttributeDependencies()
   * @see #getAttributeCapabilities()
   */
  public void enableAllAttributeDependencies() {
    for (Capability cap: Capability.values()) {
      if (cap.isAttribute())
	enableDependency(cap);
    }
  }
  
  /**
   * enables all attribute and class types (including dependencies)
   */
  public void enableAll() {
    enableAllAttributes();
    enableAllAttributeDependencies();
    enableAllClasses();
    enableAllClassDependencies();
    enable(Capability.MISSING_VALUES);
    enable(Capability.MISSING_CLASS_VALUES);
  }

  /**
   * disables the given capability
   * Disabling NOMINAL_ATTRIBUTES also disables BINARY_ATTRIBUTES, 
   * UNARY_ATTRIBUTES and EMPTY_NOMINAL_ATTRIBUTES. 
   * Disabling BINARY_ATTRIBUTES also disables UNARY_ATTRIBUTES and 
   * EMPTY_NOMINAL_ATTRIBUTES. 
   * Disabling UNARY_ATTRIBUTES also disables EMPTY_NOMINAL_ATTRIBUTES.
   * The same hierarchy applies to the class capabilities.
   *
   * @param c     the capability to disable
   */
  public void disable(Capability c) {
    // attributes
    if (c == Capability.NOMINAL_ATTRIBUTES) {
      disable(Capability.BINARY_ATTRIBUTES);
    }
    else if (c == Capability.BINARY_ATTRIBUTES) {
      disable(Capability.UNARY_ATTRIBUTES);
    }
    else if (c == Capability.UNARY_ATTRIBUTES) {
      disable(Capability.EMPTY_NOMINAL_ATTRIBUTES);
    }
    // class
    else if (c == Capability.NOMINAL_CLASS) {
      disable(Capability.BINARY_CLASS);
    }
    else if (c == Capability.BINARY_CLASS) {
      disable(Capability.UNARY_CLASS);
    }
    else if (c == Capability.UNARY_CLASS) {
      disable(Capability.EMPTY_NOMINAL_CLASS);
    }
    //RANKING BEGIN
    
    else if (c == Capability.RANKING){
    	//disable(Capability.NOMINAL_CLASS);
    }
    /*
    else if (c == Capability.PREFERENCE_ATTRIBUTE){
    	disable(Capability.RANKING);
    }*/
    //RANKING END

    m_Capabilities.remove(c);
  }

  /**
   * disables the dependency of the given capability
   * Disabling NOMINAL_ATTRIBUTES also disables BINARY_ATTRIBUTES, 
   * UNARY_ATTRIBUTES and EMPTY_NOMINAL_ATTRIBUTES. 
   * Disabling BINARY_ATTRIBUTES also disables UNARY_ATTRIBUTES and 
   * EMPTY_NOMINAL_ATTRIBUTES. 
   * Disabling UNARY_ATTRIBUTES also disables EMPTY_NOMINAL_ATTRIBUTES.
   * The same hierarchy applies to the class capabilities.
   *
   * @param c     the capability to disable the dependency flag for
   */
  public void disableDependency(Capability c) {
    // attributes
    if (c == Capability.NOMINAL_ATTRIBUTES) {
      disableDependency(Capability.BINARY_ATTRIBUTES);
    }
    else if (c == Capability.BINARY_ATTRIBUTES) {
      disableDependency(Capability.UNARY_ATTRIBUTES);
    }
    else if (c == Capability.UNARY_ATTRIBUTES) {
      disableDependency(Capability.EMPTY_NOMINAL_ATTRIBUTES);
    }
    // class
    else if (c == Capability.NOMINAL_CLASS) {
      disableDependency(Capability.BINARY_CLASS);
    }
    else if (c == Capability.BINARY_CLASS) {
      disableDependency(Capability.UNARY_CLASS);
    }
    else if (c == Capability.UNARY_CLASS) {
      disableDependency(Capability.EMPTY_NOMINAL_CLASS);
    }

    m_Dependencies.remove(c);
  }
  
  /**
   * disables all class types
   * 
   * @see #enableAllClasses()
   * @see #getClassCapabilities()
   */
  public void disableAllClasses() {
    for (Capability cap: Capability.values()) {
      if (cap.isClass())
	disable(cap);
    }
  }
  
  /**
   * disables all class type dependencies
   * 
   * @see #enableAllClassDependencies()
   * @see #getClassCapabilities()
   */
  public void disableAllClassDependencies() {
    for (Capability cap: Capability.values()) {
      if (cap.isClass())
	disableDependency(cap);
    }
  }
  
  /**
   * disables all attribute types
   * 
   * @see #enableAllAttributes()
   * @see #getAttributeCapabilities()
   */
  public void disableAllAttributes() {
    for (Capability cap: Capability.values()) {
      if (cap.isAttribute())
	disable(cap);
    }
  }
  
  /**
   * disables all attribute type dependencies
   * 
   * @see #enableAllAttributeDependencies()
   * @see #getAttributeCapabilities()
   */
  public void disableAllAttributeDependencies() {
    for (Capability cap: Capability.values()) {
      if (cap.isAttribute())
	disableDependency(cap);
    }
  }
  
  /**
   * disables all attribute and class types (including dependencies)
   */
  public void disableAll() {
    disableAllAttributes();
    disableAllAttributeDependencies();
    disableAllClasses();
    disableAllClassDependencies();
    disable(Capability.MISSING_VALUES);
    disable(Capability.MISSING_CLASS_VALUES);
    disable(Capability.NO_CLASS);
    //RANKING BEGIN
    disable(Capability.RANKING);
    disable(Capability.PREFERENCE_ATTRIBUTE);
    //RANKING END
  }
  
  /**
   * returns all class capabilities
   * 
   * @return		all capabilities regarding the class
   * @see #enableAllClasses()
   * @see #disableAllClasses()
   */
  public Capabilities getClassCapabilities() {
    Capabilities	result;
    
    result = new Capabilities(getOwner());
    
    for (Capability cap: Capability.values()) {
      if (cap.isClassCapability()) {
	if (handles(cap))
	  result.m_Capabilities.add(cap);
      }
    }
    
    return result;
  }
  
  /**
   * returns all attribute capabilities
   * 
   * @return		all capabilities regarding attributes
   * @see #enableAllAttributes()
   * @see #disableAllAttributes()
   */
  public Capabilities getAttributeCapabilities() {
    Capabilities	result;
    
    result = new Capabilities(getOwner());
    
    for (Capability cap: Capability.values()) {
      if (cap.isAttributeCapability()) {
	if (handles(cap))
	  result.m_Capabilities.add(cap);
      }
    }
    
    return result;
  }
  
  /**
   * returns all other capabilities, besides class and attribute related ones
   * 
   * @return		all other capabilities, besides class and attribute 
   * 			related ones
   */
  public Capabilities getOtherCapabilities() {
    Capabilities	result;
    
    result = new Capabilities(getOwner());
    
    for (Capability cap: Capability.values()) {
      if (cap.isOtherCapability()) {
	if (handles(cap))
	  result.m_Capabilities.add(cap);
      }
    }
    
    return result;
  }

  /**
   * returns true if the classifier handler has the specified capability
   *
   * @param c     the capability to test
   * @return      true if the classifier handler has the capability
   */
  public boolean handles(Capability c) {
    return m_Capabilities.contains(c);
  }

  /**
   * returns true if the classifier handler has a dependency for the specified 
   * capability
   *
   * @param c     the capability to test
   * @return      true if the classifier handler has a dependency for the 
   *               capability
   */
  public boolean hasDependency(Capability c) {
    return m_Dependencies.contains(c);
  }
  
  /**
   * Checks whether there are any dependencies at all
   * 
   * @return true if there is at least one dependency for a capability
   */
  public boolean hasDependencies() {
    return (m_Dependencies.size() > 0);
  }

  /**
   * returns the reason why the tests failed, is null if tests succeeded
   * 
   * @return		the reason why the tests failed
   */
  public Exception getFailReason() {
    return m_FailReason;
  }
  
  /**
   * Generates the message for, e.g., an exception. Adds the classname before the
   * actual message and returns that string.
   * 
   * @param msg		the actual content of the message, e.g., exception
   * @return		the new message
   */
  protected String createMessage(String msg) {
    String	result;
    
    result = "";
    
    if (getOwner() != null)
      result = getOwner().getClass().getName();
    else
      result = "<anonymous>";
      
    result += ": " + msg;
    
    return result;
  }
  
  /**
   * Test the given attribute, whether it can be processed by the handler,
   * given its capabilities. The method assumes that the specified attribute
   * is not the class attribute.
   * 
   * @param att		the attribute to test
   * @return		true if all the tests succeeded
   * @see		#test(Attribute, boolean)
   */
  public boolean test(Attribute att) {
    return test(att, false);
  }
  
  /**
   * Test the given attribute, whether it can be processed by the handler,
   * given its capabilities.
   * 
   * @param att		the attribute to test
   * @param isClass	whether this attribute is the class attribute
   * @return		true if all the tests succeeded
   * @see		#m_AttributeTest
   */
  public boolean test(Attribute att, boolean isClass) {
    boolean		result;
    Capability		cap;
    Capability		capBinary;
    Capability		capUnary;
    Capability		capEmpty;
    Capability      capRanking;
    String		errorStr;
    
    result = true;
    
    // shall we test the data?
    if (!m_AttributeTest)
      return result;

    // for exception
    if (isClass)
      errorStr  = "class";
    else
      errorStr  = "attributes";
    
    switch (att.type()) {
    
    //RANKING BEGIN
    case PreferenceAttribute.RANKING:
    	if (isClass) {
    	  cap       = Capability.RANKING;
    	  capBinary = Capability.BINARY_CLASS;
    	  capUnary  = Capability.UNARY_CLASS;
    	  capEmpty  = Capability.EMPTY_NOMINAL_CLASS;
    	}
    	else {
    	  cap       = Capability.PREFERENCE_ATTRIBUTE;
    	  capBinary = Capability.BINARY_ATTRIBUTES;
    	  capUnary  = Capability.UNARY_ATTRIBUTES;
    	  capEmpty  = Capability.EMPTY_NOMINAL_ATTRIBUTES;
    	}
    	
            if (handles(cap) && (att.numValues() > 2))
              break;
            else if (handles(capBinary) && (att.numValues() == 2))
              break;
            else if (handles(capUnary) && (att.numValues() == 1))
              break;
            else if (handles(capEmpty) && (att.numValues() == 0))
              break;

            if (att.numValues() == 0) {
              m_FailReason = new UnsupportedAttributeTypeException(
                  createMessage("Cannot handle empty nominal " + errorStr + "!"));
              result = false;
            }
            if (att.numValues() == 1) {
              m_FailReason = new UnsupportedAttributeTypeException(
                  createMessage("Cannot handle unary " + errorStr + "!"));
              result = false;
            }
            else if (att.numValues() == 2) {
              m_FailReason = new UnsupportedAttributeTypeException(
                  createMessage("Cannot handle binary " + errorStr + "!"));
              result = false;
            }
            else {
              m_FailReason = new UnsupportedAttributeTypeException(
                  createMessage("Cannot handle multi-valued nominal " + errorStr + "!"));
              result = false;
            }
            break;
    //RANKING END
    
      
      case Attribute.NOMINAL:
	if (isClass) {
	  cap       = Capability.NOMINAL_CLASS;
	  capBinary = Capability.BINARY_CLASS;
	  capUnary  = Capability.UNARY_CLASS;
	  capEmpty  = Capability.EMPTY_NOMINAL_CLASS;
	}
	else {
	  cap       = Capability.NOMINAL_ATTRIBUTES;
	  capBinary = Capability.BINARY_ATTRIBUTES;
	  capUnary  = Capability.UNARY_ATTRIBUTES;
	  capEmpty  = Capability.EMPTY_NOMINAL_ATTRIBUTES;
	}
	
        if (handles(cap) && (att.numValues() > 2))
          break;
        else if (handles(capBinary) && (att.numValues() == 2))
          break;
        else if (handles(capUnary) && (att.numValues() == 1))
          break;
        else if (handles(capEmpty) && (att.numValues() == 0))
          break;

        if (att.numValues() == 0) {
          m_FailReason = new UnsupportedAttributeTypeException(
              createMessage("Cannot handle empty nominal " + errorStr + "!"));
          result = false;
        }
        if (att.numValues() == 1) {
          m_FailReason = new UnsupportedAttributeTypeException(
              createMessage("Cannot handle unary " + errorStr + "!"));
          result = false;
        }
        else if (att.numValues() == 2) {
          m_FailReason = new UnsupportedAttributeTypeException(
              createMessage("Cannot handle binary " + errorStr + "!"));
          result = false;
        }
        else {
          m_FailReason = new UnsupportedAttributeTypeException(
              createMessage("Cannot handle multi-valued nominal " + errorStr + "!"));
          result = false;
        }
        break;

      case Attribute.NUMERIC:
	if (isClass)
	  cap = Capability.NUMERIC_CLASS;
	else
	  cap = Capability.NUMERIC_ATTRIBUTES;
	
        if (!handles(cap)) {
          m_FailReason = new UnsupportedAttributeTypeException(
                              createMessage("Cannot handle numeric " + errorStr + "!"));
          result = false;
        }
        break;

      case Attribute.DATE:
	if (isClass)
	  cap = Capability.DATE_CLASS;
	else
	  cap = Capability.DATE_ATTRIBUTES;
	
        if (!handles(cap)) {
          m_FailReason = new UnsupportedAttributeTypeException(
                              createMessage("Cannot handle date " + errorStr + "!"));
          result = false;
        }
        break;

      case Attribute.STRING:
	if (isClass)
	  cap = Capability.STRING_CLASS;
	else
	  cap = Capability.STRING_ATTRIBUTES;
	
        if (!handles(cap)) {
          m_FailReason = new UnsupportedAttributeTypeException(
                              createMessage("Cannot handle string " + errorStr + "!"));
          result = false;
        }
        break;

      case Attribute.RELATIONAL:
	if (isClass)
	  cap = Capability.RELATIONAL_CLASS;
	else
	  cap = Capability.RELATIONAL_ATTRIBUTES;
	
        if (!handles(cap)) {
          m_FailReason = new UnsupportedAttributeTypeException(
                              createMessage("Cannot handle relational " + errorStr + "!"));
          result = false;
        }
        // attributes in the relation of this attribute must be tested
        // separately with a different Capabilites object
        break;

      default:
        m_FailReason = new UnsupportedAttributeTypeException(
                            createMessage("Cannot handle unknown attribute type '" 
                                        + att.type() + "'!"));
        result = false;
    }
    
    return result;
  }
  
  /**
   * Tests the given data, whether it can be processed by the handler,
   * given its capabilities. Classifiers implementing the 
   * <code>MultiInstanceCapabilitiesHandler</code> interface are checked 
   * automatically for their multi-instance Capabilities (if no bags, then
   * only the bag-structure, otherwise only the first bag).
   *
   * @param data 	the data to test
   * @return		true if all the tests succeeded
   * @see 		#test(Instances, int, int)
   */
  public boolean test(Instances data) {
    return test(data, 0, data.numAttributes() - 1);
  }
  
  /**
   * Tests a certain range of attributes of the given data, whether it can be 
   * processed by the handler, given its capabilities. Classifiers 
   * implementing the <code>MultiInstanceCapabilitiesHandler</code> interface 
   * are checked automatically for their multi-instance Capabilities (if no 
   * bags, then only the bag-structure, otherwise only the first bag).
   *
   * @param data 	the data to test
   * @param fromIndex	the range of attributes - start (incl.)
   * @param toIndex	the range of attributes - end (incl.)
   * @return		true if all the tests succeeded
   * @see 		MultiInstanceCapabilitiesHandler
   * @see 		#m_InstancesTest
   * @see		#m_MissingValuesTest
   * @see		#m_MissingClassValuesTest
   * @see		#m_MinimumNumberInstancesTest
   */
  public boolean test(Instances data, int fromIndex, int toIndex) {
    int         	i;
    int         	n;
    int			m;
    Attribute   	att;
    Instance    	inst;
    boolean		testClass;
    Capabilities	cap;
    boolean		missing;
    Iterator		iter;
    
    // shall we test the data?
    if (!m_InstancesTest)
      return true;
    
    // no Capabilities? -> warning
    if (    (m_Capabilities.size() == 0) 
	 || ((m_Capabilities.size() == 1) && handles(Capability.NO_CLASS)) )
      System.err.println(createMessage("No capabilities set!"));
    
    // any attributes?
    if (toIndex - fromIndex < 0) {
      m_FailReason = new WekaException(
                          createMessage("No attributes!"));
      return false;
    }

    // do wee need to test the class attribute, i.e., is the class attribute
    // within the range of attributes?
    testClass =    (data.classIndex() > -1) 
    		&& (data.classIndex() >= fromIndex)
    		&& (data.classIndex() <= toIndex);
    
    // attributes
    for (i = fromIndex; i <= toIndex; i++) {
      att = data.attribute(i);
      
      // class is handled separately
      if (i == data.classIndex())
        continue;
      
      // check attribute types
      if (!test(att))
	return false;
    }

    // class
    if (!handles(Capability.NO_CLASS) && (data.classIndex() == -1)) {
      m_FailReason = new UnassignedClassException(
	  createMessage("Class attribute not set!"));
      return false;
    }
      
    // special case: no class attribute can be handled
    if (handles(Capability.NO_CLASS) && (data.classIndex() > -1)) {
      cap  = getClassCapabilities();
      cap.disable(Capability.NO_CLASS);
      iter = cap.capabilities();
      if (!iter.hasNext()) {
	m_FailReason = new WekaException(
	    createMessage("Cannot handle any class attribute!"));
	return false;
      }
    }
      
    if (testClass && !handles(Capability.NO_CLASS)) {
      att = data.classAttribute();
      if (!test(att, true))
	return false;

      // special handling of RELATIONAL class
      // TODO: store additional Capabilities for this case
      
      // missing class labels
      if (m_MissingClassValuesTest) {
	if (!handles(Capability.MISSING_CLASS_VALUES)) {
	  for (i = 0; i < data.numInstances(); i++) {
	    if (data.instance(i).classIsMissing()) {
	      m_FailReason = new WekaException(
		  createMessage("Cannot handle missing class values!"));
	      return false;
	    }
	  }
	}
	else {
	  if (m_MinimumNumberInstancesTest) {
	    int hasClass = 0;
	    
	    for (i = 0; i < data.numInstances(); i++) {
	      if (!data.instance(i).classIsMissing())
		hasClass++;
	    }
	    
	    // not enough instances with class labels?
	    if (hasClass < getMinimumNumberInstances()) {
	      m_FailReason = new WekaException(
		  createMessage("Not enough training instances with class labels (required: " 
		      + getMinimumNumberInstances() 
		      + ", provided: " 
		      + hasClass + ")!"));
	      return false;
	    }
	  }
	}
      }
    }

    // missing values
    if (m_MissingValuesTest) {
      if (!handles(Capability.MISSING_VALUES)) {
	missing = false;
	for (i = 0; i < data.numInstances(); i++) {
	  inst = data.instance(i);
	  
	  if (inst instanceof SparseInstance) {
	    for (m = 0; m < inst.numValues(); m++) {
	      n = inst.index(m);
	      
	      // out of scope?
	      if (n < fromIndex)
		continue;
	      if (n > toIndex)
		break;

	      // skip class
	      if (n == inst.classIndex())
		continue;
	      

	      if (inst.isMissing(n)) {
		missing = true;
		break;
	      }
	    }
	  }
	  else {
	    for (n = fromIndex; n <= toIndex; n++) {
	      // skip class
	      if (n == inst.classIndex())
		continue;

	      if (inst.isMissing(n)) {
		missing = true;
		break;
	      }
	    }
	  }
	  
	  if (missing) {
	    m_FailReason = new NoSupportForMissingValuesException(
		createMessage("Cannot handle missing values!"));
	    return false;
	  }
	}
      }
    }
    
    // instances
    if (m_MinimumNumberInstancesTest) {
      if (data.numInstances() < getMinimumNumberInstances()) {
	m_FailReason = new WekaException(
	    createMessage("Not enough training instances (required: " 
		+ getMinimumNumberInstances() 
		+ ", provided: " 
		+ data.numInstances() + ")!"));
	return false;
      }
    }

    // Multi-Instance? -> check structure (regardless of attribute range!)
    if (handles(Capability.ONLY_MULTIINSTANCE)) {
      // number of attributes?
      if (data.numAttributes() != 3) {
        m_FailReason = new WekaException(
                            createMessage("Incorrect Multi-Instance format, must be 'bag-id, bag, class'!"));
        return false;
      }
      
      // type of attributes and position of class?
      if (    !data.attribute(0).isNominal()  || !data.attribute(0).isRanking() 
           || !data.attribute(1).isRelationValued() 
           || (data.classIndex() != data.numAttributes() - 1) ) {
        m_FailReason = new WekaException(
            createMessage("Incorrect Multi-Instance format, must be 'NOMINAL att, RELATIONAL att, CLASS att'!"));
        return false;
      }

      // check data immediately
      if (getOwner() instanceof MultiInstanceCapabilitiesHandler) {
	MultiInstanceCapabilitiesHandler handler = (MultiInstanceCapabilitiesHandler) getOwner();
	cap = handler.getMultiInstanceCapabilities();
	boolean result;
	if (data.numInstances() > 0)
	  result = cap.test(data.attribute(1).relation(0));
	else
	  result = cap.test(data.attribute(1).relation());
	
	if (!result) {
	  m_FailReason = cap.m_FailReason;
	  return false;
	}
      }
    }
    
    // passed all tests!
    return true;
  }

  /**
   * tests the given attribute by calling the test(Attribute,boolean) method 
   * and throws an exception if the test fails. The method assumes that the
   * specified attribute is not the class attribute.
   *
   * @param att        	the attribute to test
   * @throws Exception  in case the attribute doesn't pass the tests
   * @see 		#test(Attribute,boolean)
   */
  public void testWithFail(Attribute att) throws Exception {
    test(att, false);
  }

  /**
   * tests the given attribute by calling the test(Attribute,boolean) method 
   * and throws an exception if the test fails.
   *
   * @param att        	the attribute to test
   * @param isClass	whether this attribute is the class attribute
   * @throws Exception  in case the attribute doesn't pass the tests
   * @see 		#test(Attribute,boolean)
   */
  public void testWithFail(Attribute att, boolean isClass) throws Exception {
    if (!test(att, isClass))
      throw m_FailReason;
  }

  /**
   * tests the given data by calling the test(Instances,int,int) method and 
   * throws an exception if the test fails.
   *
   * @param data        the data to test
   * @param fromIndex	the range of attributes - start (incl.)
   * @param toIndex	the range of attributes - end (incl.)
   * @throws Exception  in case the data doesn't pass the tests
   * @see 		#test(Instances,int,int)
   */
  public void testWithFail(Instances data, int fromIndex, int toIndex) throws Exception {
    if (!test(data, fromIndex, toIndex))
      throw m_FailReason;
  }

  /**
   * tests the given data by calling the test(Instances) method and throws 
   * an exception if the test fails.
   *
   * @param data        the data to test
   * @throws Exception  in case the data doesn't pass the tests
   * @see 		#test(Instances)
   */
  public void testWithFail(Instances data) throws Exception {
    if (!test(data))
      throw m_FailReason;
  }
  
  /**
   * returns a string representation of the capabilities
   * 
   * @return 	a string representation of this object
   */
  public String toString() {
    Vector<Capability>		sorted;
    StringBuffer	result;
    
    result = new StringBuffer();

    // capabilities
    sorted = new Vector<Capability>(m_Capabilities);
    Collections.sort(sorted);
    result.append("Capabilities: " + sorted.toString() + "\n");

    // dependencies
    sorted = new Vector<Capability>(m_Dependencies);
    Collections.sort(sorted);
    result.append("Dependencies: " + sorted.toString() + "\n");
    
    // other stuff
    result.append("min # Instance: " + getMinimumNumberInstances() + "\n");
    
    return result.toString();
  }
  
  /**
   * turns the capabilities object into source code. The returned source code
   * is a block that creates a Capabilities object named 'objectname' and
   * enables all the capabilities of this Capabilities object.
   * 
   * @param objectname	the name of the Capabilities object being instantiated
   * @return		the generated source code
   */
  public String toSource(String objectname) {
    return toSource(objectname, 0);
  }
    
  /**
   * turns the capabilities object into source code. The returned source code
   * is a block that creates a Capabilities object named 'objectname' and
   * enables all the capabilities of this Capabilities object.
   * 
   * @param objectname	the name of the Capabilities object being instantiated
   * @param indent	the number of blanks to indent
   * @return		the generated source code
   */
  public String toSource(String objectname, int indent) {
    StringBuffer	result;
    String		capsName;
    String		capName;
    String		indentStr;
    int			i;
    
    result = new StringBuffer();

    capsName = Capabilities.class.getName();
    capName  = Capabilities.Capability.class.getName().replaceAll("\\$", ".");
    
    indentStr = "";
    for (i = 0; i < indent; i++)
      indentStr += " ";
    
    // object name
    result.append(indentStr + capsName + " " + objectname + " = new " + capsName + "(this);\n");
    
    List<Capability> capsList = new ArrayList<Capability>();
    boolean hasNominalAtt = false;
    boolean hasBinaryAtt = false;
    boolean hasUnaryAtt = false;
    boolean hasEmptyNomAtt = false;
    boolean hasNominalClass = false;
    boolean hasRankingClass = false;
    boolean hasPreferenceAtt = false;
    // capabilities
    result.append("\n");
    for (Capability cap: Capability.values()) {
      // capability
      if (handles(cap)) {
          //RANKING BEGIN
        if (cap == Capability.RANKING){
          hasRankingClass = true;
        }
        if(cap == Capability.PREFERENCE_ATTRIBUTE){
          hasPreferenceAtt = true;
        }
          //RANKING END
        if (cap == Capability.NOMINAL_ATTRIBUTES) {
          hasNominalAtt = true;
        }
        if (cap == Capability.NOMINAL_CLASS) {
          hasNominalClass = true;
        }
        if (cap == Capability.BINARY_ATTRIBUTES) {
          hasBinaryAtt = true;
        }
        if (cap == Capability.UNARY_ATTRIBUTES) {
          hasUnaryAtt = true;
        }
        if (cap == Capability.EMPTY_NOMINAL_ATTRIBUTES) {
          hasEmptyNomAtt = true;
        }

        capsList.add(cap);              
      }
    }
    
    for (Capability cap : capsList) {
      if ((cap == Capability.BINARY_ATTRIBUTES && hasNominalAtt) ||
          (cap == Capability.UNARY_ATTRIBUTES && hasBinaryAtt) ||
          (cap == Capability.EMPTY_NOMINAL_ATTRIBUTES && hasUnaryAtt) ||
          (cap == Capability.BINARY_CLASS && hasNominalClass) ||
          (cap == Capability.RANKING && hasRankingClass) ||
          (cap == Capability.PREFERENCE_ATTRIBUTE && hasPreferenceAtt)){
        continue;
      }
      result.append(
          indentStr + objectname + ".enable(" + capName + "." + cap.name() + ");\n");
      // dependency
      if (hasDependency(cap))
        result.append(
            indentStr + objectname + ".enableDependency(" + capName + "." + cap.name() + ");\n");
    }
    
    // capabilities
    result.append("\n");

    // other
    result.append("\n");
    result.append(
	indentStr + objectname + ".setMinimumNumberInstances(" 
	+ getMinimumNumberInstances() + ");\n");

    result.append("\n");
    
    return result.toString();
  }
  
  /**
   * returns a Capabilities object specific for this data. The multi-instance
   * capability is not checked as well as the minimum number of instances
   * is not set.
   * 
   * @param data	the data to base the capabilities on
   * @return		a data-specific capabilities object
   * @throws Exception	in case an error occurrs, e.g., an unknown attribute 
   * 			type
   */
  public static Capabilities forInstances(Instances data) throws Exception {
    return forInstances(data, false);
  }
  
  /**
   * returns a Capabilities object specific for this data. The minimum number 
   * of instances is not set, the check for multi-instance data is optional.
   * 
   * @param data	the data to base the capabilities on
   * @param multi	if true then the structure is checked, too
   * @return		a data-specific capabilities object
   * @throws Exception	in case an error occurrs, e.g., an unknown attribute 
   * 			type
   */
  public static Capabilities forInstances(Instances data, boolean multi) throws Exception {
    Capabilities	result;
    Capabilities	multiInstance;
    int			i;
    int			n;
    int			m;
    Instance		inst;
    boolean		missing;
    
    result = new Capabilities(null);
    
    // class
    if (data.classIndex() == -1) {
      result.enable(Capability.NO_CLASS);
    }
    else {
      switch (data.classAttribute().type()) {
      
    case PreferenceAttribute.RANKING:
    	result.enable(Capability.RANKING);
    	//result.enable(Capability.NOMINAL_CLASS);
    	break;
      
	case Attribute.NOMINAL:
	  if (data.classAttribute().numValues() == 1)
	    result.enable(Capability.UNARY_CLASS);
	  else if (data.classAttribute().numValues() == 2)
	    result.enable(Capability.BINARY_CLASS);
	  else
	    result.enable(Capability.NOMINAL_CLASS);
	  break;
	 	  
	case Attribute.NUMERIC:
	  result.enable(Capability.NUMERIC_CLASS);
	  break;
	  
	case Attribute.STRING:
	  result.enable(Capability.STRING_CLASS);
	  break;
	  
	case Attribute.DATE:
	  result.enable(Capability.DATE_CLASS);
	  break;
	  
	case Attribute.RELATIONAL:
	  result.enable(Capability.RELATIONAL_CLASS);
	  break;
	  
	
	  
	default:
	  throw new UnsupportedAttributeTypeException(
	      "Unknown class attribute type '" + data.classAttribute() + "'!");
      }
      
      // missing class values
      for (i = 0; i < data.numInstances(); i++) {
	if (data.instance(i).classIsMissing()) {
	  result.enable(Capability.MISSING_CLASS_VALUES);
	  break;
	}
      }
    }
    
    // attributes
    for (i = 0; i < data.numAttributes(); i++) {
      // skip class
      if (i == data.classIndex())
	continue;

      switch (data.attribute(i).type()) {
      //RANKING BEGIN
    case PreferenceAttribute.RANKING:
    	result.enable(Capability.PREFERENCE_ATTRIBUTE);
    	break;
      //RANKING END
	case Attribute.NOMINAL:
	  result.enable(Capability.UNARY_ATTRIBUTES);
	  if (data.attribute(i).numValues() == 2)
	    result.enable(Capability.BINARY_ATTRIBUTES);
	  else if (data.attribute(i).numValues() > 2)
	    result.enable(Capability.NOMINAL_ATTRIBUTES);
	  break;
	
	case Attribute.NUMERIC:
	  result.enable(Capability.NUMERIC_ATTRIBUTES);
	  break;
		
	case Attribute.DATE:
	  result.enable(Capability.DATE_ATTRIBUTES);
	  break;

	case Attribute.STRING:
	  result.enable(Capability.STRING_ATTRIBUTES);
	  break;
	  
	case Attribute.RELATIONAL:
	  result.enable(Capability.RELATIONAL_ATTRIBUTES);
	  break;
	
	
	  
	default:
	  throw new UnsupportedAttributeTypeException(
	      "Unknown attribute type '" + data.attribute(i).type() + "'!");
      }
    }
    
    // missing values
    missing = false;
    for (i = 0; i < data.numInstances(); i++) {
      inst = data.instance(i);

      if (inst instanceof SparseInstance) {
	for (m = 0; m < inst.numValues(); m++) {
	  n = inst.index(m);

	  // skip class
	  if (n == inst.classIndex())
	    continue;

	  if (inst.isMissing(n)) {
	    missing = true;
	    break;
	  }
	}
      }
      else {
	for (n = 0; n < data.numAttributes(); n++) {
	  // skip class
	  if (n == inst.classIndex())
	    continue;

	  if (inst.isMissing(n)) {
	    missing = true;
	    break;
	  }
	}
      }

      if (missing) {
	result.enable(Capability.MISSING_VALUES);
	break;
      }
    }

    // multi-instance data?
    if (multi) {
      if (    (data.numAttributes() == 3)
	   && ((data.attribute(0).isNominal()) || data.attribute(0).isRanking())		// bag-id
	   && (data.attribute(1).isRelationValued()) 	// bag
	   && (data.classIndex() == data.numAttributes() - 1) ) {
	multiInstance = new Capabilities(null);
	multiInstance.or(result.getClassCapabilities());
	multiInstance.enable(Capability.NOMINAL_ATTRIBUTES);
	multiInstance.enable(Capability.RELATIONAL_ATTRIBUTES);
	multiInstance.enable(Capability.ONLY_MULTIINSTANCE);
	result.assign(multiInstance);
      }
    }
    
    return result;
  }
  
  /**
   * loads the given dataset and prints the Capabilities necessary to 
   * process it. <p/>
   * 
   * Valid parameters: <p/>
   * 
   * -file filename <br/>
   *  the file to load
   *  
   * -c index
   *  the explicit index of the class attribute (default: none)
   * 
   * @param args	the commandline arguments
   * @throws Exception	if something goes wrong
   */
  public static void main(String[] args) throws Exception {
    String 		tmpStr;
    String		filename;
    DataSource 		source;
    Instances 		data;
    int 		classIndex;
    Capabilities 	cap;
    Iterator		iter;

    if (args.length == 0) {
      System.out.println(
	  "\nUsage: " + Capabilities.class.getName() 
	  + " -file <dataset> [-c <class index>]\n");
      return;
    }
    
    // get parameters
    tmpStr = Utils.getOption("file", args);
    if (tmpStr.length() == 0)
      throw new Exception("No file provided with option '-file'!");
    else
      filename = tmpStr;

    tmpStr = Utils.getOption("c", args);
    if (tmpStr.length() != 0) {
      if (tmpStr.equals("first"))
	classIndex = 0;
      else if (tmpStr.equals("last"))
	classIndex = -2;  // last
      else
	classIndex = Integer.parseInt(tmpStr) - 1;
    }
    else {
      classIndex = -3;  // not set
    }
    
    // load data
    source = new DataSource(filename);
    if (classIndex == -3)
      data = source.getDataSet();
    else if (classIndex == -2)
      data = source.getDataSet(source.getStructure().numAttributes() - 1);
    else
      data = source.getDataSet(classIndex);

    // determine and print capabilities
    cap = forInstances(data);
    System.out.println("File: " + filename);
    System.out.println("Class index: " + ((data.classIndex() == -1) ? "not set" : "" + (data.classIndex() + 1)));
    System.out.println("Capabilities:");
    iter = cap.capabilities();
    while (iter.hasNext())
      System.out.println("- " + iter.next());
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6447 $");
  }
  
}
