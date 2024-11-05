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
 * XMLSerialization.java
 * Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core.xml;

import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.core.Version;

import java.beans.BeanInfo;
import java.beans.Introspector;
import java.beans.PropertyDescriptor;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStream;
import java.io.Reader;
import java.io.Writer;
import java.lang.reflect.Array;
import java.lang.reflect.Constructor;
import java.lang.reflect.Method;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Vector;

import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * With this class objects can be serialized to XML instead into a binary 
 * format. It uses introspection (cf. beans) to retrieve the data from the
 * given object, i.e. it can only access beans-conform fields automatically.
 * <p>
 * The generic approach of writing data as XML can be overriden by adding 
 * custom methods for reading/writing in a derived class
 * (cf. <code>m_Properties</code>, <code>m_CustomMethods</code>).<br>
 * Custom read and write methods must have the same signature (and also be 
 * <code>public</code>!) as the <code>readFromXML</code> and <code>writeToXML</code>
 * methods. Methods that apply to the naming rule <code>read + property name</code>
 * are added automatically to the list of methods by the method 
 * <code>XMLSerializationMethodHandler.addMethods(...)</code>.  
 * <p>
 * Other properties that are not conform the bean set/get-methods have to be 
 * processed manually in a derived class (cf. <code>readPostProcess(Object)</code>, 
 * <code>writePostProcess(Object)</code>).
 * <p>
 * For a complete XML serialization/deserialization have a look at the 
 * <code>KOML</code> class.
 * <p>
 * If a stored class has a constructor that takes a String to initialize
 * (e.g. String or Double) then the content of the tag will used for the
 * constructor, e.g. from 
 * <pre>&lt;object name="name" class="String" primitive="no"&gt;Smith&lt;/object&gt;</pre>
 * "Smith" will be used to instantiate a String object as constructor argument.
 * <p>   
 * 
 * @see KOML
 * @see #fromXML(Document)
 * @see #toXML(Object)
 * @see #m_Properties
 * @see #m_CustomMethods
 * @see #readPostProcess(Object)
 * @see #writePostProcess(Object)
 * @see #readFromXML(Element)
 * @see #writeToXML(Element, Object, String)
 * 
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5953 $ 
 */
public class XMLSerialization
   implements RevisionHandler {
  
   /** for debugging purposes only */
   protected static boolean DEBUG = false;
  
   /** the node that is currently processed, in case of writing the parent node
    * (something might go wrong writing the new child) and in case of reading 
    * the actual node that is tried to process */
   protected Element m_CurrentNode = null;
   
   /** the tag for an object */
   public final static String TAG_OBJECT = "object";
   
   /** the version attribute */
   public final static String ATT_VERSION = XMLDocument.ATT_VERSION;
  
   /** the tag for the name */
   public final static String ATT_NAME = XMLDocument.ATT_NAME;
   
   /** the tag for the class */
   public final static String ATT_CLASS = "class";
   
   /** the tag whether primitive or not (yes/no) */
   public final static String ATT_PRIMITIVE = "primitive";
   
   /** the tag whether array or not (yes/no) */
   public final static String ATT_ARRAY = "array";
   
   /** the tag whether null or not (yes/no) */
   public final static String ATT_NULL = "null";
   
   /** the value "yes" for the primitive and array attribute */
   public final static String VAL_YES = XMLDocument.VAL_YES;
   
   /** the value "no" for the primitive and array attribute */
   public final static String VAL_NO = XMLDocument.VAL_NO;
   
   /** the value of the name for the root node */
   public final static String VAL_ROOT = "__root__";
   
   /** the root node of the XML document */
   public final static String ROOT_NODE = TAG_OBJECT; 
   
   /** default value for attribute ATT_PRIMITIVE
    * @see #ATT_PRIMITIVE */
   public final static String ATT_PRIMITIVE_DEFAULT = VAL_NO;
   
   /** default value for attribute ATT_ARRAY
    * @see #ATT_ARRAY */
   public final static String ATT_ARRAY_DEFAULT = VAL_NO;
   
   /** default value for attribute ATT_NULL
    * @see #ATT_NULL */
   public final static String ATT_NULL_DEFAULT = VAL_NO;
   
   /** the DOCTYPE for the serialization */
   public final static String DOCTYPE = 
        "<!" + XMLDocument.DTD_DOCTYPE + " " + ROOT_NODE + "\n"
      + "[\n"
      + "   <!" + XMLDocument.DTD_ELEMENT + " " + TAG_OBJECT + " (" + XMLDocument.DTD_PCDATA + XMLDocument.DTD_SEPARATOR + TAG_OBJECT + ")" + XMLDocument.DTD_ZERO_OR_MORE + ">\n"
      + "   <!" + XMLDocument.DTD_ATTLIST + " " + TAG_OBJECT + " " + ATT_NAME + "      " + XMLDocument.DTD_CDATA + " " + XMLDocument.DTD_REQUIRED + ">\n"
      + "   <!" + XMLDocument.DTD_ATTLIST + " " + TAG_OBJECT + " " + ATT_CLASS + "     " + XMLDocument.DTD_CDATA + " " + XMLDocument.DTD_REQUIRED + ">\n"
      + "   <!" + XMLDocument.DTD_ATTLIST + " " + TAG_OBJECT + " " + ATT_PRIMITIVE + " " + XMLDocument.DTD_CDATA + " \"" + ATT_PRIMITIVE_DEFAULT + "\">\n"
      + "   <!" + XMLDocument.DTD_ATTLIST + " " + TAG_OBJECT + " " + ATT_ARRAY + "     " + XMLDocument.DTD_CDATA + " \"" + ATT_ARRAY_DEFAULT + "\">   <!-- the dimensions of the array; no=0, yes=1 -->\n"
      + "   <!" + XMLDocument.DTD_ATTLIST + " " + TAG_OBJECT + " " + ATT_NULL + "      " + XMLDocument.DTD_CDATA + " \"" + ATT_NULL_DEFAULT + "\">\n"
      + "   <!" + XMLDocument.DTD_ATTLIST + " " + TAG_OBJECT + " " + ATT_VERSION + "   " + XMLDocument.DTD_CDATA + " \"" + Version.VERSION + "\">\n"
      + "]\n"
      + ">";
   
   /** the XMLDocument that performs the transformation to and fro XML */
   protected XMLDocument m_Document = null;
   
   /** for handling properties (ignored/allowed) */
   protected PropertyHandler m_Properties = null;
   
   /** for handling custom read/write methods */
   protected XMLSerializationMethodHandler m_CustomMethods = null;

   /** for overriding class names (Class &lt;-&gt; Classname (String)) 
    * @see #overrideClassname(Object) */
   protected Hashtable<Class,String> m_ClassnameOverride = null;
   
   /**
    * initializes the serialization
    * 
    * @throws Exception if initialization fails
    */
   public XMLSerialization() throws Exception {
      super();
      clear();
   }
   
   /**
    * used for debugging purposes, i.e. only if DEBUG is set to true.
    * needs a newly generated Throwable instance to get the method/line from
    * @param t      a throwable instance, generated in the calling method
    * @param msg    a message to pring
    * @see          #DEBUG
    */
   protected void trace(Throwable t, String msg) {
     if ( (DEBUG) && (t.getStackTrace().length > 0) ) {
       System.out.println("trace: " + t.getStackTrace()[0] + ": " + msg);
     }
   }
   
   /**
    * generates internally a new XML document and clears also the IgnoreList and
    * the mappings for the Read/Write-Methods
    * 
    * @throws Exception	if something goes wrong
    */
   public void clear() throws Exception {
      m_Document = new XMLDocument();
      m_Document.setValidating(true);
      m_Document.newDocument(DOCTYPE, ROOT_NODE);
      
      m_Properties        = new PropertyHandler();
      m_CustomMethods     = new XMLSerializationMethodHandler(this);

      m_ClassnameOverride = new Hashtable<Class,String>();
      // java.io.File is sometimes represented as another class:
      // - Win32: sun.awt.shell.Win32ShellFolder2 
      // - Linux: sun.awt.shell.DefaultShellFolder
      // -> we set it to "java.io.File"
      m_ClassnameOverride.put(java.io.File.class, java.io.File.class.getName());
      
      setVersion(Version.VERSION); 
      
      m_CurrentNode = null;
   }
   
   /**
    * sets the given version string in the XML document
    * 
    * @param version	the new version string
    */
   private void setVersion(String version) {
      Document     doc;
      
      doc = m_Document.getDocument();
      doc.getDocumentElement().setAttribute(ATT_VERSION, version);
   }
   
   /**
    * returns the WEKA version with which the serialized object was created
    * 
    * @return		the current version
    * @see Version 
    */
   public String getVersion() {
      Document     doc;
      String       result;
      
      doc    = m_Document.getDocument();
      result = doc.getDocumentElement().getAttribute(ATT_VERSION);
      
      return result;
   }
   
   /**
    * Checks the version in the current Document with the one of the current
    * release. If the version differ, a warning is printed.
    */
   private void checkVersion() {
      String            versionStr;
      Version           version;
      
      version    = new Version();
      versionStr = getVersion();
      if (versionStr.equals(""))
         System.out.println("WARNING: has no version!");
      else if (version.isOlder(versionStr))
         System.out.println("WARNING: loading a newer version (" + versionStr + " > " + Version.VERSION + ")!");
      else if (version.isNewer(versionStr))
         System.out.println("NOTE: loading an older version (" + versionStr + " < " + Version.VERSION + ")!");
   }
   
   /**
    * returns a hashtable with PropertyDescriptors that have "get" and "set" 
    * methods indexed by the property name.
    * 
    * @see java.beans.PropertyDescriptor
    * @param o the object to retrieve the descriptors from
    * @return the PropertyDescriptors indexed by name of the property
    * @throws Exception if the introspection fails
    */
   protected Hashtable getDescriptors(Object o) throws Exception {
      BeanInfo                   info;
      PropertyDescriptor[]       desc;
      int                        i;
      Hashtable<String,PropertyDescriptor>                  result;
      
      result = new Hashtable<String,PropertyDescriptor>();

      info = Introspector.getBeanInfo(o.getClass());
      desc = info.getPropertyDescriptors();
      for (i = 0; i < desc.length; i++) {
         // get AND set method?
         if ( (desc[i].getReadMethod() != null) && (desc[i].getWriteMethod() != null) ) {
            // in ignore list, i.e. a general ignore without complete path?
            if (m_Properties.isIgnored(desc[i].getDisplayName()))
               continue;
            
            // in ignore list of the class?
            if (m_Properties.isIgnored(o, desc[i].getDisplayName()))
               continue;
            
            // not an allowed property
            if (!m_Properties.isAllowed(o, desc[i].getDisplayName()))
               continue;
            
            result.put(desc[i].getDisplayName(), desc[i]);
         }
      }
      
      return result;
   }

   /**
    * returns the path of the "name" attribute from the root down to this node
    * (including it).
    * 
    * @param node the node to get the path for
    * @return the complete "name" path of this node
    */
   protected String getPath(Element node) {
      String            result;
      
      result = node.getAttribute(ATT_NAME);

      while (node.getParentNode() != node.getOwnerDocument()) {
         node   = (Element) node.getParentNode(); 
         result = node.getAttribute(ATT_NAME) + "." + result;
      }
      
      return result;
   }
   
   /**
    * returns either <code>VAL_YES</code> or <code>VAL_NO</code> depending 
    * on the value of <code>b</code>
    * 
    * @param b the boolean to turn into a string
    * @return the value in string representation
    */
   protected String booleanToString(boolean b) {
      if (b)
         return VAL_YES;
      else
         return VAL_NO;
   }
   
   /**
    * turns the given string into a boolean, if a positive number is given, 
    * then zero is considered FALSE, every other number TRUE; the empty string 
    * is also considered being FALSE
    * 
    * @param s the string to turn into a boolean
    * @return the string as boolean
    */
   protected boolean stringToBoolean(String s) {
      if (s.equals(""))
         return false;
      else if (s.equals(VAL_YES))
         return true;
      else if (s.equalsIgnoreCase("true"))
         return true;
      else if (s.replaceAll("[0-9]*", "").equals(""))
         return (Integer.parseInt(s) != 0);
      else
         return false;
   }
   
   /**
    * appends a new node to the parent with the given parameters (a non-array)
    * 
    * @param parent the parent of this node. if it is <code>null</code> the 
    *        document root element is used
    * @param name the name of the node
    * @param classname the classname for this node
    * @param primitive whether it is a primitve data type or not (i.e. an object)
    * @return the generated node 
    */
   protected Element addElement(Element parent, String name, String classname, boolean primitive) {
     return addElement(parent, name, classname, primitive, 0);
   }
   
   /**
    * appends a new node to the parent with the given parameters
    * 
    * @param parent the parent of this node. if it is <code>null</code> the 
    *        document root element is used
    * @param name the name of the node
    * @param classname the classname for this node
    * @param primitive whether it is a primitve data type or not (i.e. an object)
    * @param array the dimensions of the array (0 if not an array)
    * @return the generated node 
    */
   protected Element addElement(Element parent, String name, String classname, boolean primitive, int array) {
     return addElement(parent, name, classname, primitive, array, false);
   }
   
   /**
    * appends a new node to the parent with the given parameters
    * 
    * @param parent the parent of this node. if it is <code>null</code> the 
    *        document root element is used
    * @param name the name of the node
    * @param classname the classname for this node
    * @param primitive whether it is a primitve data type or not (i.e. an object)
    * @param array the dimensions of the array (0 if not an array)
    * @param isnull whether it is null
    * @return the generated node 
    */
   protected Element addElement(Element parent, String name, String classname, boolean primitive, int array, boolean isnull) {
      Element           result;

      if (parent == null)
         result = m_Document.getDocument().getDocumentElement();
      else
         result = (Element) parent.appendChild(m_Document.getDocument().createElement(TAG_OBJECT));
      
      // attributes
      // mandatory attributes:
      result.setAttribute(ATT_NAME, name);
      result.setAttribute(ATT_CLASS, classname);
      
      // add following attributes only if necessary, i.e., different from default:
      if (!booleanToString(primitive).equals(ATT_PRIMITIVE_DEFAULT))
        result.setAttribute(ATT_PRIMITIVE, booleanToString(primitive));

      // multi-dimensional array?
      if (array > 1) {
        result.setAttribute(ATT_ARRAY, Integer.toString(array));
      }
      // backwards compatible: 0 -> no array ("no"), 1 -> 1-dim. array ("yes")
      else {
        if (!booleanToString(array == 1).equals(ATT_ARRAY_DEFAULT))
          result.setAttribute(ATT_ARRAY, booleanToString(array == 1));
      }
      
      if (!booleanToString(isnull).equals(ATT_NULL_DEFAULT))
        result.setAttribute(ATT_NULL, booleanToString(isnull));
      
      return result;
   }
   
   /**
    * if the class of the given object (or one of its ancestors) is stored in 
    * the classname override hashtable, then the override name is returned 
    * otherwise the classname of the given object.
    * 
    * @param o          the object to check for overriding its classname
    * @return           if overridden then the classname stored in the hashtable,
    *                   otherwise the classname of the given object
    * @see              #m_ClassnameOverride
    */
   protected String overrideClassname(Object o) {
      Enumeration    enm;
      String         result;
      Class          currentCls;
     
      result = o.getClass().getName();

      // check overrides
      enm    = m_ClassnameOverride.keys();
      while (enm.hasMoreElements()) {
         currentCls = (Class) enm.nextElement();
         if (currentCls.isInstance(o)) {
           result = (String) m_ClassnameOverride.get(currentCls);
           break;
         }
      }
      
      return result;
   }
   
   /**
    * if the given classname is stored in the classname override hashtable, 
    * then the override name is returned otherwise the given classname.
    * <b>Note:</b> in contrast to <code>overrideClassname(Object)</code> does
    * this method only look for exact name matches. The other method checks
    * whether the class of the given object is a subclass of any of the stored
    * overrides.  
    * 
    * @param classname  the classname to check for overriding
    * @return           if overridden then the classname stored in the hashtable,
    *                   otherwise the given classname
    * @see              #m_ClassnameOverride
    * @see              #overrideClassname(Object)
    */
   protected String overrideClassname(String classname) {
      Enumeration    enm;
      String         result;
      Class          currentCls;
     
      result = classname;

      // check overrides
      enm    = m_ClassnameOverride.keys();
      while (enm.hasMoreElements()) {
         currentCls = (Class) enm.nextElement();
         if (currentCls.getName().equals(classname)) {
           result = (String) m_ClassnameOverride.get(currentCls);
           break;
         }
      }
      
      return result;
   }
   
   /**
    * returns a property descriptor if possible, otherwise <code>null</code>
    * 
    * @param className the name of the class to get the descriptor for
    * @param displayName the name of the property
    * @return the descriptor if available, otherwise <code>null</code>
    */
   protected PropertyDescriptor determineDescriptor(String className, String displayName) {
      PropertyDescriptor      result;
      
      result = null;
      
      try {
         result = new PropertyDescriptor(displayName, Class.forName(className));
      }
      catch (Exception e) {
         result = null;
      }
      
      return result;
   }
   
   /**
    * adds the given primitive to the DOM structure.
    * @param parent the parent of this object, e.g. the class this object is a member of
    * @param o the primitive to describe in XML
    * @param name the name of the primitive
    * @return the node that was created
    * @throws Exception if the DOM creation fails
    */
   protected Element writeBooleanToXML(Element parent, boolean o, String name) throws Exception {
     Element      node;

     // for debugging only
     if (DEBUG)
        trace(new Throwable(), name);
     
     m_CurrentNode = parent;
     
     node = addElement(parent, name, Boolean.TYPE.getName(), true);
     node.appendChild(node.getOwnerDocument().createTextNode(new Boolean(o).toString()));
     
     return node;
   }
   
   /**
    * adds the given primitive to the DOM structure.
    * @param parent the parent of this object, e.g. the class this object is a member of
    * @param o the primitive to describe in XML
    * @param name the name of the primitive
    * @return the node that was created
    * @throws Exception if the DOM creation fails
    */
   protected Element writeByteToXML(Element parent, byte o, String name) throws Exception {
     Element      node;
     
     // for debugging only
     if (DEBUG)
        trace(new Throwable(), name);
     
     m_CurrentNode = parent;
     
     node = addElement(parent, name, Byte.TYPE.getName(), true);
     node.appendChild(node.getOwnerDocument().createTextNode(new Byte(o).toString()));
     
     return node;
   }
   
   /**
    * adds the given primitive to the DOM structure.
    * @param parent the parent of this object, e.g. the class this object is a member of
    * @param o the primitive to describe in XML
    * @param name the name of the primitive
    * @return the node that was created
    * @throws Exception if the DOM creation fails
    */
   protected Element writeCharToXML(Element parent, char o, String name) throws Exception {
     Element      node;
     
     // for debugging only
     if (DEBUG)
        trace(new Throwable(), name);
     
     m_CurrentNode = parent;
     
     node = addElement(parent, name, Character.TYPE.getName(), true);
     node.appendChild(node.getOwnerDocument().createTextNode(new Character(o).toString()));
     
     return node;
   }
   
   /**
    * adds the given primitive to the DOM structure.
    * @param parent the parent of this object, e.g. the class this object is a member of
    * @param o the primitive to describe in XML
    * @param name the name of the primitive
    * @return the node that was created
    * @throws Exception if the DOM creation fails
    */
   protected Element writeDoubleToXML(Element parent, double o, String name) throws Exception {
     Element      node;
     
     // for debugging only
     if (DEBUG)
        trace(new Throwable(), name);
     
     m_CurrentNode = parent;
     
     node = addElement(parent, name, Double.TYPE.getName(), true);
     node.appendChild(node.getOwnerDocument().createTextNode(new Double(o).toString()));
     
     return node;
   }
   
   /**
    * adds the given primitive to the DOM structure.
    * @param parent the parent of this object, e.g. the class this object is a member of
    * @param o the primitive to describe in XML
    * @param name the name of the primitive
    * @return the node that was created
    * @throws Exception if the DOM creation fails
    */
   protected Element writeFloatToXML(Element parent, float o, String name) throws Exception {
     Element      node;
     
     // for debugging only
     if (DEBUG)
        trace(new Throwable(), name);
     
     m_CurrentNode = parent;
     
     node = addElement(parent, name, Float.TYPE.getName(), true);
     node.appendChild(node.getOwnerDocument().createTextNode(new Float(o).toString()));
     
     return node;
   }
   
   /**
    * adds the given primitive to the DOM structure.
    * @param parent the parent of this object, e.g. the class this object is a member of
    * @param o the primitive to describe in XML
    * @param name the name of the primitive
    * @return the node that was created
    * @throws Exception if the DOM creation fails
    */
   protected Element writeIntToXML(Element parent, int o, String name) throws Exception {
     Element      node;
     
     // for debugging only
     if (DEBUG)
        trace(new Throwable(), name);
     
     m_CurrentNode = parent;
     
     node = addElement(parent, name, Integer.TYPE.getName(), true);
     node.appendChild(node.getOwnerDocument().createTextNode(new Integer(o).toString()));
     
     return node;
   }
   
   /**
    * adds the given primitive to the DOM structure.
    * @param parent the parent of this object, e.g. the class this object is a member of
    * @param o the primitive to describe in XML
    * @param name the name of the primitive
    * @return the node that was created
    * @throws Exception if the DOM creation fails
    */
   protected Element writeLongToXML(Element parent, long o, String name) throws Exception {
     Element      node;
     
     // for debugging only
     if (DEBUG)
        trace(new Throwable(), name);
     
     m_CurrentNode = parent;

     node = addElement(parent, name, Long.TYPE.getName(), true);
     node.appendChild(node.getOwnerDocument().createTextNode(new Long(o).toString()));
     
     return node;
   }
   
   /**
    * adds the given primitive to the DOM structure.
    * @param parent the parent of this object, e.g. the class this object is a member of
    * @param o the primitive to describe in XML
    * @param name the name of the primitive
    * @return the node that was created
    * @throws Exception if the DOM creation fails
    */
   protected Element writeShortToXML(Element parent, short o, String name) throws Exception {
     Element      node;
     
     // for debugging only
     if (DEBUG)
        trace(new Throwable(), name);
     
     m_CurrentNode = parent;
     
     node = addElement(parent, name, Short.TYPE.getName(), true);
     node.appendChild(node.getOwnerDocument().createTextNode(new Short(o).toString()));
     
     return node;
   }
   
   /**
    * checks whether the innermost class is a primitive class (handles 
    * multi-dimensional arrays)
    * @param c        the array class to inspect
    * @return         whether the array consists of primitive elements
    */
   protected boolean isPrimitiveArray(Class c) {
     if (c.getComponentType().isArray())
       return isPrimitiveArray(c.getComponentType());
    else
       return c.getComponentType().isPrimitive();
   }
   
   /**
    * adds the given Object to a DOM structure. 
    * (only public due to reflection).<br>
    * <b>Note:</b> <code>overrideClassname(Object)</code> is not invoked in case of
    * arrays, since the array class could be a superclass, whereas the elements of
    * the array can be specialized subclasses. In case of an array the method 
    * <code>overrideClassname(String)</code> is invoked, which searches for an 
    * exact match of the classname in the override hashtable.
    * 
    * @param parent the parent of this object, e.g. the class this object is a member of
    * @param o the Object to describe in XML
    * @param name the name of the object
    * @return the node that was created
    * @throws Exception if the DOM creation fails
    * @see #overrideClassname(Object)
    * @see #overrideClassname(String)
    * @see #m_ClassnameOverride
    */
   public Element writeToXML(Element parent, Object o, String name) throws Exception {
      String               classname;
      Element              node;
      Hashtable            memberlist;
      Enumeration          enm;
      Object               member;
      String               memberName;
      Method               method;
      PropertyDescriptor   desc;
      boolean              primitive;
      int                  array;
      int                  i;
      Object               obj;
      String               tmpStr;

      node = null;
      
      // for debugging only
      if (DEBUG)
         trace(new Throwable(), name);

      // special handling of null-objects
      if (o == null) {
        node = addElement(parent, name, "" + null, false, 0, true);
        return node;
      }
      
      // used for overriding the classname
      obj = null;
      
      // get information about object
      array = 0;
      if (o.getClass().isArray())
        array = Utils.getArrayDimensions(o);
      if (array > 0) {
        classname = Utils.getArrayClass(o.getClass()).getName();
        primitive = isPrimitiveArray(o.getClass()); 
      }
      else {
         // try to get property descriptor to determine real class
         // (for primitives the getClass() method returns the corresponding Object-Class!)
         desc = null;
         if (parent != null)
            desc = determineDescriptor(parent.getAttribute(ATT_CLASS), name);
         
         if (desc != null)
            primitive = desc.getPropertyType().isPrimitive(); 
         else
            primitive = o.getClass().isPrimitive(); 

         // for primitives: retrieve primitive type, otherwise the object's real 
         // class. For non-primitives we can't use the descriptor, since that
         // might only return an interface as class!
         if (primitive) {
            classname = desc.getPropertyType().getName();
         }
         else {
            obj       = o;
            classname = o.getClass().getName();
         }
      }
      
      // fix class/primitive if parent is array of primitives, thanks to 
      // reflection the elements of the array are objects and not primitives!
      if (    (parent != null) 
           && (!parent.getAttribute(ATT_ARRAY).equals(""))
           && (!parent.getAttribute(ATT_ARRAY).equals(VAL_NO))
           && (stringToBoolean(parent.getAttribute(ATT_PRIMITIVE))) ) {
         primitive = true;
         classname = parent.getAttribute(ATT_CLASS);
         obj       = null;
      }

      // perhaps we need to override the classname
      if (obj != null)
        classname = overrideClassname(obj);         // for non-arrays
      else
        classname = overrideClassname(classname);   // for arrays
      
      // create node for current object
      node = addElement(parent, name, classname, primitive, array);
      
      // array? -> save as child with 'name="<index>"'
      if (array > 0) {
         for (i = 0; i < Array.getLength(o); i++) {
            invokeWriteToXML(node, Array.get(o, i), Integer.toString(i));
         }
      }
      // non-array
      else {
         // primitive? -> only toString()
         if (primitive) {
            node.appendChild(node.getOwnerDocument().createTextNode(o.toString()));
         }
         // object
         else {
            // process recursively members of this object 
            memberlist = getDescriptors(o);
            // if no get/set methods -> we assume it has String-Constructor
            if (memberlist.size() == 0) {
              if (!o.toString().equals("")) {
        	tmpStr = o.toString();
        	// these five entities are recognized by every XML processor
        	// see http://www.xml.com/pub/a/2001/03/14/trxml10.html
        	tmpStr = tmpStr.replaceAll("&", "&amp;")
        	               .replaceAll("\"", "&quot;")
        	               .replaceAll("'", "&apos;")
        	               .replaceAll("<", "&lt;")
        	               .replaceAll(">", "&gt;");
        	// in addition, replace some other entities as well
        	tmpStr = tmpStr.replaceAll("\n", "&#10;")
        	               .replaceAll("\r", "&#13;")
        	               .replaceAll("\t", "&#9;");
        	node.appendChild(node.getOwnerDocument().createTextNode(tmpStr));
              }
            }
            else {
               enm = memberlist.keys();
               while (enm.hasMoreElements()) {
                  memberName = enm.nextElement().toString();
                  
                  // in ignore list?
                  if (    (m_Properties.isIgnored(memberName))
                       || (m_Properties.isIgnored(getPath(node) + "." + memberName))
                       || (m_Properties.isIgnored(o, getPath(node) + "." + memberName)) )
                     continue;

                  // is it allowed?
                  if (!m_Properties.isAllowed(o, memberName))
                     continue;
                  
                  desc   = (PropertyDescriptor) memberlist.get(memberName);
                  method = desc.getReadMethod();
                  member = method.invoke(o, (Object[]) null);
                  invokeWriteToXML(node, member, memberName);
               }
            }
         }
      }
      
      return node;
   }
   
   /**
    * either invokes a custom method to write a specific property/class or the standard
    * method <code>writeToXML(Element,Object,String)</code>
    * 
    * @param parent the parent XML node
    * @param o the object's content will be added as children to the given parent node
    * @param name the name of the object
    * @return the node that was created
    * @throws Exception if invocation or turning into XML fails
    */
   protected Element invokeWriteToXML(Element parent, Object o, String name) throws Exception {
      Method         method;
      Class[]        methodClasses;
      Object[]       methodArgs;
      boolean        array;
      Element        node;
      boolean        useDefault;
      
      node       = null;
      method     = null;
      useDefault = false;

      m_CurrentNode = parent;
      
      // default, if null
      if (o == null)
         useDefault = true;
      
      try {
         if (!useDefault) {
            array = o.getClass().isArray();
           
            // display name?
            if (m_CustomMethods.write().contains(name))
               method = (Method) m_CustomMethods.write().get(o.getClass());
            else
            // class?
            if ( (!array) && (m_CustomMethods.write().contains(o.getClass())) )
               method = (Method) m_CustomMethods.write().get(o.getClass());
            else
               method = null;
            
            useDefault = (method == null);
         }

         // custom
         if (!useDefault) {
             methodClasses    = new Class[3];
             methodClasses[0] = Element.class;
             methodClasses[1] = Object.class;
             methodClasses[2] = String.class;
             methodArgs       = new Object[3];
             methodArgs[0]    = parent;
             methodArgs[1]    = o;
             methodArgs[2]    = name;
             node = (Element) method.invoke(this, methodArgs);
         }
         // standard
         else {
            node = writeToXML(parent, o, name);
         }
      }
      catch (Exception e) {
         if (DEBUG)
            e.printStackTrace();
         
         if (m_CurrentNode != null) {
           System.out.println("Happened near: " + getPath(m_CurrentNode));
           // print it only once!
           m_CurrentNode = null;
         }
         System.out.println("PROBLEM (write): " + name);

         throw (Exception) e.fillInStackTrace();
      }
      
      return node;
   }
   
   /**
    * enables derived classes to due some pre-processing on the objects, that's
    * about to be serialized. Right now it only returns the object.
    * 
    * @param o the object that is serialized into XML
    * @return the possibly altered object
    * @throws Exception if post-processing fails
    */
   protected Object writePreProcess(Object o) throws Exception {
     return o;
   }
   
   /**
    * enables derived classes to add other properties to the DOM tree, e.g.
    * ones that do not apply to the get/set convention of beans. only implemented
    * with empty method body.
    * 
    * @param o the object that is serialized into XML
    * @throws Exception if post-processing fails
    */
   protected void writePostProcess(Object o) throws Exception {
   }
   
   /**
    * extracts all accesible properties from the given object
    * 
    * @param o the object to turn into an XML representation
    * @return the generated DOM document 
    * @throws Exception if XML generation fails 
    */
   public XMLDocument toXML(Object o) throws Exception {
      clear();
      invokeWriteToXML(null, writePreProcess(o), VAL_ROOT);
      writePostProcess(o);
      return m_Document;
   }
   
   /**
    * returns a descriptor for a given objet by providing the name
    * 
    * @param o the object the get the descriptor for
    * @param name the display name of the descriptor
    * @return the Descriptor, if found, otherwise <code>null</code>
    * @throws Exception if introsepction fails 
    */
   protected PropertyDescriptor getDescriptorByName(Object o, String name) throws Exception {
      PropertyDescriptor      result;
      PropertyDescriptor[]    desc;
      int                     i;
      
      result = null;
      
      desc   = Introspector.getBeanInfo(o.getClass()).getPropertyDescriptors();
      for (i = 0; i < desc.length; i++) {
         if (desc[i].getDisplayName().equals(name)) {
            result = desc[i];
            break;
         }
      }
      
      return result;
   }
   
   /**
    * returns the associated class for the given name
    * 
    * @param name the name of the class to return a Class object for
    * @return the class if  it could be retrieved
    * @throws Exception if it class retrieval fails
    */
   protected Class determineClass(String name) throws Exception {
      Class       result;
      
      if (name.equals(Boolean.TYPE.getName()))
         result = Boolean.TYPE;
      else
      if (name.equals(Byte.TYPE.getName()))
         result = Byte.TYPE;
      else
      if (name.equals(Character.TYPE.getName()))
         result = Character.TYPE;
      else
      if (name.equals(Double.TYPE.getName()))
         result = Double.TYPE;
      else
      if (name.equals(Float.TYPE.getName()))
         result = Float.TYPE;
      else
      if (name.equals(Integer.TYPE.getName()))
         result = Integer.TYPE;
      else
      if (name.equals(Long.TYPE.getName()))
         result = Long.TYPE;
      else
      if (name.equals(Short.TYPE.getName()))
         result = Short.TYPE;
      else
         result = Class.forName(name);
      
      return result;
   }
   
   /**
    * returns an Object representing the primitive described by the given node.
    * Here we use a trick to return an object even though its a primitive: by 
    * creating a primitive array with reflection of length 1, setting the 
    * primtive value as real object and then returning the "object" at 
    * position 1 of the array.
    * 
    * @param node the node to return the value as "primitive" object
    * @return the primitive as "pseudo" object
    * @throws Exception if the instantiation of the array fails or any of the
    *         String conversions fails
    */
   protected Object getPrimitive(Element node) throws Exception {
      Object            result;
      Object            tmpResult;
      Class             cls;
      
      cls       = determineClass(node.getAttribute(ATT_CLASS));
      tmpResult = Array.newInstance(cls, 1);
      
      if (cls == Boolean.TYPE)
         Array.set(tmpResult, 0, new Boolean(XMLDocument.getContent(node)));
      else
      if (cls == Byte.TYPE)
         Array.set(tmpResult, 0, new Byte(XMLDocument.getContent(node)));
      else
      if (cls == Character.TYPE)
         Array.set(tmpResult, 0, new Character(XMLDocument.getContent(node).charAt(0)));
      else
      if (cls == Double.TYPE)
         Array.set(tmpResult, 0, new Double(XMLDocument.getContent(node)));
      else
      if (cls == Float.TYPE)
         Array.set(tmpResult, 0, new Float(XMLDocument.getContent(node)));
      else
      if (cls == Integer.TYPE)
         Array.set(tmpResult, 0, new Integer(XMLDocument.getContent(node)));
      else
      if (cls == Long.TYPE)
         Array.set(tmpResult, 0, new Long(XMLDocument.getContent(node)));
      else
      if (cls == Short.TYPE)
         Array.set(tmpResult, 0, new Short(XMLDocument.getContent(node)));
      else
         throw new Exception("Cannot get primitive for class '" + cls.getName() + "'!");
      
      result = Array.get(tmpResult, 0);
      
      return result;
   }
   
   /**
    * builds the primitive from the given DOM node. 
    * 
    * @param node the associated XML node
    * @return the primitive created from the XML description
    * @throws Exception if instantiation fails 
    */
   public boolean readBooleanFromXML(Element node) throws Exception {
     // for debugging only
     if (DEBUG)
        trace(new Throwable(), node.getAttribute(ATT_NAME));

     m_CurrentNode = node;
     
     return ((Boolean) getPrimitive(node)).booleanValue();
   }
   
   /**
    * builds the primitive from the given DOM node. 
    * 
    * @param node the associated XML node
    * @return the primitive created from the XML description
    * @throws Exception if instantiation fails 
    */
   public byte readByteFromXML(Element node) throws Exception {
     // for debugging only
     if (DEBUG)
        trace(new Throwable(), node.getAttribute(ATT_NAME));

     m_CurrentNode = node;
     
     return ((Byte) getPrimitive(node)).byteValue();
   }
   
   /**
    * builds the primitive from the given DOM node. 
    * 
    * @param node the associated XML node
    * @return the primitive created from the XML description
    * @throws Exception if instantiation fails 
    */
   public char readCharFromXML(Element node) throws Exception {
     // for debugging only
     if (DEBUG)
        trace(new Throwable(), node.getAttribute(ATT_NAME));

     m_CurrentNode = node;
     
     return ((Character) getPrimitive(node)).charValue();
   }
   
   /**
    * builds the primitive from the given DOM node. 
    * 
    * @param node the associated XML node
    * @return the primitive created from the XML description
    * @throws Exception if instantiation fails 
    */
   public double readDoubleFromXML(Element node) throws Exception {
     // for debugging only
     if (DEBUG)
        trace(new Throwable(), node.getAttribute(ATT_NAME));

     m_CurrentNode = node;
     
     return ((Double) getPrimitive(node)).doubleValue();
   }
   
   /**
    * builds the primitive from the given DOM node. 
    * 
    * @param node the associated XML node
    * @return the primitive created from the XML description
    * @throws Exception if instantiation fails 
    */
   public float readFloatFromXML(Element node) throws Exception {
     // for debugging only
     if (DEBUG)
        trace(new Throwable(), node.getAttribute(ATT_NAME));

     m_CurrentNode = node;
     
     return ((Float) getPrimitive(node)).floatValue();
   }
   
   /**
    * builds the primitive from the given DOM node. 
    * 
    * @param node the associated XML node
    * @return the primitive created from the XML description
    * @throws Exception if instantiation fails 
    */
   public int readIntFromXML(Element node) throws Exception {
     // for debugging only
     if (DEBUG)
        trace(new Throwable(), node.getAttribute(ATT_NAME));

     m_CurrentNode = node;
     
     return ((Integer) getPrimitive(node)).intValue();
   }
   
   /**
    * builds the primitive from the given DOM node. 
    * 
    * @param node the associated XML node
    * @return the primitive created from the XML description
    * @throws Exception if instantiation fails 
    */
   public long readLongFromXML(Element node) throws Exception {
     // for debugging only
     if (DEBUG)
        trace(new Throwable(), node.getAttribute(ATT_NAME));

     m_CurrentNode = node;
     
     return ((Long) getPrimitive(node)).longValue();
   }
   
   /**
    * builds the primitive from the given DOM node. 
    * 
    * @param node the associated XML node
    * @return the primitive created from the XML description
    * @throws Exception if instantiation fails 
    */
   public short readShortFromXML(Element node) throws Exception {
     // for debugging only
     if (DEBUG)
        trace(new Throwable(), node.getAttribute(ATT_NAME));

     m_CurrentNode = node;
     
     return ((Short) getPrimitive(node)).shortValue();
   }
   
   /**
    * adds the specific node to the object via a set method
    * 
    * @param o            the object to set a property
    * @param name         the name of the object for which to set a property
    *                     (only for information reasons)
    * @param child        the value of the property to add
    * @return             the provided object, but augmented by the child
    * @throws Exception   if something goes wrong
    */
   public Object readFromXML(Object o, String name, Element child) throws Exception {
      Object               result;
      Hashtable            descriptors;
      PropertyDescriptor   descriptor;
      String               methodName;
      Method               method;
      Object[]             methodArgs;
      Object               tmpResult;
      Class                paramClass;
     
      result      = o;
      descriptors = getDescriptors(result);
      methodName  = child.getAttribute(ATT_NAME);

      // in ignore list?
      if (m_Properties.isIgnored(getPath(child)))
         return result;
      
      // in ignore list of class?
      if (m_Properties.isIgnored(result, getPath(child)))
        return result;
      
      // is it allowed?
      if (!m_Properties.isAllowed(result, methodName))
        return result;
      
      descriptor = (PropertyDescriptor) descriptors.get(methodName);

      // unknown property?
      if (descriptor == null) {
         if (!m_CustomMethods.read().contains(methodName))
            System.out.println("WARNING: unknown property '" + name + "." + methodName + "'!");
         return result;
      }
      
      method     = descriptor.getWriteMethod();
      methodArgs = new Object[1];
      tmpResult  = invokeReadFromXML(child);
      paramClass = method.getParameterTypes()[0];
      
      // array?
      if (paramClass.isArray()) {
         // no data?
         if (Array.getLength(tmpResult) == 0)
           return result;
         methodArgs[0] = (Object[]) tmpResult;
      }
      // non-array
      else {
         methodArgs[0] = tmpResult;
      }

      method.invoke(result, methodArgs);
     
      return result;
   }
   
   /**
    * returns an array with the dimensions of the array stored in XML
    * @param node the node to determine the dimensions for
    * @return the dimensions of the array
    */
   protected int[] getArrayDimensions(Element node) {
     Vector<Element>         children;
     Vector<Integer>         tmpVector;      
     int[]          tmp;
     int[]          result;
     int            i;
     
     // have we reached the innermost dimension?
     if (stringToBoolean(node.getAttribute(ATT_ARRAY)))
       children = XMLDocument.getChildTags(node);
     else
       children = null;
     
     if (children != null) {
       tmpVector = new Vector<Integer>();

       if (children.size() > 0) {
         // are children also arrays?
         tmp = getArrayDimensions((Element) children.get(0));
         
         // further dimensions
         if (tmp != null) {
           for (i = tmp.length - 1; i >= 0; i--)
             tmpVector.add(new Integer(tmp[i]));
         }
  
         // add current dimension
         tmpVector.add(0, new Integer(children.size()));
       }
       else {
         tmpVector.add(new Integer(0));
       }
       
       // generate result
       result = new int[tmpVector.size()];
       for (i = 0; i < result.length; i++)
         result[i] = ((Integer) tmpVector.get(tmpVector.size() - i - 1)).intValue();
     }
     else {
       result = null;
     }
     
     return result;
   }
   
   /**
    * builds the object from the given DOM node. 
    * (only public due to reflection) 
    * 
    * @param node the associated XML node
    * @return the instance created from the XML description
    * @throws Exception if instantiation fails 
    */
   public Object readFromXML(Element node) throws Exception {
      String               classname;
      String               name;
      boolean              primitive;
      boolean              array;
      boolean              isnull;
      Class<?>                cls;
      Vector<Element>               children;
      Object               result;
      int                  i;
      Constructor          constructor;
      Class[]              methodClasses;
      Object[]             methodArgs;
      Element              child;
           
      // for debugging only
      if (DEBUG)
         trace(new Throwable(), node.getAttribute(ATT_NAME));

      m_CurrentNode = node;
      
      result    = null;
      
      name      = node.getAttribute(ATT_NAME);
      classname = node.getAttribute(ATT_CLASS);
      primitive = stringToBoolean(node.getAttribute(ATT_PRIMITIVE));
      array     = stringToBoolean(node.getAttribute(ATT_ARRAY));
      isnull    = stringToBoolean(node.getAttribute(ATT_NULL));

      // special handling of null
      if (isnull)
        return result;

      children  = XMLDocument.getChildTags(node);
      cls       = determineClass(classname);
      
      // array
      if (array) {
         result = Array.newInstance(cls, getArrayDimensions(node));
         for (i = 0; i < children.size(); i++) {
            child = (Element) children.get(i);
            Array.set(result, Integer.parseInt(child.getAttribute(ATT_NAME)), invokeReadFromXML(child));
         }
      }
      // non-array
      else {
         // primitive/String-constructor
         if (children.size() == 0) {
            // primitive
            if (primitive) {
               result = getPrimitive(node);
            }
            // assumed String-constructor
            else {
               methodClasses    = new Class[1];
               methodClasses[0] = String.class;
               methodArgs       = new Object[1];
               methodArgs[0]    = XMLDocument.getContent(node);
               try {
                  constructor   = cls.getConstructor(methodClasses);
                  result        = constructor.newInstance(methodArgs);
               }
               catch (Exception e) {
                  // if it's not a class with String constructor, let's try standard constructor
                  try {
                     result = cls.newInstance();
                  }
                  catch (Exception e2) {
                     // sorry, can't instantiate!
                     result = null;
                     System.out.println("ERROR: Can't instantiate '" + classname + "'!");
                  }
               }
            }
         }
         // normal get/set methods
         else {
            result = cls.newInstance();
            for (i = 0; i < children.size(); i++)
              result = readFromXML(result, name, (Element) children.get(i));
         }
      }
            
      return result;
   }
   
   /**
    * either invokes a custom method to read a specific property/class or the standard
    * method <code>readFromXML(Element)</code>
    * 
    * @param node the associated XML node
    * @return the instance created from the XML description
    * @throws Exception if instantiation fails
    */
   protected Object invokeReadFromXML(Element node) throws Exception {
      Method         method;
      Class[]        methodClasses;
      Object[]       methodArgs;
      boolean        array;
      boolean        useDefault;

      useDefault = false;
      method     = null;
      m_CurrentNode = node;
      
      try {
         // special handling of null values
         if (stringToBoolean(node.getAttribute(ATT_NULL)))
           useDefault = true;
        
         if (!useDefault) {
            array = stringToBoolean(node.getAttribute(ATT_ARRAY));
           
            // display name?
            if (m_CustomMethods.read().contains(node.getAttribute(ATT_NAME)))
               method = (Method) m_CustomMethods.read().get(node.getAttribute(ATT_NAME));
            else
            // class name?
            if ( (!array) && (m_CustomMethods.read().contains(determineClass(node.getAttribute(ATT_CLASS)))) )
               method = (Method) m_CustomMethods.read().get(determineClass(node.getAttribute(ATT_CLASS)));
            else
               method = null;
            
            useDefault = (method == null);
         }

         // custom method
         if (!useDefault) {
            methodClasses    = new Class[1];
            methodClasses[0] = Element.class;
            methodArgs       = new Object[1];
            methodArgs[0]    = node;
            return method.invoke(this, methodArgs);
         }
         // standard
         else {
            return readFromXML(node);
         }
      }
      catch (Exception e) {
         if (DEBUG)
            e.printStackTrace();
         
         if (m_CurrentNode != null) {
           System.out.println("Happened near: " + getPath(m_CurrentNode));
           // print it only once!
           m_CurrentNode = null;
         }
         System.out.println("PROBLEM (read): " + node.getAttribute("name"));

         throw (Exception) e.fillInStackTrace();
      }
   }
   
   /**
    * additional pre-processing can happen in derived classes before the 
    * actual reading from XML (working on the raw XML). right now it does 
    * nothing with the document.
    * 
    * @param document 	the document to pre-process
    * @return the processed object
    * @throws Exception if post-processing fails
    */
   protected Document readPreProcess(Document document) throws Exception {
      return document;
   }
   
   /**
    * additional post-processing can happen in derived classes after reading 
    * from XML. right now it only returns the object as it is.
    * 
    * @param o the object to perform some additional processing on
    * @return the processed object
    * @throws Exception if post-processing fails
    */
   protected Object readPostProcess(Object o) throws Exception {
      return o;
   }

   /**
    * returns the given DOM document as an instance of the specified class
    * 
    * @param document the parsed DOM document representing the object
    * @return the XML as object 
    * @throws Exception if object instantiation fails
    */
   public Object fromXML(Document document) throws Exception {
      if (!document.getDocumentElement().getNodeName().equals(ROOT_NODE))
         throw new Exception("Expected '" + ROOT_NODE + "' as root element, but found '" + document.getDocumentElement().getNodeName() + "'!");
      m_Document.setDocument(readPreProcess(document));
      checkVersion();
      return readPostProcess(invokeReadFromXML(m_Document.getDocument().getDocumentElement()));
   }
   
   /**
    * parses the given XML string (can be XML or a filename) and returns an
    * Object generated from the representation
    * 
    * @param xml the xml to parse (if "<?xml" is not found then it is considered a file)
    * @return the generated instance
    * @throws Exception if something goes wrong with the parsing
    */
   public Object read(String xml) throws Exception {
      return fromXML(m_Document.read(xml));
   }
   
   /**
    * parses the given file and returns a DOM document
    * 
    * @param file the XML file to parse
    * @return the parsed DOM document
    * @throws Exception if something goes wrong with the parsing
    */
   public Object read(File file) throws Exception {
      return fromXML(m_Document.read(file));
   }
   
   /**
    * parses the given stream and returns a DOM document
    * 
    * @param stream the XML stream to parse
    * @return the parsed DOM document
    * @throws Exception if something goes wrong with the parsing
    */
   public Object read(InputStream stream) throws Exception {
      return fromXML(m_Document.read(stream));
   }
   
   /**
    * parses the given reader and returns a DOM document
    * 
    * @param reader the XML reader to parse
    * @return the parsed DOM document
    * @throws Exception if something goes wrong with the parsing
    */
   public Object read(Reader reader) throws Exception {
      return fromXML(m_Document.read(reader));
   }
   
   
   /**
    * writes the given object into the file
    * 
    * @param file the filename to write to
    * @param o the object to serialize as XML
    * @throws Exception if something goes wrong with the parsing
    */
   public void write(String file, Object o) throws Exception {
      toXML(o).write(file);
   }
   
   /**
    * writes the given object into the file
    * 
    * @param file the filename to write to
    * @param o the object to serialize as XML
    * @throws Exception if something goes wrong with the parsing
    */
   public void write(File file, Object o) throws Exception {
      toXML(o).write(file);
   }
   
   /**
    * writes the given object into the stream
    * 
    * @param stream the filename to write to
    * @param o the object to serialize as XML
    * @throws Exception if something goes wrong with the parsing
    */
   public void write(OutputStream stream, Object o) throws Exception {
      toXML(o).write(stream);
   }
   
   /**
    * writes the given object into the writer
    * 
    * @param writer the filename to write to
    * @param o the object to serialize as XML
    * @throws Exception if something goes wrong with the parsing
    */
   public void write(Writer writer, Object o) throws Exception {
      toXML(o).write(writer);
   }
   
   /**
    * for testing only. if the first argument is a filename with ".xml"
    * as extension it tries to generate an instance from the XML description
    * and does a <code>toString()</code> of the generated object.
    */
   public static void main(String[] args) throws Exception {
      if (args.length > 0) {
         // read xml and print
         if (args[0].toLowerCase().endsWith(".xml")) {
            System.out.println(new XMLSerialization().read(args[0]).toString());
         }
         // read binary and print generated XML
         else {
            // read
            FileInputStream fi = new FileInputStream(args[0]);
            ObjectInputStream oi = new ObjectInputStream(
                                   new BufferedInputStream(fi));
            Object o = oi.readObject();
            oi.close();
            // print to stdout
            //new XMLSerialization().write(System.out, o);
            new XMLSerialization().write(new BufferedOutputStream(new FileOutputStream(args[0] + ".xml")), o);
            // print to binary file
            FileOutputStream fo = new FileOutputStream(args[0] + ".exp");
            ObjectOutputStream oo = new ObjectOutputStream(
                                   new BufferedOutputStream(fo));
            oo.writeObject(o);
            oo.close();
         }
      }
   }
   
   /**
    * Returns the revision string.
    * 
    * @return		the revision
    */
   public String getRevision() {
     return RevisionUtils.extract("$Revision: 5953 $");
   }
}
