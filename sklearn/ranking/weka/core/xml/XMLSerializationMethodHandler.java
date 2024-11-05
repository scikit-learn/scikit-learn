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
 * XMLSerializationMethodHandler.java
 * Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core.xml;

import weka.core.RevisionHandler;
import weka.core.RevisionUtils;

import java.lang.reflect.Method;

import org.w3c.dom.Element;


/**
 * This class handles relationships between display names of properties 
 * (or classes) and Methods that are associated with them. It differentiates 
 * between read and write methods. It automatically stores public methods that 
 * have the same signature as the <code>readFromXML()</code> and 
 * <code>writeToXML()</code> methods in the <code>XMLSerialization</code>
 * class.  
 *  
 * @see MethodHandler
 * @see XMLSerialization
 * 
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5953 $ 
 */
public class XMLSerializationMethodHandler
   implements RevisionHandler {
  
   /** for storing read methods */
   protected MethodHandler m_ReadMethods = null;
   
   /** for storing write methods */
   protected MethodHandler m_WriteMethods = null;
   
   /** the object to retrieve the methods from */
   protected Object owner = null;
   
   /**
    * initializes the method handling, executes also <code>clear()</code>, which
    * adds initial methods automatically.
    * 
    * @param owner the owner to retrieve the methods from
    * @throws Exception if initialization fails
    * @see #clear() 
    */
   public XMLSerializationMethodHandler(Object owner) throws Exception {
      super();
      
      this.owner     = owner;
      m_ReadMethods  = new MethodHandler();
      m_WriteMethods = new MethodHandler();
      
      clear();
   }
   
   /**
    * adds all methods that are like <code>template</code> to the method list
    * 
    * @param handler the list to add fitting methods to
    * @param template the signature to check the given methods against
    * @param methods the methods to check
    */
   protected void addMethods(MethodHandler handler, Method template, Method[] methods) {
      int            i;
      int            n;
      Method         method;
      boolean        equal;
      String         name;
      
      for (i = 0; i < methods.length; i++) {
         method = methods[i];
         
         // is it template?
         if (template.equals(method))
            continue;
         
         // tests
         // 1. return type
         if (!template.getReturnType().equals(method.getReturnType()))
            continue;
            
         // 2. signature
         if (template.getParameterTypes().length != method.getParameterTypes().length)
            continue;
         
         equal = true;
         for (n = 0; n < template.getParameterTypes().length; n++) {
            if (!template.getParameterTypes()[n].equals(method.getParameterTypes()[n])) {
               equal = false;
               break;
            }
         }
            
         // add to list
         if (equal) {
            name = method.getName();
            name = name.replaceAll("read|write", "");
            name = name.substring(0, 1).toLowerCase() + name.substring(1);
            handler.add(name, method);
         }
      }
   }
   
   /**
    * automatically adds all fitting methods to the custom read/write lists,
    * it excludes only the generic ones. it is automatically called in 
    * <code>clear()</code>
    * It only work with methods that apply to the naming rule "read" + property
    * name (same for "write")
    *
    * @throws Exception if retrieving of methods fails
    * @see #clear()
    */
   protected void addMethods() throws Exception {
      Method         method;
      Class[]        params;
      
      // read
      params    = new Class[1];
      params[0] = Element.class;
      method    = owner.getClass().getMethod("readFromXML", params);
      addMethods(m_ReadMethods, method, owner.getClass().getMethods());
      
      // write
      params    = new Class[3];
      params[0] = Element.class;
      params[1] = Object.class;
      params[2] = String.class;
      method    = owner.getClass().getMethod("writeToXML", params);
      addMethods(m_WriteMethods, method, owner.getClass().getMethods());
   }
   
   /**
    * returns the method with the given name that has the same signature as
    * <code>readFromXML()</code> of the <code>XMLSerialiation</code> class. 
    * simplifies the adding of custom methods.
    * 
    * @param o the object to inspect
    * @param name the name of the method to return
    * @return either <code>null</code> if no method was found or a reference 
    * @see XMLSerialization#readFromXML(Element)
    */
   public static Method findReadMethod(Object o, String name) {
      Class[]        params;
      Method         result;
      
      result    = null;
      
      params    = new Class[1];
      params[0] = Element.class; 
      try {
         result = o.getClass().getMethod(name, params);
      }
      catch (Exception e) {
         result = null;  
      }
      
      return result;
   }
   
   /**
    * returns the method with the given name that has the same signature as
    * <code>writeToXML()</code> of the <code>XMLSerialiation</code> class. 
    * simplifies the adding of custom methods.
    * 
    * @param o the object to inspect
    * @param name the name of the method to return
    * @return either <code>null</code> if no method was found or a reference 
    * @see XMLSerialization#writeToXML(Element, Object, String)
    */
   public static Method findWriteMethod(Object o, String name) {
      Class[]        params;
      Method         result;
      
      result    = null;
      
      params    = new Class[3];
      params[0] = Element.class; 
      params[1] = Object.class; 
      params[2] = String.class; 
      try {
         result = o.getClass().getMethod(name, params);
      }
      catch (Exception e) {
         result = null;  
      }
      
      return result;
   }
   
   /**
    * removes all current methods and adds the methods according to the 
    *  
    */
   public void clear() {
      m_ReadMethods.clear();
      m_WriteMethods.clear();
      
      try {
         addMethods();
      }
      catch (Exception e) {
         e.printStackTrace();
      }
   }
   
   /**
    * returns the handler for read methods
    * 
    * @return the methodhandler for read methods
    */
   public MethodHandler read() {
      return m_ReadMethods;
   }
   
   /**
    * returns the handler for write methods
    * 
    * @return the methodhandler for read methods
    */
   public MethodHandler write() {
      return m_WriteMethods;
   }
   
   /**
    * adds read and write methods for the given class, i.e., read&;lt;name&gt;
    * and write&lt;name&gt; ("name" is prefixed by read and write)
    * 
    * @param handler  the handler class that contains the read and write method
    * @param cls      the class to register the read and write method for
    * @param name     the suffix of the read and write method
    */
   public void register(Object handler, Class cls, String name) {
     read().add(cls, XMLSerializationMethodHandler.findReadMethod(handler, "read" + name));
     write().add(cls, XMLSerializationMethodHandler.findWriteMethod(handler, "write" + name));
   }
   
   /**
    * returns the read and write method handlers as string
    * 
    * @return the read/write method handlers as string 
    */
   public String toString() {
      return "Read Methods:\n" + read() + "\n\n" + "Write Methods:\n" + write();
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
