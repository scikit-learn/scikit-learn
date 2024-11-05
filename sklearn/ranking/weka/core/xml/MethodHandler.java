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
 * MethodHandler.java
 * Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 */

package weka.core.xml;

import weka.core.RevisionHandler;
import weka.core.RevisionUtils;

import java.lang.reflect.Method;
import java.util.Enumeration;
import java.util.Hashtable;

/**
 * This class handles relationships between display names of properties 
 * (or classes) and Methods that are associated with them. 
 * 
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5953 $ 
 */
public class MethodHandler
   implements RevisionHandler {
  
   /** 
    * stores the properties/class - Method relationship
    * 
    * @see #keys()
    * @see #add(Class, Method)
    * @see #add(String, Method)
    * @see #remove(Class)
    * @see #remove(String)
    * @see #get(Class)
    * @see #get(String)
    * @see #contains(Class)
    * @see #contains(String)  
    */
  protected Hashtable<Object,Method> m_Methods = null;
   
   /**
    * initializes the handler
    */
   public MethodHandler() {
      super();
      m_Methods  = new Hashtable<Object,Method>();
   }
   
   /**
    * returns an enumeration over all currently stored custom methods, i.e. it
    * returns the display names/classes in the enumeration.
    * 
    * @return the currently stored methods
    * @see #m_Methods
    */
   public Enumeration keys() {
      return m_Methods.keys();
   }
   
   /**
    * adds the specified method for the property with the given displayname
    * to its internal list.
    * 
    * @param displayName the display name of the property to handle manually
    * @param method the method, which will be invoked by reflection to handle
    *        the property manually 
    * @see #m_Methods
    */
   public void add(String displayName, Method method) {
      if (method != null)
         m_Methods.put(displayName, method);
   }
   
   /**
    * adds the specified method for the given class to its internal list.
    * 
    * @param c the class to handle manually
    * @param method the method, which will be invoked by reflection to handle
    *        the property manually 
    * @see #m_Methods
    */
   public void add(Class c, Method method) {
      if (method != null)
         m_Methods.put(c, method);
   }
   
   /**
    * removes the method for the property specified by the display name
    * from its internal list. 
    * 
    * @param displayName the display name of the propery to remove the custom
    *        method for
    * @return whether the method was stored in the list at all 
    * @see #m_Methods
    */
   public boolean remove(String displayName) {
      return (m_Methods.remove(displayName) != null);
   }
   
   /**
    * removes the method for the specified class from its internal list. 
    * 
    * @param c the class to remove the custom method for
    * @return whether the method was stored in the list at all 
    * @see #m_Methods
    */
   public boolean remove(Class c) {
      return (m_Methods.remove(c) != null);
   }
   
   /**
    * checks whether a method is stored for the given property
    * 
    * @param displayName the display name of the property to check for a method
    * @return whether a method is currently stored
    * @see #m_Methods
    */
   public boolean contains(String displayName) {
      return m_Methods.containsKey(displayName);
   }
   
   /**
    * checks whether a method is stored for the given class
    * 
    * @param c the class to check for a method
    * @return whether a method is currently stored
    * @see #m_Methods
    */
   public boolean contains(Class c) {
      return m_Methods.containsKey(c);
   }
   
   /**
    * returns the stored method for the given property
    * 
    * @param displayName the display name of the property to retrieve the 
    *        method for
    * @return the method associated with the display name, can be <code>null</code>
    * @see #m_Methods
    */
   public Method get(String displayName) {
      return (Method) m_Methods.get(displayName);
   }
   
   /**
    * returns the stored method for the given class
    * 
    * @param c the class to retrieve the method for
    * @return the method associated with the class, can be <code>null</code>
    * @see #m_Methods
    */
   public Method get(Class c) {
      return (Method) m_Methods.get(c);
   }
   
   /**
    * returns the number of currently stored Methods
    * 
    * @return the nummber of methods
    */
   public int size() {
      return m_Methods.size();
   }
   
   /**
    * removes all mappings 
    */
   public void clear() {
      m_Methods.clear();
   }
   
   /**
    * returns the internal Hashtable (propety/class - method relationship) in
    * a string representation
    * 
    * @return the object as string 
    */
   public String toString() {
      return m_Methods.toString();
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
