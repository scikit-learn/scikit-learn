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
 * PropertyHandler.java
 * Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 */

package weka.core.xml;

import weka.core.RevisionHandler;
import weka.core.RevisionUtils;

import java.util.Enumeration;
import java.util.HashSet;
import java.util.Hashtable;

/**
 * This class stores information about properties to ignore or properties
 * that are allowed for a certain class.
 * 
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5953 $ 
 */
public class PropertyHandler
   implements RevisionHandler {
  
   /** 
    * contains display names of properties to ignore in the serialization
    * process
    * 
    * @see #ignored()
    * @see #addIgnored(String)
    * @see #removeIgnored(String)
    * @see #isIgnored(String)
    */
  protected Hashtable<Object,HashSet<String>> m_Ignored = null;
   
   /**
    * lists for a class the properties allowed to use for setting and getting.
    * if a class is not listed, then all get/set-methods are allowed.<br>
    * Mapping: classname (String) - properties (HashSet, containing the Property-Names)
    * 
    * @see #allowed()
    * @see #addAllowed(Class,String)
    * @see #removeAllowed(Class,String)
    * @see #isAllowed(Class,String)
    */
  protected Hashtable<Object,HashSet<String>> m_Allowed = null;

   /**
    * initializes the handling 
    */
   public PropertyHandler() {
      super();

      m_Ignored = new Hashtable<Object,HashSet<String>>();
      m_Allowed = new Hashtable<Object,HashSet<String>>();
   }
   
   /**
    * returns an enumeration of the stored display names and classes of 
    * properties to ignore.<br> 
    * <b>NOTE:</b> String and Class Objects are mixed in this enumeration, depending
    * whether it is a global property to ignore or just one for a certain class!
    * 
    * @return the display names and classes
    * @see #m_Ignored
    */
   public Enumeration ignored() {
      return m_Ignored.keys();
   }
   
   /**
    * adds the given display name of a property to the ignore list. Can either
    * be a complete path (e.g. <code>__root__.options</code>) or only a 
    * property name (e.g. <code>options</code>). In the latter case it matches 
    * all occurences of this display name.
    * 
    * @param displayName the property to ignore
    * @see #m_Ignored 
    */
   public void addIgnored(String displayName) {
      HashSet<String>        list;
      
      list = new HashSet<String>();
      list.add(displayName);
      
      m_Ignored.put(displayName, list);
   }
   
   /**
    * adds the given class with the display name of a property to the ignore list.
    * I.e. this property is only ignored for this class. 
    * 
    * @param c the class for which a property is to be ignored
    * @param displayName the property to ignore
    * @see #m_Ignored 
    */
   public void addIgnored(Class c, String displayName) {
      HashSet<String>        list;
      
      // retrieve list
      if (m_Ignored.contains(c)) {
         list = (HashSet<String>) m_Ignored.get(c);
      }
      else {
         list = new HashSet<String>();
         m_Ignored.put(c, list);
      }
      
      list.add(displayName);
   }
   
   /**
    * removes the given display name from the ignore list. returns whether the 
    * removing was succesful, i.e. whether the display name was in the list.
    * 
    * @param displayName the property to remove from the ignore list
    * @return whether the ignore list contained the specified property
    * @see #m_Ignored
    */
   public boolean removeIgnored(String displayName) {
      return (m_Ignored.remove(displayName) != null);
   }
   
   /**
    * removes the given display name from the ignore list of the class. 
    * returns whether the removing was succesful, i.e. whether the display 
    * name was in the list.
    * 
    * @param c the class to remove the property from
    * @param displayName the property to remove from the ignore list
    * @return whether the ignore list contained the specified property
    * @see #m_Ignored
    */
   public boolean removeIgnored(Class c, String displayName) {
      HashSet        list;
      
      // retrieve list
      if (m_Ignored.contains(c))
         list = (HashSet) m_Ignored.get(c);
      else
         list = new HashSet();
      
      return list.remove(displayName);
   }
   
   /**
    * checks whether the given display name is an ignored property
    * 
    * @param displayName the property to check whether it is on the ignore list
    * @return whether the property is in the ignored list
    * @see #m_Ignored
    */
   public boolean isIgnored(String displayName) {
      return m_Ignored.containsKey(displayName);
   }
   
   /**
    * checks whether the given display name of a certain class is an ignored 
    * property. It only checks for this certain class and no derivative classes.
    * If you also want to check for derivative classes, use 
    * <code>isIgnored(Object,String)</code>. 
    * 
    * @param c the class to check the property for
    * @param displayName the property to check whether it is on the ignore list
    * @return whether the property is in the ignored list
    * @see #m_Ignored
    * @see #isIgnored(Object, String)
    */
   public boolean isIgnored(Class c, String displayName) {
      HashSet        list;
      
      // retrieve list
      if (m_Ignored.containsKey(c))
         list = (HashSet) m_Ignored.get(c);
      else
         list = new HashSet();

      return list.contains(displayName);
   }
   
   /**
    * checks whether the given display name of a given object is an ignored 
    * property. The object is checked for each stored class whether it is an 
    * <code>instanceof</code>. If the class is not stored then it will default
    * to <code>false</code>, since there are no restrictions for this class. 
    * 
    * @param o the object to check the property for
    * @param displayName the property to check whether it is on the ignore list
    * @return whether the property is in the ignored list
    * @see #m_Ignored
    */
   public boolean isIgnored(Object o, String displayName) {
      Enumeration    enm;
      Class          c;
      Object         element;
      boolean        result;
      HashSet        list;
      
      result = false;
      
      enm = ignored();
      while (enm.hasMoreElements()) {
         element = enm.nextElement();
         
         // has to be class! not a display name
         if (!(element instanceof Class))
            continue;
         
         c = (Class) element;
         
         // is it an instance of this class?
         if (c.isInstance(o)) {
            list   = (HashSet) m_Ignored.get(c);
            result = list.contains(displayName); 
            break;
         }
      }
      
      return result;
   }
   
   /**
    * returns an enumeration of the classnames for which only certain properties
    * (display names) are allowed
    * 
    * @return the classnames with restriction to properties
    */
   public Enumeration allowed() {
      return m_Allowed.keys();
   }
   
   /**
    * adds the given property (display name) to the list of allowed properties 
    * for the specified class.
    * 
    * @param c the class to add a property for
    * @param displayName the property to allow for the class
    * @see #m_Allowed
    */
   public void addAllowed(Class c, String displayName) {
      HashSet<String>        list;
      
      // retrieve list
      list = (HashSet<String>) m_Allowed.get(c);
      if (list == null) {
         list = new HashSet<String>();
         m_Allowed.put(c, list);
      }
      
      // add property
      list.add(displayName);
   }
   
   /**
    * removes the given property (display name) for the specified class from 
    * the list of allowed properties.
    * 
    * @param c the class to remove the property for
    * @param displayName the property to remove
    * @return whether the property was found
    * @see #m_Allowed
    */
   public boolean removeAllowed(Class c, String displayName) {
      boolean           result;
      HashSet           list;
      
      result = false;
      
      // retrieve list
      list = (HashSet) m_Allowed.get(c);
      
      // remove property
      if (list != null)
         result = list.remove(displayName);
      
      return result;
   }
   
   /**
    * returns whether the given property (display name) is allowed for the
    * given class. It only checks for this certain class and no derivative classes.
    * If you also want to check for derivative classes, use 
    * <code>isAllowed(Object,String)</code>.
    * 
    * @param c the class to check the property for
    * @param displayName the property (display name) to check
    * @return whether the property is allowed in that context
    * @see #m_Allowed
    * @see #isAllowed(Object, String) 
    */
   public boolean isAllowed(Class c, String displayName) {
      boolean        result;
      HashSet        list;
      
      result = true;
      
      // retrieve list
      list = (HashSet) m_Allowed.get(c);
      
      // check list
      if (list != null)
         result = list.contains(displayName);
      
      return result;
   }
   
   /**
    * returns whether the given property (display name) is allowed for the given
    * object . The object is checked for each stored class whether it is an 
    * <code>instanceof</code>. If the class is not stored then it will default
    * to <code>true</code>, since there are no restrictions for this class.
    * 
    * @param o the object to check the property for
    * @param displayName the property (display name) to check
    * @return whether the property is allowed in that context 
    */
   public boolean isAllowed(Object o, String displayName) {
      Enumeration    enm;
      Class          c;
      boolean        result;
      HashSet        list;
      
      result = true;
      
      enm = allowed();
      while (enm.hasMoreElements()) {
         c = (Class) enm.nextElement();
         
         // is it an instance of this class?
         if (c.isInstance(o)) {
            list   = (HashSet) m_Allowed.get(c);
            result = list.contains(displayName); 
            break;
         }
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
}
