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
 * PropertyPath.java
 * Copyright (C) 2006 University of Waikato, Hamilton, New Zealand
 */

package weka.core;

import java.beans.PropertyDescriptor;
import java.lang.reflect.Array;
import java.lang.reflect.Method;
import java.util.StringTokenizer;
import java.util.Vector;

/**
 * A helper class for accessing properties in nested objects, e.g., accessing
 * the "getRidge" method of a LinearRegression classifier part of 
 * MultipleClassifierCombiner, e.g., Vote. For doing so, one needs to 
 * supply the object to work on and a property path. The property path is a
 * dot delimited path of property names ("getFoo()" and "setFoo(int)" have 
 * "foo" as property name), indices of arrays are 0-based. E.g.: <p/>
 * 
 * <code>getPropertyDescriptor(vote, "classifiers[1].ridge")</code> will return
 * the second classifier (which should be our LinearRegression) of the given
 * Vote meta-classifier and there the property descriptor of the "ridge" 
 * property. <code>getValue(...)</code> will return the actual value of the
 * ridge parameter and <code>setValue(...)</code> will set it.
 * 
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5953 $
 */
public class PropertyPath
  implements RevisionHandler {

  /**
   * Represents a single element of a property path
   * 
   * @author  fracpete (fracpete at waikato dot ac dot nz)
   * @version $Revision: 5953 $
   */
  public static class PathElement
    implements Cloneable, RevisionHandler {
    
    /** the property */
    protected String m_Name;
    
    /** the index of the array (-1 for none) */
    protected int m_Index;
    
    /**
     * initializes the path element with the given property
     * 
     * @param property	the property to initialize with
     */
    public PathElement(String property) {
      super();
      
      if (property.indexOf("[") > -1) {
	m_Name  = property.replaceAll("\\[.*$", "");
	m_Index = Integer.parseInt(
    		     property.replaceAll(".*\\[", "").replaceAll("\\].*", ""));
      }
      else {
	m_Name   = property;
	m_Index = -1;
      }
    }

    /**
     * returns a clone of the current object
     * 
     * @return		the clone of the current state
     */
    public Object clone() {
      return new PathElement(this.toString());
    }
    
    /**
     * returns the name of the property
     * 
     * @return		the name of the property
     */
    public String getName() {
      return m_Name;
    }
    
    /**
     * returns whether the property is an index-based one
     * 
     * @return		true if the property has an index
     */
    public boolean hasIndex() {
      return (getIndex() > -1);
    }
    
    /**
     * returns the index of the property, -1 if the property is not an
     * index-based one
     * 
     * @return		the index of the property
     */
    public int getIndex() {
      return m_Index;
    }
    
    /**
     * returns the element once again as string
     * 
     * @return		the property as string
     */
    public String toString() {
      String	result;
      
      result = getName();
      if (hasIndex())
	result += "[" + getIndex() + "]";
      
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
  
  /**
   * Contains a (property) path structure
   * 
   * @author  fracpete (fracpete at waikato dot ac dot nz)
   * @version $Revision: 5953 $
   */
  public static class Path
    implements RevisionHandler {
    
    /** the structure */
    protected Vector<PathElement> m_Elements;
    
    /**
     * default constructor, only used internally
     */
    protected Path() {
      super();
      
      m_Elements = new Vector<PathElement>();
    }
    
    /**
     * uses the given dot-path
     * 
     * @param path	path in dot-notation
     */
    public Path(String path) {
      this();
      
      m_Elements = breakUp(path);
    }
    
    /**
     * uses the vector with PathElement objects to initialize with
     * 
     * @param elements	the PathElements to use
     */
    public Path(Vector<PathElement> elements) {
      this();
      
      for (int i = 0; i < elements.size(); i++)
	m_Elements.add((PathElement) elements.get(i).clone());
    }
    
    /**
     * uses the given array as elements for the path
     * 
     * @param elements	the path elements to use
     */
    public Path(String[] elements) {
      this();
      
      for (int i = 0; i < elements.length; i++)
	m_Elements.add(new PathElement(elements[i]));
    }
    
    /**
     * breaks up the given path and returns it as vector
     * 
     * @param path	the path to break up
     * @return		the single elements of the path
     */
    protected Vector<PathElement> breakUp(String path) {
      Vector<PathElement>		result;
      StringTokenizer	tok;
      
      result = new Vector<PathElement>();
      
      tok = new StringTokenizer(path, ".");
      while (tok.hasMoreTokens())
        result.add(new PathElement(tok.nextToken()));
      
      return result;
    }

    /**
     * returns the element at the given index
     * 
     * @param index	the index of the element to return
     * @return		the specified element
     */
    public PathElement get(int index) {
      return (PathElement) m_Elements.get(index);
    }

    /**
     * returns the number of path elements of this structure
     * 
     * @return		the number of path elements
     */
    public int size() {
      return m_Elements.size();
    }
    
    /**
     * returns a path object based on the given path string
     * 
     * @param path	path to work on
     * @return		the path structure
     */
    public static Path parsePath(String path) {
      return new Path(path);
    }

    /**
     * returns a subpath of the current structure, starting with the specified
     * element index up to the end
     * 
     * @param startIndex	the first element of the subpath
     * @return			the new subpath
     */
    public Path subpath(int startIndex) {
      return subpath(startIndex, size());
    }

    /**
     * returns a subpath of the current structure, starting with the specified
     * element index up. The endIndex specifies the element that is not part
     * of the new subpath. In other words, the new path contains the elements
     * from "startIndex" up to "(endIndex-1)".
     * 
     * @param startIndex	the first element of the subpath
     * @param endIndex		the element that is after the last added element
     * @return			the new subpath
     */
    public Path subpath(int startIndex, int endIndex) {
      Vector<PathElement>	list;
      int	i;
      
      list = new Vector<PathElement>();
      for (i = startIndex; i < endIndex; i++)
	list.add(get(i));
      
      return new Path(list);
    }
    
    /**
     * returns the structure again as a dot-path
     * 
     * @return		the path structure as dot-path
     */
    public String toString() {
      String	result;
      int	i;
      
      result = "";
      
      for (i = 0; i < m_Elements.size(); i++) {
	if (i > 0)
	  result += ".";
	result += m_Elements.get(i);
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

  /**
   * A helper class that stores Object and PropertyDescriptor together.
   * 
   * @author  fracpete (fracpete at waikato dot ac dot nz)
   * @version $Revision: 5953 $
   */
  protected static class PropertyContainer
    implements RevisionHandler {
    
    /** the descriptor */
    protected PropertyDescriptor m_Descriptor;
    
    /** the associated object */
    protected Object m_Object;
    
    /**
     * initializes the container
     * 
     * @param desc	the property descriptor
     * @param obj	the associated object
     */
    public PropertyContainer(PropertyDescriptor desc, Object obj) {
      super();
      
      m_Descriptor = desc;
      m_Object     = obj;
    }
    
    /**
     * returns the stored descriptor
     * 
     * @return		the stored descriptor
     */
    public PropertyDescriptor getDescriptor() {
      return m_Descriptor;
    }
    
    /**
     * returns the stored object
     * 
     * @return		the stored object
     */
    public Object getObject() {
      return m_Object;
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
  
  /**
   * returns the property and object associated with the given path, null if 
   * a problem occurred.
   * 
   * @param src		the object to start from
   * @param path	the path to follow
   * @return		not null, if the property could be found
   */
  public static PropertyContainer find(Object src, Path path) {
    PropertyContainer	result;
    PropertyDescriptor	desc;
    Object		newSrc;
    PathElement		part;
    Method		method;
    Object		methodResult;

    // get descriptor
    part = path.get(0);
    try {
      desc = new PropertyDescriptor(part.getName(), src.getClass());
    }
    catch (Exception e) {
      desc = null;
      e.printStackTrace();
    }

    // problem occurred? -> stop
    if (desc == null)
      return null;
    
    // end of path reached?
    if (path.size() == 1) {
      result = new PropertyContainer(desc, src);
    }
    // recurse further
    else {
      try {
	method       = desc.getReadMethod();
	methodResult = method.invoke(src, (Object[]) null);
	if (part.hasIndex())
	  newSrc = Array.get(methodResult, part.getIndex());
	else
	  newSrc = methodResult;
	result = find(newSrc, path.subpath(1));
      }
      catch (Exception e) {
	result = null;
	e.printStackTrace();
      }
    }
    
    return result;
  }
  
  /**
   * returns the property associated with the given path, null if a problem
   * occurred.
   * 
   * @param src		the object to start from
   * @param path	the path to follow
   * @return		not null, if the property could be found
   */
  public static PropertyDescriptor getPropertyDescriptor(Object src, Path path) {
    PropertyContainer	cont;
    
    cont = find(src, path);
    
    if (cont == null)
      return null;
    else
      return cont.getDescriptor();
  }
  
  /**
   * returns the property associated with the given path
   * 
   * @param src		the object to start from
   * @param path	the path to follow
   * @return		not null, if the property could be found
   */
  public static PropertyDescriptor getPropertyDescriptor(Object src, String path) {
    return getPropertyDescriptor(src, new Path(path));
  }
  
  /**
   * returns the value specified by the given path from the object
   * 
   * @param src		the object to work on
   * @param path	the retrieval path
   * @return		the value, null if an error occurred
   */
  public static Object getValue(Object src, Path path) {
    Object		result;
    PropertyContainer	cont;
    Method		method;
    Object		methodResult;
    PathElement		part;
    
    result = null;
    
    cont = find(src, path);
    // problem?
    if (cont == null)
      return null;
    
    // retrieve the value
    try {
      part         = path.get(path.size() - 1);
      method       = cont.getDescriptor().getReadMethod();
      methodResult = method.invoke(cont.getObject(), (Object[]) null);
      if (part.hasIndex())
	result = Array.get(methodResult, part.getIndex());
      else
	result = methodResult;
    }
    catch (Exception e) {
      result = null;
      e.printStackTrace();
    }
    
    return result;
  }
  
  /**
   * returns the value specified by the given path from the object
   * 
   * @param src		the object to work on
   * @param path	the retrieval path
   * @return		the value, null if an error occurred
   */
  public static Object getValue(Object src, String path) {
    return getValue(src, new Path(path));
  }
  
  /**
   * set the given value specified by the given path in the object
   * 
   * @param src		the object to work on
   * @param path	the retrieval path
   * @param value	the value to set
   * @return		true if the value could be set
   */
  public static boolean setValue(Object src, Path path, Object value) {
    boolean		result;
    PropertyContainer	cont;
    Method		methodRead;
    Method		methodWrite;
    Object		methodResult;
    PathElement		part;
    
    result = false;
    
    cont = find(src, path);
    // problem?
    if (cont == null)
      return result;
    
    // set the value
    try {
      part         = path.get(path.size() - 1);
      methodRead   = cont.getDescriptor().getReadMethod();
      methodWrite  = cont.getDescriptor().getWriteMethod();
      if (part.hasIndex()) {
	methodResult = methodRead.invoke(cont.getObject(), (Object[]) null);
	Array.set(methodResult, part.getIndex(), value);
	methodWrite.invoke(cont.getObject(), new Object[]{methodResult});
      }
      else {
	methodWrite.invoke(cont.getObject(), new Object[]{value});
      }
      result = true;
    }
    catch (Exception e) {
      result = false;
      e.printStackTrace();
    }
    
    return result;
  }
  
  /**
   * set the given value specified by the given path in the object
   * 
   * @param src		the object to work on
   * @param path	the retrieval path
   * @param value	the value to set
   */
  public static void setValue(Object src, String path, Object value) {
    setValue(src, new Path(path), value);
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
   * for testing only
   * 
   * @param args	the commandline options - ignored
   * @throws Exception	if something goes wrong
   */
  public static void main(String[] args) throws Exception {
    // Path
    Path path = new Path("hello.world[2].nothing");
    System.out.println("Path: " + path);
    System.out.println(" -size: " + path.size());
    System.out.println(" -elements:");
    for (int i = 0; i < path.size(); i++)
      System.out.println(
	  "  " + i + ". " + path.get(i).getName() 
	  + " -> " + path.get(i).getIndex());
    
    /*
    // retrieving ridge with path
    weka.classifiers.meta.Vote vote = new weka.classifiers.meta.Vote();
    vote.setClassifiers(
	new weka.classifiers.Classifier[]{
	    new weka.classifiers.trees.J48(),
	    new weka.classifiers.functions.LinearRegression()});
    path = new Path("classifiers[1].ridge");
    System.out.println("path: " + path + " -> " + getValue(vote, path));
    
    // setting ridge with path and retrieving it again
    setValue(vote, path.toString(), new Double(0.1));
    System.out.println("path: " + path + " -> " + getValue(vote, path));
    */
  }
}
