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
 * XMLBasicSerialization.java
 * Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core.xml;

import weka.core.RevisionUtils;
import weka.core.Utils;

import java.io.StringReader;
import java.io.StringWriter;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map;
import java.util.Properties;
import java.util.Stack;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Vector;

import javax.swing.DefaultListModel;

import org.w3c.dom.Element;

/**
 * This serializer contains some read/write methods for common classes that
 * are not beans-conform. Currently supported are:
 * <ul>
 *    <li>java.util.HashMap</li>
 *    <li>java.util.HashSet</li>
 *    <li>java.util.Hashtable</li>
 *    <li>java.util.LinkedList</li>
 *    <li>java.util.Properties</li>
 *    <li>java.util.Stack</li>
 *    <li>java.util.TreeMap</li>
 *    <li>java.util.TreeSet</li>
 *    <li>java.util.Vector</li>
 *    <li>javax.swing.DefaultListModel</li>
 * </ul>
 *
 * Weka classes:
 * <ul>
 *    <li>weka.core.Matrix</li>
 *    <li>weka.core.matrix.Matrix</li>
 * </ul>
 * 
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5953 $ 
 */
public class XMLBasicSerialization
   extends XMLSerialization {

   /** the value for mapping, e.g., Maps */
   public final static String VAL_MAPPING = "mapping";

   /** the value for a mapping-key, e.g., Maps */
   public final static String VAL_KEY = "key";

   /** the value for mapping-value, e.g., Maps */
   public final static String VAL_VALUE = "value";

   /** the matrix cells */
   public final static String VAL_CELLS = "cells";

   /**
    * initializes the serialization
    * 
    * @throws Exception if initialization fails
    */
   public XMLBasicSerialization() throws Exception {
      super();
   }
   
   /**
    * generates internally a new XML document and clears also the IgnoreList
    * and the mappings for the Read/Write-Methods
    * 
    * @throws Exception if initializing fails
    */
   public void clear() throws Exception {
      super.clear();
      
      // Java classes
      m_CustomMethods.register(this, DefaultListModel.class, "DefaultListModel");
      m_CustomMethods.register(this, HashMap.class, "Map");
      m_CustomMethods.register(this, HashSet.class, "Collection");
      m_CustomMethods.register(this, Hashtable.class, "Map");
      m_CustomMethods.register(this, LinkedList.class, "Collection");
      m_CustomMethods.register(this, Properties.class, "Map");
      m_CustomMethods.register(this, Stack.class, "Collection");
      m_CustomMethods.register(this, TreeMap.class, "Map");
      m_CustomMethods.register(this, TreeSet.class, "Collection");
      m_CustomMethods.register(this, Vector.class, "Collection");

      // Weka classes
      m_CustomMethods.register(this, weka.core.matrix.Matrix.class, "Matrix");
      m_CustomMethods.register(this, weka.core.Matrix.class, "MatrixOld");
      m_CustomMethods.register(this, weka.classifiers.CostMatrix.class, "CostMatrixOld");
   }
   
   /**
    * adds the given DefaultListModel to a DOM structure. 
    * 
    * @param parent the parent of this object, e.g. the class this object is a
    * member of
    * @param o the Object to describe in XML
    * @param name the name of the object
    * @return the node that was created
    * @throws Exception if the DOM creation fails
    * @see javax.swing.DefaultListModel
    */
   public Element writeDefaultListModel(Element parent, Object o, String name) 
      throws Exception {

      Element              node;
      int                  i;
      DefaultListModel     model;

      // for debugging only
      if (DEBUG)
         trace(new Throwable(), name);
      
      m_CurrentNode = parent;
      
      model = (DefaultListModel) o;
      node = addElement(parent, name, o.getClass().getName(), false);

      for (i = 0; i < model.getSize(); i++)
         invokeWriteToXML(node, model.get(i), Integer.toString(i));
      
      return node;
   }

   /**
    * builds the DefaultListModel from the given DOM node. 
    * 
    * @param node the associated XML node
    * @return the instance created from the XML description
    * @throws Exception if instantiation fails 
    * @see javax.swing.DefaultListModel
    */
   public Object readDefaultListModel(Element node) throws Exception {
      DefaultListModel     model;
      Vector               children;
      Element              child;
      int                  i;
      int                  index;
      int                  currIndex;

      // for debugging only
      if (DEBUG)
         trace(new Throwable(), node.getAttribute(ATT_NAME));

      m_CurrentNode = node;
      
      children = XMLDocument.getChildTags(node); 
      model    = new DefaultListModel();
      
      // determine highest index for size
      index    = children.size() - 1;
      for (i = 0; i < children.size(); i++) {
        child     = (Element) children.get(i);
        currIndex = Integer.parseInt(child.getAttribute(ATT_NAME));
        if (currIndex > index)
          index = currIndex;
      }
      model.setSize(index + 1);

      // set values
      for (i = 0; i < children.size(); i++) {
         child = (Element) children.get(i);
         model.set(
             Integer.parseInt(child.getAttribute(ATT_NAME)), 
             invokeReadFromXML(child));
      }
      
      return model;
   }
   
   /**
    * adds the given Collection to a DOM structure. 
    * 
    * @param parent the parent of this object, e.g. the class this object is a
    * member of
    * @param o the Object to describe in XML
    * @param name the name of the object
    * @return the node that was created
    * @throws Exception if the DOM creation fails
    * @see java.util.Collection
    */
   public Element writeCollection(Element parent, Object o, String name) 
      throws Exception {

      Element         node;
      Iterator        iter;
      int             i;

      // for debugging only
      if (DEBUG)
         trace(new Throwable(), name);
      
      m_CurrentNode = parent;
      
      iter = ((Collection) o).iterator();
      node = addElement(parent, name, o.getClass().getName(), false);

      i = 0;
      while (iter.hasNext()) {
         invokeWriteToXML(node, iter.next(), Integer.toString(i));
         i++;
      }
      
      return node;
   }

   /**
    * builds the Collection from the given DOM node. 
    * 
    * @param node the associated XML node
    * @return the instance created from the XML description
    * @throws Exception if instantiation fails 
    * @see java.util.Collection
    */
   public Object readCollection(Element node) throws Exception {
      Collection<Object>           coll;
      Vector<Object>               v;
      Vector<Element>               children;
      Element              child;
      int                  i;
      int                  index;
      int                  currIndex;

      // for debugging only
      if (DEBUG)
         trace(new Throwable(), node.getAttribute(ATT_NAME));

      m_CurrentNode = node;
      
      children = XMLDocument.getChildTags(node); 
      v        = new Vector<Object>();

      // determine highest index for size
      index    = children.size() - 1;
      for (i = 0; i < children.size(); i++) {
        child     = (Element) children.get(i);
        currIndex = Integer.parseInt(child.getAttribute(ATT_NAME));
        if (currIndex > index)
          index = currIndex;
      }
      v.setSize(index + 1);


      // put the children in the vector to sort them according their index
      for (i = 0; i < children.size(); i++) {
         child = (Element) children.get(i);
         v.set(
               Integer.parseInt(child.getAttribute(ATT_NAME)), 
               invokeReadFromXML(child));
      }
      
      // populate collection
      coll = Utils.cast(Class.forName(node.getAttribute(ATT_CLASS)).
                        newInstance());
      coll.addAll(v);
      
      return coll;
   }
   
   /**
    * adds the given Map to a DOM structure. 
    * 
    * @param parent the parent of this object, e.g. the class this object is a
    * member of
    * @param o the Object to describe in XML
    * @param name the name of the object
    * @return the node that was created
    * @throws Exception if the DOM creation fails
    * @see java.util.Map
    */
   public Element writeMap(Element parent, Object o, String name) 
      throws Exception {

      Map            map;
      Object         key;
      Element        node;
      Element        child;
      Iterator       iter;

      // for debugging only
      if (DEBUG)
         trace(new Throwable(), name);
      
      m_CurrentNode = parent;
      
      map  = (Map) o;
      iter = map.keySet().iterator();
      node = addElement(parent, name, o.getClass().getName(), false);

      while (iter.hasNext()) {
         key   = iter.next();
         child = addElement(
                     node, VAL_MAPPING, Object.class.getName(), false);
         invokeWriteToXML(child, key,          VAL_KEY);
         invokeWriteToXML(child, map.get(key), VAL_VALUE);
      }
      
      return node;
   }

   /**
    * builds the Map from the given DOM node. 
    * 
    * @param node the associated XML node
    * @return the instance created from the XML description
    * @throws Exception if instantiation fails 
    * @see java.util.Map
    */
   public Object readMap(Element node) throws Exception {
     Map<Object,Object>                  map;
      Object               key;
      Object               value;
      Vector               children;
      Vector               cchildren;
      Element              child;
      Element              cchild;
      int                  i;
      int                  n;
      String               name;

      // for debugging only
      if (DEBUG)
         trace(new Throwable(), node.getAttribute(ATT_NAME));

      m_CurrentNode = node;
      
      map      = Utils.cast(Class.forName(node.getAttribute(ATT_CLASS)).
                            newInstance());
      children = XMLDocument.getChildTags(node); 

      for (i = 0; i < children.size(); i++) {
         child     = (Element) children.get(i);
         cchildren = XMLDocument.getChildTags(child);
         key       = null;
         value     = null;
         
         for (n = 0; n < cchildren.size(); n++) {
            cchild = (Element) cchildren.get(n);
            name   = cchild.getAttribute(ATT_NAME);
            if (name.equals(VAL_KEY))
               key = invokeReadFromXML(cchild);
            else if (name.equals(VAL_VALUE))
               value = invokeReadFromXML(cchild);
            else
               System.out.println("WARNING: '" 
                     + name + "' is not a recognized name for maps!");
         }
         
         map.put(key, value);
      }
      
      return map;
   }
   
   /**
    * adds the given Matrix to a DOM structure. 
    * 
    * @param parent the parent of this object, e.g. the class this object is a
    * member of
    * @param o the Object to describe in XML
    * @param name the name of the object
    * @return the node that was created
    * @throws Exception if the DOM creation fails
    * @see weka.core.matrix.Matrix
    */
   public Element writeMatrix(Element parent, Object o, String name) 
      throws Exception {

      weka.core.matrix.Matrix    matrix;
      Element                    node;

      // for debugging only
      if (DEBUG)
         trace(new Throwable(), name);
      
      m_CurrentNode = parent;
      
      matrix = (weka.core.matrix.Matrix) o;
      node   = addElement(parent, name, o.getClass().getName(), false);

      invokeWriteToXML(node, matrix.getArray(), VAL_CELLS);
      
      return node;
   }

   /**
    * builds the Matrix from the given DOM node. 
    * 
    * @param node the associated XML node
    * @return the instance created from the XML description
    * @throws Exception if instantiation fails 
    * @see weka.core.matrix.Matrix
    */
   public Object readMatrix(Element node) throws Exception {
      weka.core.matrix.Matrix    matrix;
      Vector                     children;
      Element                    child;
      int                        i;
      String                     name;
      Object                     o;

      // for debugging only
      if (DEBUG)
         trace(new Throwable(), node.getAttribute(ATT_NAME));

      m_CurrentNode = node;
      
      matrix   = null;
      children = XMLDocument.getChildTags(node); 
      for (i = 0; i < children.size(); i++) {
         child = (Element) children.get(i);
         name  = child.getAttribute(ATT_NAME);

         if (name.equals(VAL_CELLS)) {
            o = invokeReadFromXML(child);
            matrix = new weka.core.matrix.Matrix(
                        (double[][]) o);
         }
      }
      
      return matrix;
   }
   
   /**
    * adds the given Matrix (old) to a DOM structure. 
    * 
    * @param parent the parent of this object, e.g. the class this object is a
    * member of
    * @param o the Object to describe in XML
    * @param name the name of the object
    * @return the node that was created
    * @throws Exception if the DOM creation fails
    * @see weka.core.Matrix
    */
   public Element writeMatrixOld(Element parent, Object o, String name) 
      throws Exception {

      weka.core.Matrix  matrix;
      Element           node;
      double[][]        array;
      int               i;

      // for debugging only
      if (DEBUG)
         trace(new Throwable(), name);
      
      m_CurrentNode = parent;
      
      matrix = (weka.core.Matrix) o;
      node   = addElement(parent, name, o.getClass().getName(), false);

      array = new double[matrix.numRows()][];
      for (i = 0; i < array.length; i++)
         array[i] = matrix.getRow(i);
      invokeWriteToXML(node, array, VAL_CELLS);
      
      return node;
   }

   /**
    * builds the Matrix (old) from the given DOM node. 
    * 
    * @param node the associated XML node
    * @return the instance created from the XML description
    * @throws Exception if instantiation fails 
    * @see weka.core.Matrix
    */
   public Object readMatrixOld(Element node) throws Exception {
      weka.core.Matrix           matrix;
      weka.core.matrix.Matrix    matrixNew;

      // for debugging only
      if (DEBUG)
         trace(new Throwable(), node.getAttribute(ATT_NAME));

      m_CurrentNode = node;
      
      matrixNew = (weka.core.matrix.Matrix) readMatrix(node);
      matrix    = new weka.core.Matrix(matrixNew.getArrayCopy());
      
      return matrix;
   }
   
   /**
    * adds the given CostMatrix (old) to a DOM structure. 
    * 
    * @param parent the parent of this object, e.g. the class this object is a
    * member of
    * @param o the Object to describe in XML
    * @param name the name of the object
    * @return the node that was created
    * @throws Exception if the DOM creation fails
    * @see weka.classifiers.CostMatrix
    */
   public Element writeCostMatrixOld(Element parent, Object o, String name) 
      throws Exception {

      // for debugging only
      if (DEBUG)
         trace(new Throwable(), name);
      
      m_CurrentNode = parent;
      
      return writeMatrixOld(parent, o, name);
   }

   /**
    * builds the Matrix (old) from the given DOM node. 
    * 
    * @param node the associated XML node
    * @return the instance created from the XML description
    * @throws Exception if instantiation fails 
    * @see weka.classifiers.CostMatrix
    */
   public Object readCostMatrixOld(Element node) throws Exception {
      weka.classifiers.CostMatrix   matrix;
      weka.core.matrix.Matrix       matrixNew;
      StringWriter                  writer;

      // for debugging only
      if (DEBUG)
         trace(new Throwable(), node.getAttribute(ATT_NAME));

      m_CurrentNode = node;
      
      matrixNew = (weka.core.matrix.Matrix) readMatrix(node);
      writer    = new StringWriter();
      matrixNew.write(writer);
      matrix    = new weka.classifiers.CostMatrix(new StringReader(writer.toString()));
      
      return matrix;
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
