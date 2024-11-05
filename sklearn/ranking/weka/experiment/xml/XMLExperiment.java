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
 * XMLExperiment.java
 * Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 */

package weka.experiment.xml;

import java.beans.PropertyDescriptor;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.Vector;

import org.w3c.dom.Element;

import weka.core.RevisionUtils;
import weka.core.xml.XMLBasicSerialization;
import weka.core.xml.XMLDocument;
import weka.experiment.Experiment;
import weka.experiment.PropertyNode;

/**
 * This class serializes and deserializes an Experiment instance to and
 * fro XML.<br>
 * It omits the <code>options</code> from the Experiment, since these are handled
 * by the get/set-methods. For the <code>Classifier</code> class with all its 
 * derivative classes it stores only <code>debug</code> and <code>options</code>.
 * For <code>SplitEvaluator</code> and <code>ResultProducer</code> only the
 * options are retrieved. The <code>PropertyNode</code> is done manually since
 * it has no get/set-methods for its public fields.<br>
 * Since there's no read-method for <code>m_ClassFirst</code> we always save it
 * as <code>false</code>.
 * 
 * @see Experiment#m_ClassFirst
 * 
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 1.6 $ 
 */
public class XMLExperiment
   extends XMLBasicSerialization {
  
   /** the name of the classFirst property */
   public final static String NAME_CLASSFIRST = "classFirst";

   /** PropertyNode member */
   public final static String NAME_PROPERTYNODE_VALUE = "value";

   /** PropertyNode member */
   public final static String NAME_PROPERTYNODE_PARENTCLASS = "parentClass";

   /** PropertyNode member */
   public final static String NAME_PROPERTYNODE_PROPERTY = "property";
   
   /**
    * initializes the serialization
    * 
    * @throws Exception if initialization fails
    */
   public XMLExperiment() throws Exception {
      super();
   }
   
   /**
    * generates internally a new XML document and clears also the IgnoreList and
    * the mappings for the Read/Write-Methods
    * 
    * @throws Exception if initializing fails
    */
   public void clear() throws Exception {
      super.clear();

      // ignore
      m_Properties.addIgnored(VAL_ROOT + ".options");
      m_Properties.addIgnored(Experiment.class, "options");
      
      // allow
      m_Properties.addAllowed(weka.classifiers.Classifier.class, "debug");
      m_Properties.addAllowed(weka.classifiers.Classifier.class, "options");
      // we assume that classes implementing SplitEvaluator also implement OptionHandler
      m_Properties.addAllowed(weka.experiment.SplitEvaluator.class, "options");
      // we assume that classes implementing ResultProducer also implement OptionHandler
      m_Properties.addAllowed(weka.experiment.ResultProducer.class, "options");
      
      // read/write methods
      m_CustomMethods.register(this, PropertyNode.class, "PropertyNode");
   }
   
   /**
    * enables derived classes to add other properties to the DOM tree, e.g.
    * ones that do not apply to the get/set convention of beans. 
    * 
    * @param o the object that is serialized into XML
    * @throws Exception if post-processing fails
    */
   protected void writePostProcess(Object o) throws Exception {
      Element              node;
      Experiment           exp;

      exp = (Experiment) o;
      
      // classfirst
      node = addElement(m_Document.getDocument().getDocumentElement(), NAME_CLASSFIRST, Boolean.class.getName(), false);
      node.appendChild(node.getOwnerDocument().createTextNode(new Boolean(false).toString()));   // TODO: get-Method for classFirst in Experiment???
   }
   
   /**
    * additional post-processing can happen in derived classes after reading 
    * from XML. 
    * 
    * @param o the object to perform some additional processing on
    * @return the processed object
    * @throws Exception if post-processing fails
    */
   protected Object readPostProcess(Object o) throws Exception {
      Element              node;
      Experiment           exp;
      int                  i;
      Vector               children;

      exp = (Experiment) o;
      
      // classfirst
      children = XMLDocument.getChildTags(m_Document.getDocument().getDocumentElement());
      for (i = 0; i < children.size(); i++) {
         node = (Element) children.get(i);
         if (node.getAttribute(ATT_NAME).equals(NAME_CLASSFIRST)) {
            exp.classFirst(new Boolean(XMLDocument.getContent(node)).booleanValue());
            break;
         }
      }
      
      return o;
   }
   
   /**
    * adds the given PropertyNode to a DOM structure. 
    * 
    * @param parent the parent of this object, e.g. the class this object is a member of
    * @param o the Object to describe in XML
    * @param name the name of the object
    * @return the node that was created
    * @throws Exception if the DOM creation fails
    */
   public Element writePropertyNode(Element parent, Object o, String name) throws Exception {
      Element              node;
      PropertyNode         pnode;
      Vector               children;
      int                  i;
      Element              child;

      // for debugging only
      if (DEBUG)
         trace(new Throwable(), name);
      
      m_CurrentNode = parent;
      
      pnode = (PropertyNode) o;
      node  = (Element) parent.appendChild(m_Document.getDocument().createElement(TAG_OBJECT));
      node.setAttribute(ATT_NAME, name);
      node.setAttribute(ATT_CLASS, pnode.getClass().getName());
      node.setAttribute(ATT_PRIMITIVE, VAL_NO);
      node.setAttribute(ATT_ARRAY, VAL_NO);
      
      if (pnode.value != null)
         invokeWriteToXML(node, pnode.value, NAME_PROPERTYNODE_VALUE);
      if (pnode.parentClass != null)
         invokeWriteToXML(node, pnode.parentClass.getName(), NAME_PROPERTYNODE_PARENTCLASS);
      if (pnode.property != null)
         invokeWriteToXML(node, pnode.property.getDisplayName(), NAME_PROPERTYNODE_PROPERTY);
      
      // fix primitive values
      if (    (pnode.value != null) 
           && (pnode.property != null) 
           && (pnode.property.getPropertyType().isPrimitive())) {
         children = XMLDocument.getChildTags(node);
         for (i = 0; i < children.size(); i++) {
            child = (Element) children.get(i);
            if (!child.getAttribute(ATT_NAME).equals(NAME_PROPERTYNODE_VALUE))
               continue;
            child.setAttribute(ATT_CLASS, pnode.property.getPropertyType().getName());
            child.setAttribute(ATT_PRIMITIVE, VAL_YES);
         }
      }
      
      return node;
   }

   /**
    * builds the PropertyNode from the given DOM node. 
    * 
    * @param node the associated XML node
    * @return the instance created from the XML description
    * @throws Exception if instantiation fails 
    * @see javax.swing.DefaultListModel
    */
   public Object readPropertyNode(Element node) throws Exception {
      Object               result;
      Object               value;
      String               parentClass;
      String               property;
      Vector               children;
      Element              child;
      int                  i;
      Class                cls;

      // for debugging only
      if (DEBUG)
         trace(new Throwable(), node.getAttribute(ATT_NAME));

      m_CurrentNode = node;
      
      result      = null;

      children    = XMLDocument.getChildTags(node);
      value       = null;
      parentClass = null;
      property    = null;
      
      for (i = 0; i < children.size(); i++) {
         child = (Element) children.get(i);
         
         if (child.getAttribute(ATT_NAME).equals(NAME_PROPERTYNODE_VALUE)) {
            if (stringToBoolean(child.getAttribute(ATT_PRIMITIVE)))
               value = getPrimitive(child);
            else
               value = invokeReadFromXML(child);
         }
         if (child.getAttribute(ATT_NAME).equals(NAME_PROPERTYNODE_PARENTCLASS))
            parentClass = XMLDocument.getContent(child);
         if (child.getAttribute(ATT_NAME).equals(NAME_PROPERTYNODE_PROPERTY))
            property = XMLDocument.getContent(child);
      }
      
      if (parentClass != null)
         cls = Class.forName(parentClass);
      else
         cls = null;
      
      if (cls != null)
         result = new PropertyNode(value, new PropertyDescriptor(property, cls), cls);
      else
         result = new PropertyNode(value);
      
      return result;
   }
   
   /**
    * Returns the revision string.
    * 
    * @return		the revision
    */
   public String getRevision() {
     return RevisionUtils.extract("$Revision: 1.6 $");
   }

   /**
    * for testing only. if the first argument is a filename with ".xml"
    * as extension it tries to generate an instance from the XML description
    * and does a <code>toString()</code> of the generated object.
    * Otherwise it loads the binary file, saves the XML representation in a 
    * file with the original filename appended by ".xml" and once again in a
    * binary file with the original filename appended by ".exp".
    * 
    * @param args 	the commandline arguments
    * @throws Exception	if something goes wrong, e.g., file not found
    */
   public static void main(String[] args) throws Exception {
      if (args.length > 0) {
         // read xml and print
         if (args[0].toLowerCase().endsWith(".xml")) {
            System.out.println(new XMLExperiment().read(args[0]).toString());
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
            //new XMLExperiment().write(System.out, o);
            // write to XML file
            new XMLExperiment().write(new BufferedOutputStream(new FileOutputStream(args[0] + ".xml")), o);
            // print to binary file
            FileOutputStream fo = new FileOutputStream(args[0] + ".exp");
            ObjectOutputStream oo = new ObjectOutputStream(
                                   new BufferedOutputStream(fo));
            oo.writeObject(o);
            oo.close();
         }
      }
   }
}
