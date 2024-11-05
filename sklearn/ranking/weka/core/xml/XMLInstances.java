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
 * XMLInstances.java
 * Copyright (C) 2006 University of Waikato, Hamilton, New Zealand
 */

package weka.core.xml;

import weka.core.Attribute;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.ProtectedProperties;
import weka.core.RevisionUtils;
import weka.core.DenseInstance;
import weka.core.SparseInstance;
import weka.core.Utils;
import weka.core.Version;
import weka.core.labelranking.PreferenceAttribute;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.io.Serializable;
import java.util.Enumeration;
import java.util.Properties;
import java.util.Vector;
import java.util.zip.GZIPInputStream;
import java.util.ArrayList;

import org.w3c.dom.Element;

/**
 * XML representation of the Instances class.
 * 
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5987 $
 */
public class XMLInstances
  extends XMLDocument
  implements Serializable {

  /** for serialization */
  private static final long serialVersionUID = 3626821327547416099L;
  
  /** The filename extension that should be used for xrff files */
  public static String FILE_EXTENSION = ".xrff";
  
  // tags
  /** the root element */
  public final static String TAG_DATASET = "dataset";

  /** the header element */
  public final static String TAG_HEADER = "header";

  /** the body element */
  public final static String TAG_BODY = "body";

  /** the notes element */
  public final static String TAG_NOTES = "notes";

  /** the attributes element */
  public final static String TAG_ATTRIBUTES = "attributes";

  /** the attribute element */
  public final static String TAG_ATTRIBUTE = "attribute";

  /** the labels element */
  public final static String TAG_LABELS = "labels";

  /** the label element */
  public final static String TAG_LABEL = "label";

  /** the meta-data element */
  public final static String TAG_METADATA = "metadata";

  /** the property element */
  public final static String TAG_PROPERTY = "property";

  /** the data element */
  public final static String TAG_INSTANCES = "instances";

  /** the instance element */
  public final static String TAG_INSTANCE = "instance";

  /** the value element */
  public final static String TAG_VALUE = "value";
  
  // attributes
  /** the version attribute */
  public final static String ATT_VERSION = "version";
  
  /** the type attribute */
  public final static String ATT_TYPE = "type";
  
  /** the format attribute (for date attributes) */
  public final static String ATT_FORMAT = "format";
  
  /** the class attribute */
  public final static String ATT_CLASS = "class";
  
  /** the index attribute */
  public final static String ATT_INDEX = "index";
  
  /** the weight attribute */
  public final static String ATT_WEIGHT = "weight";
  
  /** the missing attribute */
  public final static String ATT_MISSING = "missing";
  
  // values
  /** the value for numeric */
  public final static String VAL_NUMERIC = "numeric";
  
  /** the value for date */
  public final static String VAL_DATE = "date";
  
  /** the value for nominal */
  public final static String VAL_NOMINAL = "nominal";
  
  /** the value for string */
  public final static String VAL_STRING = "string";
  
  /** the value for relational */
  public final static String VAL_RELATIONAL = "relational";
  
  /** the value for normal */
  public final static String VAL_NORMAL = "normal";
  
  /** the value for sparse */
  public final static String VAL_SPARSE = "sparse";
  
  /** the DTD */
  public final static String DOCTYPE = 
      "<!" + DTD_DOCTYPE + " " + TAG_DATASET + "\n"
    + "[\n"
    + "   <!" + DTD_ELEMENT + " " + TAG_DATASET + " (" + TAG_HEADER + "," + TAG_BODY + ")" + ">\n"
    + "   <!" + DTD_ATTLIST + " " + TAG_DATASET + " " + ATT_NAME + " " + DTD_CDATA + " " + DTD_REQUIRED + ">\n"
    + "   <!" + DTD_ATTLIST + " " + TAG_DATASET + " " + ATT_VERSION + " " + DTD_CDATA + " \"" + Version.VERSION + "\">\n"
    + "\n"
    + "   <!" + DTD_ELEMENT + " " + TAG_HEADER + " (" + TAG_NOTES + DTD_OPTIONAL + "," + TAG_ATTRIBUTES + ")" + ">\n"
    + "   <!" + DTD_ELEMENT + " " + TAG_BODY   + " (" + TAG_INSTANCES  + ")" + ">\n"
    + "   <!" + DTD_ELEMENT + " " + TAG_NOTES + " " + DTD_ANY + ">   <!--  comments, information, copyright, etc. -->\n"
    + "\n"
    + "   <!" + DTD_ELEMENT + " " + TAG_ATTRIBUTES + " (" + TAG_ATTRIBUTE + DTD_AT_LEAST_ONE + ")" + ">\n"
    + "   <!" + DTD_ELEMENT + " " + TAG_ATTRIBUTE + " (" + TAG_LABELS + DTD_OPTIONAL + "," + TAG_METADATA + DTD_OPTIONAL + "," + TAG_ATTRIBUTES + DTD_OPTIONAL + ")" + ">\n"
    + "   <!" + DTD_ATTLIST + " " + TAG_ATTRIBUTE + " " + ATT_NAME + " " + DTD_CDATA + " " + DTD_REQUIRED + ">\n"
    + "   <!" + DTD_ATTLIST + " " + TAG_ATTRIBUTE + " " + ATT_TYPE + " (" + VAL_NUMERIC + DTD_SEPARATOR + VAL_DATE + DTD_SEPARATOR + VAL_NOMINAL + DTD_SEPARATOR + VAL_STRING + DTD_SEPARATOR + VAL_RELATIONAL + ") " + DTD_REQUIRED + ">\n"
    + "   <!" + DTD_ATTLIST + " " + TAG_ATTRIBUTE + " " + ATT_FORMAT + " " + DTD_CDATA + " " + DTD_IMPLIED + ">\n"
    + "   <!" + DTD_ATTLIST + " " + TAG_ATTRIBUTE + " " + ATT_CLASS + " (" + VAL_YES + DTD_SEPARATOR + VAL_NO + ") \"" + VAL_NO + "\"" + ">\n"
    + "   <!" + DTD_ELEMENT + " " + TAG_LABELS + " (" + TAG_LABEL + DTD_ZERO_OR_MORE + ")" + ">   <!-- only for type \"nominal\" -->\n"
    + "   <!" + DTD_ELEMENT + " " + TAG_LABEL + " " + DTD_ANY + ">\n"
    + "   <!" + DTD_ELEMENT + " " + TAG_METADATA + " (" + TAG_PROPERTY + DTD_ZERO_OR_MORE + ")" + ">\n"
    + "   <!" + DTD_ELEMENT + " " + TAG_PROPERTY + " " + DTD_ANY + ">\n"
    + "   <!" + DTD_ATTLIST + " " + TAG_PROPERTY + " " + ATT_NAME + " " + DTD_CDATA + " " + DTD_REQUIRED + ">\n"
    + "\n"
    + "   <!" + DTD_ELEMENT + " " + TAG_INSTANCES + " (" + TAG_INSTANCE + DTD_ZERO_OR_MORE + ")" + ">\n"
    + "   <!" + DTD_ELEMENT + " " + TAG_INSTANCE + " (" + TAG_VALUE + DTD_ZERO_OR_MORE + ")" + ">\n"
    + "   <!" + DTD_ATTLIST + " " + TAG_INSTANCE + " " + ATT_TYPE + " (" + VAL_NORMAL + DTD_SEPARATOR + VAL_SPARSE + ") \"" + VAL_NORMAL + "\"" + ">\n"
    + "   <!" + DTD_ATTLIST + " " + TAG_INSTANCE + " " + ATT_WEIGHT + " " + DTD_CDATA + " " + DTD_IMPLIED + ">\n"
    + "   <!" + DTD_ELEMENT + " " + TAG_VALUE + " (" + DTD_PCDATA + DTD_SEPARATOR + TAG_INSTANCES + ")" + DTD_ZERO_OR_MORE + ">\n"
    + "   <!" + DTD_ATTLIST + " " + TAG_VALUE + " " + ATT_INDEX + " " + DTD_CDATA + " " + DTD_IMPLIED + ">   <!-- 1-based index (only used for instance format \"sparse\") -->\n"
    + "   <!" + DTD_ATTLIST + " " + TAG_VALUE + " " + ATT_MISSING + " (" + VAL_YES + DTD_SEPARATOR + VAL_NO + ") \"" + VAL_NO + "\"" + ">\n"
    + "]\n"
    + ">";

  /** the precision for numbers */
  protected int m_Precision = 6;
  
  /** the underlying Instances */
  protected Instances m_Instances;
  
  /**
   * the default constructor
   * 
   * @throws Exception	if XML initialization fails
   */
  public XMLInstances() throws Exception {
    super();
    
    m_Instances = null;

    setDocType(DOCTYPE);
    setRootNode(TAG_DATASET);
    setValidating(true);
  }
  
  /**
   * generates the XML structure based on the given data
   * 
   * @param data	the data to build the XML structure from
   * @throws Exception	if initialization/generation fails
   */
  public XMLInstances(Instances data) throws Exception {
    this();
    
    setInstances(data);
  }
  
  /**
   * generates the Instances directly from the reader containing the
   * XML data.
   * 
   * @param reader	the reader for the XML data
   * @throws Exception	if something goes wrong
   */
  public XMLInstances(Reader reader) throws Exception {
    this();
    
    setXML(reader);
  }
  
  /**
   * adds the attribute to the XML structure
   * 
   * @param parent	the parent node to add the attribute node as child
   * @param att		the attribute to add
   */
  protected void addAttribute(Element parent, Attribute att) {
    Element		node;
    Element		child;
    Element		property;
    Element		label;
    String		tmpStr;
    Enumeration		enm;
    int			i;
    
    node = m_Document.createElement(TAG_ATTRIBUTE);
    parent.appendChild(node);
    
    // XML attributes
    // name
    node.setAttribute(ATT_NAME, validContent(att.name()));
    
    // type
    switch (att.type()) {
      case Attribute.NUMERIC:
	node.setAttribute(ATT_TYPE, VAL_NUMERIC);
	break;
	
      case Attribute.DATE:
	node.setAttribute(ATT_TYPE, VAL_DATE);
	break;
	
      case Attribute.NOMINAL:
      case PreferenceAttribute.RANKING:
	node.setAttribute(ATT_TYPE, VAL_NOMINAL);
	break;
	
      case Attribute.STRING:
	node.setAttribute(ATT_TYPE, VAL_STRING);
	break;
	
      case Attribute.RELATIONAL:
	node.setAttribute(ATT_TYPE, VAL_RELATIONAL);
	break;
	
      default:
	node.setAttribute(ATT_TYPE, "???");
    }
    
    // labels
    if (att.isNominal() || att.isRanking()) {
      child = m_Document.createElement(TAG_LABELS);
      node.appendChild(child);
      enm = att.enumerateValues();
      while (enm.hasMoreElements()) {
	tmpStr = enm.nextElement().toString();
	label = m_Document.createElement(TAG_LABEL);
	child.appendChild(label);
	label.appendChild(m_Document.createTextNode(validContent(tmpStr)));
      }
    }
    
    // format
    if (att.isDate())
      node.setAttribute(ATT_FORMAT, validContent(att.getDateFormat()));
    
    // class
    if (m_Instances.classIndex() > -1) {
      if (att == m_Instances.classAttribute())
	node.setAttribute(ATT_CLASS, VAL_YES);
    }
    
    // add meta-data
    if ( (att.getMetadata() != null) && (att.getMetadata().size() > 0) ) {
      child = m_Document.createElement(TAG_METADATA);
      node.appendChild(child);
      enm = att.getMetadata().propertyNames();
      while (enm.hasMoreElements()) {
	tmpStr = enm.nextElement().toString();
	property = m_Document.createElement(TAG_PROPERTY);
	child.appendChild(property);
	property.setAttribute(ATT_NAME, validContent(tmpStr));
	property.appendChild(m_Document.createTextNode(validContent(att.getMetadata().getProperty(tmpStr, ""))));
      }
    }
    
    // relational attribute?
    if (att.isRelationValued()) {
      child = m_Document.createElement(TAG_ATTRIBUTES);
      node.appendChild(child);
      for (i = 0; i < att.relation().numAttributes(); i++)
	addAttribute(child, att.relation().attribute(i));
    }
  }

  /**
   * turns all &lt;, &gt; and &amp;into character entities and returns that 
   * string. Necessary for TextNodes.
   * 
   * @param content	string to convert
   * @return		the valid content string
   */
  protected String validContent(String content) {
    String	result;
    
    result = content;
    
    // these five entities are recognized by every XML processor
    // see http://www.xml.com/pub/a/2001/03/14/trxml10.html
    result = result.replaceAll("&", "&amp;")
                   .replaceAll("\"", "&quot;")
                   .replaceAll("'", "&apos;")
                   .replaceAll("<", "&lt;")
                   .replaceAll(">", "&gt;");
    // in addition, replace some other entities as well
    result = result.replaceAll("\n", "&#10;")
                   .replaceAll("\r", "&#13;")
                   .replaceAll("\t", "&#9;");
    
    return result;
  }
  
  /**
   * adds the instance to the XML structure
   * 
   * @param parent	the parent node to add the instance node as child
   * @param inst	the instance to add
   */
  protected void addInstance(Element parent, Instance inst) {
    Element		node;
    Element		value;
    Element		child;
    boolean		sparse;
    int			i;
    int			n;
    int			index;
    
    node = m_Document.createElement(TAG_INSTANCE);
    parent.appendChild(node);
    
    // sparse?
    sparse = (inst instanceof SparseInstance);
    if (sparse)
      node.setAttribute(ATT_TYPE, VAL_SPARSE);
    
    // weight
    if (inst.weight() != 1.0)
      node.setAttribute(ATT_WEIGHT, Utils.doubleToString(inst.weight(), m_Precision));
    
    // values
    for (i = 0; i < inst.numValues(); i++) {
      index = inst.index(i);
      
      value = m_Document.createElement(TAG_VALUE);
      node.appendChild(value);

      if (inst.isMissing(index)) {
	value.setAttribute(ATT_MISSING, VAL_YES);
      }
      else {
	if (inst.attribute(index).isRelationValued()) {
	  child = m_Document.createElement(TAG_INSTANCES);
	  value.appendChild(child);
	  for (n = 0; n < inst.relationalValue(i).numInstances(); n++)
	    addInstance(child, inst.relationalValue(i).instance(n));
	}
	else {
	  if (inst.attribute(index).type() == Attribute.NUMERIC)
	    value.appendChild(m_Document.createTextNode(Utils.doubleToString(inst.value(index), m_Precision)));
	  else
	    value.appendChild(m_Document.createTextNode(validContent(inst.stringValue(index))));
	}
      }
      
      if (sparse)
	value.setAttribute(ATT_INDEX, "" + (index+1));
    }
  }
  
  /**
   * generates the XML structure for the header
   */
  protected void headerToXML() {
    Element	root;
    Element	node;
    Element	child;
    int		i;
    
    root = m_Document.getDocumentElement();
    root.setAttribute(ATT_NAME, validContent(m_Instances.relationName()));
    root.setAttribute(ATT_VERSION, Version.VERSION);
    
    // create "header" node
    node = m_Document.createElement(TAG_HEADER);
    root.appendChild(node);

    // add all attributes
    child = m_Document.createElement(TAG_ATTRIBUTES);
    node.appendChild(child);
    for (i = 0; i < m_Instances.numAttributes(); i++)
      addAttribute(child, m_Instances.attribute(i));
  }
  
  /**
   * generates the XML structure from the rows
   */
  protected void dataToXML() {
    Element	root;
    Element	node;
    Element	child;
    int		i;
    
    root = m_Document.getDocumentElement();
    
    // create "body" node
    node = m_Document.createElement(TAG_BODY);
    root.appendChild(node);

    // add all instances
    child = m_Document.createElement(TAG_INSTANCES);
    node.appendChild(child);
    for (i = 0; i < m_Instances.numInstances(); i++)
      addInstance(child, m_Instances.instance(i));
  }
  
  /**
   * builds up the XML structure based on the given data
   * 
   * @param data	data to generate the XML from
   */
  public void setInstances(Instances data) {
    m_Instances = new Instances(data);
    clear();
    headerToXML();
    dataToXML();
  }
  
  /**
   * returns the current instances, either the ones that were set or the ones
   * that were generated from the XML structure.
   * 
   * @return		the current instances
   */
  public Instances getInstances() {
    return m_Instances;
  }

  /**
   * returns the metadata, if any available underneath this node, otherwise
   * just null
   * 
   * @param parent	the attribute node
   * @return		the metadata, or null if none found
   * @throws Exception	if generation fails
   */
  protected ProtectedProperties createMetadata(Element parent) throws Exception {
    ProtectedProperties	result;
    Properties		props;
    Vector		list;
    Element		node;
    Element		metanode;
    int			i;
    
    result = null;
    
    // find metadata node directly underneath this attribute, but not in
    // deeper nested attributes (e.g., within relational attributes)
    metanode = null;
    list     = getChildTags(parent, TAG_METADATA);
    if (list.size() > 0)
      metanode = (Element) list.get(0);
    
    // generate properties
    if (metanode != null) {
      props = new Properties();
      list  = getChildTags(metanode, TAG_PROPERTY);
      for (i = 0; i < list.size(); i++) {
	node = (Element) list.get(i);
	props.setProperty(node.getAttribute(ATT_NAME), getContent(node));
      }
      result = new ProtectedProperties(props);
    }
    
    return result;
  }

  /**
   * returns the labels listed underneath this (nominal) attribute in a 
   * ArrayList
   * 
   * @param parent	the (nominal) attribute node
   * @return		the label vector
   * @throws Exception	if generation fails
   */
  protected ArrayList<String> createLabels(Element parent) throws Exception {
    ArrayList<String>		result;
    Vector		list;
    Element		node;
    Element		labelsnode;
    int			i;
    
    result = new ArrayList<String>();
    
    // find labels node directly underneath this attribute, but not in
    // deeper nested attributes (e.g., within relational attributes)
    labelsnode = null;
    list     = getChildTags(parent, TAG_LABELS);
    if (list.size() > 0)
      labelsnode = (Element) list.get(0);
    
    // retrieve all labels
    if (labelsnode != null) {
      list  = getChildTags(labelsnode, TAG_LABEL);
      for (i = 0; i < list.size(); i++) {
	node = (Element) list.get(i);
	result.add(getContent(node));
      }
    }
    
    return result;
  }
  
  /**
   * creates an Attribute from the given XML node
   * 
   * @param node	the node with the setup
   * @return		the configured Attribute
   * @throws Exception	if generation fails, e.g., due to unknown attribute type
   */
  protected Attribute createAttribute(Element node) throws Exception {
    String		typeStr;
    String		name;
    int			type;
    Attribute		result;
    ArrayList<String>		values;
    ProtectedProperties	metadata;
    Vector		list;
    ArrayList<Attribute>		atts;
    
    result = null;
    
    // name
    name = node.getAttribute(ATT_NAME);

    // type
    typeStr = node.getAttribute(ATT_TYPE);
    if (typeStr.equals(VAL_NUMERIC))
      type = Attribute.NUMERIC;
    else if (typeStr.equals(VAL_DATE))
      type = Attribute.DATE;
    else if (typeStr.equals(VAL_NOMINAL))
      type = Attribute.NOMINAL;
    else if (typeStr.equals(VAL_STRING))
      type = Attribute.STRING;
    else if (typeStr.equals(VAL_RELATIONAL))
      type = Attribute.RELATIONAL;
    else
      throw new Exception(
	  "Attribute type '" + typeStr + "' is not supported!");

    // metadata
    metadata = createMetadata(node);
    
    switch (type) {
      case Attribute.NUMERIC:
	if (metadata == null)
	  result = new Attribute(name);
	else
	  result = new Attribute(name, metadata);
	break;

      case Attribute.DATE:
	if (metadata == null)
	  result = new Attribute(name, node.getAttribute(ATT_FORMAT));
	else
	  result = new Attribute(name, node.getAttribute(ATT_FORMAT), metadata);
	break;
	
      case Attribute.NOMINAL:
      case PreferenceAttribute.RANKING:
	values = createLabels(node);
	if (metadata == null)
	  result = new Attribute(name, values);
	else
	  result = new Attribute(name, values, metadata);
	break;
	
      case Attribute.STRING:
	if (metadata == null)
	  result = new Attribute(name, (ArrayList<String>) null);
	else
	  result = new Attribute(name, (ArrayList<String>) null, metadata);
	break;
	
      case Attribute.RELATIONAL:
	list = getChildTags(node, TAG_ATTRIBUTES);
	node = (Element) list.get(0);
	atts = createAttributes(node, new int[1]);
	if (metadata == null)
	  result = new Attribute(name, new Instances(name, atts, 0));
	else
	  result = new Attribute(name, new Instances(name, atts, 0), metadata);
	break;
    }
    
    return result;
  }
  
  /**
   * returns a list of generated attributes
   * 
   * @param parent	the attributes node
   * @param classIndex	array of length 1 to return the class index, if any
   * @return		the vector with the generated attributes
   * @throws Exception	if generation fails, e.g., due to unknown attribute type
   */
  protected ArrayList<Attribute> createAttributes(Element parent, int[] classIndex) throws Exception {
    Vector	list;
    ArrayList<Attribute>	result;
    int		i;
    Element	node;
    Attribute	att;

    result        = new ArrayList<Attribute>();
    classIndex[0] = -1;
    
    list = getChildTags(parent, TAG_ATTRIBUTE);
    for (i = 0; i < list.size(); i++) {
      node = (Element) list.get(i);
      att = createAttribute(node);
      if (node.getAttribute(ATT_CLASS).equals(VAL_YES))
	classIndex[0] = i;
      result.add(att);
    }
    
    return result;
  }
  
  /**
   * creates an Instance from the given XML node
   * 
   * @param header	the data this instance will belong to
   * @param parent	the instance node
   * @return		the configured Instance
   * @throws Exception	if generation fails, e.g., due to unknown attribute type
   */
  protected Instance createInstance(Instances header, Element parent) throws Exception {
    Instance	result;
    Element	node;
    Element	child;
    boolean	sparse;
    int		i;
    int		index;
    Vector	list;
    Vector	subList;
    double[]	values;
    String	content;
    double	weight;
    Instances	data;
    
    result = null;

    // sparse?
    sparse = (parent.getAttribute(ATT_TYPE).equals(VAL_SPARSE));
    values = new double[header.numAttributes()];
    
    // weight
    if (parent.getAttribute(ATT_WEIGHT).length() != 0)
      weight = Double.parseDouble(parent.getAttribute(ATT_WEIGHT));
    else
      weight = 1.0;
    
    list = getChildTags(parent, TAG_VALUE);
    for (i = 0; i < list.size(); i++) {
      node = (Element) list.get(i);
      
      // determine index
      if (sparse)
	index = Integer.parseInt(node.getAttribute(ATT_INDEX)) - 1;
      else
	index = i;

      // set value
      if (node.getAttribute(ATT_MISSING).equals(VAL_YES)) {
	values[index] = Utils.missingValue();
      }
      else {
	content = getContent(node);
	switch (header.attribute(index).type()) {
	  case Attribute.NUMERIC:
	    values[index] = Double.parseDouble(content);
	    break;
	    
	  case Attribute.DATE:
	    values[index] = header.attribute(index).parseDate(content);
	    break;
	    
	  case Attribute.NOMINAL:
	  case PreferenceAttribute.RANKING:
	    values[index] = header.attribute(index).indexOfValue(content);
	    break;
	    
	  case Attribute.STRING:
	    values[index] = header.attribute(index).addStringValue(content);
	    break;
	    
	  case Attribute.RELATIONAL:
	    subList       = getChildTags(node, TAG_INSTANCES);
	    child         = (Element) subList.get(0);
	    data          = createInstances(header.attribute(index).relation(), child);
	    values[index] = header.attribute(index).addRelation(data);
	    break;
	    
	  default:
	    throw new Exception(
		"Attribute type " + header.attribute(index).type() 
		+ " is not supported!");  
	}
      }
    }
    
    // create instance
    if (sparse)
      result = new SparseInstance(weight, values);
    else
      result = new DenseInstance(weight, values);
    
    return result;
  }
  
  /**
   * creates Instances from the given XML node
   * 
   * @param header	the header of this data
   * @param parent	the instances node
   * @return		the generated Instances
   * @throws Exception	if generation fails, e.g., due to unknown attribute type
   */
  protected Instances createInstances(Instances header, Element parent) throws Exception {
    Instances	result;
    Vector	list;
    int		i;
    
    result = new Instances(header, 0);
    
    list = getChildTags(parent, TAG_INSTANCE);
    for (i = 0; i < list.size(); i++)
      result.add(createInstance(result, (Element) list.get(i)));
    
    return result;
  }
  
  /**
   * generates the header from the XML document
   * 
   * @return		the generated header
   * @throws Exception	if generation fails
   */
  protected Instances headerFromXML() throws Exception {
    Instances	result;
    Element	root;
    Element	node;
    Vector	list;
    ArrayList<Attribute>	atts;
    Version	version;
    int[]	classIndex;

    root = m_Document.getDocumentElement();
    
    // check version
    version = new Version();
    if (version.isOlder(root.getAttribute(ATT_VERSION)))
      System.out.println(
	  "WARNING: loading data of version " + root.getAttribute(ATT_VERSION)
	  + " with version " + Version.VERSION);
    
    // attributes
    list       = getChildTags(root, TAG_HEADER);
    node       = (Element) list.get(0);
    list       = getChildTags(node, TAG_ATTRIBUTES);
    node       = (Element) list.get(0);
    classIndex = new int[1];
    atts       = createAttributes(node, classIndex);

    // generate header
    result = new Instances(root.getAttribute(ATT_NAME), atts, 0);
    result.setClassIndex(classIndex[0]);
    
    return result;
  }
  
  /**
   * generates the complete dataset from the XML document
   * 
   * @param header	the header structure
   * @return		the complete dataset
   * @throws Exception	if generation fails
   */
  protected Instances dataFromXML(Instances header) throws Exception {
    Instances	result;
    Element	node;
    Vector	list;

    list   = getChildTags(m_Document.getDocumentElement(), TAG_BODY);
    node   = (Element) list.get(0);
    list   = getChildTags(node, TAG_INSTANCES);
    node   = (Element) list.get(0);
    result = createInstances(header, node);
    
    return result;
  }
  
  /**
   * reads the XML structure from the given reader
   * 
   * @param reader	the reader to get the XML from
   * @throws Exception	if 
   */
  public void setXML(Reader reader) throws Exception {
    read(reader);
    
    // interprete XML structure
    m_Instances = dataFromXML(headerFromXML());
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5987 $");
  }
  
  /**
   * takes an XML document as first argument and then outputs the Instances
   * statistics
   * 
   * @param args	the commandline options
   */
  public static void main(String[] args) {
    try {
      Reader r = null;
      if (args.length != 1) {
	throw (new Exception("Usage: XMLInstances <filename>"));
      }
      else {
	InputStream in = new FileInputStream(args[0]);
	// compressed file?
	if (args[0].endsWith(".gz"))
	  in = new GZIPInputStream(in);
        r = new BufferedReader(new InputStreamReader(in));
     }
      
      if (args[0].endsWith(Instances.FILE_EXTENSION)) {
	XMLInstances i = new XMLInstances(new Instances(r));
	System.out.println(i.toString());
      }
      else {
	Instances i = new XMLInstances(r).getInstances();
	System.out.println(i.toSummaryString());
      }
    }
    catch (Exception ex) {
      ex.printStackTrace();
      System.err.println(ex.getMessage());
    }
  }
}
