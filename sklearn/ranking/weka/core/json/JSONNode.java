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
 * JSONObject.java
 * Copyright (C) 2009 University of Waikato, Hamilton, New Zealand
 */

package weka.core.json;

import java.awt.BorderLayout;
import java.io.Reader;

import java_cup.runtime.DefaultSymbolFactory;
import java_cup.runtime.SymbolFactory;

import javax.swing.JFrame;
import javax.swing.JScrollPane;
import javax.swing.JTree;
import javax.swing.tree.DefaultMutableTreeNode;

/**
 * Container class for storing a <a href="http://www.json.org/" target="_blank">JSON</a> 
 * data structure.
 * 
 * @author  FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5786 $
 */
public class JSONNode
  extends DefaultMutableTreeNode {

  /** for serialization. */
  private static final long serialVersionUID = -3047440914507883491L;

  /**
   * The type of a node.
   * 
   * @author  FracPete (fracpete at waikato dot ac dot nz)
   * @version $Revision: 5786 $
   */
  public static enum NodeType {
    /** a primitive. */
    PRIMITIVE,
    /** an object with nested key-value pairs. */
    OBJECT,
    /** an array. */
    ARRAY
  }
 
  /** the name of the node. */
  protected String m_Name;
  
  /** the value of the node. */
  protected Object m_Value;
  
  /** the type of the node. */
  protected NodeType m_NodeType;
  
  /**
   * Initializes the root container.
   */
  public JSONNode() {
    this(null, NodeType.OBJECT);
  }
  
  /**
   * Initializes the primitive container.
   * 
   * @param name	the name
   * @param value	the primitive value
   */
  public JSONNode(String name, Boolean value) {
    this(name, value, NodeType.PRIMITIVE);
  }
  
  /**
   * Initializes the primitive container.
   * 
   * @param name	the name
   * @param value	the primitive value
   */
  public JSONNode(String name, Integer value) {
    this(name, value, NodeType.PRIMITIVE);
  }
  
  /**
   * Initializes the primitive container.
   * 
   * @param name	the name
   * @param value	the primitive value
   */
  public JSONNode(String name, Double value) {
    this(name, value, NodeType.PRIMITIVE);
  }
  
  /**
   * Initializes the primitive container.
   * 
   * @param name	the name
   * @param value	the primitive value
   */
  public JSONNode(String name, String value) {
    this(name, value, NodeType.PRIMITIVE);
  }
  
  /**
   * Initializes the object container with null value.
   * 
   * @param name	the name
   * @param type	the node type
   */
  protected JSONNode(String name, NodeType type) {
    this(name, null, type);
  }
  
  /**
   * Initializes the container.
   * 
   * @param name	the name
   * @param value	the primitive value
   * @param type	the type of the node, null for primitives
   */
  protected JSONNode(String name, Object value, NodeType type) {
    super();
    
    m_Name     = name;
    m_Value    = value;
    m_NodeType = type;
  }
  
  /**
   * Checks whether the node is anonymous.
   * 
   * @return		true if no name available
   */
  public boolean isAnonymous() {
    return (m_Name == null);
  }
  
  /**
   * Returns the name of the node.
   * 
   * @return		the name, null for anonymous nodes
   */
  public String getName() {
    return m_Name;
  }
  
  /**
   * Returns the stored value.
   * 
   * @return		the stored value, can be null
   */
  public Object getValue() {
    return getValue(null);
  }
  
  /**
   * Returns the stored value.
   * 
   * @param defValue	the default value, if value is null
   * @return		the stored value, can be null
   */
  public Object getValue(Object defValue) {
    if (m_Value == null)
      return defValue;
    else
      return m_Value;
  }
  
  /**
   * Returns whether the node stores a primitive value or a an array/object.
   * 
   * @return		true if a primitive, false in case of an array/object
   */
  public boolean isPrimitive() {
    return (m_NodeType == NodeType.PRIMITIVE);
  }
  
  /**
   * Returns wether the node is an array.
   * 
   * @return		true if the node is array container
   */
  public boolean isArray() {
    return (m_NodeType == NodeType.ARRAY);
  }
  
  /**
   * Returns wether the node is an object.
   * 
   * @return		true if the node is object container
   */
  public boolean isObject() {
    return (m_NodeType == NodeType.OBJECT);
  }
  
  /**
   * Returns the type of the container.
   * 
   * @return		the type
   */
  public NodeType getNodeType() {
    return m_NodeType;
  }
  
  /**
   * Adds a "null" child to the object.
   * 
   * @param name	the name of the null value
   * @return		the new node, or null if none added
   */
  public JSONNode addNull(String name) {
    return add(name, null, NodeType.PRIMITIVE);
  }
  
  /**
   * Adds a key-value child to the object.
   * 
   * @param name	the name of the pair
   * @param value	the value
   * @return		the new node, or null if none added
   */
  public JSONNode addPrimitive(String name, Boolean value) {
    return add(name, value, NodeType.PRIMITIVE);
  }
  
  /**
   * Adds a key-value child to the object.
   * 
   * @param name	the name of the pair
   * @param value	the value
   * @return		the new node, or null if none added
   */
  public JSONNode addPrimitive(String name, Integer value) {
    return add(name, value, NodeType.PRIMITIVE);
  }
  
  /**
   * Adds a key-value child to the object.
   * 
   * @param name	the name of the pair
   * @param value	the value
   * @return		the new node, or null if none added
   */
  public JSONNode addPrimitive(String name, Double value) {
    return add(name, value, NodeType.PRIMITIVE);
  }
  
  /**
   * Adds a key-value child to the object.
   * 
   * @param name	the name of the pair
   * @param value	the value
   * @return		the new node, or null if none added
   */
  public JSONNode addPrimitive(String name, String value) {
    return add(name, value, NodeType.PRIMITIVE);
  }
  
  /**
   * Adds an array child to the object.
   * 
   * @param name	the name of the pair
   * @return		the new node, or null if none added
   */
  public JSONNode addArray(String name) {
    return add(name, null, NodeType.ARRAY);
  }
  
  /**
   * Adds an array element child to the array.
   * 
   * @param value	the value of the element array
   * @return		the new node, or null if none added
   */
  public JSONNode addArrayElement(Object value) {
    NodeType	type;

    if (getNodeType() != NodeType.ARRAY)
      return null;
    
    type = null;
    
    if (value != null) {
      if (value instanceof Boolean)
	type = NodeType.PRIMITIVE;
      else if (value instanceof Integer)
	type = NodeType.PRIMITIVE;
      else if (value instanceof Double)
	type = NodeType.PRIMITIVE;
      else if (value instanceof String)
	type = NodeType.PRIMITIVE;
      else if (value.getClass().isArray())
	type = NodeType.ARRAY;
      else
	type = NodeType.OBJECT;
    }
      
    return add(null, value, type);
  }
  
  /**
   * Adds an object child to the object.
   * 
   * @param name	the name of the pair
   * @return		the new node, or null if none added
   */
  public JSONNode addObject(String name) {
    return add(name, null, NodeType.OBJECT);
  }
  
  /**
   * Adds a key-value child to the object.
   * 
   * @param name	the name of the pair
   * @param value	the value
   * @param type	the node type, null for primitives
   * @return		the new node, or null if none added
   */
  protected JSONNode add(String name, Object value, NodeType type) {
    JSONNode	child;
    
    if (isPrimitive())
      return null;
    
    child = new JSONNode(name, value, type);
    add(child);
    
    return child;
  }
  
  /**
   * Checks whether the node has a child with the given name.
   * 
   * @param name	the name of the child
   * @return		true if child with that name is available
   */
  public boolean hasChild(String name) {
    return (getChild(name) != null);
  }
  
  /**
   * Returns the child with the given name.
   * 
   * @param name	the name of the child
   * @return		the child if available, null otherwise
   */
  public JSONNode getChild(String name) {
    JSONNode	result;
    JSONNode	node;
    int		i;
    
    result = null;
    
    for (i = 0; i < getChildCount(); i++) {
      node = (JSONNode) getChildAt(i);
      if (!node.isAnonymous() && node.getName().equals(name)) {
	result = node;
	break;
      }
    }
    
    return result;
  }
  
  /**
   * Generates the indentation string.
   * 
   * @param level	the level
   * @return		the indentation string (tabs)
   */
  protected String getIndentation(int level) {
    StringBuffer	result;
    int			i;
    
    result = new StringBuffer();
    for (i = 0; i < level; i++)
      result.append("\t");
    
    return result.toString();
  }
  
  /**
   * Escapes ", \, /, \b, \f, \n, \r, \t in strings.
   * 
   * @param o		the object to process (only strings get processed)
   * @return		the processed object
   */
  protected Object escape(Object o) {
    if (o instanceof String)
      return escape((String) o);
    else
      return o;
  }
  
  /**
   * Escapes ", /, \b, \f, \n, \r, \t.
   * 
   * @param s		the string to process
   * @return		the processed
   */
  protected String escape(String s) {
    StringBuffer	result;
    int			i;
    char		c;
    
    if (    (s.indexOf('\"') > -1)
         || (s.indexOf('\\') > -1) 
         || (s.indexOf('\b') > -1) 
         || (s.indexOf('\f') > -1) 
         || (s.indexOf('\n') > -1) 
         || (s.indexOf('\r') > -1) 
         || (s.indexOf('\t') > -1) ) {
      result = new StringBuffer();
      for (i = 0; i < s.length(); i++) {
	c = s.charAt(i);
	if (c == '\"')
	  result.append("\\\"");
	else if (c == '\\')
	  result.append("\\\\");
	else if (c == '\b')
	  result.append("\\b");
	else if (c == '\f')
	  result.append("\\f");
	else if (c == '\n')
	  result.append("\\n");
	else if (c == '\r')
	  result.append("\\r");
	else if (c == '\t')
	  result.append("\\t");
	else
	  result.append(c);
      }
    }
    else {
      result = new StringBuffer(s);
    }
    
    return result.toString();
  }
  
  /**
   * Dumps the node structure into JSON format.
   * 
   * @param buffer	the buffer to add the data to
   */
  public void toString(StringBuffer buffer) {
    int		level;
    boolean	isLast;
    String	indent;
    int		i;
    
    level  = getLevel();
    isLast = (getNextSibling() == null);
    indent = getIndentation(level);
    
    buffer.append(indent);
    if (m_Name != null) {
      buffer.append("\"");
      buffer.append(escape(m_Name));
      buffer.append("\" : ");
    }
    
    if (isObject()) {
      buffer.append("{\n");
      for (i = 0; i < getChildCount(); i++)
	((JSONNode) getChildAt(i)).toString(buffer);
      buffer.append(indent);
      buffer.append("}");
    }
    else if (isArray()) {
      buffer.append("[\n");
      for (i = 0; i < getChildCount(); i++)
	((JSONNode) getChildAt(i)).toString(buffer);
      buffer.append(indent);
      buffer.append("]");
    }
    else {
      if (m_Value == null) {
	buffer.append("null");
      }
      else if (m_Value instanceof String) {
	buffer.append("\"");
	buffer.append(escape((String) m_Value));
	buffer.append("\"");
      }
      else {
	buffer.append(m_Value.toString());
      }
    }
    
    if (!isLast)
      buffer.append(",");
    buffer.append("\n");
  }

  /**
   * Returns a string representation of the node.
   * 
   * @return		the string representation
   */
  public String toString() {
    String	result;
    
    result = null;
    
    if (isObject()) {
      if (isRoot())
	result = "JSON";
      else if (m_Name == null)
	result = "<object>";
      else
	result = escape(m_Name) + " (Object)";
    }
    else if (isArray()) {
      if (m_Name == null)
	result = "<array>";
      else
	result = escape(m_Name) + " (Array)";
    }
    else {
      if (m_Name != null)
	result = escape(m_Name) + ": " + escape(m_Value);
      else
	result = "" + m_Value;
    }
    
    return result;
  }
  
  /**
   * Reads the JSON object from the given reader.
   * 
   * @param reader	the reader to read the JSON object from
   * @return		the generated JSON object
   * @throws Exception	if parsing fails
   */
  public static JSONNode read(Reader reader) throws Exception {
    SymbolFactory 	sf;
    Parser 		parser;
    
    sf     = new DefaultSymbolFactory();
    parser = new Parser(new Scanner(reader, sf), sf);
    parser.parse();
    
    return parser.getResult();
  }
  
  /**
   * Only for testing. Generates a simple JSON object and displays it.
   * 
   * @param args	ignored
   * @throws Exception	if something goes wrong
   */
  public static void main(String[] args) throws Exception {
    // generates the example listed here:
    // http://en.wikipedia.org/wiki/JSON
    JSONNode person = new JSONNode();
    person.addPrimitive("firstName", "John");
    person.addPrimitive("lastName", "Smith");
    JSONNode address = person.addObject("address");
    address.addPrimitive("streetAddress", "21 2nd Street");
    address.addPrimitive("city", "New York");
    address.addPrimitive("state", "NY");
    address.addPrimitive("postalCode", 10021);
    JSONNode phonenumbers = person.addArray("phoneNumbers");
    phonenumbers.addArrayElement("212 555-1234");
    phonenumbers.addArrayElement("646 555-4567");
    
    // output in console
    StringBuffer buffer = new StringBuffer();
    person.toString(buffer);
    System.out.println(buffer.toString());
    
    // display GUI
    JTree tree = new JTree(person);
    JFrame frame = new JFrame("JSON");
    frame.setSize(800, 600);
    frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    frame.getContentPane().setLayout(new BorderLayout());
    frame.getContentPane().add(new JScrollPane(tree), BorderLayout.CENTER);
    frame.setLocationRelativeTo(null);
    frame.setVisible(true);
  }
}
