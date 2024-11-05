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
 * XMLOptions.java
 * Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core.xml;

import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.core.Utils;

import java.io.File;
import java.io.InputStream;
import java.io.Reader;
import java.util.Vector;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

/**
 * A class for transforming options listed in XML to a regular WEKA command
 * line string.<p>
 *
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5953 $
 */
public class XMLOptions
  implements RevisionHandler {
  
  /** tag for a single option. */
  public final static String TAG_OPTION = "option";

  /** tag for a list of options. */
  public final static String TAG_OPTIONS = "options";

  /** the name attribute. */
  public final static String ATT_NAME = "name";

  /** the type attribute. */
  public final static String ATT_TYPE = "type";

  /** the value attribute. */
  public final static String ATT_VALUE = "value";

  /** a value of the type attribute. */
  public final static String VAL_TYPE_FLAG = "flag";

  /** a value of the type attribute. */
  public final static String VAL_TYPE_SINGLE = "single";

  /** a value of the type attribute. */
  public final static String VAL_TYPE_HYPHENS = "hyphens";

  /** a value of the type attribute. */
  public final static String VAL_TYPE_QUOTES = "quotes";

  /** a value of the type attribute. */
  public final static String VAL_TYPE_CLASSIFIER = "classifier";

  /** a value of the type attribute. */
  public final static String VAL_TYPE_OPTIONHANDLER = "optionhandler";

  /** the root node. */
  public final static String ROOT_NODE = TAG_OPTIONS;

  /** the DTD for the XML file. */
  public final static String DOCTYPE = 
    "<!DOCTYPE " + ROOT_NODE + "\n"
    + "[\n"
    + "   <!ELEMENT " + TAG_OPTIONS + " (" + TAG_OPTION + ")*>\n"
    + "   <!ATTLIST " + TAG_OPTIONS + " " + ATT_TYPE + " CDATA \"" + VAL_TYPE_OPTIONHANDLER + "\">\n"
    + "   <!ATTLIST " + TAG_OPTIONS + " " + ATT_VALUE + " CDATA \"\">\n"
    + "   <!ELEMENT " + TAG_OPTION + " (#PCDATA | " + TAG_OPTIONS + ")*>\n"
    + "   <!ATTLIST " + TAG_OPTION + " " + ATT_NAME + " CDATA #REQUIRED>\n"
    + "   <!ATTLIST " + TAG_OPTION + " " + ATT_TYPE + " (" + VAL_TYPE_FLAG + " | " + VAL_TYPE_SINGLE + " | " + VAL_TYPE_HYPHENS + " | " + VAL_TYPE_QUOTES + ") \"" + VAL_TYPE_SINGLE + "\">\n"
    + "]\n"
    + ">";

  /** the XML document. */
  protected XMLDocument m_XMLDocument = null;

  /** 
   * Creates a new instance of XMLOptions.
   * 
   * @throws Exception 	if the construction of the DocumentBuilder fails
   * @see 		#setValidating(boolean)
   */
  public XMLOptions() throws Exception {
    m_XMLDocument = new XMLDocument(); 
    m_XMLDocument.setRootNode(ROOT_NODE);
    m_XMLDocument.setDocType(DOCTYPE);
    setValidating(true);
  }

  /** 
   * Creates a new instance of XMLOptions.
   *  
   * @param xml 	the xml to parse (if "<?xml" is not found then it is considered a file)
   * @throws Exception 	if the construction of the DocumentBuilder fails
   * @see 		#setValidating(boolean)
   */
  public XMLOptions(String xml) throws Exception {
    this();
    getXMLDocument().read(xml);
  }

  /** 
   * Creates a new instance of XMLOptions.
   * 
   * @param file 	the XML file to parse
   * @throws Exception 	if the construction of the DocumentBuilder fails
   * @see 		#setValidating(boolean)
   */
  public XMLOptions(File file) throws Exception {
    this();
    getXMLDocument().read(file);
  }

  /** 
   * Creates a new instance of XMLOptions.
   * 
   * @param stream 	the XML stream to parse
   * @throws Exception 	if the construction of the DocumentBuilder fails
   * @see 		#setValidating(boolean)
   */
  public XMLOptions(InputStream stream) throws Exception {
    this();
    getXMLDocument().read(stream);
  }

  /** 
   * Creates a new instance of XMLOptions.
   * 
   * @param reader 	the XML reader to parse
   * @throws Exception 	if the construction of the DocumentBuilder fails
   * @see 		#setValidating(boolean)
   */
  public XMLOptions(Reader reader) throws Exception {
    this();
    getXMLDocument().read(reader);
  }

  /**
   * returns whether a validating parser is used.
   * 
   * @return 		whether a validating parser is used
   */
  public boolean getValidating() {
    return m_XMLDocument.getValidating();
  }

  /**
   * sets whether to use a validating parser or not. <br>
   * Note: this does clear the current DOM document! 
   * 
   * @param validating 	whether to use a validating parser
   * @throws Exception 	if the instantiating of the DocumentBuilder fails
   */
  public void setValidating(boolean validating) throws Exception {
    m_XMLDocument.setValidating(validating);
  }

  /**
   * returns the parsed DOM document.
   * 
   * @return 		the parsed DOM document
   */
  public Document getDocument() {
    fixHyphens();
    return m_XMLDocument.getDocument();
  }

  /**
   * returns the handler of the XML document. the internal DOM document can 
   * be accessed via the <code>getDocument()</code> method.
   * 
   * @return 		the object handling the XML document
   * @see 		#getDocument()
   */
  public XMLDocument getXMLDocument() {
    fixHyphens();
    return m_XMLDocument;
  }

  /**
   * pushes any options with type VAL_TYPE_HYPHENS to the end, i.e., the "--" 
   * are really added at the end.
   * 
   * @see		#VAL_TYPE_HYPHENS
   */
  protected void fixHyphens() {
    NodeList	list;
    Vector<Element>	hyphens;
    int		i;
    Node	node;
    Node	tmpNode;
    boolean	isLast;

    // get all option tags
    list = m_XMLDocument.findNodes("//" + TAG_OPTION);

    // get all hyphen tags
    hyphens = new Vector<Element>();
    for (i = 0; i < list.getLength(); i++) {
      if (((Element) list.item(i)).getAttribute(ATT_TYPE).equals(VAL_TYPE_HYPHENS))
	hyphens.add((Element)list.item(i));
    }

    // check all hyphen tags whether they are at the end, if not fix it
    for (i = 0; i < hyphens.size(); i++) {
      node = (Node) hyphens.get(i);

      // at the end?
      isLast  = true;
      tmpNode = node;
      while (tmpNode.getNextSibling() != null) {
	// normal tag?
	if (tmpNode.getNextSibling().getNodeType() == Node.ELEMENT_NODE) {
	  isLast = false;
	  break;
	}
	tmpNode = tmpNode.getNextSibling();
      }

      // move
      if (!isLast) {
	tmpNode = node.getParentNode();
	tmpNode.removeChild(node);
	tmpNode.appendChild(node);
      }
    }
  }

  /**
   * converts the given node into a command line representation and returns it.
   * 
   * @param parent 	the node to convert to command line
   * @return 		the new command line
   */
  protected String toCommandLine(Element parent) {
    Vector<String>	result;
    Vector		list;
    Vector		subList;
    NodeList		subNodeList;
    String[]		params;
    int			i;
    int			n;
    String		tmpStr;
    
    result = new Vector<String>();
    
    // "options" tag
    if (parent.getNodeName().equals(TAG_OPTIONS)) {
      // children
      list = XMLDocument.getChildTags(parent);

      if (parent.getAttribute(ATT_TYPE).equals(VAL_TYPE_CLASSIFIER)) {
	System.err.println(
	    "Type '" + VAL_TYPE_CLASSIFIER + "' is deprecated, "
	    + "use '" + VAL_TYPE_OPTIONHANDLER + "' instead!");
	parent.setAttribute(ATT_TYPE, VAL_TYPE_OPTIONHANDLER);
      }
      
      if (parent.getAttribute(ATT_TYPE).equals(VAL_TYPE_OPTIONHANDLER)) {
	result.add(parent.getAttribute(ATT_VALUE));

	// hyphens?
	if (    (list.size() > 0)
	     && (parent.getParentNode() != null) 
	     && (parent.getParentNode() instanceof Element) 
	     && (((Element) parent.getParentNode()).getNodeName().equals(TAG_OPTION)) 
	     && (((Element) parent.getParentNode()).getAttribute(ATT_TYPE).equals(VAL_TYPE_HYPHENS)) )
	  result.add("--");
      }

      // process children
      for (i = 0; i < list.size(); i++) {
	tmpStr = toCommandLine((Element) list.get(i));
	try {
	  params = Utils.splitOptions(tmpStr);
	  for (n = 0; n < params.length; n++)
	    result.add(params[n]);
	}
	catch (Exception e) {
	  System.err.println("Error splitting: " + tmpStr);
	  e.printStackTrace();
	}
      }
    }
    // "option" tag
    else if (parent.getNodeName().equals(TAG_OPTION)) {
      subList     = XMLDocument.getChildTags(parent);
      subNodeList = parent.getChildNodes();

      result.add("-" + parent.getAttribute(ATT_NAME));

      // single argument
      if (parent.getAttribute(ATT_TYPE).equals(VAL_TYPE_SINGLE)) {
	if (    (subNodeList.getLength() > 0) 
	     && (subNodeList.item(0).getNodeValue().trim().length() > 0) )
          result.add(subNodeList.item(0).getNodeValue());
      }
      // compound argument surrounded by quotes
      else if (parent.getAttribute(ATT_TYPE).equals(VAL_TYPE_QUOTES)) {
	result.add(toCommandLine((Element) subList.get(0)));
      }
      // classname + further options after "--"
      else if (parent.getAttribute(ATT_TYPE).equals(VAL_TYPE_HYPHENS)) {
	tmpStr = toCommandLine((Element) subList.get(0));
	try {
	  params = Utils.splitOptions(tmpStr);
	  for (n = 0; n < params.length; n++)
	    result.add(params[n]);
	}
	catch (Exception e) {
	  System.err.println("Error splitting: " + tmpStr);
	  e.printStackTrace();
	}
      }
    }
    // other tag
    else {
      System.err.println("Unsupported tag '" + parent.getNodeName() + "' - skipped!");
    }
    
    return Utils.joinOptions(result.toArray(new String[result.size()]));
  }

  /**
   * returns the given DOM document as command line.
   * 
   * @return 		the document as command line
   * @throws Exception 	if anything goes wrong initializing the parsing
   */
  public String toCommandLine() throws Exception {
    return toCommandLine(getDocument().getDocumentElement());
  }

  /**
   * returns the current DOM document as string array.
   * 
   * @return 		the document as string array
   * @throws Exception 	if anything goes wrong initializing the parsing
   */
  public String[] toArray() throws Exception {
    return Utils.splitOptions(toCommandLine());
  }

  /**
   * returns the object in a string representation (as indented XML output).
   * 
   * @return 		the object in a string representation
   */
  public String toString() {
    return getXMLDocument().toString();
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
   * for testing only. prints the given XML file, the resulting commandline and
   * the string array.
   * 
   * @param args	the commandline options.
   * @throws Exception	if something goes wrong
   */
  public static void main(String[] args) throws Exception {
    if (args.length > 0) {
      System.out.println("\nXML:\n\n" + new XMLOptions(args[0]).toString()); 

      System.out.println("\nCommandline:\n\n" + new XMLOptions(args[0]).toCommandLine());

      System.out.println("\nString array:\n");
      String[] options = new XMLOptions(args[0]).toArray();
      for (int i = 0; i < options.length; i++)
	System.out.println(options[i]);
    }
  }
}
