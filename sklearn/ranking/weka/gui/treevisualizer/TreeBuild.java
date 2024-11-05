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
 *    Tree_build.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.treevisualizer;

import java.util.*;
import java.io.*;
import java.awt.*;

/**
 * This class will parse a dotty file and construct a tree structure from it 
 * with Edge's and Node's
 *
 * @author Malcolm Ware (mfw4@cs.waikato.ac.nz)
 * @version $Revision: 1.4 $
 */
public class TreeBuild {
  //this class will parse the tree into relevant strings
  //into info objects 
  //from there it will create the nodes and edges from the info objects

  /** The name of the tree, Not in use. */
  private String m_graphName;

  /** An array with all the nodes initially constructed into it. */
  private Vector m_aNodes;

  /** An array with all the edges initially constructed into it. */
  private Vector m_aEdges;

  /** An array containing a structure that describes the node without 
   * actually creating it. */
  private Vector m_nodes;

  /** An arry containing a structure that describes the edge without 
   * actually creating it. */
  private Vector m_edges;

  /** An object setup to take graph data. */
  private InfoObject m_grObj;

  /** An object setup to take node data. */
  private InfoObject m_noObj;

  /** An object setup to take edge data. */
  private InfoObject m_edObj;

  /** true if it is a digraph. (note that this can't build digraphs). */
  private boolean m_digraph;
  
  /** The stream to parse. */
  private StreamTokenizer m_st;

  /** A table containing all the colors. */
  private Hashtable m_colorTable;
  
  /**
   * Upon construction this will only setup the color table for quick 
   * reference of a color.
   */
  public TreeBuild() {
    m_colorTable = new Hashtable();
    
    Colors ab = new Colors();
    for (int noa = 0;noa < ab.m_cols.length;noa++) {
      m_colorTable.put(ab.m_cols[noa].m_name,ab.m_cols[noa].m_col);
    }
  }  
    
  /**
   * This will build A node structure from the dotty format passed. Don't 
   * send a dotty format with multiple parents
   * per node, and ensure that there is 1 and only 1 node with no parent.
   *
   * @param t The reader with the dotty string to be read.
   * @return The top node of the tree structure (the last node with no parent).
   */
  public Node create(Reader t) {
    m_nodes = new Vector(50,50);
    m_edges = new Vector(50,50);
    m_grObj = new InfoObject("graph");
    m_noObj = new InfoObject("node");
    m_edObj = new InfoObject("edge");
    m_digraph = false;
    
    m_st = new StreamTokenizer(new BufferedReader(t));
    setSyntax();

    graph();

    Node top = generateStructures();
    
    return top;
  }
  
  /**
   * This will go through all the found Nodes and Edges build instances of 
   * these
   * and link them together.
   *
   * @return The node with no parent (the top of the tree).
   */
  private Node generateStructures() {
    String id,label,source,target;
    Integer shape,style;
    int fontsize;
    Color fontcolor,color;

    InfoObject t;
    m_aNodes = new Vector(50,50);
    m_aEdges = new Vector(50,50);
    for (int noa = 0;noa < m_nodes.size();noa++) {
      t = (InfoObject)m_nodes.elementAt(noa);
      id = t.m_id;
      
      if (t.m_label == null) {
	if (m_noObj.m_label == null) {
	  label = "";
	}
	else {
	  label = m_noObj.m_label;
	}
      }
      else {
	label = t.m_label;
      }
      
      if (t.m_shape == null) {
	if (m_noObj.m_shape == null) {
	  shape = new Integer(2);
	}
	else {
	  shape = getShape(m_noObj.m_shape);
	}
      }
      else {
	shape = getShape(t.m_shape);
      }
      if (shape == null) {
	shape = new Integer(2);
      }
      
      if (t.m_style == null) {
	if (m_noObj.m_style == null) {
	  style = new Integer(1);
	}
	else {
	  style = getStyle(m_noObj.m_style);
	}
      }
      else {
	style = getStyle(t.m_style);
      }
      if (style == null) {
	style = new Integer(1);
      }
      
      if (t.m_fontSize == null) {
	if (m_noObj.m_fontSize == null) {
	  fontsize = 12;
	}
	else {
	  fontsize = Integer.valueOf(m_noObj.m_fontSize).intValue();
	}
      }
      else {
	fontsize = Integer.valueOf(t.m_fontSize).intValue();
      }
      
      if (t.m_fontColor == null) {
	if (m_noObj.m_fontColor == null) {
	  fontcolor = Color.black;
	}
	else {
	  fontcolor = (Color)m_colorTable
	    .get(m_noObj.m_fontColor.toLowerCase());
	}
      }
      else {
	fontcolor = (Color)m_colorTable.get(t.m_fontColor.toLowerCase());
      }
      if (fontcolor == null) {
	fontcolor = Color.black;
      }
      
      if (t.m_color == null) {
	if (m_noObj.m_color == null) {
	  color = Color.gray;
	}
	else {
	  color = (Color)m_colorTable.get(m_noObj.m_color.toLowerCase());
	}
      }
      else {
	color = (Color)m_colorTable.get(t.m_color.toLowerCase());
      }
      if (color == null) {
	color = Color.gray;
      }
      
      m_aNodes.addElement(new Node(label,id,style.intValue(),shape.intValue(),
				  color,t.m_data));
    }


    for (int noa = 0;noa < m_edges.size();noa++) {
      t = (InfoObject)m_edges.elementAt(noa);
      id = t.m_id;
      
      if (t.m_label == null) {
	if (m_noObj.m_label == null) {
	  label = "";
	}
	else {
	  label = m_noObj.m_label;
	}
      }
      else {
	label = t.m_label;
      }
      
      if (t.m_shape == null) {
	if (m_noObj.m_shape == null) {
	  shape = new Integer(2);
	}
	else {
	  shape = getShape(m_noObj.m_shape);
	}
      }
      else {
	shape = getShape(t.m_shape);
      }
      if (shape == null) {
	shape = new Integer(2);
      }
      
      if (t.m_style == null) {
	if (m_noObj.m_style == null) {
	  style = new Integer(1);
	}
	else {
	  style = getStyle(m_noObj.m_style);
	}
      }
      else {
	style = getStyle(t.m_style);
      }
      if (style == null) {
	style = new Integer(1);
      }
      
      if (t.m_fontSize == null) {
	if (m_noObj.m_fontSize == null) {
	  fontsize = 12;
	}
	else {
	  fontsize = Integer.valueOf(m_noObj.m_fontSize).intValue();
	}
      }
      else {
	fontsize = Integer.valueOf(t.m_fontSize).intValue();
      }
      
      if (t.m_fontColor == null) {
	if (m_noObj.m_fontColor == null) {
	  fontcolor = Color.black;
	}
	else {
	  fontcolor = (Color)m_colorTable
	    .get(m_noObj.m_fontColor.toLowerCase());
	}
      }
      else {
	fontcolor = (Color)m_colorTable.get(t.m_fontColor.toLowerCase());
      }
      if (fontcolor == null) {
	fontcolor = Color.black;
      }
      
      if (t.m_color == null) {
	if (m_noObj.m_color == null) {
	  color = Color.white;
	}
	else {
	  color = (Color)m_colorTable.get(m_noObj.m_color.toLowerCase());
	}
      }
      else {
	color = (Color)m_colorTable.get(t.m_color.toLowerCase());
      }
      if (color == null) {
	color = Color.white;
      }
      
      m_aEdges.addElement(new Edge(label,t.m_source,t.m_target));
    }
        
    boolean f_set,s_set;
    Node x,sour = null,targ = null;
    Edge y;
    for (int noa = 0;noa < m_aEdges.size();noa++) {
      f_set = false;
      s_set = false;
      y = (Edge)m_aEdges.elementAt(noa);
      for (int nob = 0;nob < m_aNodes.size();nob++) {
	x = (Node)m_aNodes.elementAt(nob);
	if (x.getRefer().equals(y.getRtarget())) {
	  f_set = true;
	  targ = x;
	}
	if (x.getRefer().equals(y.getRsource())) {
	  s_set = true;
	  sour = x;
	}
	if (f_set == true && s_set == true) {
	  break;
	}
      }
      if (targ != sour) {
	y.setTarget(targ);
	y.setSource(sour);
      }
      else {
	System.out.println("logic error");
      }
    }
    
    for (int noa = 0;noa < m_aNodes.size();noa++) {
      if (((Node)m_aNodes.elementAt(noa)).getParent(0) == null) {
	sour = (Node)m_aNodes.elementAt(noa);
      }
    }

    return sour;
  }
  
  /**
   * This will convert the shape string to an int representing that shape.
   *
   * @param sh The name of the shape.
   * @return An Integer representing the shape.
   */
  private Integer getShape(String sh) {
    if (sh.equalsIgnoreCase("box") || sh.equalsIgnoreCase("rectangle")) {
      return new Integer(1);
    }
    else if (sh.equalsIgnoreCase("oval")) {
      return new Integer(2);
    }
    else if (sh.equalsIgnoreCase("diamond")) {
      return new Integer(3);
    }
    else {
      return null;
    }
  }
  
  /**
   * Converts the string representing the fill style int oa number 
   * representing it.
   *
   * @param sty The name of the style.
   * @return An Integer representing the shape.
   */
  private Integer getStyle(String sty) {
    if (sty.equalsIgnoreCase("filled")) {
      return new Integer(1);
    }
    else {
      return null;
    }
  }

  /**
   * This will setup the syntax for the tokenizer so that it parses the 
   * string properly.
   *
   */
  private void setSyntax() {
    m_st.resetSyntax();
    m_st.eolIsSignificant(false);
    m_st.slashStarComments(true);
    m_st.slashSlashComments(true);
    //System.out.println("slash");
    m_st.whitespaceChars(0,' ');
    m_st.wordChars(' '+1,'\u00ff');
    m_st.ordinaryChar('[');
    m_st.ordinaryChar(']');
    m_st.ordinaryChar('{');
    m_st.ordinaryChar('}');
    m_st.ordinaryChar('-');
    m_st.ordinaryChar('>');
    m_st.ordinaryChar('/');
    m_st.ordinaryChar('*');
    m_st.quoteChar('"');
    m_st.whitespaceChars(';',';');
    m_st.ordinaryChar('=');
  }

  /**
   * This is the alternative syntax for the tokenizer.
   */
  private void alterSyntax() {
    m_st.resetSyntax();
    m_st.wordChars('\u0000', '\u00ff');
    m_st.slashStarComments(false);
    m_st.slashSlashComments(false);
    m_st.ordinaryChar('\n');
    m_st.ordinaryChar('\r');
  }

  /**
   * This will parse the next token out of the stream and check for certain 
   * conditions.
   *
   * @param r The error string to print out if something goes wrong.
   */
  private void nextToken(String r) {
    int t = 0;
    try {
      t = m_st.nextToken();
    } catch(IOException e) {
    }
    
    if (t == m_st.TT_EOF) {
      System.out.println("eof , " + r);
    }
    else if (t == m_st.TT_NUMBER) {
      System.out.println("got a number , " + r);
    }
  }
  
  /**
   * Parses the top of the dotty stream that has the graph information.
   *
   */
  private void graph() {
    boolean flag = true;
    String s;
    //expect digraph
    int t;
    nextToken("expected 'digraph'");
    
    if (m_st.sval.equalsIgnoreCase("digraph")) {
      m_digraph = true;
    }
    else {
      System.out.println("expected 'digraph'");
    }
    
    nextToken("expected a Graph Name");
    if (m_st.sval != null) {
      m_graphName = m_st.sval;
    }
    else {
      System.out.println("expected a Graph Name");
    }
    
    nextToken("expected '{'");
    
    if (m_st.ttype == '{') {
      stmtList();
    }
    else {
      System.out.println("expected '{'");
    }
  }

  /**
   * This is one of the states, this one is where new items can be defined 
   * or the structure can end.
   *
   */
  private void stmtList() {
    boolean flag = true;
    String s;
    int t;
    while(flag) {
      nextToken("expects a STMT_LIST item or '}'");
      if (m_st.ttype == '}') {
	flag = false;
      }
      else if (m_st.sval.equalsIgnoreCase("graph") ||
	       m_st.sval.equalsIgnoreCase("node") ||
	       m_st.sval.equalsIgnoreCase("edge")) {
	m_st.pushBack();
	attrStmt();
      }
      else if (m_st.sval != null) {
	nodeId(m_st.sval,0);
      }
      else {
	System.out.println("expects a STMT_LIST item or '}'");
      }
    }
  }

  /**
   * This will deal specifically with a new object such as graph , node , edge.
   *
   */
  private void attrStmt() {
    
    nextToken("expected 'graph' or 'node' or 'edge'");
    
    if (m_st.sval.equalsIgnoreCase("graph")) {
      nextToken("expected a '['");
      if (m_st.ttype == '[') {
	attrList(m_grObj);
      }
      else {
	System.out.println("expected a '['");
      }
    }
    else if (m_st.sval.equalsIgnoreCase("node")) {
      nextToken("expected a '['");
      if (m_st.ttype == '[') {
	attrList(m_noObj);
      }
      else {
	System.out.println("expected a '['");
      }
    }
    else if (m_st.sval.equalsIgnoreCase("edge")) {
      nextToken("expected a '['");
      if (m_st.ttype == '[') {
	attrList(m_edObj);
      }
      else {
	System.out.println("expected a '['");
      }
      
    }
    else {
      System.out.println("expected 'graph' or 'node' or 'edge'"); 
    }
  }

  /**
   * Generates a new InfoObject with the specified name and either does 
   * further processing if applicable
   * Otherwise it is an edge and will deal with that.
   *
   * @param s The ID string.
   * @param t Not sure!.
   */
  private void nodeId(String s,int t) {
    
    nextToken("error occurred in node_id");

    if (m_st.ttype == '}') {
      //creates a node if t is zero
      if (t == 0) {
	m_nodes.addElement(new InfoObject(s));
      }
      m_st.pushBack();
    }
    else if (m_st.ttype == '-') {
      nextToken("error occurred checking for an edge");
      if (m_st.ttype == '>') {
	edgeStmt(s);
      }
      else {
	System.out.println("error occurred checking for an edge");
      }
    }
    else if (m_st.ttype == '[') {
      //creates a node if t is zero and sends it to attr
      if (t == 0) {
	m_nodes.addElement(new InfoObject(s));
	attrList((InfoObject)m_nodes.lastElement());
      }
      else {
	attrList((InfoObject)m_edges.lastElement());
      }
    }
    else if (m_st.sval != null) {
      //creates a node if t is zero
      if (t == 0) {
	m_nodes.addElement(new InfoObject(s));
      }
      m_st.pushBack();
    }
    else {
      System.out.println("error occurred in node_id");
    }
  }

  /**
   * This will get the target of the edge.
   *
   * @param i The source of the edge.
   */
  private void edgeStmt(String i) {
    nextToken("error getting target of edge");
    
    if (m_st.sval != null) {
      m_edges.addElement(new InfoObject("an edge ,no id"));
      ((InfoObject)m_edges.lastElement()).m_source = i;
      ((InfoObject)m_edges.lastElement()).m_target = m_st.sval;
      nodeId(m_st.sval,1);
    }
    else {
      System.out.println("error getting target of edge");
    }
  }

  /**
   * This will parse all the items in the attrib list for an object.
   *
   * @param a The object that the attribs apply to.
   */
  private void attrList(InfoObject a) {
    boolean flag = true;
    
    while (flag) {
      nextToken("error in attr_list");
      //System.out.println(st.sval);
      if (m_st.ttype == ']') {
	flag = false;
      }
      else if (m_st.sval.equalsIgnoreCase("color")) {
	nextToken("error getting color");
	if (m_st.ttype == '=') {
	  nextToken("error getting color");
	  if (m_st.sval != null) {
	    a.m_color = m_st.sval;
	  }
	  else {
	    System.out.println("error getting color");
	  }
	}
	else {
	  System.out.println("error getting color");
	}
      }
      else if (m_st.sval.equalsIgnoreCase("fontcolor")) {
	nextToken("error getting font color");
	if (m_st.ttype == '=') {
	  nextToken("error getting font color");
	  if (m_st.sval != null) {
	    a.m_fontColor = m_st.sval;
	  }
	  else {
	    System.out.println("error getting font color");
	  }
	}
	else {
	  System.out.println("error getting font color");
	}
      }
      else if (m_st.sval.equalsIgnoreCase("fontsize")) {
	nextToken("error getting font size");
	if (m_st.ttype == '=') {
	  nextToken("error getting font size");
	  if (m_st.sval != null) {
	    a.m_fontSize = m_st.sval;
	  }
	  else {
	    System.out.println("error getting font size");
	  }
	}
	else {
	  System.out.println("error getting font size");
	}
      }
      else if (m_st.sval.equalsIgnoreCase("label")) {
	nextToken("error getting label");
	if (m_st.ttype == '=') {
	  nextToken("error getting label");
	  if (m_st.sval != null) {
	    a.m_label = m_st.sval;
	  }
	  else {
	    System.out.println("error getting label");
	  }
	}
	else {
	  System.out.println("error getting label");
	}
      }
      else if (m_st.sval.equalsIgnoreCase("shape")) {
	nextToken("error getting shape");
	if (m_st.ttype == '=') {
	  nextToken("error getting shape");
	  if (m_st.sval != null) {
	    a.m_shape = m_st.sval;
	  }
	  else {
	    System.out.println("error getting shape");
	  }
	}
	else {
	  System.out.println("error getting shape");
	}
      }
      else if (m_st.sval.equalsIgnoreCase("style")) {
	nextToken("error getting style");
	if (m_st.ttype == '=') {
	  nextToken("error getting style");
	  if (m_st.sval != null) {
	    a.m_style = m_st.sval;
	  }
	  else {
	    System.out.println("error getting style");
	  }
	}
	else {
	  System.out.println("error getting style");
	}
      }
      else if (m_st.sval.equalsIgnoreCase("data")) {
	nextToken("error getting data");
	if (m_st.ttype == '=') {
	  //data has a special data string that can have anything
	  //this is delimited by a single comma on an otherwise empty line
	  alterSyntax();
	  a.m_data = new String("");
	  
	  while (true) {
	    nextToken("error getting data");
	    if (m_st.sval != null && a.m_data 
		!= null && m_st.sval.equals(",")) {
	      break;
	    }
	    else if (m_st.sval != null) {
	      a.m_data = a.m_data.concat(m_st.sval);
	    }
	    else if (m_st.ttype == '\r') {
	      a.m_data = a.m_data.concat("\r");
	    }
	    else if (m_st.ttype == '\n') {
	      a.m_data = a.m_data.concat("\n");
	    }
	    else {
	      System.out.println("error getting data");
	    }
	  }
	  setSyntax();
	}
	else {
	  System.out.println("error getting data");
	}
      }
    }
  }

  //special class for use in creating the tree

  /**
   * This is an internal class used to keep track of the info for the objects 
   * before they are 
   * actually created.
   */
  private class InfoObject {
    /** The ID string for th object. */
    public String m_id;

    /** The color name for the object. */
    public String m_color;

    /** The font color for the object. not in use. */
    public String m_fontColor;

    /** The fontsize for the object. not in use. */
    public String m_fontSize;

    /** The label for the object. */
    public String m_label;

    /** The shape name of for the object. */
    public String m_shape;

    /** The backstyle name for the object. */
    public String m_style;

    /** The source ID of the object. */
    public String m_source;

    /** The target ID of the object. */
    public String m_target;

    /** The data for this object. */
    public String m_data;

    /**
     * This will construct a new InfoObject with the specified ID string.
     */
    public InfoObject(String i) {
      m_id = i;
      m_color = null;
      m_fontColor = null;
      m_fontSize = null;
      m_label = null;
      m_shape = null;
      m_style = null;
      m_source = null;
      m_target = null;
      m_data = null;
    }
  }
}
