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
 *    DotParser.java
 *    Copyright (C) 2003 University of Waikato, Hamilton, New Zealand
 *
 */
package weka.gui.graphvisualizer;

import java.io.Reader;
import java.io.StreamTokenizer;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.FileWriter;

import weka.core.FastVector;
import weka.gui.graphvisualizer.GraphNode;
import weka.gui.graphvisualizer.GraphEdge;



/**
 * This class parses input in DOT format, and
 * builds the datastructures that are passed to it.
 * It is NOT 100% compatible with the DOT format. The
 * GraphNode and GraphEdge classes do not have any provision
 * for dealing with different colours or  shapes of nodes,
 * there can however, be a different label and ID for a
 * node. It  also does not do anything for labels for
 * edges. However, this class won't crash or throw an
 * exception if it encounters any of the above
 * attributes of an edge or a node. This class however,
 * won't be able to deal with things like subgraphs and
 * grouping of nodes.
 *
 * @author Ashraf M. Kibriya (amk14@cs.waikato.ac.nz)
 * @version $Revision: 1.4 $ - 23 Apr 2003 - Initial version (Ashraf M. Kibriya)
 */
public class DotParser implements GraphConstants  {
  
  /** These holds the nodes and edges of the graph */
  protected FastVector m_nodes, m_edges;
  /** This is the input containing DOT stream  to be parsed */
  protected Reader m_input;
  /**  This holds the name of the graph if there is any otherwise it is null */
  protected String m_graphName;
  
  /**
   *
   *  Dot parser Constructor
   *
   * @param input - The input, if passing in a string then
   *                encapsulate that in String reader object
   * @param nodes - Vector to put in GraphNode objects,
   *                corresponding to the nodes parsed in from
   *                the input
   * @param edges - Vector to put in GraphEdge objects,
   *                corresponding to the edges parsed in from
   *                the input
   */
  public DotParser(Reader input, FastVector nodes, FastVector edges) {
    m_nodes = nodes; m_edges = edges;
    m_input = input;
  }
  
  
  /**
   * This method parses the string or the InputStream that we
   * passed in through the constructor and builds up the
   * m_nodes and m_edges vectors
   *
   * @return - returns the name of the graph
   */
  public String parse() {
    StreamTokenizer tk = new StreamTokenizer(new BufferedReader(m_input));
    setSyntax(tk);
    
    graph(tk);
    
    return m_graphName;
  }
  
  
  /**
   * This method sets the syntax of the StreamTokenizer.
   * i.e. set the whitespace, comment and delimit chars.
   *
   */
  protected void setSyntax(StreamTokenizer tk) {
    tk.resetSyntax();
    tk.eolIsSignificant(false);
    tk.slashStarComments(true);
    tk.slashSlashComments(true);
    tk.whitespaceChars(0,' ');
    tk.wordChars(' '+1,'\u00ff');
    tk.ordinaryChar('[');
    tk.ordinaryChar(']');
    tk.ordinaryChar('{');
    tk.ordinaryChar('}');
    tk.ordinaryChar('-');
    tk.ordinaryChar('>');
    tk.ordinaryChar('/');
    tk.ordinaryChar('*');
    tk.quoteChar('"');
    tk.whitespaceChars(';',';');
    tk.ordinaryChar('=');
  }
  
  /*************************************************************
   *
   * Following methods parse the DOT input and mimic the DOT
   * language's grammar in their structure
   *
   *************************************************************
   */
  protected void graph(StreamTokenizer tk) {
    try {
      tk.nextToken();
      
      if(tk.ttype==tk.TT_WORD) {
        if(tk.sval.equalsIgnoreCase("digraph")) {
          tk.nextToken();
          if(tk.ttype==tk.TT_WORD) {
            m_graphName = tk.sval;
            tk.nextToken();
          }
          
          while(tk.ttype!='{') {
            System.err.println("Error at line "+tk.lineno()+" ignoring token "+
            tk.sval);
            tk.nextToken();
            if(tk.ttype==tk.TT_EOF)
              return;
          }
          stmtList(tk);
        }
        else if(tk.sval.equalsIgnoreCase("graph"))
          System.err.println("Error. Undirected graphs cannot be used");
        else
          System.err.println("Error. Expected graph or digraph at line "+
          tk.lineno());
      }
    }
    catch(Exception ex) { ex.printStackTrace(); }
    
    
    //int tmpMatrix[][] = new int[m_nodes.size()][m_nodes.size()];
    //for(int i=0; i<m_edges.size(); i++)
    //    tmpMatrix[((GraphEdge)m_edges.elementAt(i)).src]
    //	     [((GraphEdge)m_edges.elementAt(i)).dest] =
    //                                   ((GraphEdge)m_edges.elementAt(i)).type;
    //for(int i=0; i<m_nodes.size(); i++) {
    //    GraphNode n = (GraphNode)m_nodes.elementAt(i);
    //    n.edges = tmpMatrix[i];
    //}
    
    //Adding parents, and those edges to a node which are coming out from it
    int tmpEdges[], noOfEdgesOfNode[]=new int[m_nodes.size()];
    int noOfPrntsOfNode[]=new int[m_nodes.size()];
    for(int i=0; i<m_edges.size(); i++) {
      GraphEdge e = (GraphEdge)m_edges.elementAt(i);
      noOfEdgesOfNode[e.src]++;
      noOfPrntsOfNode[e.dest]++;
    }
    for(int i=0; i<m_edges.size(); i++) {
      GraphEdge e  = (GraphEdge)m_edges.elementAt(i);
      GraphNode n  = (GraphNode)m_nodes.elementAt(e.src);
      GraphNode n2 = (GraphNode)m_nodes.elementAt(e.dest);
      if(n.edges==null) {
        n.edges = new int[noOfEdgesOfNode[e.src]][2];
        for(int k=0; k<n.edges.length; k++)
          n.edges[k][1]=0;
      }
      if(n2.prnts==null) {
        n2.prnts = new int[noOfPrntsOfNode[e.dest]];
        for(int k=0; k<n2.prnts.length; k++)
          n2.prnts[k]=-1;
      }
      int k=0;
      while(n.edges[k][1]!=0) k++;
      n.edges[k][0] = e.dest;
      n.edges[k][1] = e.type;
      
      k=0;
      while(n2.prnts[k]!=-1) k++;
      n2.prnts[k] = e.src;
    }
  }
  
  
  protected void stmtList(StreamTokenizer tk) throws Exception{
    tk.nextToken();
    if(tk.ttype=='}' || tk.ttype==tk.TT_EOF)
      return;
    else {
      stmt(tk);
      stmtList(tk);
    }
  }
  
  
  protected void stmt(StreamTokenizer tk) {
    //tk.nextToken();
    
    if(tk.sval.equalsIgnoreCase("graph") || tk.sval.equalsIgnoreCase("node") ||
    tk.sval.equalsIgnoreCase("edge") )
      ; //attribStmt(k);
    else {
      try {
        nodeID(tk);
        int nodeindex= m_nodes.indexOf(new GraphNode(tk.sval, null));
        tk.nextToken();
        
        if(tk.ttype == '[')
          nodeStmt(tk, nodeindex);
        else if(tk.ttype == '-')
          edgeStmt(tk, nodeindex);
        else
          System.err.println("error at lineno "+tk.lineno()+" in stmt");
      }
      catch(Exception ex) {
        System.err.println("error at lineno "+tk.lineno()+" in stmtException");
        ex.printStackTrace();
      }
    }
  }
  
  
  protected void nodeID(StreamTokenizer tk) throws Exception{
    
    if(tk.ttype=='"' || tk.ttype==tk.TT_WORD || (tk.ttype>='a' && tk.ttype<='z')
    || (tk.ttype>='A' && tk.ttype<='Z')) {
      if(m_nodes!=null && !(m_nodes.contains( new GraphNode(tk.sval, null))) ) {
        m_nodes.addElement( new GraphNode(tk.sval, tk.sval) );
        //System.out.println("Added node >"+tk.sval+"<");
      }
    }
    else
    { throw new Exception(); }
    
    //tk.nextToken();
  }
  
  
  protected void nodeStmt(StreamTokenizer tk, final int nindex)
  throws Exception {
    tk.nextToken();
    
    GraphNode temp = (GraphNode) m_nodes.elementAt(nindex);
    
    if(tk.ttype==']' || tk.ttype==tk.TT_EOF)
      return;
    else if(tk.ttype==tk.TT_WORD) {
      
      if(tk.sval.equalsIgnoreCase("label")) {
        
        tk.nextToken();
        if(tk.ttype=='=') {
          tk.nextToken();
          if(tk.ttype==tk.TT_WORD || tk.ttype=='"')
            temp.lbl = tk.sval;
          else {
            System.err.println("couldn't find label at line "+tk.lineno());
            tk.pushBack();
          }
        }
        else {
          System.err.println("couldn't find label at line "+tk.lineno());
          tk.pushBack();
        }
      }
      
      else if(tk.sval.equalsIgnoreCase("color")){
        
        tk.nextToken();
        if(tk.ttype=='=') {
          tk.nextToken();
          if(tk.ttype==tk.TT_WORD || tk.ttype=='"')
            ;
          else {
            System.err.println("couldn't find color at line "+tk.lineno());
            tk.pushBack();
          }
        }
        else {
          System.err.println("couldn't find color at line "+tk.lineno());
          tk.pushBack();
        }
      }
      
      else if(tk.sval.equalsIgnoreCase("style")) {
        
        tk.nextToken();
        if(tk.ttype=='=') {
          tk.nextToken();
          if(tk.ttype==tk.TT_WORD || tk.ttype=='"')
            ;
          else {
            System.err.println("couldn't find style at line "+tk.lineno());
            tk.pushBack();
          }
        }
        else {
          System.err.println("couldn't find style at line "+tk.lineno());
          tk.pushBack();
        }
      }
    }
    nodeStmt(tk, nindex);
  }
  
  
  protected void edgeStmt(StreamTokenizer tk, final int nindex)
  throws Exception {
    tk.nextToken();
    
    GraphEdge e=null;
    if(tk.ttype=='>') {
      tk.nextToken();
      if(tk.ttype=='{') {
        while(true) {
          tk.nextToken();
          if(tk.ttype=='}')
            break;
          else {
            nodeID(tk);
            e = new GraphEdge(nindex,
            m_nodes.indexOf( new GraphNode(tk.sval, null) ),
            DIRECTED);
            if( m_edges!=null && !(m_edges.contains(e)) ) {
              m_edges.addElement( e );
              //System.out.println("Added edge from "+
              //                  ((GraphNode)(m_nodes.elementAt(nindex))).ID+
              //                  " to "+
              //	        ((GraphNode)(m_nodes.elementAt(e.dest))).ID);
            }
          }
        }
      }
      else {
        nodeID(tk);
        e = new GraphEdge(nindex,
        m_nodes.indexOf( new GraphNode(tk.sval, null) ),
        DIRECTED);
        if( m_edges!=null && !(m_edges.contains(e)) ) {
          m_edges.addElement( e );
          //System.out.println("Added edge from "+
          //                 ((GraphNode)(m_nodes.elementAt(nindex))).ID+" to "+
          //		     ((GraphNode)(m_nodes.elementAt(e.dest))).ID);
        }
      }
    }
    else if(tk.ttype=='-') {
      System.err.println("Error at line "+tk.lineno()+
      ". Cannot deal with undirected edges");
      if(tk.ttype==tk.TT_WORD)
        tk.pushBack();
      return;
    }
    else {
      System.err.println("Error at line "+tk.lineno()+" in edgeStmt");
      if(tk.ttype==tk.TT_WORD)
        tk.pushBack();
      return;
    }
    
    tk.nextToken();
    
    if(tk.ttype=='[')
      edgeAttrib(tk, e);
    else
      tk.pushBack();
  }
  
  
  protected void edgeAttrib(StreamTokenizer tk, final GraphEdge e)
  throws Exception {
    tk.nextToken();
    
    if(tk.ttype==']' || tk.ttype==tk.TT_EOF)
      return;
    else if(tk.ttype==tk.TT_WORD) {
      
      if(tk.sval.equalsIgnoreCase("label")) {
        
        tk.nextToken();
        if(tk.ttype=='=') {
          tk.nextToken();
          if(tk.ttype==tk.TT_WORD || tk.ttype=='"')
            System.err.println("found label "+tk.sval);//e.lbl = tk.sval;
          else {
            System.err.println("couldn't find label at line "+tk.lineno());
            tk.pushBack();
          }
        }
        else {
          System.err.println("couldn't find label at line "+tk.lineno());
          tk.pushBack();
        }
      }
      else if(tk.sval.equalsIgnoreCase("color")) {
        
        tk.nextToken();
        if(tk.ttype=='=') {
          tk.nextToken();
          if(tk.ttype==tk.TT_WORD || tk.ttype=='"')
            ;
          else {
            System.err.println("couldn't find color at line "+tk.lineno());
            tk.pushBack();
          }
        }
        else {
          System.err.println("couldn't find color at line "+tk.lineno());
          tk.pushBack();
        }
      }
      
      else if(tk.sval.equalsIgnoreCase("style")) {
        
        tk.nextToken();
        if(tk.ttype=='=') {
          tk.nextToken();
          if(tk.ttype==tk.TT_WORD || tk.ttype=='"')
            ;
          else {
            System.err.println("couldn't find style at line "+tk.lineno());
            tk.pushBack();
          }
        }
        else {
          System.err.println("couldn't find style at line "+tk.lineno());
          tk.pushBack();
        }
      }
    }
    edgeAttrib(tk, e);
  }
  
  
  /**
   *
   * This method saves a graph in a file in DOT format.
   * However, if reloaded in GraphVisualizer we would need
   * to layout the graph again.
   *
   * @param filename  - The name of the file to write in. (will overwrite)
   * @param graphName - The name of the graph
   * @param nodes     - Vector containing all the nodes
   * @param edges     - Vector containing all the edges
   */
  public static void writeDOT(String filename, String graphName,
  FastVector nodes, FastVector edges) {
    try {
      FileWriter os = new FileWriter(filename);
      os.write("digraph ", 0, ("digraph ").length());
      if(graphName!=null)
        os.write(graphName+" ", 0, graphName.length()+1);
      os.write("{\n", 0, ("{\n").length());
      
      GraphEdge e;
      for(int i=0; i<edges.size(); i++) {
        e = (GraphEdge) edges.elementAt(i);
        os.write(((GraphNode)nodes.elementAt(e.src)).ID, 0,
        ((GraphNode)nodes.elementAt(e.src)).ID.length());
        os.write("->", 0, ("->").length() );
        os.write(((GraphNode)nodes.elementAt(e.dest)).ID+"\n",
        0,
        ((GraphNode)nodes.elementAt(e.dest)).ID.length()+1);
      }
      os.write("}\n", 0, ("}\n").length());
      os.close();
    }
    catch(IOException ex) { ex.printStackTrace(); }
  }
  
} // DotParser
