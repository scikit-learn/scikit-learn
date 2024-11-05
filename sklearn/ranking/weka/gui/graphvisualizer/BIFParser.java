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
 *    BIFParser.java
 *    Copyright (C) 2003 University of Waikato, Hamilton, New Zealand
 *
 */
package weka.gui.graphvisualizer;

import java.io.InputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.StringReader;
import java.util.StringTokenizer;

import weka.core.FastVector;

import org.w3c.dom.Document;
import org.w3c.dom.NodeList;
import org.w3c.dom.Element;


/**
 * This class parses an inputstream or a string in
 * XMLBIF ver. 0.3 format, and builds the datastructures
 * that are passed to it through the constructor.
 *
 * @author Ashraf M. Kibriya (amk14@cs.waikato.ac.nz)
 * @version $Revision: 1.7 $ - 24 Apr 2003 - Initial version (Ashraf M. Kibriya)
 */
public class BIFParser implements GraphConstants {
  
  /** These holds the nodes and edges of the graph */
  protected FastVector m_nodes, m_edges;
  /**  This holds the name of the graph (i.e. the name of network tag in XMLBIF
   * input)
   */
  protected String graphName;
  /** This holds the string to be parsed */
  protected String inString;
  /** This holds the InputStream to be parsed */
  protected InputStream inStream;
  
  
  /**
   * Constructor (if our input is a String)
   *
   * @param input the string to be parsed (should not be null)
   * @param nodes vector containing GraphNode objects (should be empty)
   * @param edges vector containing GraphEdge objects (should be empty)
   */
  public BIFParser(String input, FastVector nodes, FastVector edges) {
    m_nodes = nodes; m_edges = edges; inString = input;
  }
  
  
  /**
   * Constructor (if our input is an InputStream)
   *
   * @param instream the InputStream to be parsed (should not be null)
   * @param nodes vector containing GraphNode objects (should be empty)
   * @param edges vector containing GraphEdge objects (should be empty)
   */
  public BIFParser(InputStream instream, FastVector nodes, FastVector edges) {
    m_nodes = nodes; m_edges = edges; inStream = instream;
  }
  
  
  /**
   * This method parses the string or the InputStream that we
   * passed in through the constructor and builds up the
   * m_nodes and m_edges vectors
   * @exception Exception if both the inString and inStream are
   *              null, i.e. no input has been provided
   * @exception BIFFormatException if there is format of the
   *              input is not correct. The format should conform to
   *              XMLBIF version 0.3
   * @exception NumberFormatException if there is an invalid
   *              char in the probability table of a node.
   * @return    returns the name of the graph
   */
  public String parse() throws Exception {
    Document dc=null;
    
    javax.xml.parsers.DocumentBuilderFactory dbf =
    javax.xml.parsers.DocumentBuilderFactory.newInstance();
    dbf.setIgnoringElementContentWhitespace(true);
    javax.xml.parsers.DocumentBuilder db = dbf.newDocumentBuilder();
    
    if(inStream!=null)
      dc = db.parse(inStream);
    else if(inString!=null)
      dc = db.parse(new org.xml.sax.InputSource(new StringReader(inString)));
    else
    { throw new Exception("No input given"); }
    
    NodeList nl = dc.getElementsByTagName( "NETWORK" );
    
    if(nl.getLength()==0) {
      throw new BIFFormatException( "NETWORK tag not found" );
    }
    
    //take only the first network node
    NodeList templist = ((Element)nl.item(0)).getElementsByTagName( "NAME" );
    graphName = templist.item(0).getFirstChild().getNodeValue();
    //System.out.println("The name of the network is "+
    //templist.item(0).getFirstChild().getNodeValue());
    
    //Get all the variables
    nl = dc.getElementsByTagName("VARIABLE");
    for(int i=0; i<nl.getLength(); i++) {
      
      templist = ((Element)nl.item(i)).getElementsByTagName("NAME");
      if(templist.getLength()>1)
        throw new BIFFormatException("More than one name tags found for "+
        "variable no. "+(i+1));
      
      String nodename =
      ((org.w3c.dom.Node)templist.item(0)).getFirstChild().getNodeValue();
      GraphNode n = new GraphNode( nodename, nodename, GraphNode.NORMAL );
      m_nodes.addElement(n);
      //getting nodes position
      templist = ((Element)nl.item(i)).getElementsByTagName("PROPERTY");
      for(int j=0; j<templist.getLength(); j++) {
        if( ((org.w3c.dom.Node)templist.item(j)).getFirstChild()
        .getNodeValue().startsWith("position") ) {
          String xy = templist.item(j).getFirstChild().getNodeValue();
          //System.out.println("x: "+
          //                   xy.substring(xy.indexOf('(')+1, xy.indexOf(','))+
          //                   " y: "+
          //                   xy.substring(xy.indexOf(',')+1, xy.indexOf(')'))
          //                  );
          n.x = Integer.parseInt( xy.substring(xy.indexOf('(')+
          1, xy.indexOf(',')).trim() );
          n.y = Integer.parseInt( xy.substring(xy.indexOf(',')+
          1, xy.indexOf(')')).trim() );
          break;
        }
      }
      //getting all the outcomes of the node
      templist = ((Element)nl.item(i)).getElementsByTagName("OUTCOME");
      n.outcomes = new String[templist.getLength()];
      for(int j=0; j<templist.getLength(); j++) {
        n.outcomes[j] = templist.item(j).getFirstChild().getNodeValue();
        //System.out.println("Outcome["+j+"]: "+n.outcomes[j]);
      }
    } //end for (for variables)
    
    //Get all the edges and probability tables by getting all the definitions
    nl = dc.getElementsByTagName("DEFINITION");
    for(int i=0; i<nl.getLength(); i++) {
      
      templist = ((Element)nl.item(i)).getElementsByTagName("FOR");
      //the Label of the node the edges are coming into
      String nid = templist.item(0).getFirstChild().getNodeValue();
      
      //getting the GraphNode object with the above label
      GraphNode n = (GraphNode)m_nodes.elementAt(0);
      for(int j=1; j<m_nodes.size() && !n.ID.equals(nid); j++)
        n = (GraphNode)m_nodes.elementAt(j);
      
      
      templist = ((Element)nl.item(i)).getElementsByTagName("GIVEN");
      int parntOutcomes = 1; //for creating the probability table later on
      //creating all the edges coming into the node
      for(int j=0; j<templist.getLength(); j++) {
        nid = templist.item(j).getFirstChild().getNodeValue();
        
        GraphNode n2 = (GraphNode)m_nodes.elementAt(0);
        for(int k=1; k<m_nodes.size() && !n2.ID.equals(nid); k++)
          n2 = (GraphNode)m_nodes.elementAt(k);
        m_edges.addElement( new GraphEdge(m_nodes.indexOf(n2),
        m_nodes.indexOf(n), 1) );
        
        parntOutcomes *= n2.outcomes.length;
      }
      
      //creating the probability table for the node
      templist = ((Element)nl.item(i)).getElementsByTagName("TABLE");
      if(templist.getLength()>1)
        throw new BIFFormatException("More than one Probability Table for "+
        n.ID);
      
      String probs = templist.item(0).getFirstChild().getNodeValue();
      StringTokenizer tk = new StringTokenizer(probs, " \n\t");
      
      if(parntOutcomes*n.outcomes.length > tk.countTokens())
        throw new BIFFormatException("Probability Table for "+n.ID+
        " contains more values than it should");
      else if(parntOutcomes*n.outcomes.length < tk.countTokens())
        throw new BIFFormatException("Probability Table for "+n.ID+
        " contains less values than it should");
      else {
        n.probs = new double[parntOutcomes][n.outcomes.length];
        for(int r=0; r<parntOutcomes; r++)     //row
          for(int c=0; c<n.outcomes.length; c++) //column
            try {
              n.probs[r][c] =  Double.parseDouble( tk.nextToken() );
            }
            catch(NumberFormatException ne) { throw ne; }
      } // end of creating probability table
    } //endfor (for edges)
    
    
    //int tmpMatrix[][] = new int[m_nodes.size()][m_nodes.size()];
    //for(int i=0; i<m_edges.size(); i++)
    //    tmpMatrix[((GraphEdge)m_edges.elementAt(i)).src]
    //	       [((GraphEdge)m_edges.elementAt(i)).dest] =
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
          n.edges[k][0]=-1;
      }
      if(n2.prnts==null) {
        n2.prnts = new int[noOfPrntsOfNode[e.dest]];
        for(int k=0; k<n2.prnts.length; k++)
          n2.prnts[k]=-1;
      }
      
      int k=0;
      while(n.edges[k][0]!=-1) k++;
      n.edges[k][0] = e.dest;
      n.edges[k][1] = e.type;
      
      k=0;
      while(n2.prnts[k]!=-1) k++;
      n2.prnts[k] = e.src;
    }
    
    //processGraph();
    //setAppropriateSize();
    return graphName;
  } //end readBIF
  
  
  /**
   * This method writes a graph in XMLBIF ver. 0.3 format to a file.
   * However, if is reloaded in GraphVisualizer we would need to layout
   * the graph again to display it correctly.
   *
   * @param filename  The name of the file to write in. (will overwrite)
   * @param graphName The name of the graph. (will be the name of network
   * tag in XMLBIF)
   * @param nodes     Vector containing all the nodes
   * @param edges     Vector containing all the edges
   */
  public static void writeXMLBIF03(String filename, String graphName,
  FastVector nodes, FastVector edges) {
    try {
      FileWriter outfile = new FileWriter(filename);
      
      StringBuffer text = new StringBuffer();
      
      text.append("<?xml version=\"1.0\"?>\n");
      text.append("<!-- DTD for the XMLBIF 0.3 format -->\n");
      text.append("<!DOCTYPE BIF [\n");
      text.append("	<!ELEMENT BIF ( NETWORK )*>\n");
      text.append("	      <!ATTLIST BIF VERSION CDATA #REQUIRED>\n");
      text.append("	<!ELEMENT NETWORK ( NAME, ( PROPERTY | VARIABLE | DEFI"+
      "NITION )* )>\n");
      text.append("	<!ELEMENT NAME (#PCDATA)>\n");
      text.append("	<!ELEMENT VARIABLE ( NAME, ( OUTCOME |  PROPERTY )* )"+
      " >\n");
      text.append("	      <!ATTLIST VARIABLE TYPE (nature|decision|utility"+
      ") \"nature\">\n");
      text.append("	<!ELEMENT OUTCOME (#PCDATA)>\n");
      text.append("	<!ELEMENT DEFINITION ( FOR | GIVEN | TABLE | PROPERTY"+
      " )* >\n");
      text.append("	<!ELEMENT FOR (#PCDATA)>\n");
      text.append("	<!ELEMENT GIVEN (#PCDATA)>\n");
      text.append("	<!ELEMENT TABLE (#PCDATA)>\n");
      text.append("	<!ELEMENT PROPERTY (#PCDATA)>\n");
      text.append("]>\n");
      text.append("\n");
      text.append("\n");
      text.append("<BIF VERSION=\"0.3\">\n");
      text.append("<NETWORK>\n");
      text.append("<NAME>" + XMLNormalize(graphName)  + "</NAME>\n");
      
      //Writing all the node names and their outcomes
      //If outcome is null(ie if the graph was loaded from DOT file) then
      //simply write TRUE
      for(int nodeidx=0; nodeidx<nodes.size(); nodeidx++) {
        GraphNode n = (GraphNode)nodes.elementAt(nodeidx);
        if(n.nodeType!=GraphNode.NORMAL)
          continue;
        
        text.append("<VARIABLE TYPE=\"nature\">\n");
        text.append("\t<NAME>" + XMLNormalize(n.ID) + "</NAME>\n");
        
        if(n.outcomes!=null)
          for(int outidx=0; outidx<n.outcomes.length; outidx++)
            text.append("\t<OUTCOME>" + XMLNormalize(n.outcomes[outidx])+
            "</OUTCOME>\n");
        else
          text.append("\t<OUTCOME>true</OUTCOME>\n");
        
        text.append("\t<PROPERTY>position = ("+n.x+","+n.y+")</PROPERTY>\n");
        text.append("</VARIABLE>\n");
      }
      
      //Writing all the nodes definitions and their probability tables
      //If probability table is null then simply write 1 for all
      //the posible outcomes of the parents
      for (int nodeidx=0; nodeidx<nodes.size(); nodeidx++) {
        GraphNode n = (GraphNode) nodes.elementAt(nodeidx);
        if(n.nodeType!=GraphNode.NORMAL)
          continue;
        
        text.append("<DEFINITION>\n");
        text.append("<FOR>" + XMLNormalize(n.ID) + "</FOR>\n");
        int parntOutcomes = 1;
        if(n.prnts!=null)
          for(int pidx=0; pidx<n.prnts.length; pidx++) {
            GraphNode prnt = (GraphNode)nodes.elementAt(n.prnts[pidx]);
            text.append("\t<GIVEN>" + XMLNormalize(prnt.ID) + "</GIVEN>\n");
            if(prnt.outcomes!=null)
              parntOutcomes *= prnt.outcomes.length;
          }
        
        text.append("<TABLE>\n");
        for(int i=0; i<parntOutcomes; i++) {
          if(n.outcomes!=null)
            for(int outidx=0; outidx<n.outcomes.length; outidx++){
              text.append(n.probs[i][outidx]+" ");
            }
          else
            text.append("1");
          text.append('\n');
        }
        text.append("</TABLE>\n");
        text.append("</DEFINITION>\n");
      }
      
      text.append("</NETWORK>\n");
      text.append("</BIF>\n");
      
      outfile.write(text.toString());
      outfile.close();
    }
    catch(IOException ex) { ex.printStackTrace(); }
  } // writeXMLBIF
  
  /** XMLNormalize converts the five standard XML entities in a string
   * g.e. the string V&D's is returned as V&amp;D&apos;s
   * @author Remco Bouckaert (rrb@xm.co.nz)
   * @param sStr string to normalize
   * @return normalized string
   */
  private static String XMLNormalize(String sStr) {
    StringBuffer sStr2 = new StringBuffer();
    for (int iStr = 0; iStr < sStr.length(); iStr++) {
      char c = sStr.charAt(iStr);
      switch (c) {
        case '&': sStr2.append("&amp;"); break;
        case '\'': sStr2.append("&apos;"); break;
        case '\"': sStr2.append("&quot;"); break;
        case '<': sStr2.append("&lt;"); break;
        case '>': sStr2.append("&gt;"); break;
        default:
          sStr2.append(c);
      }
    }
    return sStr2.toString();
  } // XMLNormalize
  
  
} // BIFParser
