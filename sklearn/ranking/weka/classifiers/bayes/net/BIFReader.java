/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/*
 * BIFReader.java
 * Copyright (C) 2003 University of Waikato, Hamilton, New Zealand
 * 
 */

package weka.classifiers.bayes.net;

import weka.classifiers.bayes.BayesNet;
import weka.classifiers.bayes.net.estimate.DiscreteEstimatorBayes;
import weka.core.FastVector;
import weka.core.Instances;
import weka.core.RevisionUtils;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformation.Type;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformationHandler;
import weka.estimators.Estimator;

import java.io.File;
import java.io.StringReader;
import java.util.StringTokenizer;

import javax.xml.parsers.DocumentBuilderFactory;

import org.w3c.dom.CharacterData;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

/**
 <!-- globalinfo-start -->
 * Builds a description of a Bayes Net classifier stored in XML BIF 0.3 format.<br/>
 * <br/>
 * For more details on XML BIF see:<br/>
 * <br/>
 * Fabio Cozman, Marek Druzdzel, Daniel Garcia (1998). XML BIF version 0.3. URL http://www-2.cs.cmu.edu/~fgcozman/Research/InterchangeFormat/.
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;misc{Cozman1998,
 *    author = {Fabio Cozman and Marek Druzdzel and Daniel Garcia},
 *    title = {XML BIF version 0.3},
 *    year = {1998},
 *    URL = {http://www-2.cs.cmu.edu/\~fgcozman/Research/InterchangeFormat/}
 * }
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -D
 *  Do not use ADTree data structure
 * </pre>
 * 
 * <pre> -B &lt;BIF file&gt;
 *  BIF file to compare with
 * </pre>
 * 
 * <pre> -Q weka.classifiers.bayes.net.search.SearchAlgorithm
 *  Search algorithm
 * </pre>
 * 
 * <pre> -E weka.classifiers.bayes.net.estimate.SimpleEstimator
 *  Estimator algorithm
 * </pre>
 * 
 <!-- options-end -->
 *
 * @author Remco Bouckaert (rrb@xm.co.nz)
 * @version $Revision: 1.15 $
 */
public class BIFReader 
    extends BayesNet
    implements TechnicalInformationHandler {
  
    protected int [] m_nPositionX;
    protected int [] m_nPositionY;
    private int [] m_order;
    
    /** for serialization */
    static final long serialVersionUID = -8358864680379881429L;

    /**
     * This will return a string describing the classifier.
     * @return The string.
     */
    public String globalInfo() {
        return 
            "Builds a description of a Bayes Net classifier stored in XML "
        + "BIF 0.3 format.\n\n"
        + "For more details on XML BIF see:\n\n"
        + getTechnicalInformation().toString();
    }

	/** processFile reads a BIFXML file and initializes a Bayes Net
	 * @param sFile name of the file to parse
	 * @return the BIFReader
	 * @throws Exception if processing fails
	 */
	public BIFReader processFile(String sFile) throws Exception {
		m_sFile = sFile;
        DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
        factory.setValidating(true);
        Document doc = factory.newDocumentBuilder().parse(new File(sFile));
        doc.normalize();

        buildInstances(doc, sFile);
        buildStructure(doc);
        return this;
	} // processFile

	public BIFReader processString(String sStr) throws Exception {
        DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
        factory.setValidating(true);
		Document doc = factory.newDocumentBuilder().parse(new org.xml.sax.InputSource(new StringReader(sStr)));
        doc.normalize();
        buildInstances(doc, "from-string");
        buildStructure(doc);
        return this;
	} // processString
	
	
	/** the current filename */
	String m_sFile;
	
	/**
	 * returns the current filename
	 * 
	 * @return the current filename
	 */
	public String getFileName() {
	  return m_sFile;
	}
	
	
	/**
	 * Returns an instance of a TechnicalInformation object, containing 
	 * detailed information about the technical background of this class,
	 * e.g., paper reference or book this class is based on.
	 * 
	 * @return the technical information about this class
	 */
	public TechnicalInformation getTechnicalInformation() {
	  TechnicalInformation 	result;
	  
	  result = new TechnicalInformation(Type.MISC);
	  result.setValue(Field.AUTHOR, "Fabio Cozman and Marek Druzdzel and Daniel Garcia");
	  result.setValue(Field.YEAR, "1998");
	  result.setValue(Field.TITLE, "XML BIF version 0.3");
	  result.setValue(Field.URL, "http://www-2.cs.cmu.edu/~fgcozman/Research/InterchangeFormat/");
	  
	  return result;
	}
	
	/** buildStructure parses the BIF document in the DOM tree contained
	 * in the doc parameter and specifies the the network structure and 
	 * probability tables.
	 * It assumes that buildInstances has been called before
	 * @param doc DOM document containing BIF document in DOM tree
	 * @throws Exception if building of structure fails
	 */
    void buildStructure(Document doc)  throws Exception {
        // Get the name of the network
		// initialize conditional distribution tables
		m_Distributions = new Estimator[m_Instances.numAttributes()][];
        for (int iNode = 0; iNode < m_Instances.numAttributes(); iNode++) {
        	// find definition that goes with this node
        	String sName = m_Instances.attribute(iNode).name();
			Element definition = getDefinition(doc, sName);
/*
	        if (nodelist.getLength() == 0) {
	        	throw new Exception("No definition found for node " + sName);
	        }
	        if (nodelist.getLength() > 1) {
	        	System.err.println("More than one definition found for node " + sName + ". Using first definition.");
	        }
	        Element definition = (Element) nodelist.item(0);
*/	        
	        
	        // get the parents for this node
	        // resolve structure
	        FastVector nodelist = getParentNodes(definition);
	        for (int iParent = 0; iParent < nodelist.size(); iParent++) {
	        	Node parentName = ((Node) nodelist.elementAt(iParent)).getFirstChild();
	        	String sParentName = ((CharacterData) (parentName)).getData();
	        	int nParent = getNode(sParentName);
	        	m_ParentSets[iNode].addParent(nParent, m_Instances);
	        }
	        // resolve conditional probability table
		        int nCardinality = m_ParentSets[iNode].getCardinalityOfParents();
	        int nValues = m_Instances.attribute(iNode).numValues();
	        m_Distributions[iNode] = new Estimator[nCardinality];
			for (int i = 0; i < nCardinality; i++) {
				m_Distributions[iNode][i] = new DiscreteEstimatorBayes(nValues, 0.0f);
			}

/*
	        StringBuffer sTable = new StringBuffer();
	        for (int iText = 0; iText < nodelist.getLength(); iText++) {
	        	sTable.append(((CharacterData) (nodelist.item(iText))).getData());
	        	sTable.append(' ');
	        }
	        StringTokenizer st = new StringTokenizer(sTable.toString());
*/
	        String sTable = getTable(definition);
			StringTokenizer st = new StringTokenizer(sTable.toString());
	        
	        
			for (int i = 0; i < nCardinality; i++) {
				DiscreteEstimatorBayes d = (DiscreteEstimatorBayes) m_Distributions[iNode][i];
				for (int iValue = 0; iValue < nValues; iValue++) {
					String sWeight = st.nextToken();
					d.addValue(iValue, new Double(sWeight).doubleValue());
				}
			}
         }
    } // buildStructure

    /** synchronizes the node ordering of this Bayes network with
     * those in the other network (if possible).
     * @param other Bayes network to synchronize with
     * @throws Exception if nr of attributes differs or not all of the variables have the same name.
     */
    public void Sync(BayesNet other) throws Exception {
    	int nAtts = m_Instances.numAttributes();
    	if (nAtts != other.m_Instances.numAttributes()) {
    		throw new Exception ("Cannot synchronize networks: different number of attributes.");
    	}
        m_order = new int[nAtts];
        for (int iNode = 0; iNode < nAtts; iNode++) {
        	String sName = other.getNodeName(iNode);
        	m_order[getNode(sName)] = iNode;
        }
    } // Sync


    /**
     * Returns all TEXT children of the given node in one string. Between
     * the node values new lines are inserted.
     * 
     * @param node the node to return the content for
     * @return the content of the node
     */
    public String getContent(Element node) {
      NodeList       list;
      Node           item;
      int            i;
      String         result;
      
      result = "";
      list   = node.getChildNodes();
      
      for (i = 0; i < list.getLength(); i++) {
         item = list.item(i);
         if (item.getNodeType() == Node.TEXT_NODE)
            result += "\n" + item.getNodeValue();
      }
         
      return result;
    }


	/** buildInstances parses the BIF document and creates a Bayes Net with its 
	 * nodes specified, but leaves the network structure and probability tables empty.
	 * @param doc DOM document containing BIF document in DOM tree
	 * @param sName default name to give to the Bayes Net. Will be overridden if specified in the BIF document.
	 * @throws Exception if building fails
	 */
	void buildInstances(Document doc, String sName) throws Exception {
		NodeList nodelist;
        // Get the name of the network
        nodelist = selectAllNames(doc);
        if (nodelist.getLength() > 0) {
        	sName = ((CharacterData) (nodelist.item(0).getFirstChild())).getData();
        }

        // Process variables
        nodelist = selectAllVariables(doc);
		int nNodes = nodelist.getLength();
		// initialize structure
		FastVector attInfo = new FastVector(nNodes);

        // Initialize
        m_nPositionX = new int[nodelist.getLength()];
        m_nPositionY = new int[nodelist.getLength()];

        // Process variables
        for (int iNode = 0; iNode < nodelist.getLength(); iNode++) {
            // Get element
			FastVector valueslist;
	        // Get the name of the network
    	    valueslist = selectOutCome(nodelist.item(iNode));

			int nValues = valueslist.size();
			// generate value strings
	        FastVector nomStrings = new FastVector(nValues + 1);
	        for (int iValue = 0; iValue < nValues; iValue++) {
	        	Node node = ((Node) valueslist.elementAt(iValue)).getFirstChild();
	        	String sValue = ((CharacterData) (node)).getData();
	        	if (sValue == null) {
	        		sValue = "Value" + (iValue + 1);
	        	}
				nomStrings.addElement(sValue);
	        }
			FastVector nodelist2;
	        // Get the name of the network
    	    nodelist2 = selectName(nodelist.item(iNode));
    	    if (nodelist2.size() == 0) {
    	    	throw new Exception ("No name specified for variable");
    	    }
    	    String sNodeName = ((CharacterData) (((Node) nodelist2.elementAt(0)).getFirstChild())).getData();

			weka.core.Attribute att = new weka.core.Attribute(sNodeName, nomStrings);
			attInfo.addElement(att);

    	    valueslist = selectProperty(nodelist.item(iNode));
			nValues = valueslist.size();
			// generate value strings
	        for (int iValue = 0; iValue < nValues; iValue++) {
                // parsing for strings of the form "position = (73, 165)"
	        	Node node = ((Node)valueslist.elementAt(iValue)).getFirstChild();
	        	String sValue = ((CharacterData) (node)).getData();
                if (sValue.startsWith("position")) {
                    int i0 = sValue.indexOf('(');
                    int i1 = sValue.indexOf(',');
                    int i2 = sValue.indexOf(')');
                    String sX = sValue.substring(i0 + 1, i1).trim();
                    String sY = sValue.substring(i1 + 1, i2).trim();
                    try {
                    	m_nPositionX[iNode] = (int) Integer.parseInt(sX);
                    	m_nPositionY[iNode] = (int) Integer.parseInt(sY);
                    } catch (NumberFormatException e) {
                    	System.err.println("Wrong number format in position :(" + sX + "," + sY +")");
                   	    m_nPositionX[iNode] = 0;
                   	    m_nPositionY[iNode] = 0;
                    }
                }
            }

        }
        
 		m_Instances = new Instances(sName, attInfo, 100);
		m_Instances.setClassIndex(nNodes - 1);
		setUseADTree(false);
		initStructure();
	} // buildInstances

//	/** selectNodeList selects list of nodes from document specified in XPath expression
//	 * @param doc : document (or node) to query
//	 * @param sXPath : XPath expression
//	 * @return list of nodes conforming to XPath expression in doc
//	 * @throws Exception
//	 */
//	private NodeList selectNodeList(Node doc, String sXPath) throws Exception {
//		NodeList nodelist = org.apache.xpath.XPathAPI.selectNodeList(doc, sXPath);
//		return nodelist;
//	} // selectNodeList

	NodeList selectAllNames(Document doc) throws Exception {
		//NodeList nodelist = selectNodeList(doc, "//NAME");
		NodeList nodelist = doc.getElementsByTagName("NAME");
		return nodelist;
	} // selectAllNames

	NodeList selectAllVariables(Document doc) throws Exception {
		//NodeList nodelist = selectNodeList(doc, "//VARIABLE");
		NodeList nodelist = doc.getElementsByTagName("VARIABLE");
		return nodelist;
	} // selectAllVariables

	Element getDefinition(Document doc, String sName) throws Exception {
		//NodeList nodelist = selectNodeList(doc, "//DEFINITION[normalize-space(FOR/text())=\"" + sName + "\"]");

		NodeList nodelist = doc.getElementsByTagName("DEFINITION");
		for (int iNode = 0; iNode < nodelist.getLength(); iNode++) {
			Node node = nodelist.item(iNode);
			FastVector list = selectElements(node, "FOR");
			if (list.size() > 0) {
				Node forNode = (Node) list.elementAt(0);
				if (getContent((Element) forNode).trim().equals(sName)) {
					return (Element) node;
				}
			}
		}
		throw new Exception("Could not find definition for ((" + sName + "))");
	} // getDefinition

	FastVector getParentNodes(Node definition) throws Exception {
		//NodeList nodelist = selectNodeList(definition, "GIVEN");
		FastVector nodelist = selectElements(definition, "GIVEN");
		return nodelist;
	} // getParentNodes

	String getTable(Node definition) throws Exception {
		//NodeList nodelist = selectNodeList(definition, "TABLE/text()");
		FastVector nodelist = selectElements(definition, "TABLE");
		String sTable = getContent((Element) nodelist.elementAt(0));
		sTable = sTable.replaceAll("\\n"," ");
		return sTable;
	} // getTable

	FastVector selectOutCome(Node item) throws Exception {
		//NodeList nodelist = selectNodeList(item, "OUTCOME");
		FastVector nodelist = selectElements(item, "OUTCOME");
		return nodelist;
	} // selectOutCome

	FastVector selectName(Node item) throws Exception {
	   //NodeList nodelist = selectNodeList(item, "NAME");
	   FastVector nodelist = selectElements(item, "NAME");
	   return nodelist;
   } // selectName

   FastVector selectProperty(Node item) throws Exception {
	  // NodeList nodelist = selectNodeList(item, "PROPERTY");
	  FastVector nodelist = selectElements(item, "PROPERTY");
	  return nodelist;
   } // selectProperty

	FastVector selectElements(Node item, String sElement) throws Exception {
	  NodeList children = item.getChildNodes();
	  FastVector nodelist = new FastVector();
	  for (int iNode = 0; iNode < children.getLength(); iNode++) {
		Node node = children.item(iNode);
		if ((node.getNodeType() == Node.ELEMENT_NODE) && node.getNodeName().equals(sElement)) {
			nodelist.addElement(node);
		}
	  }
	  return nodelist;
  } // selectElements
	/** Count nr of arcs missing from other network compared to current network
	 * Note that an arc is not 'missing' if it is reversed.
	 * @param other network to compare with
	 * @return nr of missing arcs
	 */
	public int missingArcs(BayesNet other) {
		try {
			Sync(other);
			int nMissing = 0;
			for (int iAttribute = 0; iAttribute < m_Instances.numAttributes(); iAttribute++) {
				for (int iParent = 0; iParent < m_ParentSets[iAttribute].getNrOfParents(); iParent++) {
					int nParent = m_ParentSets[iAttribute].getParent(iParent);
					if (!other.getParentSet(m_order[iAttribute]).contains(m_order[nParent]) && !other.getParentSet(m_order[nParent]).contains(m_order[iAttribute])) {
						nMissing++;
					}
				}
			}
			return nMissing;
		} catch (Exception e) {
			System.err.println(e.getMessage());
			return 0;
		}
	} // missingArcs

	/** Count nr of exta arcs  from other network compared to current network
	 * Note that an arc is not 'extra' if it is reversed.
	 * @param other network to compare with
	 * @return nr of missing arcs
	 */
	public int extraArcs(BayesNet other) {
		try {
			Sync(other);
			int nExtra = 0;
			for (int iAttribute = 0; iAttribute < m_Instances.numAttributes(); iAttribute++) {
				for (int iParent = 0; iParent < other.getParentSet(m_order[iAttribute]).getNrOfParents(); iParent++) {
					int nParent = m_order[other.getParentSet(m_order[iAttribute]).getParent(iParent)];
					if (!m_ParentSets[iAttribute].contains(nParent) && !m_ParentSets[nParent].contains(iAttribute)) {
						nExtra++;
					}
				}
			}
			return nExtra;
		} catch (Exception e) {
			System.err.println(e.getMessage());
			return 0;
		}
	} // extraArcs


	/** calculates the divergence between the probability distribution
	 * represented by this network and that of another, that is,
	 * \sum_{x\in X} P(x)log P(x)/Q(x)
	 * where X is the set of values the nodes in the network can take,
	 * P(x) the probability of this network for configuration x
	 * Q(x) the probability of the other network for configuration x
	 * @param other network to compare with
	 * @return divergence between this and other Bayes Network
	 */
	public double divergence(BayesNet other) {
		try {
			Sync(other);
			// D: divergence
			double D = 0.0;
			int nNodes = m_Instances.numAttributes();
			int [] nCard = new int[nNodes];
			for (int iNode = 0; iNode < nNodes; iNode++) {
				nCard[iNode] = m_Instances.attribute(iNode).numValues();
			}
			// x: holds current configuration of nodes
			int [] x = new int[nNodes];
			// simply sum over all configurations to calc divergence D
			int i = 0;
			while (i < nNodes) {
				// update configuration
				x[i]++;
				while (i < nNodes && x[i] == m_Instances.attribute(i).numValues()) {
					x[i] = 0;
					i++;
					if (i < nNodes){
						x[i]++;
					}
				}
				if (i < nNodes) {
					i = 0;
					// calc P(x) and Q(x)
					double P = 1.0;
					for (int iNode = 0; iNode < nNodes; iNode++) {
						int iCPT = 0;
						for (int iParent = 0; iParent < m_ParentSets[iNode].getNrOfParents(); iParent++) {
					    	int nParent = m_ParentSets[iNode].getParent(iParent);
						    iCPT = iCPT * nCard[nParent] + x[nParent];
						} 
						P = P * m_Distributions[iNode][iCPT].getProbability(x[iNode]);
					}
	
					double Q = 1.0;
					for (int iNode = 0; iNode < nNodes; iNode++) {
						int iCPT = 0;
						for (int iParent = 0; iParent < other.getParentSet(m_order[iNode]).getNrOfParents(); iParent++) {
					    	int nParent = m_order[other.getParentSet(m_order[iNode]).getParent(iParent)];
						    iCPT = iCPT * nCard[nParent] + x[nParent];
						} 
						Q = Q * other.m_Distributions[m_order[iNode]][iCPT].getProbability(x[iNode]);
					}
	
					// update divergence if probabilities are positive
					if (P > 0.0 && Q > 0.0) {
						D = D + P * Math.log(Q / P);
					}
				}
			}
			return D;
		} catch (Exception e) {
			System.err.println(e.getMessage());
			return 0;
		}
	} // divergence

	/** Count nr of reversed arcs from other network compared to current network
	 * @param other network to compare with
	 * @return nr of missing arcs
	 */
	public int reversedArcs(BayesNet other) {
		try {
			Sync(other);
			int nReversed = 0;
		    for (int iAttribute = 0; iAttribute < m_Instances.numAttributes(); iAttribute++) {
				for (int iParent = 0; iParent < m_ParentSets[iAttribute].getNrOfParents(); iParent++) {
					int nParent = m_ParentSets[iAttribute].getParent(iParent);
					if (!other.getParentSet(m_order[iAttribute]).contains(m_order[nParent]) && other.getParentSet(m_order[nParent]).contains(m_order[iAttribute])) {
						nReversed++;
					}
				}
			}
			return nReversed;
		} catch (Exception e) {
			System.err.println(e.getMessage());
			return 0;
		}
	} // reversedArcs
	/** getNode finds the index of the node with name sNodeName
	 * and throws an exception if no such node can be found.
	 * @param sNodeName name of the node to get the index from
	 * @return index of the node with name sNodeName
	 * @throws Exception if node cannot be found
	 */
    public int getNode(String sNodeName) throws Exception {
		int iNode = 0;
		while (iNode < m_Instances.numAttributes()) {
			if (m_Instances.attribute(iNode).name().equals(sNodeName)) {
				return iNode;
			}
			iNode++;
		}
   		throw new Exception("Could not find node [[" + sNodeName + "]]");
    } // getNode

    /**
     * the default constructor
     */
    public BIFReader() {
    }
    
    /**
     * Returns the revision string.
     * 
     * @return		the revision
     */
    public String getRevision() {
      return RevisionUtils.extract("$Revision: 1.15 $");
    }

    /**
     * Loads the file specified as first parameter and prints it to stdout.
     * 
     * @param args the command line parameters
     */
    public static void main(String[] args) {
        try {
            BIFReader br = new BIFReader();
            br.processFile(args[0]);
	    System.out.println(br.toString());
        
        }
        catch (Throwable t) {
            t.printStackTrace();
        }
    } // main
} // class BIFReader

