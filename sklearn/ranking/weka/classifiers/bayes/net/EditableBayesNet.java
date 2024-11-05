package weka.classifiers.bayes.net;
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
 * EditableBayesNet.java
 *
 */

import java.io.Serializable;
import java.io.StringReader;
import java.util.StringTokenizer;

import javax.xml.parsers.DocumentBuilderFactory;

import org.w3c.dom.CharacterData;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import weka.classifiers.bayes.BayesNet;
import weka.classifiers.bayes.net.estimate.DiscreteEstimatorBayes;
import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instances;
import weka.core.RevisionUtils;
import weka.core.SerializedObject;
import weka.estimators.Estimator;
import weka.filters.Filter;
import weka.filters.unsupervised.attribute.Reorder;


/**
 <!-- globalinfo-start -->
 * Bayes Network learning using various search algorithms and quality measures.<br/>
 * Base class for a Bayes Network classifier. Provides datastructures (network structure, conditional probability distributions, etc.) and facilities common to Bayes Network learning algorithms like K2 and B.<br/>
 * <br/>
 * For more information see:<br/>
 * <br/>
 * http://www.cs.waikato.ac.nz/~remco/weka.pdf
 * <p/>
 <!-- globalinfo-end -->
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
 * @version $Revision: 4899 $
 */

public class EditableBayesNet extends BayesNet {
	/** for serialization */
	static final long serialVersionUID = 746037443258735954L;

	/** location of nodes, used for graph drawing * */
	protected FastVector m_nPositionX;

	protected FastVector m_nPositionY;

	/** marginal distributions * */
	protected FastVector m_fMarginP;

	/** evidence values, used for evidence propagation * */
	protected FastVector m_nEvidence;

	/** standard constructor * */
	public EditableBayesNet() {
		super();
		m_nEvidence = new FastVector(0);
		m_fMarginP = new FastVector(0);
		m_nPositionX = new FastVector();
		m_nPositionY = new FastVector();
		clearUndoStack();
	} // c'tor

	/** constructor, creates empty network with nodes based on the attributes in a data set */
	public EditableBayesNet(Instances instances) {
		try {
			if (instances.classIndex() < 0) {
				instances.setClassIndex(instances.numAttributes() - 1);
			}
			m_Instances = normalizeDataSet(instances);
		} catch (Exception e) {
			e.printStackTrace();
		}

		int nNodes = getNrOfNodes();
		m_ParentSets = new ParentSet[nNodes];
		for (int i = 0; i < nNodes; i++) {
			m_ParentSets[i] = new ParentSet();
		}
		m_Distributions = new Estimator[nNodes][];
		for (int iNode = 0; iNode < nNodes; iNode++) {
			m_Distributions[iNode] = new Estimator[1];
			m_Distributions[iNode][0] = new DiscreteEstimatorBayes(getCardinality(iNode), 0.5);
		}

		m_nEvidence = new FastVector(nNodes);
		for (int i = 0; i < nNodes; i++) {
			m_nEvidence.addElement(-1);
		}
		m_fMarginP = new FastVector(nNodes);
		for (int i = 0; i < nNodes; i++) {
			double[] P = new double[getCardinality(i)];
			m_fMarginP.addElement(P);
		}

		m_nPositionX = new FastVector(nNodes);
		m_nPositionY = new FastVector(nNodes);
		for (int iNode = 0; iNode < nNodes; iNode++) {
			m_nPositionX.addElement(iNode%10 * 50);
			m_nPositionY.addElement(((int)(iNode/10)) * 50);
		}

	} // c'tor

	/** constructor, copies Bayesian network structure from a Bayesian network
	 * encapsulated in a BIFReader
	 */
	public EditableBayesNet(BIFReader other) {
		m_Instances = other.m_Instances;
		m_ParentSets = other.getParentSets();
		m_Distributions = other.getDistributions();

		int nNodes = getNrOfNodes();
		m_nPositionX = new FastVector(nNodes);
		m_nPositionY = new FastVector(nNodes);
		for (int i = 0; i < nNodes; i++) {
			m_nPositionX.addElement(other.m_nPositionX[i]);
			m_nPositionY.addElement(other.m_nPositionY[i]);
		}
		m_nEvidence = new FastVector(nNodes);
		for (int i = 0; i < nNodes; i++) {
			m_nEvidence.addElement(-1);
		}
		m_fMarginP = new FastVector(nNodes);
		for (int i = 0; i < nNodes; i++) {
			double[] P = new double[getCardinality(i)];
			m_fMarginP.addElement(P);
		}
		clearUndoStack();
	} // c'tor

	/**
	 * constructor that potentially initializes instances as well
	 *
	 * @param bSetInstances
	 *            flag indicating whether to initialize instances or not
	 */
	public EditableBayesNet(boolean bSetInstances) {
		super();
		m_nEvidence = new FastVector(0);
		m_fMarginP = new FastVector(0);
		m_nPositionX = new FastVector();
		m_nPositionY = new FastVector();
		clearUndoStack();
		if (bSetInstances) {
			m_Instances = new Instances("New Network", new FastVector(0), 0);
		}
	} // c'tor


	/** Assuming a network structure is defined and we want to learn from data,
	 * the data set must be put if correct order first and possibly discretized/missing
	 * values filled in before proceeding to CPT learning.
	 * @param instances data set to learn from
	 * @exception Exception when data sets are not compatible, e.g., a variable is missing
	 * or a variable has different nr of values.
	 */
	public void setData(Instances instances) throws Exception {
		// sync order of variables
		int [] order = new int [getNrOfNodes()];
		for (int iNode = 0; iNode < getNrOfNodes(); iNode++) {
			String sName = getNodeName(iNode);
			int nNode = 0;
			while (nNode < getNrOfNodes() && !sName.equals(instances.attribute(nNode).name())) {
				nNode++;
			}
			if (nNode >= getNrOfNodes()) {
				throw new Exception("Cannot find node named [[[" + sName + "]]] in the data");
			}
			order[iNode] = nNode;
		}
		Reorder reorderFilter = new Reorder();
		reorderFilter.setAttributeIndicesArray(order);
		reorderFilter.setInputFormat(instances);
		instances = Filter.useFilter(instances, reorderFilter);
		// filter using discretization/missing values filter
		Instances newInstances = new Instances(m_Instances, 0);
		if (m_DiscretizeFilter == null && m_MissingValuesFilter == null) {
			newInstances = normalizeDataSet(instances);
		} else {
			for (int iInstance = 0; iInstance < instances.numInstances(); iInstance++) {
				newInstances.add(normalizeInstance(instances.instance(iInstance)));
			}
		}
		//sanity check
		for (int iNode = 0; iNode < getNrOfNodes(); iNode++) {
			if (newInstances.attribute(iNode).numValues() != getCardinality(iNode)) {
				throw new Exception("Number of values of node [[[" + getNodeName(iNode) + "]]] differs in (discretized) dataset." );
			}
		}
		// if we got this far, all is ok with the data set and
		// we can replace data set of Bayes net
		m_Instances = newInstances;
	} // setData

	/** returns index of node with given name, or -1 if no such node exists
	 * @param sNodeName name of the node to get index for
	 */
	public int getNode2(String sNodeName) {
		int iNode = 0;
		while (iNode < m_Instances.numAttributes()) {
			if (m_Instances.attribute(iNode).name().equals(sNodeName)) {
				return iNode;
			}
			iNode++;
		}
		return -1;
	} // getNode2

	/** returns index of node with given name. Throws exception if no such node exists
	 * @param sNodeName name of the node to get index for
	 */
	public int getNode(String sNodeName) throws Exception {
		int iNode = getNode2(sNodeName);
		if (iNode < 0) {
			throw new Exception("Could not find node [[" + sNodeName + "]]");
		}
		return iNode;
	} // getNode

	/**
	 * Add new node to the network, initializing instances, parentsets,
	 * distributions. Used for manual manipulation of the Bayesian network.
	 *
	 * @param sName
	 *            name of the node. If the name already exists, an x is appended
	 *            to the name
	 * @param nCardinality
	 *            number of values for this node
	 * @throws Exception
	 */
	public void addNode(String sName, int nCardinality) throws Exception {
		addNode(sName, nCardinality, 100 + getNrOfNodes() * 10, 100 + getNrOfNodes() * 10);
	} // addNode

	/** Add node to network at a given position, initializing instances, parentsets,
	 * distributions. Used for manual manipulation of the Bayesian network.
	 *
	 * @param sName
	 *            name of the node. If the name already exists, an x is appended
	 *            to the name
	 * @param nCardinality
	 *            number of values for this node
	 * @param nPosX x-coordiate of the position to place this node
	 * @param nPosY y-coordiate of the position to place this node
	 * @throws Exception
	 */
	public void addNode(String sName, int nCardinality, int nPosX, int nPosY) throws Exception {
		if (getNode2(sName) >= 0) {
			addNode(sName + "x", nCardinality);
			return ;
		}
		// update instances
		FastVector values = new FastVector(nCardinality);
		for (int iValue = 0; iValue < nCardinality; iValue++) {
			values.addElement("Value" + (iValue + 1));
		}
		Attribute att = new Attribute(sName, values);
		m_Instances.insertAttributeAt(att, m_Instances.numAttributes());
		int nAtts = m_Instances.numAttributes();
		// update parentsets
		ParentSet[] parentSets = new ParentSet[nAtts];
		for (int iParentSet = 0; iParentSet < nAtts - 1; iParentSet++) {
			parentSets[iParentSet] = m_ParentSets[iParentSet];
		}
		parentSets[nAtts - 1] = new ParentSet();
		m_ParentSets = parentSets;
		// update distributions
		Estimator[][] distributions = new Estimator[nAtts][];
		for (int iNode = 0; iNode < nAtts - 1; iNode++) {
			distributions[iNode] = m_Distributions[iNode];
		}
		distributions[nAtts - 1] = new Estimator[1];
		distributions[nAtts - 1][0] = new DiscreteEstimatorBayes(nCardinality, 0.5);
		m_Distributions = distributions;
		// update positions
		m_nPositionX.addElement(nPosX);
		m_nPositionY.addElement(nPosY);
		// update evidence & margins
		m_nEvidence.addElement(-1);
		double[] fMarginP = new double[nCardinality];
		for (int iValue = 0; iValue < nCardinality; iValue++) {
			fMarginP[iValue] = 1.0 / nCardinality;
		}
		m_fMarginP.addElement(fMarginP);
		// update undo stack
		if (m_bNeedsUndoAction) {
			addUndoAction(new AddNodeAction(sName, nCardinality, nPosX, nPosY));
		}
	} // addNode

	/**
	 * Delete node from the network, updating instances, parentsets,
	 * distributions Conditional distributions are condensed by taking the
	 * values for the target node to be its first value. Used for manual
	 * manipulation of the Bayesian network.
	 *
	 * @param sName
	 *            name of the node. If the name does not exists an exception is
	 *            thrown
	 * @throws Exception
	 */
	public void deleteNode(String sName) throws Exception {
		int nTargetNode = getNode(sName);
		deleteNode(nTargetNode);
	} // deleteNode

	/**
	 * Delete node from the network, updating instances, parentsets,
	 * distributions Conditional distributions are condensed by taking the
	 * values for the target node to be its first value. Used for manual
	 * manipulation of the Bayesian network.
	 *
	 * @param nTargetNode
	 *            index of the node to delete.
	 * @throws Exception
	 */
	public void deleteNode(int nTargetNode) throws Exception {
		// update undo stack
		if (m_bNeedsUndoAction) {
			addUndoAction(new DeleteNodeAction(nTargetNode));
		}
		int nAtts = m_Instances.numAttributes() - 1;
		int nTargetCard = m_Instances.attribute(nTargetNode).numValues();
		// update distributions
		Estimator[][] distributions = new Estimator[nAtts][];
		for (int iNode = 0; iNode < nAtts; iNode++) {
			int iNode2 = iNode;
			if (iNode >= nTargetNode) {
				iNode2++;
			}
			Estimator[] distribution = m_Distributions[iNode2];
			if (m_ParentSets[iNode2].contains(nTargetNode)) {
				// condense distribution, use values for targetnode = 0
				int nParentCard = m_ParentSets[iNode2].getCardinalityOfParents();
				nParentCard = nParentCard / nTargetCard;
				Estimator[] distribution2 = new Estimator[nParentCard];
				for (int iParent = 0; iParent < nParentCard; iParent++) {
					distribution2[iParent] = distribution[iParent];
				}
				distribution = distribution2;
			}
			distributions[iNode] = distribution;
		}
		m_Distributions = distributions;
		// update parentsets
		ParentSet[] parentSets = new ParentSet[nAtts];
		for (int iParentSet = 0; iParentSet < nAtts; iParentSet++) {
			int iParentSet2 = iParentSet;
			if (iParentSet >= nTargetNode) {
				iParentSet2++;
			}
			ParentSet parentset = m_ParentSets[iParentSet2];
			parentset.deleteParent(nTargetNode, m_Instances);
			for (int iParent = 0; iParent < parentset.getNrOfParents(); iParent++) {
				int nParent = parentset.getParent(iParent);
				if (nParent > nTargetNode) {
					parentset.SetParent(iParent, nParent - 1);
				}
			}
			parentSets[iParentSet] = parentset;
		}
		m_ParentSets = parentSets;
		// update instances
		m_Instances.setClassIndex(-1);
		m_Instances.deleteAttributeAt(nTargetNode);
		m_Instances.setClassIndex(nAtts - 1);

		// update positions
		m_nPositionX.removeElementAt(nTargetNode);
		m_nPositionY.removeElementAt(nTargetNode);
		// update evidence & margins
		m_nEvidence.removeElementAt(nTargetNode);
		m_fMarginP.removeElementAt(nTargetNode);
	} // deleteNode

	/**
	 * Delete nodes with indexes in selection from the network, updating instances, parentsets,
	 * distributions Conditional distributions are condensed by taking the
	 * values for the target node to be its first value. Used for manual
	 * manipulation of the Bayesian network.
	 *
	 * @param nodes
	 *            array of indexes of nodes to delete.
	 * @throws Exception
	 */
	public void deleteSelection(FastVector nodes) {
		// sort before proceeding
		for (int i = 0; i < nodes.size(); i++) {
			for (int j = i + 1; j < nodes.size(); j++) {
				if ((Integer) nodes.elementAt(i) > (Integer) nodes.elementAt(j)) {
					int h = (Integer) nodes.elementAt(i);
					nodes.setElementAt(nodes.elementAt(j), i);
					nodes.setElementAt(h, j);
				}
			}
		}
		// update undo stack
		if (m_bNeedsUndoAction) {
			addUndoAction(new DeleteSelectionAction(nodes));
		}
		boolean bNeedsUndoAction = m_bNeedsUndoAction;
		m_bNeedsUndoAction = false;
		try {
			for (int iNode = nodes.size() - 1; iNode >= 0; iNode--) {
				deleteNode((Integer) nodes.elementAt(iNode));
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		m_bNeedsUndoAction = bNeedsUndoAction;
	} // deleteSelection

	/** XML helper function for selecting elements under a node with a given name
	 * @param item XMLNode to select items from
	 * @param sElement name of the element to return
	 */
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

	/**
	 * XML helper function. Returns all TEXT children of the given node in one string. Between the
	 * node values new lines are inserted.
	 *
	 * @param node
	 *            the node to return the content for
	 * @return the content of the node
	 */
	public String getContent(Element node) {
		NodeList list;
		Node item;
		int i;
		String result;

		result = "";
		list = node.getChildNodes();

		for (i = 0; i < list.getLength(); i++) {
			item = list.item(i);
			if (item.getNodeType() == Node.TEXT_NODE)
				result += "\n" + item.getNodeValue();
		}

		return result;
	}

	/** XML helper function that returns DEFINITION element from a XMLBIF document
	 * for a node with a given name.
	 * @param doc XMLBIF document
	 * @param sName name of the node to get the definition for
	 */
	Element getDefinition(Document doc, String sName) throws Exception {
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


	/** Paste modes. This allows for verifying that a past action does not cause
	 * any problems before actually performing the paste operation.
	 */
	final static int TEST = 0;
	final static int EXECUTE = 1;

	/** Apply paste operation with XMLBIF fragment. This adds nodes in the XMLBIF fragment
	 * to the network, together with its parents. First, paste in test mode to verify
	 * no problems occur, then execute paste operation. If a problem occurs (e.g. parent
	 * does not exist) then a exception is thrown.
	 * @param sXML XMLBIF fragment to paste into the network
	 */
	public void paste(String sXML) throws Exception {
		try {
			paste(sXML, TEST);
		} catch (Exception e) {
			throw e;
		}
		paste(sXML, EXECUTE);
	} // paste

	/** Apply paste operation with XMLBIF fragment. Depending on the paste mode, the
	 * nodes are actually added to the network or it is just tested that the nodes can
	 * be added to the network.
	 * @param sXML XMLBIF fragment to paste into the network
	 * @param mode paste mode TEST or EXECUTE
	 */
	void paste(String sXML, int mode) throws Exception {
		DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
		factory.setValidating(true);
		Document doc = factory.newDocumentBuilder().parse(new org.xml.sax.InputSource(new StringReader(sXML)));
		doc.normalize();

		// create nodes first
		NodeList nodelist = doc.getElementsByTagName("VARIABLE");
		FastVector sBaseNames = new FastVector();
		Instances instances = new Instances(m_Instances, 0);
		int nBase = instances.numAttributes();
		for (int iNode = 0; iNode < nodelist.getLength(); iNode++) {
			// Get element
			FastVector valueslist;
			// Get the name of the node
			valueslist = selectElements(nodelist.item(iNode), "OUTCOME");

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
			nodelist2 = selectElements(nodelist.item(iNode), "NAME");
			if (nodelist2.size() == 0) {
				throw new Exception("No name specified for variable");
			}
			String sBaseName = ((CharacterData) (((Node) nodelist2.elementAt(0)).getFirstChild())).getData();
			sBaseNames.addElement(sBaseName);
			String sNodeName = sBaseName;
			if (getNode2(sNodeName) >= 0) {
				sNodeName = "Copy of " + sBaseName;
			}
			int iAttempt = 2;
			while (getNode2(sNodeName) >= 0) {
				sNodeName = "Copy (" + iAttempt + ") of " + sBaseName;
				iAttempt++;
			}

			Attribute att = new Attribute(sNodeName, nomStrings);
			instances.insertAttributeAt(att, instances.numAttributes());

			valueslist = selectElements(nodelist.item(iNode), "PROPERTY");
			nValues = valueslist.size();
			// generate value strings
			int nPosX = iAttempt * 10;
			int nPosY = iAttempt * 10;
			for (int iValue = 0; iValue < nValues; iValue++) {
				// parsing for strings of the form "position = (73, 165)"
				Node node = ((Node) valueslist.elementAt(iValue)).getFirstChild();
				String sValue = ((CharacterData) (node)).getData();
				if (sValue.startsWith("position")) {
					int i0 = sValue.indexOf('(');
					int i1 = sValue.indexOf(',');
					int i2 = sValue.indexOf(')');
					String sX = sValue.substring(i0 + 1, i1).trim();
					String sY = sValue.substring(i1 + 1, i2).trim();
					try {
						nPosX = (Integer.parseInt(sX) + iAttempt * 10);
						nPosY = (Integer.parseInt(sY) + iAttempt * 10);
					} catch (NumberFormatException e) {
						System.err.println("Wrong number format in position :(" + sX + "," + sY + ")");
					}
				}
			}
			if (mode == EXECUTE) {
				m_nPositionX.addElement(nPosX);
				m_nPositionY.addElement(nPosY);
			}

		}

		FastVector nodelist2;
		Estimator[][] distributions = new Estimator[nBase + sBaseNames.size()][];
		ParentSet[] parentsets = new ParentSet[nBase + sBaseNames.size()];
		for (int iNode = 0; iNode < nBase; iNode++) {
			distributions[iNode] = m_Distributions[iNode];
			parentsets[iNode] = m_ParentSets[iNode];
		}
		if (mode == EXECUTE) {
			m_Instances = instances;
		}
		// create arrows & create distributions
		for (int iNode = 0; iNode < sBaseNames.size(); iNode++) {
			// find definition that goes with this node
			String sName = (String) sBaseNames.elementAt(iNode);
			Element definition = getDefinition(doc, sName);
			parentsets[nBase + iNode] = new ParentSet();

			// get the parents for this node
			// resolve structure
			nodelist2 = selectElements(definition, "GIVEN");
			for (int iParent = 0; iParent < nodelist2.size(); iParent++) {
				Node parentName = ((Node) nodelist2.elementAt(iParent)).getFirstChild();
				String sParentName = ((CharacterData) (parentName)).getData();
				int nParent = -1;
				for (int iBase = 0; iBase < sBaseNames.size(); iBase++) {
					if (sParentName.equals((String) sBaseNames.elementAt(iBase))) {
						nParent = nBase + iBase;
					}
				}
				if (nParent < 0) {
					nParent = getNode(sParentName);
				}
				parentsets[nBase + iNode].addParent(nParent, instances);
			}
			// resolve conditional probability table
			int nCardinality = parentsets[nBase + iNode].getCardinalityOfParents();
			int nValues = instances.attribute(nBase + iNode).numValues();
			distributions[nBase + iNode] = new Estimator[nCardinality];
			for (int i = 0; i < nCardinality; i++) {
				distributions[nBase + iNode][i] = new DiscreteEstimatorBayes(nValues, 0.0f);
			}

			String sTable = getContent((Element) selectElements(definition, "TABLE").elementAt(0));
			sTable = sTable.replaceAll("\\n", " ");
			StringTokenizer st = new StringTokenizer(sTable.toString());

			for (int i = 0; i < nCardinality; i++) {
				DiscreteEstimatorBayes d = (DiscreteEstimatorBayes) distributions[nBase + iNode][i];
				for (int iValue = 0; iValue < nValues; iValue++) {
					String sWeight = st.nextToken();
					d.addValue(iValue, new Double(sWeight).doubleValue());
				}
			}
			if (mode == EXECUTE) {
				m_nEvidence.insertElementAt(-1, nBase + iNode);
				m_fMarginP.insertElementAt(new double[getCardinality(nBase + iNode)], nBase + iNode);
			}
		}
		if (mode == EXECUTE) {
			m_Distributions = distributions;
			m_ParentSets = parentsets;
		}
		// update undo stack
		if (mode == EXECUTE && m_bNeedsUndoAction) {
			addUndoAction(new PasteAction(sXML, nBase));
		}
	} // paste

	/**
	 * Add arc between two nodes Distributions are updated by duplication for
	 * every value of the parent node.
	 *
	 * @param sParent
	 *            name of the parent node
	 * @param sChild
	 *            name of the child node
	 * @throws Exception
	 *             if parent or child cannot be found in network
	 */
	public void addArc(String sParent, String sChild) throws Exception {
		int nParent = getNode(sParent);
		int nChild = getNode(sChild);
		addArc(nParent, nChild);
	} // addArc

	/**
	 * Add arc between two nodes Distributions are updated by duplication for
	 * every value of the parent node.
	 *
	 * @param nParent
	 *            index of the parent node
	 * @param nChild
	 *            index of the child node
	 * @throws Exception
	 */
	public void addArc(int nParent, int nChild) throws Exception {
		// update undo stack
		if (m_bNeedsUndoAction) {
			addUndoAction(new AddArcAction(nParent, nChild));
		}
		int nOldCard = m_ParentSets[nChild].getCardinalityOfParents();
		// update parentsets
		m_ParentSets[nChild].addParent(nParent, m_Instances);
		// update distributions
		int nNewCard = m_ParentSets[nChild].getCardinalityOfParents();
		Estimator[] ds = new Estimator[nNewCard];
		for (int iParent = 0; iParent < nNewCard; iParent++) {
			ds[iParent] = Estimator.clone(m_Distributions[nChild][iParent % nOldCard]);
		}
		m_Distributions[nChild] = ds;
	} // addArc

	/**
	 * Add arc between parent node and each of the nodes in a given list.
	 * Distributions are updated as above.
	 *
	 * @param sParent
	 *            name of the parent node
	 * @param nodes
	 *            array of indexes of child nodes
	 * @throws Exception
	 */
	public void addArc(String sParent, FastVector nodes) throws Exception {
		int nParent = getNode(sParent);
		// update undo stack
		if (m_bNeedsUndoAction) {
			addUndoAction(new AddArcAction(nParent, nodes));
		}
		boolean bNeedsUndoAction = m_bNeedsUndoAction;
		m_bNeedsUndoAction = false;
		for (int iNode = 0; iNode < nodes.size(); iNode++) {
			int nNode = (Integer) nodes.elementAt(iNode);
			addArc(nParent, nNode);
		}
		m_bNeedsUndoAction = bNeedsUndoAction;
	} // addArc

	/**
	 * Delete arc between two nodes. Distributions are updated by condensing for
	 * the parent node taking its first value.
	 *
	 * @param sParent
	 *            name of the parent node
	 * @param sChild
	 *            name of the child node
	 * @throws Exception
	 *             if parent or child cannot be found in network
	 */
	public void deleteArc(String sParent, String sChild) throws Exception {
		int nParent = getNode(sParent);
		int nChild = getNode(sChild);
		deleteArc(nParent, nChild);
	} // deleteArc

	/**
	 * Delete arc between two nodes. Distributions are updated by condensing for
	 * the parent node taking its first value.
	 *
	 * @param nParent
	 *            index of the parent node
	 * @param nChild
	 *            index of the child node
	 * @throws Exception
	 */
	public void deleteArc(int nParent, int nChild) throws Exception {
		// update undo stack
		if (m_bNeedsUndoAction) {
			addUndoAction(new DeleteArcAction(nParent, nChild));
		}
		// update distributions
		// condense distribution, use values for targetnode = 0
		int nParentCard = m_ParentSets[nChild].getCardinalityOfParents();
		int nTargetCard = m_Instances.attribute(nChild).numValues();
		nParentCard = nParentCard / nTargetCard;
		Estimator[] distribution2 = new Estimator[nParentCard];
		for (int iParent = 0; iParent < nParentCard; iParent++) {
			distribution2[iParent] = m_Distributions[nChild][iParent];
		}
		m_Distributions[nChild] = distribution2;
		// update parentsets
		m_ParentSets[nChild].deleteParent(nParent, m_Instances);
	} // deleteArc


	/** specify distribution of a node
	 * @param sName name of the node to specify distribution for
	 * @param P matrix representing distribution with P[i][j] = P(node = j | parent configuration = i)
	 * @throws Exception
	 *             if parent or child cannot be found in network
	 */
	public void setDistribution(String sName, double[][] P) throws Exception {
		int nTargetNode = getNode(sName);
		setDistribution(nTargetNode, P);
	} // setDistribution

	/** specify distribution of a node
	 * @param nTargetNode index of the node to specify distribution for
	 * @param P matrix representing distribution with P[i][j] = P(node = j | parent configuration = i)
	 * @throws Exception
	 *             if parent or child cannot be found in network
	 */
	public void setDistribution(int nTargetNode, double[][] P) throws Exception {
		// update undo stack
		if (m_bNeedsUndoAction) {
			addUndoAction(new SetDistributionAction(nTargetNode, P));
		}
		Estimator[] distributions = m_Distributions[nTargetNode];
		for (int iParent = 0; iParent < distributions.length; iParent++) {
			DiscreteEstimatorBayes distribution = new DiscreteEstimatorBayes(P[0].length, 0);
			for (int iValue = 0; iValue < distribution.getNumSymbols(); iValue++) {
				distribution.addValue(iValue, P[iParent][iValue]);
			}
			distributions[iParent] = distribution;
		}
		// m_Distributions[nTargetNode] = distributions;
	} // setDistribution

	/** returns distribution of a node in matrix form with matrix representing distribution
	 * with P[i][j] = P(node = j | parent configuration = i)
	 * @param sName name of the node to get distribution from
	 */
	public double[][] getDistribution(String sName) {
		int nTargetNode = getNode2(sName);
		return getDistribution(nTargetNode);
	} // getDistribution

	/** returns distribution of a node in matrix form with matrix representing distribution
	 * with P[i][j] = P(node = j | parent configuration = i)
	 * @param nTargetNode index of the node to get distribution from
	 */
	public double[][] getDistribution(int nTargetNode) {
		int nParentCard = m_ParentSets[nTargetNode].getCardinalityOfParents();
		int nCard = m_Instances.attribute(nTargetNode).numValues();
		double[][] P = new double[nParentCard][nCard];
		for (int iParent = 0; iParent < nParentCard; iParent++) {
			for (int iValue = 0; iValue < nCard; iValue++) {
				P[iParent][iValue] = m_Distributions[nTargetNode][iParent].getProbability(iValue);
			}
		}
		return P;
	} // getDistribution

	/** returns array of values of a node
	 * @param sName name of the node to get values from
	 */
	public String[] getValues(String sName) {
		int nTargetNode = getNode2(sName);
		return getValues(nTargetNode);
	} // getValues

	/** returns array of values of a node
	 * @param nTargetNode index of the node to get values from
	 */
	public String[] getValues(int nTargetNode) {
		String[] values = new String[getCardinality(nTargetNode)];
		for (int iValue = 0; iValue < values.length; iValue++) {
			values[iValue] = m_Instances.attribute(nTargetNode).value(iValue);
		}
		return values;
	} // getValues

	/** returns value of a node
	 * @param nTargetNode index of the node to get values from
	 * @param iValue index of the value
	 */
	public String getValueName(int nTargetNode, int iValue) {
		return m_Instances.attribute(nTargetNode).value(iValue);
	} // getNodeValue

	/** change the name of a node
	 * @param nTargetNode index of the node to set name for
	 * @param sName new name to assign
	 */
	public void setNodeName(int nTargetNode, String sName) {
		// update undo stack
		if (m_bNeedsUndoAction) {
			addUndoAction(new RenameAction(nTargetNode, getNodeName(nTargetNode), sName));
		}
		Attribute att = m_Instances.attribute(nTargetNode);
		int nCardinality = att.numValues();
		FastVector values = new FastVector(nCardinality);
		for (int iValue = 0; iValue < nCardinality; iValue++) {
			values.addElement(att.value(iValue));
		}
		replaceAtt(nTargetNode, sName, values);
	} // setNodeName

	/** change the name of a value of a node
	 * @param nTargetNode index of the node to set name for
	 * @param sValue current name of the value
	 * @param sNewValue new name of the value
	 */
	public void renameNodeValue(int nTargetNode, String sValue, String sNewValue) {
		// update undo stack
		if (m_bNeedsUndoAction) {
			addUndoAction(new RenameValueAction(nTargetNode, sValue, sNewValue));
		}
		Attribute att = m_Instances.attribute(nTargetNode);
		int nCardinality = att.numValues();
		FastVector values = new FastVector(nCardinality);
		for (int iValue = 0; iValue < nCardinality; iValue++) {
			if (att.value(iValue).equals(sValue)) {
				values.addElement(sNewValue);
			} else {
				values.addElement(att.value(iValue));
			}
		}
		replaceAtt(nTargetNode, att.name(), values);
	} // renameNodeValue


	/** Add node value to a node. Distributions for the node assign zero probability
	 * to the new value. Child nodes duplicate CPT conditioned on the new value.
	 * @param nTargetNode index of the node to add value for
	 * @param sNewValue name of the value
	 */
	public void addNodeValue(int nTargetNode, String sNewValue) {
		// update undo stack
		if (m_bNeedsUndoAction) {
			addUndoAction(new AddValueAction(nTargetNode, sNewValue));
		}
		Attribute att = m_Instances.attribute(nTargetNode);
		int nCardinality = att.numValues();
		FastVector values = new FastVector(nCardinality);
		for (int iValue = 0; iValue < nCardinality; iValue++) {
			values.addElement(att.value(iValue));
		}
		values.addElement(sNewValue);
		replaceAtt(nTargetNode, att.name(), values);

		// update distributions of this node
		Estimator[] distributions = m_Distributions[nTargetNode];
		int nNewCard = values.size();
		for (int iParent = 0; iParent < distributions.length; iParent++) {
			DiscreteEstimatorBayes distribution = new DiscreteEstimatorBayes(nNewCard, 0);
			for (int iValue = 0; iValue < nNewCard - 1; iValue++) {
				distribution.addValue(iValue, distributions[iParent].getProbability(iValue));
			}
			distributions[iParent] = distribution;
		}

		// update distributions of all children
		for (int iNode = 0; iNode < getNrOfNodes(); iNode++) {
			if (m_ParentSets[iNode].contains(nTargetNode)) {
				distributions = m_Distributions[iNode];
				ParentSet parentSet = m_ParentSets[iNode];
				int nParentCard = parentSet.getFreshCardinalityOfParents(m_Instances);
				Estimator[] newDistributions = new Estimator[nParentCard];
				int nCard = getCardinality(iNode);
				int nParents = parentSet.getNrOfParents();
				int[] values2 = new int[nParents];
				int iOldPos = 0;
				int iTargetNode = 0;
				while (parentSet.getParent(iTargetNode) != nTargetNode) {
					iTargetNode++;
				}
				for (int iPos = 0; iPos < nParentCard; iPos++) {
					DiscreteEstimatorBayes distribution = new DiscreteEstimatorBayes(nCard, 0);
					for (int iValue = 0; iValue < nCard; iValue++) {
						distribution.addValue(iValue, distributions[iOldPos].getProbability(iValue));
					}
					newDistributions[iPos] = distribution;
					// update values
					int i = 0;
					values2[i]++;
					while (i < nParents && values2[i] == getCardinality(parentSet.getParent(i))) {
						values2[i] = 0;
						i++;
						if (i < nParents) {
							values2[i]++;
						}
					}
					if (values2[iTargetNode] != nNewCard - 1) {
						iOldPos++;
					}
				}
				m_Distributions[iNode] = newDistributions;
			}
		}
	} // addNodeValue


	/** Delete node value from a node. Distributions for the node are scaled
	 * up proportional to existing distribution
	 * (or made uniform if zero probability is assigned to remainder of values).
	.* Child nodes delete CPTs conditioned on the new value.
	 * @param nTargetNode index of the node to delete value from
	 * @param sValue name of the value to delete
	 */
	public void delNodeValue(int nTargetNode, String sValue) throws Exception {
		// update undo stack
		if (m_bNeedsUndoAction) {
			addUndoAction(new DelValueAction(nTargetNode, sValue));
		}
		Attribute att = m_Instances.attribute(nTargetNode);
		int nCardinality = att.numValues();
		FastVector values = new FastVector(nCardinality);
		int nValue = -1;
		for (int iValue = 0; iValue < nCardinality; iValue++) {
			if (att.value(iValue).equals(sValue)) {
				nValue = iValue;
			} else {
				values.addElement(att.value(iValue));
			}
		}
		if (nValue < 0) {
			// could not find value
			throw new Exception("Node " + nTargetNode + " does not have value (" + sValue + ")");
		}
		replaceAtt(nTargetNode, att.name(), values);

		// update distributions
		Estimator[] distributions = m_Distributions[nTargetNode];
		int nCard = values.size();
		for (int iParent = 0; iParent < distributions.length; iParent++) {
			DiscreteEstimatorBayes distribution = new DiscreteEstimatorBayes(nCard, 0);
			double sum = 0;
			for (int iValue = 0; iValue < nCard; iValue++) {
				sum += distributions[iParent].getProbability(iValue);
			}
			if (sum > 0) {
				for (int iValue = 0; iValue < nCard; iValue++) {
					distribution.addValue(iValue, distributions[iParent].getProbability(iValue) / sum);
				}
			} else {
				for (int iValue = 0; iValue < nCard; iValue++) {
					distribution.addValue(iValue, 1.0 / nCard);
				}
			}
			distributions[iParent] = distribution;
		}

		// update distributions of all children
		for (int iNode = 0; iNode < getNrOfNodes(); iNode++) {
			if (m_ParentSets[iNode].contains(nTargetNode)) {
				ParentSet parentSet = m_ParentSets[iNode];
				distributions = m_Distributions[iNode];
				Estimator[] newDistributions = new Estimator[distributions.length * nCard / (nCard + 1)];
				int iCurrentDist = 0;

				int nParents = parentSet.getNrOfParents();
				int[] values2 = new int[nParents];
				// fill in the values
				int nParentCard = parentSet.getFreshCardinalityOfParents(m_Instances) * (nCard + 1) / nCard;
				int iTargetNode = 0;
				while (parentSet.getParent(iTargetNode) != nTargetNode) {
					iTargetNode++;
				}
				int[] nCards = new int[nParents];
				for (int iParent = 0; iParent < nParents; iParent++) {
					nCards[iParent] = getCardinality(parentSet.getParent(iParent));
				}
				nCards[iTargetNode]++;
				for (int iPos = 0; iPos < nParentCard; iPos++) {
					if (values2[iTargetNode] != nValue) {
						newDistributions[iCurrentDist++] = distributions[iPos];
					}
					// update values
					int i = 0;
					values2[i]++;
					while (i < nParents && values2[i] == nCards[i]) {
						values2[i] = 0;
						i++;
						if (i < nParents) {
							values2[i]++;
						}
					}
				}

				m_Distributions[iNode] = newDistributions;
			}
		}
		// update evidence
		if (getEvidence(nTargetNode) > nValue) {
			setEvidence(nTargetNode, getEvidence(nTargetNode) - 1);
		}
	} // delNodeValue

	/** set position of node
	 * @param iNode index of node to set position for
	 * @param nX x position of new position
	 * @param nY y position of new position
	 */
	public void setPosition(int iNode, int nX, int nY) {
		// update undo stack
		if (m_bNeedsUndoAction) {
			boolean isUpdate = false;
			UndoAction undoAction = null;
			try {
				if (m_undoStack.size() > 0) {
					undoAction = (UndoAction) m_undoStack.elementAt(m_undoStack.size() - 1);
					SetPositionAction posAction = (SetPositionAction) undoAction;
					if (posAction.m_nTargetNode == iNode) {
						isUpdate = true;
						posAction.setUndoPosition(nX, nY);
					}
				}
			} catch (Exception e) {
				// ignore. it's not a SetPositionAction
			}
			if (!isUpdate) {
				addUndoAction(new SetPositionAction(iNode, nX, nY));
			}
		}
		m_nPositionX.setElementAt(nX, iNode);
		m_nPositionY.setElementAt(nY, iNode);
	} // setPosition

	/** Set position of node. Move set of nodes with the same displacement
	 * as a specified node.
	 * @param nNode index of node to set position for
	 * @param nX x position of new position
	 * @param nY y position of new position
	 * @param nodes array of indexes of nodes to move
	 */
	public void setPosition(int nNode, int nX, int nY, FastVector nodes) {
		int dX = nX - getPositionX(nNode);
		int dY = nY - getPositionY(nNode);
		// update undo stack
		if (m_bNeedsUndoAction) {
			boolean isUpdate = false;
			try {
				UndoAction undoAction = null;
				if (m_undoStack.size() > 0) {
					undoAction = (UndoAction) m_undoStack.elementAt(m_undoStack.size() - 1);
						SetGroupPositionAction posAction = (SetGroupPositionAction) undoAction;
						isUpdate = true;
						int iNode = 0;
						while (isUpdate && iNode < posAction.m_nodes.size()) {
							if ((Integer)posAction.m_nodes.elementAt(iNode) != (Integer) nodes.elementAt(iNode)) {
								isUpdate = false;
							}
							iNode++;
						}
						if (isUpdate == true) {
							posAction.setUndoPosition(dX, dY);
						}
				}
			} catch (Exception e) {
				// ignore. it's not a SetPositionAction
			}
			if (!isUpdate) {
				addUndoAction(new SetGroupPositionAction(nodes, dX, dY));
			}
		}
		for (int iNode = 0; iNode < nodes.size(); iNode++) {
			nNode = (Integer) nodes.elementAt(iNode);
			m_nPositionX.setElementAt(getPositionX(nNode) + dX, nNode);
			m_nPositionY.setElementAt(getPositionY(nNode) + dY, nNode);
		}
	} // setPosition

	/** set positions of all nodes
	 * @param nPosX new x positions for all nodes
	 * @param nPosY new y positions for all nodes
	 */
	public void layoutGraph(FastVector nPosX, FastVector nPosY) {
		if (m_bNeedsUndoAction) {
			addUndoAction(new LayoutGraphAction(nPosX, nPosY));
		}
		m_nPositionX = nPosX;
		m_nPositionY = nPosY;
	} // layoutGraph

	/** get x position of a node
	 * @param iNode index of node of interest
	 */
	public int getPositionX(int iNode) {
		return (Integer) (m_nPositionX.elementAt(iNode));
	}

	/** get y position of a node
	 * @param iNode index of node of interest
	 */
	public int getPositionY(int iNode) {
		return (Integer) (m_nPositionY.elementAt(iNode));
	}

	/** align set of nodes with the left most node in the list
	 * @param nodes list of indexes of nodes to align
	 */
	public void alignLeft(FastVector nodes) {
		// update undo stack
		if (m_bNeedsUndoAction) {
			addUndoAction(new alignLeftAction(nodes));
		}
		int nMinX = -1;
		for (int iNode = 0; iNode < nodes.size(); iNode++) {
			int nX = getPositionX((Integer) nodes.elementAt(iNode));
			if (nX < nMinX || iNode == 0) {
				nMinX = nX;
			}
		}
		for (int iNode = 0; iNode < nodes.size(); iNode++) {
			int nNode = (Integer) nodes.elementAt(iNode);
			m_nPositionX.setElementAt(nMinX, nNode);
		}
	} // alignLeft

	/** align set of nodes with the right most node in the list
	 * @param nodes list of indexes of nodes to align
	 */
	public void alignRight(FastVector nodes) {
		// update undo stack
		if (m_bNeedsUndoAction) {
			addUndoAction(new alignRightAction(nodes));
		}
		int nMaxX = -1;
		for (int iNode = 0; iNode < nodes.size(); iNode++) {
			int nX = getPositionX((Integer) nodes.elementAt(iNode));
			if (nX > nMaxX || iNode == 0) {
				nMaxX = nX;
			}
		}
		for (int iNode = 0; iNode < nodes.size(); iNode++) {
			int nNode = (Integer) nodes.elementAt(iNode);
			m_nPositionX.setElementAt(nMaxX, nNode);
		}
	} // alignRight

	/** align set of nodes with the top most node in the list
	 * @param nodes list of indexes of nodes to align
	 */
	public void alignTop(FastVector nodes) {
		// update undo stack
		if (m_bNeedsUndoAction) {
			addUndoAction(new alignTopAction(nodes));
		}
		int nMinY = -1;
		for (int iNode = 0; iNode < nodes.size(); iNode++) {
			int nY = getPositionY((Integer) nodes.elementAt(iNode));
			if (nY < nMinY || iNode == 0) {
				nMinY = nY;
			}
		}
		for (int iNode = 0; iNode < nodes.size(); iNode++) {
			int nNode = (Integer) nodes.elementAt(iNode);
			m_nPositionY.setElementAt(nMinY, nNode);
		}
	} // alignTop

	/** align set of nodes with the bottom most node in the list
	 * @param nodes list of indexes of nodes to align
	 */
	public void alignBottom(FastVector nodes) {
		// update undo stack
		if (m_bNeedsUndoAction) {
			addUndoAction(new alignBottomAction(nodes));
		}
		int nMaxY = -1;
		for (int iNode = 0; iNode < nodes.size(); iNode++) {
			int nY = getPositionY((Integer) nodes.elementAt(iNode));
			if (nY > nMaxY || iNode == 0) {
				nMaxY = nY;
			}
		}
		for (int iNode = 0; iNode < nodes.size(); iNode++) {
			int nNode = (Integer) nodes.elementAt(iNode);
			m_nPositionY.setElementAt(nMaxY, nNode);
		}
	} // alignBottom

	/** center set of nodes half way between left and right most node in the list
	 * @param nodes list of indexes of nodes to center
	 */
	public void centerHorizontal(FastVector nodes) {
		// update undo stack
		if (m_bNeedsUndoAction) {
			addUndoAction(new centerHorizontalAction(nodes));
		}
		int nMinY = -1;
		int nMaxY = -1;
		for (int iNode = 0; iNode < nodes.size(); iNode++) {
			int nY = getPositionY((Integer) nodes.elementAt(iNode));
			if (nY < nMinY || iNode == 0) {
				nMinY = nY;
			}
			if (nY > nMaxY || iNode == 0) {
				nMaxY = nY;
			}
		}
		for (int iNode = 0; iNode < nodes.size(); iNode++) {
			int nNode = (Integer) nodes.elementAt(iNode);
			m_nPositionY.setElementAt((nMinY + nMaxY) / 2, nNode);
		}
	} // centerHorizontal

	/** center set of nodes half way between top and bottom most node in the list
	 * @param nodes list of indexes of nodes to center
	 */
	public void centerVertical(FastVector nodes) {
		// update undo stack
		if (m_bNeedsUndoAction) {
			addUndoAction(new centerVerticalAction(nodes));
		}
		int nMinX = -1;
		int nMaxX = -1;
		for (int iNode = 0; iNode < nodes.size(); iNode++) {
			int nX = getPositionX((Integer) nodes.elementAt(iNode));
			if (nX < nMinX || iNode == 0) {
				nMinX = nX;
			}
			if (nX > nMaxX || iNode == 0) {
				nMaxX = nX;
			}
		}
		for (int iNode = 0; iNode < nodes.size(); iNode++) {
			int nNode = (Integer) nodes.elementAt(iNode);
			m_nPositionX.setElementAt((nMinX + nMaxX) / 2, nNode);
		}
	} // centerVertical

	/** space out set of nodes evenly between left and right most node in the list
	 * @param nodes list of indexes of nodes to space out
	 */
	public void spaceHorizontal(FastVector nodes) {
		// update undo stack
		if (m_bNeedsUndoAction) {
			addUndoAction(new spaceHorizontalAction(nodes));
		}
		int nMinX = -1;
		int nMaxX = -1;
		for (int iNode = 0; iNode < nodes.size(); iNode++) {
			int nX = getPositionX((Integer) nodes.elementAt(iNode));
			if (nX < nMinX || iNode == 0) {
				nMinX = nX;
			}
			if (nX > nMaxX || iNode == 0) {
				nMaxX = nX;
			}
		}
		for (int iNode = 0; iNode < nodes.size(); iNode++) {
			int nNode = (Integer) nodes.elementAt(iNode);
			m_nPositionX.setElementAt((int) (nMinX + iNode * (nMaxX - nMinX) / (nodes.size() - 1.0)), nNode);
		}
	} // spaceHorizontal

	/** space out set of nodes evenly between top and bottom most node in the list
	 * @param nodes list of indexes of nodes to space out
	 */
	public void spaceVertical(FastVector nodes) {
		// update undo stack
		if (m_bNeedsUndoAction) {
			addUndoAction(new spaceVerticalAction(nodes));
		}
		int nMinY = -1;
		int nMaxY = -1;
		for (int iNode = 0; iNode < nodes.size(); iNode++) {
			int nY = getPositionY((Integer) nodes.elementAt(iNode));
			if (nY < nMinY || iNode == 0) {
				nMinY = nY;
			}
			if (nY > nMaxY || iNode == 0) {
				nMaxY = nY;
			}
		}
		for (int iNode = 0; iNode < nodes.size(); iNode++) {
			int nNode = (Integer) nodes.elementAt(iNode);
			m_nPositionY.setElementAt((int) (nMinY + iNode * (nMaxY - nMinY) / (nodes.size() - 1.0)), nNode);
		}
	} // spaceVertical


	/** replace attribute with specified name and values
	 * @param nTargetNode index of node the replace specification for
	 * @param sName new name of the node
	 * @param values array of values of the node
	 */
	void replaceAtt(int nTargetNode, String sName, FastVector values) {
		Attribute newAtt = new Attribute(sName, values);
		if (m_Instances.classIndex() == nTargetNode) {
			m_Instances.setClassIndex(-1);
			m_Instances.insertAttributeAt(newAtt, nTargetNode);
			m_Instances.deleteAttributeAt(nTargetNode + 1);
			m_Instances.setClassIndex(nTargetNode);
		} else {
			m_Instances.insertAttributeAt(newAtt, nTargetNode);
			m_Instances.deleteAttributeAt(nTargetNode + 1);
		}
	} // replaceAtt

	/** return marginal distibution for a node
	 * @param iNode index of node of interest
	 */
	public double[] getMargin(int iNode) {
		return (double[]) m_fMarginP.elementAt(iNode);
	};

	/** set marginal distibution for a node
	 * @param iNode index of node to set marginal distribution for
	 * @param fMarginP marginal distribution
	 */
	public void setMargin(int iNode, double[] fMarginP) {
		m_fMarginP.setElementAt(fMarginP, iNode);
	}

	/** get evidence state of a node. -1 represents no evidence set, otherwise
	 * the index of a value of the node
	 * @param iNode index of node of interest
	 */
	public int getEvidence(int iNode) {
		return (Integer) m_nEvidence.elementAt(iNode);
	}

	/** set evidence state of a node. -1 represents no evidence set, otherwise
	 * the index of a value of the node
 	 * @param iNode index of node of interest
	 * @param iValue evidence value to set
	 */
	public void setEvidence(int iNode, int iValue) {
		m_nEvidence.setElementAt(iValue, iNode);
	}

	/** return list of children of a node
	 * @param nTargetNode index of node of interest
	 */
	public FastVector getChildren(int nTargetNode) {
		FastVector children = new FastVector();
		for (int iNode = 0; iNode < getNrOfNodes(); iNode++) {
			if (m_ParentSets[iNode].contains(nTargetNode)) {
				children.addElement(iNode);
			}
		}
		return children;
	} // getChildren

	/** returns network in XMLBIF format
	*/
	public String toXMLBIF03() {
		if (m_Instances == null) {
			return ("<!--No model built yet-->");
		}

		StringBuffer text = new StringBuffer();
		text.append(getBIFHeader());
		text.append("\n");
		text.append("\n");
		text.append("<BIF VERSION=\"0.3\">\n");
		text.append("<NETWORK>\n");
		text.append("<NAME>" + XMLNormalize(m_Instances.relationName()) + "</NAME>\n");
		for (int iAttribute = 0; iAttribute < m_Instances.numAttributes(); iAttribute++) {
			text.append("<VARIABLE TYPE=\"nature\">\n");
			text.append("<NAME>" + XMLNormalize(m_Instances.attribute(iAttribute).name()) + "</NAME>\n");
			for (int iValue = 0; iValue < m_Instances.attribute(iAttribute).numValues(); iValue++) {
				text.append("<OUTCOME>" + XMLNormalize(m_Instances.attribute(iAttribute).value(iValue))
						+ "</OUTCOME>\n");
			}
			text.append("<PROPERTY>position = (" + getPositionX(iAttribute) + "," + getPositionY(iAttribute)
					+ ")</PROPERTY>\n");
			text.append("</VARIABLE>\n");
		}

		for (int iAttribute = 0; iAttribute < m_Instances.numAttributes(); iAttribute++) {
			text.append("<DEFINITION>\n");
			text.append("<FOR>" + XMLNormalize(m_Instances.attribute(iAttribute).name()) + "</FOR>\n");
			for (int iParent = 0; iParent < m_ParentSets[iAttribute].getNrOfParents(); iParent++) {
				text.append("<GIVEN>"
						+ XMLNormalize(m_Instances.attribute(m_ParentSets[iAttribute].getParent(iParent)).name())
						+ "</GIVEN>\n");
			}
			text.append("<TABLE>\n");
			for (int iParent = 0; iParent < m_ParentSets[iAttribute].getCardinalityOfParents(); iParent++) {
				for (int iValue = 0; iValue < m_Instances.attribute(iAttribute).numValues(); iValue++) {
					text.append(m_Distributions[iAttribute][iParent].getProbability(iValue));
					text.append(' ');
				}
				text.append('\n');
			}
			text.append("</TABLE>\n");
			text.append("</DEFINITION>\n");
		}
		text.append("</NETWORK>\n");
		text.append("</BIF>\n");
		return text.toString();
	} // toXMLBIF03

	/** return fragment of network in XMLBIF format
	 * @param nodes array of indexes of nodes that should be in the fragment
	 */
	public String toXMLBIF03(FastVector nodes) {
		StringBuffer text = new StringBuffer();
		text.append(getBIFHeader());
		text.append("\n");
		text.append("\n");
		text.append("<BIF VERSION=\"0.3\">\n");
		text.append("<NETWORK>\n");
		text.append("<NAME>" + XMLNormalize(m_Instances.relationName()) + "</NAME>\n");
		for (int iNode = 0; iNode < nodes.size(); iNode++) {
			int nNode = (Integer) nodes.elementAt(iNode);
			text.append("<VARIABLE TYPE=\"nature\">\n");
			text.append("<NAME>" + XMLNormalize(m_Instances.attribute(nNode).name()) + "</NAME>\n");
			for (int iValue = 0; iValue < m_Instances.attribute(nNode).numValues(); iValue++) {
				text.append("<OUTCOME>" + XMLNormalize(m_Instances.attribute(nNode).value(iValue)) + "</OUTCOME>\n");
			}
			text.append("<PROPERTY>position = (" + getPositionX(nNode) + "," + getPositionY(nNode) + ")</PROPERTY>\n");
			text.append("</VARIABLE>\n");
		}

		for (int iNode = 0; iNode < nodes.size(); iNode++) {
			int nNode = (Integer) nodes.elementAt(iNode);
			text.append("<DEFINITION>\n");
			text.append("<FOR>" + XMLNormalize(m_Instances.attribute(nNode).name()) + "</FOR>\n");
			for (int iParent = 0; iParent < m_ParentSets[nNode].getNrOfParents(); iParent++) {
				text.append("<GIVEN>"
						+ XMLNormalize(m_Instances.attribute(m_ParentSets[nNode].getParent(iParent)).name())
						+ "</GIVEN>\n");
			}
			text.append("<TABLE>\n");
			for (int iParent = 0; iParent < m_ParentSets[nNode].getCardinalityOfParents(); iParent++) {
				for (int iValue = 0; iValue < m_Instances.attribute(nNode).numValues(); iValue++) {
					text.append(m_Distributions[nNode][iParent].getProbability(iValue));
					text.append(' ');
				}
				text.append('\n');
			}
			text.append("</TABLE>\n");
			text.append("</DEFINITION>\n");
		}
		text.append("</NETWORK>\n");
		text.append("</BIF>\n");
		return text.toString();
	} // toXMLBIF03

	/** undo stack for undoin edit actions, or redo edit actions */
	FastVector m_undoStack = new FastVector();

	/** current action in undo stack */
	int m_nCurrentEditAction = -1;

	/** action that the network is saved */
	int m_nSavedPointer = -1;

	/***************************************************************************
	 * flag to indicate whether an edit action needs to introduce an undo
	 * action. This is only false when an undo or redo action is performed.
	 **************************************************************************/
	boolean m_bNeedsUndoAction = true;

	/** return whether there is something on the undo stack that can be performed */
	public boolean canUndo() {
		return m_nCurrentEditAction >= 0;
	}

	/** return whether there is something on the undo stack that can be performed */
	public boolean canRedo() {
		return m_nCurrentEditAction < m_undoStack.size() - 1;
	}

	/** return true when current state differs from the state the network was last saved */
	public boolean isChanged() {
		return m_nCurrentEditAction != m_nSavedPointer;
	}

	/** indicate the network state was saved */
	public void isSaved() {
		m_nSavedPointer = m_nCurrentEditAction;
	}

	/** get message representing the last action performed on the network */
	public String lastActionMsg() {
		if (m_undoStack.size() == 0) {
			return "";
		}
		return ((UndoAction) m_undoStack.lastElement()).getRedoMsg();
	} // lastActionMsg


	/** undo the last edit action performed on the network.
	 * returns message representing the action performed.
	 */
	public String undo() {
		if (!canUndo()) {
			return "";
		}
		UndoAction undoAction = (UndoAction) m_undoStack.elementAt(m_nCurrentEditAction);
		m_bNeedsUndoAction = false;
		undoAction.undo();
		m_bNeedsUndoAction = true;
		m_nCurrentEditAction--;

		// undo stack debugging
		/*
		if (m_nCurrentEditAction>0) {
			String sXML = (String) m_sXMLStack.elementAt(m_nCurrentEditAction);
			String sXMLCurrent = toXMLBIF03();
			if (!sXML.equals(sXMLCurrent)) {
				String sDiff = "";
				String sDiff2 = "";
				for (int i = 0; i < sXML.length() && sDiff.length() < 80; i++) {
					if (sXML.charAt(i) != sXMLCurrent.charAt(i)) {
						sDiff += sXML.charAt(i);
						sDiff2 += sXMLCurrent.charAt(i);
					}
				}

				JOptionPane.showMessageDialog(null,"Undo error\n" + sDiff + " \n" + sDiff2);
			}
		}
		*/
		return undoAction.getUndoMsg();
	} // undo

	/** redo the last edit action performed on the network.
	 * returns message representing the action performed.
	 */
	public String redo() {
		if (!canRedo()) {
			return "";
		}
		m_nCurrentEditAction++;
		UndoAction undoAction = (UndoAction) m_undoStack.elementAt(m_nCurrentEditAction);
		m_bNeedsUndoAction = false;
		undoAction.redo();
		m_bNeedsUndoAction = true;

		// undo stack debugging
		/*
		if (m_nCurrentEditAction < m_sXMLStack.size()) {
			String sXML = (String) m_sXMLStack.elementAt(m_nCurrentEditAction);
			String sXMLCurrent = toXMLBIF03();
			if (!sXML.equals(sXMLCurrent)) {
				String sDiff = "";
				String sDiff2 = "";
				for (int i = 0; i < sXML.length() && sDiff.length() < 80; i++) {
					if (sXML.charAt(i) != sXMLCurrent.charAt(i)) {
						sDiff += sXML.charAt(i);
						sDiff2 += sXMLCurrent.charAt(i);
					}
				}

				JOptionPane.showMessageDialog(null,"redo error\n" + sDiff + " \n" + sDiff2);
			}
		}
		*/
		return undoAction.getRedoMsg();
	} // redo

	/** add undo action to the undo stack.
	 * @param action operation that needs to be added to the undo stack
	 */
	void addUndoAction(UndoAction action) {
		int iAction = m_undoStack.size() - 1;
		while (iAction > m_nCurrentEditAction) {
			m_undoStack.removeElementAt(iAction--);
		}
		if (m_nSavedPointer > m_nCurrentEditAction) {
			m_nSavedPointer = -2;
		}
		m_undoStack.addElement(action);
		//m_sXMLStack.addElement(toXMLBIF03());
		m_nCurrentEditAction++;
	} // addUndoAction

	/** remove all actions from the undo stack */
	public void clearUndoStack() {
		m_undoStack = new FastVector();
		//m_sXMLStack = new FastVector();
		m_nCurrentEditAction = -1;
		m_nSavedPointer = -1;
	} // clearUndoStack

	/** base class for actions representing operations on the Bayesian network
	 * that can be undone/redone
	 */
	class UndoAction implements Serializable {
		/** for serialization */
		static final long serialVersionUID = 1;
		public void undo() {
		}

		public void redo() {
		}

		public String getUndoMsg() {
			return getMsg();
		}

		public String getRedoMsg() {
			return getMsg();
		}
		String getMsg() {
			String sStr = toString();
			int iStart = sStr.indexOf('$');
			int iEnd = sStr.indexOf('@');
			StringBuffer sBuffer = new StringBuffer();
			for(int i= iStart + 1; i < iEnd; i++) {
				char c = sStr.charAt(i);
				if (Character.isUpperCase(c)) {
					sBuffer.append(' ');
				}
				sBuffer.append(sStr.charAt(i));
			}
			return sBuffer.toString();
		} // getMsg
	} // class UndoAction

	class AddNodeAction extends UndoAction {
		/** for serialization */
		static final long serialVersionUID = 1;
		String m_sName;
		int m_nPosX;
		int m_nPosY;

		int m_nCardinality;

		AddNodeAction(String sName, int nCardinality, int nPosX, int nPosY) {
			m_sName = sName;
			m_nCardinality = nCardinality;
			m_nPosX = nPosX;
			m_nPosY = nPosY;
		} // c'tor

		public void undo() {
			try {
				deleteNode(getNrOfNodes() - 1);
			} catch (Exception e) {
				e.printStackTrace();
			}
		} // undo

		public void redo() {
			try {
				addNode(m_sName, m_nCardinality, m_nPosX, m_nPosY);
			} catch (Exception e) {
				e.printStackTrace();
			}
		} // redo
	} // class AddNodeAction

	class DeleteNodeAction extends UndoAction {
		/** for serialization */
		static final long serialVersionUID = 1;
		int m_nTargetNode;

		Attribute m_att;

		Estimator[] m_CPT;

		ParentSet m_ParentSet;

		FastVector m_deleteArcActions;

		int m_nPosX;

		int m_nPosY;

		DeleteNodeAction(int nTargetNode) {
			m_nTargetNode = nTargetNode;
			m_att = m_Instances.attribute(nTargetNode);
			try {
				SerializedObject so = new SerializedObject(m_Distributions[nTargetNode]);
				m_CPT = (Estimator[]) so.getObject();
				;
				so = new SerializedObject(m_ParentSets[nTargetNode]);
				m_ParentSet = (ParentSet) so.getObject();
			} catch (Exception e) {
				e.printStackTrace();
			}
			m_deleteArcActions = new FastVector();
			for (int iNode = 0; iNode < getNrOfNodes(); iNode++) {
				if (m_ParentSets[iNode].contains(nTargetNode)) {
					m_deleteArcActions.addElement(new DeleteArcAction(nTargetNode, iNode));
				}
			}
			m_nPosX = getPositionX(m_nTargetNode);
			m_nPosY = getPositionY(m_nTargetNode);
		} // c'tor

		public void undo() {
			try {
				m_Instances.insertAttributeAt(m_att, m_nTargetNode);
				int nAtts = m_Instances.numAttributes();
				// update parentsets
				ParentSet[] parentSets = new ParentSet[nAtts];
				int nX = 0;
				for (int iParentSet = 0; iParentSet < nAtts; iParentSet++) {
					if (iParentSet == m_nTargetNode) {
						SerializedObject so = new SerializedObject(m_ParentSet);
						parentSets[iParentSet] = (ParentSet) so.getObject();
						nX = 1;
					} else {
						parentSets[iParentSet] = m_ParentSets[iParentSet - nX];
						for (int iParent = 0; iParent < parentSets[iParentSet].getNrOfParents(); iParent++) {
							int nParent = parentSets[iParentSet].getParent(iParent);
							if (nParent >= m_nTargetNode) {
								parentSets[iParentSet].SetParent(iParent, nParent + 1);
							}
						}
					}
				}
				m_ParentSets = parentSets;
				// update distributions
				Estimator[][] distributions = new Estimator[nAtts][];
				nX = 0;
				for (int iNode = 0; iNode < nAtts; iNode++) {
					if (iNode == m_nTargetNode) {
						SerializedObject so = new SerializedObject(m_CPT);
						distributions[iNode] = (Estimator[]) so.getObject();
						nX = 1;
					} else {
						distributions[iNode] = m_Distributions[iNode - nX];
					}
				}
				m_Distributions = distributions;

				for (int deletedArc = 0; deletedArc < m_deleteArcActions.size(); deletedArc++) {
					DeleteArcAction action = (DeleteArcAction) m_deleteArcActions.elementAt(deletedArc);
					action.undo();
				}
				m_nPositionX.insertElementAt(m_nPosX, m_nTargetNode);
				m_nPositionY.insertElementAt(m_nPosY, m_nTargetNode);
				m_nEvidence.insertElementAt(-1, m_nTargetNode);
				m_fMarginP.insertElementAt(new double[getCardinality(m_nTargetNode)], m_nTargetNode);
			} catch (Exception e) {
				e.printStackTrace();
			}
		} // undo

		public void redo() {
			try {
				deleteNode(m_nTargetNode);
			} catch (Exception e) {
				e.printStackTrace();
			}
		} // redo
	} // class DeleteNodeAction

	class DeleteSelectionAction extends UndoAction {
		/** for serialization */
		static final long serialVersionUID = 1;
		FastVector m_nodes;

		Attribute[] m_att;

		Estimator[][] m_CPT;

		ParentSet[] m_ParentSet;

		FastVector m_deleteArcActions;

		int[] m_nPosX;

		int[] m_nPosY;

		public DeleteSelectionAction(FastVector nodes) {
			m_nodes = new FastVector();
			int nNodes = nodes.size();
			m_att = new Attribute[nNodes];
			m_CPT = new Estimator[nNodes][];
			m_ParentSet = new ParentSet[nNodes];
			m_nPosX = new int[nNodes];
			m_nPosY = new int[nNodes];
			m_deleteArcActions = new FastVector();
			for (int iNode = 0; iNode < nodes.size(); iNode++) {
				int nTargetNode = (Integer) nodes.elementAt(iNode);
				m_nodes.addElement(nTargetNode);
				m_att[iNode] = m_Instances.attribute(nTargetNode);
				try {
					SerializedObject so = new SerializedObject(m_Distributions[nTargetNode]);
					m_CPT[iNode] = (Estimator[]) so.getObject();
					;
					so = new SerializedObject(m_ParentSets[nTargetNode]);
					m_ParentSet[iNode] = (ParentSet) so.getObject();
				} catch (Exception e) {
					e.printStackTrace();
				}
				m_nPosX[iNode] = getPositionX(nTargetNode);
				m_nPosY[iNode] = getPositionY(nTargetNode);
				for (int iNode2 = 0; iNode2 < getNrOfNodes(); iNode2++) {
					if (!nodes.contains(iNode2) && m_ParentSets[iNode2].contains(nTargetNode)) {
						m_deleteArcActions.addElement(new DeleteArcAction(nTargetNode, iNode2));
					}
				}
			}
		} // c'tor

		public void undo() {
			try {
				for (int iNode = 0; iNode < m_nodes.size(); iNode++) {
					int nTargetNode = (Integer) m_nodes.elementAt(iNode);
					m_Instances.insertAttributeAt(m_att[iNode], nTargetNode);
				}
				int nAtts = m_Instances.numAttributes();
				// update parentsets
				ParentSet[] parentSets = new ParentSet[nAtts];
				int[] offset = new int[nAtts];
				for (int iNode = 0; iNode < nAtts; iNode++) {
					offset[iNode] = iNode;
				}
				for (int iNode = m_nodes.size() - 1; iNode >= 0; iNode--) {
					int nTargetNode = (Integer) m_nodes.elementAt(iNode);
					for (int i = nTargetNode; i < nAtts - 1; i++) {
						offset[i] = offset[i + 1];
					}
				}

				int iTargetNode = 0;
				for (int iParentSet = 0; iParentSet < nAtts; iParentSet++) {
					if (iTargetNode < m_nodes.size()
							&& (Integer) m_nodes.elementAt(iTargetNode) == (Integer) iParentSet) {
						SerializedObject so = new SerializedObject(m_ParentSet[iTargetNode]);
						parentSets[iParentSet] = (ParentSet) so.getObject();
						iTargetNode++;
					} else {
						parentSets[iParentSet] = m_ParentSets[iParentSet - iTargetNode];
						for (int iParent = 0; iParent < parentSets[iParentSet].getNrOfParents(); iParent++) {
							int nParent = parentSets[iParentSet].getParent(iParent);
							parentSets[iParentSet].SetParent(iParent, offset[nParent]);
						}
					}
				}
				m_ParentSets = parentSets;
				// update distributions
				Estimator[][] distributions = new Estimator[nAtts][];
				iTargetNode = 0;
				for (int iNode = 0; iNode < nAtts; iNode++) {
					if (iTargetNode < m_nodes.size() && (Integer) m_nodes.elementAt(iTargetNode) == (Integer) iNode) {
						SerializedObject so = new SerializedObject(m_CPT[iTargetNode]);
						distributions[iNode] = (Estimator[]) so.getObject();
						iTargetNode++;
					} else {
						distributions[iNode] = m_Distributions[iNode - iTargetNode];
					}
				}
				m_Distributions = distributions;

				for (int iNode = 0; iNode < m_nodes.size(); iNode++) {
					int nTargetNode = (Integer) m_nodes.elementAt(iNode);
					m_nPositionX.insertElementAt(m_nPosX[iNode], nTargetNode);
					m_nPositionY.insertElementAt(m_nPosY[iNode], nTargetNode);
					m_nEvidence.insertElementAt(-1, nTargetNode);
					m_fMarginP.insertElementAt(new double[getCardinality(nTargetNode)], nTargetNode);
				}
				for (int deletedArc = 0; deletedArc < m_deleteArcActions.size(); deletedArc++) {
					DeleteArcAction action = (DeleteArcAction) m_deleteArcActions.elementAt(deletedArc);
					action.undo();
				}
			} catch (Exception e) {
				e.printStackTrace();
			}
		} // undo

		public void redo() {
			try {
				for (int iNode = m_nodes.size() - 1; iNode >= 0; iNode--) {
					int nNode = (Integer) m_nodes.elementAt(iNode);
					deleteNode(nNode);
				}
			} catch (Exception e) {
				e.printStackTrace();
			}
		} // redo
	} // class DeleteSelectionAction

	class AddArcAction extends UndoAction {
		/** for serialization */
		static final long serialVersionUID = 1;
		//int m_nChild;
		FastVector m_children;

		int m_nParent;

		Estimator[][] m_CPT;

		AddArcAction(int nParent, int nChild) {
			try {
				m_nParent = nParent;
				m_children = new FastVector();
				m_children.addElement(nChild);
				//m_nChild = nChild;
				SerializedObject so = new SerializedObject(m_Distributions[nChild]);
				m_CPT = new Estimator[1][];
				m_CPT[0] = (Estimator[]) so.getObject();
				;
			} catch (Exception e) {
				e.printStackTrace();
			}
		} // c'tor

		AddArcAction(int nParent, FastVector children) {
			try {
				m_nParent = nParent;
				m_children = new FastVector();
				m_CPT = new Estimator[children.size()][];
				for (int iChild = 0; iChild < children.size(); iChild++) {
					int nChild = (Integer) children.elementAt(iChild);
					m_children.addElement(nChild);
					SerializedObject so = new SerializedObject(m_Distributions[nChild]);
					m_CPT[iChild] = (Estimator[]) so.getObject();
				}
			} catch (Exception e) {
				e.printStackTrace();
			}
		} // c'tor

		public void undo() {
			try {
				for (int iChild = 0; iChild < m_children.size(); iChild++) {
					int nChild = (Integer) m_children.elementAt(iChild);
					deleteArc(m_nParent, nChild);
					SerializedObject so = new SerializedObject(m_CPT[iChild]);
					m_Distributions[nChild] = (Estimator[]) so.getObject();
				}
			} catch (Exception e) {
				e.printStackTrace();
			}
		} // undo

		public void redo() {
			try {
				for (int iChild = 0; iChild < m_children.size(); iChild++) {
					int nChild = (Integer) m_children.elementAt(iChild);
					addArc(m_nParent, nChild);
				}
			} catch (Exception e) {
				e.printStackTrace();
			}
		} // redo
	} // class AddArcAction

	class DeleteArcAction extends UndoAction {
		/** for serialization */
		static final long serialVersionUID = 1;
		int[] m_nParents;
		int m_nChild;
		int m_nParent;
		Estimator[] m_CPT;

		DeleteArcAction(int nParent, int nChild) {
			try {
			m_nChild = nChild;
			m_nParent = nParent;
			m_nParents = new int[getNrOfParents(nChild)];
			for (int iParent = 0; iParent < m_nParents.length; iParent++) {
				m_nParents[iParent] = getParent(nChild, iParent);
			}
			SerializedObject so = new SerializedObject(m_Distributions[nChild]);
			m_CPT = (Estimator[]) so.getObject();
			} catch (Exception e) {
				e.printStackTrace();
			}
		} // c'tor

		public void undo() {
			try {
				SerializedObject so = new SerializedObject(m_CPT);
				m_Distributions[m_nChild] = (Estimator[]) so.getObject();
				ParentSet parentSet = new ParentSet();
				for (int iParent = 0; iParent < m_nParents.length; iParent++) {
					parentSet.addParent(m_nParents[iParent], m_Instances);
				}
				m_ParentSets[m_nChild] = parentSet;
			} catch (Exception e) {
				e.printStackTrace();
			}
		} // undo

		public void redo() {
			try {
				deleteArc(m_nParent, m_nChild);
			} catch (Exception e) {
				e.printStackTrace();
			}
		} // redo
	} // class DeleteArcAction

	class SetDistributionAction extends UndoAction {
		/** for serialization */
		static final long serialVersionUID = 1;
		int m_nTargetNode;

		Estimator[] m_CPT;

		double[][] m_P;

		SetDistributionAction(int nTargetNode, double[][] P) {
			try {
				m_nTargetNode = nTargetNode;
				SerializedObject so = new SerializedObject(m_Distributions[nTargetNode]);
				m_CPT = (Estimator[]) so.getObject();
				;
				m_P = P;
			} catch (Exception e) {
				e.printStackTrace();
			}
		} // c'tor

		public void undo() {
			try {
				SerializedObject so = new SerializedObject(m_CPT);
				m_Distributions[m_nTargetNode] = (Estimator[]) so.getObject();
			} catch (Exception e) {
				e.printStackTrace();
			}
		} // undo

		public void redo() {
			try {
				setDistribution(m_nTargetNode, m_P);
			} catch (Exception e) {
				e.printStackTrace();
			}
		} // redo

		public String getUndoMsg() {
			return "Distribution of node " + getNodeName(m_nTargetNode) + " changed";
		}

		public String getRedoMsg() {
			return "Distribution of node " + getNodeName(m_nTargetNode) + " changed";
		}
	} // class SetDistributionAction

	class RenameAction extends UndoAction {
		/** for serialization */
		static final long serialVersionUID = 1;
		int m_nTargetNode;

		String m_sNewName;

		String m_sOldName;

		RenameAction(int nTargetNode, String sOldName, String sNewName) {
			m_nTargetNode = nTargetNode;
			m_sNewName = sNewName;
			m_sOldName = sOldName;
		} // c'tor

		public void undo() {
			setNodeName(m_nTargetNode, m_sOldName);
		} // undo

		public void redo() {
			setNodeName(m_nTargetNode, m_sNewName);
		} // redo
	} // class RenameAction

	class RenameValueAction extends RenameAction {
		/** for serialization */
		static final long serialVersionUID = 1;
		RenameValueAction(int nTargetNode, String sOldName, String sNewName) {
			super(nTargetNode, sOldName, sNewName);
		} // c'tor

		public void undo() {
			renameNodeValue(m_nTargetNode, m_sNewName, m_sOldName);
		} // undo

		public void redo() {
			renameNodeValue(m_nTargetNode, m_sOldName, m_sNewName);
		} // redo

		public String getUndoMsg() {
			return "Value of node " + getNodeName(m_nTargetNode) + " changed from " + m_sNewName + " to " + m_sOldName;
		}

		public String getRedoMsg() {
			return "Value of node " + getNodeName(m_nTargetNode) + " changed from " + m_sOldName + " to " + m_sNewName;
		}
	} // class RenameValueAction

	class AddValueAction extends UndoAction {
		/** for serialization */
		static final long serialVersionUID = 1;
		int m_nTargetNode;

		String m_sValue;

		AddValueAction(int nTargetNode, String sValue) {
			m_nTargetNode = nTargetNode;
			m_sValue = sValue;
		} // c'tor

		public void undo() {
			try {
				delNodeValue(m_nTargetNode, m_sValue);
			} catch (Exception e) {
				e.printStackTrace();
			}
		} // undo

		public void redo() {
			addNodeValue(m_nTargetNode, m_sValue);
		} // redo

		public String getUndoMsg() {
			return "Value " + m_sValue + " removed from node " + getNodeName(m_nTargetNode);
		}

		public String getRedoMsg() {
			return "Value " + m_sValue + " added to node " + getNodeName(m_nTargetNode);
		}
	} // class AddValueAction

	class DelValueAction extends UndoAction {
		/** for serialization */
		static final long serialVersionUID = 1;
		int m_nTargetNode;

		String m_sValue;

		Estimator[] m_CPT;

		FastVector m_children;

		Estimator[][] m_childAtts;

		Attribute m_att;

		DelValueAction(int nTargetNode, String sValue) {
			try {
				m_nTargetNode = nTargetNode;
				m_sValue = sValue;
				m_att = m_Instances.attribute(nTargetNode);
				SerializedObject so = new SerializedObject(m_Distributions[nTargetNode]);
				m_CPT = (Estimator[]) so.getObject();
				;
				m_children = new FastVector();
				for (int iNode = 0; iNode < getNrOfNodes(); iNode++) {
					if (m_ParentSets[iNode].contains(nTargetNode)) {
						m_children.addElement(iNode);
					}
				}
				m_childAtts = new Estimator[m_children.size()][];
				for (int iChild = 0; iChild < m_children.size(); iChild++) {
					int nChild = (Integer) m_children.elementAt(iChild);
					m_childAtts[iChild] = m_Distributions[nChild];
				}
			} catch (Exception e) {
				e.printStackTrace();
			}
		} // c'tor

		public void undo() {
			try {
				m_Instances.insertAttributeAt(m_att, m_nTargetNode);
				SerializedObject so = new SerializedObject(m_CPT);
				m_Distributions[m_nTargetNode] = (Estimator[]) so.getObject();
				for (int iChild = 0; iChild < m_children.size(); iChild++) {
					int nChild = (Integer) m_children.elementAt(iChild);
					m_Instances.insertAttributeAt(m_att, m_nTargetNode);
					so = new SerializedObject(m_childAtts[iChild]);
					m_Distributions[nChild] = (Estimator[]) so.getObject();
				}
			} catch (Exception e) {
				e.printStackTrace();
			}
		} // undo

		public void redo() {
			try {
				delNodeValue(m_nTargetNode, m_sValue);
			} catch (Exception e) {
				e.printStackTrace();
			}
		} // redo

		public String getUndoMsg() {
			return "Value " + m_sValue + " added to node " + getNodeName(m_nTargetNode);
		}

		public String getRedoMsg() {
			return "Value " + m_sValue + " removed from node " + getNodeName(m_nTargetNode);
		}
	} // class DelValueAction

	class alignAction extends UndoAction {
		/** for serialization */
		static final long serialVersionUID = 1;
		FastVector m_nodes;

		FastVector m_posX;

		FastVector m_posY;

		alignAction(FastVector nodes) {
			m_nodes = new FastVector(nodes.size());
			m_posX = new FastVector(nodes.size());
			m_posY = new FastVector(nodes.size());
			for (int iNode = 0; iNode < nodes.size(); iNode++) {
				int nNode = (Integer) nodes.elementAt(iNode);
				m_nodes.addElement(nNode);
				m_posX.addElement(getPositionX(nNode));
				m_posY.addElement(getPositionY(nNode));
			}
		} // c'tor

		public void undo() {
			try {
				for (int iNode = 0; iNode < m_nodes.size(); iNode++) {
					int nNode = (Integer) m_nodes.elementAt(iNode);
					setPosition(nNode, (Integer) m_posX.elementAt(iNode), (Integer) m_posY.elementAt(iNode));
				}
			} catch (Exception e) {
				e.printStackTrace();
			}
		} // undo
	} // class alignAction

	class alignLeftAction extends alignAction {
		/** for serialization */
		static final long serialVersionUID = 1;
		public alignLeftAction(FastVector nodes) {
			super(nodes);
		} // c'tor

		public void redo() {
			try {
				alignLeft(m_nodes);
			} catch (Exception e) {
				e.printStackTrace();
			}
		} // redo

		public String getUndoMsg() {
			return "Returning " + m_nodes.size() + " from aliging nodes to the left.";
		}

		public String getRedoMsg() {
			return "Aligning " + m_nodes.size() + " nodes to the left.";
		}
	} // class alignLeftAction

	class alignRightAction extends alignAction {
		/** for serialization */
		static final long serialVersionUID = 1;
		public alignRightAction(FastVector nodes) {
			super(nodes);
		} // c'tor

		public void redo() {
			try {
				alignRight(m_nodes);
			} catch (Exception e) {
				e.printStackTrace();
			}
		} // redo

		public String getUndoMsg() {
			return "Returning " + m_nodes.size() + " from aliging nodes to the right.";
		}

		public String getRedoMsg() {
			return "Aligning " + m_nodes.size() + " nodes to the right.";
		}
	} // class alignLeftAction

	class alignTopAction extends alignAction {
		/** for serialization */
		static final long serialVersionUID = 1;
		public alignTopAction(FastVector nodes) {
			super(nodes);
		} // c'tor

		public void redo() {
			try {
				alignTop(m_nodes);
			} catch (Exception e) {
				e.printStackTrace();
			}
		} // redo

		public String getUndoMsg() {
			return "Returning " + m_nodes.size() + " from aliging nodes to the top.";
		}

		public String getRedoMsg() {
			return "Aligning " + m_nodes.size() + " nodes to the top.";
		}
	} // class alignTopAction

	class alignBottomAction extends alignAction {
		/** for serialization */
		static final long serialVersionUID = 1;
		public alignBottomAction(FastVector nodes) {
			super(nodes);
		} // c'tor

		public void redo() {
			try {
				alignBottom(m_nodes);
			} catch (Exception e) {
				e.printStackTrace();
			}
		} // redo

		public String getUndoMsg() {
			return "Returning " + m_nodes.size() + " from aliging nodes to the bottom.";
		}

		public String getRedoMsg() {
			return "Aligning " + m_nodes.size() + " nodes to the bottom.";
		}
	} // class alignBottomAction

	class centerHorizontalAction extends alignAction {
		/** for serialization */
		static final long serialVersionUID = 1;
		public centerHorizontalAction(FastVector nodes) {
			super(nodes);
		} // c'tor

		public void redo() {
			try {
				centerHorizontal(m_nodes);
			} catch (Exception e) {
				e.printStackTrace();
			}
		} // redo

		public String getUndoMsg() {
			return "Returning " + m_nodes.size() + " from centering horizontally.";
		}

		public String getRedoMsg() {
			return "Centering " + m_nodes.size() + " nodes horizontally.";
		}
	} // class centerHorizontalAction

	class centerVerticalAction extends alignAction {
		/** for serialization */
		static final long serialVersionUID = 1;
		public centerVerticalAction(FastVector nodes) {
			super(nodes);
		} // c'tor

		public void redo() {
			try {
				centerVertical(m_nodes);
			} catch (Exception e) {
				e.printStackTrace();
			}
		} // redo

		public String getUndoMsg() {
			return "Returning " + m_nodes.size() + " from centering vertically.";
		}

		public String getRedoMsg() {
			return "Centering " + m_nodes.size() + " nodes vertically.";
		}
	} // class centerVerticalAction

	class spaceHorizontalAction extends alignAction {
		/** for serialization */
		static final long serialVersionUID = 1;
		public spaceHorizontalAction(FastVector nodes) {
			super(nodes);
		} // c'tor

		public void redo() {
			try {
				spaceHorizontal(m_nodes);
			} catch (Exception e) {
				e.printStackTrace();
			}
		} // redo

		public String getUndoMsg() {
			return "Returning " + m_nodes.size() + " from spaceing horizontally.";
		}

		public String getRedoMsg() {
			return "spaceing " + m_nodes.size() + " nodes horizontally.";
		}
	} // class spaceHorizontalAction

	class spaceVerticalAction extends alignAction {
		/** for serialization */
		static final long serialVersionUID = 1;
		public spaceVerticalAction(FastVector nodes) {
			super(nodes);
		} // c'tor

		public void redo() {
			try {
				spaceVertical(m_nodes);
			} catch (Exception e) {
				e.printStackTrace();
			}
		} // redo

		public String getUndoMsg() {
			return "Returning " + m_nodes.size() + " from spaceng vertically.";
		}

		public String getRedoMsg() {
			return "Spaceng " + m_nodes.size() + " nodes vertically.";
		}
	} // class spaceVerticalAction

	class SetPositionAction extends UndoAction {
		/** for serialization */
		static final long serialVersionUID = 1;
		int m_nTargetNode;

		int m_nX;

		int m_nY;

		int m_nX2;

		int m_nY2;

		SetPositionAction(int nTargetNode, int nX, int nY) {
			m_nTargetNode = nTargetNode;
			m_nX2 = nX;
			m_nY2 = nY;
			m_nX = getPositionX(nTargetNode);
			m_nY = getPositionY(nTargetNode);
		} // c'tor

		public void undo() {
			setPosition(m_nTargetNode, m_nX, m_nY);
		} // undo

		public void redo() {
			setPosition(m_nTargetNode, m_nX2, m_nY2);
		} // redo

		public void setUndoPosition(int nX, int nY) {
			m_nX2 = nX;
			m_nY2 = nY;
		} // setPosition
	} // class SetPositionAction

	class SetGroupPositionAction extends UndoAction {
		/** for serialization */
		static final long serialVersionUID = 1;
		FastVector m_nodes;
		int m_dX;
		int m_dY;

		SetGroupPositionAction(FastVector nodes, int dX, int dY) {
			m_nodes = new FastVector(nodes.size());
			for (int iNode = 0; iNode < nodes.size(); iNode++) {
				m_nodes.addElement(nodes.elementAt(iNode));
			}
			m_dX = dX;
			m_dY = dY;
		} // c'tor

		public void undo() {
			for (int iNode = 0; iNode < m_nodes.size(); iNode++) {
				int nNode = (Integer) m_nodes.elementAt(iNode);
				setPosition(nNode, getPositionX(nNode) - m_dX,  getPositionY(nNode) - m_dY);
			}
		} // undo

		public void redo() {
			for (int iNode = 0; iNode < m_nodes.size(); iNode++) {
				int nNode = (Integer) m_nodes.elementAt(iNode);
				setPosition(nNode, getPositionX(nNode) + m_dX,  getPositionY(nNode) + m_dY);
			}
		} // redo
		public void setUndoPosition(int dX, int dY) {
			m_dX += dX;
			m_dY += dY;
		} // setPosition
	} // class SetGroupPositionAction

	class LayoutGraphAction extends UndoAction {
		/** for serialization */
		static final long serialVersionUID = 1;
		FastVector m_nPosX;
		FastVector m_nPosY;
		FastVector m_nPosX2;
		FastVector m_nPosY2;

		LayoutGraphAction(FastVector nPosX, FastVector nPosY) {
			m_nPosX = new FastVector(nPosX.size());
			m_nPosY = new FastVector(nPosX.size());
			m_nPosX2 = new FastVector(nPosX.size());
			m_nPosY2 = new FastVector(nPosX.size());
			for (int iNode = 0; iNode < nPosX.size(); iNode++) {
				m_nPosX.addElement(m_nPositionX.elementAt(iNode));
				m_nPosY.addElement(m_nPositionY.elementAt(iNode));
				m_nPosX2.addElement(nPosX.elementAt(iNode));
				m_nPosY2.addElement(nPosY.elementAt(iNode));
			}
		} // c'tor

		public void undo() {
			for (int iNode = 0; iNode < m_nPosX.size(); iNode++) {
				setPosition(iNode, (Integer) m_nPosX.elementAt(iNode), (Integer) m_nPosY.elementAt(iNode));
			}
		} // undo

		public void redo() {
			for (int iNode = 0; iNode < m_nPosX.size(); iNode++) {
				setPosition(iNode, (Integer) m_nPosX2.elementAt(iNode), (Integer) m_nPosY2.elementAt(iNode));
			}
		} // redo
	} // class LayoutGraphAction

	class PasteAction extends UndoAction {
		/** for serialization */
		static final long serialVersionUID = 1;
		int m_nBase;

		String m_sXML;

		PasteAction(String sXML, int nBase) {
			m_sXML = sXML;
			m_nBase = nBase;
		} // c'tor

		public void undo() {
			try {
				int iNode = getNrOfNodes() - 1;
				while (iNode >= m_nBase) {
					deleteNode(iNode);
					iNode--;
				}
			} catch (Exception e) {
				e.printStackTrace();
			}
		} // undo

		public void redo() {
			try {
				paste(m_sXML, EXECUTE);
			} catch (Exception e) {
				e.printStackTrace();
			}
		} // redo
	} // class PasteAction
	  
	  /**
	   * Returns the revision string.
	   * 
	   * @return		the revision
	   */
	  public String getRevision() {
	    return RevisionUtils.extract("$Revision: 4899 $");
	  }

	/**
	 * @param args
	 */
	public static void main(String[] args) {
	} // main
} // class EditableBayesNet

