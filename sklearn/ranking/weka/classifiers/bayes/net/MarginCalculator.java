package weka.classifiers.bayes.net;

import weka.classifiers.bayes.BayesNet;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;

import java.io.Serializable;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.Vector;


public class MarginCalculator implements Serializable, RevisionHandler {
	  /** for serialization */
	  private static final long serialVersionUID = 650278019241175534L;

	  boolean m_debug = false;
	  public JunctionTreeNode m_root = null;
	JunctionTreeNode [] jtNodes;

	public int getNode(String sNodeName) {
    	int iNode = 0;
    	while (iNode < m_root.m_bayesNet.m_Instances.numAttributes()) {
    		if (m_root.m_bayesNet.m_Instances.attribute(iNode).name().equals(sNodeName)) {
    			return iNode;
    		}
	    	iNode++; 
    	}
    	//throw new Exception("Could not find node [[" + sNodeName + "]]");
    	return -1;
	}
	public String toXMLBIF03() {return m_root.m_bayesNet.toXMLBIF03();}
	
	/**
	 * Calc marginal distributions of nodes in Bayesian network
	 *	 Note that a connected network is assumed. 
	 *	 Unconnected networks may give unexpected results.
	 * @param bayesNet
	 */
	public void calcMargins(BayesNet bayesNet) throws Exception {
		//System.out.println(bayesNet.toString());
		boolean[][] bAdjacencyMatrix = moralize(bayesNet);
		process(bAdjacencyMatrix, bayesNet);
	} // calcMargins

	public void calcFullMargins(BayesNet bayesNet) throws Exception {
		//System.out.println(bayesNet.toString());
		int nNodes = bayesNet.getNrOfNodes();
		boolean[][] bAdjacencyMatrix = new boolean[nNodes][nNodes];
		for (int iNode = 0; iNode < nNodes; iNode++) {
			for (int iNode2 = 0; iNode2 < nNodes; iNode2++) {
				bAdjacencyMatrix[iNode][iNode2] = true;
			}
		}
		process(bAdjacencyMatrix, bayesNet);
	} // calcMargins
	
	
	public void process(boolean[][] bAdjacencyMatrix, BayesNet bayesNet) throws Exception {
		int[] order = getMaxCardOrder(bAdjacencyMatrix);
		bAdjacencyMatrix = fillIn(order, bAdjacencyMatrix);
		order = getMaxCardOrder(bAdjacencyMatrix);
		Set [] cliques = getCliques(order, bAdjacencyMatrix);
		Set [] separators = getSeparators(order, cliques);
		int [] parentCliques = getCliqueTree(order, cliques, separators);
		// report cliques
		int nNodes = bAdjacencyMatrix.length;
		if (m_debug) {
		for (int i = 0; i < nNodes; i++) {
			int iNode = order[i];
			if (cliques[iNode] != null) {
				System.out.print("Clique " + iNode + " (");
				Iterator nodes = cliques[iNode].iterator();
				while (nodes.hasNext()) {
					int iNode2 = (Integer) nodes.next();
					System.out.print(iNode2 + " " + bayesNet.getNodeName(iNode2));
					if (nodes.hasNext()) {
						System.out.print(",");
					}
				}
				System.out.print(") S(");
				nodes = separators[iNode].iterator();
				while (nodes.hasNext()) {
					int iNode2 = (Integer) nodes.next();
					System.out.print(iNode2 + " " + bayesNet.getNodeName(iNode2));
					if (nodes.hasNext()) {
						System.out.print(",");
					}
				}
				System.out.println(") parent clique " + parentCliques[iNode]);
			}		
		}
		}
				
		jtNodes = getJunctionTree(cliques, separators, parentCliques, order, bayesNet);
		m_root = null;
		for (int iNode = 0; iNode < nNodes; iNode++) {
			if (parentCliques[iNode] < 0 && jtNodes[iNode] != null) {
				m_root = jtNodes[iNode];
				break;
			}
		}
		m_Margins = new double[nNodes][];
		initialize(jtNodes, order, cliques, separators, parentCliques);
		
		// sanity check
		for (int i = 0; i < nNodes; i++) {
			int iNode = order[i];
			if (cliques[iNode] != null) {
				if (parentCliques[iNode] == -1 && separators[iNode].size() > 0) {
					throw new Exception("Something wrong in clique tree");
				}
			}
		}
		if (m_debug) {
			//System.out.println(m_root.toString());
		}
	} // process
		
	void initialize(JunctionTreeNode [] jtNodes, int [] order, Set [] cliques, Set [] separators, int [] parentCliques) {
		int nNodes = order.length;
		for (int i = nNodes - 1; i >= 0; i--) {
			int iNode = order[i];
			if (jtNodes[iNode]!=null) {
				jtNodes[iNode].initializeUp();
			}
		}	
		for (int i = 0; i < nNodes; i++) {
			int iNode = order[i];
			if (jtNodes[iNode]!=null) {
				jtNodes[iNode].initializeDown(false);
			}
		}	
	} // initialize
	
	JunctionTreeNode [] getJunctionTree(Set [] cliques, Set [] separators, int [] parentCliques, int [] order, BayesNet bayesNet) {
		int nNodes = order.length;
		JunctionTreeNode root = null;
		JunctionTreeNode [] jtns = new JunctionTreeNode[nNodes]; 
		boolean [] bDone = new boolean[nNodes];
		// create junction tree nodes
		for (int i = 0; i < nNodes; i++) {
			int iNode = order[i];
			if (cliques[iNode] != null) {
				jtns[iNode] = new JunctionTreeNode(cliques[iNode], bayesNet, bDone);
			}
		}
		// create junction tree separators
		for (int i = 0; i < nNodes; i++) {
			int iNode = order[i];
			if (cliques[iNode] != null) {
				JunctionTreeNode parent = null;
				if (parentCliques[iNode] > 0) {
					parent = jtns[parentCliques[iNode]];
					JunctionTreeSeparator jts = new JunctionTreeSeparator(separators[iNode], bayesNet, jtns[iNode], parent);
					jtns[iNode].setParentSeparator(jts);
					jtns[parentCliques[iNode]].addChildClique(jtns[iNode]);
				} else {
					root = jtns[iNode];	
				}
			}
		}
		return jtns;
	} // getJunctionTree
	
	public class JunctionTreeSeparator implements Serializable, RevisionHandler {
	  
		  private static final long serialVersionUID = 6502780192411755343L;
		int [] m_nNodes;
		int m_nCardinality;
		double [] m_fiParent;
		double [] m_fiChild;
		JunctionTreeNode m_parentNode;
		JunctionTreeNode m_childNode;
		BayesNet m_bayesNet;
		
		JunctionTreeSeparator(Set separator, BayesNet bayesNet, JunctionTreeNode childNode, JunctionTreeNode parentNode) {
			//////////////////////
			// initialize node set
			m_nNodes = new int[separator.size()];
			int iPos = 0;
			m_nCardinality = 1;
			for(Iterator nodes = separator.iterator(); nodes.hasNext();) {
				int iNode = (Integer) nodes.next();
				m_nNodes[iPos++] = iNode;
				m_nCardinality *= bayesNet.getCardinality(iNode);
			}
			m_parentNode = parentNode;
			m_childNode = childNode;
			m_bayesNet = bayesNet;
		} // c'tor
		
		/** marginalize junciontTreeNode node over all nodes outside the separator set
		 * of the parent clique
		 *
		 */
		public void updateFromParent() {
			double [] fis = update(m_parentNode); 
			if (fis == null) {
				m_fiParent = null;
			} else {
				m_fiParent = fis;
				// normalize
				double sum = 0;
				for (int iPos = 0; iPos < m_nCardinality; iPos++) {
					sum += m_fiParent[iPos];
				}
				for (int iPos = 0; iPos < m_nCardinality; iPos++) {
					m_fiParent[iPos] /= sum;
				}
			}
		} // updateFromParent

		/** marginalize junciontTreeNode node over all nodes outside the separator set
		 * of the child clique
		 *
		 */
		public void updateFromChild() {
			double [] fis = update(m_childNode); 
			if (fis == null) {
				m_fiChild = null;
			} else {
				m_fiChild = fis;
				// normalize
				double sum = 0;
				for (int iPos = 0; iPos < m_nCardinality; iPos++) {
					sum += m_fiChild[iPos];
				}
				for (int iPos = 0; iPos < m_nCardinality; iPos++) {
					m_fiChild[iPos] /= sum;
				}
			}
		} // updateFromChild
		
		/** marginalize junciontTreeNode node over all nodes outside the separator set
		 * 
		 * @param node one of the neighboring junciont tree nodes of this separator
		 */
		public double [] update(JunctionTreeNode node) {
			if (node.m_P == null) {
				return null;
			}
			double [] fi = new double[m_nCardinality];

			int [] values = new int[node.m_nNodes.length];
			int [] order = new int[m_bayesNet.getNrOfNodes()];
			for (int iNode = 0; iNode < node.m_nNodes.length; iNode++) {
				order[node.m_nNodes[iNode]] = iNode;
			}
			// fill in the values
			for (int iPos = 0; iPos < node.m_nCardinality; iPos++) {
				int iNodeCPT = getCPT(node.m_nNodes, node.m_nNodes.length, values, order, m_bayesNet);
				int iSepCPT =  getCPT(m_nNodes, m_nNodes.length, values, order, m_bayesNet);
				fi[iSepCPT] += node.m_P[iNodeCPT];
				// update values
				int i = 0;
				values[i]++;
				while (i < node.m_nNodes.length && values[i] == m_bayesNet.getCardinality(node.m_nNodes[i])) {
					values[i] = 0;
					i++;
					if (i < node.m_nNodes.length) {
						values[i]++;
					}
				}
			}
			return fi;
		} // update
		  
		/**
		 * Returns the revision string.
		 * 
		 * @return		the revision
		 */
		public String getRevision() {
		  return RevisionUtils.extract("$Revision: 4899 $");
		}

	} // class JunctionTreeSeparator

	public class JunctionTreeNode implements Serializable, RevisionHandler {
	  
		  private static final long serialVersionUID = 650278019241175536L;
		/** reference Bayes net for information about variables like name, cardinality, etc.
		 * but not for relations between nodes **/
		BayesNet m_bayesNet;
		/** nodes of the Bayes net in this junction node **/
		public int [] m_nNodes;
		/** cardinality of the instances of variables in this junction node **/
		int m_nCardinality;
		/** potentials for first network **/
		double [] m_fi;

		/** distribution over this junction node according to first Bayes network **/
		double [] m_P;


		double [][] m_MarginalP;		

		
		JunctionTreeSeparator m_parentSeparator;
		public void setParentSeparator(JunctionTreeSeparator parentSeparator) {m_parentSeparator = parentSeparator;}
		public Vector m_children;
		public void addChildClique(JunctionTreeNode child) {m_children.add(child);}

		public void initializeUp() {
			m_P = new double[m_nCardinality];
			for (int iPos = 0; iPos < m_nCardinality; iPos++) {
				m_P[iPos] = m_fi[iPos];
			}
			int [] values = new int[m_nNodes.length];
			int [] order = new int[m_bayesNet.getNrOfNodes()];
			for (int iNode = 0; iNode < m_nNodes.length; iNode++) {
				order[m_nNodes[iNode]] = iNode;
			}
			for (Iterator child = m_children.iterator(); child.hasNext(); ) {
				JunctionTreeNode childNode = (JunctionTreeNode) child.next();
				JunctionTreeSeparator separator = childNode.m_parentSeparator;
			// Update the values
			for (int iPos = 0; iPos < m_nCardinality; iPos++) {
				int iSepCPT = getCPT(separator.m_nNodes, separator.m_nNodes.length, values, order, m_bayesNet);
				int iNodeCPT =  getCPT(m_nNodes, m_nNodes.length, values, order, m_bayesNet);
					m_P[iNodeCPT] *= separator.m_fiChild[iSepCPT];					
				// update values
				int i = 0;
				values[i]++;
				while (i < m_nNodes.length && values[i] == m_bayesNet.getCardinality(m_nNodes[i])) {
					values[i] = 0;
					i++;
					if (i < m_nNodes.length) {
						values[i]++;
					}
				}
			}
			}
			// normalize
			double sum = 0;
			for (int iPos = 0; iPos < m_nCardinality; iPos++) {
				sum += m_P[iPos];
			}
			for (int iPos = 0; iPos < m_nCardinality; iPos++) {
				m_P[iPos] /= sum;
			}

			if (m_parentSeparator != null) { // not a root node
				m_parentSeparator.updateFromChild();
			}
		} // initializeUp

		public void initializeDown(boolean recursively) {
			if (m_parentSeparator == null) { // a root node
				calcMarginalProbabilities();
			} else {
			m_parentSeparator.updateFromParent();
				int [] values = new int[m_nNodes.length];
				int [] order = new int[m_bayesNet.getNrOfNodes()];
				for (int iNode = 0; iNode < m_nNodes.length; iNode++) {
					order[m_nNodes[iNode]] = iNode;
				}

				
				// Update the values
				for (int iPos = 0; iPos < m_nCardinality; iPos++) {
					int iSepCPT = getCPT(m_parentSeparator.m_nNodes, m_parentSeparator.m_nNodes.length, values, order, m_bayesNet);
					int iNodeCPT =  getCPT(m_nNodes, m_nNodes.length, values, order, m_bayesNet);
					if ( m_parentSeparator.m_fiChild[iSepCPT] > 0) {
						m_P[iNodeCPT] *= m_parentSeparator.m_fiParent[iSepCPT] / m_parentSeparator.m_fiChild[iSepCPT];
					} else {
						m_P[iNodeCPT] = 0;
					}
					// update values
					int i = 0;
					values[i]++;
					while (i < m_nNodes.length && values[i] == m_bayesNet.getCardinality(m_nNodes[i])) {
						values[i] = 0;
						i++;
						if (i < m_nNodes.length) {
							values[i]++;
						}
					}
				}
				// normalize
				double sum = 0;
				for (int iPos = 0; iPos < m_nCardinality; iPos++) {
					sum += m_P[iPos];
				}
				for (int iPos = 0; iPos < m_nCardinality; iPos++) {
					m_P[iPos] /= sum;
				}
				m_parentSeparator.updateFromChild();
				calcMarginalProbabilities();
			}
			if (recursively) {
				for (Iterator child = m_children.iterator(); child.hasNext(); ) {
					JunctionTreeNode childNode = (JunctionTreeNode) child.next();
					childNode.initializeDown(true);
				}			
			}
		} // initializeDown
		
		
		/** calculate marginal probabilities for the individual nodes in the clique.
		 * Store results in m_MarginalP 
		 */
		void calcMarginalProbabilities() {			
			// calculate marginal probabilities
			int [] values = new int[m_nNodes.length];
			int [] order = new int[m_bayesNet.getNrOfNodes()];
			m_MarginalP = new double[m_nNodes.length][];
			for (int iNode = 0; iNode < m_nNodes.length; iNode++) {
				order[m_nNodes[iNode]] = iNode;
				m_MarginalP[iNode]=new double[m_bayesNet.getCardinality(m_nNodes[iNode])];
			}
			for (int iPos = 0; iPos < m_nCardinality; iPos++) {
				int iNodeCPT =  getCPT(m_nNodes, m_nNodes.length, values, order, m_bayesNet);
				for (int iNode = 0; iNode < m_nNodes.length; iNode++) {
					m_MarginalP[iNode][values[iNode]] += m_P[iNodeCPT];
				}
				// update values
				int i = 0;
				values[i]++;
				while (i < m_nNodes.length && values[i] == m_bayesNet.getCardinality(m_nNodes[i])) {
					values[i] = 0;
					i++;
					if (i < m_nNodes.length) {
						values[i]++;
					}
				}
			}
			
			for (int iNode = 0; iNode < m_nNodes.length; iNode++) {
				m_Margins[m_nNodes[iNode]] = m_MarginalP[iNode]; 
			}
		} // calcMarginalProbabilities
		
		public String toString() {
			StringBuffer buf = new StringBuffer();
			for (int iNode = 0; iNode < m_nNodes.length; iNode++) {
				buf.append(m_bayesNet.getNodeName(m_nNodes[iNode]) + ": ");
				for (int iValue = 0; iValue < m_MarginalP[iNode].length; iValue++) {
					buf.append(m_MarginalP[iNode][iValue] + " ");
				}
				buf.append('\n');
			}
			for (Iterator child = m_children.iterator(); child.hasNext(); ) {
				JunctionTreeNode childNode = (JunctionTreeNode) child.next();
				buf.append("----------------\n");
				buf.append(childNode.toString());
			}			
			return buf.toString();
		} // toString
		
		void calculatePotentials(BayesNet bayesNet, Set clique, boolean [] bDone) {
			m_fi = new double[m_nCardinality];
			
			int [] values = new int[m_nNodes.length];
			int [] order = new int[bayesNet.getNrOfNodes()];
			for (int iNode = 0; iNode < m_nNodes.length; iNode++) {
				order[m_nNodes[iNode]] = iNode;
			}
			// find conditional probabilities that need to be taken in account
			boolean [] bIsContained = new boolean[m_nNodes.length];
			for (int iNode = 0; iNode < m_nNodes.length; iNode++) {
				int nNode = m_nNodes[iNode];
				bIsContained[iNode] = !bDone[nNode];
				for (int iParent = 0; iParent < bayesNet.getNrOfParents(nNode); iParent++) {
					int nParent = bayesNet.getParent(nNode, iParent);
					if (!clique.contains(nParent)) {
						bIsContained[iNode] = false;
					}
				}
				if (bIsContained[iNode]) {
					bDone[nNode] = true;
					if (m_debug) {
						System.out.println("adding node " +nNode);
					}
				}
			}			

			// fill in the values
			for (int iPos = 0; iPos < m_nCardinality; iPos++) {
				int iCPT = getCPT(m_nNodes, m_nNodes.length, values, order, bayesNet);
				m_fi[iCPT] = 1.0;
				for (int iNode = 0; iNode < m_nNodes.length; iNode++) {
					if (bIsContained[iNode]) {
						int nNode = m_nNodes[iNode];
						int [] nNodes = bayesNet.getParentSet(nNode).getParents();
						int iCPT2 = getCPT(nNodes, bayesNet.getNrOfParents(nNode), values, order, bayesNet);
						double f = bayesNet.getDistributions()[nNode][iCPT2].getProbability(values[iNode]);
						m_fi[iCPT] *= f;
					}
				}
				
				// update values
				int i = 0;
				values[i]++;
				while (i < m_nNodes.length && values[i] == bayesNet.getCardinality(m_nNodes[i])) {
					values[i] = 0;
					i++;
					if (i < m_nNodes.length) {
						values[i]++;
					}
				}
			}
		} // calculatePotentials

		JunctionTreeNode(Set clique, BayesNet bayesNet, boolean [] bDone) {
			m_bayesNet = bayesNet;
			m_children = new Vector();
			//////////////////////
			// initialize node set
			m_nNodes = new int[clique.size()];
			int iPos = 0;
			m_nCardinality = 1;
			for(Iterator nodes = clique.iterator(); nodes.hasNext();) {
				int iNode = (Integer) nodes.next();
				m_nNodes[iPos++] = iNode;
				m_nCardinality *= bayesNet.getCardinality(iNode);
			}
			////////////////////////////////
			// initialize potential function
			calculatePotentials(bayesNet, clique, bDone);
       } // JunctionTreeNode c'tor

		/* check whether this junciton tree node contains node nNode
		 * 
		 */
		boolean contains(int nNode) {
			for (int iNode = 0; iNode < m_nNodes.length; iNode++) {
				if (m_nNodes[iNode]== nNode){
					return true;
				}
			}
			return false;
		} // contains
		
		public void setEvidence(int nNode, int iValue) throws Exception {
			int [] values = new int[m_nNodes.length];
			int [] order = new int[m_bayesNet.getNrOfNodes()];

			int nNodeIdx = -1;
			for (int iNode = 0; iNode < m_nNodes.length; iNode++) {
				order[m_nNodes[iNode]] = iNode;
				if (m_nNodes[iNode] == nNode) {
					nNodeIdx = iNode;
				}
			}
			if (nNodeIdx < 0) {
				throw new Exception("setEvidence: Node " + nNode + " not found in this clique");
			}
			for (int iPos = 0; iPos < m_nCardinality; iPos++) {
				if (values[nNodeIdx] != iValue) {
					int iNodeCPT =  getCPT(m_nNodes, m_nNodes.length, values, order, m_bayesNet);
					m_P[iNodeCPT] = 0;
				}
				// update values
				int i = 0;
				values[i]++;
				while (i < m_nNodes.length && values[i] == m_bayesNet.getCardinality(m_nNodes[i])) {
					values[i] = 0;
					i++;
					if (i < m_nNodes.length) {
						values[i]++;
					}
				}
			}		
			// normalize
			double sum = 0;
			for (int iPos = 0; iPos < m_nCardinality; iPos++) {
				sum += m_P[iPos];
			}
			for (int iPos = 0; iPos < m_nCardinality; iPos++) {
				m_P[iPos] /= sum;
			}
			calcMarginalProbabilities();
			updateEvidence(this);
		} // setEvidence

		void updateEvidence(JunctionTreeNode source) {
			if (source != this) {
				int [] values = new int[m_nNodes.length];
				int [] order = new int[m_bayesNet.getNrOfNodes()];
				for (int iNode = 0; iNode < m_nNodes.length; iNode++) {
					order[m_nNodes[iNode]] = iNode;
				}
				int [] nChildNodes = source.m_parentSeparator.m_nNodes;
				int nNumChildNodes = nChildNodes.length; 
				for (int iPos = 0; iPos < m_nCardinality; iPos++) {
					int iNodeCPT =  getCPT(m_nNodes, m_nNodes.length, values, order, m_bayesNet);
					int iChildCPT =  getCPT(nChildNodes, nNumChildNodes, values, order, m_bayesNet);
					if (source.m_parentSeparator.m_fiParent[iChildCPT] != 0) {
						m_P[iNodeCPT] *= source.m_parentSeparator.m_fiChild[iChildCPT]/source.m_parentSeparator.m_fiParent[iChildCPT];
					} else {
						m_P[iNodeCPT] = 0;
					}
					// update values
					int i = 0;
					values[i]++;
					while (i < m_nNodes.length && values[i] == m_bayesNet.getCardinality(m_nNodes[i])) {
						values[i] = 0;
						i++;
						if (i < m_nNodes.length) {
							values[i]++;
						}
					}
				}		
				// normalize
				double sum = 0;
				for (int iPos = 0; iPos < m_nCardinality; iPos++) {
					sum += m_P[iPos];
				}
				for (int iPos = 0; iPos < m_nCardinality; iPos++) {
					m_P[iPos] /= sum;
				}
				calcMarginalProbabilities();
			}
			for (Iterator child = m_children.iterator(); child.hasNext(); ) {
				JunctionTreeNode childNode = (JunctionTreeNode) child.next();
				if (childNode != source) {
					childNode.initializeDown(true);
				}
			}			
			if (m_parentSeparator != null) {
				m_parentSeparator.updateFromChild();
				m_parentSeparator.m_parentNode.updateEvidence(this);
				m_parentSeparator.updateFromParent();
			}
		} // updateEvidence

		/**
		 * Returns the revision string.
		 * 
		 * @return		the revision
		 */
		public String getRevision() {
		  return RevisionUtils.extract("$Revision: 4899 $");
		}
		
	} // class JunctionTreeNode

	int getCPT(int [] nodeSet, int nNodes, int[] values, int[] order, BayesNet bayesNet) {
		int iCPTnew = 0;
		for (int iNode = 0; iNode < nNodes; iNode++) {
			int nNode = nodeSet[iNode];
			iCPTnew = iCPTnew * bayesNet.getCardinality(nNode);
			iCPTnew += values[order[nNode]];
		}
		return iCPTnew;
	} // getCPT

	int [] getCliqueTree(int [] order, Set [] cliques, Set [] separators) {
		int nNodes = order.length;
		int [] parentCliques = new int[nNodes];
		//for (int i = nNodes - 1; i >= 0; i--) {
		for (int i = 0; i < nNodes; i++) {
			int iNode = order[i];
			parentCliques[iNode] = -1;
			if (cliques[iNode] != null && separators[iNode].size() > 0) {
				//for (int j = nNodes - 1; j > i; j--) {
				for (int j = 0; j < nNodes; j++) {
					int iNode2 = order[j];
					if (iNode!= iNode2 && cliques[iNode2] != null && cliques[iNode2].containsAll(separators[iNode])) {
						parentCliques[iNode] = iNode2;
						j = i;
						j = 0;
						j = nNodes;
					}
				}
				
			}
		}
		return parentCliques;
	} // getCliqueTree
	
	/** calculate separator sets in clique tree
	 * 
	 * @param order: maximum cardinality ordering of the graph
	 * @param cliques: set of cliques
	 * @return set of separator sets
	 */
	Set [] getSeparators(int [] order, Set [] cliques) {
		int nNodes = order.length;
		Set [] separators = new HashSet[nNodes];
		Set processedNodes = new HashSet(); 
		//for (int i = nNodes - 1; i >= 0; i--) {
		for (int i = 0; i < nNodes; i++) {
			int iNode = order[i];
			if (cliques[iNode] != null) {
				Set separator = new HashSet();
				separator.addAll(cliques[iNode]);
				separator.retainAll(processedNodes);
				separators[iNode] = separator;
				processedNodes.addAll(cliques[iNode]);
			}
		}
		return separators;
	} // getSeparators
	
	/**
	 * get cliques in a decomposable graph represented by an adjacency matrix
	 * 
	 * @param order: maximum cardinality ordering of the graph
	 * @param bAdjacencyMatrix: decomposable graph
	 * @return set of cliques
	 */
	Set [] getCliques(int[] order, boolean[][] bAdjacencyMatrix) throws Exception {
		int nNodes = bAdjacencyMatrix.length;
		Set [] cliques = new HashSet[nNodes];
		//int[] inverseOrder = new int[nNodes];
		//for (int iNode = 0; iNode < nNodes; iNode++) {
			//inverseOrder[order[iNode]] = iNode;
		//}
		// consult nodes in reverse order
		for (int i = nNodes - 1; i >= 0; i--) {
			int iNode = order[i];
			if (iNode == 22) {
				int h = 3;
				h ++;
			}
			Set clique = new HashSet();
			clique.add(iNode);
			for (int j = 0; j < i; j++) {
				int iNode2 = order[j];
				if (bAdjacencyMatrix[iNode][iNode2]) {
					clique.add(iNode2);
				}
			}
			
			//for (int iNode2 = 0; iNode2 < nNodes; iNode2++) {
				//if (bAdjacencyMatrix[iNode][iNode2] && inverseOrder[iNode2] < inverseOrder[iNode]) {
					//clique.add(iNode2);
				//}
			//}
			cliques[iNode] = clique;
		}
		for (int iNode = 0; iNode < nNodes; iNode++) {
			for (int iNode2 = 0; iNode2 < nNodes; iNode2++) {
				if (iNode != iNode2 && cliques[iNode]!= null && cliques[iNode2]!= null && cliques[iNode].containsAll(cliques[iNode2])) {
					cliques[iNode2] = null;
				}
			}
		}		
		// sanity check
		if (m_debug) {
		int [] nNodeSet = new int[nNodes];
		for (int iNode = 0; iNode < nNodes; iNode++) {
			if (cliques[iNode] != null) {
				Iterator it = cliques[iNode].iterator();
				int k = 0;
				while (it.hasNext()) {
					nNodeSet[k++] = (Integer) it.next();
				}
				for (int i = 0; i < cliques[iNode].size(); i++) {
					for (int j = 0; j < cliques[iNode].size(); j++) {
						if (i!=j && !bAdjacencyMatrix[nNodeSet[i]][nNodeSet[j]]) {
							throw new Exception("Non clique" + i + " " + j);
						}
					}
				}
			}
		}
		}
		return cliques;
	} // getCliques

	/**
	 * moralize DAG and calculate
	 * adjacency matrix representation for a Bayes Network, effecively
	 * converting the directed acyclic graph to an undirected graph.
	 * 
	 * @param bayesNet
	 *            Bayes Network to process
	 * @return adjacencies in boolean matrix format
	 */
	public boolean[][] moralize(BayesNet bayesNet) {
		int nNodes = bayesNet.getNrOfNodes();
		boolean[][] bAdjacencyMatrix = new boolean[nNodes][nNodes];
		for (int iNode = 0; iNode < nNodes; iNode++) {
			ParentSet parents = bayesNet.getParentSets()[iNode];
			moralizeNode(parents, iNode, bAdjacencyMatrix);
		}
		return bAdjacencyMatrix;
	} // moralize

	private void moralizeNode(ParentSet parents, int iNode, boolean[][] bAdjacencyMatrix) {
		for (int iParent = 0; iParent < parents.getNrOfParents(); iParent++) {
			int nParent = parents.getParent(iParent);
			if ( m_debug && !bAdjacencyMatrix[iNode][nParent])
				System.out.println("Insert " + iNode + "--" + nParent);
			bAdjacencyMatrix[iNode][nParent] = true;
			bAdjacencyMatrix[nParent][iNode] = true;
			for (int iParent2 = iParent + 1; iParent2 < parents.getNrOfParents(); iParent2++) {
				int nParent2 = parents.getParent(iParent2);
				if (m_debug && !bAdjacencyMatrix[nParent2][nParent])
					System.out.println("Mary " + nParent + "--" + nParent2);
				bAdjacencyMatrix[nParent2][nParent] = true;
				bAdjacencyMatrix[nParent][nParent2] = true;
			}
		}	
	} // moralizeNode
	
	/**
	 * Apply Tarjan and Yannakakis (1984) fill in algorithm for graph
	 * triangulation. In reverse order, insert edges between any non-adjacent
	 * neighbors that are lower numbered in the ordering.
	 * 
	 * Side effect: input matrix is used as output
	 * 
	 * @param order
	 *            node ordering
	 * @param bAdjacencyMatrix
	 *            boolean matrix representing the graph
	 * @return boolean matrix representing the graph with fill ins
	 */
	public boolean[][] fillIn(int[] order, boolean[][] bAdjacencyMatrix) {
		int nNodes = bAdjacencyMatrix.length;
		int[] inverseOrder = new int[nNodes];
		for (int iNode = 0; iNode < nNodes; iNode++) {
			inverseOrder[order[iNode]] = iNode;
		}
		// consult nodes in reverse order
		for (int i = nNodes - 1; i >= 0; i--) {
			int iNode = order[i];
			// find pairs of neighbors with lower order
			for (int j = 0; j < i; j++) {
				int iNode2 = order[j];
				if (bAdjacencyMatrix[iNode][iNode2]) {
					for (int k = j+1; k < i; k++) {
						int iNode3 = order[k];
						if (bAdjacencyMatrix[iNode][iNode3]) {
							// fill in
							if (m_debug && (!bAdjacencyMatrix[iNode2][iNode3] || !bAdjacencyMatrix[iNode3][iNode2]) )
								System.out.println("Fill in " + iNode2 + "--" + iNode3);
							bAdjacencyMatrix[iNode2][iNode3] = true;
							bAdjacencyMatrix[iNode3][iNode2] = true;
						}
					}
				}
			}
		}
		return bAdjacencyMatrix;
	} // fillIn

	/**
	 * calculate maximum cardinality ordering; start with first node add node
	 * that has most neighbors already ordered till all nodes are in the
	 * ordering
	 * 
	 * This implementation does not assume the graph is connected
	 * 
	 * @param bAdjacencyMatrix:
	 *            n by n matrix with adjacencies in graph of n nodes
	 * @return maximum cardinality ordering
	 */
	int[] getMaxCardOrder(boolean[][] bAdjacencyMatrix) {
		int nNodes = bAdjacencyMatrix.length;
		int[] order = new int[nNodes];
		if (nNodes==0) {return order;}
		boolean[] bDone = new boolean[nNodes];
		// start with node 0
		order[0] = 0;
		bDone[0] = true;
		// order remaining nodes
		for (int iNode = 1; iNode < nNodes; iNode++) {
			int nMaxCard = -1;
			int iBestNode = -1;
			// find node with higest cardinality of previously ordered nodes
			for (int iNode2 = 0; iNode2 < nNodes; iNode2++) {
				if (!bDone[iNode2]) {
					int nCard = 0;
					// calculate cardinality for node iNode2
					for (int iNode3 = 0; iNode3 < nNodes; iNode3++) {
						if (bAdjacencyMatrix[iNode2][iNode3] && bDone[iNode3]) {
							nCard++;
						}
					}
					if (nCard > nMaxCard) {
						nMaxCard = nCard;
						iBestNode = iNode2;
					}
				}
			}
			order[iNode] = iBestNode;
			bDone[iBestNode] = true;
		}
		return order;
	} // getMaxCardOrder

	public void setEvidence(int nNode, int iValue) throws Exception {
		if (m_root == null) {
			throw new Exception("Junction tree not initialize yet");
		}
		int iJtNode = 0;
		while (iJtNode < jtNodes.length && (jtNodes[iJtNode] == null ||!jtNodes[iJtNode].contains(nNode))) {
			iJtNode++;
		}
		if (jtNodes.length == iJtNode) {
			throw new Exception("Could not find node " + nNode + " in junction tree");
		}
		jtNodes[iJtNode].setEvidence(nNode, iValue);
	} // setEvidence
	
	public String toString() {
		return m_root.toString();
	} // toString

	double [][] m_Margins;
	public double [] getMargin(int iNode) {
		return m_Margins[iNode];
	} // getMargin

	/**
	 * Returns the revision string.
	 * 
	 * @return		the revision
	 */
	public String getRevision() {
	  return RevisionUtils.extract("$Revision: 4899 $");
	}
	
	public static void main(String[] args) {
		try {
			BIFReader bayesNet = new BIFReader();
			bayesNet.processFile(args[0]);

			MarginCalculator dc = new MarginCalculator();
			dc.calcMargins(bayesNet);
			int iNode = 2;
			int iValue = 0;
			int iNode2 = 4;
			int iValue2 = 0;
			dc.setEvidence(iNode, iValue);
			dc.setEvidence(iNode2, iValue2);
			System.out.print(dc.toString());


			dc.calcFullMargins(bayesNet);
			dc.setEvidence(iNode, iValue);
			dc.setEvidence(iNode2, iValue2);
			System.out.println("==============");
			System.out.print(dc.toString());
			
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	} // main

} // class MarginCalculator
