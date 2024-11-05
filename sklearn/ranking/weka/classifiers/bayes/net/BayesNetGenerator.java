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
 * BayesNet.java
 * Copyright (C) 2003 University of Waikato, Hamilton, New Zealand
 * 
 */

package weka.classifiers.bayes.net;

import weka.classifiers.bayes.net.estimate.DiscreteEstimatorBayes;
import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.DenseInstance;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.estimators.Estimator;

import java.util.Enumeration;
import java.util.Random;
import java.util.Vector;

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
 * <pre> -B
 *  Generate network (instead of instances)
 * </pre>
 * 
 * <pre> -N &lt;integer&gt;
 *  Nr of nodes
 * </pre>
 * 
 * <pre> -A &lt;integer&gt;
 *  Nr of arcs
 * </pre>
 * 
 * <pre> -M &lt;integer&gt;
 *  Nr of instances
 * </pre>
 * 
 * <pre> -C &lt;integer&gt;
 *  Cardinality of the variables
 * </pre>
 * 
 * <pre> -S &lt;integer&gt;
 *  Seed for random number generator
 * </pre>
 * 
 * <pre> -F &lt;file&gt;
 *  The BIF file to obtain the structure from.
 * </pre>
 * 
 <!-- options-end -->
 *
 * @author Remco Bouckaert (rrb@xm.co.nz)
 * @version $Revision: 5987 $
 */
public class BayesNetGenerator extends EditableBayesNet {
    /** the seed value */
    int m_nSeed = 1;
    
    /** the random number generator */
    Random random;
    
    /** for serialization */
    static final long serialVersionUID = -7462571170596157720L;

	/**
	 * Constructor for BayesNetGenerator.
	 */
	public BayesNetGenerator() {
		super();
	} // c'tor

	/** 
	 * Generate random connected Bayesian network with discrete nodes
	 * having all the same cardinality.
	 * 
	 * @throws Exception if something goes wrong
	 */
	public void generateRandomNetwork () throws Exception {
		if (m_otherBayesNet == null) {
			// generate from scratch
			Init(m_nNrOfNodes, m_nCardinality);
			generateRandomNetworkStructure(m_nNrOfNodes, m_nNrOfArcs);
			generateRandomDistributions(m_nNrOfNodes, m_nCardinality);
		} else {
			// read from file, just copy parent sets and distributions
			m_nNrOfNodes = m_otherBayesNet.getNrOfNodes();
			m_ParentSets = m_otherBayesNet.getParentSets();
			m_Distributions = m_otherBayesNet.getDistributions();


			random = new Random(m_nSeed);
			// initialize m_Instances
			FastVector attInfo = new FastVector(m_nNrOfNodes);
			// generate value strings

			for (int iNode = 0; iNode < m_nNrOfNodes; iNode++) {
				int nValues = m_otherBayesNet.getCardinality(iNode);
				FastVector nomStrings = new FastVector(nValues + 1);
				for (int iValue = 0; iValue < nValues; iValue++) {
					nomStrings.addElement(m_otherBayesNet.getNodeValue(iNode, iValue));
				}
				Attribute att = new Attribute(m_otherBayesNet.getNodeName(iNode), nomStrings);
				attInfo.addElement(att);
			}

			m_Instances = new Instances(m_otherBayesNet.getName(), attInfo, 100);
			m_Instances.setClassIndex(m_nNrOfNodes - 1);
		}
	} // GenerateRandomNetwork

	/** 
	 * Init defines a minimal Bayes net with no arcs
	 * @param nNodes number of nodes in the Bayes net 
	 * @param nValues number of values each of the nodes can take
	 * @throws Exception if something goes wrong
	 */
	public void Init(int nNodes, int nValues) throws Exception {
		random = new Random(m_nSeed);
		// initialize structure
		FastVector attInfo = new FastVector(nNodes);
		// generate value strings
        FastVector nomStrings = new FastVector(nValues + 1);
        for (int iValue = 0; iValue < nValues; iValue++) {
			nomStrings.addElement("Value" + (iValue + 1));
        }

		for (int iNode = 0; iNode < nNodes; iNode++) {
			Attribute att = new Attribute("Node" + (iNode + 1), nomStrings);
			attInfo.addElement(att);
		}
 		m_Instances = new Instances("RandomNet", attInfo, 100);
 		m_Instances.setClassIndex(nNodes - 1);
		setUseADTree(false);
// 		m_bInitAsNaiveBayes = false;
// 		m_bMarkovBlanketClassifier = false;
		initStructure();
		
		// initialize conditional distribution tables
		m_Distributions = new Estimator[nNodes][1];
		for (int iNode = 0; iNode < nNodes; iNode++) {
			m_Distributions[iNode][0] = 
			  new DiscreteEstimatorBayes(nValues, getEstimator().getAlpha());
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
	} // DefineNodes

	/** 
	 * GenerateRandomNetworkStructure generate random connected Bayesian network 
	 * @param nNodes number of nodes in the Bayes net to generate
	 * @param nArcs number of arcs to generate. Must be between nNodes - 1 and nNodes * (nNodes-1) / 2
	 * @throws Exception if number of arcs is incorrect
	 */
	public void generateRandomNetworkStructure(int nNodes, int nArcs) 
		throws Exception
	{
		if (nArcs < nNodes - 1) {
			throw new Exception("Number of arcs should be at least (nNodes - 1) = " + (nNodes - 1) + " instead of " + nArcs);
		}
		if (nArcs > nNodes * (nNodes - 1) / 2) {
			throw new Exception("Number of arcs should be at most nNodes * (nNodes - 1) / 2 = "+ (nNodes * (nNodes - 1) / 2) + " instead of " + nArcs);
		}
		if (nArcs == 0) {return;} // deal with  patalogical case for nNodes = 1

	    // first generate tree connecting all nodes
	    generateTree(nNodes);
	    // The tree contains nNodes - 1 arcs, so there are 
	    // nArcs - (nNodes-1) to add at random.
	    // All arcs point from lower to higher ordered nodes
	    // so that acyclicity is ensured.
	    for (int iArc = nNodes - 1; iArc < nArcs; iArc++) {
	    	boolean bDone = false;
	    	while (!bDone) {
				int nNode1 = random.nextInt(nNodes);
				int nNode2 = random.nextInt(nNodes);
				if (nNode1 == nNode2) {nNode2 = (nNode1 + 1) % nNodes;}
				if (nNode2 < nNode1) {int h = nNode1; nNode1 = nNode2; nNode2 = h;}
				if (!m_ParentSets[nNode2].contains(nNode1)) {
					m_ParentSets[nNode2].addParent(nNode1, m_Instances);
					bDone = true;
				}
	    	}
	    }

	} // GenerateRandomNetworkStructure
	
	/** 
	 * GenerateTree creates a tree-like network structure (actually a
	 * forest) by starting with a randomly selected pair of nodes, add 
	 * an arc between. Then keep on selecting one of the connected nodes 
	 * and one of the unconnected ones and add an arrow between them, 
	 * till all nodes are connected.
	 * @param nNodes number of nodes in the Bayes net to generate
	 */
	void generateTree(int nNodes) {
        boolean [] bConnected = new boolean [nNodes];
        // start adding an arc at random
		int nNode1 = random.nextInt(nNodes);
		int nNode2 = random.nextInt(nNodes);
		if (nNode1 == nNode2) {nNode2 = (nNode1 + 1) % nNodes;}
		if (nNode2 < nNode1) {int h = nNode1; nNode1 = nNode2; nNode2 = h;}
		m_ParentSets[nNode2].addParent(nNode1, m_Instances);
		bConnected[nNode1] = true;
		bConnected[nNode2] = true;
		// Repeatedly, select one of the connected nodes, and one of 
		// the unconnected nodes and add an arc.
	    // All arcs point from lower to higher ordered nodes
	    // so that acyclicity is ensured.
		for (int iArc = 2; iArc < nNodes; iArc++ ) {
			int nNode = random.nextInt(nNodes);
			nNode1 = 0; //  one of the connected nodes
			while (nNode >= 0) {
				nNode1 = (nNode1 + 1) % nNodes;
				while (!bConnected[nNode1]) {
					nNode1 = (nNode1 + 1) % nNodes;
				}
				nNode--;
			}
			nNode = random.nextInt(nNodes);
			nNode2 = 0; //  one of the unconnected nodes
			while (nNode >= 0) {
				nNode2 = (nNode2 + 1) % nNodes;
				while (bConnected[nNode2]) {
					nNode2 = (nNode2 + 1) % nNodes;
				}
				nNode--;
			}
			if (nNode2 < nNode1) {int h = nNode1; nNode1 = nNode2; nNode2 = h;}
			m_ParentSets[nNode2].addParent(nNode1, m_Instances);
			bConnected[nNode1] = true;
			bConnected[nNode2] = true;
		}
	} // GenerateTree
	
	/** 
	 * GenerateRandomDistributions generates discrete conditional distribution tables
	 * for all nodes of a Bayes network once a network structure has been determined.
	 * @param nNodes number of nodes in the Bayes net 
	 * @param nValues number of values each of the nodes can take
	 */
    void generateRandomDistributions(int nNodes, int nValues) {
	    // Reserve space for CPTs
    	int nMaxParentCardinality = 1;
	    for (int iAttribute = 0; iAttribute < nNodes; iAttribute++) {
            if (m_ParentSets[iAttribute].getCardinalityOfParents() > nMaxParentCardinality) {
	             nMaxParentCardinality = m_ParentSets[iAttribute].getCardinalityOfParents();
            } 
        } 

        // Reserve plenty of memory
        m_Distributions = new Estimator[m_Instances.numAttributes()][nMaxParentCardinality];

        // estimate CPTs
        for (int iAttribute = 0; iAttribute < nNodes; iAttribute++) {
        	int [] nPs = new int [nValues + 1];
        	nPs[0] = 0;
        	nPs[nValues] = 1000;
            for (int iParent = 0; iParent < m_ParentSets[iAttribute].getCardinalityOfParents(); iParent++) {
            	// fill array with random nr's
            	for (int iValue = 1; iValue < nValues; iValue++)  {
            		nPs[iValue] = random.nextInt(1000);
            	}
            	// sort
            	for (int iValue = 1; iValue < nValues; iValue++)  {
	            	for (int iValue2 = iValue + 1; iValue2 < nValues; iValue2++)  {
	            		if (nPs[iValue2] < nPs[iValue]) {
	            			int h = nPs[iValue2]; nPs[iValue2] = nPs[iValue]; nPs[iValue] = h;
	            		}
	            	}
            	}
            	// assign to probability tables
            	DiscreteEstimatorBayes d = new DiscreteEstimatorBayes(nValues, getEstimator().getAlpha());
            	for (int iValue = 0; iValue < nValues; iValue++)  {
            		d.addValue(iValue, nPs[iValue + 1] - nPs[iValue]);
            	}
	            m_Distributions[iAttribute][iParent] = d;
            } 
        } 
    } // GenerateRandomDistributions
    
	/**
	 * GenerateInstances generates random instances sampling from the
	 * distribution represented by the Bayes network structure. It assumes
	 * a Bayes network structure has been initialized
	 * 
	 * @throws Exception if something goes wrong
	 */
	public void generateInstances () throws Exception {
	    int [] order = getOrder();
		for (int iInstance = 0; iInstance < m_nNrOfInstances; iInstance++) {
		    int nNrOfAtts = m_Instances.numAttributes();
			Instance instance = new DenseInstance(nNrOfAtts);
			instance.setDataset(m_Instances);
			for (int iAtt2 = 0; iAtt2 < nNrOfAtts; iAtt2++) {
			    int iAtt = order[iAtt2];

				double iCPT = 0;

				for (int iParent = 0; iParent < m_ParentSets[iAtt].getNrOfParents(); iParent++) {
				  int nParent = m_ParentSets[iAtt].getParent(iParent);
				  iCPT = iCPT * m_Instances.attribute(nParent).numValues() + instance.value(nParent);
				} 
	
				double fRandom = random.nextInt(1000) / 1000.0f;
				int iValue = 0;
				while (fRandom > m_Distributions[iAtt][(int) iCPT].getProbability(iValue)) {
					fRandom = fRandom - m_Distributions[iAtt][(int) iCPT].getProbability(iValue);
					iValue++ ;
				}
				instance.setValue(iAtt, iValue);
			}
			m_Instances.add(instance);
		}
	} // GenerateInstances

    /**
     * @throws Exception if there's a cycle in the graph
     */	
    int [] getOrder() throws Exception {
	int nNrOfAtts = m_Instances.numAttributes();
	int [] order = new int[nNrOfAtts];
	boolean [] bDone = new boolean[nNrOfAtts];
	for (int iAtt = 0; iAtt < nNrOfAtts; iAtt++) {
	    int iAtt2 = 0; 
	    boolean allParentsDone = false;
	    while (!allParentsDone && iAtt2 < nNrOfAtts) {
		if (!bDone[iAtt2]) {
		    allParentsDone = true;
		    int iParent = 0;
		    while (allParentsDone && iParent < m_ParentSets[iAtt2].getNrOfParents()) {
			allParentsDone = bDone[m_ParentSets[iAtt2].getParent(iParent++)];
		    }
		    if (allParentsDone && iParent == m_ParentSets[iAtt2].getNrOfParents()) {
			order[iAtt] = iAtt2;
			bDone[iAtt2] = true;
		    } else {
			iAtt2++;
		    }
		} else {
		    iAtt2++;
		}
	    }
	    if (!allParentsDone && iAtt2 == nNrOfAtts) {
		throw new Exception("There appears to be a cycle in the graph");
	    }
	}
	return order;
    } // getOrder

    	/**
    	 * Returns either the net (if BIF format) or the generated instances
    	 * 
    	 * @return either the net or the generated instances
    	 */
  	public String toString() {
  	  if (m_bGenerateNet) {
  	    return toXMLBIF03();
  	  }
  	  return m_Instances.toString();
  	} // toString
  	

	boolean m_bGenerateNet = false;
	int m_nNrOfNodes = 10;
	int m_nNrOfArcs = 10;
	int m_nNrOfInstances = 10;
	int m_nCardinality = 2;
	String m_sBIFFile = "";

	void setNrOfNodes(int nNrOfNodes) {m_nNrOfNodes = nNrOfNodes;}
	void setNrOfArcs(int nNrOfArcs) {m_nNrOfArcs = nNrOfArcs;}
	void setNrOfInstances(int nNrOfInstances) {m_nNrOfInstances = nNrOfInstances;}
	void setCardinality(int nCardinality) {m_nCardinality = nCardinality;}
	void setSeed(int nSeed) {m_nSeed = nSeed;}

	/**
	 * Returns an enumeration describing the available options
	 * 
	 * @return an enumeration of all the available options
	 */
	public Enumeration listOptions() {
		Vector newVector = new Vector(6);

		newVector.addElement(new Option("\tGenerate network (instead of instances)\n", "B", 0, "-B"));
		newVector.addElement(new Option("\tNr of nodes\n", "N", 1, "-N <integer>"));
		newVector.addElement(new Option("\tNr of arcs\n", "A", 1, "-A <integer>"));
		newVector.addElement(new Option("\tNr of instances\n", "M", 1, "-M <integer>"));
		newVector.addElement(new Option("\tCardinality of the variables\n", "C", 1, "-C <integer>"));
		newVector.addElement(new Option("\tSeed for random number generator\n", "S", 1, "-S <integer>"));
		newVector.addElement(new Option("\tThe BIF file to obtain the structure from.\n", "F", 1, "-F <file>"));

		return newVector.elements();
	} // listOptions

	/**
	 * Parses a given list of options. <p/>
	 * 
	 <!-- options-start -->
	 * Valid options are: <p/>
	 * 
	 * <pre> -B
	 *  Generate network (instead of instances)
	 * </pre>
	 * 
	 * <pre> -N &lt;integer&gt;
	 *  Nr of nodes
	 * </pre>
	 * 
	 * <pre> -A &lt;integer&gt;
	 *  Nr of arcs
	 * </pre>
	 * 
	 * <pre> -M &lt;integer&gt;
	 *  Nr of instances
	 * </pre>
	 * 
	 * <pre> -C &lt;integer&gt;
	 *  Cardinality of the variables
	 * </pre>
	 * 
	 * <pre> -S &lt;integer&gt;
	 *  Seed for random number generator
	 * </pre>
	 * 
	 * <pre> -F &lt;file&gt;
	 *  The BIF file to obtain the structure from.
	 * </pre>
	 * 
	 <!-- options-end -->
	 *
	 * @param options the list of options as an array of strings
	 * @exception Exception if an option is not supported
	 */
	public void setOptions(String[] options) throws Exception {
		m_bGenerateNet = Utils.getFlag('B', options);

		String sNrOfNodes = Utils.getOption('N', options);
		if (sNrOfNodes.length() != 0) {
		  setNrOfNodes(Integer.parseInt(sNrOfNodes));
		} else {
			setNrOfNodes(10);
		} 

		String sNrOfArcs = Utils.getOption('A', options);
		if (sNrOfArcs.length() != 0) {
		  setNrOfArcs(Integer.parseInt(sNrOfArcs));
		} else {
			setNrOfArcs(10);
		} 

		String sNrOfInstances = Utils.getOption('M', options);
		if (sNrOfInstances.length() != 0) {
		  setNrOfInstances(Integer.parseInt(sNrOfInstances));
		} else {
			setNrOfInstances(10);
		} 

		String sCardinality = Utils.getOption('C', options);
		if (sCardinality.length() != 0) {
		  setCardinality(Integer.parseInt(sCardinality));
		} else {
			setCardinality(2);
		} 

		String sSeed = Utils.getOption('S', options);
		if (sSeed.length() != 0) {
		  setSeed(Integer.parseInt(sSeed));
		} else {
			setSeed(1);
		} 

		String sBIFFile = Utils.getOption('F', options);
		if ((sBIFFile != null) && (sBIFFile != "")) {
			setBIFFile(sBIFFile);
		}
	} // setOptions

	/**
	 * Gets the current settings of the classifier.
	 * 
	 * @return an array of strings suitable for passing to setOptions
	 */
	public String[] getOptions() {
		String[] options = new String[13];
		int current = 0;
		if (m_bGenerateNet) {
		  options[current++] = "-B";
		} 

		options[current++] = "-N";
		options[current++] = "" + m_nNrOfNodes;

		options[current++] = "-A";
		options[current++] = "" + m_nNrOfArcs;

		options[current++] = "-M";
		options[current++] = "" + m_nNrOfInstances;

		options[current++] = "-C";
		options[current++] = "" + m_nCardinality;

		options[current++] = "-S";
		options[current++] = "" + m_nSeed;

                if (m_sBIFFile.length() != 0) {
                  options[current++] = "-F";
                  options[current++] = "" + m_sBIFFile;
                }

		// Fill up rest with empty strings, not nulls!
		while (current < options.length) {
			options[current++] = "";
		}

		return options;
	} // getOptions

    /**
     * prints all the options to stdout
     */
    protected static void printOptions(OptionHandler o) {
      Enumeration enm = o.listOptions();
      
      System.out.println("Options for " + o.getClass().getName() + ":\n");
      
      while (enm.hasMoreElements()) {
        Option option = (Option) enm.nextElement();
        System.out.println(option.synopsis());
        System.out.println(option.description());
      }
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
     * Main method
     * 
     * @param args the commandline parameters
     */
    static public void main(String [] args) {
		BayesNetGenerator b = new BayesNetGenerator();
    	try {
		if ( (args.length == 0) || (Utils.getFlag('h', args)) ) {
                        printOptions(b);
                        return;
		}
	    	b.setOptions(args);
	    	
	    	b.generateRandomNetwork();
	    	if (!b.m_bGenerateNet) { // skip if not required
				b.generateInstances();
	    	}
	    	System.out.println(b.toString());
    	} catch (Exception e) {
    		e.printStackTrace();
    		printOptions(b);
    	}
    } // main
    
} // class BayesNetGenerator
