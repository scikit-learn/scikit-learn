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
 * ICSSearchAlgorithm.java
 * Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 * 
 */


package weka.classifiers.bayes.net.search.ci;

import weka.classifiers.bayes.BayesNet;
import weka.classifiers.bayes.net.ParentSet;
import weka.core.Instances;
import weka.core.Option;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.core.Utils;

import java.io.FileReader;
import java.util.Enumeration;
import java.util.Vector;

/** 
 <!-- globalinfo-start -->
 * This Bayes Network learning algorithm uses conditional independence tests to find a skeleton, finds V-nodes and applies a set of rules to find the directions of the remaining arrows.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -cardinality &lt;num&gt;
 *  When determining whether an edge exists a search is performed 
 *  for a set Z that separates the nodes. MaxCardinality determines 
 *  the maximum size of the set Z. This greatly influences the 
 *  length of the search. (default 2)</pre>
 * 
 * <pre> -mbc
 *  Applies a Markov Blanket correction to the network structure, 
 *  after a network structure is learned. This ensures that all 
 *  nodes in the network are part of the Markov blanket of the 
 *  classifier node.</pre>
 * 
 * <pre> -S [BAYES|MDL|ENTROPY|AIC|CROSS_CLASSIC|CROSS_BAYES]
 *  Score type (BAYES, BDeu, MDL, ENTROPY and AIC)</pre>
 * 
 <!-- options-end -->
 * 
 * @author Remco Bouckaert
 * @version $Revision: 1.8 $
 */ 
public class ICSSearchAlgorithm 
    extends CISearchAlgorithm {

    /** for serialization */
    static final long serialVersionUID = -2510985917284798576L;
  
    /**
     * returns the name of the attribute with the given index
     * 
     * @param iAttribute the index of the attribute
     * @return the name of the attribute
     */
    String name(int iAttribute) {
      return m_instances.attribute(iAttribute).name();
    }
    
    /**
     * returns the number of attributes
     * 
     * @return the number of attributes
     */
    int maxn() {
      return m_instances.numAttributes();
    }
    
    /** maximum size of separating set **/
    private int m_nMaxCardinality = 2; 

    /**
     * sets the cardinality
     * 
     * @param nMaxCardinality the max cardinality
     */
    public void setMaxCardinality(int nMaxCardinality) {
      m_nMaxCardinality = nMaxCardinality;
    }
    
    /**
     * returns the max cardinality
     * 
     * @return the max cardinality
     */
    public int getMaxCardinality() {
      return m_nMaxCardinality;
    }
	

	
    class SeparationSet
        implements RevisionHandler {
      
        public int [] m_set;
        
        /**
         * constructor
         */
        public SeparationSet() {
			m_set= new int [getMaxCardinality() + 1];
        } // c'tor

        public boolean contains(int nItem) {
        	for (int iItem = 0; iItem < getMaxCardinality() && m_set[iItem] != -1; iItem++) {
        		if (m_set[iItem] == nItem) {
					return true;
				}
        	}
			return false;
    	} // contains
        
        /**
         * Returns the revision string.
         * 
         * @return		the revision
         */
        public String getRevision() {
          return RevisionUtils.extract("$Revision: 1.8 $");
        }

    } // class sepset

	/**
	 * Search for Bayes network structure using ICS algorithm
	 * @param bayesNet datastructure to build network structure for
	 * @param instances data set to learn from
	 * @throws Exception if something goes wrong
	 */
	protected void search(BayesNet bayesNet, Instances instances) throws Exception {
        // init
        m_BayesNet = bayesNet;
        m_instances = instances;

        boolean edges[][] = new boolean [maxn() + 1][];
		boolean [] [] arrows = new boolean [maxn() + 1][];
        SeparationSet [] [] sepsets = new SeparationSet [maxn() + 1][];
        for (int iNode = 0 ; iNode < maxn() + 1; iNode++) {
            edges[iNode] = new boolean[maxn()];
			arrows[iNode] = new boolean[maxn()]; 
            sepsets[iNode] = new SeparationSet[maxn()];
        }

        calcDependencyGraph(edges, sepsets);
        calcVeeNodes(edges, arrows, sepsets);
        calcArcDirections(edges, arrows);

		// transfrom into BayesNet datastructure 
		for (int iNode = 0; iNode < maxn(); iNode++) {
			// clear parent set of AttributeX
			ParentSet oParentSet = m_BayesNet.getParentSet(iNode);
			while (oParentSet.getNrOfParents() > 0) {
				oParentSet.deleteLastParent(m_instances);
			}
			for (int iParent = 0; iParent < maxn(); iParent++) {
				if (arrows[iParent][iNode]) {
					oParentSet.addParent(iParent, m_instances);
				}
			}
		}
	} // search
	
	
	/** CalcDependencyGraph determines the skeleton of the BayesNetwork by
	 * starting with a complete graph and removing edges (a--b) if it can
	 * find a set Z such that a and b conditionally independent given Z.
	 * The set Z is found by trying all possible subsets of nodes adjacent 
	 * to a and b, first of size 0, then of size 1, etc. up to size 
	 * m_nMaxCardinality
	 * @param edges boolean matrix representing the edges
	 * @param sepsets set of separating sets
	 */
	void calcDependencyGraph(boolean[][] edges, SeparationSet[][] sepsets) {
		/*calc undirected graph a-b iff D(a,S,b) for all S)*/
		SeparationSet oSepSet;

		for (int iNode1 = 0; iNode1 < maxn(); iNode1++) { 
			/*start with complete graph*/
			for (int iNode2 = 0; iNode2 < maxn(); iNode2++) {
				edges[iNode1][iNode2] = true;
			}
		}
		for (int iNode1 = 0; iNode1 < maxn(); iNode1++) {
			edges[iNode1][iNode1] = false;
		}

		for (int iCardinality = 0; iCardinality <= getMaxCardinality(); iCardinality++) {
			for (int iNode1 = 0; iNode1 <= maxn() - 2; iNode1++) {
				for (int iNode2 = iNode1 + 1; iNode2 < maxn(); iNode2++) {
					if (edges[iNode1][iNode2]) {
						oSepSet = existsSepSet(iNode1, iNode2, iCardinality, edges);
						if (oSepSet != null) {
							edges[iNode1][iNode2] = false;
							edges[iNode2][iNode1] = false;
							sepsets[iNode1][iNode2] = oSepSet;
							sepsets[iNode2][iNode1] = oSepSet;
							// report separating set
							System.err.print("I(" + name(iNode1) + ", {");
							for (int iNode3 = 0; iNode3 < iCardinality; iNode3++) {
								System.err.print(name(oSepSet.m_set[iNode3]) + " ");
							}
							System.err.print("} ," + name(iNode2) + ")\n");
						}
					}
				}
			}
			// report current state of dependency graph
			System.err.print(iCardinality + " ");
			for (int iNode1 = 0; iNode1 < maxn(); iNode1++) {
				System.err.print(name(iNode1) + " ");
			}
			System.err.print('\n');
			for (int iNode1 = 0; iNode1 < maxn(); iNode1++) {
				for (int iNode2 = 0; iNode2 < maxn(); iNode2++) {
					if (edges[iNode1][iNode2])
						System.err.print("X ");
					else
						System.err.print(". ");
				}
				System.err.print(name(iNode1) + " ");
				System.err.print('\n');
			}
		}
	} /*CalcDependencyGraph*/

	/** ExistsSepSet tests if a separating set Z of node a and b exists of given 
	 * cardiniality exists. 
	 * The set Z is found by trying all possible subsets of nodes adjacent 
	 * to both a and b of the requested cardinality.
	 * @param iNode1 index of first node a
	 * @param iNode2 index of second node b
	 * @param nCardinality size of the separating set Z
	 * @param edges
	 * @return SeparationSet containing set that separates iNode1 and iNode2 or null if no such set exists
	 */
    SeparationSet existsSepSet(int iNode1, int iNode2, int nCardinality, boolean [] [] edges)
    {
        /*Test if a separating set of node d and e exists of cardiniality k*/
//        int iNode1_, iNode2_;
        int iNode3, iZ;
		SeparationSet Z = new SeparationSet();
		Z.m_set[nCardinality] = -1;

//        iNode1_ = iNode1;
//        iNode2_ = iNode2;

		// find first candidate separating set Z
        if (nCardinality > 0) {
            Z.m_set[0] = next(-1, iNode1, iNode2, edges);
            iNode3 = 1;
            while (iNode3 < nCardinality) {
              Z.m_set[iNode3] = next(Z.m_set[iNode3 - 1], iNode1, iNode2, edges);
              iNode3++;
            }
        }

        if (nCardinality > 0) {
	        iZ = maxn() - Z.m_set[nCardinality - 1] - 1;
        } else {
    	    iZ = 0;
        }
        

        while (iZ >= 0)
        {  
        	//check if candidate separating set makes iNode2_ and iNode1_ independent
            if (isConditionalIndependent(iNode2, iNode1, Z.m_set, nCardinality))	{
                return Z;
            }
			// calc next candidate separating set
            if (nCardinality > 0) {
                Z.m_set[nCardinality - 1] = next(Z.m_set[nCardinality - 1], iNode1, iNode2, edges);
            }
            iZ = nCardinality - 1;   
            while (iZ >= 0 && Z.m_set[iZ] >= maxn()) {
                iZ = nCardinality - 1;
                while (iZ >= 0 && Z.m_set[iZ] >= maxn()) {
                	iZ--;
                }
                if (iZ < 0) {
                    break;
                }
                Z.m_set[iZ] = next(Z.m_set[iZ], iNode1, iNode2, edges);
                for (iNode3 = iZ + 1; iNode3 < nCardinality; iNode3++) {
                    Z.m_set[iNode3] = next(Z.m_set[iNode3 - 1], iNode1, iNode2, edges);
                }
                iZ = nCardinality - 1;
            }
        }

        return null;
    }  /*ExistsSepSet*/

	/** 
	 * determine index of node that makes next candidate separating set
	 * adjacent to iNode1 and iNode2, but not iNode2 itself
	 * @param x index of current node
	 * @param iNode1 first node
	 * @param iNode2 second node (must be larger than iNode1)
	 * @param edges skeleton so far
	 * @return int index of next node adjacent to iNode1 after x
	 */
	int next(int x, int iNode1, int iNode2, boolean [] [] edges)
	{
		x++;
		while (x < maxn() && (!edges[iNode1][x] || !edges[iNode2][x] ||x == iNode2)) {
			x++;
		}
		return x;
	}  /*next*/


	/** CalcVeeNodes tries to find V-nodes, i.e. nodes a,b,c such that
	 * a->c<-b and a-/-b. These nodes are identified by finding nodes
	 * a,b,c in the skeleton such that a--c, c--b and a-/-b and furthermore
	 * c is not in the set Z that separates a and b
	 * @param edges skeleton
	 * @param arrows resulting partially directed skeleton after all V-nodes 
	 * have been identified
	 * @param sepsets separating sets
	 */
	void calcVeeNodes(
		boolean[][] edges,
		boolean[][] arrows,
		SeparationSet[][] sepsets) {

		// start with complete empty graph
		for (int iNode1 = 0; iNode1 < maxn(); iNode1++) {
			for (int iNode2 = 0; iNode2 < maxn(); iNode2++) {
				arrows[iNode1][iNode2] = false;
			}
		}

		for (int iNode1 = 0; iNode1 < maxn() - 1; iNode1++) {
			for (int iNode2 = iNode1 + 1; iNode2 < maxn(); iNode2++) {
				if (!edges[iNode1][iNode2]) { /*i nonadj j*/
					for (int iNode3 = 0; iNode3 < maxn(); iNode3++) {
						if ((iNode3 != iNode1
							&& iNode3 != iNode2
							&& edges[iNode1][iNode3]
							&& edges[iNode2][iNode3])
							& (!sepsets[iNode1][iNode2].contains(iNode3))) {
							arrows[iNode1][iNode3] = true; /*add arc i->k*/
							arrows[iNode2][iNode3] = true; /*add arc j->k*/
						}
					}
				}
			}
		}
	} // CalcVeeNodes


	/** CalcArcDirections assigns directions to edges that remain after V-nodes have
	 * been identified. The arcs are directed using the following rules:
	   Rule 1: i->j--k & i-/-k => j->k
	   Rule 2: i->j->k & i--k => i->k
	   Rule 3  m
			 /|\
			i | k  => m->j
	i->j<-k  \|/
			  j
	
	   Rule 4  m
			 / \
			i---k  => i->m & k->m
	  i->j   \ /
			  j
	   Rule 5: if no edges are directed then take a random one (first we can find)
	 * @param edges skeleton
	 * @param arrows resulting fully directed DAG
	 */
	void calcArcDirections(boolean[][] edges, boolean[][] arrows) {
		/*give direction to remaining arcs*/
		int i, j, k, m;
		boolean bFound;

		do {
			bFound = false;

			/*Rule 1: i->j--k & i-/-k => j->k*/

			for (i = 0; i < maxn(); i++) {
				for (j = 0; j < maxn(); j++) {
					if (i != j && arrows[i][j]) {
						for (k = 0; k < maxn(); k++) {
							if (i != k
								&& j != k
								&& edges[j][k]
								&& !edges[i][k]
								&& !arrows[j][k]
								&& !arrows[k][j]) {
								arrows[j][k] = true;
								bFound = true;
							}
						}
					}
				}
			}

			/*Rule 2: i->j->k & i--k => i->k*/

			for (i = 0; i < maxn(); i++) {
				for (j = 0; j < maxn(); j++) {
					if (i != j && arrows[i][j]) {
						for (k = 0; k < maxn(); k++) {
							if (i != k
								&& j != k
								&& edges[i][k]
								&& arrows[j][k]
								&& !arrows[i][k]
								&& !arrows[k][i]) {
								arrows[i][k] = true;
								bFound = true;
							}
						}
					}
				}
			}

			/* Rule 3  m
			         /|\
			        i | k  => m->j
			i->j<-k  \|/
			          j
			*/
			for (i = 0; i < maxn(); i++) {
				for (j = 0; j < maxn(); j++) {
					if (i != j && arrows[i][j]) {
						for (k = 0; k < maxn(); k++) {
							if (k != i
								&& k != j
								&& arrows[k][j]
								&& !edges[k][i]) {
								for (m = 0; m < maxn(); m++) {
									if (m != i
										&& m != j
										&& m != k
										&& edges[m][i]
										&& !arrows[m][i]
										&& !arrows[i][m]
										&& edges[m][j]
										&& !arrows[m][j]
										&& !arrows[j][m]
										&& edges[m][k]
										&& !arrows[m][k]
										&& !arrows[k][m]) {
										arrows[m][j] = true;
										bFound = true;
									}
								}
							}
						}
					}
				}
			}

			/* Rule 4  m
			         / \
			        i---k  => i->m & k->m
			  i->j   \ /
			          j
			*/
			for (i = 0; i < maxn(); i++) {
				for (j = 0; j < maxn(); j++) {
					if (i != j && arrows[j][i]) {
						for (k = 0; k < maxn(); k++) {
							if (k != i
								&& k != j
								&& edges[k][j]
								&& !arrows[k][j]
								&& !arrows[j][k]
								&& edges[k][i]
								&& !arrows[k][i]
								&& !arrows[i][k]) {
								for (m = 0; m < maxn(); m++) {
									if (m != i
										&& m != j
										&& m != k
										&& edges[m][i]
										&& !arrows[m][i]
										&& !arrows[i][m]
										&& edges[m][k]
										&& !arrows[m][k]
										&& !arrows[k][m]) {
										arrows[i][m] = true;
										arrows[k][m] = true;
										bFound = true;
									}
								}
							}
						}
					}
				}
			}

			/*Rule 5: if no edges are directed then take a random one (first we can find)*/

			if (!bFound) {
				i = 0;
				while (!bFound && i < maxn()) {
					j = 0;
					while (!bFound && j < maxn()) {
						if (edges[i][j]
							&& !arrows[i][j]
							&& !arrows[j][i]) {
							arrows[i][j] = true;
							bFound = true;
						}
						j++;
					}
					i++;
				}
			}

		}
		while (bFound);

	} // CalcArcDirections

	/**
	 * Returns an enumeration describing the available options.
	 *
	 * @return an enumeration of all the available options.
	 */
	public Enumeration listOptions() {
	  Vector result = new Vector();
	  
	  result.addElement(new Option(
                "\tWhen determining whether an edge exists a search is performed \n"
              + "\tfor a set Z that separates the nodes. MaxCardinality determines \n"
              + "\tthe maximum size of the set Z. This greatly influences the \n"
              + "\tlength of the search. (default 2)",
	      "cardinality", 1, "-cardinality <num>"));
	  
	  Enumeration en = super.listOptions();
	  while (en.hasMoreElements())
	    result.addElement(en.nextElement());
	  
	  return result.elements();
	} // listOption
	
	/**
	 * Parses a given list of options. <p/>
	 *
	 <!-- options-start -->
	 * Valid options are: <p/>
	 * 
	 * <pre> -cardinality &lt;num&gt;
	 *  When determining whether an edge exists a search is performed 
	 *  for a set Z that separates the nodes. MaxCardinality determines 
	 *  the maximum size of the set Z. This greatly influences the 
	 *  length of the search. (default 2)</pre>
	 * 
	 * <pre> -mbc
	 *  Applies a Markov Blanket correction to the network structure, 
	 *  after a network structure is learned. This ensures that all 
	 *  nodes in the network are part of the Markov blanket of the 
	 *  classifier node.</pre>
	 * 
	 * <pre> -S [BAYES|MDL|ENTROPY|AIC|CROSS_CLASSIC|CROSS_BAYES]
	 *  Score type (BAYES, BDeu, MDL, ENTROPY and AIC)</pre>
	 * 
	 <!-- options-end -->
	 * 
	 * @param options the list of options as an array of strings
	 * @throws Exception if an option is not supported
	 */
	public void setOptions(String[] options) throws Exception {
	  String        tmpStr;
	  
	  tmpStr = Utils.getOption("cardinality", options);
	  if (tmpStr.length() != 0)
	    setMaxCardinality(Integer.parseInt(tmpStr));
	  else
            setMaxCardinality(2);
            
          super.setOptions(options);
	} // setOptions
	
	/**
	 * Gets the current settings of the Classifier.
	 *
	 * @return an array of strings suitable for passing to setOptions
	 */
	public String[] getOptions() {
	  Vector        result;
	  String[]      options;
	  int           i;
	  
	  result  = new Vector();
	  options = super.getOptions();
	  for (i = 0; i < options.length; i++)
	    result.add(options[i]);
	  
	  result.add("-cardinality");
	  result.add("" + getMaxCardinality());
	  
	  return (String[]) result.toArray(new String[result.size()]);
	} // getOptions
	

	/**
	 * @return a string to describe the MaxCardinality option.
	 */
	public String maxCardinalityTipText() {
	  return "When determining whether an edge exists a search is performed for a set Z "+
	  "that separates the nodes. MaxCardinality determines the maximum size of the set Z. " +
	  "This greatly influences the length of the search. Default value is 2.";
	} // maxCardinalityTipText

	/**
	 * This will return a string describing the search algorithm.
	 * @return The string.
	 */
	public String globalInfo() {
	  return "This Bayes Network learning algorithm uses conditional independence tests " +
	  "to find a skeleton, finds V-nodes and applies a set of rules to find the directions " +
	  "of the remaining arrows.";
	}

	/**
	 * Returns the revision string.
	 * 
	 * @return		the revision
	 */
	public String getRevision() {
	  return RevisionUtils.extract("$Revision: 1.8 $");
	}

	/**
	 * for testing the class
	 * 
	 * @param argv the commandline parameters
	 */
	static public void main(String [] argv) {
		try {
			BayesNet b = new BayesNet();
			b.setSearchAlgorithm( new ICSSearchAlgorithm());
			Instances instances = new Instances(new FileReader("C:\\eclipse\\workspace\\weka\\data\\contact-lenses.arff"));
			instances.setClassIndex(instances.numAttributes() - 1);
			b.buildClassifier(instances);
			System.out.println(b.toString());
		} catch (Exception e) {
			e.printStackTrace();
		}
	} // main

} // class ICSSearchAlgorithm
