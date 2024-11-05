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
 * GeneticSearch.java
 * Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 * 
 */
 
package weka.classifiers.bayes.net.search.global;

import weka.classifiers.bayes.BayesNet;
import weka.classifiers.bayes.net.ParentSet;
import weka.core.Instances;
import weka.core.Option;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.core.Utils;

import java.util.Enumeration;
import java.util.Random;
import java.util.Vector;

/** 
 <!-- globalinfo-start -->
 * This Bayes Network learning algorithm uses genetic search for finding a well scoring Bayes network structure. Genetic search works by having a population of Bayes network structures and allow them to mutate and apply cross over to get offspring. The best network structure found during the process is returned.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -L &lt;integer&gt;
 *  Population size</pre>
 * 
 * <pre> -A &lt;integer&gt;
 *  Descendant population size</pre>
 * 
 * <pre> -U &lt;integer&gt;
 *  Number of runs</pre>
 * 
 * <pre> -M
 *  Use mutation.
 *  (default true)</pre>
 * 
 * <pre> -C
 *  Use cross-over.
 *  (default true)</pre>
 * 
 * <pre> -O
 *  Use tournament selection (true) or maximum subpopulatin (false).
 *  (default false)</pre>
 * 
 * <pre> -R &lt;seed&gt;
 *  Random number seed</pre>
 * 
 * <pre> -mbc
 *  Applies a Markov Blanket correction to the network structure, 
 *  after a network structure is learned. This ensures that all 
 *  nodes in the network are part of the Markov blanket of the 
 *  classifier node.</pre>
 * 
 * <pre> -S [LOO-CV|k-Fold-CV|Cumulative-CV]
 *  Score type (LOO-CV,k-Fold-CV,Cumulative-CV)</pre>
 * 
 * <pre> -Q
 *  Use probabilistic or 0/1 scoring.
 *  (default probabilistic scoring)</pre>
 * 
 <!-- options-end -->
 * 
 * @author Remco Bouckaert (rrb@xm.co.nz)
 * @version $Revision: 1.5 $
 */
public class GeneticSearch 
    extends GlobalScoreSearchAlgorithm {

    /** for serialization */
    static final long serialVersionUID = 4236165533882462203L;
  
    /** number of runs **/
    int m_nRuns = 10;

	/** size of population **/
	int m_nPopulationSize = 10;

	/** size of descendant population **/
	int m_nDescendantPopulationSize = 100;

	/** use cross-over? **/
	boolean m_bUseCrossOver = true;

	/** use mutation? **/
	boolean m_bUseMutation = true;
	
	/** use tournament selection or take best sub-population **/
	boolean m_bUseTournamentSelection = false;	
	
	/** random number seed **/
	int m_nSeed = 1;
	
	/** random number generator **/
	Random m_random = null;


	/** used in BayesNetRepresentation for efficiently determining
	 * whether a number is square  
	 */
	static boolean [] g_bIsSquare;
	
	class BayesNetRepresentation implements RevisionHandler {
		/** number of nodes in network **/		
		int m_nNodes = 0;

		/** bit representation of parent sets 
		 * m_bits[iTail + iHead * m_nNodes] represents arc iTail->iHead
		 */
		boolean [] m_bits;
		
		/** score of represented network structure **/
		double m_fScore = 0.0f;
		
		/** 
		 * return score of represented network structure
		 * 
		 * @return the score
		 */
		public double getScore() {
			return m_fScore;
		} // getScore

		/** 
		 * c'tor 
		 * 
		 * @param nNodes the number of nodes
		 */
		BayesNetRepresentation (int nNodes) {
			m_nNodes = nNodes;
		} // c'tor
		
		/** initialize with a random structure by randomly placing
		 * m_nNodes arcs.
		 */
		public void randomInit() {
			do {
				m_bits = new boolean [m_nNodes * m_nNodes];
				for (int i = 0; i < m_nNodes; i++) {
					int iPos;
					do {
						iPos = m_random.nextInt(m_nNodes * m_nNodes);
					} while (isSquare(iPos));
					m_bits[iPos] = true;
				}
			} while (hasCycles());
			calcGlobalScore();
		}

		/** calculate score of current network representation
		 * As a side effect, the parent sets are set
		 */
		void calcGlobalScore() {
			// clear current network
			for (int iNode = 0; iNode < m_nNodes; iNode++) {
				ParentSet parentSet = m_BayesNet.getParentSet(iNode);
				while (parentSet.getNrOfParents() > 0) {
					parentSet.deleteLastParent(m_BayesNet.m_Instances);
				}
			}
			// insert arrows
			for (int iNode = 0; iNode < m_nNodes; iNode++) {
				ParentSet parentSet = m_BayesNet.getParentSet(iNode);
				for (int iNode2 = 0; iNode2 < m_nNodes; iNode2++) {
					if (m_bits[iNode2 + iNode * m_nNodes]) {
						parentSet.addParent(iNode2, m_BayesNet.m_Instances);
					}
				}
			}
			// calc score
			try {
				m_fScore = calcScore(m_BayesNet);
			} catch (Exception e) {
				// ignore
			}
		} // calcScore

		/** check whether there are cycles in the network
 		* 
 		* @return true if a cycle is found, false otherwise
		 */
		public boolean hasCycles() {
			// check for cycles
			boolean[] bDone = new boolean[m_nNodes];
			for (int iNode = 0; iNode < m_nNodes; iNode++) {

				// find a node for which all parents are 'done'
				boolean bFound = false;

				for (int iNode2 = 0; !bFound && iNode2 < m_nNodes; iNode2++) {
					if (!bDone[iNode2]) {
						boolean bHasNoParents = true;
						for (int iParent = 0; iParent < m_nNodes; iParent++) {
							if (m_bits[iParent + iNode2 * m_nNodes] && !bDone[iParent]) {
								bHasNoParents = false;
							}
						}
						if (bHasNoParents) {
							bDone[iNode2] = true;
							bFound = true;
						}
					}
				}
				if (!bFound) {
					return true;
				}
			}
			return false;
		} // hasCycles

		/** create clone of current object 
		 * @return cloned object
		 */
		BayesNetRepresentation copy() {
			BayesNetRepresentation b = new BayesNetRepresentation(m_nNodes);
			b.m_bits = new boolean [m_bits.length];
			for (int i = 0; i < m_nNodes * m_nNodes; i++) {
				b.m_bits[i] = m_bits[i];
			}
			b.m_fScore = m_fScore;
			return b;		
		} // copy

		/** Apply mutation operation to BayesNet
		 * Calculate score and as a side effect sets BayesNet parent sets.
		 */
		void mutate() {
			// flip a bit
			do {				
				int iBit;
				do {
					iBit = m_random.nextInt(m_nNodes * m_nNodes);
				} while (isSquare(iBit));
				
				m_bits[iBit] = !m_bits[iBit];
			} while (hasCycles());

			calcGlobalScore();
		} // mutate

		/** Apply cross-over operation to BayesNet 
		 * Calculate score and as a side effect sets BayesNet parent sets.
		 * @param other BayesNetRepresentation to cross over with
		 */
		void crossOver(BayesNetRepresentation other) {
			boolean [] bits = new boolean [m_bits.length];
			for (int i = 0; i < m_bits.length; i++) {
				bits[i] = m_bits[i];
			}
			int iCrossOverPoint = m_bits.length;
			do {
				// restore to original state
				for (int i = iCrossOverPoint; i < m_bits.length; i++) {
					m_bits[i] = bits[i];
				}
				// take all bits from cross-over points onwards
				iCrossOverPoint = m_random.nextInt(m_bits.length);
				for (int i = iCrossOverPoint; i < m_bits.length; i++) {
					m_bits[i] = other.m_bits[i];
				}
			} while (hasCycles());
			calcGlobalScore();
		} // crossOver
				
		/** check if number is square and initialize g_bIsSquare structure
		 * if necessary
		 * @param nNum number to check (should be below m_nNodes * m_nNodes)
		 * @return true if number is square
		 */
		boolean isSquare(int nNum) {
			if (g_bIsSquare == null || g_bIsSquare.length < nNum) {
				g_bIsSquare = new boolean [m_nNodes * m_nNodes];
				for (int i = 0; i < m_nNodes; i++) {
					g_bIsSquare[i * m_nNodes + i] = true;
				}
			}
			return g_bIsSquare[nNum];
		} // isSquare

		/**
		 * Returns the revision string.
		 * 
		 * @return		the revision
		 */
		public String getRevision() {
		  return RevisionUtils.extract("$Revision: 1.5 $");
		}
	} // class BayesNetRepresentation 
	    	
	/**
	 * search determines the network structure/graph of the network
	 * with a genetic search algorithm.
	 * 
	 * @param bayesNet the network to search
	 * @param instances the instances to use
	 * @throws Exception if population size doesn fit or neither cross-over or mutation was chosen
	 */
	protected void search(BayesNet bayesNet, Instances instances) throws Exception {
		// sanity check
		if (getDescendantPopulationSize() < getPopulationSize()) {
			throw new Exception ("Descendant PopulationSize should be at least Population Size");
		}
		if (!getUseCrossOver() && !getUseMutation()) {
			throw new Exception ("At least one of mutation or cross-over should be used");
		}
		
		m_random = new Random(m_nSeed);

		// keeps track of best structure found so far 
		BayesNet bestBayesNet;
		// keeps track of score pf best structure found so far 
		double fBestScore = calcScore(bayesNet);	

		// initialize bestBayesNet
		bestBayesNet = new BayesNet();
		bestBayesNet.m_Instances = instances;
		bestBayesNet.initStructure();
		copyParentSets(bestBayesNet, bayesNet);
		
                
        // initialize population        
		BayesNetRepresentation  [] population = new BayesNetRepresentation [getPopulationSize()];
        for (int i = 0; i < getPopulationSize(); i++) {
        	population[i] = new BayesNetRepresentation (instances.numAttributes());
			population[i].randomInit();
			if (population[i].getScore() > fBestScore) {
				copyParentSets(bestBayesNet, bayesNet);
				fBestScore = population[i].getScore();
				
			}
        }
        
        // go do the search        
        for (int iRun = 0; iRun < m_nRuns; iRun++) {
        	// create descendants
			BayesNetRepresentation  [] descendantPopulation = new BayesNetRepresentation  [getDescendantPopulationSize()];
			for (int i = 0; i < getDescendantPopulationSize(); i++) {
				descendantPopulation[i] = population[m_random.nextInt(getPopulationSize())].copy();
				if (getUseMutation()) {
					if (getUseCrossOver() && m_random.nextBoolean()) {
						descendantPopulation[i].crossOver(population[m_random.nextInt(getPopulationSize())]);						
					} else {
						descendantPopulation[i].mutate();								
					}
				} else {
					// use crossover
					descendantPopulation[i].crossOver(population[m_random.nextInt(getPopulationSize())]);
				}

				if (descendantPopulation[i].getScore() > fBestScore) {
					copyParentSets(bestBayesNet, bayesNet);
					fBestScore = descendantPopulation[i].getScore();
				}
			}
			// select new population
			boolean [] bSelected = new boolean [getDescendantPopulationSize()];
			for (int i = 0; i < getPopulationSize(); i++) {
				int iSelected = 0;
				if (m_bUseTournamentSelection) {
					// use tournament selection
					iSelected = m_random.nextInt(getDescendantPopulationSize());
					while (bSelected[iSelected]) {
						iSelected = (iSelected + 1) % getDescendantPopulationSize();
					}
					int iSelected2 =  m_random.nextInt(getDescendantPopulationSize());
					while (bSelected[iSelected2]) {
						iSelected2 = (iSelected2 + 1) % getDescendantPopulationSize();
					}
					if (descendantPopulation[iSelected2].getScore() > descendantPopulation[iSelected].getScore()) {
						iSelected = iSelected2;
					}
				} else {
					// find best scoring network in population
					while (bSelected[iSelected]) {
						iSelected++;
					}
					double fScore = descendantPopulation[iSelected].getScore();
					for (int j = 0; j < getDescendantPopulationSize(); j++) {
						if (!bSelected[j] && descendantPopulation[j].getScore() > fScore) {
							fScore = descendantPopulation[j].getScore();
							iSelected = j;
						}
					}
				}
				population[i] = descendantPopulation[iSelected];
				bSelected[iSelected] = true;
			}
        }
        
        // restore current network to best network
		copyParentSets(bayesNet, bestBayesNet);
		
		// free up memory
		bestBayesNet = null;
    } // search


	/** copyParentSets copies parent sets of source to dest BayesNet
	 * @param dest destination network
	 * @param source source network
	 */
	void copyParentSets(BayesNet dest, BayesNet source) {
		int nNodes = source.getNrOfNodes();
		// clear parent set first
		for (int iNode = 0; iNode < nNodes; iNode++) {
			dest.getParentSet(iNode).copy(source.getParentSet(iNode));
		}		
	} // CopyParentSets

    /**
    * @return number of runs
    */
    public int getRuns() {
        return m_nRuns;
    } // getRuns

    /**
     * Sets the number of runs
     * @param nRuns The number of runs to set
     */
    public void setRuns(int nRuns) {
        m_nRuns = nRuns;
    } // setRuns

	/**
	 * Returns an enumeration describing the available options.
	 *
	 * @return an enumeration of all the available options.
	 */
	public Enumeration listOptions() {
		Vector newVector = new Vector(7);

		newVector.addElement(new Option("\tPopulation size", "L", 1, "-L <integer>"));
		newVector.addElement(new Option("\tDescendant population size", "A", 1, "-A <integer>"));
		newVector.addElement(new Option("\tNumber of runs", "U", 1, "-U <integer>"));
		newVector.addElement(new Option("\tUse mutation.\n\t(default true)", "M", 0, "-M"));
		newVector.addElement(new Option("\tUse cross-over.\n\t(default true)", "C", 0, "-C"));
		newVector.addElement(new Option("\tUse tournament selection (true) or maximum subpopulatin (false).\n\t(default false)", "O", 0, "-O"));
		newVector.addElement(new Option("\tRandom number seed", "R", 1, "-R <seed>"));

		Enumeration enu = super.listOptions();
		while (enu.hasMoreElements()) {
			newVector.addElement(enu.nextElement());
		}
		return newVector.elements();
	} // listOptions

	/**
	 * Parses a given list of options. <p/>
	 *
	 <!-- options-start -->
	 * Valid options are: <p/>
	 * 
	 * <pre> -L &lt;integer&gt;
	 *  Population size</pre>
	 * 
	 * <pre> -A &lt;integer&gt;
	 *  Descendant population size</pre>
	 * 
	 * <pre> -U &lt;integer&gt;
	 *  Number of runs</pre>
	 * 
	 * <pre> -M
	 *  Use mutation.
	 *  (default true)</pre>
	 * 
	 * <pre> -C
	 *  Use cross-over.
	 *  (default true)</pre>
	 * 
	 * <pre> -O
	 *  Use tournament selection (true) or maximum subpopulatin (false).
	 *  (default false)</pre>
	 * 
	 * <pre> -R &lt;seed&gt;
	 *  Random number seed</pre>
	 * 
	 * <pre> -mbc
	 *  Applies a Markov Blanket correction to the network structure, 
	 *  after a network structure is learned. This ensures that all 
	 *  nodes in the network are part of the Markov blanket of the 
	 *  classifier node.</pre>
	 * 
	 * <pre> -S [LOO-CV|k-Fold-CV|Cumulative-CV]
	 *  Score type (LOO-CV,k-Fold-CV,Cumulative-CV)</pre>
	 * 
	 * <pre> -Q
	 *  Use probabilistic or 0/1 scoring.
	 *  (default probabilistic scoring)</pre>
	 * 
	 <!-- options-end -->
	 *
	 * @param options the list of options as an array of strings
	 * @throws Exception if an option is not supported
	 */
	public void setOptions(String[] options) throws Exception {
		String sPopulationSize = Utils.getOption('L', options);
		if (sPopulationSize.length() != 0) {
			setPopulationSize(Integer.parseInt(sPopulationSize));
		}
		String sDescendantPopulationSize = Utils.getOption('A', options);
		if (sDescendantPopulationSize.length() != 0) {
			setDescendantPopulationSize(Integer.parseInt(sDescendantPopulationSize));
		}
		String sRuns = Utils.getOption('U', options);
		if (sRuns.length() != 0) {
			setRuns(Integer.parseInt(sRuns));
		}
		String sSeed = Utils.getOption('R', options);
		if (sSeed.length() != 0) {
			setSeed(Integer.parseInt(sSeed));
		}
		setUseMutation(Utils.getFlag('M', options));
		setUseCrossOver(Utils.getFlag('C', options));
		setUseTournamentSelection(Utils.getFlag('O', options));
		
		super.setOptions(options);
	} // setOptions

	/**
	 * Gets the current settings of the search algorithm.
	 *
	 * @return an array of strings suitable for passing to setOptions
	 */
	public String[] getOptions() {
		String[] superOptions = super.getOptions();
		String[] options = new String[11 + superOptions.length];
		int current = 0;
		
		options[current++] = "-L";
		options[current++] = "" + getPopulationSize();

		options[current++] = "-A";
		options[current++] = "" + getDescendantPopulationSize();

		options[current++] = "-U";
		options[current++] = "" + getRuns();

		options[current++] = "-R";
		options[current++] = "" + getSeed();

		if (getUseMutation()) {
		  options[current++] = "-M";
		}
		if (getUseCrossOver()) {
		  options[current++] = "-C";
		}
		if (getUseTournamentSelection()) {
		  options[current++] = "-O";
		}

		// insert options from parent class
		for (int iOption = 0; iOption < superOptions.length; iOption++) {
			options[current++] = superOptions[iOption];
		}

		// Fill up rest with empty strings, not nulls!
		while (current < options.length) {
			options[current++] = "";
		}
		return options;
	} // getOptions

	/**
	 * @return whether cross-over is used
	 */
	public boolean getUseCrossOver() {
		return m_bUseCrossOver;
	}

	/**
	 * @return whether mutation is used
	 */
	public boolean getUseMutation() {
		return m_bUseMutation;
	}

	/**
	 * @return descendant population size
	 */
	public int getDescendantPopulationSize() {
		return m_nDescendantPopulationSize;
	}

	/**
	 * @return population size
	 */
	public int getPopulationSize() {
		return m_nPopulationSize;
	}

	/**
	 * @param bUseCrossOver sets whether cross-over is used
	 */
	public void setUseCrossOver(boolean bUseCrossOver) {
		m_bUseCrossOver = bUseCrossOver;
	}

	/**
	 * @param bUseMutation sets whether mutation is used
	 */
	public void setUseMutation(boolean bUseMutation) {
		m_bUseMutation = bUseMutation;
	}

	/**
	 * @return whether Tournament Selection (true) or Maximum Sub-Population (false) should be used
	 */
	public boolean getUseTournamentSelection() {
		return m_bUseTournamentSelection;
	}

	/**
	 * @param bUseTournamentSelection sets whether Tournament Selection or Maximum Sub-Population should be used
	 */
	public void setUseTournamentSelection(boolean bUseTournamentSelection) {
		m_bUseTournamentSelection = bUseTournamentSelection;
	}

	/**
	 * @param iDescendantPopulationSize sets descendant population size
	 */
	public void setDescendantPopulationSize(int iDescendantPopulationSize) {
		m_nDescendantPopulationSize = iDescendantPopulationSize;
	}

	/**
	 * @param iPopulationSize sets population size
	 */
	public void setPopulationSize(int iPopulationSize) {
		m_nPopulationSize = iPopulationSize;
	}

	/**
	* @return random number seed
	*/
	public int getSeed() {
		return m_nSeed;
	} // getSeed

	/**
	 * Sets the random number seed
	 * @param nSeed The number of the seed to set
	 */
	public void setSeed(int nSeed) {
		m_nSeed = nSeed;
	} // setSeed

	/**
	 * This will return a string describing the classifier.
	 * @return The string.
	 */
	public String globalInfo() {
		return "This Bayes Network learning algorithm uses genetic search for finding a well scoring " +
		"Bayes network structure. Genetic search works by having a population of Bayes network structures " +
		"and allow them to mutate and apply cross over to get offspring. The best network structure " +
		"found during the process is returned.";
	} // globalInfo
	
	/**
	 * @return a string to describe the Runs option.
	 */
	public String runsTipText() {
	  return "Sets the number of generations of Bayes network structure populations.";
	} // runsTipText
	
	/**
	 * @return a string to describe the Seed option.
	 */
	public String seedTipText() {
	  return "Initialization value for random number generator." +
	  " Setting the seed allows replicability of experiments.";
	} // seedTipText

	/**
	 * @return a string to describe the Population Size option.
	 */
	public String populationSizeTipText() {
	  return "Sets the size of the population of network structures that is selected each generation.";
	} // populationSizeTipText

	/**
	 * @return a string to describe the Descendant Population Size option.
	 */
	public String descendantPopulationSizeTipText() {
	  return "Sets the size of the population of descendants that is created each generation.";
	} // descendantPopulationSizeTipText

	/**
	 * @return a string to describe the Use Mutation option.
	 */
	public String useMutationTipText() {
		return "Determines whether mutation is allowed. Mutation flips a bit in the bit " +
			"representation of the network structure. At least one of mutation or cross-over " +
			"should be used.";
	} // useMutationTipText

	/**
	 * @return a string to describe the Use Cross-Over option.
	 */
	public String useCrossOverTipText() {
		return "Determines whether cross-over is allowed. Cross over combined the bit " +
			"representations of network structure by taking a random first k bits of one" +
			"and adding the remainder of the other. At least one of mutation or cross-over " +
			"should be used.";
	} // useCrossOverTipText

	/**
	 * @return a string to describe the Use Tournament Selection option.
	 */
	public String useTournamentSelectionTipText() {
		return "Determines the method of selecting a population. When set to true, tournament " +
			"selection is used (pick two at random and the highest is allowed to continue). " +
			"When set to false, the top scoring network structures are selected.";
	} // useTournamentSelectionTipText

	/**
	 * Returns the revision string.
	 * 
	 * @return		the revision
	 */
	public String getRevision() {
	  return RevisionUtils.extract("$Revision: 1.5 $");
	}
} // GeneticSearch
