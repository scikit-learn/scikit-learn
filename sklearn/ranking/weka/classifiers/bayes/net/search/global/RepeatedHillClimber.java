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
 * RepeatedHillClimber.java
 * Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 * 
 */
 
package weka.classifiers.bayes.net.search.global;

import weka.classifiers.bayes.BayesNet;
import weka.classifiers.bayes.net.ParentSet;
import weka.core.Instances;
import weka.core.Option;
import weka.core.RevisionUtils;
import weka.core.Utils;

import java.util.Enumeration;
import java.util.Random;
import java.util.Vector;

/** 
 <!-- globalinfo-start -->
 * This Bayes Network learning algorithm repeatedly uses hill climbing starting with a randomly generated network structure and return the best structure of the various runs.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -U &lt;integer&gt;
 *  Number of runs</pre>
 * 
 * <pre> -A &lt;seed&gt;
 *  Random number seed</pre>
 * 
 * <pre> -P &lt;nr of parents&gt;
 *  Maximum number of parents</pre>
 * 
 * <pre> -R
 *  Use arc reversal operation.
 *  (default false)</pre>
 * 
 * <pre> -N
 *  Initial structure is empty (instead of Naive Bayes)</pre>
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
 * @version $Revision: 1.6 $
 */
public class RepeatedHillClimber 
    extends HillClimber {

    /** for serialization */
    static final long serialVersionUID = -7359197180460703069L;
  
    /** number of runs **/
    int m_nRuns = 10;
    /** random number seed **/
    int m_nSeed = 1;
    /** random number generator **/
    Random m_random;

	/**
	* search determines the network structure/graph of the network
	* with the repeated hill climbing.
	* 
	* @param bayesNet the network to use
	* @param instances the data to use
	* @throws Exception if something goes wrong
	**/
	protected void search(BayesNet bayesNet, Instances instances) throws Exception {
		m_random = new Random(getSeed());
		// keeps track of score pf best structure found so far 
		double fBestScore;	
		double fCurrentScore = calcScore(bayesNet);

		// keeps track of best structure found so far 
		BayesNet bestBayesNet;

		// initialize bestBayesNet
		fBestScore = fCurrentScore;
		bestBayesNet = new BayesNet();
		bestBayesNet.m_Instances = instances;
		bestBayesNet.initStructure();
		copyParentSets(bestBayesNet, bayesNet);
		
                
        // go do the search        
        for (int iRun = 0; iRun < m_nRuns; iRun++) {
        	// generate random nework
        	generateRandomNet(bayesNet, instances);

        	// search
        	super.search(bayesNet, instances);

			// calculate score
			fCurrentScore = calcScore(bayesNet);

			// keep track of best network seen so far
			if (fCurrentScore > fBestScore) {
				fBestScore = fCurrentScore;
				copyParentSets(bestBayesNet, bayesNet);
			}
        }
        
        // restore current network to best network
		copyParentSets(bayesNet, bestBayesNet);
		
		// free up memory
		bestBayesNet = null;
    } // search

	/**
	 * 
	 * @param bayesNet
	 * @param instances
	 */
	void generateRandomNet(BayesNet bayesNet, Instances instances) {
		int nNodes = instances.numAttributes();
		// clear network
		for (int iNode = 0; iNode < nNodes; iNode++) {
			ParentSet parentSet = bayesNet.getParentSet(iNode);
			while (parentSet.getNrOfParents() > 0) {
				parentSet.deleteLastParent(instances);
			}
		}
		
		// initialize as naive Bayes?
		if (getInitAsNaiveBayes()) {
			int iClass = instances.classIndex();
			// initialize parent sets to have arrow from classifier node to
			// each of the other nodes
			for (int iNode = 0; iNode < nNodes; iNode++) {
				if (iNode != iClass) {
					bayesNet.getParentSet(iNode).addParent(iClass, instances);
				}
			}
		}

		// insert random arcs
		int nNrOfAttempts = m_random.nextInt(nNodes * nNodes);
		for (int iAttempt = 0; iAttempt < nNrOfAttempts; iAttempt++) {
			int iTail = m_random.nextInt(nNodes);
			int iHead = m_random.nextInt(nNodes);
			if (bayesNet.getParentSet(iHead).getNrOfParents() < getMaxNrOfParents() &&
			    addArcMakesSense(bayesNet, instances, iHead, iTail)) {
					bayesNet.getParentSet(iHead).addParent(iTail, instances);
			}
		}
	} // generateRandomNet

	/** 
	 * copyParentSets copies parent sets of source to dest BayesNet
	 * 
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
     * Returns the number of runs
     * 
     * @return number of runs
     */
    public int getRuns() {
        return m_nRuns;
    } // getRuns

    /**
     * Sets the number of runs
     * 
     * @param nRuns The number of runs to set
     */
    public void setRuns(int nRuns) {
        m_nRuns = nRuns;
    } // setRuns

	/**
	 * Returns the random seed
	 * 
	 * @return random number seed
	 */
	public int getSeed() {
		return m_nSeed;
	} // getSeed

	/**
	 * Sets the random number seed
	 * 
	 * @param nSeed The number of the seed to set
	 */
	public void setSeed(int nSeed) {
		m_nSeed = nSeed;
	} // setSeed

	/**
	 * Returns an enumeration describing the available options.
	 *
	 * @return an enumeration of all the available options.
	 */
	public Enumeration listOptions() {
		Vector newVector = new Vector(4);

		newVector.addElement(new Option("\tNumber of runs", "U", 1, "-U <integer>"));
		newVector.addElement(new Option("\tRandom number seed", "A", 1, "-A <seed>"));

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
	 * <pre> -U &lt;integer&gt;
	 *  Number of runs</pre>
	 * 
	 * <pre> -A &lt;seed&gt;
	 *  Random number seed</pre>
	 * 
	 * <pre> -P &lt;nr of parents&gt;
	 *  Maximum number of parents</pre>
	 * 
	 * <pre> -R
	 *  Use arc reversal operation.
	 *  (default false)</pre>
	 * 
	 * <pre> -N
	 *  Initial structure is empty (instead of Naive Bayes)</pre>
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
		String sRuns = Utils.getOption('U', options);
		if (sRuns.length() != 0) {
			setRuns(Integer.parseInt(sRuns));
		}
		
		String sSeed = Utils.getOption('A', options);
		if (sSeed.length() != 0) {
			setSeed(Integer.parseInt(sSeed));
		}

		super.setOptions(options);
	} // setOptions

	/**
	 * Gets the current settings of the search algorithm.
	 *
	 * @return an array of strings suitable for passing to setOptions
	 */
	public String[] getOptions() {
		String[] superOptions = super.getOptions();
		String[] options = new String[7 + superOptions.length];
		int current = 0;

		options[current++] = "-U";
		options[current++] = "" + getRuns();

		options[current++] = "-A";
		options[current++] = "" + getSeed();

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
	 * This will return a string describing the classifier.
	 * 
	 * @return The string.
	 */
	public String globalInfo() {
		return "This Bayes Network learning algorithm repeatedly uses hill climbing starting " +
		"with a randomly generated network structure and return the best structure of the " +
		"various runs.";
	} // globalInfo
	
	/**
	 * @return a string to describe the Runs option.
	 */
	public String runsTipText() {
	  return "Sets the number of times hill climbing is performed.";
	} // runsTipText

	/**
	 * @return a string to describe the Seed option.
	 */
	public String seedTipText() {
	  return "Initialization value for random number generator." +
	  " Setting the seed allows replicability of experiments.";
	} // seedTipText

	/**
	 * Returns the revision string.
	 * 
	 * @return		the revision
	 */
	public String getRevision() {
	  return RevisionUtils.extract("$Revision: 1.6 $");
	}
}
