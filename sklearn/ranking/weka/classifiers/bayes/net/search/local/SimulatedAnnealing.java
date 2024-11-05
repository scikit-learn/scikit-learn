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
 * SimulatedAnnealing.java
 * Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 * 
 */
 
package weka.classifiers.bayes.net.search.local;

import weka.classifiers.bayes.BayesNet;
import weka.core.Instances;
import weka.core.Option;
import weka.core.RevisionUtils;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformation.Type;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformationHandler;
import weka.core.Utils;

import java.util.Enumeration;
import java.util.Random;
import java.util.Vector;

/** 
 <!-- globalinfo-start -->
 * This Bayes Network learning algorithm uses the general purpose search method of simulated annealing to find a well scoring network structure.<br/>
 * <br/>
 * For more information see:<br/>
 * <br/>
 * R.R. Bouckaert (1995). Bayesian Belief Networks: from Construction to Inference. Utrecht, Netherlands.
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;phdthesis{Bouckaert1995,
 *    address = {Utrecht, Netherlands},
 *    author = {R.R. Bouckaert},
 *    institution = {University of Utrecht},
 *    title = {Bayesian Belief Networks: from Construction to Inference},
 *    year = {1995}
 * }
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -A &lt;float&gt;
 *  Start temperature</pre>
 * 
 * <pre> -U &lt;integer&gt;
 *  Number of runs</pre>
 * 
 * <pre> -D &lt;float&gt;
 *  Delta temperature</pre>
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
 * <pre> -S [BAYES|MDL|ENTROPY|AIC|CROSS_CLASSIC|CROSS_BAYES]
 *  Score type (BAYES, BDeu, MDL, ENTROPY and AIC)</pre>
 * 
 <!-- options-end -->
 * 
 * @author Remco Bouckaert (rrb@xm.co.nz)
 * @version $Revision: 1.6 $
 */
public class SimulatedAnnealing 
	extends LocalScoreSearchAlgorithm
	implements TechnicalInformationHandler {
  
  	/** for serialization */
  	static final long serialVersionUID = 6951955606060513191L;

    	/** start temperature **/
	double m_fTStart = 10;

	/** change in temperature at every run **/
	double m_fDelta = 0.999;

	/** number of runs **/
	int m_nRuns = 10000;

	/** use the arc reversal operator **/
	boolean m_bUseArcReversal = false;

	/** random number seed **/
	int m_nSeed = 1;

	/** random number generator **/
	Random m_random;

	/**
	 * Returns an instance of a TechnicalInformation object, containing 
	 * detailed information about the technical background of this class,
	 * e.g., paper reference or book this class is based on.
	 * 
	 * @return the technical information about this class
	 */
	public TechnicalInformation getTechnicalInformation() {
	  TechnicalInformation 	result;
	  
	  result = new TechnicalInformation(Type.PHDTHESIS);
	  result.setValue(Field.AUTHOR, "R.R. Bouckaert");
	  result.setValue(Field.YEAR, "1995");
	  result.setValue(Field.TITLE, "Bayesian Belief Networks: from Construction to Inference");
	  result.setValue(Field.INSTITUTION, "University of Utrecht");
	  result.setValue(Field.ADDRESS, "Utrecht, Netherlands");
	  
	  return result;
	}

    /**
     * 
     * @param bayesNet the network
     * @param instances the data to use
     * @throws Exception if something goes wrong
     */
    public void search (BayesNet bayesNet, Instances instances) throws Exception {
		m_random = new Random(m_nSeed);
		
        // determine base scores
        double [] fBaseScores = new double [instances.numAttributes()];
		double fCurrentScore = 0;
        for (int iAttribute = 0; iAttribute < instances.numAttributes(); iAttribute++) {
            fBaseScores[iAttribute] = calcNodeScore(iAttribute);
			fCurrentScore += fBaseScores[iAttribute];
        }

		// keep track of best scoring network
		double fBestScore = fCurrentScore;
		BayesNet bestBayesNet = new BayesNet();
		bestBayesNet.m_Instances = instances;
		bestBayesNet.initStructure();
		copyParentSets(bestBayesNet, bayesNet);

        double fTemp = m_fTStart;
        for (int iRun = 0; iRun < m_nRuns; iRun++) {
            boolean bRunSucces = false;
            double fDeltaScore = 0.0;
            while (!bRunSucces) {
	            // pick two nodes at random
	            int iTailNode = Math.abs(m_random.nextInt()) % instances.numAttributes();
	            int iHeadNode = Math.abs(m_random.nextInt()) % instances.numAttributes();
	            while (iTailNode == iHeadNode) {
		            iHeadNode = Math.abs(m_random.nextInt()) % instances.numAttributes();
	            }
	            if (isArc(bayesNet, iHeadNode, iTailNode)) {
                    bRunSucces = true;
	                // either try a delete
                    bayesNet.getParentSet(iHeadNode).deleteParent(iTailNode, instances);
                    double fScore = calcNodeScore(iHeadNode);
                    fDeltaScore = fScore - fBaseScores[iHeadNode];
//System.out.println("Try delete " + iTailNode + "->" + iHeadNode + " dScore = " + fDeltaScore);                    
                    if (fTemp * Math.log((Math.abs(m_random.nextInt()) % 10000)/10000.0  + 1e-100) < fDeltaScore) {
//System.out.println("success!!!");                    
						fCurrentScore += fDeltaScore;
                        fBaseScores[iHeadNode] = fScore;
                    } else {
                        // roll back
                        bayesNet.getParentSet(iHeadNode).addParent(iTailNode, instances);
                    }
	            } else {
	                // try to add an arc
	                if (addArcMakesSense(bayesNet, instances, iHeadNode, iTailNode)) {
                        bRunSucces = true;
                        double fScore = calcScoreWithExtraParent(iHeadNode, iTailNode);
                        fDeltaScore = fScore - fBaseScores[iHeadNode];
//System.out.println("Try add " + iTailNode + "->" + iHeadNode + " dScore = " + fDeltaScore);                    
                        if (fTemp * Math.log((Math.abs(m_random.nextInt()) % 10000)/10000.0  + 1e-100) < fDeltaScore) {
//System.out.println("success!!!");                    
                            bayesNet.getParentSet(iHeadNode).addParent(iTailNode, instances);
                            fBaseScores[iHeadNode] = fScore;
							fCurrentScore += fDeltaScore;
                        }
	                }
	            }
            }
			if (fCurrentScore > fBestScore) {
				copyParentSets(bestBayesNet, bayesNet);				
			}
            fTemp = fTemp * m_fDelta;
        }

		copyParentSets(bayesNet, bestBayesNet);
    } // buildStructure 
	
	/** CopyParentSets copies parent sets of source to dest BayesNet
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
     * @return double
     */
    public double getDelta() {
        return m_fDelta;
    }

    /**
     * @return double
     */
    public double getTStart() {
        return m_fTStart;
    }

    /**
     * @return int
     */
    public int getRuns() {
        return m_nRuns;
    }

    /**
     * Sets the m_fDelta.
     * @param fDelta The m_fDelta to set
     */
    public void setDelta(double fDelta) {
        m_fDelta = fDelta;
    }

    /**
     * Sets the m_fTStart.
     * @param fTStart The m_fTStart to set
     */
    public void setTStart(double fTStart) {
        m_fTStart = fTStart;
    }

    /**
     * Sets the m_nRuns.
     * @param nRuns The m_nRuns to set
     */
    public void setRuns(int nRuns) {
        m_nRuns = nRuns;
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
	 * Returns an enumeration describing the available options.
	 *
	 * @return an enumeration of all the available options.
	 */
	public Enumeration listOptions() {
		Vector newVector = new Vector(3);

		newVector.addElement(new Option("\tStart temperature", "A", 1, "-A <float>"));
		newVector.addElement(new Option("\tNumber of runs", "U", 1, "-U <integer>"));
		newVector.addElement(new Option("\tDelta temperature", "D", 1, "-D <float>"));
		newVector.addElement(new Option("\tRandom number seed", "R", 1, "-R <seed>"));

		Enumeration enu = super.listOptions();
		while (enu.hasMoreElements()) {
			newVector.addElement(enu.nextElement());
		}
		return newVector.elements();
	}

	/**
	 * Parses a given list of options. <p/>
	 *
	 <!-- options-start -->
	 * Valid options are: <p/>
	 * 
	 * <pre> -A &lt;float&gt;
	 *  Start temperature</pre>
	 * 
	 * <pre> -U &lt;integer&gt;
	 *  Number of runs</pre>
	 * 
	 * <pre> -D &lt;float&gt;
	 *  Delta temperature</pre>
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
	 * <pre> -S [BAYES|MDL|ENTROPY|AIC|CROSS_CLASSIC|CROSS_BAYES]
	 *  Score type (BAYES, BDeu, MDL, ENTROPY and AIC)</pre>
	 * 
	 <!-- options-end -->
	 *
	 * @param options the list of options as an array of strings
	 * @throws Exception if an option is not supported
	 */
	public void setOptions(String[] options) throws Exception {
		String sTStart = Utils.getOption('A', options);
		if (sTStart.length() != 0) {
			setTStart(Double.parseDouble(sTStart));
		}
		String sRuns = Utils.getOption('U', options);
		if (sRuns.length() != 0) {
			setRuns(Integer.parseInt(sRuns));
		}
		String sDelta = Utils.getOption('D', options);
		if (sDelta.length() != 0) {
			setDelta(Double.parseDouble(sDelta));
		}
		String sSeed = Utils.getOption('R', options);
		if (sSeed.length() != 0) {
			setSeed(Integer.parseInt(sSeed));
		}
		super.setOptions(options);
	}

	/**
	 * Gets the current settings of the search algorithm.
	 *
	 * @return an array of strings suitable for passing to setOptions
	 */
	public String[] getOptions() {
		String[] superOptions = super.getOptions();
		String[] options = new String[8 + superOptions.length];
		int current = 0;
		options[current++] = "-A";
		options[current++] = "" + getTStart();

		options[current++] = "-U";
		options[current++] = "" + getRuns();

		options[current++] = "-D";
		options[current++] = "" + getDelta();

		options[current++] = "-R";
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
	}

	/**
	 * This will return a string describing the classifier.
	 * @return The string.
	 */
	public String globalInfo() {
		return 
		    "This Bayes Network learning algorithm uses the general purpose search method "
		  + "of simulated annealing to find a well scoring network structure.\n\n"
		  + "For more information see:\n\n"
		  + getTechnicalInformation().toString();
	} // globalInfo
	
	/**
	 * @return a string to describe the TStart option.
	 */
	public String TStartTipText() {
	  return "Sets the start temperature of the simulated annealing search. "+
	  "The start temperature determines the probability that a step in the 'wrong' direction in the " +
	  "search space is accepted. The higher the temperature, the higher the probability of acceptance.";
	} // TStartTipText

	/**
	 * @return a string to describe the Runs option.
	 */
	public String runsTipText() {
	  return "Sets the number of iterations to be performed by the simulated annealing search.";
	} // runsTipText
	
	/**
	 * @return a string to describe the Delta option.
	 */
	public String deltaTipText() {
	  return "Sets the factor with which the temperature (and thus the acceptance probability of " +
	  	"steps in the wrong direction in the search space) is decreased in each iteration.";
	} // deltaTipText

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

} // SimulatedAnnealing
