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
 * LocalScoreSearchAlgorithm.java
 * Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 * 
 */
 
package weka.classifiers.bayes.net.search.local;

import weka.classifiers.bayes.BayesNet;
import weka.classifiers.bayes.net.ParentSet;
import weka.classifiers.bayes.net.search.SearchAlgorithm;
import weka.core.Instances;
import weka.core.Instance;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.core.Statistics;
import weka.core.Tag;
import weka.core.Option;
import weka.core.SelectedTag;

import java.util.Vector;
import java.util.Enumeration;

/** 
 <!-- globalinfo-start -->
 * The ScoreBasedSearchAlgorithm class supports Bayes net structure search algorithms that are based on maximizing scores (as opposed to for example conditional independence based search algorithms).
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
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
 * @version $Revision: 5195 $
 */
public class LocalScoreSearchAlgorithm 
	extends SearchAlgorithm {

  	/** for serialization */
  	static final long serialVersionUID = 3325995552474190374L;
  	
	/** points to Bayes network for which a structure is searched for **/
	BayesNet m_BayesNet;
	
	/**
	 * default constructor
	 */
	public LocalScoreSearchAlgorithm() {
	} // c'tor
	
	/**
	 * constructor
	 * 
	 * @param bayesNet the network
	 * @param instances the data
	 */
	public LocalScoreSearchAlgorithm(BayesNet bayesNet, Instances instances) {
		m_BayesNet = bayesNet;
//		m_Instances = instances;
	} // c'tor
	
	/**
	 * Holds prior on count
	 */
	double m_fAlpha = 0.5;

	/** the score types */
	public static final Tag[] TAGS_SCORE_TYPE = {
	  new Tag(Scoreable.BAYES, "BAYES"),
	  new Tag(Scoreable.BDeu, "BDeu"),
	  new Tag(Scoreable.MDL, "MDL"),
	  new Tag(Scoreable.ENTROPY, "ENTROPY"),
	  new Tag(Scoreable.AIC, "AIC")
	};

	/**
	 * Holds the score type used to measure quality of network
	 */
	int m_nScoreType = Scoreable.BAYES;

	/**
	 * logScore returns the log of the quality of a network
	 * (e.g. the posterior probability of the network, or the MDL
	 * value).
	 * @param nType score type (Bayes, MDL, etc) to calculate score with
	 * @return log score.
	 */
    public double logScore(int nType) {
    	if (m_BayesNet.m_Distributions == null) {return 0;}
        if (nType < 0) {
            nType = m_nScoreType;
        }

        double fLogScore = 0.0;
        
        Instances instances = m_BayesNet.m_Instances;

        for (int iAttribute = 0; iAttribute < instances.numAttributes(); iAttribute++) {
        	int nCardinality = m_BayesNet.getParentSet(iAttribute).getCardinalityOfParents();
            for (int iParent = 0; iParent < nCardinality; iParent++) {
                fLogScore += ((Scoreable) m_BayesNet.m_Distributions[iAttribute][iParent]).logScore(nType, nCardinality);
            }

            switch (nType) {
                case (Scoreable.MDL) :
                    {
                        fLogScore -= 0.5
                            * m_BayesNet.getParentSet(iAttribute).getCardinalityOfParents()
                            * (instances.attribute(iAttribute).numValues() - 1)
                            * Math.log(instances.numInstances());
                    }
                    break;
                case (Scoreable.AIC) :
                    {
                        fLogScore -= m_BayesNet.getParentSet(iAttribute).getCardinalityOfParents()
                            * (instances.attribute(iAttribute).numValues() - 1);
                    }
                    break;
            }
        }

        return fLogScore;
    } // logScore

	/**
	* buildStructure determines the network structure/graph of the network
	* with the K2 algorithm, restricted by its initial structure (which can
	* be an empty graph, or a Naive Bayes graph.
	* 
	* @param bayesNet the network
	* @param instances the data to use
	* @throws Exception if something goes wrong
	*/
	public void buildStructure (BayesNet bayesNet, Instances instances) throws Exception {
		m_BayesNet = bayesNet;
		super.buildStructure(bayesNet, instances);
	} // buildStructure


	/**
	 * Calc Node Score for given parent set
	 * 
	 * @param nNode node for which the score is calculate
	 * @return log score
	 */
	public double calcNodeScore(int nNode) {
		if (m_BayesNet.getUseADTree() && m_BayesNet.getADTree() != null) {
			return calcNodeScoreADTree(nNode);
		} else {
			return calcNodeScorePlain(nNode);
		}
	}

	/**
	 * helper function for CalcNodeScore above using the ADTree data structure
	 * 
	 * @param nNode node for which the score is calculate
	 * @return log score
	 */
	private double calcNodeScoreADTree(int nNode) {
		Instances instances = m_BayesNet.m_Instances;
		ParentSet oParentSet = m_BayesNet.getParentSet(nNode);
		// get set of parents, insert iNode
		int nNrOfParents = oParentSet.getNrOfParents();
		int[] nNodes = new int[nNrOfParents + 1];
		for (int iParent = 0; iParent < nNrOfParents; iParent++) {
			nNodes[iParent] = oParentSet.getParent(iParent);
		}
		nNodes[nNrOfParents] = nNode;

		// calculate offsets
		int[] nOffsets = new int[nNrOfParents + 1];
		int nOffset = 1;
		nOffsets[nNrOfParents] = 1;
		nOffset *= instances.attribute(nNode).numValues();
		for (int iNode = nNrOfParents - 1; iNode >= 0; iNode--) {
			nOffsets[iNode] = nOffset;
			nOffset *= instances.attribute(nNodes[iNode]).numValues();
		}

		// sort nNodes & offsets
		for (int iNode = 1; iNode < nNodes.length; iNode++) {
			int iNode2 = iNode;
			while (iNode2 > 0 && nNodes[iNode2] < nNodes[iNode2 - 1]) {
				int h = nNodes[iNode2];
				nNodes[iNode2] = nNodes[iNode2 - 1];
				nNodes[iNode2 - 1] = h;
				h = nOffsets[iNode2];
				nOffsets[iNode2] = nOffsets[iNode2 - 1];
				nOffsets[iNode2 - 1] = h;
				iNode2--;
			}
		}

		// get counts from ADTree
		int nCardinality = oParentSet.getCardinalityOfParents();
		int numValues = instances.attribute(nNode).numValues();
		int[] nCounts = new int[nCardinality * numValues];
		//if (nNrOfParents > 1) {

		m_BayesNet.getADTree().getCounts(nCounts, nNodes, nOffsets, 0, 0, false);

		return calcScoreOfCounts(nCounts, nCardinality, numValues, instances);
	} // CalcNodeScore

	private double calcNodeScorePlain(int nNode) {
		Instances instances = m_BayesNet.m_Instances;
		ParentSet oParentSet = m_BayesNet.getParentSet(nNode);

		// determine cardinality of parent set & reserve space for frequency counts
		int nCardinality = oParentSet.getCardinalityOfParents();
		int numValues = instances.attribute(nNode).numValues();
		int[] nCounts = new int[nCardinality * numValues];

		// initialize (don't need this?)
		for (int iParent = 0; iParent < nCardinality * numValues; iParent++) {
			nCounts[iParent] = 0;
		}

		// estimate distributions
		Enumeration enumInsts = instances.enumerateInstances();

		while (enumInsts.hasMoreElements()) {
			Instance instance = (Instance) enumInsts.nextElement();

			// updateClassifier;
			double iCPT = 0;

			for (int iParent = 0; iParent < oParentSet.getNrOfParents(); iParent++) {
				int nParent = oParentSet.getParent(iParent);

				iCPT = iCPT * instances.attribute(nParent).numValues() + instance.value(nParent);
			}

			nCounts[numValues * ((int) iCPT) + (int) instance.value(nNode)]++;
		}

		return calcScoreOfCounts(nCounts, nCardinality, numValues, instances);
	} // CalcNodeScore

	/**
	 * utility function used by CalcScore and CalcNodeScore to determine the score
	 * based on observed frequencies.
	 * 
	 * @param nCounts array with observed frequencies
	 * @param nCardinality ardinality of parent set
	 * @param numValues number of values a node can take
	 * @param instances to calc score with
	 * @return log score
	 */
	protected double calcScoreOfCounts(int[] nCounts, int nCardinality, int numValues, Instances instances) {

		// calculate scores using the distributions
		double fLogScore = 0.0;

		for (int iParent = 0; iParent < nCardinality; iParent++) {
			switch (m_nScoreType) {

				case (Scoreable.BAYES) :
					{
						double nSumOfCounts = 0;

						for (int iSymbol = 0; iSymbol < numValues; iSymbol++) {
							if (m_fAlpha + nCounts[iParent * numValues + iSymbol] != 0) {
								fLogScore += Statistics.lnGamma(m_fAlpha + nCounts[iParent * numValues + iSymbol]);
								nSumOfCounts += m_fAlpha + nCounts[iParent * numValues + iSymbol];
							}
						}

						if (nSumOfCounts != 0) {
							fLogScore -= Statistics.lnGamma(nSumOfCounts);
						}

						if (m_fAlpha != 0) {
							fLogScore -= numValues * Statistics.lnGamma(m_fAlpha);
							fLogScore += Statistics.lnGamma(numValues * m_fAlpha);
						}
					}

					break;
                case (Scoreable.BDeu) :
                {
                    double nSumOfCounts = 0;

                    for (int iSymbol = 0; iSymbol < numValues; iSymbol++) {
                        if (m_fAlpha + nCounts[iParent * numValues + iSymbol] != 0) {
                            fLogScore += Statistics.lnGamma(1.0/(numValues * nCardinality) + nCounts[iParent * numValues + iSymbol]);
                            nSumOfCounts += 1.0/(numValues * nCardinality) + nCounts[iParent * numValues + iSymbol];
                        }
                    }
                    fLogScore -= Statistics.lnGamma(nSumOfCounts);

                    fLogScore -= numValues * Statistics.lnGamma(1.0/(numValues * nCardinality));
                    fLogScore += Statistics.lnGamma(1.0/nCardinality);
                }
	                break;

				case (Scoreable.MDL) :

				case (Scoreable.AIC) :

				case (Scoreable.ENTROPY) :
					{
						double nSumOfCounts = 0;

						for (int iSymbol = 0; iSymbol < numValues; iSymbol++) {
							nSumOfCounts += nCounts[iParent * numValues + iSymbol];
						}

						for (int iSymbol = 0; iSymbol < numValues; iSymbol++) {
							if (nCounts[iParent * numValues + iSymbol] > 0) {
								fLogScore += nCounts[iParent * numValues
									+ iSymbol] * Math.log(nCounts[iParent * numValues + iSymbol] / nSumOfCounts);
							}
						}
					}

					break;

				default :
					{
					}
			}
		}

		switch (m_nScoreType) {

			case (Scoreable.MDL) :
				{
					fLogScore -= 0.5 * nCardinality * (numValues - 1) * Math.log(instances.numInstances());

					// it seems safe to assume that numInstances>0 here
				}

				break;

			case (Scoreable.AIC) :
				{
					fLogScore -= nCardinality * (numValues - 1);
				}

				break;
		}

		return fLogScore;
	} // CalcNodeScore

	protected double calcScoreOfCounts2(int[][] nCounts, int nCardinality, int numValues, Instances instances) {

		// calculate scores using the distributions
		double fLogScore = 0.0;

		for (int iParent = 0; iParent < nCardinality; iParent++) {
			switch (m_nScoreType) {

				case (Scoreable.BAYES) :
					{
						double nSumOfCounts = 0;

						for (int iSymbol = 0; iSymbol < numValues; iSymbol++) {
							if (m_fAlpha + nCounts[iParent][iSymbol] != 0) {
								fLogScore += Statistics.lnGamma(m_fAlpha + nCounts[iParent][iSymbol]);
								nSumOfCounts += m_fAlpha + nCounts[iParent][iSymbol];
							}
						}

						if (nSumOfCounts != 0) {
							fLogScore -= Statistics.lnGamma(nSumOfCounts);
						}

						if (m_fAlpha != 0) {
							fLogScore -= numValues * Statistics.lnGamma(m_fAlpha);
							fLogScore += Statistics.lnGamma(numValues * m_fAlpha);
						}
					}

					break;

				case (Scoreable.BDeu) :
				{
					double nSumOfCounts = 0;

					for (int iSymbol = 0; iSymbol < numValues; iSymbol++) {
						if (m_fAlpha + nCounts[iParent][iSymbol] != 0) {
							fLogScore += Statistics.lnGamma(1.0/(numValues * nCardinality) + nCounts[iParent][iSymbol]);
							nSumOfCounts += 1.0/(numValues * nCardinality) + nCounts[iParent][iSymbol];
						}
					}
					fLogScore -= Statistics.lnGamma(nSumOfCounts);

					fLogScore -= numValues * Statistics.lnGamma(1.0/(nCardinality * numValues));
					fLogScore += Statistics.lnGamma(1.0/ nCardinality);
				}
					break;

				case (Scoreable.MDL) :

				case (Scoreable.AIC) :

				case (Scoreable.ENTROPY) :
					{
						double nSumOfCounts = 0;

						for (int iSymbol = 0; iSymbol < numValues; iSymbol++) {
							nSumOfCounts += nCounts[iParent][iSymbol];
						}

						for (int iSymbol = 0; iSymbol < numValues; iSymbol++) {
							if (nCounts[iParent][iSymbol] > 0) {
								fLogScore += nCounts[iParent][iSymbol]
									* Math.log(nCounts[iParent][iSymbol] / nSumOfCounts);
							}
						}
					}

					break;

				default :
					{
					}
			}
		}

		switch (m_nScoreType) {

			case (Scoreable.MDL) :
				{
					fLogScore -= 0.5 * nCardinality * (numValues - 1) * Math.log(instances.numInstances());

					// it seems safe to assume that numInstances>0 here
				}

				break;

			case (Scoreable.AIC) :
				{
					fLogScore -= nCardinality * (numValues - 1);
				}

				break;
		}

		return fLogScore;
	} // CalcNodeScore


	/**
	 * Calc Node Score With AddedParent
	 * 
	 * @param nNode node for which the score is calculate
	 * @param nCandidateParent candidate parent to add to the existing parent set
	 * @return log score
	 */
	public double calcScoreWithExtraParent(int nNode, int nCandidateParent) {
		ParentSet oParentSet = m_BayesNet.getParentSet(nNode);

		// sanity check: nCandidateParent should not be in parent set already
		if (oParentSet.contains(nCandidateParent)) {
				return -1e100;
		}

		// set up candidate parent
		oParentSet.addParent(nCandidateParent, m_BayesNet.m_Instances);

		// calculate the score
		double logScore = calcNodeScore(nNode);

		// delete temporarily added parent
		oParentSet.deleteLastParent(m_BayesNet.m_Instances);

		return logScore;
	} // CalcScoreWithExtraParent


	/**
	 * Calc Node Score With Parent Deleted
	 * 
	 * @param nNode node for which the score is calculate
	 * @param nCandidateParent candidate parent to delete from the existing parent set
	 * @return log score
	 */
	public double calcScoreWithMissingParent(int nNode, int nCandidateParent) {
		ParentSet oParentSet = m_BayesNet.getParentSet(nNode);

		// sanity check: nCandidateParent should be in parent set already
		if (!oParentSet.contains( nCandidateParent)) {
				return -1e100;
		}

		// set up candidate parent
		int iParent = oParentSet.deleteParent(nCandidateParent, m_BayesNet.m_Instances);

		// calculate the score
		double logScore = calcNodeScore(nNode);

		// restore temporarily deleted parent
		oParentSet.addParent(nCandidateParent, iParent, m_BayesNet.m_Instances);

		return logScore;
	} // CalcScoreWithMissingParent

	/**
	 * set quality measure to be used in searching for networks.
	 * 
	 * @param newScoreType the new score type
	 */
	public void setScoreType(SelectedTag newScoreType) {
		if (newScoreType.getTags() == TAGS_SCORE_TYPE) {
			m_nScoreType = newScoreType.getSelectedTag().getID();
		}
	}

	/**
	 * get quality measure to be used in searching for networks.
	 * @return quality measure
	 */
	public SelectedTag getScoreType() {
		return new SelectedTag(m_nScoreType, TAGS_SCORE_TYPE);
	}

	/**
	 * 
	 * @param bMarkovBlanketClassifier
	 */
	public void setMarkovBlanketClassifier(boolean bMarkovBlanketClassifier) {
	  super.setMarkovBlanketClassifier(bMarkovBlanketClassifier);
	}
	
	/**
	 * 
	 * @return
	 */
	public boolean getMarkovBlanketClassifier() {
	  return super.getMarkovBlanketClassifier();
	}

	/**
	 * Returns an enumeration describing the available options
	 * 
	 * @return an enumeration of all the available options
	 */
	public Enumeration listOptions() {
		Vector newVector = new Vector();

		newVector.addElement(new Option(
		    "\tApplies a Markov Blanket correction to the network structure, \n"
		  + "\tafter a network structure is learned. This ensures that all \n"
		  + "\tnodes in the network are part of the Markov blanket of the \n"
		  + "\tclassifier node.",
		  "mbc", 0, "-mbc"));
      
		newVector.addElement(
			new Option(
				"\tScore type (BAYES, BDeu, MDL, ENTROPY and AIC)",
				"S",
				1,
				"-S [BAYES|MDL|ENTROPY|AIC|CROSS_CLASSIC|CROSS_BAYES]"));

		return newVector.elements();
	} // listOptions

	/**
	 * Parses a given list of options. <p/>
	 *
	 <!-- options-start -->
	 * Valid options are: <p/>
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

	  	setMarkovBlanketClassifier(Utils.getFlag("mbc", options));

		String sScore = Utils.getOption('S', options);

		if (sScore.compareTo("BAYES") == 0) {
			setScoreType(new SelectedTag(Scoreable.BAYES, TAGS_SCORE_TYPE));
		}
		if (sScore.compareTo("BDeu") == 0) {
			setScoreType(new SelectedTag(Scoreable.BDeu, TAGS_SCORE_TYPE));
		}
		if (sScore.compareTo("MDL") == 0) {
			setScoreType(new SelectedTag(Scoreable.MDL, TAGS_SCORE_TYPE));
		}
		if (sScore.compareTo("ENTROPY") == 0) {
			setScoreType(new SelectedTag(Scoreable.ENTROPY, TAGS_SCORE_TYPE));
		}
		if (sScore.compareTo("AIC") == 0) {
			setScoreType(new SelectedTag(Scoreable.AIC, TAGS_SCORE_TYPE));
		}
	} // setOptions

	/**
	 * Gets the current settings of the search algorithm.
	 *
	 * @return an array of strings suitable for passing to setOptions
	 */
	public String[] getOptions() {
                String[] superOptions = super.getOptions();
		String[] options = new String[3 + superOptions.length];
		int current = 0;

		if (getMarkovBlanketClassifier())
		  options[current++] = "-mbc";

		options[current++] = "-S";

		switch (m_nScoreType) {

			case (Scoreable.BAYES) :
				options[current++] = "BAYES";
				break;

			case (Scoreable.BDeu) :
				options[current++] = "BDeu";
				break;

			case (Scoreable.MDL) :
				options[current++] = "MDL";
				break;

			case (Scoreable.ENTROPY) :
				options[current++] = "ENTROPY";

				break;

			case (Scoreable.AIC) :
				options[current++] = "AIC";
				break;
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
	 * @return a string to describe the ScoreType option.
	 */
	public String scoreTypeTipText() {
		return "The score type determines the measure used to judge the quality of a"
			+ " network structure. It can be one of Bayes, BDeu, Minimum Description Length (MDL),"
			+ " Akaike Information Criterion (AIC), and Entropy.";
	}
	
	/**
	 * @return a string to describe the MarkovBlanketClassifier option.
	 */
	public String markovBlanketClassifierTipText() {
	  return super.markovBlanketClassifierTipText();
	}

	/**
	 * This will return a string describing the search algorithm.
	 * @return The string.
	 */
	public String globalInfo() {
	  return 
	      "The ScoreBasedSearchAlgorithm class supports Bayes net "
	    + "structure search algorithms that are based on maximizing "
	    + "scores (as opposed to for example conditional independence "
	    + "based search algorithms).";
	} // globalInfo

	/**
	 * Returns the revision string.
	 * 
	 * @return		the revision
	 */
	public String getRevision() {
	  return RevisionUtils.extract("$Revision: 5195 $");
	}
}
