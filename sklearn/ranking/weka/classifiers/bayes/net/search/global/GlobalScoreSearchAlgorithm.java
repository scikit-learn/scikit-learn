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
 * GlobalScoreSearchAlgorithm.java
 * Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 * 
 */

package weka.classifiers.bayes.net.search.global;

import weka.classifiers.bayes.BayesNet;
import weka.classifiers.bayes.net.ParentSet;
import weka.classifiers.bayes.net.search.SearchAlgorithm;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.RevisionUtils;
import weka.core.SelectedTag;
import weka.core.Tag;
import weka.core.Utils;

import java.util.Enumeration;
import java.util.Vector;

/** 
 <!-- globalinfo-start -->
 * This Bayes Network learning algorithm uses cross validation to estimate classification accuracy.
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
 * <pre> -S [LOO-CV|k-Fold-CV|Cumulative-CV]
 *  Score type (LOO-CV,k-Fold-CV,Cumulative-CV)</pre>
 * 
 * <pre> -Q
 *  Use probabilistic or 0/1 scoring.
 *  (default probabilistic scoring)</pre>
 * 
 <!-- options-end -->
 * 
 * @author Remco Bouckaert
 * @version $Revision: 1.10 $
 */
public class GlobalScoreSearchAlgorithm 
	extends SearchAlgorithm {

  	/** for serialization */
  	static final long serialVersionUID = 7341389867906199781L;
	
	/** points to Bayes network for which a structure is searched for **/
	BayesNet m_BayesNet;
	
	/** toggle between scoring using accuracy = 0-1 loss (when false) or class probabilities (when true) **/
	boolean m_bUseProb = true;
	
	/** number of folds for k-fold cross validation **/
	int m_nNrOfFolds = 10;

	/** constant for score type: LOO-CV */
	final static int LOOCV = 0;
	/** constant for score type: k-fold-CV */
	final static int KFOLDCV = 1;
	/** constant for score type: Cumulative-CV */
	final static int CUMCV = 2;

	/** the score types **/
	public static final Tag[] TAGS_CV_TYPE =
		{
			new Tag(LOOCV, "LOO-CV"),
			new Tag(KFOLDCV, "k-Fold-CV"),
			new Tag(CUMCV, "Cumulative-CV")
		};
	/**
	 * Holds the cross validation strategy used to measure quality of network
	 */
	int m_nCVType = LOOCV;

	/**
	 * performCV returns the accuracy calculated using cross validation.  
	 * The dataset used is m_Instances associated with the Bayes Network.
	 * 
	 * @param bayesNet : Bayes Network containing structure to evaluate
	 * @return accuracy (in interval 0..1) measured using cv.
	 * @throws Exception whn m_nCVType is invalided + exceptions passed on by updateClassifier
	 */
	public double calcScore(BayesNet bayesNet) throws Exception {
		switch (m_nCVType) {
			case LOOCV: 
				return leaveOneOutCV(bayesNet);
			case CUMCV: 
				return cumulativeCV(bayesNet);
			case KFOLDCV: 
				return kFoldCV(bayesNet, m_nNrOfFolds);
			default:
				throw new Exception("Unrecognized cross validation type encountered: " + m_nCVType);
		}
	} // calcScore

	/**
	 * Calc Node Score With Added Parent
	 * 
	 * @param nNode node for which the score is calculate
	 * @param nCandidateParent candidate parent to add to the existing parent set
	 * @return log score
	 * @throws Exception if something goes wrong
	 */
	public double calcScoreWithExtraParent(int nNode, int nCandidateParent) throws Exception {
		ParentSet oParentSet = m_BayesNet.getParentSet(nNode);
		Instances instances = m_BayesNet.m_Instances;

		// sanity check: nCandidateParent should not be in parent set already
		for (int iParent = 0; iParent < oParentSet.getNrOfParents(); iParent++) {
			if (oParentSet.getParent(iParent) == nCandidateParent) {
				return -1e100;
			}
		}

		// set up candidate parent
		oParentSet.addParent(nCandidateParent, instances);

		// calculate the score
		double fAccuracy = calcScore(m_BayesNet);

		// delete temporarily added parent
		oParentSet.deleteLastParent(instances);

		return fAccuracy;
	} // calcScoreWithExtraParent


	/**
	 * Calc Node Score With Parent Deleted
	 * 
	 * @param nNode node for which the score is calculate
	 * @param nCandidateParent candidate parent to delete from the existing parent set
	 * @return log score
	 * @throws Exception if something goes wrong
	 */
	public double calcScoreWithMissingParent(int nNode, int nCandidateParent) throws Exception {
		ParentSet oParentSet = m_BayesNet.getParentSet(nNode);
		Instances instances = m_BayesNet.m_Instances;

		// sanity check: nCandidateParent should be in parent set already
		if (!oParentSet.contains( nCandidateParent)) {
				return -1e100;
		}

		// set up candidate parent
		int iParent = oParentSet.deleteParent(nCandidateParent, instances);

		// calculate the score
		double fAccuracy = calcScore(m_BayesNet);

		// reinsert temporarily deleted parent
		oParentSet.addParent(nCandidateParent, iParent, instances);

		return fAccuracy;
	} // calcScoreWithMissingParent

	/**
	 * Calc Node Score With Arrow reversed
	 * 
	 * @param nNode node for which the score is calculate
	 * @param nCandidateParent candidate parent to delete from the existing parent set
	 * @return log score
	 * @throws Exception if something goes wrong
	 */
	public double calcScoreWithReversedParent(int nNode, int nCandidateParent) throws Exception {
		ParentSet oParentSet = m_BayesNet.getParentSet(nNode);
		ParentSet oParentSet2 = m_BayesNet.getParentSet(nCandidateParent);
		Instances instances = m_BayesNet.m_Instances;

		// sanity check: nCandidateParent should be in parent set already
		if (!oParentSet.contains( nCandidateParent)) {
				return -1e100;
		}

		// set up candidate parent
		int iParent = oParentSet.deleteParent(nCandidateParent, instances);
		oParentSet2.addParent(nNode, instances);

		// calculate the score
		double fAccuracy = calcScore(m_BayesNet);

		// restate temporarily reversed arrow
		oParentSet2.deleteLastParent(instances);
		oParentSet.addParent(nCandidateParent, iParent, instances);

		return fAccuracy;
	} // calcScoreWithReversedParent

	/**
	 * LeaveOneOutCV returns the accuracy calculated using Leave One Out
	 * cross validation. The dataset used is m_Instances associated with
	 * the Bayes Network.
	 * @param bayesNet : Bayes Network containing structure to evaluate
	 * @return accuracy (in interval 0..1) measured using leave one out cv.
	 * @throws Exception passed on by updateClassifier
	 */
	public double leaveOneOutCV(BayesNet bayesNet) throws Exception {
		m_BayesNet = bayesNet;
		double fAccuracy = 0.0;
		double fWeight = 0.0;
		Instances instances = bayesNet.m_Instances;
		bayesNet.estimateCPTs();
		for (int iInstance = 0; iInstance < instances.numInstances(); iInstance++) {
			Instance instance = instances.instance(iInstance);
			instance.setWeight(-instance.weight());
			bayesNet.updateClassifier(instance);
			fAccuracy += accuracyIncrease(instance);
			fWeight += instance.weight();
			instance.setWeight(-instance.weight());
			bayesNet.updateClassifier(instance);
		}
		return fAccuracy / fWeight;
	} // LeaveOneOutCV

	/**
	 * CumulativeCV returns the accuracy calculated using cumulative
	 * cross validation. The idea is to run through the data set and
	 * try to classify each of the instances based on the previously
	 * seen data.
	 * The data set used is m_Instances associated with the Bayes Network.
	 * @param bayesNet : Bayes Network containing structure to evaluate
	 * @return accuracy (in interval 0..1) measured using leave one out cv.
	 * @throws Exception passed on by updateClassifier
	 */
	public double cumulativeCV(BayesNet bayesNet) throws Exception {
		m_BayesNet = bayesNet;
		double fAccuracy = 0.0;
		double fWeight = 0.0;
		Instances instances = bayesNet.m_Instances;
		bayesNet.initCPTs();
		for (int iInstance = 0; iInstance < instances.numInstances(); iInstance++) {
			Instance instance = instances.instance(iInstance);
			fAccuracy += accuracyIncrease(instance);
			bayesNet.updateClassifier(instance);
			fWeight += instance.weight();
		}
		return fAccuracy / fWeight;
	} // LeaveOneOutCV
	
	/**
	 * kFoldCV uses k-fold cross validation to measure the accuracy of a Bayes
	 * network classifier.
	 * @param bayesNet : Bayes Network containing structure to evaluate
	 * @param nNrOfFolds : the number of folds k to perform k-fold cv
	 * @return accuracy (in interval 0..1) measured using leave one out cv.
	 * @throws Exception passed on by updateClassifier
	 */
	public double kFoldCV(BayesNet bayesNet, int nNrOfFolds) throws Exception {
		m_BayesNet = bayesNet;
		double fAccuracy = 0.0;
		double fWeight = 0.0;
		Instances instances = bayesNet.m_Instances;
		// estimate CPTs based on complete data set
		bayesNet.estimateCPTs();
		int nFoldStart = 0;
		int nFoldEnd = instances.numInstances() / nNrOfFolds;
		int iFold = 1;
		while (nFoldStart < instances.numInstances()) {
			// remove influence of fold iFold from the probability distribution
			for (int iInstance = nFoldStart; iInstance < nFoldEnd; iInstance++) {
				Instance instance = instances.instance(iInstance);
				instance.setWeight(-instance.weight());
				bayesNet.updateClassifier(instance);
			}
			
			// measure accuracy on fold iFold
			for (int iInstance = nFoldStart; iInstance < nFoldEnd; iInstance++) {
				Instance instance = instances.instance(iInstance);
				instance.setWeight(-instance.weight());
				fAccuracy += accuracyIncrease(instance);
				instance.setWeight(-instance.weight());
				fWeight += instance.weight();
			}

			// restore influence of fold iFold from the probability distribution
			for (int iInstance = nFoldStart; iInstance < nFoldEnd; iInstance++) {
				Instance instance = instances.instance(iInstance);
				instance.setWeight(-instance.weight());
				bayesNet.updateClassifier(instance);
			}

			// go to next fold
			nFoldStart = nFoldEnd;
			iFold++;
			nFoldEnd = iFold * instances.numInstances() / nNrOfFolds;
		}
		return fAccuracy / fWeight;
	} // kFoldCV
	
	/** accuracyIncrease determines how much the accuracy estimate should
	 * be increased due to the contribution of a single given instance. 
	 * 
	 * @param instance : instance for which to calculate the accuracy increase.
	 * @return increase in accuracy due to given instance. 
	 * @throws Exception passed on by distributionForInstance and classifyInstance
	 */
	double accuracyIncrease(Instance instance) throws Exception {
		if (m_bUseProb) {
			double [] fProb = m_BayesNet.distributionForInstance(instance);
			return fProb[(int) instance.classValue()] * instance.weight();
		} else {
			if (m_BayesNet.classifyInstance(instance) == instance.classValue()) {
				return instance.weight();
			}
		}
		return 0;
	} // accuracyIncrease

	/**
	 * @return use probabilities or not in accuracy estimate
	 */
	public boolean getUseProb() {
		return m_bUseProb;
	} // getUseProb

	/**
	 * @param useProb : use probabilities or not in accuracy estimate
	 */
	public void setUseProb(boolean useProb) {
		m_bUseProb = useProb;
	} // setUseProb
	
	/**
	 * set cross validation strategy to be used in searching for networks.
	 * @param newCVType : cross validation strategy 
	 */
	public void setCVType(SelectedTag newCVType) {
		if (newCVType.getTags() == TAGS_CV_TYPE) {
			m_nCVType = newCVType.getSelectedTag().getID();
		}
	} // setCVType

	/**
	 * get cross validation strategy to be used in searching for networks.
	 * @return cross validation strategy 
	 */
	public SelectedTag getCVType() {
		return new SelectedTag(m_nCVType, TAGS_CV_TYPE);
	} // getCVType

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
				"\tScore type (LOO-CV,k-Fold-CV,Cumulative-CV)",
				"S",
				1,
				"-S [LOO-CV|k-Fold-CV|Cumulative-CV]"));

		newVector.addElement(new Option("\tUse probabilistic or 0/1 scoring.\n\t(default probabilistic scoring)", "Q", 0, "-Q"));

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

	  	setMarkovBlanketClassifier(Utils.getFlag("mbc", options));

		String sScore = Utils.getOption('S', options);

		if (sScore.compareTo("LOO-CV") == 0) {
			setCVType(new SelectedTag(LOOCV, TAGS_CV_TYPE));
		}
		if (sScore.compareTo("k-Fold-CV") == 0) {
			setCVType(new SelectedTag(KFOLDCV, TAGS_CV_TYPE));
		}
		if (sScore.compareTo("Cumulative-CV") == 0) {
			setCVType(new SelectedTag(CUMCV, TAGS_CV_TYPE));
		}
		setUseProb(!Utils.getFlag('Q', options));		
		super.setOptions(options);
	} // setOptions

	/**
	 * Gets the current settings of the search algorithm.
	 *
	 * @return an array of strings suitable for passing to setOptions
	 */
	public String[] getOptions() {
		String[] superOptions = super.getOptions();
		String[] options = new String[4 + superOptions.length];
		int current = 0;

		if (getMarkovBlanketClassifier())
		  options[current++] = "-mbc";

		options[current++] = "-S";

		switch (m_nCVType) {
			case (LOOCV) :
				options[current++] = "LOO-CV";
				break;
			case (KFOLDCV) :
				options[current++] = "k-Fold-CV";
				break;
			case (CUMCV) :
				options[current++] = "Cumulative-CV";
				break;
		}
		
		if (!getUseProb()) {
		  options[current++] = "-Q";
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
	 * @return a string to describe the CVType option.
	 */
	public String CVTypeTipText() {
	  return "Select cross validation strategy to be used in searching for networks." +
	  "LOO-CV = Leave one out cross validation\n" +
	  "k-Fold-CV = k fold cross validation\n" +
	  "Cumulative-CV = cumulative cross validation."
	  ;
	} // CVTypeTipText

	/**
	 * @return a string to describe the UseProb option.
	 */
	public String useProbTipText() {
	  return "If set to true, the probability of the class if returned in the estimate of the "+
	  "accuracy. If set to false, the accuracy estimate is only increased if the classifier returns " +
	  "exactly the correct class.";
	} // useProbTipText

	/**
	 * This will return a string describing the search algorithm.
	 * @return The string.
	 */
	public String globalInfo() {
	  return "This Bayes Network learning algorithm uses cross validation to estimate " +
	  "classification accuracy.";
	} // globalInfo
	
	/**
	 * @return a string to describe the MarkovBlanketClassifier option.
	 */
	public String markovBlanketClassifierTipText() {
	  return super.markovBlanketClassifierTipText();
	}

	/**
	 * Returns the revision string.
	 * 
	 * @return		the revision
	 */
	public String getRevision() {
	  return RevisionUtils.extract("$Revision: 1.10 $");
	}
}
