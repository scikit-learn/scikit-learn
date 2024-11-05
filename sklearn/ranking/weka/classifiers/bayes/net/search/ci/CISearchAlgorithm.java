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
 * CISearchAlgorithm.java
 * Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 * 
 */

package weka.classifiers.bayes.net.search.ci;

import weka.classifiers.bayes.BayesNet;
import weka.classifiers.bayes.net.ParentSet;
import weka.classifiers.bayes.net.search.local.LocalScoreSearchAlgorithm;
import weka.core.Instances;
import weka.core.RevisionUtils;

/** 
 <!-- globalinfo-start -->
 * The CISearchAlgorithm class supports Bayes net structure search algorithms that are based on conditional independence test (as opposed to for example score based of cross validation based search algorithms).
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
 * @author Remco Bouckaert (rrb@xm.co.nz)
 * @version $Revision: 1.7 $
 */
public class CISearchAlgorithm 
	extends LocalScoreSearchAlgorithm {
  	
  	/** for serialization */
  	static final long serialVersionUID = 3165802334119704560L;
  	
	BayesNet  m_BayesNet;
	Instances m_instances;
	    
	/**
	 * Returns a string describing this object
	 * @return a description of the classifier suitable for
	 * displaying in the explorer/experimenter gui
	 */
	public String globalInfo() {
	  return 
	      "The CISearchAlgorithm class supports Bayes net structure "
	    + "search algorithms that are based on conditional independence "
	    + "test (as opposed to for example score based of cross validation "
	    + "based search algorithms).";
	}

	/** IsConditionalIndependent tests whether two nodes X and Y are independent
	 *  given a set of variables Z. The test compares the score of the Bayes network
	 * with and without arrow Y->X where all nodes in Z are parents of X.
	 * @param iAttributeX - index of attribute representing variable X
	 * @param iAttributeY - index of attribute representing variable Y
 	 * @param iAttributesZ - array of integers representing indices of attributes in set Z
 	 * @param nAttributesZ - cardinality of Z
 	 * @return true if X and Y conditionally independent given Z 
	 */
	protected boolean isConditionalIndependent(
		int iAttributeX, 
		int iAttributeY, 
		int [] iAttributesZ, 
		int nAttributesZ) {
		ParentSet oParentSetX = m_BayesNet.getParentSet(iAttributeX);
		// clear parent set of AttributeX
		while (oParentSetX.getNrOfParents() > 0) {
			oParentSetX.deleteLastParent(m_instances);
		}
		
		// insert parents in iAttributeZ
		for (int iAttributeZ = 0; iAttributeZ < nAttributesZ; iAttributeZ++) {
			oParentSetX.addParent( iAttributesZ[iAttributeZ], m_instances);
		}
		
		double fScoreZ = calcNodeScore(iAttributeX);
		double fScoreZY = calcScoreWithExtraParent(iAttributeX, iAttributeY);
		if (fScoreZY <= fScoreZ) {
			// the score does not improve by adding Y to the parent set of X
			// so we conclude that nodes X and Y are conditionally independent
			// given the set of variables Z
			return true;
		}
		return false;
	} // IsConditionalIndependent

	/**
	 * Returns the revision string.
	 * 
	 * @return		the revision
	 */
	public String getRevision() {
	  return RevisionUtils.extract("$Revision: 1.7 $");
	}
} // class CISearchAlgorithm
