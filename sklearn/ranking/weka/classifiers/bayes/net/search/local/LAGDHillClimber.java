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
 * LAGDHillClimber.java
 * Copyright (C) 2005 Manuel Neubach
 * 
 */

package weka.classifiers.bayes.net.search.local;

import weka.classifiers.bayes.BayesNet;
import weka.core.Instances;
import weka.core.Option;
import weka.core.RevisionUtils;
import weka.core.Utils;

import java.util.Enumeration;
import java.util.Vector;

/** 
 <!-- globalinfo-start -->
 * This Bayes Network learning algorithm uses a Look Ahead Hill Climbing algorithm called LAGD Hill Climbing. Unlike Greedy Hill Climbing it doesn't calculate a best greedy operation (adding, deleting or reversing an arc) but a sequence of nrOfLookAheadSteps operations, which leads to a network structure whose score is most likely higher in comparison to the network obtained by performing a sequence of nrOfLookAheadSteps greedy operations. The search is not restricted by an order on the variables (unlike K2). The difference with B and B2 is that this hill climber also considers arrows part of the naive Bayes structure for deletion.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -L &lt;nr of look ahead steps&gt;
 *  Look Ahead Depth</pre>
 * 
 * <pre> -G &lt;nr of good operations&gt;
 *  Nr of Good Operations</pre>
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
 * <pre> -S [BAYES|MDL|ENTROPY|AIC|CROSS_CLASSIC|CROSS_BAYES]
 *  Score type (BAYES, BDeu, MDL, ENTROPY and AIC)</pre>
 * 
 <!-- options-end -->
 * 
 * @author Manuel Neubach
 * @version $Revision: 1.7 $
 */
public class LAGDHillClimber 
    extends HillClimber {
  
    /** for serialization */
    static final long serialVersionUID = 7217437499439184344L;

    /** Number of Look Ahead Steps **/
    int m_nNrOfLookAheadSteps = 2;

    /** Number of Good Operations per Step **/
    int m_nNrOfGoodOperations = 5;

   /**
     * search determines the network structure/graph of the network
     * 
     * @param bayesNet the network
     * @param instances the data to use
     * @throws Exception if something goes wrong
     */
   protected void search(BayesNet bayesNet, Instances instances) throws Exception {
        int k=m_nNrOfLookAheadSteps;  // Number of Look Ahead Steps
        int l=m_nNrOfGoodOperations; // Number of Good Operations per step
        lookAheadInGoodDirectionsSearch(bayesNet, instances, k, l);
   } // search


   /**
    * lookAheadInGoodDirectionsSearch determines the network structure/graph of the network
    * with best score according to LAGD Hill Climbing
    * 
    * @param bayesNet the network
    * @param instances the data to use
    * @param nrOfLookAheadSteps
    * @param nrOfGoodOperations
    * @throws Exception if something goes wrong
    */
    protected void lookAheadInGoodDirectionsSearch(BayesNet bayesNet, Instances instances, int nrOfLookAheadSteps, int nrOfGoodOperations) throws Exception {
         System.out.println("Initializing Cache");
         initCache(bayesNet, instances);

         while (nrOfLookAheadSteps>1) {         
            System.out.println("Look Ahead Depth: "+nrOfLookAheadSteps);
            boolean legalSequence = true;
            double sequenceDeltaScore = 0;
            Operation [] bestOperation=new Operation [nrOfLookAheadSteps];
         
            bestOperation = getOptimalOperations(bayesNet, instances, nrOfLookAheadSteps, nrOfGoodOperations);
            for (int i = 0; i < nrOfLookAheadSteps; i++) {
               if (bestOperation [i] == null) {
                  legalSequence=false;
               } else {
                  sequenceDeltaScore += bestOperation [i].m_fDeltaScore;
               }
            }
            while (legalSequence && sequenceDeltaScore > 0) {
               System.out.println("Next Iteration..........................");
               for (int i = 0; i < nrOfLookAheadSteps; i++) {
                  performOperation(bayesNet, instances,bestOperation [i]);
               }
               bestOperation = getOptimalOperations(bayesNet, instances, nrOfLookAheadSteps, nrOfGoodOperations);
               sequenceDeltaScore = 0;
               for (int i = 0; i < nrOfLookAheadSteps; i++) {
                  if (bestOperation [i] != null) {
                     System.out.println(bestOperation [i].m_nOperation + " " + bestOperation [i].m_nHead + " " + bestOperation [i].m_nTail);
                     sequenceDeltaScore += bestOperation [i].m_fDeltaScore;
                  } else {
                     legalSequence = false;

                  }
                  System.out.println("DeltaScore: "+sequenceDeltaScore);
               }
            }
            --nrOfLookAheadSteps;
         }

         /** last steps with greedy HC **/          
         Operation oOperation = getOptimalOperation(bayesNet, instances);
         while ((oOperation != null) && (oOperation.m_fDeltaScore > 0)) {
	    performOperation(bayesNet, instances, oOperation);
            System.out.println("Performing last greedy steps");
            oOperation = getOptimalOperation(bayesNet, instances);
         }               
	 // free up memory
	 m_Cache = null;
    } // lookAheadInGoodDirectionsSearch

    /**
      * getAntiOperation determines the Operation, which is needed to cancel oOperation
      * 
      * @param oOperation Operation to cancel
      * @return antiOperation to oOperation
      * @throws Exception if something goes wrong
      */
    protected Operation getAntiOperation(Operation oOperation) throws Exception {
        if (oOperation.m_nOperation == Operation.OPERATION_ADD)
           return (new Operation (oOperation.m_nTail, oOperation.m_nHead, Operation.OPERATION_DEL));
        else {
           if (oOperation.m_nOperation == Operation.OPERATION_DEL)
              return (new Operation (oOperation.m_nTail, oOperation.m_nHead, Operation.OPERATION_ADD));
           else {
              return (new Operation (oOperation.m_nHead, oOperation.m_nTail, Operation.OPERATION_REVERSE));
           }
         }
    } // getAntiOperation


    /**
      * getGoodOperations determines the nrOfGoodOperations best Operations, which are considered for
      * the calculation of an optimal operationsequence
      * @param bayesNet Bayes network to apply operation on
      * @param instances data set to learn from
      * @param nrOfGoodOperations number of good operations to consider
      * @return good operations to consider
      * @throws Exception if something goes wrong
      **/
    protected Operation [] getGoodOperations(BayesNet bayesNet, Instances instances, int nrOfGoodOperations) throws Exception {
		Operation [] goodOperations=new Operation [nrOfGoodOperations];
       		for (int i = 0; i < nrOfGoodOperations; i++) {
                   goodOperations [i] = getOptimalOperation(bayesNet, instances);
                   if (goodOperations[i] != null) {
                      m_Cache.put(goodOperations [i], -1E100);
                   } else i=nrOfGoodOperations;
                }
                for (int i = 0; i < nrOfGoodOperations; i++) {
                   if (goodOperations[i] != null) {
                      if (goodOperations [i].m_nOperation!=Operation.OPERATION_REVERSE) {
                         m_Cache.put(goodOperations [i], goodOperations [i].m_fDeltaScore);
                      } else {
                         m_Cache.put(goodOperations [i], goodOperations [i].m_fDeltaScore - m_Cache.m_fDeltaScoreAdd[goodOperations[i].m_nHead] [goodOperations [i].m_nTail]);
                      }
                   } else i=nrOfGoodOperations;
                }
                return goodOperations;
    } // getGoodOperations

    /**
      * getOptimalOperations determines an optimal operationsequence in respect of the parameters 
      * nrOfLookAheadSteps and nrOfGoodOperations
      * @param bayesNet Bayes network to apply operation on
      * @param instances data set to learn from
      * @param nrOfLookAheadSteps number of lood ahead steps to use
      * @param nrOfGoodOperations number of good operations to consider
      * @return optimal sequence of operations in respect to nrOfLookAheadSteps and nrOfGoodOperations
      * @throws Exception if something goes wrong
      **/
    protected Operation [] getOptimalOperations(BayesNet bayesNet, Instances instances, int nrOfLookAheadSteps, int nrOfGoodOperations) throws Exception {
       if (nrOfLookAheadSteps == 1) { // Abbruch der Rekursion
          Operation [] bestOperation = new Operation [1];
          bestOperation [0] = getOptimalOperation(bayesNet, instances);
          return(bestOperation);  // Abbruch der Rekursion
       } else {
          double bestDeltaScore = 0;
          double currentDeltaScore = 0;
          Operation [] bestOperation = new Operation [nrOfLookAheadSteps];
          Operation [] goodOperations = new Operation [nrOfGoodOperations];
          Operation [] tempOperation = new Operation [nrOfLookAheadSteps-1];
          goodOperations = getGoodOperations(bayesNet, instances, nrOfGoodOperations);
          for (int i = 0; i < nrOfGoodOperations; i++) {
           if (goodOperations[i] != null) {
             performOperation(bayesNet, instances, goodOperations [i]);
             tempOperation = getOptimalOperations(bayesNet, instances, nrOfLookAheadSteps-1, nrOfGoodOperations); // rekursiver Abstieg
             currentDeltaScore = goodOperations [i].m_fDeltaScore;
             for (int j = 0; j < nrOfLookAheadSteps-1; j++) {
                if (tempOperation [j] != null) {
                   currentDeltaScore += tempOperation [j].m_fDeltaScore;
                }
             }
             performOperation(bayesNet, instances, getAntiOperation(goodOperations [i]));
                      if (currentDeltaScore > bestDeltaScore) {
                        bestDeltaScore = currentDeltaScore;
                        bestOperation [0] = goodOperations [i];
                        for (int j = 1; j < nrOfLookAheadSteps; j++) {
                          bestOperation [j] = tempOperation [j-1];
                        }
                    }
            } else i=nrOfGoodOperations;
           }
           return(bestOperation);
       }
    } // getOptimalOperations


	/**
	 * Sets the max number of parents
	 *
	 * @param nMaxNrOfParents the max number of parents
	 */
	public void setMaxNrOfParents(int nMaxNrOfParents) {
	  m_nMaxNrOfParents = nMaxNrOfParents;
	} 

	/**
	 * Gets the max number of parents.
	 *
	 * @return the max number of parents
	 */
	public int getMaxNrOfParents() {
	  return m_nMaxNrOfParents;
	} 

	/**
	 * Sets the number of look-ahead steps
	 *
	 * @param nNrOfLookAheadSteps the number of look-ahead steps
	 */
	public void setNrOfLookAheadSteps(int nNrOfLookAheadSteps) {
	  m_nNrOfLookAheadSteps = nNrOfLookAheadSteps;
	} 

	/**
	 * Gets the number of look-ahead steps
	 *
	 * @return the number of look-ahead step
	 */
	public int getNrOfLookAheadSteps() {
	  return m_nNrOfLookAheadSteps;
	} 

	/**
	 * Sets the number of "good operations"
	 *
	 * @param nNrOfGoodOperations the number of "good operations"
	 */
	public void setNrOfGoodOperations(int nNrOfGoodOperations) {
	  m_nNrOfGoodOperations = nNrOfGoodOperations;
	} 

	/**
	 * Gets the number of "good operations"
	 *
	 * @return the number of "good operations"
	 */
	public int getNrOfGoodOperations() {
	  return m_nNrOfGoodOperations;
	} 


	/**
	 * Returns an enumeration describing the available options.
	 *
	 * @return an enumeration of all the available options.
	 */
	public Enumeration listOptions() {
		Vector newVector = new Vector();

		newVector.addElement(new Option("\tLook Ahead Depth", "L", 2, "-L <nr of look ahead steps>"));
		newVector.addElement(new Option("\tNr of Good Operations", "G", 5, "-G <nr of good operations>"));

		Enumeration enm = super.listOptions();
		while (enm.hasMoreElements()) {
			newVector.addElement(enm.nextElement());
		}
		return newVector.elements();
	} // listOptions

	/**
	 * Parses a given list of options. Valid options are:<p>
	 *
	 <!-- options-start -->
	 * Valid options are: <p/>
	 * 
	 * <pre> -L &lt;nr of look ahead steps&gt;
	 *  Look Ahead Depth</pre>
	 * 
	 * <pre> -G &lt;nr of good operations&gt;
	 *  Nr of Good Operations</pre>
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
	 * <pre> -S [BAYES|MDL|ENTROPY|AIC|CROSS_CLASSIC|CROSS_BAYES]
	 *  Score type (BAYES, BDeu, MDL, ENTROPY and AIC)</pre>
	 * 
	 <!-- options-end -->
	 *
	 * @param options the list of options as an array of strings
	 * @throws Exception if an option is not supported
	 */
	public void setOptions(String[] options) throws Exception {
	  	String sNrOfLookAheadSteps = Utils.getOption('L', options);
		if (sNrOfLookAheadSteps.length() != 0) {
		  setNrOfLookAheadSteps(Integer.parseInt(sNrOfLookAheadSteps));
		} else {
		  setNrOfLookAheadSteps(2);
		}

                String sNrOfGoodOperations = Utils.getOption('G', options);
		if (sNrOfGoodOperations.length() != 0) {
		  setNrOfGoodOperations(Integer.parseInt(sNrOfGoodOperations));
		} else {
		  setNrOfGoodOperations(5);
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
		String[] options = new String[9 + superOptions.length];
		int current = 0;

		options[current++] = "-L";
		options[current++] = "" + m_nNrOfLookAheadSteps;
		
		options[current++] = "-G";
		options[current++] = "" + m_nNrOfGoodOperations;

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
	 * This will return a string describing the search algorithm.
	 * @return The string.
	 */
	public String globalInfo() {
	  return "This Bayes Network learning algorithm uses a Look Ahead Hill Climbing algorithm called LAGD Hill Climbing." +
	  " Unlike Greedy Hill Climbing it doesn't calculate a best greedy operation (adding, deleting or reversing an arc) " + 
	  "but a sequence of nrOfLookAheadSteps operations, which leads to a network structure whose score is most likely " +
	  "higher in comparison to the network obtained by performing a sequence of nrOfLookAheadSteps greedy operations. " +
	  "The search is not restricted by an order " +
	  "on the variables (unlike K2). The difference with B and B2 is that this hill " +
          "climber also considers arrows part of the naive Bayes structure for deletion.";
	} // globalInfo

	/**
	 * @return a string to describe the Number of Look Ahead Steps option.
	 */
	public String nrOfLookAheadStepsTipText() {
	  return "Sets the Number of Look Ahead Steps. 'nrOfLookAheadSteps = 2' means that all network structures in a " +
	  "distance of 2 (from the current network structure) are taken into account for the decision which arcs to add, " +
	  "remove or reverse. 'nrOfLookAheadSteps = 1' results in Greedy Hill Climbing." ;
	} // nrOfLookAheadStepsTipText

	/**
	 * @return a string to describe the Number of Good Operations option.
	 */
	public String nrOfGoodOperationsTipText() {
	  return "Sets the Number of Good Operations per Look Ahead Step. 'nrOfGoodOperations = 5' means that for the next " +
	  "Look Ahead Step only the 5 best Operations (adding, deleting or reversing an arc) are taken into account for the " +
	  "calculation of the best sequence consisting of nrOfLookAheadSteps operations." ;
	} // nrOfGoodOperationsTipText

	/**
	 * Returns the revision string.
	 * 
	 * @return		the revision
	 */
	public String getRevision() {
	  return RevisionUtils.extract("$Revision: 1.7 $");
	}

} // LAGDHillClimber
