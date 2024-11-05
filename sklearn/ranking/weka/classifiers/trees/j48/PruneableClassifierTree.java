/*
 *    This program is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program; if not, write to the Free Software
 *    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/*
 *    PruneableClassifierTree.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.trees.j48;

import weka.core.Capabilities;
import weka.core.Instances;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.core.Capabilities.Capability;

import java.util.Random;

/**
 * Class for handling a tree structure that can
 * be pruned using a pruning set. 
 *
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version $Revision: 5533 $
 */
public class PruneableClassifierTree 
  extends ClassifierTree {
  
  /** for serialization */
  static final long serialVersionUID = -555775736857600201L;

  /** True if the tree is to be pruned. */
  private boolean pruneTheTree = false;

  /** How many subsets of equal size? One used for pruning, the rest for training. */
  private int numSets = 3;

  /** Cleanup after the tree has been built. */
  private boolean m_cleanup = true;

  /** The random number seed. */
  private int m_seed = 1;

  /**
   * Constructor for pruneable tree structure. Stores reference
   * to associated training data at each node.
   *
   * @param toSelectLocModel selection method for local splitting model
   * @param pruneTree true if the tree is to be pruned
   * @param num number of subsets of equal size
   * @param cleanup
   * @param seed the seed value to use
   * @throws Exception if something goes wrong
   */
  public PruneableClassifierTree(ModelSelection toSelectLocModel,
				 boolean pruneTree, int num, boolean cleanup,
				 int seed)
       throws Exception {

    super(toSelectLocModel);

    pruneTheTree = pruneTree;
    numSets = num;
    m_cleanup = cleanup;
    m_seed = seed;
  }

  /**
   * Returns default capabilities of the classifier tree.
   *
   * @return      the capabilities of this classifier tree
   */
  public Capabilities getCapabilities() {
    Capabilities result = super.getCapabilities();
    result.disableAll();

    // attributes
    result.enable(Capability.NOMINAL_ATTRIBUTES);
    result.enable(Capability.NUMERIC_ATTRIBUTES);
    result.enable(Capability.DATE_ATTRIBUTES);
    result.enable(Capability.MISSING_VALUES);

    // class
    result.enable(Capability.NOMINAL_CLASS);
    result.enable(Capability.MISSING_CLASS_VALUES);

    // instances
    result.setMinimumNumberInstances(0);
    
    return result;
  }

  /**
   * Method for building a pruneable classifier tree.
   *
   * @param data the data to build the tree from 
   * @throws Exception if tree can't be built successfully
   */
  public void buildClassifier(Instances data) 
       throws Exception {

    // can classifier tree handle the data?
    getCapabilities().testWithFail(data);

    // remove instances with missing class
    data = new Instances(data);
    data.deleteWithMissingClass();
    
   Random random = new Random(m_seed);
   data.stratify(numSets);
   buildTree(data.trainCV(numSets, numSets - 1, random),
	     data.testCV(numSets, numSets - 1), false);
   if (pruneTheTree) {
     prune();
   }
   if (m_cleanup) {
     cleanup(new Instances(data, 0));
   }
  }

  /**
   * Prunes a tree.
   *
   * @throws Exception if tree can't be pruned successfully
   */
  public void prune() throws Exception {
  
    if (!m_isLeaf) {
      
      // Prune all subtrees.
      for (int i = 0; i < m_sons.length; i++)
	son(i).prune();
      
      // Decide if leaf is best choice.
      if (Utils.smOrEq(errorsForLeaf(),errorsForTree())) {
	
	// Free son Trees
	m_sons = null;
	m_isLeaf = true;
	
	// Get NoSplit Model for node.
	m_localModel = new NoSplit(localModel().distribution());
      }
    }
  }

  /**
   * Returns a newly created tree.
   *
   * @param train the training data
   * @param test the test data
   * @return the generated tree
   * @throws Exception if something goes wrong
   */
  protected ClassifierTree getNewTree(Instances train, Instances test) 
       throws Exception {

    PruneableClassifierTree newTree = 
      new PruneableClassifierTree(m_toSelectModel, pruneTheTree, numSets, m_cleanup,
				  m_seed);
    newTree.buildTree(train, test, false);
    return newTree;
  }

  /**
   * Computes estimated errors for tree.
   *
   * @return the estimated errors
   * @throws Exception if error estimate can't be computed
   */
  private double errorsForTree() throws Exception {

    double errors = 0;

    if (m_isLeaf)
      return errorsForLeaf();
    else{
      for (int i = 0; i < m_sons.length; i++)
	if (Utils.eq(localModel().distribution().perBag(i), 0)) {
	  errors += m_test.perBag(i)-
	    m_test.perClassPerBag(i,localModel().distribution().
				maxClass());
	} else
	  errors += son(i).errorsForTree();

      return errors;
    }
  }

  /**
   * Computes estimated errors for leaf.
   *
   * @return the estimated errors
   * @throws Exception if error estimate can't be computed
   */
  private double errorsForLeaf() throws Exception {

    return m_test.total()-
      m_test.perClass(localModel().distribution().maxClass());
  }

  /**
   * Method just exists to make program easier to read.
   */
  private ClassifierSplitModel localModel() {
    
    return (ClassifierSplitModel)m_localModel;
  }

  /**
   * Method just exists to make program easier to read.
   */
  private PruneableClassifierTree son(int index) {

    return (PruneableClassifierTree)m_sons[index];
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5533 $");
  }
}
