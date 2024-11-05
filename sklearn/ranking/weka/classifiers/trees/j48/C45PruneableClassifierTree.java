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
 *    C45PruneableClassifierTree.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.trees.j48;

import weka.core.Capabilities;
import weka.core.Instances;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.core.Capabilities.Capability;

/**
 * Class for handling a tree structure that can
 * be pruned using C4.5 procedures.
 *
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version $Revision: 6073 $
 */

public class C45PruneableClassifierTree 
  extends ClassifierTree {

  /** for serialization */
  static final long serialVersionUID = -4813820170260388194L;
  
  /** True if the tree is to be pruned. */
  boolean m_pruneTheTree = false;
  
  /** True if the tree is to be collapsed. */
  boolean m_collapseTheTree = false;

  /** The confidence factor for pruning. */
  float m_CF = 0.25f;

  /** Is subtree raising to be performed? */
  boolean m_subtreeRaising = true;

  /** Cleanup after the tree has been built. */
  boolean m_cleanup = true;

  /**
   * Constructor for pruneable tree structure. Stores reference
   * to associated training data at each node.
   *
   * @param toSelectLocModel selection method for local splitting model
   * @param pruneTree true if the tree is to be pruned
   * @param cf the confidence factor for pruning
   * @param raiseTree
   * @param cleanup
   * @throws Exception if something goes wrong
   */
  public C45PruneableClassifierTree(ModelSelection toSelectLocModel,
				    boolean pruneTree,float cf,
				    boolean raiseTree,
				    boolean cleanup,
                                    boolean collapseTree)
       throws Exception {

    super(toSelectLocModel);

    m_pruneTheTree = pruneTree;
    m_CF = cf;
    m_subtreeRaising = raiseTree;
    m_cleanup = cleanup;
    m_collapseTheTree = collapseTree;
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
   * @param data the data for building the tree
   * @throws Exception if something goes wrong
   */
  public void buildClassifier(Instances data) throws Exception {

    // can classifier tree handle the data?
    getCapabilities().testWithFail(data);

    // remove instances with missing class
    data = new Instances(data);
    data.deleteWithMissingClass();
    
   buildTree(data, m_subtreeRaising);
   if (m_collapseTheTree) {
     collapse();
   }
   if (m_pruneTheTree) {
     prune();
   }
   if (m_cleanup) {
     cleanup(new Instances(data, 0));
   }
  }

  /**
   * Collapses a tree to a node if training error doesn't increase.
   */
  public final void collapse(){

    double errorsOfSubtree;
    double errorsOfTree;
    int i;

    if (!m_isLeaf){
      errorsOfSubtree = getTrainingErrors();
      errorsOfTree = localModel().distribution().numIncorrect();
      if (errorsOfSubtree >= errorsOfTree-1E-3){

	// Free adjacent trees
	m_sons = null;
	m_isLeaf = true;
			
	// Get NoSplit Model for tree.
	m_localModel = new NoSplit(localModel().distribution());
      }else
	for (i=0;i<m_sons.length;i++)
	  son(i).collapse();
    }
  }

  /**
   * Prunes a tree using C4.5's pruning procedure.
   *
   * @throws Exception if something goes wrong
   */
  public void prune() throws Exception {

    double errorsLargestBranch;
    double errorsLeaf;
    double errorsTree;
    int indexOfLargestBranch;
    C45PruneableClassifierTree largestBranch;
    int i;

    if (!m_isLeaf){

      // Prune all subtrees.
      for (i=0;i<m_sons.length;i++)
	son(i).prune();

      // Compute error for largest branch
      indexOfLargestBranch = localModel().distribution().maxBag();
      if (m_subtreeRaising) {
	errorsLargestBranch = son(indexOfLargestBranch).
	  getEstimatedErrorsForBranch((Instances)m_train);
      } else {
	errorsLargestBranch = Double.MAX_VALUE;
      }

      // Compute error if this Tree would be leaf
      errorsLeaf = 
	getEstimatedErrorsForDistribution(localModel().distribution());

      // Compute error for the whole subtree
      errorsTree = getEstimatedErrors();

      // Decide if leaf is best choice.
      if (Utils.smOrEq(errorsLeaf,errorsTree+0.1) &&
	  Utils.smOrEq(errorsLeaf,errorsLargestBranch+0.1)){

	// Free son Trees
	m_sons = null;
	m_isLeaf = true;
		
	// Get NoSplit Model for node.
	m_localModel = new NoSplit(localModel().distribution());
	return;
      }

      // Decide if largest branch is better choice
      // than whole subtree.
      if (Utils.smOrEq(errorsLargestBranch,errorsTree+0.1)){
	largestBranch = son(indexOfLargestBranch);
	m_sons = largestBranch.m_sons;
	m_localModel = largestBranch.localModel();
	m_isLeaf = largestBranch.m_isLeaf;
	newDistribution(m_train);
	prune();
      }
    }
  }

  /**
   * Returns a newly created tree.
   *
   * @param data the data to work with
   * @return the new tree
   * @throws Exception if something goes wrong
   */
  protected ClassifierTree getNewTree(Instances data) throws Exception {
    
    C45PruneableClassifierTree newTree = 
      new C45PruneableClassifierTree(m_toSelectModel, m_pruneTheTree, m_CF,
				     m_subtreeRaising, m_cleanup, m_collapseTheTree);
    newTree.buildTree((Instances)data, m_subtreeRaising);

    return newTree;
  }

  /**
   * Computes estimated errors for tree.
   * 
   * @return the estimated errors
   */
  private double getEstimatedErrors(){

    double errors = 0;
    int i;

    if (m_isLeaf)
      return getEstimatedErrorsForDistribution(localModel().distribution());
    else{
      for (i=0;i<m_sons.length;i++)
	errors = errors+son(i).getEstimatedErrors();
      return errors;
    }
  }
  
  /**
   * Computes estimated errors for one branch.
   *
   * @param data the data to work with
   * @return the estimated errors
   * @throws Exception if something goes wrong
   */
  private double getEstimatedErrorsForBranch(Instances data) 
       throws Exception {

    Instances [] localInstances;
    double errors = 0;
    int i;

    if (m_isLeaf)
      return getEstimatedErrorsForDistribution(new Distribution(data));
    else{
      Distribution savedDist = localModel().m_distribution;
      localModel().resetDistribution(data);
      localInstances = (Instances[])localModel().split(data);
      localModel().m_distribution = savedDist;
      for (i=0;i<m_sons.length;i++)
	errors = errors+
	  son(i).getEstimatedErrorsForBranch(localInstances[i]);
      return errors;
    }
  }

  /**
   * Computes estimated errors for leaf.
   * 
   * @param theDistribution the distribution to use
   * @return the estimated errors
   */
  private double getEstimatedErrorsForDistribution(Distribution 
						   theDistribution){

    if (Utils.eq(theDistribution.total(),0))
      return 0;
    else
      return theDistribution.numIncorrect()+
	Stats.addErrs(theDistribution.total(),
		      theDistribution.numIncorrect(),m_CF);
  }

  /**
   * Computes errors of tree on training data.
   * 
   * @return the training errors
   */
  private double getTrainingErrors(){

    double errors = 0;
    int i;

    if (m_isLeaf)
      return localModel().distribution().numIncorrect();
    else{
      for (i=0;i<m_sons.length;i++)
	errors = errors+son(i).getTrainingErrors();
      return errors;
    }
  }

  /**
   * Method just exists to make program easier to read.
   * 
   * @return the local split model
   */
  private ClassifierSplitModel localModel(){
    
    return (ClassifierSplitModel)m_localModel;
  }

  /**
   * Computes new distributions of instances for nodes
   * in tree.
   *
   * @param data the data to compute the distributions for
   * @throws Exception if something goes wrong
   */
  private void newDistribution(Instances data) throws Exception {

    Instances [] localInstances;

    localModel().resetDistribution(data);
    m_train = data;
    if (!m_isLeaf){
      localInstances = 
	(Instances [])localModel().split(data);
      for (int i = 0; i < m_sons.length; i++)
	son(i).newDistribution(localInstances[i]);
    } else {

      // Check whether there are some instances at the leaf now!
      if (!Utils.eq(data.sumOfWeights(), 0)) {
	m_isEmpty = false;
      }
    }
  }

  /**
   * Method just exists to make program easier to read.
   */
  private C45PruneableClassifierTree son(int index){

    return (C45PruneableClassifierTree)m_sons[index];
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6073 $");
  }
}
