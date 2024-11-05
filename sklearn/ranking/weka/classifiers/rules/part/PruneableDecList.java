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
 *    PruneableDecList.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.rules.part;

import weka.classifiers.trees.j48.Distribution;
import weka.classifiers.trees.j48.ModelSelection;
import weka.classifiers.trees.j48.NoSplit;
import weka.core.Instances;
import weka.core.RevisionUtils;
import weka.core.Utils;

/**
 * Class for handling a partial tree structure that
 * can be pruned using a pruning set.
 *
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version $Revision: 1.10 $
 */
public class PruneableDecList
  extends ClassifierDecList {

  /** for serialization */
  private static final long serialVersionUID = -7228103346297172921L;
  
  /**
   * Constructor for pruneable partial tree structure. 
   *
   * @param toSelectLocModel selection method for local splitting model
   * @param minNum minimum number of objects in leaf
   */
  public PruneableDecList(ModelSelection toSelectLocModel,
			  int minNum) {
			       
    super(toSelectLocModel, minNum);
  }
  
  /**
   * Method for building a pruned partial tree.
   *
   * @throws Exception if tree can't be built successfully
   */
  public void buildRule(Instances train,
			Instances test) throws Exception { 
    
    buildDecList(train, test, false);

    cleanup(new Instances(train, 0));
  }

  /**
   * Builds the partial tree with hold out set
   *
   * @throws Exception if something goes wrong
   */
  public void buildDecList(Instances train, Instances test, 
			   boolean leaf) throws Exception {
    
    Instances [] localTrain,localTest;
    int index,ind;
    int i,j;
    double sumOfWeights;
    NoSplit noSplit;
    
    m_train = null;
    m_isLeaf = false;
    m_isEmpty = false;
    m_sons = null;
    indeX = 0;
    sumOfWeights = train.sumOfWeights();
    noSplit = new NoSplit (new Distribution((Instances)train));
    if (leaf)
      m_localModel = noSplit;
    else
      m_localModel = m_toSelectModel.selectModel(train, test);
    m_test = new Distribution(test, m_localModel);
    if (m_localModel.numSubsets() > 1) {
      localTrain = m_localModel.split(train);
      localTest = m_localModel.split(test);
      train = null;
      test = null;
      m_sons = new ClassifierDecList [m_localModel.numSubsets()];
      i = 0;
      do {
	i++;
	ind = chooseIndex();
	if (ind == -1) {
	  for (j = 0; j < m_sons.length; j++) 
	    if (m_sons[j] == null)
	      m_sons[j] = getNewDecList(localTrain[j],localTest[j],true);
	  if (i < 2) {
	    m_localModel = noSplit;
	    m_isLeaf = true;
	    m_sons = null;
	    if (Utils.eq(sumOfWeights,0))
	      m_isEmpty = true;
	    return;
	  }
	  ind = 0;
	  break;
	} else 
	  m_sons[ind] = getNewDecList(localTrain[ind],localTest[ind],false);
      } while ((i < m_sons.length) && (m_sons[ind].m_isLeaf));
      
      // Check if all successors are leaves
      for (j = 0; j < m_sons.length; j++) 
	if ((m_sons[j] == null) || (!m_sons[j].m_isLeaf))
	  break;
      if (j == m_sons.length) {
	pruneEnd();
	if (!m_isLeaf) 
	  indeX = chooseLastIndex();
      }else 
	indeX = chooseLastIndex();
    }else{
      m_isLeaf = true;
      if (Utils.eq(sumOfWeights, 0))
	m_isEmpty = true;
    }
  }
  
  /**
   * Returns a newly created tree.
   *
   * @param train train data
   * @param test test data
   * @param leaf
   * @throws Exception if something goes wrong
   */
  protected ClassifierDecList getNewDecList(Instances train, Instances test, 
					    boolean leaf) throws Exception {
	 
    PruneableDecList newDecList = 
      new PruneableDecList(m_toSelectModel, m_minNumObj);
    
    newDecList.buildDecList((Instances)train, test, leaf);
    
    return newDecList;
  }

  /**
   * Prunes the end of the rule.
   */
  protected void pruneEnd() throws Exception {
    
    double errorsLeaf, errorsTree;
    
    errorsTree = errorsForTree();
    errorsLeaf = errorsForLeaf();
    if (Utils.smOrEq(errorsLeaf,errorsTree)){ 
      m_isLeaf = true;
      m_sons = null;
      m_localModel = new NoSplit(localModel().distribution());
    }
  }

  /**
   * Computes error estimate for tree.
   */
  private double errorsForTree() throws Exception {

    Distribution test;

    if (m_isLeaf)
      return errorsForLeaf();
    else {
      double error = 0;
      for (int i = 0; i < m_sons.length; i++) 
	if (Utils.eq(son(i).localModel().distribution().total(),0)) {
	  error += m_test.perBag(i)-
	    m_test.perClassPerBag(i,localModel().distribution().
				maxClass());
	} else
	  error += ((PruneableDecList)son(i)).errorsForTree();

      return error;
    }
  }

  /**
   * Computes estimated errors for leaf.
   */
  private double errorsForLeaf() throws Exception {

    return m_test.total()-
	    m_test.perClass(localModel().distribution().maxClass());
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
