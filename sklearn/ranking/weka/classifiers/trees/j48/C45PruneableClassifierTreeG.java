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
 *    C45PruneableClassifierTreeG.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *    Copyright (C) 2007 Geoff Webb & Janice Boughton
 *
 */

package weka.classifiers.trees.j48;

import weka.core.Capabilities;
import weka.core.Instances;
import weka.core.Instance;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import java.util.ArrayList;
import java.util.Collections;

/**
 * Class for handling a tree structure that can
 * be pruned using C4.5 procedures and have nodes grafted on.
 *
 * @author Janice Boughton (based on code by Eibe Frank)
 * @version $Revision: 5532 $
 */

public class C45PruneableClassifierTreeG extends ClassifierTree{

  /** for serialization */
  static final long serialVersionUID = 66981207374331964L;

  /** True if the tree is to be pruned. */
  boolean m_pruneTheTree = false;

  /** The confidence factor for pruning. */
  float m_CF = 0.25f;

  /** Is subtree raising to be performed? */
  boolean m_subtreeRaising = true;

  /** Cleanup after the tree has been built. */
  boolean m_cleanup = true;

  /** flag for using relabelling when grafting */
  boolean m_relabel = false;

  /** binomial probability critical value */
  double m_BiProbCrit = 1.64;

  boolean m_Debug = false;

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
  public C45PruneableClassifierTreeG(ModelSelection toSelectLocModel,
				    boolean pruneTree,float cf,
				    boolean raiseTree,
				    boolean relabel, boolean cleanup)
       throws Exception {

    super(toSelectLocModel);

    m_pruneTheTree = pruneTree;
    m_CF = cf;
    m_subtreeRaising = raiseTree;
    m_cleanup = cleanup;
    m_relabel = relabel;
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
    result.enable(Capability.MISSING_VALUES);

    // class
    result.enable(Capability.NOMINAL_CLASS);
    result.enable(Capability.MISSING_CLASS_VALUES);

    // instances
    result.setMinimumNumberInstances(0);

    return result;
  }

  /**
   * Constructor for pruneable tree structure. Used to create new nodes
   * in the tree during grafting.
   *
   * @param toSelectLocModel selection method for local splitting model
   * @param data the dta used to produce split model
   * @param gs the split model
   * @param prune true if the tree is to be pruned
   * @param cf the confidence factor for pruning
   * @param raise
   * @param isLeaf if this node is a leaf or not
   * @param relabel whether relabeling occured
   * @param cleanup
   * @throws Exception if something goes wrong
   */
  public C45PruneableClassifierTreeG(ModelSelection toSelectLocModel, 
                                    Instances data, ClassifierSplitModel gs, 
                                    boolean prune, float cf, boolean raise,
                                    boolean isLeaf, boolean relabel, 
                                    boolean cleanup) {

    super(toSelectLocModel);
    m_relabel = relabel;
    m_cleanup = cleanup;
    m_localModel = gs;
    m_train = data;
    m_test = null;
    m_isLeaf = isLeaf;
    if(gs.distribution().total() > 0)
       m_isEmpty = false;
    else
       m_isEmpty = true;

    m_pruneTheTree = prune;
    m_CF = cf;
    m_subtreeRaising = raise;
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
    collapse();
    if (m_pruneTheTree) {
      prune();
    }
    doGrafting(data);
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
    C45PruneableClassifierTreeG largestBranch;
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
    
    C45PruneableClassifierTreeG newTree = 
      new C45PruneableClassifierTreeG(m_toSelectModel, m_pruneTheTree, m_CF,
	     m_subtreeRaising, m_relabel, m_cleanup);
	// ATBOP Modification     // m_subtreeRaising, m_cleanup);

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
  private C45PruneableClassifierTreeG son(int index){
    return (C45PruneableClassifierTreeG)m_sons[index];
  }


  /**
   * Initializes variables for grafting.
   * sets up limits array (for numeric attributes) and calls 
   * the recursive function traverseTree.
   *
   * @param data the data for the tree
   * @throws Exception if anything goes wrong
   */
  public void doGrafting(Instances data) throws Exception {

    // 2d array for the limits
    double [][] limits = new double[data.numAttributes()][2];
    // 2nd dimension: index 0 == lower limit, index 1 == upper limit
    // initialise to no limit
    for(int i = 0; i < data.numAttributes(); i++) {
       limits[i][0] = Double.NEGATIVE_INFINITY;
       limits[i][1] = Double.POSITIVE_INFINITY;
    }

    // use an index instead of creating new Insances objects all the time
    // instanceIndex[0] == array for weights at leaf
    // instanceIndex[1] == array for weights in atbop
    double [][] instanceIndex = new double[2][data.numInstances()];
    // initialize the weight for each instance
    for(int x = 0; x < data.numInstances(); x++) {
        instanceIndex[0][x] = 1;
        instanceIndex[1][x] = 1;  // leaf instances are in atbop
    }

    // first call to graft
    traverseTree(data, instanceIndex, limits, this, 0, -1);
  }


  /**
   * recursive function.
   * if this node is a leaf then calls findGraft, otherwise sorts 
   * the two sets of instances (tracked in iindex array) and calls
   * sortInstances for each of the child nodes (which then calls
   * this method).
   *
   * @param fulldata all instances
   * @param iindex array the tracks the weight of each instance in
   *        the atbop and at the leaf (0.0 if not present)
   * @param limits array specifying current upper/lower limits for numeric atts
   * @param parent the node immediately before the current one
   * @param pL laplace for node, as calculated by parent (in case leaf is empty)
   * @param nodeClass class of node, determined by parent (in case leaf empty)
   */
  private void traverseTree(Instances fulldata, double [][] iindex, 
     double[][] limits, C45PruneableClassifierTreeG parent, 
     double pL, int nodeClass) throws Exception {
    
    if(m_isLeaf) {

       findGraft(fulldata, iindex, limits, 
                 (ClassifierTree)parent, pL, nodeClass);

    } else {

       // traverse each branch
       for(int i = 0; i < localModel().numSubsets(); i++) {

          double [][] newiindex = new double[2][fulldata.numInstances()];
          for(int x = 0; x < 2; x++)
             System.arraycopy(iindex[x], 0, newiindex[x], 0, iindex[x].length);
          sortInstances(fulldata, newiindex, limits, i);
       }
    }
  }

  /**
   * sorts/deletes instances into/from node and atbop according to 
   * the test for subset, then calls traverseTree for subset's node.
   *
   * @param fulldata all instances
   * @param iindex array the tracks the weight of each instance in
   *        the atbop and at the leaf (0.0 if not present)
   * @param limits array specifying current upper/lower limits for numeric atts
   * @param subset the subset for which to sort instances into inode & iatbop
   */
  private void sortInstances(Instances fulldata, double [][] iindex, 
                   double [][] limits, int subset) throws Exception {

    C45Split test = (C45Split)localModel();

    // update the instances index for subset
    double knownCases = 0;
    double thisSubsetCount = 0;
    for(int x = 0; x < iindex[0].length; x++) {
       if(iindex[0][x] == 0 && iindex[1][x] == 0) // skip "discarded" instances
          continue;
       if(!fulldata.instance(x).isMissing(test.attIndex())) {
          knownCases += iindex[0][x];
          if(test.whichSubset(fulldata.instance(x)) != subset) {
             if(iindex[0][x] > 0) {
                // move to atbop, delete from leaf
                iindex[1][x] = iindex[0][x];
                iindex[0][x] = 0;
             } else {
                if(iindex[1][x] > 0) {
                   // instance is now "discarded"
                   iindex[1][x] = 0;
                }
             }
          } else {
             thisSubsetCount += iindex[0][x];
          }
       }
    }

    // work out proportions of weight for missing values for leaf and atbop
    double lprop = (knownCases == 0) ? (1.0 / (double)test.numSubsets()) 
                                : (thisSubsetCount / (double)knownCases);

    // add in the instances that have missing value for attIndex 
    for(int x = 0; x < iindex[0].length; x++) {
       if(iindex[0][x] == 0 && iindex[1][x] == 0)
          continue;     // skip "discarded" instances
       if(fulldata.instance(x).isMissing(test.attIndex())) {
          iindex[1][x] -= (iindex[1][x] - iindex[0][x]) * (1-lprop);
          iindex[0][x] *= lprop;
       }
    }

    int nodeClass = localModel().distribution().maxClass(subset);
    double pL = (localModel().distribution().perClass(nodeClass) + 1.0)
               / (localModel().distribution().total() + 2.0);

    // call traerseTree method for the child node
    son(subset).traverseTree(fulldata, iindex,
          test.minsAndMaxs(fulldata, limits, subset), this, pL, nodeClass);
  }

  /**
   * finds new nodes that improve accuracy and grafts them onto the tree
   *
   * @param fulldata the instances in whole trainset
   * @param iindex records num tests each instance has failed up to this node
   * @param limits the upper/lower limits for numeric attributes
   * @param parent the node immediately before the current one
   * @param pLaplace laplace for leaf, calculated by parent (in case leaf empty)
   * @param pLeafClass class of leaf, determined by parent (in case leaf empty)
   */
  private void findGraft(Instances fulldata, double [][] iindex, 
   double [][] limits, ClassifierTree parent, double pLaplace, 
   int pLeafClass) throws Exception {

    // get the class for this leaf
    int leafClass = (m_isEmpty)
                       ? pLeafClass
                       :  localModel().distribution().maxClass();

    // get the laplace value for this leaf
    double leafLaplace = (m_isEmpty)
                            ? pLaplace
                            : laplaceLeaf(leafClass);

    // sort the instances into those at the leaf, those in atbop, and discarded
    Instances l = new Instances(fulldata, fulldata.numInstances());
    Instances n = new Instances(fulldata, fulldata.numInstances());
    int lcount = 0;
    int acount = 0;
    for(int x = 0; x < fulldata.numInstances(); x++) {
       if(iindex[0][x] <= 0 && iindex[1][x] <= 0)
          continue;
       if(iindex[0][x] != 0) {
          l.add(fulldata.instance(x));
          l.instance(lcount).setWeight(iindex[0][x]);
          // move instance's weight in iindex to same index as in l
          iindex[0][lcount++] = iindex[0][x];
       }
       if(iindex[1][x] > 0) {
          n.add(fulldata.instance(x));
          n.instance(acount).setWeight(iindex[1][x]);
          // move instance's weight in iindex to same index as in n
          iindex[1][acount++] = iindex[1][x];
       }
    }

    boolean graftPossible = false;
    double [] classDist = new double[n.numClasses()];
    for(int x = 0; x < n.numInstances(); x++) {
       if(iindex[1][x] > 0 && !n.instance(x).classIsMissing())
          classDist[(int)n.instance(x).classValue()] += iindex[1][x];
    }

    for(int cVal = 0; cVal < n.numClasses(); cVal++) {
       double theLaplace = (classDist[cVal] + 1.0) / (classDist[cVal] + 2.0);
       if(cVal != leafClass && (theLaplace > leafLaplace) && 
        (biprob(classDist[cVal], classDist[cVal], leafLaplace)
         > m_BiProbCrit)) {
          graftPossible = true;
          break;
       }
    }

    if(!graftPossible) {
       return;
    }

    // 1. Initialize to {} a set of tuples t containing potential tests
    ArrayList t = new ArrayList();

    // go through each attribute
    for(int a = 0; a < n.numAttributes(); a++) {
       if(a == n.classIndex())
          continue;   // skip the class

       // sort instances in atbop by $a
       int [] sorted = sortByAttribute(n, a);

       // 2. For each continuous attribute $a:
       if(n.attribute(a).isNumeric()) {

          // find min and max values for this attribute at the leaf
          boolean prohibited = false;
          double minLeaf = Double.POSITIVE_INFINITY;
          double maxLeaf = Double.NEGATIVE_INFINITY;
          for(int i = 0; i < l.numInstances(); i++) {
             if(l.instance(i).isMissing(a)) {
                if(l.instance(i).classValue() == leafClass) {
                   prohibited = true;
                   break;
                }
             }
             double value = l.instance(i).value(a);
             if(!m_relabel || l.instance(i).classValue() == leafClass) {
                if(value < minLeaf)
                   minLeaf = value;
                if(value > maxLeaf)
                   maxLeaf = value;
             }
          }
          if(prohibited) {
             continue;
	  }

          // (a) find values of
          //    $n: instances in atbop (already have that, actually)
          //    $v: a value for $a that exists for a case in the atbop, where
          //       $v is < the min value for $a for a case at the leaf which
          //       has the class $c, and $v is > the lowerlimit of $a at
          //       the leaf.
          //       (note: error in original paper stated that $v must be
          //       smaller OR EQUAL TO the min value).
          //    $k: $k is a class
          //  that maximize L' = Laplace({$x: $x contained in cases($n)
          //    & value($a,$x) <= $v & value($a,$x) > lowerlim($l,$a)}, $k).
          double minBestClass = Double.NaN;
          double minBestLaplace = leafLaplace;
          double minBestVal = Double.NaN;
          double minBestPos = Double.NaN;
          double minBestTotal = Double.NaN;
          double [][] minBestCounts = null;
          double [][] counts = new double[2][n.numClasses()];
          for(int x = 0; x < n.numInstances(); x++) {
             if(n.instance(sorted[x]).isMissing(a))
                break;   // missing are sorted to end: no more valid vals

             double theval = n.instance(sorted[x]).value(a);
             if(m_Debug)
                System.out.println("\t " + theval);

             if(theval <= limits[a][0]) {
                if(m_Debug)
                   System.out.println("\t  <= lowerlim: continuing...");
                continue;
             }
             // note: error in paper would have this read "theVal > minLeaf)
             if(theval >= minLeaf) {
                if(m_Debug)
                   System.out.println("\t  >= minLeaf; breaking...");
                break;
             }
             counts[0][(int)n.instance(sorted[x]).classValue()]
                += iindex[1][sorted[x]];

             if(x != n.numInstances() - 1) {
                int z = x + 1;
                while(z < n.numInstances()
                 && n.instance(sorted[z]).value(a) == theval) {
                   z++; x++;
                   counts[0][(int)n.instance(sorted[x]).classValue()] 
                    += iindex[1][sorted[x]];
                }
             }

             // work out the best laplace/class (for <= theval)
             double total = Utils.sum(counts[0]);
             for(int c = 0; c < n.numClasses(); c++) {
                double temp = (counts[0][c]+1.0)/(total+2.0);
                if(temp > minBestLaplace) {
                   minBestPos = counts[0][c];
                   minBestTotal = total;
                   minBestLaplace = temp;
                   minBestClass = c;
                   minBestCounts = copyCounts(counts);

                   minBestVal = (x == n.numInstances()-1) 
                      ? theval
                      : ((theval + n.instance(sorted[x+1]).value(a)) / 2.0);
                }
             }
          }

          // (b) add to t tuple <n,a,v,k,L',"<=">
          if(!Double.isNaN(minBestVal)
             && biprob(minBestPos, minBestTotal, leafLaplace) > m_BiProbCrit) {
             GraftSplit gsplit = null;
             try {
                gsplit = new GraftSplit(a, minBestVal, 0,
                                        leafClass, minBestCounts);
             } catch (Exception e) {
                System.err.println("graftsplit error: "+e.getMessage());
                System.exit(1);
             }
             t.add(gsplit);
	  }
          // free space
          minBestCounts = null;

          // (c) find values of
          //    n: instances in atbop (already have that, actually)
          //    $v: a value for $a that exists for a case in the atbop, where
          //       $v is > the max value for $a for a case at the leaf which
          //       has the class $c, and $v is <= the upperlimit of $a at
          //       the leaf.
          //    k: k is a class
          //   that maximize L' = Laplace({x: x contained in cases(n)
          //       & value(a,x) > v & value(a,x) <= upperlim(l,a)}, k).
          double maxBestClass = -1;
          double maxBestLaplace = leafLaplace;
          double maxBestVal = Double.NaN;
          double maxBestPos = Double.NaN;
          double maxBestTotal = Double.NaN;
          double [][] maxBestCounts = null;
          for(int c = 0; c < n.numClasses(); c++) {  // zero the counts
             counts[0][c] = 0;
             counts[1][c] = 0;  // shouldn't need to do this ...
          }

          // check smallest val for a in atbop is < upper limit
          if(n.numInstances() >= 1
           && n.instance(sorted[0]).value(a) < limits[a][1]) {
             for(int x = n.numInstances() - 1; x >= 0; x--) {
                if(n.instance(sorted[x]).isMissing(a))
                   continue;

                double theval = n.instance(sorted[x]).value(a);
                if(m_Debug)
                   System.out.println("\t " + theval);

                if(theval > limits[a][1]) {
                   if(m_Debug)
                      System.out.println("\t  >= upperlim; continuing...");
                   continue;
                }
                if(theval <= maxLeaf) {
                   if(m_Debug)
                      System.out.println("\t  < maxLeaf; breaking...");
                   break;
                }

                // increment counts
                counts[1][(int)n.instance(sorted[x]).classValue()] 
                   += iindex[1][sorted[x]];

                if(x != 0 && !n.instance(sorted[x-1]).isMissing(a)) {
                   int z = x - 1;
                   while(z >= 0 && n.instance(sorted[z]).value(a) == theval) {
                      z--; x--;
                      counts[1][(int)n.instance(sorted[x]).classValue()]
                         += iindex[1][sorted[x]];
                   }
                }

                // work out best laplace for > theval
                double total = Utils.sum(counts[1]);
                for(int c = 0; c < n.numClasses(); c++) {
                   double temp = (counts[1][c]+1.0)/(total+2.0);
                   if(temp > maxBestLaplace ) {
                      maxBestPos = counts[1][c];
                      maxBestTotal = total;
                      maxBestLaplace = temp;
                      maxBestClass = c;
                      maxBestCounts = copyCounts(counts);
                      maxBestVal = (x == 0) 
                        ? theval
                        : ((theval + n.instance(sorted[x-1]).value(a)) / 2.0);
                   }
                }
             }

             // (d) add to t tuple <n,a,v,k,L',">">
             if(!Double.isNaN(maxBestVal)
               && biprob(maxBestPos,maxBestTotal,leafLaplace) > m_BiProbCrit) {
                GraftSplit gsplit = null;
                try {
                   gsplit = new GraftSplit(a, maxBestVal, 1,
                      leafClass, maxBestCounts);
                } catch (Exception e) {
                   System.err.println("graftsplit error:" + e.getMessage());
                   System.exit(1);
                }
                t.add(gsplit);
             }
          }
       } else {    // must be a nominal attribute

          // 3. for each discrete attribute a for which there is no
          //    test at an ancestor of l

          // skip if this attribute has already been used
          if(limits[a][1] == 1) {
             continue;
          }

          boolean [] prohibit = new boolean[l.attribute(a).numValues()];
          for(int aval = 0; aval < n.attribute(a).numValues(); aval++) {
             for(int x = 0; x < l.numInstances(); x++) {
                if((l.instance(x).isMissing(a)
                    || l.instance(x).value(a) == aval) 
                 && (!m_relabel || (l.instance(x).classValue() == leafClass))) {
                   prohibit[aval] = true;
                   break;
                }
             }
          }

          // (a) find values of
          //       $n: instances in atbop (already have that, actually)
          //       $v: $v is a value for $a
          //       $k: $k is a class
          //     that maximize L' = Laplace({$x: $x contained in cases($n)
          //           & value($a,$x) = $v}, $k).
          double bestVal = Double.NaN;
          double bestClass = Double.NaN;
          double bestLaplace = leafLaplace;
          double [][] bestCounts = null;
          double [][] counts = new double[2][n.numClasses()];

          for(int x = 0; x < n.numInstances(); x++) {
             if(n.instance(sorted[x]).isMissing(a))
                continue;

             // zero the counts
             for(int c = 0; c < n.numClasses(); c++)
                counts[0][c] = 0;

             double theval = n.instance(sorted[x]).value(a);
             counts[0][(int)n.instance(sorted[x]).classValue()] 
               += iindex[1][sorted[x]];

             if(x != n.numInstances() - 1) {
                int z = x + 1;
                while(z < n.numInstances() 
                 && n.instance(sorted[z]).value(a) == theval) {
                   z++; x++;
                   counts[0][(int)n.instance(sorted[x]).classValue()]
                      += iindex[1][sorted[x]];
                }
             }

             if(!prohibit[(int)theval]) {
                // work out best laplace for > theval
                double total = Utils.sum(counts[0]);
                bestLaplace = leafLaplace;
                bestClass = Double.NaN;
                for(int c = 0; c < n.numClasses(); c++) {
                   double temp = (counts[0][c]+1.0)/(total+2.0);
                   if(temp > bestLaplace
                    && biprob(counts[0][c],total,leafLaplace) > m_BiProbCrit) {
                      bestLaplace = temp;
                      bestClass = c;
                      bestVal = theval;
                      bestCounts = copyCounts(counts);
                   }
                }
		// add to graft list
                if(!Double.isNaN(bestClass)) {
                   GraftSplit gsplit = null;
                   try {
                      gsplit = new GraftSplit(a, bestVal, 2,
                         leafClass, bestCounts);
                   } catch (Exception e) {
                     System.err.println("graftsplit error: "+e.getMessage());
                     System.exit(1);
                   }
                   t.add(gsplit);
                }
             }
          }
          // (b) add to t tuple <n,a,v,k,L',"=">
          // done this already
       }
    }

    // 4. remove from t all tuples <n,a,v,c,L,x> such that L <=
    //    Laplace(cases(l),c) or prob(x,n,Laplace(cases(l),c) <= 0.05
    //      -- checked this constraint prior to adding a tuple --

    // *** step six done before step five for efficiency ***
    // 6. for each <n,a,v,k,L,x> in t ordered on L from highest to lowest
    // order the tuples from highest to lowest laplace
    // (this actually orders lowest to highest)
    Collections.sort(t);

    // 5. remove from t all tuples <n,a,v,c,L,x> such that there is
    //    no tuple <n',a',v',k',L',x'> such that k' != c & L' < L.
    for(int x = 0; x < t.size(); x++) {
       GraftSplit gs = (GraftSplit)t.get(x);
       if(gs.maxClassForSubsetOfInterest() != leafClass) {
          break; // reached a graft with class != leafClass, so stop deleting
       } else {
          t.remove(x);
          x--;
       }
    }

    // if no potential grafts were found, do nothing and return
    if(t.size() < 1) {
       return;
    }

    // create the distributions for each graft
    for(int x = t.size()-1; x >= 0; x--) {
       GraftSplit gs = (GraftSplit)t.get(x);
       try {
          gs.buildClassifier(l);
          gs.deleteGraftedCases(l); // so they don't go down the other branch
       } catch (Exception e) {
          System.err.println("graftsplit build error: " + e.getMessage());
       }
    }

    // add this stuff to the tree
    ((C45PruneableClassifierTreeG)parent).setDescendents(t, this);
  }

  /**
   * sorts the int array in ascending order by attribute indexed 
   * by a in dataset data.  
   * @param the data the indices represent
   * @param the index of the attribute to sort by
   * @return array of sorted indicies
   */
  private int [] sortByAttribute(Instances data, int a) {

    double [] attList = data.attributeToDoubleArray(a);
    int [] temp = Utils.sort(attList);
    return temp;
  }

  /**
   * deep copy the 2d array of counts
   *
   * @param src the array to copy
   * @return a copy of src
   */
  private double [][] copyCounts(double [][] src) {

    double [][] newArr = new double[src.length][0];
    for(int x = 0; x < src.length; x++) {
       newArr[x] = new double[src[x].length];
       for(int y = 0; y < src[x].length; y++) {
          newArr[x][y] = src[x][y];
       }
    }
    return newArr;
  }
  

  /**
   * Help method for computing class probabilities of
   * a given instance.
   *
   * @throws Exception if something goes wrong
   */
  private double getProbsLaplace(int classIndex, Instance instance, double weight)
       throws Exception {

    double [] weights;
    double prob = 0;
    int treeIndex;
    int i,j;

    if (m_isLeaf) {
       return weight * localModel().classProbLaplace(classIndex, instance, -1);
    } else {
       treeIndex = localModel().whichSubset(instance);

       if (treeIndex == -1) {
          weights = localModel().weights(instance);
          for (i = 0; i < m_sons.length; i++) {
             if (!son(i).m_isEmpty) {
                if (!son(i).m_isLeaf) {
                   prob += son(i).getProbsLaplace(classIndex, instance,
                                                  weights[i] * weight);
                } else {
                   prob += weight * weights[i] *
                     localModel().classProbLaplace(classIndex, instance, i);
                }
             }
          }
          return prob;
       } else {

          if (son(treeIndex).m_isLeaf) {
             return weight * localModel().classProbLaplace(classIndex, instance,
                                                           treeIndex);
          } else {
             return son(treeIndex).getProbsLaplace(classIndex,instance,weight);
          }
       }
    }
  }


  /**
   * Help method for computing class probabilities of
   * a given instance.
   *
   * @throws Exception if something goes wrong
   */
  private double getProbs(int classIndex, Instance instance, double weight)
      throws Exception {

    double [] weights;
    double prob = 0;
    int treeIndex;
    int i,j;

    if (m_isLeaf) {
       return weight * localModel().classProb(classIndex, instance, -1);
    } else {
       treeIndex = localModel().whichSubset(instance);
       if (treeIndex == -1) {
          weights = localModel().weights(instance);
          for (i = 0; i < m_sons.length; i++) {
             if (!son(i).m_isEmpty) {
                prob += son(i).getProbs(classIndex, instance,
                                 weights[i] * weight);
             }
          }
          return prob;
       } else {

          if (son(treeIndex).m_isEmpty) {
             return weight * localModel().classProb(classIndex, instance,
                                                    treeIndex);
          } else {
             return son(treeIndex).getProbs(classIndex, instance, weight);
          }
       }
    }
  }



  /**
   * add the grafted nodes at originalLeaf's position in tree.
   * a recursive function that terminates when t is empty.
   * 
   * @param t the list of nodes to graft
   * @param originalLeaf the leaf that the grafts are replacing
   */
  public void setDescendents(ArrayList t, 
                             C45PruneableClassifierTreeG originalLeaf) {

    Instances headerInfo = new Instances(m_train, 0);

    boolean end = false;
    ClassifierSplitModel splitmod = null;
    C45PruneableClassifierTreeG newNode;
    if(t.size() > 0) {
       splitmod = (ClassifierSplitModel)t.remove(t.size() - 1);
       newNode = new C45PruneableClassifierTreeG(m_toSelectModel, headerInfo,
                           splitmod, m_pruneTheTree, m_CF, m_subtreeRaising,
                           false, m_relabel, m_cleanup);
    } else {
       // get the leaf for one of newNode's children
       NoSplit kLeaf = ((GraftSplit)localModel()).getOtherLeaf();
       newNode = 
             new C45PruneableClassifierTreeG(m_toSelectModel, headerInfo,
                           kLeaf, m_pruneTheTree, m_CF, m_subtreeRaising,
                           true, m_relabel, m_cleanup);
       end = true;
    }

    // behave differently for parent of original leaf, since we don't
    // want to destroy any of its other branches
    if(m_sons != null) {
       for(int x = 0; x < m_sons.length; x++) {
          if(son(x).equals(originalLeaf)) {
             m_sons[x] = newNode;  // replace originalLeaf with newNode
          }
       }
    } else {

       // allocate space for the children
       m_sons = new C45PruneableClassifierTreeG[localModel().numSubsets()];
 
       // get the leaf for one of newNode's children
       NoSplit kLeaf = ((GraftSplit)localModel()).getLeaf();
       C45PruneableClassifierTreeG kNode = 
                 new C45PruneableClassifierTreeG(m_toSelectModel, headerInfo,
                               kLeaf, m_pruneTheTree, m_CF, m_subtreeRaising,
                               true, m_relabel, m_cleanup);
 
       // figure where to put the new node
       if(((GraftSplit)localModel()).subsetOfInterest() == 0) {
          m_sons[0] = kNode;
          m_sons[1] = newNode;
       } else {
          m_sons[0] = newNode;
          m_sons[1] = kNode;
       }
    }
    if(!end)
       ((C45PruneableClassifierTreeG)newNode).setDescendents
                  (t, (C45PruneableClassifierTreeG)originalLeaf);
  }


  /**
   *  class prob with laplace correction (assumes binary class)
   */
  private double laplaceLeaf(double classIndex) {
    double l =  (localModel().distribution().perClass((int)classIndex) + 1.0)
               / (localModel().distribution().total() + 2.0);
    return l;
  }


  /**
   * Significance test
   * @param x 
   * @param n
   * @param r
   * @return returns the probability of obtaining x or MORE out of n
   * if r proportion of n are positive.
   *
   * z for normal estimation of binomial probability of obtaining x 
   * or more out of n, if r proportion of n are positive
   */
  public double biprob(double x, double n, double r) throws Exception {

    return ((((x) - 0.5) - (n) * (r)) / Math.sqrt((n) * (r) * (1.0 - (r))));
  }

  /**
   * Prints tree structure.
   */
  public String toString() {

    try {
       StringBuffer text = new StringBuffer();

       if(m_isLeaf) {
          text.append(": ");
          if(m_localModel instanceof GraftSplit)
             text.append(((GraftSplit)m_localModel).dumpLabelG(0,m_train));
          else
             text.append(m_localModel.dumpLabel(0,m_train));
       } else
          dumpTree(0,text);
       text.append("\n\nNumber of Leaves  : \t"+numLeaves()+"\n");
       text.append("\nSize of the tree : \t"+numNodes()+"\n");

       return text.toString();
    } catch (Exception e) {
       return "Can't print classification tree.";
    }
  }

  /**
   * Help method for printing tree structure.
   *
   * @throws Exception if something goes wrong
   */
  protected void dumpTree(int depth,StringBuffer text) throws Exception {

    int i,j;

    for(i=0;i<m_sons.length;i++) {
       text.append("\n");;
       for(j=0;j<depth;j++)
          text.append("|   ");
       text.append(m_localModel.leftSide(m_train));
       text.append(m_localModel.rightSide(i, m_train));
       if(m_sons[i].m_isLeaf) {
          text.append(": ");
          if(m_localModel instanceof GraftSplit)
             text.append(((GraftSplit)m_localModel).dumpLabelG(i,m_train));
          else
             text.append(m_localModel.dumpLabel(i,m_train));
       } else
          ((C45PruneableClassifierTreeG)m_sons[i]).dumpTree(depth+1,text);
     }
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5532 $");
  }
}
