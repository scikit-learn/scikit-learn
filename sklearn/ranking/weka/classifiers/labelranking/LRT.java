package weka.classifiers.labelranking;

import java.io.File;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.util.Vector;
import weka.classifiers.AbstractClassifier;
import weka.classifiers.Evaluation;
import weka.core.Capabilities;
import weka.core.DenseInstance;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.core.converters.XArffLoader;
import weka.core.labelranking.PreferenceDenseInstance;
import weka.core.labelranking.RankUtilities;



/**
 * This class implements a Label Ranking learner by generating decision trees.
 * For more information, see:
 * "Decision Tree and Instance-Based Learning for Label Ranking"
 * by Eyke Hüllermeier, Weiwei Cheng and Jens Hühn.
 * 
 * @author Jens Hühn, Alexander Balz
 *
 */
@SuppressWarnings({ "deprecation", "rawtypes", "unchecked" })
//This class makes massive use of FastVectors. They are deprecated,
//so many annoying warnings would occur if not suppressed..
public class LRT extends AbstractClassifier implements OptionHandler{
  /**
	 * 
	 */
	private static final long serialVersionUID = -9156584845054553462L;
			
  /**The amount of attributes.*/
  int numVars = 0;
  /**The amount of classes or labels in a ranking.*/
  int numLabels = 0;
  /**This seed is used for a special test scenario described in the paper.*/
  static int seed=0;
  
  LRT tl;
  public LRTree lrTr;
  
  
  int m_minSizeForBranch = 2;
  double m_measureStopCriterion = 0.99;
  
  private static Instances instData;
  
  
  /**
   * Constructor: initially setting the options of 
   * parent class.
   * @throws Exception
   */
  public LRT() throws Exception{
		String[] options = new String[]{};
		
		super.setOptions(options);
		setOptions(options);
  }

  /**
   * Contains the elements of a LRTree like nodes
   * and their evaluation measures.
   *
   */
  public class LRTree{
    Node root; 
    int depth = 0;
    int numInnerNodes = 0;
    
    class Node{
      Node child1, child2, parent;
      double[] ranking = null;
      double measure;
      double splitPoint = Double.NaN;
      int splitVar;


	FastVector trainingData = null;
      private double rightEdges;
      private double wrongEdges;
      
      
      
      /**
       * In here the possible splits are evaluated.
       * Afterwards, the split with best information gain 
       * will be performed, if the stop criterion was not fulfilled.
       * @param data
       * @param p
       */
      private void growTree(FastVector data, Node p){
	parent = p;
	trainingData = data;
	int numInst = data.size();
	
//	if (measure > m_measureStopCriterion)
//	  return;
	
	//Here, all the possible splits are evaluated, trying to find
	//the one with highest information gain.
	double maxInfoGain = 0;
	int maxInfoGainVar = -1;
	double maxInfoGainSplitValue = Double.NaN;
	DataEvaluationMeasure[] maxInfGainDems = null;
	
	for (int var = 0; var < numVars; var++){
	  HashSet<Double> testedSplitValues = new HashSet<Double>();
	  for (int splitIndex = 0; splitIndex < numInst; splitIndex++){
	    double splitValue = ((double[])data.elementAt(splitIndex))[var];
	    if(testedSplitValues.add(splitValue)){
	      DataEvaluationMeasure[] dems = evaluateSplit(data,var,splitValue);
	      if(dems != null){
//		double infoGain = dems[0].size * Utils.log2(dems[0].bestResult) + dems[1].size  * Utils.log2(dems[1].bestResult) - Utils.log2(this.measure);
		double infoGain = dems[0].size * Utils.log2(dems[0].bestResult) + dems[1].size  * Utils.log2(dems[1].bestResult) - Utils.log2(this.measure);
		if (Utils.gr(infoGain,maxInfoGain)){
		  maxInfoGain = infoGain;
		  maxInfoGainVar = var;
		  maxInfoGainSplitValue = ((double[])data.elementAt(splitIndex))[var];
		  maxInfGainDems = dems;
		}
	      }
	    }
	  }
	  
	}
	
	//stop-condition
	if (Utils.smOrEq(maxInfoGain,0))
	  return;
	
	// test whether the child nodes return the same ranking. If yes, remove both and make this node to a leaf
//	boolean splitIsStupid = true;
//	for (int i = 0; i < numLabels; i++){
//	  if (maxInfGainDems[0].bestResultRanking[i] != maxInfGainDems[1].bestResultRanking[i]){
//	    splitIsStupid = false;
//	    break;
//	  }
//	}
//	if (splitIsStupid)
//	  System.out.println("STUPID SPLIT");
	
		
	this.splitPoint = maxInfoGainSplitValue;
	this.splitVar = maxInfoGainVar;


	// split once again
	//Performing split with optimal information gain.
	FastVector dataCopy1 = new FastVector();
	FastVector dataCopy2 = new FastVector();
	for (int i = 0; i < data.size(); i++){
		
//	if(((double[])data.elementAt(i))[maxInfoGainVar]== maxInfoGainSplitValue)	
		if(!instData.attribute(maxInfoGainVar).isNumeric()){
			if(((double[])data.elementAt(i))[maxInfoGainVar]== maxInfoGainSplitValue)
				dataCopy1.addElement(((double[])data.elementAt(i)));
			else
				dataCopy2.addElement(((double[])data.elementAt(i)));
		}
		else{
			if(((double[])data.elementAt(i))[maxInfoGainVar]<= maxInfoGainSplitValue)
				dataCopy1.addElement(((double[])data.elementAt(i)));
			else
				dataCopy2.addElement(((double[])data.elementAt(i)));
		}
	}
	
	if (dataCopy1.size() < m_minSizeForBranch || dataCopy2.size() < m_minSizeForBranch ){
	  splitPoint = Double.NaN;
	  return;
	}
	

	  child1 = new Node();
	  child2 = new Node();
	  child1.ranking = maxInfGainDems[0].bestResultRanking;
	  child2.ranking = maxInfGainDems[1].bestResultRanking;
	  child1.measure = maxInfGainDems[0].bestResult;
	  child2.measure = maxInfGainDems[1].bestResult;
	  child1.rightEdges = maxInfGainDems[0].rightEdges;
	  child1.wrongEdges = maxInfGainDems[0].wrongEdges;
	  child2.rightEdges = maxInfGainDems[1].rightEdges;
	  child2.wrongEdges = maxInfGainDems[1].wrongEdges;
	  child1.growTree(dataCopy1, this);
	  child2.growTree(dataCopy2, this);

      }
      
      /**
       * Pruning the tree, if
       * dem.bestResult > (dataCopy1.size()/data.size()) * result1 + (dataCopy2.size()/data.size()) * result2)
       * 
       * @param data
       * @return
       *
      private double pruneTree(FastVector data){
	DataEvaluationMeasure dem = new DataEvaluationMeasure(data);
	double result = dem.bestResult;
	
	if(!isLeaf()){
	  FastVector dataCopy1 = new FastVector();
	  FastVector dataCopy2 = new FastVector();
	  for (int i = 0; i < data.size(); i++){
	    if(((double[])data.elementAt(i))[splitVar]<= splitPoint)
	      dataCopy1.addElement(((double[])data.elementAt(i)));
	    else
	      dataCopy2.addElement(((double[])data.elementAt(i)));
	  }
	  double result1 = child1.pruneTree(dataCopy1);
	  double result2 = child2.pruneTree(dataCopy2);
	  
	  
	  if (Utils.gr(dem.bestResult,(dataCopy1.size()/(double)data.size()) * result1 + (dataCopy2.size()/(double)data.size()) * result2)){
	    //prune
	    child1 = null;
	    child2 = null;
	    splitPoint = Double.NaN;
	  }else{
	    result = (dataCopy1.size()/(double)data.size()) * result1 + (dataCopy2.size()/(double)data.size()) * result2;
	  }
	  
	}
	
	return result;
	
      }*/
      
      
      
      
      /**
       * Find the goodness-of-split measure i.e. the weighted average
       * of within-leaf variances.
       * @param dataF
       * @param splitIndex
       * @return
       */
      private DataEvaluationMeasure[] evaluateSplit(FastVector dataIn, int var, double splitValue) {
	FastVector data1 = new FastVector();
	FastVector data2 = new FastVector();
	for (int i = 0; i < dataIn.size(); i++){
		

		if(!instData.attribute(var).isNumeric()){
			if(((double[])dataIn.elementAt(i))[var] == splitValue)
				 data1.addElement(((double[])dataIn.elementAt(i)));
			  else
			    data2.addElement(((double[])dataIn.elementAt(i)));
		}
		else{
	  if (((double[])dataIn.elementAt(i))[var] <= splitValue)
	    data1.addElement(((double[])dataIn.elementAt(i)));
	  else
	    data2.addElement(((double[])dataIn.elementAt(i)));
		}
	}
	
	if (data1.size() < m_minSizeForBranch || data2.size() < m_minSizeForBranch)
	  return null;
		
	DataEvaluationMeasure dem1 = new DataEvaluationMeasure(data1);
	
	
//	if((rE+wE) > 0 && (rE+wE) !=  (dem1.rightEdges+dem1.wrongEdges)  + (dem2.rightEdges+dem2.wrongEdges))
//	System.out.println((rE+wE) + "\t" + (dem1.rightEdges+dem1.wrongEdges) +"\t" + (dem2.rightEdges+dem2.wrongEdges));
		
	// If first split is so bad that it is not possible to beat the parents performance then break immediately
	double theoreticalMaxInformationGain = 
	  (dem1.rightEdges+dem1.wrongEdges)/(rightEdges+wrongEdges) * Utils.log2(dem1.bestResult) - /*rE/(rE+wE) */ Utils.log2(measure);
	if (rightEdges+wrongEdges > 0 && Utils.smOrEq(theoreticalMaxInformationGain, 0))
	  return null;
	
	
	DataEvaluationMeasure dem2 = new DataEvaluationMeasure(data2);
	dem1.size = (dem1.rightEdges+dem1.wrongEdges) / (dem1.rightEdges+dem1.wrongEdges + dem2.rightEdges+dem2.wrongEdges);
	dem2.size = (dem2.rightEdges+dem2.wrongEdges) / (dem1.rightEdges+dem1.wrongEdges + dem2.rightEdges+dem2.wrongEdges);
		
	return new DataEvaluationMeasure[]{dem1,dem2};
      }

      /**
       * 
       * @author jens
       *
       */
      class DataEvaluationMeasure{
	public double bestResult = 0;
	public double[] bestResultRanking = null; 
	public double size = 0;
	private double rightEdges;
	private double wrongEdges;

	public DataEvaluationMeasure(FastVector data) {
	  boolean[] labelStillAvailable = new boolean[numLabels]; 
	  
	  double[][] graph = new double[numLabels][numLabels];
	  for (int inst = 0; inst < data.size(); inst++){
	    for (int i = numVars; i < numVars+numLabels; i++){
	      if (((double[])data.elementAt(inst))[i] > 0){
		if (!labelStillAvailable[i-numVars]) labelStillAvailable[i-numVars] = true;
		for (int j = i+1; j < numVars+numLabels; j++){
		  // Do not count missing labels (pos 0)
		 
		  if (((double[])data.elementAt(inst))[j] > 0){
		    if(((double[])data.elementAt(inst))[i] < ((double[])data.elementAt(inst))[j]){
		      graph[i-numVars][j-numVars]++;
		    }else{
		      graph[j-numVars][i-numVars]++;
		    }
//		    if (!labelStillAvailable[j-numVars]) labelStillAvailable[j-numVars] = true; 
		  }
		}
	      }
	    }
	  }
	  	  
	  double[][] graphBackup = new double[graph.length][graph.length];
	  for (int i = 0; i < graph.length; i++){
	    for (int j = 0; j < graph.length; j++){
	      graphBackup[i][j] = graph[i][j];
	    }
	  }
	  
	  FastVector s1 = new FastVector();
	  FastVector s2 = new FastVector();
	  
	  int countOfNotAvailableLabels = 0 ;
	  
	  HashSet<Double> lSet = new HashSet<Double>();
	  double[] ranking = new double[numLabels];
	  for (int i = 0; i < numLabels; i++){
	    //only add label to available set of labels if it is still here
	    if(labelStillAvailable[i])
	      lSet.add((double)i);
	    else{
	      countOfNotAvailableLabels++;
	      ranking[i] = -1;
	    }
	  }
	  
	  
	  while(lSet.size() > 0){
	    //remove all sinks
	    boolean tryFindingSink = true;
	    while(tryFindingSink){
	      tryFindingSink = false;
	      for (int i = 0; i < numLabels; i++){
		if (Utils.sum(graph[i]) == 0 && lSet.contains(i)){
		  lSet.remove(i);
		  FastVector temp = new FastVector();
		  temp.addElement(i);
		  temp.appendElements(s2);
		  s2 = temp;
		  for (int j = 0; j < numLabels; j++){
		    graph[j][i] = 0; 
		  }
		  tryFindingSink = true;
		  break;
		}
	      }
	    }
	    
	    // remove all sources
	    boolean tryFindingSource = true;
	    while(tryFindingSource){
	      tryFindingSource = false;
	      for (int i = 0; i < numLabels; i++){
		double inBound = 0;
		for (int j = 0; j < numLabels; j++){
		  inBound += graph[j][i];
		}		  
		if (inBound == 0 && lSet.contains(i)){
		  lSet.remove(i);
		  s1.addElement(i);
		  for (int j = 0; j < numLabels; j++){
		    graph[i][j] = 0; 
		  }
		  tryFindingSource = true;
		  break;
		}
	      }
	    }
	    
	    if (lSet.size() > 0){
	      double[] deltas = new double[numLabels];
	      for (int i = 0; i < numLabels; i++){
		for (int j = 0; j < numLabels; j++){
		  deltas[i] += graph[i][j] - graph[j][i];
		}
	      }
	      double maxDeltaU;
	      boolean onlyZeros = true;
	      for (int i = 0; i < deltas.length; i++){
		if (deltas[i] != 0.0){
		  onlyZeros = false;
		  break;
		}
	      }
	      if (onlyZeros){
	    	  //int seed=m_Seed;
	   // 	  int sd = (int)seed;

		Random rand = new Random(seed);
		maxDeltaU = (Double)(lSet.toArray()[rand.nextInt(lSet.size())]);
		if (parent != null){
		  double tRank = parent.getFirstRankedLabel(lSet);
		  if (tRank > -1){
		    maxDeltaU = tRank;
		  }
		}
	      }else{
		maxDeltaU = Utils.maxIndex(deltas);
	      }
			      
	      s1.addElement(maxDeltaU);
	      lSet.remove(maxDeltaU);


	      for (int i = 0; i < numLabels; i++){
		for (int j = 0; j < numLabels; j++){
		  if (i == maxDeltaU || j == maxDeltaU)
		    graph[i][j] = 0;
		}
	      }  
	    }
	  }
	  
	  s1.appendElements(s2);
	  
	  for (int r = 0; r < numLabels-countOfNotAvailableLabels; r++){
	    ranking[((Double)s1.elementAt(r)).intValue()] = (double)r;
	  }
	  
	  graph = graphBackup;

	  
	  double result = -1;
	  double rightDir = 0;
	  double wrongDir = 0;
	  
	  for (int i = 0; i < ranking.length; i++){
	      for (int j = i+1; j < ranking.length; j++){
		if (labelStillAvailable[i] && labelStillAvailable[j]){
		  if(ranking[i] < ranking[j]){
		    rightDir += graph[i][j];
		    wrongDir += graph[j][i];
		  }else{
		    rightDir += graph[j][i];
		    wrongDir += graph[i][j];
		  }
		}
	      }
	    }

	    result = rightDir/(rightDir+wrongDir);
	    if (result > bestResult || (result == bestResult && bestResultRanking != null && ranking.length > bestResultRanking.length)){
	      bestResult = result;
	      bestResultRanking = ranking;
	      rightEdges = rightDir;
	      wrongEdges = wrongDir;
	    }
	}
      }
      
      /**
       * 
       * @return
       */
      private boolean isLeaf(){
	if(child1 == null && child2 == null){
	  return true;
	}else{
	  return false;
	}
      }
      
      
      public double getFirstRankedLabel(HashSet<Double> lset) {
	// if not all desired labels are available, directly go up
	for (int i = 0; i < lset.size(); i++){
	  if (ranking[((Double)lset.toArray()[i]).intValue()] == -1.0){
	    return parent.getFirstRankedLabel(lset);
	  }
	}
	      
	double firstRankedLabelPos = Double.MAX_VALUE;
	double firstRankedLabel = -1;
	for (int i = 0; i < lset.size(); i++){
	  if (ranking[((Double)lset.toArray()[i]).intValue()] >= 0 && ranking[((Double)lset.toArray()[i]).intValue()] < firstRankedLabelPos){
	    firstRankedLabelPos = ranking[((Double)lset.toArray()[i]).intValue()];
	    firstRankedLabel = (Double)lset.toArray()[i]; 
	  }

	}
	if (firstRankedLabelPos == Double.MAX_VALUE){
	  if (parent != null)
	    return parent.getFirstRankedLabel(lset);
	  else return -1;
	}else
	  return firstRankedLabel;
      }

      /**
       * Contains the visualization returned by dumpTree
       * @return label ranking tree visualization.
       */
      public String toString(){
	StringBuffer text = new StringBuffer();
	try {
	  dumpTree(0,text);
        } catch (Exception e) {
	  e.printStackTrace();
        }
	return text.toString();
      }

      /**
       * Returning a visualization of the tree.
       * @param depth
       * @param text
       * @throws Exception
       */
      private void dumpTree(int depth, StringBuffer text) 
      throws Exception {

	int j;

	
	text.append("\n");;
	for (j=0;j<depth;j++)
	  text.append("|   ");
	if (child1 == null){
	  text.append(dumpLabel(0));
	}else{
	  text.append(leftSide());
	  text.append(rightSide(0));
	  if (child1.isLeaf()) {
	    text.append(": ");
	    text.append(dumpLabel(0));
	  }else
	    child1.dumpTree(depth+1,text);
	}
	
	
	text.append("\n");;
	for (j=0;j<depth;j++)
	  text.append("|   ");
	if (child2 == null) {
	  text.append(dumpLabel(1));
	}else{
	  text.append(leftSide());
	  text.append(rightSide(1));
	  if (child2.isLeaf()) {
	    text.append(": ");
	    text.append(dumpLabel(1));
	  }else
	    child2.dumpTree(depth+1,text);
	}
	
      }
      
      /**
       * Returns information for each node, like accuracy
       * and ranking.
       * @param index
       * @return
       * @throws Exception
       */
      public final String dumpLabel(int index) throws Exception {
	StringBuffer text;

	Node nodeToOutput = this;
	if (child1 != null || child2 != null){
	  if (index == 0){
	    nodeToOutput = child1;
	  }else{
	    nodeToOutput = child2;
	  }
	}
	
	
	text = new StringBuffer();
	boolean labelOutputtedYet = false;
		for (int i = 0; i < nodeToOutput.ranking.length; i++){
		  if (nodeToOutput.ranking[i] > -1){
			  if (i>0 && nodeToOutput.ranking[i-1] == -1 && labelOutputtedYet)
				  text.append(">");
			  text.append(((int)nodeToOutput.ranking[i])+1);
			  labelOutputtedYet = true;
			  if (i<nodeToOutput.ranking.length-1 && nodeToOutput.ranking[i+1] > -1)
				  text.append(">");
			 }
		}
	
	
	


	  int correct = 0;
	  int incorrect = 0;

	  for (int i = 0; i < nodeToOutput.trainingData.size(); i++){
	    for (int j = numVars; j < numVars+numLabels; j++){
	      if(((double[])nodeToOutput.trainingData.elementAt(i))[j] == nodeToOutput.ranking[j-numVars]){
		incorrect++;
		break;
	      }
	    }
	    correct++;
	  }
	  text.append(" (Acc: " + correct +"/" + incorrect + ")");

	return text.toString();
      }
      
      
      /**
       * Prints left side of condition..
       *
       * @param data training set.
       */
      public final String leftSide() {
        return "Att " + (splitVar+1);
      }

      /**
       * Prints the condition satisfied by instances in a subset.
       *
       * @param index of subset 
       * @param data training set.
       */
      public final String rightSide(int index) {

	StringBuffer text;

	text = new StringBuffer();

	if(!instData.attribute(splitVar).isNumeric()){
			if (index == 0){
				text.append("  { "+
						instData.attribute(splitVar).value((int)splitPoint)+"}");
			}
			else {			
				text.append("  { ");
				for(int i=0; i<instData.attribute(splitVar).numValues(); i++){
					if(i!=(int)splitPoint){
						text.append(instData.attribute(splitVar).value(i)+",");
					}
				}
				text.append(" }");
			}
	}
	else{
		if (index == 0)
			  text.append(" <= "+
			      Utils.doubleToString(splitPoint,6));
			else
			  text.append(" > "+
			      Utils.doubleToString(splitPoint,6));
	}
		
		text.append(" maxInfoMeasure: " + Utils.doubleToString(measure,4) + " ");
		return text.toString();
	      }
	
	

      public void createStats(int d) {
	numInnerNodes++;	
	if (depth<d){
	  depth = d;
	}
	
	//prevents crash if informationGain has been 0.
	if(child1==null && child2==null){}
	else{
		if (!child1.isLeaf())
			child1.createStats(d+1);

		if (!child2.isLeaf())
			child2.createStats(d+1);
      }
    }


      /**
       * Returns a complete ranking for features
       * stored in r.
       * @param r
       * @return
       * @throws Exception
       */
      private double[] completeRanking(double[] r) throws Exception {
	double[] comparisonRanking = ranking;
	if (numMissingLabels() > 0)
	  r = parent.completeRanking(r);
	
	// Test whether ranking is complete
	HashSet<Double> missingLabels = new HashSet<Double>();
	for (int i = 0; i < numLabels; i++){
	  if (r[i] == -1){
	    missingLabels.add((double)i);
	  }
	}
	while (missingLabels.size() > 0){
	  double bestMeasure = -1;
	  double[] bestResultRanking = r;
	  int posToStartInserting = 0;
	  double firstOfMissingLabels  = -1;
	  if (missingLabels.size() == 1){
	    firstOfMissingLabels = (Double)missingLabels.toArray()[0]; 
	  }else{
	    firstOfMissingLabels = getFirstRankedLabel(missingLabels);
	  }
	  double[] firstParentRankingContainingFOML = comparisonRanking;
	  missingLabels.remove(firstOfMissingLabels);
	  // try to insert in a optimal way
	  for (int j = posToStartInserting; j < numLabels-missingLabels.size(); j++){
	    double[] potentialRanking = new double[numLabels];
	    for (int k = 0; k < numLabels; k++){
	      if (k == (int)firstOfMissingLabels){
		potentialRanking[k] = j;
	      }else{
		if (r[k] > -1)
		  if (r[k] < j)
		    potentialRanking[k] = r[k];
		  else 
		    potentialRanking[k] = r[k]+1;
		else potentialRanking[k] = -1;
	      }
	    }

	    double[] temp1R = new double[numLabels];
	    double[] temp2R = new double[numLabels];
	    for (int k = 0; k < numLabels; k++){
	      if (potentialRanking[k] == -1){
		temp1R[(int)(k)] = 0;
		temp2R[(int)(k)] = 0;
	      }else{
		temp1R[(int)(k)] = potentialRanking[k];
		temp2R[(int)(k)] = firstParentRankingContainingFOML[k];
	      }
	    }

	    double measure = TreeLabelRankingEvaluation.kendall(temp1R, temp2R);
	    if (measure > bestMeasure){
	      bestMeasure = measure;
	      bestResultRanking = potentialRanking;
	      posToStartInserting = j+1;
	    }
	  }
	  r = bestResultRanking;


	}


	return r;
      }
    
      /**
       * Counting missing labels inside of a ranking.
       * @return
       */
    private double numMissingLabels(){
      double missingCount = 0;
      for (int i = 0; i < ranking.length;i++){
	if (ranking[i] == -1){
	  missingCount++;
	}
      }
      return missingCount;
    }

    private boolean hasMissingLabels() {
      for (int i = 0; i < ranking.length;i++){
	if (ranking[i] == -1){
	  return true;
	}
      }
      return false;
    }

      
    }
    
    /**
     * Initialization and generation of our tree.
     * @param data
     */
    private void buildModel(FastVector data){
      root = new Node();
      Node.DataEvaluationMeasure dem = root.new DataEvaluationMeasure(data);
      root.ranking = dem.bestResultRanking;
      root.measure = dem.bestResult;
           
      root.growTree(data,null);
      root.createStats(1);
      toString();
    }
    
    /**
     * Returns the predicted ranking for an input
     * instance. Uses the completeRanking method
     * in order to do so.
     * @param inst
     * @return
     * @throws Exception
     */
    public double[] getLabelRanking(double[] inst) throws Exception {
      Node iterNode = root;      
      
      if(!instData.attribute(iterNode.splitVar).isNumeric()){
    	  while(!iterNode.isLeaf()){
    		  if (inst[iterNode.splitVar] == iterNode.splitPoint){
    			  iterNode = iterNode.child1;
    		  }else{
    			  iterNode = iterNode.child2;
    		  }
    	  }
      }
      else{     
    	  while(!iterNode.isLeaf()){
    		  if (inst[iterNode.splitVar] <= iterNode.splitPoint){
    			  iterNode = iterNode.child1;
    		  }else{
    			  iterNode = iterNode.child2;
    		  }
    	  }
      }
      
      double[] resultRanking = iterNode.ranking; 
      if(iterNode.hasMissingLabels()){
	resultRanking = iterNode.completeRanking(resultRanking);
      }
      
      double[] rankingFormatted = new double[numLabels];
      for (int i = 0; i < resultRanking.length;i++){
	rankingFormatted[i] = resultRanking[i]+1;
      }
      
      return rankingFormatted;
    }

 
    /**
     * This toString method returns general information
     * of the tree and its visualization.
     * @return information of the tree and its visualization.
     */
    @Override
    public String toString() {
      StringBuffer text = new StringBuffer();
      text.append("Num Inner Nodes: " + numInnerNodes + "\n");
      text.append("Depth: " + depth + "\n\n");
      text.append(root.toString());
           
      return text.toString();
    }       
  }

  /**
   * A special main method executing evaluations described in the paper
   * of the TreeLabelRanking algorithm.
   */
  public static void mainDeleteLabels(String[] args) throws Exception{
	  long startTime = System.currentTimeMillis();

	    double avgtau = 0;
	    double testCounter = 0;
	    
	    	    
	    double[] missingRatesToTest = new double[]{0,0.3,0.6};
	    	    
	    // Setting for paper
	    final int numtrials = 5;
	    final int numfolders = 10;
	    int numDS = 16;
	    	    
	    String dataset = new String();
	    String subfolder = new String();
	    String option = new String();
	    
	    double[][][] totalResults = new double[numDS][missingRatesToTest.length][4];
	    FastVector dsNames = new FastVector(); 
	    
	    for (int mr = 0; mr < missingRatesToTest.length;mr++){
	      double missingRate = missingRatesToTest[mr]; // 1 = all labels missing, 0 = no label missing
	      System.out.println(missingRate);
	    
	    String folder="/home/alex/Desktop/RankingDatasets/";
	    
	      for (int ds = 0; ds < numDS; ds++){
		switch(ds){
		
		case 0:
		  dataset = "analcatdata-authorship_dense.xarff";//args[0];
		  subfolder = folder;
		  break;		  
		case 1:
		  dataset = "bodyfat_dense.xarff";//args[0];
		  subfolder = folder;
		  break;
		case 2:
		  dataset = "calhousing_dense.xarff";//args[0];
		  subfolder = folder;
		  break;
		case 3:
		  dataset = "cpu-small_dense.xarff";//args[0];
		  subfolder = folder;
		  break;
		case 4:
		  dataset = "elevators_dense.xarff";//args[0];
		  subfolder = folder;
		  break;
		case 5:
		  dataset = "fried_dense.xarff";//args[0];
		  subfolder = folder;
		  break;
		case 6:
		  dataset = "glass_dense.xarff";//args[0];
		  subfolder = folder;
		  break;
		case 7:
		  dataset = "housing_dense.xarff";//args[0];
		  subfolder = folder;
		  break;
		case 8:
		  dataset = "iris_dense.xarff";//args[0];
		  subfolder = folder;
		  break;
		case 9:
		  dataset = "pendigits_dense.xarff";//args[0];
		  subfolder = folder;
		  break;
		case 10:
		  dataset = "segment_dense.xarff";//args[0];
		  subfolder = folder;
		  break;
		case 11:
		  dataset = "stock_dense.xarff";//args[0];
		  subfolder = folder;
		  break;
		case 12:
		  dataset = "vehicle_dense.xarff";//args[0];
		  subfolder = folder;
		  break;
		case 13:
		  dataset = "vowel_dense.xarff";//args[0];
		  subfolder = folder;
		  break;
		case 14:
		  dataset = "wine_dense.xarff";//args[0];
		  subfolder = folder;
		  break;
		case 15:
		  dataset = "wisconsin_dense.xarff";//args[0];
		  subfolder = folder;
		  break;
		  default:
			  continue;
		}

		if (mr == 0){
		  dsNames.addElement(dataset.subSequence(0, dataset.indexOf("_")));
		}
		
		if (ds!=15)
		  continue;
	//
//		if (subset){
//		  if (ds ==0 || ds == 3 || ds == 4 || ds == 5 || ds == 6)
//		    continue;
////		  if (ds != 1)
////		  continue;
//		}


		double tau_sum = 0;
		double rho_sum = 0;
		double[] taus_total = new double[numfolders*numtrials];
		double[] rhos_total = new double[numfolders*numtrials];

		for(int cvindex = 0; cvindex<numtrials; cvindex++){

		  double[] taus = new double[numfolders];
		  double[] rhos = new double[numfolders];
		  String[] path = new String[numfolders];

		  //Global parameters setting starts
		  //Be careful of the meaning of labels

		  for(int i=0; i<path.length; i++){
		    path[i] = subfolder + dataset;
		  }


		  option = "soft";

		  //Global parameters setting overs

		  String testfolder = new String();
		  for(int f=numfolders-1; f>=0; f--){
		    testfolder = path[f];
		    path[f] = path[numfolders-1];
		    path[numfolders-1] = testfolder;



		    LRT tlr = new LRT();
		    FastVector trainingData = new FastVector();;
		    FastVector fv = new FastVector();
		    double[] feat = new double[5];
		    for(int i=0; i<numfolders-1; i++){
		    	fv = new FastVector();
		    	XArffLoader xarff = new XArffLoader();
	    		File file = new File(path[i]);
	    		xarff.setFile(file);
	    		Instances inst = xarff.getDataSet();
	    		inst.setClassIndex(inst.numAttributes()-1);
	    		tlr.numLabels = inst.getLabels().size();
    			tlr.numVars = inst.numAttributes()-1;
	    		
	    		
	    		for(int j=0; j<inst.numInstances(); j++){
	    			
	    			feat = new double[inst.numAttributes()-1+inst.getLabels().size()];
	    			double[] labels = new double[inst.getLabels().size()];
	    			int ft=0;
	    			PreferenceDenseInstance pdi = (PreferenceDenseInstance)inst.get(j);
	    			for(int x=0; ; x++){
	    				
	    				if(pdi.getHashMap().get(x)!=null){	    			
	    					RankUtilities.noPartialRanking(pdi.getHashMap().get(x));
	    					labels = RankUtilities.triangleLabels;
	    					break;
	    				}
	    				
	    			}
	    			
	    			for(int k=0; k<feat.length; k++){
	    				if(k<inst.numAttributes()-1){
	    				  feat[k] = inst.get(j).value(k);
	    				}
	    				else{
	    					feat[k] = labels[ft]+1;
	    					ft++;
	    				}
	    			}
	    			
	  		      fv.add(feat);	    			
	    		}
	    		
	    	  trainingData.appendElements(fv);
		      
		    }


		    trainingData = labeldelete(trainingData,missingRate, tlr.numVars, seed + trainingData.hashCode());


		    LRTree lrTree = tlr.new LRTree();
		    lrTree.buildModel(trainingData);

		    //
		    // Prediction 
		    //   
		    FastVector testingData = new FastVector();;
		    for(int i=0; i<numfolders-1; i++){
		        //Putting predictions in Instances.
		    	fv = new FastVector();
		    	XArffLoader xarff = new XArffLoader();
	    		File file = new File(path[numfolders-1]);
	    		xarff.setFile(file);
	    		Instances inst = xarff.getDataSet();
	    		inst.setClassIndex(inst.numAttributes()-1);
	    		
	    		
	    		for(int j=0; j<inst.numInstances(); j++){
	    			
	    			feat = new double[inst.numAttributes()-1+inst.getLabels().size()];
	    			double[] labels = new double[inst.getLabels().size()];
	    			int ft=0;
	    			PreferenceDenseInstance pdi = (PreferenceDenseInstance)inst.get(j);
	    			for(int x=0; ; x++){
	    				
	    				if(pdi.getHashMap().get(x)!=null){	    			
	    					RankUtilities.noPartialRanking(pdi.getHashMap().get(x));
	    					labels = RankUtilities.triangleLabels;
	    					break;
	    				}
	    				
	    			}
	    			
	    			for(int k=0; k<feat.length; k++){
	    				if(k<inst.numAttributes()-1){
	    				  feat[k] = inst.get(j).value(k);
	    				}
	    				else{
	    					feat[k] = labels[ft]+1;
	    					ft++;
	    				}
	    			}
	    			
	  		      fv.add(feat);
	    			
	    		}
	    		
	    	  testingData.appendElements(fv);
		    	
		    	
		    	//testingData.appendElements(tlr.readTxTDataFile(path[numfolders-1]));
		    }




		    ArrayList<double[]> features = new ArrayList<double[]>();
		    ArrayList<double[]> labels = new ArrayList<double[]>();

		    for(int i=0; i<testingData.size(); i++){
		    	
		    //tlr.numVars = testingData.get(i).
		      double[] feature = new double[tlr.numVars];
		      double[] label = new double[tlr.numLabels];
		      for(int j=0; j<tlr.numVars; j++){
			feature[j] = ((double[])testingData.elementAt(i))[j];
		      }   
		      features.add(feature);

		      for(int j=tlr.numVars; j<tlr.numVars+tlr.numLabels; j++){                    
			label[j-tlr.numVars] = ((double[])testingData.elementAt(i))[j];
		      }
		      // from order to rank if needed (depends on the dataset. if not simply add array label)
		      //  Start here
//		      double[] label_trans = new double[label.length];
//		      for(int k=0; k<label.length; k++){
//		      label_trans[(int)label[k]-1] = k+1;
//		      }
		      // End here
		      labels.add(label);
		    }

		    //get predicted labels.
		    ArrayList<double[]> labelrankings = new ArrayList<double[]>(); 
		    for(int i=0; i<testingData.size(); i++){
		      labelrankings.add(lrTree.getLabelRanking(features.get(i)));
		    }

		    //compute the correlation coefficients
		    double tau = 0;
		    double rho = 0;
		    ArrayList<Double> plabels;
		    ArrayList<Double> tlabels;
		    for(int i=0; i<labels.size(); i++){
		      plabels = new ArrayList<Double>();
		      tlabels = new ArrayList<Double>();
		      for(int j=0; j<labels.get(i).length; j++){
			plabels.add(labelrankings.get(i)[j]-1);
			tlabels.add(labels.get(i)[j]-1);
		      }
		      tau = tau + TreeLabelRankingEvaluation.kendall(plabels, tlabels);
		      rho = rho + TreeLabelRankingEvaluation.spearman(plabels, tlabels);
		    }

		    taus[f] = tau/labels.size();
		    rhos[f] = rho/labels.size();

		  }

		  for(int i=0; i<numfolders; i++){
		    tau_sum = tau_sum + taus[i];
		    taus_total[numfolders*cvindex+i] = taus[i];
		    rho_sum = rho_sum + rhos[i];
		    rhos_total[numfolders*cvindex+i] = rhos[i];
		  }

		}

		System.out.println(numtrials + " times " + numfolders + " folder cv.\t" + dataset + ".\t" + option + "-voting." +
		    "\ntau: " +  tau_sum/(numfolders*numtrials) +
		    "\tsigma: " +  TreeLabelRankingEvaluation.standarddeviation(taus_total) +
		    "\nrho: " + rho_sum/(numfolders*numtrials) +
		    "\tsigma: " +  TreeLabelRankingEvaluation.standarddeviation(rhos_total));
		
		totalResults[ds][mr][0] = tau_sum/(numfolders*numtrials);
		totalResults[ds][mr][1] = TreeLabelRankingEvaluation.standarddeviation(taus_total);
		totalResults[ds][mr][2] = rho_sum/(numfolders*numtrials);
		totalResults[ds][mr][3] = TreeLabelRankingEvaluation.standarddeviation(rhos_total);
		
		avgtau += tau_sum/(numfolders*numtrials);
		testCounter++;
	      }
	      System.out.println("AVG TAU: " + avgtau/testCounter);
	    }
	    
	    for (int h = 0; h < 4; h++){
	      switch(h){
	      case 0:
		System.out.print("tau\r");
		break;
	      case 1:
		System.out.print("tau_sigma\r");
		break;
	      case 2:
		System.out.print("rho\r");
		break;
	      case 3:
		System.out.print("tau_sigma\r");
		break;
	      }
	      
	      System.out.print("data set\t");
	      for (int i = 0; i < missingRatesToTest.length; i++){
		System.out.print(missingRatesToTest[i] + "\t");
	      }
	      System.out.print("\r");

	      for (int i = 0; i < totalResults.length; i++){
		System.out.print(dsNames.elementAt(i) + "\t");
		for (int j = 0; j < missingRatesToTest.length; j++){
		  System.out.print(totalResults[i][j][h] + "\t");
		}
		System.out.print("\r");
	      }
	    }

	    System.out.println("Time: " + (System.currentTimeMillis() - startTime));
  }
  
  /**
   * Automatic evaluations of multiple data sets, calculating Kendall's
   * tau and Spearman rank correlation.
   * @param args
   * @throws Exception
   */
  	    public static void main(String[] args)throws Exception{
	    	//Making an evaluation like described on the paper.
	    	Evaluation eval;
	    	Instances inst=null;
	    	String path;
	    	LRT tlr;
	    	double missing=0.3;
    	
	    	    path = args[0];
//	    		path = "C:\\Users\\Alex\\Desktop\\RankingDatasets\\";		
	    		tlr = new LRT();	
	    			    		
	    		for(int i=0; i<16; i++){
	    			switch(i){
	    			 
	    			case 0:
	    				inst = tlr.init(path, tlr, "analcatdata-authorship_dense.xarff");		    			
	    				System.out.println("Using authorship data set:");
	    				break;	 	    						
	    			case 1:
	    				inst = tlr.init(path, tlr, "bodyfat_dense.xarff");		
	    	    		System.out.println("Using bodyfat data set:");	
	    	    		break;	    		
	    			case 2:
	    				inst = tlr.init(path, tlr, "calhousing_dense.xarff");		
	    	    		System.out.println("Using calhousing data set:");	
	    				break;
	    			case 3:
	    				inst = tlr.init(path, tlr, "cpu-small_dense.xarff");		
	    	    		System.out.println("Using cpu-small data set:");	
	    	    		break;
	    			case 4:
	    				inst = tlr.init(path, tlr, "elevators_dense.xarff");		
	    	    		System.out.println("Using elevators data set:");	
	    	    		break;
	    			case 5:
	    				inst = tlr.init(path, tlr, "fried_dense.xarff");		
	    	    		System.out.println("Using fried data set:");	
	    				break;
	    			case 6:
	    				inst = tlr.init(path, tlr, "glass_dense.xarff");		
	    	    		System.out.println("Using glass data set:");	
	    	    		break;
	    			case 7:
	    				inst = tlr.init(path, tlr, "housing_dense.xarff");		
	    	    		System.out.println("Using housing data set:");
	    	    		break;
	    			case 8:
	    				inst = tlr.init(path, tlr, "iris_dense.xarff");		
	    	    		System.out.println("Using iris_dense data set:");	
	    	    		break;
	    			case 9:
	    				inst = tlr.init(path, tlr, "pendigits_dense.xarff");		
	    	    		System.out.println("Using pendigits data set:");	
	    	    		break;
	    			case 10:
	    				inst = tlr.init(path, tlr, "segment_dense.xarff");		
	    	    		System.out.println("Using segment data set:");	
	    	    		break;
	    			case 11:
	    				inst = tlr.init(path, tlr, "stock_dense.xarff");		
	    	    		System.out.println("Using stock data set:");	
	    	    		break;
	    			case 12:
	    				inst = tlr.init(path, tlr,"vehicle_dense.xarff");	
	    	    		System.out.println("Using vehicle data set:");
	    	    		break;
	    			case 13:
	    				inst = tlr.init(path, tlr,"vowel_dense.xarff");		
	    	    		System.out.println("Using vowel data set:");	
	    	    		break;
	    			case 14:
	    				inst = tlr.init(path,  tlr,"wine_dense.xarff");	
	    	    		System.out.println("Using wine data set:");	
	    				break;
	    			case 15:
	    				inst = tlr.init(path, tlr,"wisconsin_dense.xarff");		
	    	    		System.out.println("Using wisconsin data set:");	
	    	    		break;
	    	    	default: 
	    	    		continue;
	    			}
	    			eval = new Evaluation(inst);
	    			tlr.numVars = inst.numAttributes()-1;
	    			tlr.numLabels = inst.getLabels().size();
	    			double[] res = tlr.createStats(tlr,eval, inst,missing);
		    		System.out.println("kendall's tau: "+res[0]+"\t spearman correlation: "+res[1]+"\n\n");
	    			
	    		}
	    		    	
	    	}
	    	
	    	/**
	    	 * Computing kendall's tau and spearman correlation with 5 repeats.
	    	 * This evaluation corresponds to the one mentioned inside of
	    	 * the paper mentioned at the beginning of this file.
	    	 * @param LRT
	    	 * @param Evaluation
	    	 * @param Instances inst
	    	 * @param double missing
	    	 * @return double[]
	    	 * @throws Exception
	    	 */
	    	public double[] createStats(LRT tlr, Evaluation eval, Instances inst, double missing) throws Exception{

	    		int thisSeed=5;
	    		Instances tmp = new Instances(inst);
	    		
			   
	    		for(int i=1; i<thisSeed+1; i++){	    			   	    									    
	    				eval.crossValidateModel("weka.classifiers.labelranking.LRT", tmp, 10, tlr.getOptions(), new Random(i));	    					    			
	    		}
	    		return new double[]{eval.tau(), eval.spearman()[1]};
	
	    	}
	    	
	    	/**
	    	 * Setting up the instances and our classifier.
	    	 * @param d
	    	 * @param filename
	    	 * @throws Exception
	    	 */
	    	public Instances init(String path, LRT tlr, String filename) throws Exception{
	    		File file;
	    		XArffLoader xarff = new XArffLoader();
	    		file = new File(path+filename);
	    		xarff.setFile(file);
	    		Instances inst = xarff.getDataSet();
	    		inst.setClassIndex(inst.numAttributes()-1);

	    		tlr.buildClassifier(inst);
	    		
	    		return inst;
	    	}
	    	

 
  
  
  /**
   * Randomly delete labels based on a given missing rate.
   * The deleted labels are replaced by Double.NaN.
   * @param labels
   * @param missingrate
   * @throws Exception
   */
  public static FastVector labeldelete(FastVector labels, double missingrate, int numVars, long seed) throws Exception{
      FastVector labels_new = new FastVector();
      double[] label_new;
      Random r = new Random(seed);
     
      for(int i=0; i<labels.size(); i++){
          label_new = new double[((double[])(labels.elementAt(i))).length];
          double addedVars = 0;
          for(int j=0; j<label_new.length; j++){
            if (j < numVars)
              label_new[j] = ((double[])labels.elementAt(i))[j];
            else
              if(r.nextFloat()>=missingrate){
        	addedVars++;
                  label_new[j] = ((double[])labels.elementAt(i))[j];
              }else{
                  label_new[j] = -1;
              }
          }
          if (addedVars > 1)
            labels_new.addElement(label_new);
      }
      return labels_new;
  }
  
  /**
   * Randomly delete label based on a given missing rate.
   * The deleted labels are replaced by Double.NaN.
   * This method is an adaption of the original labeldelete
   * function included in here. It is used for evaluation. 
   * @param labels
   * @param missingrate
   * @throws Exception
   */
  public static Instances labeldeleteInst(Instances labels, double missingrate, int numVars, int seeding) throws Exception{
	  Instances labels_new = new Instances(labels);
	  labels_new.clear();
	  int key=0;
	  int count=0;
      double[] label_new;
      double[] readRank;
      double[] rank = new double[labels.getLabels().size()];
      
      Random r = new Random(seeding);
      int[][] pairList = new int[rank.length][rank.length];
      for(int i=0; i<labels.size(); i++){
    	     	  
    	  PreferenceDenseInstance pd = (PreferenceDenseInstance) labels.get(i);
    	  for(int z=0; ;z++){
    		  if(pd.getHashMap().get(z)!=null){
    			  RankUtilities.noPartialRanking(pd.getHashMap().get(z));
    			  readRank = RankUtilities.triangleLabels;
    			  key = z;
    			  break;
    		  }
    	  }
    	  
    	  
          label_new = new double[labels.get(i).numAttributes()];
          double addedVars = 0;
          pairList = new int[rank.length][rank.length];
          count=0;
          for(int j=0; j<label_new.length+labels.getLabels().size(); j++){
            if (j < numVars+1)
              label_new[j] = labels.get(i).value(j);
            else     
              if(r.nextFloat()>=missingrate){
            	  rank[count]=readRank[count];
            	  count++;
            	  addedVars++;       		  
              }else{
            		 rank[count] = -1;
            		 count++;
            	  }
              }
                       
          
          for(int y=0; y<pairList.length; y++){
          	for(int z=0; z<pairList.length; z++){
          		if(y==z)continue;
          		if(rank[y]!=-1&&rank[z]!=-1){
          			if(rank[y]<rank[z])pairList[y][z]=1;
          		}
          	}
          }
          
          if (addedVars > 1){
        	  HashMap<Integer,int[][]> hm = new HashMap<Integer,int[][]>();
        	  hm.put(key, pairList);

        	  //If the first label has been deleted, we must determine the new high order label
        	  //for compatibility reasons.
        	  RankUtilities.noPartialRanking(pairList);
        	  double[] tmp = RankUtilities.triangleLabels;
        	  label_new[label_new.length-1]=tmp[0];
           	  PreferenceDenseInstance pdi = new PreferenceDenseInstance(1, label_new, hm);
           	  labels_new.add(pdi);
          }
          
      }
      return labels_new;
  }


@Override
/**
 * Initializing the classifier by setting necessary standard values.
 * @param the data set.
 */
public void buildClassifier(Instances data) throws Exception {
	FastVector trainingData = new FastVector();
	instData = data;
	
	for(int i=0; i<data.numInstances(); i++){
		double[] inst = new double[data.numAttributes()+data.labels.size()-1];
		for(int j=0; j<data.numAttributes()-1; j++){
				inst[j]=data.get(i).value(j);				
		}
		
		PreferenceDenseInstance pdi;
		if(data.get(i)instanceof DenseInstance){		
			pdi = new PreferenceDenseInstance((DenseInstance)data.get(i));
		}
		else
			pdi = (PreferenceDenseInstance) data.get(i);
		
		int[][] preferences;
		for(int k=0; ;k++){
			if(pdi.getHashMap().get(k)!=null){
				preferences = pdi.getHashMap().get(k);
				break;
			}
		}
		if(RankUtilities.noPartialRanking(preferences)){
			int index=0;
			for(int l=data.numAttributes()-1; l<inst.length; l++){
				//TreeLabelRanking uses labels counting from 1 on, not 0.
				inst[l]=RankUtilities.triangleLabels[index]+1;
				index++;
			}
		}
		
		trainingData.add(inst);		
	}
	
	LRT tlr = new LRT();
    tlr.numLabels = data.labels.size();
    tlr.numVars = data.numAttributes()-1;
    LRTree lrTree = tlr.new LRTree();
    
    lrTree.buildModel(trainingData);
    tl = tlr;
    lrTr = lrTree;
 
}

/**
 * Capabilities necessary for the ranker.
 */
 public Capabilities getCapabilities() {
	    Capabilities result = super.getCapabilities();   // returns the object from weka.classifiers.Classifier
	 
	    result.disableAll();
	    // attributes
	    result.enableAllAttributes();
	 
	    // class   
	    result.enable(Capability.MISSING_CLASS_VALUES);	    
	    result.enable(Capability.RANKING);
	 
	    return result;
}
 
 /**
  * Displays a short description of the classifier.
  * @return
  */
 public String globalInfo(){
	 return "Decision tree based label lanking";
 }
 
 public String toString(){
	 try {
		return "decision tree based label ranking\n\n"+lrTr.toString();
	} catch (Exception e) {
		return "decision tree based label ranking, no model built yet.\n\n";
	}
 }
 
 /**
  * Inside of the output of distributionForInstance, we store the Label Rankings
  * predicted by the Ranking algorithm, so we will not have distributions over
  * the separate labels like in nominal case.
  * @param a single Instance of a data set.
  * @return a predicted Label Ranking of that Instance.
  */
 public double[] distributionForInstance(Instance inst) throws Exception{
	 
	 PreferenceDenseInstance pdi = (PreferenceDenseInstance) inst;
	 double[] features = new double[pdi.numAttributes()-1];
	 double[] ranking;
	 int[] prediction;
	 
	 for(int i=0; i<pdi.numAttributes()-1; i++){
		 features[i]=pdi.value(i);
	 }
	 
	 ranking = lrTr.getLabelRanking(features);
	 prediction = new int[ranking.length];
	 
	 //The LabelRanking Tree's rankings are based on 1, our ranks started with
	 //0, so we have to take this into consideration here.
	 for(int j=0; j<ranking.length; j++){
		 prediction[j]=(int)ranking[j]-1;
		 ranking[j]--;
	 }
	 
	 int key=0;
	 for(int i=0; ;i++){
		 if(pdi.getHashMap().get(i)!=null){
			 if(RankUtilities.noPartialRanking(pdi.getHashMap().get(i))){
				 key=i;
			 }
			 break;
		 }
	 }
	 
	 pdi.insertPrediction(key,prediction);
	 
	 return ranking;
 }

 
 /**
  * Returns an enumeration describing the available options.
  * 
  * @return an enumeration of all the available options.
  */
 public Enumeration<Option> listOptions() {
   Vector<Option> newVector = new Vector<Option>(2);
      
   Enumeration<Option> enu = super.listOptions();
   while (enu.hasMoreElements()) {
     newVector.addElement(enu.nextElement());
   }
   
   return newVector.elements();
 }

 /**
  * Parses a given list of options. <p/>
  * Options after -- are passed to the designated classifier.<p>
  *
  * @param options the list of options as an array of strings
  * @throws Exception if an option is not supported
  */
 public void setOptions(String[] options) throws Exception {
   super.setOptions(options);
 }

 /**
  * Gets the current settings of the Classifier.
  *
  * @return an array of strings suitable for passing to setOptions
 * @throws Exception 
 * @throws NumberFormatException 
  */
 public String[] getOptions(){
 	   String[] superOptions = super.getOptions();
 	   String[] options = new String[superOptions.length + 2];
 	    
 	    int current = 0;
 	    	    	     	    	
 	    System.arraycopy(superOptions, 0, options, current, 
 	        superOptions.length);

 	    current += superOptions.length;
 	    while (current < options.length) {
 	      options[current++] = "";
 	    }

 	   try {
 		for(int x=0; x<options.length; x++){
 			if(options[x].equals("-S")){
 				LRT.seed=Integer.parseInt(options[x+1]);
 				break;
 			}
 		}

	}catch (Exception e) {
		e.printStackTrace();
	}
 	    return options;
 }

  

}

