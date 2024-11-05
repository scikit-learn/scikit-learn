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
 *    REPTree.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.trees;

import weka.classifiers.Classifier;
import weka.classifiers.AbstractClassifier;
import weka.classifiers.Sourcable;
import weka.classifiers.rules.ZeroR;
import weka.core.AdditionalMeasureProducer;
import weka.core.Attribute;
import weka.core.Capabilities;
import weka.core.ContingencyTables;
import weka.core.Drawable;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.core.WeightedInstancesHandler;
import weka.core.Capabilities.Capability;

import java.io.Serializable;
import java.util.Enumeration;
import java.util.Random;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * Fast decision tree learner. Builds a decision/regression tree using information gain/variance and prunes it using reduced-error pruning (with backfitting).  Only sorts values for numeric attributes once. Missing values are dealt with by splitting the corresponding instances into pieces (i.e. as in C4.5).
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -M &lt;minimum number of instances&gt;
 *  Set minimum number of instances per leaf (default 2).</pre>
 * 
 * <pre> -V &lt;minimum variance for split&gt;
 *  Set minimum numeric class variance proportion
 *  of train variance for split (default 1e-3).</pre>
 * 
 * <pre> -N &lt;number of folds&gt;
 *  Number of folds for reduced error pruning (default 3).</pre>
 * 
 * <pre> -S &lt;seed&gt;
 *  Seed for random data shuffling (default 1).</pre>
 * 
 * <pre> -P
 *  No pruning.</pre>
 * 
 * <pre> -L
 *  Maximum tree depth (default -1, no maximum)</pre>
 * 
 <!-- options-end -->
 *
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version $Revision: 6746 $ 
 */
public class REPTree 
  extends AbstractClassifier 
  implements OptionHandler, WeightedInstancesHandler, Drawable, 
	     AdditionalMeasureProducer, Sourcable {

  /** for serialization */
  static final long serialVersionUID = -8562443428621539458L;
  
  /** ZeroR model that is used if no attributes are present. */
  protected ZeroR m_zeroR;

  /**
   * Returns a string describing classifier
   * @return a description suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {

    return  "Fast decision tree learner. Builds a decision/regression tree using "
      + "information gain/variance and prunes it using reduced-error pruning "
      + "(with backfitting).  Only sorts values for numeric attributes "
      + "once. Missing values are dealt with by splitting the corresponding "
      + "instances into pieces (i.e. as in C4.5).";
  }

  /** An inner class for building and storing the tree structure */
  protected class Tree 
    implements Serializable, RevisionHandler {
    
    /** for serialization */
    static final long serialVersionUID = -1635481717888437935L;
    
    /** The header information (for printing the tree). */
    protected Instances m_Info = null;

    /** The subtrees of this tree. */ 
    protected Tree[] m_Successors;
    
    /** The attribute to split on. */
    protected int m_Attribute = -1;

    /** The split point. */
    protected double m_SplitPoint = Double.NaN;
    
    /** The proportions of training instances going down each branch. */
    protected double[] m_Prop = null;

    /** Class probabilities from the training data in the nominal case. 
	Holds the mean in the numeric case. */
    protected double[] m_ClassProbs = null;
    
    /** The (unnormalized) class distribution in the nominal
	case. Holds the sum of squared errors and the weight 
	in the numeric case. */
    protected double[] m_Distribution = null;
    
    /** Class distribution of hold-out set at node in the nominal case. 
	Straight sum of weights in the numeric case (i.e. array has
	only one element. */
    protected double[] m_HoldOutDist = null;
    
    /** The hold-out error of the node. The number of miss-classified
	instances in the nominal case, the sum of squared errors in the 
	numeric case. */
    protected double m_HoldOutError = 0;
  
    /**
     * Computes class distribution of an instance using the tree.
     * 
     * @param instance the instance to compute the distribution for
     * @return the distribution
     * @throws Exception if computation fails
     */
    protected double[] distributionForInstance(Instance instance) 
      throws Exception {

      double[] returnedDist = null;
      
      if (m_Attribute > -1) {
	
	// Node is not a leaf
	if (instance.isMissing(m_Attribute)) {
	  
	  // Value is missing
	  returnedDist = new double[m_Info.numClasses()];

	  // Split instance up
	  for (int i = 0; i < m_Successors.length; i++) {
	    double[] help = 
	      m_Successors[i].distributionForInstance(instance);
	    if (help != null) {
	      for (int j = 0; j < help.length; j++) {
		returnedDist[j] += m_Prop[i] * help[j];
	      }
	    }
	  }
	} else if (m_Info.attribute(m_Attribute).isNominal() || m_Info.attribute(m_Attribute).isRanking()) {
	  
	  // For nominal attributes
	  returnedDist =  m_Successors[(int)instance.value(m_Attribute)].
	    distributionForInstance(instance);
	} else {
	  
	  // For numeric attributes
	  if (instance.value(m_Attribute) < m_SplitPoint) {
	    returnedDist = 
	      m_Successors[0].distributionForInstance(instance);
	  } else {
	    returnedDist = 
	      m_Successors[1].distributionForInstance(instance);
	  }
	}
      }
      if ((m_Attribute == -1) || (returnedDist == null)) {
	
	// Node is a leaf or successor is empty
	return m_ClassProbs;
      } else {
	return returnedDist;
      }
    }

   /**
    * Returns a string containing java source code equivalent to the test
    * made at this node. The instance being tested is called "i". This
    * routine assumes to be called in the order of branching, enabling us to
    * set the >= condition test (the last one) of a numeric splitpoint 
    * to just "true" (because being there in the flow implies that the 
    * previous less-than test failed).
    *
    * @param index index of the value tested
    * @return a value of type 'String'
    */
    public final String sourceExpression(int index) {
      
      StringBuffer expr = null;
      if (index < 0) {
        return "i[" + m_Attribute + "] == null";
      }
      if (m_Info.attribute(m_Attribute).isNominal() || m_Info.attribute(m_Attribute).isRanking()) {
        expr = new StringBuffer("i[");
	expr.append(m_Attribute).append("]");
	expr.append(".equals(\"").append(m_Info.attribute(m_Attribute)
		.value(index)).append("\")");
      } else {
        expr = new StringBuffer("");
	if (index == 0) {
	  expr.append("((Double)i[")
	    .append(m_Attribute).append("]).doubleValue() < ")
	    .append(m_SplitPoint);
	} else {
	  expr.append("true");
	}
      }
      return expr.toString();
    }

   /**
    * Returns source code for the tree as if-then statements. The 
    * class is assigned to variable "p", and assumes the tested 
    * instance is named "i". The results are returned as two stringbuffers: 
    * a section of code for assignment of the class, and a section of
    * code containing support code (eg: other support methods).
    * <p/>
    * TODO: If the outputted source code encounters a missing value
    * for the evaluated attribute, it stops branching and uses the 
    * class distribution of the current node to decide the return value. 
    * This is unlike the behaviour of distributionForInstance(). 
    *
    * @param className the classname that this static classifier has
    * @param parent parent node of the current node 
    * @return an array containing two stringbuffers, the first string containing
    * assignment code, and the second containing source for support code.
    * @throws Exception if something goes wrong
    */
    public StringBuffer [] toSource(String className, Tree parent) 
      throws Exception {
    
      StringBuffer [] result = new StringBuffer[2];
      double[] currentProbs;

      if(m_ClassProbs == null)
        currentProbs = parent.m_ClassProbs;
      else
        currentProbs = m_ClassProbs;

      long printID = nextID();

      // Is this a leaf?
      if (m_Attribute == -1) {
        result[0] = new StringBuffer("	p = ");
 	if(m_Info.classAttribute().isNumeric())
	  result[0].append(currentProbs[0]);
	else {
	  result[0].append(Utils.maxIndex(currentProbs));
	}
	result[0].append(";\n");
	result[1] = new StringBuffer("");
      } else {
	StringBuffer text = new StringBuffer("");
	StringBuffer atEnd = new StringBuffer("");

	text.append("  static double N")
	  .append(Integer.toHexString(this.hashCode()) + printID)
	  .append("(Object []i) {\n")
	  .append("    double p = Double.NaN;\n");

        text.append("    /* " + m_Info.attribute(m_Attribute).name() + " */\n");
	// Missing attribute?
	text.append("    if (" + this.sourceExpression(-1) + ") {\n")
	  .append("      p = ");
	if(m_Info.classAttribute().isNumeric())
	  text.append(currentProbs[0] + ";\n");
	else
	  text.append(Utils.maxIndex(currentProbs) + ";\n");
	text.append("    } ");
	
	// Branching of the tree
	for (int i=0;i<m_Successors.length; i++) {
          text.append("else if (" + this.sourceExpression(i) + ") {\n");
	  // Is the successor a leaf?
	  if(m_Successors[i].m_Attribute == -1) {
	    double[] successorProbs = m_Successors[i].m_ClassProbs;
	    if(successorProbs == null)
	      successorProbs = m_ClassProbs;
	    text.append("      p = ");
	    if(m_Info.classAttribute().isNumeric()) {
	      text.append(successorProbs[0] + ";\n");
	    } else {
	      text.append(Utils.maxIndex(successorProbs) + ";\n");
	    }
	  } else {
	    StringBuffer [] sub = m_Successors[i].toSource(className, this);
	    text.append("" + sub[0]);
            atEnd.append("" + sub[1]);
	  }
	  text.append("    } ");
	  if (i == m_Successors.length - 1) {
	    text.append("\n");
	  }
        }

        text.append("    return p;\n  }\n");

        result[0] = new StringBuffer("    p = " + className + ".N");
        result[0].append(Integer.toHexString(this.hashCode()) + printID)
          .append("(i);\n");
        result[1] = text.append("" + atEnd);
      }
      return result;
    }

	
    /**
     * Outputs one node for graph.
     * 
     * @param text the buffer to append the output to
     * @param num the current node id
     * @param parent the parent of the nodes
     * @return the next node id
     * @throws Exception if something goes wrong
     */
    protected int toGraph(StringBuffer text, int num,
			Tree parent) throws Exception {
      
      num++;
      if (m_Attribute == -1) {
	text.append("N" + Integer.toHexString(Tree.this.hashCode()) +
		    " [label=\"" + num + leafString(parent) +"\"" +
		    "shape=box]\n");
      } else {
	text.append("N" + Integer.toHexString(Tree.this.hashCode()) +
		    " [label=\"" + num + ": " + 
		    m_Info.attribute(m_Attribute).name() + 
		    "\"]\n");
	for (int i = 0; i < m_Successors.length; i++) {
	  text.append("N" + Integer.toHexString(Tree.this.hashCode()) 
		      + "->" + 
		      "N" + 
		      Integer.toHexString(m_Successors[i].hashCode())  +
		      " [label=\"");
	  if (m_Info.attribute(m_Attribute).isNumeric()) {
	    if (i == 0) {
	      text.append(" < " +
			  Utils.doubleToString(m_SplitPoint, 2));
	    } else {
	      text.append(" >= " +
			  Utils.doubleToString(m_SplitPoint, 2));
	    }
	  } else {
	    text.append(" = " + m_Info.attribute(m_Attribute).value(i));
	  }
	  text.append("\"]\n");
	  num = m_Successors[i].toGraph(text, num, this);
	}
      }
      
      return num;
    }

    /**
     * Outputs description of a leaf node.
     * 
     * @param parent the parent of the node
     * @return the description of the node
     * @throws Exception if generation fails
     */
    protected String leafString(Tree parent) throws Exception {
    
      if (m_Info.classAttribute().isNumeric()) {
	double classMean;
	if (m_ClassProbs == null) {
	  classMean = parent.m_ClassProbs[0];
	} else {
	  classMean = m_ClassProbs[0];
	}
	StringBuffer buffer = new StringBuffer();
	buffer.append(" : " + Utils.doubleToString(classMean, 2));
	double avgError = 0;
	if (m_Distribution[1] > 0) {
	  avgError = m_Distribution[0] / m_Distribution[1];
	}
	buffer.append(" (" +
		      Utils.doubleToString(m_Distribution[1], 2) + "/" +
		      Utils.doubleToString(avgError, 2) 
		      + ")");
	avgError = 0;
	if (m_HoldOutDist[0] > 0) {
	  avgError = m_HoldOutError / m_HoldOutDist[0];
	}
	buffer.append(" [" +
		      Utils.doubleToString(m_HoldOutDist[0], 2) + "/" +
		      Utils.doubleToString(avgError, 2) 
		      + "]");
	return buffer.toString();
      } else { 
	int maxIndex;
	if (m_ClassProbs == null) {
	  maxIndex = Utils.maxIndex(parent.m_ClassProbs);
	} else {
	  maxIndex = Utils.maxIndex(m_ClassProbs);
	}
	return " : " + m_Info.classAttribute().value(maxIndex) + 
	  " (" + Utils.doubleToString(Utils.sum(m_Distribution), 2) + 
	  "/" + 
	  Utils.doubleToString((Utils.sum(m_Distribution) - 
				m_Distribution[maxIndex]), 2) + ")" +
	  " [" + Utils.doubleToString(Utils.sum(m_HoldOutDist), 2) + "/" + 
	  Utils.doubleToString((Utils.sum(m_HoldOutDist) - 
				m_HoldOutDist[maxIndex]), 2) + "]";
      }
    }
  
    /**
     * Recursively outputs the tree.
     * 
     * @param level the current level
     * @param parent the current parent
     * @return the generated substree
     */
    protected String toString(int level, Tree parent) {

      try {
	StringBuffer text = new StringBuffer();
      
	if (m_Attribute == -1) {
	
	  // Output leaf info
	  return leafString(parent);
	} else if (m_Info.attribute(m_Attribute).isNominal() || m_Info.attribute(m_Attribute).isRanking()) {
	
	  // For nominal attributes
	  for (int i = 0; i < m_Successors.length; i++) {
	    text.append("\n");
	    for (int j = 0; j < level; j++) {
	      text.append("|   ");
	    }
	    text.append(m_Info.attribute(m_Attribute).name() + " = " +
			m_Info.attribute(m_Attribute).value(i));
	    text.append(m_Successors[i].toString(level + 1, this));
	  }
	} else {
	
	  // For numeric attributes
	  text.append("\n");
	  for (int j = 0; j < level; j++) {
	    text.append("|   ");
	  }
	  text.append(m_Info.attribute(m_Attribute).name() + " < " +
		      Utils.doubleToString(m_SplitPoint, 2));
	  text.append(m_Successors[0].toString(level + 1, this));
	  text.append("\n");
	  for (int j = 0; j < level; j++) {
	    text.append("|   ");
	  }
	  text.append(m_Info.attribute(m_Attribute).name() + " >= " +
		      Utils.doubleToString(m_SplitPoint, 2));
	  text.append(m_Successors[1].toString(level + 1, this));
	}
      
	return text.toString();
      } catch (Exception e) {
	e.printStackTrace();
	return "Decision tree: tree can't be printed";
      }
    }     

    /**
     * Recursively generates a tree.
     * 
     * @param sortedIndices the sorted indices of the instances
     * @param weights the weights of the instances
     * @param data the data to work with
     * @param totalWeight
     * @param classProbs the class probabilities
     * @param header the header of the data
     * @param minNum the minimum number of instances in a leaf
     * @param minVariance
     * @param depth the current depth of the tree
     * @param maxDepth the maximum allowed depth of the tree
     * @throws Exception if generation fails
     */
    protected void buildTree(int[][][] sortedIndices, double[][][] weights,
			     Instances data, double totalWeight, 
			     double[] classProbs, Instances header,
			     double minNum, double minVariance,
			     int depth, int maxDepth) 
      throws Exception {
      
      // Store structure of dataset, set minimum number of instances
      // and make space for potential info from pruning data
      m_Info = header;
      m_HoldOutDist = new double[data.numClasses()];
	
      // Make leaf if there are no training instances
      int helpIndex = 0;
      if (data.classIndex() == 0) {
	helpIndex = 1;
      }
      if (sortedIndices[0][helpIndex].length == 0) {
	if (data.classAttribute().isNumeric()) {
	  m_Distribution = new double[2];
	} else {
	  m_Distribution = new double[data.numClasses()];
	}
	m_ClassProbs = null;
        sortedIndices[0] = null;
        weights[0] = null;
	return;
      }
      
      double priorVar = 0;
      if (data.classAttribute().isNumeric()) {

	// Compute prior variance
	double totalSum = 0, totalSumSquared = 0, totalSumOfWeights = 0; 
	for (int i = 0; i < sortedIndices[0][helpIndex].length; i++) {
	  Instance inst = data.instance(sortedIndices[0][helpIndex][i]);
	  totalSum += inst.classValue() * weights[0][helpIndex][i];
	  totalSumSquared += 
	    inst.classValue() * inst.classValue() * weights[0][helpIndex][i];
	  totalSumOfWeights += weights[0][helpIndex][i];
	}
	priorVar = singleVariance(totalSum, totalSumSquared, 
				  totalSumOfWeights);
      }

      // Check if node doesn't contain enough instances, is pure
      // or the maximum tree depth is reached
      m_ClassProbs = new double[classProbs.length];
      System.arraycopy(classProbs, 0, m_ClassProbs, 0, classProbs.length);
      if ((totalWeight < (2 * minNum)) ||

	  // Nominal case
	  ((data.classAttribute().isNominal() || data.classAttribute().isRanking()) &&
	   Utils.eq(m_ClassProbs[Utils.maxIndex(m_ClassProbs)],
		    Utils.sum(m_ClassProbs))) ||

	  // Numeric case
	  (data.classAttribute().isNumeric() && 
	   ((priorVar / totalWeight) < minVariance)) ||

	  // Check tree depth
	  ((m_MaxDepth >= 0) && (depth >= maxDepth))) {

	// Make leaf
	m_Attribute = -1;
	if (data.classAttribute().isNominal() || data.classAttribute().isRanking()) {

	  // Nominal case
	  m_Distribution = new double[m_ClassProbs.length];
	  for (int i = 0; i < m_ClassProbs.length; i++) {
	    m_Distribution[i] = m_ClassProbs[i];
	  }
	  Utils.normalize(m_ClassProbs);
	} else {

	  // Numeric case
	  m_Distribution = new double[2];
	  m_Distribution[0] = priorVar;
	  m_Distribution[1] = totalWeight;
	}
        sortedIndices[0] = null;
        weights[0] = null;
	return;
      }

      // Compute class distributions and value of splitting
      // criterion for each attribute
      double[] vals = new double[data.numAttributes()];
      double[][][] dists = new double[data.numAttributes()][0][0];
      double[][] props = new double[data.numAttributes()][0];
      double[][] totalSubsetWeights = new double[data.numAttributes()][0];
      double[] splits = new double[data.numAttributes()];
      if (data.classAttribute().isNominal() || data.classAttribute().isRanking()) { 

	// Nominal case
	for (int i = 0; i < data.numAttributes(); i++) {
	  if (i != data.classIndex()) {
	    splits[i] = distribution(props, dists, i, sortedIndices[0][i], 
				     weights[0][i], totalSubsetWeights, data);
	    vals[i] = gain(dists[i], priorVal(dists[i]));
	  }
	}
      } else {

	// Numeric case
	for (int i = 0; i < data.numAttributes(); i++) {
	  if (i != data.classIndex()) {
	    splits[i] = 
	      numericDistribution(props, dists, i, sortedIndices[0][i], 
				  weights[0][i], totalSubsetWeights, data, 
				  vals);
	  }
	}
      }

      // Find best attribute
      m_Attribute = Utils.maxIndex(vals);
      int numAttVals = dists[m_Attribute].length;

      // Check if there are at least two subsets with
      // required minimum number of instances
      int count = 0;
      for (int i = 0; i < numAttVals; i++) {
	if (totalSubsetWeights[m_Attribute][i] >= minNum) {
	  count++;
	}
	if (count > 1) {
	  break;
	}
      }

      // Any useful split found?
      if ((vals[m_Attribute] > 0) && (count > 1)) {

        // Set split point, proportions, and temp arrays
	m_SplitPoint = splits[m_Attribute];
	m_Prop = props[m_Attribute];
        double[][] attSubsetDists = dists[m_Attribute];
        double[] attTotalSubsetWeights = totalSubsetWeights[m_Attribute];

        // Release some memory before proceeding further
        vals = null;
        dists = null;
        props = null;
        totalSubsetWeights = null;
        splits = null;

	// Split data
	int[][][][] subsetIndices = 
	  new int[numAttVals][1][data.numAttributes()][0];
	double[][][][] subsetWeights = 
	  new double[numAttVals][1][data.numAttributes()][0];
	splitData(subsetIndices, subsetWeights, m_Attribute, m_SplitPoint, 
		  sortedIndices[0], weights[0], data);

        // Release memory
        sortedIndices[0] = null;
        weights[0] = null;

        // Build successors
	m_Successors = new Tree[numAttVals];
	for (int i = 0; i < numAttVals; i++) {
	  m_Successors[i] = new Tree();
	  m_Successors[i].
	    buildTree(subsetIndices[i], subsetWeights[i], 
		      data, attTotalSubsetWeights[i],
		      attSubsetDists[i], header, minNum, 
		      minVariance, depth + 1, maxDepth);

          // Release as much memory as we can
          attSubsetDists[i] = null;
	}
      } else {
      
	// Make leaf
	m_Attribute = -1;
        sortedIndices[0] = null;
        weights[0] = null;
      }

      // Normalize class counts
      if (data.classAttribute().isNominal() || data.classAttribute().isRanking()) {
	m_Distribution = new double[m_ClassProbs.length];
	for (int i = 0; i < m_ClassProbs.length; i++) {
	    m_Distribution[i] = m_ClassProbs[i];
	}
	Utils.normalize(m_ClassProbs);
      } else {
	m_Distribution = new double[2];
	m_Distribution[0] = priorVar;
	m_Distribution[1] = totalWeight;
      }
    }

    /**
     * Computes size of the tree.
     * 
     * @return the number of nodes
     */
    protected int numNodes() {
    
      if (m_Attribute == -1) {
	return 1;
      } else {
	int size = 1;
	for (int i = 0; i < m_Successors.length; i++) {
	  size += m_Successors[i].numNodes();
	}
	return size;
      }
    }

    /**
     * Splits instances into subsets.
     * 
     * @param subsetIndices the sorted indices in the subset
     * @param subsetWeights the weights of the subset
     * @param att the attribute index
     * @param splitPoint the split point for numeric attributes
     * @param sortedIndices the sorted indices of the whole set
     * @param weights the weights of the whole set
     * @param data the data to work with
     * @throws Exception if something goes wrong
     */
    protected void splitData(int[][][][] subsetIndices, 
			     double[][][][] subsetWeights,
			     int att, double splitPoint, 
			     int[][] sortedIndices, double[][] weights, 
			     Instances data) throws Exception {
    
      int j;
      int[] num;
   
      // For each attribute
      for (int i = 0; i < data.numAttributes(); i++) {
	if (i != data.classIndex()) {
	  if (data.attribute(att).isNominal() || data.attribute(att).isRanking()) {

	    // For nominal attributes
	    num = new int[data.attribute(att).numValues()];
	    for (int k = 0; k < num.length; k++) {
	      subsetIndices[k][0][i] = new int[sortedIndices[i].length];
	      subsetWeights[k][0][i] = new double[sortedIndices[i].length];
	    }
	    for (j = 0; j < sortedIndices[i].length; j++) {
	      Instance inst = data.instance(sortedIndices[i][j]);
	      if (inst.isMissing(att)) {

		// Split instance up
		for (int k = 0; k < num.length; k++) {
		  if (m_Prop[k] > 0) {
		    subsetIndices[k][0][i][num[k]] = sortedIndices[i][j];
		    subsetWeights[k][0][i][num[k]] = 
		      m_Prop[k] * weights[i][j];
		    num[k]++;
		  }
		}
	      } else {
		int subset = (int)inst.value(att);
		subsetIndices[subset][0][i][num[subset]] = 
		  sortedIndices[i][j];
		subsetWeights[subset][0][i][num[subset]] = weights[i][j];
		num[subset]++;
	      }
	    }
	  } else {

	    // For numeric attributes
	    num = new int[2];
	    for (int k = 0; k < 2; k++) {
	      subsetIndices[k][0][i] = new int[sortedIndices[i].length];
	      subsetWeights[k][0][i] = new double[weights[i].length];
	    }
	    for (j = 0; j < sortedIndices[i].length; j++) {
	      Instance inst = data.instance(sortedIndices[i][j]);
	      if (inst.isMissing(att)) {

		// Split instance up
		for (int k = 0; k < num.length; k++) {
		  if (m_Prop[k] > 0) {
		    subsetIndices[k][0][i][num[k]] = sortedIndices[i][j];
		    subsetWeights[k][0][i][num[k]] = 
		      m_Prop[k] * weights[i][j];
		    num[k]++;
		  }
		}
	      } else {
		int subset = (inst.value(att) < splitPoint) ? 0 : 1;
		subsetIndices[subset][0][i][num[subset]] = 
		  sortedIndices[i][j];
		subsetWeights[subset][0][i][num[subset]] = weights[i][j];
		num[subset]++;
	      } 
	    }
	  }
	
	  // Trim arrays
	  for (int k = 0; k < num.length; k++) {
	    int[] copy = new int[num[k]];
	    System.arraycopy(subsetIndices[k][0][i], 0, copy, 0, num[k]);
	    subsetIndices[k][0][i] = copy;
	    double[] copyWeights = new double[num[k]];
	    System.arraycopy(subsetWeights[k][0][i], 0,
			     copyWeights, 0, num[k]);
	    subsetWeights[k][0][i] = copyWeights;
	  }
	}
      }
    }

    /**
     * Computes class distribution for an attribute.
     * 
     * @param props
     * @param dists
     * @param att the attribute index
     * @param sortedIndices the sorted indices of the instances
     * @param weights the weights of the instances
     * @param subsetWeights the weights of the subset
     * @param data the data to work with
     * @return the split point
     * @throws Exception if computation fails
     */
    protected double distribution(double[][] props,
				  double[][][] dists, int att, 
				  int[] sortedIndices,
				  double[] weights, 
				  double[][] subsetWeights, 
				  Instances data) 
      throws Exception {

      double splitPoint = Double.NaN;
      Attribute attribute = data.attribute(att);
      double[][] dist = null;
      int i;

      if (attribute.isNominal() || attribute.isRanking()) {

	// For nominal attributes
	dist = new double[attribute.numValues()][data.numClasses()];
	for (i = 0; i < sortedIndices.length; i++) {
	  Instance inst = data.instance(sortedIndices[i]);
	  if (inst.isMissing(att)) {
	    break;
	  }
	  dist[(int)inst.value(att)][(int)inst.classValue()] += weights[i];
	}
      } else {

	// For numeric attributes
	double[][] currDist = new double[2][data.numClasses()];
	dist = new double[2][data.numClasses()];

	// Move all instances into second subset
	for (int j = 0; j < sortedIndices.length; j++) {
	  Instance inst = data.instance(sortedIndices[j]);
	  if (inst.isMissing(att)) {
	    break;
	  }
	  currDist[1][(int)inst.classValue()] += weights[j];
	}
	double priorVal = priorVal(currDist);
	System.arraycopy(currDist[1], 0, dist[1], 0, dist[1].length);

	// Try all possible split points
	double currSplit = data.instance(sortedIndices[0]).value(att);
	double currVal, bestVal = -Double.MAX_VALUE;
	for (i = 0; i < sortedIndices.length; i++) {
	  Instance inst = data.instance(sortedIndices[i]);
	  if (inst.isMissing(att)) {
	    break;
	  }
	  if (inst.value(att) > currSplit) {
	    currVal = gain(currDist, priorVal);
	    if (currVal > bestVal) {
	      bestVal = currVal;
	      splitPoint = (inst.value(att) + currSplit) / 2.0;
	      for (int j = 0; j < currDist.length; j++) {
		System.arraycopy(currDist[j], 0, dist[j], 0, 
				 dist[j].length);
	      }
	    } 
	  } 
	  currSplit = inst.value(att);
	  currDist[0][(int)inst.classValue()] += weights[i];
	  currDist[1][(int)inst.classValue()] -= weights[i];
	}
      }

      // Compute weights
      props[att] = new double[dist.length];
      for (int k = 0; k < props[att].length; k++) {
	props[att][k] = Utils.sum(dist[k]);
      }
      if (!(Utils.sum(props[att]) > 0)) {
	for (int k = 0; k < props[att].length; k++) {
	  props[att][k] = 1.0 / (double)props[att].length;
	}
      } else {
	Utils.normalize(props[att]);
      }
    
      // Distribute counts
      while (i < sortedIndices.length) {
	Instance inst = data.instance(sortedIndices[i]);
	for (int j = 0; j < dist.length; j++) {
	  dist[j][(int)inst.classValue()] += props[att][j] * weights[i];
	}
	i++;
      }

      // Compute subset weights
      subsetWeights[att] = new double[dist.length];
      for (int j = 0; j < dist.length; j++) {
	subsetWeights[att][j] += Utils.sum(dist[j]);
      }

      // Return distribution and split point
      dists[att] = dist;
      return splitPoint;
    }      

    /**
     * Computes class distribution for an attribute.
     * 
     * @param props
     * @param dists
     * @param att the attribute index
     * @param sortedIndices the sorted indices of the instances
     * @param weights the weights of the instances
     * @param subsetWeights the weights of the subset
     * @param data the data to work with
     * @param vals
     * @return the split point
     * @throws Exception if computation fails
     */
    protected double numericDistribution(double[][] props, 
					 double[][][] dists, int att, 
					 int[] sortedIndices,
					 double[] weights, 
					 double[][] subsetWeights, 
					 Instances data,
					 double[] vals) 
      throws Exception {

      double splitPoint = Double.NaN;
      Attribute attribute = data.attribute(att);
      double[][] dist = null;
      double[] sums = null;
      double[] sumSquared = null;
      double[] sumOfWeights = null;
      double totalSum = 0, totalSumSquared = 0, totalSumOfWeights = 0;

      int i;

      if (attribute.isNominal() || attribute.isRanking()) {

	// For nominal attributes
	sums = new double[attribute.numValues()];
        sumSquared = new double[attribute.numValues()];
	sumOfWeights = new double[attribute.numValues()];
	int attVal;
	for (i = 0; i < sortedIndices.length; i++) {
	  Instance inst = data.instance(sortedIndices[i]);
	  if (inst.isMissing(att)) {
	    break;
	  }
	  attVal = (int)inst.value(att);
	  sums[attVal] += inst.classValue() * weights[i];
	  sumSquared[attVal] += 
	    inst.classValue() * inst.classValue() * weights[i];
	  sumOfWeights[attVal] += weights[i];
	}
	totalSum = Utils.sum(sums);
	totalSumSquared = Utils.sum(sumSquared);
	totalSumOfWeights = Utils.sum(sumOfWeights);
      } else {

	// For numeric attributes
	sums = new double[2];
        sumSquared = new double[2];
	sumOfWeights = new double[2];
	double[] currSums = new double[2];
        double[] currSumSquared = new double[2];
	double[] currSumOfWeights = new double[2];

	// Move all instances into second subset
	for (int j = 0; j < sortedIndices.length; j++) {
	  Instance inst = data.instance(sortedIndices[j]);
	  if (inst.isMissing(att)) {
	    break;
	  }
	  currSums[1] += inst.classValue() * weights[j];
	  currSumSquared[1] += 
	    inst.classValue() * inst.classValue() * weights[j];
	  currSumOfWeights[1] += weights[j];
	  
	}
	totalSum = currSums[1];
	totalSumSquared = currSumSquared[1];
	totalSumOfWeights = currSumOfWeights[1];
	
	sums[1] = currSums[1];
	sumSquared[1] = currSumSquared[1];
	sumOfWeights[1] = currSumOfWeights[1];

	// Try all possible split points
	double currSplit = data.instance(sortedIndices[0]).value(att);
	double currVal, bestVal = Double.MAX_VALUE;
	for (i = 0; i < sortedIndices.length; i++) {
	  Instance inst = data.instance(sortedIndices[i]);
	  if (inst.isMissing(att)) {
	    break;
	  }
	  if (inst.value(att) > currSplit) {
	    currVal = variance(currSums, currSumSquared, currSumOfWeights);
	    if (currVal < bestVal) {
	      bestVal = currVal;
	      splitPoint = (inst.value(att) + currSplit) / 2.0;
	      for (int j = 0; j < 2; j++) {
		sums[j] = currSums[j];
		sumSquared[j] = currSumSquared[j];
		sumOfWeights[j] = currSumOfWeights[j];
	      }
	    } 
	  } 

	  currSplit = inst.value(att);

	  double classVal = inst.classValue() * weights[i];
	  double classValSquared = inst.classValue() * classVal;

	  currSums[0] += classVal;
	  currSumSquared[0] += classValSquared;
	  currSumOfWeights[0] += weights[i];

	  currSums[1] -= classVal;
	  currSumSquared[1] -= classValSquared;
	  currSumOfWeights[1] -= weights[i];
	}
      }

      // Compute weights
      props[att] = new double[sums.length];
      for (int k = 0; k < props[att].length; k++) {
	props[att][k] = sumOfWeights[k];
      }
      if (!(Utils.sum(props[att]) > 0)) {
	for (int k = 0; k < props[att].length; k++) {
	  props[att][k] = 1.0 / (double)props[att].length;
	}
      } else {
	Utils.normalize(props[att]);
      }
    
	
      // Distribute counts for missing values
      while (i < sortedIndices.length) {
	Instance inst = data.instance(sortedIndices[i]);
	for (int j = 0; j < sums.length; j++) {
	  sums[j] += props[att][j] * inst.classValue() * weights[i];
	  sumSquared[j] += props[att][j] * inst.classValue() * 
	    inst.classValue() * weights[i];
	  sumOfWeights[j] += props[att][j] * weights[i];
	}
	totalSum += inst.classValue() * weights[i];
	totalSumSquared += 
	  inst.classValue() * inst.classValue() * weights[i]; 
	totalSumOfWeights += weights[i];
	i++;
      }

      // Compute final distribution
      dist = new double[sums.length][data.numClasses()];
      for (int j = 0; j < sums.length; j++) {
	if (sumOfWeights[j] > 0) {
	  dist[j][0] = sums[j] / sumOfWeights[j];
	} else {
	  dist[j][0] = totalSum / totalSumOfWeights;
	}
      }
      
      // Compute variance gain
      double priorVar =
	singleVariance(totalSum, totalSumSquared, totalSumOfWeights);
      double var = variance(sums, sumSquared, sumOfWeights);
      double gain = priorVar - var;
      
      // Return distribution and split point
      subsetWeights[att] = sumOfWeights;
      dists[att] = dist;
      vals[att] = gain;
      return splitPoint;
    }      

    /**
     * Computes variance for subsets.
     * 
     * @param s
     * @param sS
     * @param sumOfWeights
     * @return the variance
     */
    protected double variance(double[] s, double[] sS, 
			    double[] sumOfWeights) {
      
      double var = 0;
      
      for (int i = 0; i < s.length; i++) {
	if (sumOfWeights[i] > 0) {
	  var += singleVariance(s[i], sS[i], sumOfWeights[i]);
	}
      }
      
      return var;
    }
    
    /** 
     * Computes the variance for a single set
     * 
     * @param s
     * @param sS
     * @param weight the weight
     * @return the variance
     */
    protected double singleVariance(double s, double sS, double weight) {
      
      return sS - ((s * s) / weight);
    }

    /**
     * Computes value of splitting criterion before split.
     * 
     * @param dist
     * @return the splitting criterion
     */
    protected double priorVal(double[][] dist) {

      return ContingencyTables.entropyOverColumns(dist);
    }

    /**
     * Computes value of splitting criterion after split.
     * 
     * @param dist
     * @param priorVal the splitting criterion
     * @return the gain after splitting
     */
    protected double gain(double[][] dist, double priorVal) {

      return priorVal - ContingencyTables.entropyConditionedOnRows(dist);
    }

    /**
     * Prunes the tree using the hold-out data (bottom-up).
     * 
     * @return the error
     * @throws Exception if pruning fails for some reason
     */
    protected double reducedErrorPrune() throws Exception {

      // Is node leaf ? 
      if (m_Attribute == -1) {
	return m_HoldOutError;
      }

      // Prune all sub trees
      double errorTree = 0;
      for (int i = 0; i < m_Successors.length; i++) {
	errorTree += m_Successors[i].reducedErrorPrune();
      }

      // Replace sub tree with leaf if error doesn't get worse
      if (errorTree >= m_HoldOutError) {
	m_Attribute = -1;
	m_Successors = null;
	return m_HoldOutError;
      } else {
	return errorTree;
      }
    }

    /**
     * Inserts hold-out set into tree.
     * 
     * @param data the data to insert
     * @throws Exception if something goes wrong
     */
    protected void insertHoldOutSet(Instances data) throws Exception {

      for (int i = 0; i < data.numInstances(); i++) {
	insertHoldOutInstance(data.instance(i), data.instance(i).weight(),
			      this);
      }
    }

    /**
     * Inserts an instance from the hold-out set into the tree.
     * 
     * @param inst the instance to insert
     * @param weight the weight of the instance
     * @param parent the parent of the node
     * @throws Exception if insertion fails
     */
    protected void insertHoldOutInstance(Instance inst, double weight, 
					 Tree parent) throws Exception {
      
      // Insert instance into hold-out class distribution
      if (inst.classAttribute().isNominal() || inst.classAttribute().isRanking()) {
	
	// Nominal case
	m_HoldOutDist[(int)inst.classValue()] += weight;
	int predictedClass = 0;
	if (m_ClassProbs == null) {
	  predictedClass = Utils.maxIndex(parent.m_ClassProbs);
	} else {
	  predictedClass = Utils.maxIndex(m_ClassProbs);
	}
	if (predictedClass != (int)inst.classValue()) {
	  m_HoldOutError += weight;
	}
      } else {
	
	// Numeric case
	m_HoldOutDist[0] += weight;
	double diff = 0;
	if (m_ClassProbs == null) {
	  diff = parent.m_ClassProbs[0] - inst.classValue();
	} else {
	  diff =  m_ClassProbs[0] - inst.classValue();
	}
	m_HoldOutError += diff * diff * weight;
      }	
      
      // The process is recursive
      if (m_Attribute != -1) {
	
	// If node is not a leaf
	if (inst.isMissing(m_Attribute)) {
	  
	  // Distribute instance
	  for (int i = 0; i < m_Successors.length; i++) {
	    if (m_Prop[i] > 0) {
	      m_Successors[i].insertHoldOutInstance(inst, weight * 
						    m_Prop[i], this);
	    }
	  }
	} else {
	  
	  if (m_Info.attribute(m_Attribute).isNominal() || m_Info.attribute(m_Attribute).isRanking()) {
	    
	    // Treat nominal attributes
	    m_Successors[(int)inst.value(m_Attribute)].
	      insertHoldOutInstance(inst, weight, this);
	  } else {
	    
	    // Treat numeric attributes
	    if (inst.value(m_Attribute) < m_SplitPoint) {
	      m_Successors[0].insertHoldOutInstance(inst, weight, this);
	    } else {
	      m_Successors[1].insertHoldOutInstance(inst, weight, this);
	    }
	  }
	}
      }
    }
  
    /**
     * Inserts hold-out set into tree.
     * 
     * @param data the data to insert
     * @throws Exception if insertion fails
     */
    protected void backfitHoldOutSet(Instances data) throws Exception {
      
      for (int i = 0; i < data.numInstances(); i++) {
	backfitHoldOutInstance(data.instance(i), data.instance(i).weight(),
			       this);
      }
    }
    
    /**
     * Inserts an instance from the hold-out set into the tree.
     * 
     * @param inst the instance to insert
     * @param weight the weight of the instance
     * @param parent the parent node
     * @throws Exception if insertion fails
     */
    protected void backfitHoldOutInstance(Instance inst, double weight, 
					  Tree parent) throws Exception {
      
      // Insert instance into hold-out class distribution
      if (inst.classAttribute().isNominal() || inst.classAttribute().isRanking()) {
	
	// Nominal case
	if (m_ClassProbs == null) {
	  m_ClassProbs = new double[inst.numClasses()];
	}
	System.arraycopy(m_Distribution, 0, m_ClassProbs, 0, inst.numClasses());
	m_ClassProbs[(int)inst.classValue()] += weight;
	Utils.normalize(m_ClassProbs);
      } else {
	
	// Numeric case
	if (m_ClassProbs == null) {
	  m_ClassProbs = new double[1];
	}
	m_ClassProbs[0] *= m_Distribution[1];
	m_ClassProbs[0] += weight * inst.classValue();
	m_ClassProbs[0] /= (m_Distribution[1] + weight);
      }	
      
      // The process is recursive
      if (m_Attribute != -1) {
	
	// If node is not a leaf
	if (inst.isMissing(m_Attribute)) {
	  
	  // Distribute instance
	  for (int i = 0; i < m_Successors.length; i++) {
	    if (m_Prop[i] > 0) {
	      m_Successors[i].backfitHoldOutInstance(inst, weight * 
						     m_Prop[i], this);
	    }
	  }
	} else {
	  
	  if (m_Info.attribute(m_Attribute).isNominal() || m_Info.attribute(m_Attribute).isRanking()) {
	    
	    // Treat nominal attributes
	    m_Successors[(int)inst.value(m_Attribute)].
	      backfitHoldOutInstance(inst, weight, this);
	  } else {
	    
	    // Treat numeric attributes
	    if (inst.value(m_Attribute) < m_SplitPoint) {
	      m_Successors[0].backfitHoldOutInstance(inst, weight, this);
	    } else {
	      m_Successors[1].backfitHoldOutInstance(inst, weight, this);
	    }
	  }
	}
      }
    }
    
    /**
     * Returns the revision string.
     * 
     * @return		the revision
     */
    public String getRevision() {
      return RevisionUtils.extract("$Revision: 6746 $");
    }
  }

  /** The Tree object */
  protected Tree m_Tree = null;
    
  /** Number of folds for reduced error pruning. */
  protected int m_NumFolds = 3;
    
  /** Seed for random data shuffling. */
  protected int m_Seed = 1;
    
  /** Don't prune */
  protected boolean m_NoPruning = false;

  /** The minimum number of instances per leaf. */
  protected double m_MinNum = 2;

  /** The minimum proportion of the total variance (over all the data)
      required for split. */
  protected double m_MinVarianceProp = 1e-3;

  /** Upper bound on the tree depth */
  protected int m_MaxDepth = -1;
  
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String noPruningTipText() {
    return "Whether pruning is performed.";
  }
  
  /**
   * Get the value of NoPruning.
   *
   * @return Value of NoPruning.
   */
  public boolean getNoPruning() {
    
    return m_NoPruning;
  }
  
  /**
   * Set the value of NoPruning.
   *
   * @param newNoPruning Value to assign to NoPruning.
   */
  public void setNoPruning(boolean newNoPruning) {
    
    m_NoPruning = newNoPruning;
  }
  
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String minNumTipText() {
    return "The minimum total weight of the instances in a leaf.";
  }

  /**
   * Get the value of MinNum.
   *
   * @return Value of MinNum.
   */
  public double getMinNum() {
    
    return m_MinNum;
  }
  
  /**
   * Set the value of MinNum.
   *
   * @param newMinNum Value to assign to MinNum.
   */
  public void setMinNum(double newMinNum) {
    
    m_MinNum = newMinNum;
  }
  
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String minVariancePropTipText() {
    return "The minimum proportion of the variance on all the data " +
      "that needs to be present at a node in order for splitting to " +
      "be performed in regression trees.";
  }

  /**
   * Get the value of MinVarianceProp.
   *
   * @return Value of MinVarianceProp.
   */
  public double getMinVarianceProp() {
    
    return m_MinVarianceProp;
  }
  
  /**
   * Set the value of MinVarianceProp.
   *
   * @param newMinVarianceProp Value to assign to MinVarianceProp.
   */
  public void setMinVarianceProp(double newMinVarianceProp) {
    
    m_MinVarianceProp = newMinVarianceProp;
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String seedTipText() {
    return "The seed used for randomizing the data.";
  }

  /**
   * Get the value of Seed.
   *
   * @return Value of Seed.
   */
  public int getSeed() {
    
    return m_Seed;
  }
  
  /**
   * Set the value of Seed.
   *
   * @param newSeed Value to assign to Seed.
   */
  public void setSeed(int newSeed) {
    
    m_Seed = newSeed;
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String numFoldsTipText() {
    return "Determines the amount of data used for pruning. One fold is used for "
      + "pruning, the rest for growing the rules.";
  }
  
  /**
   * Get the value of NumFolds.
   *
   * @return Value of NumFolds.
   */
  public int getNumFolds() {
    
    return m_NumFolds;
  }
  
  /**
   * Set the value of NumFolds.
   *
   * @param newNumFolds Value to assign to NumFolds.
   */
  public void setNumFolds(int newNumFolds) {
    
    m_NumFolds = newNumFolds;
  }
  
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String maxDepthTipText() {
    return "The maximum tree depth (-1 for no restriction).";
  }

  /**
   * Get the value of MaxDepth.
   *
   * @return Value of MaxDepth.
   */
  public int getMaxDepth() {
    
    return m_MaxDepth;
  }
  
  /**
   * Set the value of MaxDepth.
   *
   * @param newMaxDepth Value to assign to MaxDepth.
   */
  public void setMaxDepth(int newMaxDepth) {
    
    m_MaxDepth = newMaxDepth;
  }
  
  /**
   * Lists the command-line options for this classifier.
   * 
   * @return an enumeration over all commandline options
   */
  public Enumeration listOptions() {
    
    Vector newVector = new Vector(5);

    newVector.
      addElement(new Option("\tSet minimum number of instances per leaf " +
			    "(default 2).",
			    "M", 1, "-M <minimum number of instances>"));
    newVector.
      addElement(new Option("\tSet minimum numeric class variance proportion\n" +
			    "\tof train variance for split (default 1e-3).",
			    "V", 1, "-V <minimum variance for split>"));
    newVector.
      addElement(new Option("\tNumber of folds for reduced error pruning " +
			    "(default 3).",
			    "N", 1, "-N <number of folds>"));
    newVector.
      addElement(new Option("\tSeed for random data shuffling (default 1).",
			    "S", 1, "-S <seed>"));
    newVector.
      addElement(new Option("\tNo pruning.",
			    "P", 0, "-P"));
    newVector.
      addElement(new Option("\tMaximum tree depth (default -1, no maximum)",
			    "L", 1, "-L"));

    return newVector.elements();
  } 

  /**
   * Gets options from this classifier.
   * 
   * @return the options for the current setup
   */
  public String[] getOptions() {
    
    String [] options = new String [12];
    int current = 0;
    options[current++] = "-M"; 
    options[current++] = "" + (int)getMinNum();
    options[current++] = "-V"; 
    options[current++] = "" + getMinVarianceProp();
    options[current++] = "-N"; 
    options[current++] = "" + getNumFolds();
    options[current++] = "-S"; 
    options[current++] = "" + getSeed();
    options[current++] = "-L"; 
    options[current++] = "" + getMaxDepth();
    if (getNoPruning()) {
      options[current++] = "-P";
    }
    while (current < options.length) {
      options[current++] = "";
    }
    return options;
  }

  /**
   * Parses a given list of options. <p/>
   * 
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -M &lt;minimum number of instances&gt;
   *  Set minimum number of instances per leaf (default 2).</pre>
   * 
   * <pre> -V &lt;minimum variance for split&gt;
   *  Set minimum numeric class variance proportion
   *  of train variance for split (default 1e-3).</pre>
   * 
   * <pre> -N &lt;number of folds&gt;
   *  Number of folds for reduced error pruning (default 3).</pre>
   * 
   * <pre> -S &lt;seed&gt;
   *  Seed for random data shuffling (default 1).</pre>
   * 
   * <pre> -P
   *  No pruning.</pre>
   * 
   * <pre> -L
   *  Maximum tree depth (default -1, no maximum)</pre>
   * 
   <!-- options-end -->
   * 
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    
    String minNumString = Utils.getOption('M', options);
    if (minNumString.length() != 0) {
      m_MinNum = (double)Integer.parseInt(minNumString);
    } else {
      m_MinNum = 2;
    }
    String minVarString = Utils.getOption('V', options);
    if (minVarString.length() != 0) {
      m_MinVarianceProp = Double.parseDouble(minVarString);
    } else {
      m_MinVarianceProp = 1e-3;
    }
    String numFoldsString = Utils.getOption('N', options);
    if (numFoldsString.length() != 0) {
      m_NumFolds = Integer.parseInt(numFoldsString);
    } else {
      m_NumFolds = 3;
    }
    String seedString = Utils.getOption('S', options);
    if (seedString.length() != 0) {
      m_Seed = Integer.parseInt(seedString);
    } else {
      m_Seed = 1;
    }
    m_NoPruning = Utils.getFlag('P', options);
    String depthString = Utils.getOption('L', options);
    if (depthString.length() != 0) {
      m_MaxDepth = Integer.parseInt(depthString);
    } else {
      m_MaxDepth = -1;
    }
    Utils.checkForRemainingOptions(options);
  }
  
  /**
   * Computes size of the tree.
   * 
   * @return the number of nodes
   */
  public int numNodes() {

    return m_Tree.numNodes();
  }

  /**
   * Returns an enumeration of the additional measure names.
   *
   * @return an enumeration of the measure names
   */
  public Enumeration enumerateMeasures() {
    
    Vector newVector = new Vector(1);
    newVector.addElement("measureTreeSize");
    return newVector.elements();
  }
 
  /**
   * Returns the value of the named measure.
   *
   * @param additionalMeasureName the name of the measure to query for its value
   * @return the value of the named measure
   * @throws IllegalArgumentException if the named measure is not supported
   */
  public double getMeasure(String additionalMeasureName) {
    
    if (additionalMeasureName.equalsIgnoreCase("measureTreeSize")) {
      return (double) numNodes();
    }
    else {throw new IllegalArgumentException(additionalMeasureName 
			      + " not supported (REPTree)");
    }
  }

  /**
   * Returns default capabilities of the classifier.
   *
   * @return      the capabilities of this classifier
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
    result.enable(Capability.NUMERIC_CLASS);
    result.enable(Capability.DATE_CLASS);
    result.enable(Capability.MISSING_CLASS_VALUES);
    
    return result;
  }

  /**
   * Builds classifier.
   * 
   * @param data the data to train with
   * @throws Exception if building fails
   */
  public void buildClassifier(Instances data) throws Exception {

    // can classifier handle the data?
    getCapabilities().testWithFail(data);

    // remove instances with missing class
    data = new Instances(data);
    data.deleteWithMissingClass();
    
    Random random = new Random(m_Seed);

    m_zeroR = null;
    if (data.numAttributes() == 1) {
      m_zeroR = new ZeroR();
      m_zeroR.buildClassifier(data);
      return;
    }

    // Randomize and stratify
    data.randomize(random);
    if (data.classAttribute().isNominal() || data.classAttribute().isRanking()) {
      data.stratify(m_NumFolds);
    }

    // Split data into training and pruning set
    Instances train = null;
    Instances prune = null;
    if (!m_NoPruning) {
      train = data.trainCV(m_NumFolds, 0, random);
      prune = data.testCV(m_NumFolds, 0);
    } else {
      train = data;
    }

    // Create array of sorted indices and weights
    int[][][] sortedIndices = new int[1][train.numAttributes()][0];
    double[][][] weights = new double[1][train.numAttributes()][0];
    double[] vals = new double[train.numInstances()];
    for (int j = 0; j < train.numAttributes(); j++) {
      if (j != train.classIndex()) {
	weights[0][j] = new double[train.numInstances()];
	if (train.attribute(j).isNominal() || train.attribute(j).isRanking()) {

	  // Handling nominal attributes. Putting indices of
	  // instances with missing values at the end.
	  sortedIndices[0][j] = new int[train.numInstances()];
	  int count = 0;
	  for (int i = 0; i < train.numInstances(); i++) {
	    Instance inst = train.instance(i);
	    if (!inst.isMissing(j)) {
	      sortedIndices[0][j][count] = i;
	      weights[0][j][count] = inst.weight();
	      count++;
	    }
	  }
	  for (int i = 0; i < train.numInstances(); i++) {
	    Instance inst = train.instance(i);
	    if (inst.isMissing(j)) {
	      sortedIndices[0][j][count] = i;
	      weights[0][j][count] = inst.weight();
	      count++;
	    }
	  }
	} else {

	  // Sorted indices are computed for numeric attributes
	  for (int i = 0; i < train.numInstances(); i++) {
	    Instance inst = train.instance(i);
	    vals[i] = inst.value(j);
	  }
	  sortedIndices[0][j] = Utils.sort(vals);
	  for (int i = 0; i < train.numInstances(); i++) {
	    weights[0][j][i] = train.instance(sortedIndices[0][j][i]).weight();
	  }
	}
      }
    }

    // Compute initial class counts
    double[] classProbs = new double[train.numClasses()];
    double totalWeight = 0, totalSumSquared = 0;
    for (int i = 0; i < train.numInstances(); i++) {
      Instance inst = train.instance(i);
      if (data.classAttribute().isNominal() || data.classAttribute().isRanking()) {
	classProbs[(int)inst.classValue()] += inst.weight();
	totalWeight += inst.weight();
      } else {
	classProbs[0] += inst.classValue() * inst.weight();
	totalSumSquared += inst.classValue() * inst.classValue() * inst.weight();
	totalWeight += inst.weight();
      }
    }
    m_Tree = new Tree();
    double trainVariance = 0;
    if (data.classAttribute().isNumeric()) {
      trainVariance = m_Tree.
	singleVariance(classProbs[0], totalSumSquared, totalWeight) / totalWeight;
      classProbs[0] /= totalWeight;
    }

    // Build tree
    m_Tree.buildTree(sortedIndices, weights, train, totalWeight, classProbs,
		     new Instances(train, 0), m_MinNum, m_MinVarianceProp * 
		     trainVariance, 0, m_MaxDepth);
    
    // Insert pruning data and perform reduced error pruning
    if (!m_NoPruning) {
      m_Tree.insertHoldOutSet(prune);
      m_Tree.reducedErrorPrune();
      m_Tree.backfitHoldOutSet(prune);
    }
  }

  /**
   * Computes class distribution of an instance using the tree.
   * 
   * @param instance the instance to compute the distribution for
   * @return the computed class probabilities
   * @throws Exception if computation fails
   */
  public double[] distributionForInstance(Instance instance) 
    throws Exception {
      
      if (m_zeroR != null) {
	return m_zeroR.distributionForInstance(instance);
      } else {
	return m_Tree.distributionForInstance(instance);
      }
  }


  /** 
   * For getting a unique ID when outputting the tree source
   * (hashcode isn't guaranteed unique) 
   */
  private static long PRINTED_NODES = 0;

  /**
   * Gets the next unique node ID.
   *
   * @return the next unique node ID.
   */
  protected static long nextID() {

    return PRINTED_NODES ++;
  }

  /**
   * resets the counter for the nodes
   */
  protected static void resetID() {
    PRINTED_NODES = 0;
  }

  /**
   * Returns the tree as if-then statements.
   *
   * @param className the name for the generated class
   * @return the tree as a Java if-then type statement
   * @throws Exception if something goes wrong
   */
  public String toSource(String className) 
    throws Exception {
     
    if (m_Tree == null) {
      throw new Exception("REPTree: No model built yet.");
    } 
    StringBuffer [] source = m_Tree.toSource(className, m_Tree);
    return
    "class " + className + " {\n\n"
    +"  public static double classify(Object [] i)\n"
    +"    throws Exception {\n\n"
    +"    double p = Double.NaN;\n"
    + source[0]  // Assignment code
    +"    return p;\n"
    +"  }\n"
    + source[1]  // Support code
    +"}\n";
  }

  /**
   *  Returns the type of graph this classifier
   *  represents.
   *  @return Drawable.TREE
   */   
  public int graphType() {
      return Drawable.TREE;
  }

  /**
   * Outputs the decision tree as a graph
   * 
   * @return the tree as a graph
   * @throws Exception if generation fails
   */
  public String graph() throws Exception {

    if (m_Tree == null) {
      throw new Exception("REPTree: No model built yet.");
    } 
    StringBuffer resultBuff = new StringBuffer();
    m_Tree.toGraph(resultBuff, 0, null);
    String result = "digraph Tree {\n" + "edge [style=bold]\n" + resultBuff.toString()
      + "\n}\n";
    return result;
  }
  
  /**
   * Outputs the decision tree.
   * 
   * @return a string representation of the classifier 
   */
  public String toString() {

    if (m_zeroR != null) {
      return "No attributes other than class. Using ZeroR.\n\n" + m_zeroR.toString();
    }
    if ((m_Tree == null)) {
      return "REPTree: No model built yet.";
    } 
    return     
      "\nREPTree\n============\n" + m_Tree.toString(0, null) + "\n" +
      "\nSize of the tree : " + numNodes();
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6746 $");
  }

  /**
   * Main method for this class.
   * 
   * @param argv the commandline options
   */
  public static void main(String[] argv) {
    runClassifier(new REPTree(), argv);
  }
}
