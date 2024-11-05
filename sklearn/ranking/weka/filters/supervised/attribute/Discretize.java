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
 *    Discretize.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.filters.supervised.attribute;

import weka.core.Attribute;
import weka.core.Capabilities;
import weka.core.ContingencyTables;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.DenseInstance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.Range;
import weka.core.RevisionUtils;
import weka.core.SparseInstance;
import weka.core.SpecialFunctions;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformationHandler;
import weka.core.Utils;
import weka.core.WeightedInstancesHandler;
import weka.core.Capabilities.Capability;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;
import weka.filters.Filter;
import weka.filters.SupervisedFilter;

import java.util.Enumeration;
import java.util.Vector;

/** 
 <!-- globalinfo-start -->
 * An instance filter that discretizes a range of numeric attributes in the dataset into nominal attributes. Discretization is by Fayyad &amp; Irani's MDL method (the default).<br/>
 * <br/>
 * For more information, see:<br/>
 * <br/>
 * Usama M. Fayyad, Keki B. Irani: Multi-interval discretization of continuousvalued attributes for classification learning. In: Thirteenth International Joint Conference on Articial Intelligence, 1022-1027, 1993.<br/>
 * <br/>
 * Igor Kononenko: On Biases in Estimating Multi-Valued Attributes. In: 14th International Joint Conference on Articial Intelligence, 1034-1040, 1995.
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;inproceedings{Fayyad1993,
 *    author = {Usama M. Fayyad and Keki B. Irani},
 *    booktitle = {Thirteenth International Joint Conference on Articial Intelligence},
 *    pages = {1022-1027},
 *    publisher = {Morgan Kaufmann Publishers},
 *    title = {Multi-interval discretization of continuousvalued attributes for classification learning},
 *    volume = {2},
 *    year = {1993}
 * }
 * 
 * &#64;inproceedings{Kononenko1995,
 *    author = {Igor Kononenko},
 *    booktitle = {14th International Joint Conference on Articial Intelligence},
 *    pages = {1034-1040},
 *    title = {On Biases in Estimating Multi-Valued Attributes},
 *    year = {1995},
 *    PS = {http://ai.fri.uni-lj.si/papers/kononenko95-ijcai.ps.gz}
 * }
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -R &lt;col1,col2-col4,...&gt;
 *  Specifies list of columns to Discretize. First and last are valid indexes.
 *  (default none)</pre>
 * 
 * <pre> -V
 *  Invert matching sense of column indexes.</pre>
 * 
 * <pre> -D
 *  Output binary attributes for discretized attributes.</pre>
 * 
 * <pre> -E
 *  Use better encoding of split point for MDL.</pre>
 * 
 * <pre> -K
 *  Use Kononenko's MDL criterion.</pre>
 * 
 <!-- options-end -->
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version $Revision: 6564 $
 */
public class Discretize 
  extends Filter 
  implements SupervisedFilter, OptionHandler, WeightedInstancesHandler, 
  	     TechnicalInformationHandler {
  
  /** for serialization */
  static final long serialVersionUID = -3141006402280129097L;

  /** Stores which columns to Discretize */
  protected Range m_DiscretizeCols = new Range();

  /** Store the current cutpoints */
  protected double [][] m_CutPoints = null;

  /** Output binary attributes for discretized attributes. */
  protected boolean m_MakeBinary = false;

  /** Use better encoding of split point for MDL. */
  protected boolean m_UseBetterEncoding = false;

  /** Use Kononenko's MDL criterion instead of Fayyad et al.'s */
  protected boolean m_UseKononenko = false;

  /** Constructor - initialises the filter */
  public Discretize() {

    setAttributeIndices("first-last");
  }


  /**
   * Gets an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {

    Vector newVector = new Vector(7);

    newVector.addElement(new Option(
              "\tSpecifies list of columns to Discretize. First"
	      + " and last are valid indexes.\n"
	      + "\t(default none)",
              "R", 1, "-R <col1,col2-col4,...>"));

    newVector.addElement(new Option(
              "\tInvert matching sense of column indexes.",
              "V", 0, "-V"));

    newVector.addElement(new Option(
              "\tOutput binary attributes for discretized attributes.",
              "D", 0, "-D"));

    newVector.addElement(new Option(
              "\tUse better encoding of split point for MDL.",
              "E", 0, "-E"));

    newVector.addElement(new Option(
              "\tUse Kononenko's MDL criterion.",
              "K", 0, "-K"));

    return newVector.elements();
  }


  /**
   * Parses a given list of options. <p/>
   * 
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -R &lt;col1,col2-col4,...&gt;
   *  Specifies list of columns to Discretize. First and last are valid indexes.
   *  (default none)</pre>
   * 
   * <pre> -V
   *  Invert matching sense of column indexes.</pre>
   * 
   * <pre> -D
   *  Output binary attributes for discretized attributes.</pre>
   * 
   * <pre> -E
   *  Use better encoding of split point for MDL.</pre>
   * 
   * <pre> -K
   *  Use Kononenko's MDL criterion.</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {

    setMakeBinary(Utils.getFlag('D', options));
    setUseBetterEncoding(Utils.getFlag('E', options));
    setUseKononenko(Utils.getFlag('K', options));
    setInvertSelection(Utils.getFlag('V', options));
    
    String convertList = Utils.getOption('R', options);
    if (convertList.length() != 0) {
      setAttributeIndices(convertList);
    } else {
      setAttributeIndices("first-last");
    }

    if (getInputFormat() != null) {
      setInputFormat(getInputFormat());
    }
  }
  /**
   * Gets the current settings of the filter.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String [] getOptions() {

    String [] options = new String [12];
    int current = 0;

    if (getMakeBinary()) {
      options[current++] = "-D";
    }
    if (getUseBetterEncoding()) {
      options[current++] = "-E";
    }
    if (getUseKononenko()) {
      options[current++] = "-K";
    }
    if (getInvertSelection()) {
      options[current++] = "-V";
    }
    if (!getAttributeIndices().equals("")) {
      options[current++] = "-R"; options[current++] = getAttributeIndices();
    }
    while (current < options.length) {
      options[current++] = "";
    }
    return options;
  }

  /** 
   * Returns the Capabilities of this filter.
   *
   * @return            the capabilities of this object
   * @see               Capabilities
   */
  public Capabilities getCapabilities() {
    Capabilities result = super.getCapabilities();
    result.disableAll();

    // attributes
    result.enableAllAttributes();
    result.enable(Capability.MISSING_VALUES);
    
    // class
    result.enable(Capability.NOMINAL_CLASS);

    
    return result;
  }

  /**
   * Sets the format of the input instances.
   *
   * @param instanceInfo an Instances object containing the input instance
   * structure (any instances contained in the object are ignored - only the
   * structure is required).
   * @return true if the outputFormat may be collected immediately
   * @throws Exception if the input format can't be set successfully
   */
  public boolean setInputFormat(Instances instanceInfo) throws Exception {

    super.setInputFormat(instanceInfo);

    m_DiscretizeCols.setUpper(instanceInfo.numAttributes() - 1);
    m_CutPoints = null;
    
    // If we implement loading cutfiles, then load 
    //them here and set the output format
    return false;
  }

  

  /**
   * Input an instance for filtering. Ordinarily the instance is processed
   * and made available for output immediately. Some filters require all
   * instances be read before producing output.
   *
   * @param instance the input instance
   * @return true if the filtered instance may now be
   * collected with output().
   * @throws IllegalStateException if no input format has been defined.
   */
  public boolean input(Instance instance) {

    if (getInputFormat() == null) {
      throw new IllegalStateException("No input instance format defined");
    }
    if (m_NewBatch) {
      resetQueue();
      m_NewBatch = false;
    }
    
    if (m_CutPoints != null) {
      convertInstance(instance);
      return true;
    }

    bufferInput(instance);
    return false;
  }


  /**
   * Signifies that this batch of input to the filter is finished. If the 
   * filter requires all instances prior to filtering, output() may now 
   * be called to retrieve the filtered instances.
   *
   * @return true if there are instances pending output
   * @throws IllegalStateException if no input structure has been defined
   */
  public boolean batchFinished() {

    if (getInputFormat() == null) {
      throw new IllegalStateException("No input instance format defined");
    }
    if (m_CutPoints == null) {
      calculateCutPoints();

      setOutputFormat();

      // If we implement saving cutfiles, save the cuts here

      // Convert pending input instances
      for(int i = 0; i < getInputFormat().numInstances(); i++) {
	convertInstance(getInputFormat().instance(i));
      }
    } 
    flushInput();

    m_NewBatch = true;
    return (numPendingOutput() != 0);
  }

  /**
   * Returns a string describing this filter
   *
   * @return a description of the filter suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {

    return "An instance filter that discretizes a range of numeric"
      + " attributes in the dataset into nominal attributes."
      + " Discretization is by Fayyad & Irani's MDL method (the default).\n\n"
      + "For more information, see:\n\n"
      + getTechnicalInformation().toString();
  }

  /**
   * Returns an instance of a TechnicalInformation object, containing 
   * detailed information about the technical background of this class,
   * e.g., paper reference or book this class is based on.
   * 
   * @return the technical information about this class
   */
  public TechnicalInformation getTechnicalInformation() {
    TechnicalInformation 	result;
    TechnicalInformation 	additional;
    
    result = new TechnicalInformation(Type.INPROCEEDINGS);
    result.setValue(Field.AUTHOR, "Usama M. Fayyad and Keki B. Irani");
    result.setValue(Field.TITLE, "Multi-interval discretization of continuousvalued attributes for classification learning");
    result.setValue(Field.BOOKTITLE, "Thirteenth International Joint Conference on Articial Intelligence");
    result.setValue(Field.YEAR, "1993");
    result.setValue(Field.VOLUME, "2");
    result.setValue(Field.PAGES, "1022-1027");
    result.setValue(Field.PUBLISHER, "Morgan Kaufmann Publishers");
    
    additional = result.add(Type.INPROCEEDINGS);
    additional.setValue(Field.AUTHOR, "Igor Kononenko");
    additional.setValue(Field.TITLE, "On Biases in Estimating Multi-Valued Attributes");
    additional.setValue(Field.BOOKTITLE, "14th International Joint Conference on Articial Intelligence");
    additional.setValue(Field.YEAR, "1995");
    additional.setValue(Field.PAGES, "1034-1040");
    additional.setValue(Field.PS, "http://ai.fri.uni-lj.si/papers/kononenko95-ijcai.ps.gz");
    
    return result;
  }
  
  /**
   * Returns the tip text for this property
   *
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String makeBinaryTipText() {

    return "Make resulting attributes binary.";
  }

  /**
   * Gets whether binary attributes should be made for discretized ones.
   *
   * @return true if attributes will be binarized
   */
  public boolean getMakeBinary() {

    return m_MakeBinary;
  }

  /** 
   * Sets whether binary attributes should be made for discretized ones.
   *
   * @param makeBinary if binary attributes are to be made
   */
  public void setMakeBinary(boolean makeBinary) {

    m_MakeBinary = makeBinary;
  }

  /**
   * Returns the tip text for this property
   *
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String useKononenkoTipText() {

    return "Use Kononenko's MDL criterion. If set to false"
      + " uses the Fayyad & Irani criterion.";
  }
  
  /**
   * Gets whether Kononenko's MDL criterion is to be used.
   *
   * @return true if Kononenko's criterion will be used.
   */
  public boolean getUseKononenko() {

    return m_UseKononenko;
  }

  /** 
   * Sets whether Kononenko's MDL criterion is to be used.
   *
   * @param useKon true if Kononenko's one is to be used
   */
  public void setUseKononenko(boolean useKon) {

    m_UseKononenko = useKon;
  }

  /**
   * Returns the tip text for this property
   *
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String useBetterEncodingTipText() {

    return "Uses a more efficient split point encoding.";
  }

  /**
   * Gets whether better encoding is to be used for MDL.
   *
   * @return true if the better MDL encoding will be used
   */
  public boolean getUseBetterEncoding() {

    return m_UseBetterEncoding;
  }

  /** 
   * Sets whether better encoding is to be used for MDL.
   *
   * @param useBetterEncoding true if better encoding to be used.
   */
  public void setUseBetterEncoding(boolean useBetterEncoding) {

    m_UseBetterEncoding = useBetterEncoding;
  }

  /**
   * Returns the tip text for this property
   *
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String invertSelectionTipText() {

    return "Set attribute selection mode. If false, only selected"
      + " (numeric) attributes in the range will be discretized; if"
      + " true, only non-selected attributes will be discretized.";
  }

  /**
   * Gets whether the supplied columns are to be removed or kept
   *
   * @return true if the supplied columns will be kept
   */
  public boolean getInvertSelection() {

    return m_DiscretizeCols.getInvert();
  }

  /**
   * Sets whether selected columns should be removed or kept. If true the 
   * selected columns are kept and unselected columns are deleted. If false
   * selected columns are deleted and unselected columns are kept.
   *
   * @param invert the new invert setting
   */
  public void setInvertSelection(boolean invert) {

    m_DiscretizeCols.setInvert(invert);
  }

  /**
   * Returns the tip text for this property
   *
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String attributeIndicesTipText() {
    return "Specify range of attributes to act on."
      + " This is a comma separated list of attribute indices, with"
      + " \"first\" and \"last\" valid values. Specify an inclusive"
      + " range with \"-\". E.g: \"first-3,5,6-10,last\".";
  }

  /**
   * Gets the current range selection
   *
   * @return a string containing a comma separated list of ranges
   */
  public String getAttributeIndices() {

    return m_DiscretizeCols.getRanges();
  }

  /**
   * Sets which attributes are to be Discretized (only numeric
   * attributes among the selection will be Discretized).
   *
   * @param rangeList a string representing the list of attributes. Since
   * the string will typically come from a user, attributes are indexed from
   * 1. <br>
   * eg: first-3,5,6-last
   * @throws IllegalArgumentException if an invalid range list is supplied 
   */
  public void setAttributeIndices(String rangeList) {

    m_DiscretizeCols.setRanges(rangeList);
  }

  /**
   * Sets which attributes are to be Discretized (only numeric
   * attributes among the selection will be Discretized).
   *
   * @param attributes an array containing indexes of attributes to Discretize.
   * Since the array will typically come from a program, attributes are indexed
   * from 0.
   * @throws IllegalArgumentException if an invalid set of ranges
   * is supplied 
   */
  public void setAttributeIndicesArray(int [] attributes) {

    setAttributeIndices(Range.indicesToRangeList(attributes));
  }

  /**
   * Gets the cut points for an attribute
   *
   * @param attributeIndex the index (from 0) of the attribute to get the cut points of
   * @return an array containing the cutpoints (or null if the
   * attribute requested isn't being Discretized
   */
  public double [] getCutPoints(int attributeIndex) {

    if (m_CutPoints == null) {
      return null;
    }
    return m_CutPoints[attributeIndex];
  }

  /** Generate the cutpoints for each attribute */
  protected void calculateCutPoints() {

    Instances copy = null;

    m_CutPoints = new double [getInputFormat().numAttributes()] [];
    for(int i = getInputFormat().numAttributes() - 1; i >= 0; i--) {
      if ((m_DiscretizeCols.isInRange(i)) && 
	  (getInputFormat().attribute(i).isNumeric())) {

	// Use copy to preserve order
	if (copy == null) {
	  copy = new Instances(getInputFormat());
	}
	calculateCutPointsByMDL(i, copy);
      }
    }
  }

  /**
   * Set cutpoints for a single attribute using MDL.
   *
   * @param index the index of the attribute to set cutpoints for
   * @param data the data to work with
   */
  protected void calculateCutPointsByMDL(int index,
					 Instances data) {

    // Sort instances
    data.sort(data.attribute(index));

    // Find first instances that's missing
    int firstMissing = data.numInstances();
    for (int i = 0; i < data.numInstances(); i++) {
      if (data.instance(i).isMissing(index)) {
        firstMissing = i;
        break;
      }
    }
    m_CutPoints[index] = cutPointsForSubset(data, index, 0, firstMissing);
  }

  /** 
   * Test using Kononenko's MDL criterion.
   * 
   * @param priorCounts
   * @param bestCounts
   * @param numInstances
   * @param numCutPoints
   * @return true if the split is acceptable
   */
  private boolean KononenkosMDL(double[] priorCounts,
				double[][] bestCounts,
				double numInstances,
				int numCutPoints) {

    double distPrior, instPrior, distAfter = 0, sum, instAfter = 0;
    double before, after;
    int numClassesTotal;

    // Number of classes occuring in the set
    numClassesTotal = 0;
    for (int i = 0; i < priorCounts.length; i++) {
      if (priorCounts[i] > 0) {
	numClassesTotal++;
      }
    }

    // Encode distribution prior to split
    distPrior = SpecialFunctions.log2Binomial(numInstances 
					      + numClassesTotal - 1,
					      numClassesTotal - 1);

    // Encode instances prior to split.
    instPrior = SpecialFunctions.log2Multinomial(numInstances,
						 priorCounts);

    before = instPrior + distPrior;

    // Encode distributions and instances after split.
    for (int i = 0; i < bestCounts.length; i++) {
      sum = Utils.sum(bestCounts[i]);
      distAfter += SpecialFunctions.log2Binomial(sum + numClassesTotal - 1,
						 numClassesTotal - 1);
      instAfter += SpecialFunctions.log2Multinomial(sum,
						    bestCounts[i]);
    }

    // Coding cost after split
    after = Utils.log2(numCutPoints) + distAfter + instAfter;

    // Check if split is to be accepted
    return (before > after);
  }


  /** 
   * Test using Fayyad and Irani's MDL criterion.
   * 
   * @param priorCounts
   * @param bestCounts
   * @param numInstances
   * @param numCutPoints
   * @return true if the splits is acceptable
   */
  private boolean FayyadAndIranisMDL(double[] priorCounts,
				     double[][] bestCounts,
				     double numInstances,
				     int numCutPoints) {

    double priorEntropy, entropy, gain; 
    double entropyLeft, entropyRight, delta;
    int numClassesTotal, numClassesRight, numClassesLeft;

    // Compute entropy before split.
    priorEntropy = ContingencyTables.entropy(priorCounts);

    // Compute entropy after split.
    entropy = ContingencyTables.entropyConditionedOnRows(bestCounts);

    // Compute information gain.
    gain = priorEntropy - entropy;

    // Number of classes occuring in the set
    numClassesTotal = 0;
    for (int i = 0; i < priorCounts.length; i++) {
      if (priorCounts[i] > 0) {
	numClassesTotal++;
      }
    }

    // Number of classes occuring in the left subset
    numClassesLeft = 0;
    for (int i = 0; i < bestCounts[0].length; i++) {
      if (bestCounts[0][i] > 0) {
	numClassesLeft++;
      }
    }

    // Number of classes occuring in the right subset
    numClassesRight = 0;
    for (int i = 0; i < bestCounts[1].length; i++) {
      if (bestCounts[1][i] > 0) {
	numClassesRight++;
      }
    }

    // Entropy of the left and the right subsets
    entropyLeft = ContingencyTables.entropy(bestCounts[0]);
    entropyRight = ContingencyTables.entropy(bestCounts[1]);

    // Compute terms for MDL formula
    delta = Utils.log2(Math.pow(3, numClassesTotal) - 2) - 
      (((double) numClassesTotal * priorEntropy) - 
       (numClassesRight * entropyRight) - 
       (numClassesLeft * entropyLeft));

    // Check if split is to be accepted
    return (gain > (Utils.log2(numCutPoints) + delta) / (double)numInstances);
  }
    

  /** 
   * Selects cutpoints for sorted subset.
   * 
   * @param instances
   * @param attIndex
   * @param first
   * @param lastPlusOne
   * @return
   */
  private double[] cutPointsForSubset(Instances instances, int attIndex, 
				      int first, int lastPlusOne) { 

    double[][] counts, bestCounts;
    double[] priorCounts, left, right, cutPoints;
    double currentCutPoint = -Double.MAX_VALUE, bestCutPoint = -1, 
      currentEntropy, bestEntropy, priorEntropy, gain;
    int bestIndex = -1, numInstances = 0, numCutPoints = 0;

    // Compute number of instances in set
    if ((lastPlusOne - first) < 2) {
      return null;
    }

    // Compute class counts.
    counts = new double[2][instances.numClasses()];
    for (int i = first; i < lastPlusOne; i++) {
      numInstances += instances.instance(i).weight();
      counts[1][(int)instances.instance(i).classValue()] +=
	instances.instance(i).weight();
    }

    // Save prior counts
    priorCounts = new double[instances.numClasses()];
    System.arraycopy(counts[1], 0, priorCounts, 0, 
		     instances.numClasses());

    // Entropy of the full set
    priorEntropy = ContingencyTables.entropy(priorCounts);
    bestEntropy = priorEntropy;
    
    // Find best entropy.
    bestCounts = new double[2][instances.numClasses()];
    for (int i = first; i < (lastPlusOne - 1); i++) {
      counts[0][(int)instances.instance(i).classValue()] +=
	instances.instance(i).weight();
      counts[1][(int)instances.instance(i).classValue()] -=
	instances.instance(i).weight();
      if (instances.instance(i).value(attIndex) < 
	  instances.instance(i + 1).value(attIndex)) {
	currentCutPoint = (instances.instance(i).value(attIndex) + 
	  instances.instance(i + 1).value(attIndex)) / 2.0;
	currentEntropy = ContingencyTables.entropyConditionedOnRows(counts);
	if (currentEntropy < bestEntropy) {
	  bestCutPoint = currentCutPoint;
	  bestEntropy = currentEntropy;
	  bestIndex = i;
	  System.arraycopy(counts[0], 0, 
			   bestCounts[0], 0, instances.numClasses());
	  System.arraycopy(counts[1], 0, 
			   bestCounts[1], 0, instances.numClasses()); 
	}
	numCutPoints++;
      }
    }

    // Use worse encoding?
    if (!m_UseBetterEncoding) {
      numCutPoints = (lastPlusOne - first) - 1;
    }

    // Checks if gain is zero
    gain = priorEntropy - bestEntropy;
    if (gain <= 0) {
      return null;
    }

    // Check if split is to be accepted
    if ((m_UseKononenko && KononenkosMDL(priorCounts, bestCounts,
					 numInstances, numCutPoints)) ||
	(!m_UseKononenko && FayyadAndIranisMDL(priorCounts, bestCounts,
					       numInstances, numCutPoints))) {
      
      // Select split points for the left and right subsets
      left = cutPointsForSubset(instances, attIndex, first, bestIndex + 1);
      right = cutPointsForSubset(instances, attIndex, 
				 bestIndex + 1, lastPlusOne);
      
      // Merge cutpoints and return them
      if ((left == null) && (right) == null) {
	cutPoints = new double[1];
	cutPoints[0] = bestCutPoint;
      } else if (right == null) {
	cutPoints = new double[left.length + 1];
	System.arraycopy(left, 0, cutPoints, 0, left.length);
	cutPoints[left.length] = bestCutPoint;
      } else if (left == null) {
	cutPoints = new double[1 + right.length];
	cutPoints[0] = bestCutPoint;
	System.arraycopy(right, 0, cutPoints, 1, right.length);
      } else {
	cutPoints = new double[left.length + right.length + 1];
	System.arraycopy(left, 0, cutPoints, 0, left.length);
	cutPoints[left.length] = bestCutPoint;
	System.arraycopy(right, 0, cutPoints, left.length + 1, right.length);
      }
      
      return cutPoints;
    } else
      return null;
  }
 
  /**
   * Set the output format. Takes the currently defined cutpoints and 
   * m_InputFormat and calls setOutputFormat(Instances) appropriately.
   */
  protected void setOutputFormat() {

    if (m_CutPoints == null) {
      setOutputFormat(null);
      return;
    }
    FastVector attributes = new FastVector(getInputFormat().numAttributes());
    int classIndex = getInputFormat().classIndex();
    for(int i = 0; i < getInputFormat().numAttributes(); i++) {
      if ((m_DiscretizeCols.isInRange(i)) 
	  && (getInputFormat().attribute(i).isNumeric())) {
	if (!m_MakeBinary) {
	  FastVector attribValues = new FastVector(1);
	  if (m_CutPoints[i] == null) {
	    attribValues.addElement("'All'");
	  } else {
	    for(int j = 0; j <= m_CutPoints[i].length; j++) {
	      if (j == 0) {
		attribValues.addElement("'(-inf-"
			+ Utils.doubleToString(m_CutPoints[i][j], 6) + "]'");
	      } else if (j == m_CutPoints[i].length) {
		attribValues.addElement("'("
			+ Utils.doubleToString(m_CutPoints[i][j - 1], 6) 
					+ "-inf)'");
	      } else {
		attribValues.addElement("'("
			+ Utils.doubleToString(m_CutPoints[i][j - 1], 6) + "-"
			+ Utils.doubleToString(m_CutPoints[i][j], 6) + "]'");
	      }
	    }
	  }
	  attributes.addElement(new Attribute(getInputFormat().
					      attribute(i).name(),
					      attribValues));
	} else {
	  if (m_CutPoints[i] == null) {
	    FastVector attribValues = new FastVector(1);
	    attribValues.addElement("'All'");
	    attributes.addElement(new Attribute(getInputFormat().
						attribute(i).name(),
						attribValues));
	  } else {
	    if (i < getInputFormat().classIndex()) {
	      classIndex += m_CutPoints[i].length - 1;
	    }
	    for(int j = 0; j < m_CutPoints[i].length; j++) {
	      FastVector attribValues = new FastVector(2);
	      attribValues.addElement("'(-inf-"
		      + Utils.doubleToString(m_CutPoints[i][j], 6) + "]'");
	      attribValues.addElement("'("
		      + Utils.doubleToString(m_CutPoints[i][j], 6) + "-inf)'");
	      attributes.addElement(new Attribute(getInputFormat().
						  attribute(i).name() + "_" + (j+1),
						  attribValues));
	    }
	  }
	}
      } else {
	attributes.addElement(getInputFormat().attribute(i).copy());
      }
    }
    Instances outputFormat = 
      new Instances(getInputFormat().relationName(), attributes, 0);
    outputFormat.setClassIndex(classIndex);
    setOutputFormat(outputFormat);
  }

  /**
   * Convert a single instance over. The converted instance is added to 
   * the end of the output queue.
   *
   * @param instance the instance to convert
   */
  protected void convertInstance(Instance instance) {

    int index = 0;
    double [] vals = new double [outputFormatPeek().numAttributes()];
    // Copy and convert the values
    for(int i = 0; i < getInputFormat().numAttributes(); i++) {
      if (m_DiscretizeCols.isInRange(i) && 
	  getInputFormat().attribute(i).isNumeric()) {
	int j;
	double currentVal = instance.value(i);
	if (m_CutPoints[i] == null) {
	  if (instance.isMissing(i)) {
	    vals[index] = Utils.missingValue();
	  } else {
	    vals[index] = 0;
	  }
	  index++;
	} else {
	  if (!m_MakeBinary) {
	    if (instance.isMissing(i)) {
	      vals[index] = Utils.missingValue();
	    } else {
	      for (j = 0; j < m_CutPoints[i].length; j++) {
		if (currentVal <= m_CutPoints[i][j]) {
		  break;
		}
	      }
              vals[index] = j;
	    }
	    index++;
	  } else {
	    for (j = 0; j < m_CutPoints[i].length; j++) {
	      if (instance.isMissing(i)) {
                vals[index] = Utils.missingValue();
	      } else if (currentVal <= m_CutPoints[i][j]) {
                vals[index] = 0;
	      } else {
                vals[index] = 1;
	      }
	      index++;
	    }
	  }   
	}
      } else {
        vals[index] = instance.value(i);
	index++;
      }
    }
    
    Instance inst = null;
    if (instance instanceof SparseInstance) {
      inst = new SparseInstance(instance.weight(), vals);
    } else {
      inst = new DenseInstance(instance.weight(), vals);
    }
    inst.setDataset(getOutputFormat());
    copyValues(inst, false, instance.dataset(), getOutputFormat());
    inst.setDataset(getOutputFormat());
    push(inst);
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6564 $");
  }

  /**
   * Main method for testing this class.
   *
   * @param argv should contain arguments to the filter: use -h for help
   */
  public static void main(String [] argv) {
    runFilter(new Discretize(), argv);
  }
}
