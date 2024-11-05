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


package weka.filters.unsupervised.attribute;

import weka.core.Attribute;
import weka.core.Capabilities;
import weka.core.FastVector;
import weka.core.Instance; 
import weka.core.DenseInstance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.Range;
import weka.core.RevisionUtils;
import weka.core.SparseInstance;
import weka.core.Utils;
import weka.core.WeightedInstancesHandler;
import weka.core.Capabilities.Capability;
import weka.core.labelranking.PreferenceDenseInstance;
import weka.filters.UnsupervisedFilter;

import java.util.Enumeration;
import java.util.Vector;

/** 
 <!-- globalinfo-start -->
 * An instance filter that discretizes a range of numeric attributes in the dataset into nominal attributes. Discretization is by simple binning. Skips the class attribute if set.
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -unset-class-temporarily
 *  Unsets the class index temporarily before the filter is
 *  applied to the data.
 *  (default: no)</pre>
 * 
 * <pre> -B &lt;num&gt;
 *  Specifies the (maximum) number of bins to divide numeric attributes into.
 *  (default = 10)</pre>
 * 
 * <pre> -M &lt;num&gt;
 *  Specifies the desired weight of instances per bin for
 *  equal-frequency binning. If this is set to a positive
 *  number then the -B option will be ignored.
 *  (default = -1)</pre>
 * 
 * <pre> -F
 *  Use equal-frequency instead of equal-width discretization.</pre>
 * 
 * <pre> -O
 *  Optimize number of bins using leave-one-out estimate
 *  of estimated entropy (for equal-width discretization).
 *  If this is set then the -B option will be ignored.</pre>
 * 
 * <pre> -R &lt;col1,col2-col4,...&gt;
 *  Specifies list of columns to Discretize. First and last are valid indexes.
 *  (default: first-last)</pre>
 * 
 * <pre> -V
 *  Invert matching sense of column indexes.</pre>
 * 
 * <pre> -D
 *  Output binary attributes for discretized attributes.</pre>
 * 
 <!-- options-end -->
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version $Revision: 6567 $
 */
public class Discretize 
  extends PotentialClassIgnorer 
  implements UnsupervisedFilter, WeightedInstancesHandler {
  
  /** for serialization */
  static final long serialVersionUID = -1358531742174527279L;

  /** Stores which columns to Discretize */
  protected Range m_DiscretizeCols = new Range();

  /** The number of bins to divide the attribute into */
  protected int m_NumBins = 10;

  /** The desired weight of instances per bin */
  protected double m_DesiredWeightOfInstancesPerInterval = -1;

  /** Store the current cutpoints */
  protected double [][] m_CutPoints = null;

  /** Output binary attributes for discretized attributes. */
  protected boolean m_MakeBinary = false;

  /** Find the number of bins using cross-validated entropy. */
  protected boolean m_FindNumBins = false;

  /** Use equal-frequency binning if unsupervised discretization turned on */
  protected boolean m_UseEqualFrequency = false;

  /** The default columns to discretize */
  protected String m_DefaultCols;

  /** Constructor - initialises the filter */
  public Discretize() {

    m_DefaultCols = "first-last";
    setAttributeIndices("first-last");
  }

  /** 
   * Another constructor, sets the attribute indices immediately
   * 
   * @param cols the attribute indices
   */
  public Discretize(String cols) {

    m_DefaultCols = cols;
    setAttributeIndices(cols);
  }

  /**
   * Gets an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector result = new Vector();
    Enumeration enm = super.listOptions();
    while (enm.hasMoreElements())
      result.add(enm.nextElement());
      
    result.addElement(new Option(
	"\tSpecifies the (maximum) number of bins to divide numeric"
	+ " attributes into.\n"
	+ "\t(default = 10)",
	"B", 1, "-B <num>"));
    
    result.addElement(new Option(
	"\tSpecifies the desired weight of instances per bin for\n"
	+ "\tequal-frequency binning. If this is set to a positive\n"
	+ "\tnumber then the -B option will be ignored.\n"
	+ "\t(default = -1)",
	"M", 1, "-M <num>"));
    
    result.addElement(new Option(
	"\tUse equal-frequency instead of equal-width discretization.",
	"F", 0, "-F"));
    
    result.addElement(new Option(
	"\tOptimize number of bins using leave-one-out estimate\n"+
	"\tof estimated entropy (for equal-width discretization).\n"+
	"\tIf this is set then the -B option will be ignored.",
	"O", 0, "-O"));
    
    result.addElement(new Option(
	"\tSpecifies list of columns to Discretize. First"
	+ " and last are valid indexes.\n"
	+ "\t(default: first-last)",
	"R", 1, "-R <col1,col2-col4,...>"));
    
    result.addElement(new Option(
	"\tInvert matching sense of column indexes.",
	"V", 0, "-V"));
    
    result.addElement(new Option(
	"\tOutput binary attributes for discretized attributes.",
	"D", 0, "-D"));

    return result.elements();
  }


  /**
   * Parses a given list of options. <p/>
   * 
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -unset-class-temporarily
   *  Unsets the class index temporarily before the filter is
   *  applied to the data.
   *  (default: no)</pre>
   * 
   * <pre> -B &lt;num&gt;
   *  Specifies the (maximum) number of bins to divide numeric attributes into.
   *  (default = 10)</pre>
   * 
   * <pre> -M &lt;num&gt;
   *  Specifies the desired weight of instances per bin for
   *  equal-frequency binning. If this is set to a positive
   *  number then the -B option will be ignored.
   *  (default = -1)</pre>
   * 
   * <pre> -F
   *  Use equal-frequency instead of equal-width discretization.</pre>
   * 
   * <pre> -O
   *  Optimize number of bins using leave-one-out estimate
   *  of estimated entropy (for equal-width discretization).
   *  If this is set then the -B option will be ignored.</pre>
   * 
   * <pre> -R &lt;col1,col2-col4,...&gt;
   *  Specifies list of columns to Discretize. First and last are valid indexes.
   *  (default: first-last)</pre>
   * 
   * <pre> -V
   *  Invert matching sense of column indexes.</pre>
   * 
   * <pre> -D
   *  Output binary attributes for discretized attributes.</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {

    super.setOptions(options);

    setMakeBinary(Utils.getFlag('D', options));
    setUseEqualFrequency(Utils.getFlag('F', options));
    setFindNumBins(Utils.getFlag('O', options));
    setInvertSelection(Utils.getFlag('V', options));

    String weight = Utils.getOption('M', options);
    if (weight.length() != 0) {
      setDesiredWeightOfInstancesPerInterval((new Double(weight)).doubleValue());
    } else {
      setDesiredWeightOfInstancesPerInterval(-1);
    }

    String numBins = Utils.getOption('B', options);
    if (numBins.length() != 0) {
      setBins(Integer.parseInt(numBins));
    } else {
      setBins(10);
    }
    
    String convertList = Utils.getOption('R', options);
    if (convertList.length() != 0) {
      setAttributeIndices(convertList);
    } else {
      setAttributeIndices(m_DefaultCols);
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
    Vector        result;
    String[]      options;
    int           i;

    result = new Vector();

    options = super.getOptions();
    for (i = 0; i < options.length; i++)
      result.add(options[i]);

    if (getMakeBinary())
      result.add("-D");
    
    if (getUseEqualFrequency())
      result.add("-F");
    
    if (getFindNumBins())
      result.add("-O");
    
    if (getInvertSelection())
      result.add("-V");
    
    result.add("-B");
    result.add("" + getBins());
    
    result.add("-M");
    result.add("" + getDesiredWeightOfInstancesPerInterval());
    
    if (!getAttributeIndices().equals("")) {
      result.add("-R");
      result.add(getAttributeIndices());
    }

    return (String[]) result.toArray(new String[result.size()]);
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
    result.enableAllClasses();
    result.enable(Capability.MISSING_CLASS_VALUES);
    if (!getMakeBinary())
      result.enable(Capability.NO_CLASS);
    
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

    if (m_MakeBinary && m_IgnoreClass) {
      throw new IllegalArgumentException("Can't ignore class when " +
					 "changing the number of attributes!");
    }

    super.setInputFormat(instanceInfo);

    m_DiscretizeCols.setUpper(instanceInfo.numAttributes() - 1);
    m_CutPoints = null;
    
    if (getFindNumBins() && getUseEqualFrequency()) {
      throw new IllegalArgumentException("Bin number optimization in conjunction "+
					 "with equal-frequency binning not implemented.");
    }

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
      + " Discretization is by simple binning. Skips the class"
      + " attribute if set.";
  }
  
  /**
   * Returns the tip text for this property
   *
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String findNumBinsTipText() {

    return "Optimize number of equal-width bins using leave-one-out. Doesn't " +
      "work for equal-frequency binning";
  }

  /**
   * Get the value of FindNumBins.
   *
   * @return Value of FindNumBins.
   */
  public boolean getFindNumBins() {
    
    return m_FindNumBins;
  }
  
  /**
   * Set the value of FindNumBins.
   *
   * @param newFindNumBins Value to assign to FindNumBins.
   */
  public void setFindNumBins(boolean newFindNumBins) {
    
    m_FindNumBins = newFindNumBins;
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
  public String desiredWeightOfInstancesPerIntervalTipText() {

    return "Sets the desired weight of instances per interval for " +
      "equal-frequency binning.";
  }
  
  /**
   * Get the DesiredWeightOfInstancesPerInterval value.
   * @return the DesiredWeightOfInstancesPerInterval value.
   */
  public double getDesiredWeightOfInstancesPerInterval() {

    return m_DesiredWeightOfInstancesPerInterval;
  }

  /**
   * Set the DesiredWeightOfInstancesPerInterval value.
   * @param newDesiredNumber The new DesiredNumber value.
   */
  public void setDesiredWeightOfInstancesPerInterval(double newDesiredNumber) {
    
    m_DesiredWeightOfInstancesPerInterval = newDesiredNumber;
  }
  
  /**
   * Returns the tip text for this property
   *
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String useEqualFrequencyTipText() {

    return "If set to true, equal-frequency binning will be used instead of" +
      " equal-width binning.";
  }
  
  /**
   * Get the value of UseEqualFrequency.
   *
   * @return Value of UseEqualFrequency.
   */
  public boolean getUseEqualFrequency() {
    
    return m_UseEqualFrequency;
  }
  
  /**
   * Set the value of UseEqualFrequency.
   *
   * @param newUseEqualFrequency Value to assign to UseEqualFrequency.
   */
  public void setUseEqualFrequency(boolean newUseEqualFrequency) {
    
    m_UseEqualFrequency = newUseEqualFrequency;
  }

  /**
   * Returns the tip text for this property
   *
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String binsTipText() {

    return "Number of bins.";
  }

  /**
   * Gets the number of bins numeric attributes will be divided into
   *
   * @return the number of bins.
   */
  public int getBins() {

    return m_NumBins;
  }

  /**
   * Sets the number of bins to divide each selected numeric attribute into
   *
   * @param numBins the number of bins
   */
  public void setBins(int numBins) {

    m_NumBins = numBins;
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
   * attribute requested has been discretized into only one interval.)
   */
  public double [] getCutPoints(int attributeIndex) {

    if (m_CutPoints == null) {
      return null;
    }
    return m_CutPoints[attributeIndex];
  }

  /** Generate the cutpoints for each attribute */
  protected void calculateCutPoints() {

    m_CutPoints = new double [getInputFormat().numAttributes()] [];
    for(int i = getInputFormat().numAttributes() - 1; i >= 0; i--) {
      if ((m_DiscretizeCols.isInRange(i)) && 
	  (getInputFormat().attribute(i).isNumeric()) &&
	  (getInputFormat().classIndex() != i)) {
	if (m_FindNumBins) {
	  findNumBins(i);
	} else if (!m_UseEqualFrequency) {
	  calculateCutPointsByEqualWidthBinning(i);
	} else {
	  calculateCutPointsByEqualFrequencyBinning(i);
	}
      }
    }
  }
 
  /**
   * Set cutpoints for a single attribute.
   *
   * @param index the index of the attribute to set cutpoints for
   */
  protected void calculateCutPointsByEqualWidthBinning(int index) {

    // Scan for max and min values
    double max = 0, min = 1, currentVal;
    Instance currentInstance;
    for(int i = 0; i < getInputFormat().numInstances(); i++) {
      currentInstance = getInputFormat().instance(i);
      if (!currentInstance.isMissing(index)) {
	currentVal = currentInstance.value(index);
	if (max < min) {
	  max = min = currentVal;
	}
	if (currentVal > max) {
	  max = currentVal;
	}
	if (currentVal < min) {
	  min = currentVal;
	}
      }
    }
    double binWidth = (max - min) / m_NumBins;
    double [] cutPoints = null;
    if ((m_NumBins > 1) && (binWidth > 0)) {
      cutPoints = new double [m_NumBins - 1];
      for(int i = 1; i < m_NumBins; i++) {
	cutPoints[i - 1] = min + binWidth * i;
      }
    }
    m_CutPoints[index] = cutPoints;
  }
 
  /**
   * Set cutpoints for a single attribute.
   *
   * @param index the index of the attribute to set cutpoints for
   */
  protected void calculateCutPointsByEqualFrequencyBinning(int index) {

    // Copy data so that it can be sorted
    Instances data = new Instances(getInputFormat());

    // Sort input data
    data.sort(index);

    // Compute weight of instances without missing values
    double sumOfWeights = 0;
    for (int i = 0; i < data.numInstances(); i++) {
      if (data.instance(i).isMissing(index)) {
	break;
      } else {
	sumOfWeights += data.instance(i).weight();
      }
    }
    double freq;
    double[] cutPoints = new double[m_NumBins - 1];
    if (getDesiredWeightOfInstancesPerInterval() > 0) {
      freq = getDesiredWeightOfInstancesPerInterval();
      cutPoints = new double[(int)(sumOfWeights / freq)];
    } else {
      freq = sumOfWeights / m_NumBins;
      cutPoints = new double[m_NumBins - 1];
    }

    // Compute break points
    double counter = 0, last = 0;
    int cpindex = 0, lastIndex = -1;
    for (int i = 0; i < data.numInstances() - 1; i++) {

      // Stop if value missing
      if (data.instance(i).isMissing(index)) {
	break;
      }
      counter += data.instance(i).weight();
      sumOfWeights -= data.instance(i).weight();

      // Do we have a potential breakpoint?
      if (data.instance(i).value(index) < 
	  data.instance(i + 1).value(index)) {

	// Have we passed the ideal size?
	if (counter >= freq) {

	  // Is this break point worse than the last one?
	  if (((freq - last) < (counter - freq)) && (lastIndex != -1)) {
	    cutPoints[cpindex] = (data.instance(lastIndex).value(index) +
				  data.instance(lastIndex + 1).value(index)) / 2;
	    counter -= last;
	    last = counter;
	    lastIndex = i;
	  } else {
	    cutPoints[cpindex] = (data.instance(i).value(index) +
				  data.instance(i + 1).value(index)) / 2;
	    counter = 0;
	    last = 0;
	    lastIndex = -1;
	  }
	  cpindex++;
	  freq = (sumOfWeights + counter) / ((cutPoints.length + 1) - cpindex);
	} else {
	  lastIndex = i;
	  last = counter;
	}
      }
    }

    // Check whether there was another possibility for a cut point
    if ((cpindex < cutPoints.length) && (lastIndex != -1)) {
      cutPoints[cpindex] = (data.instance(lastIndex).value(index) +
			    data.instance(lastIndex + 1).value(index)) / 2;      
      cpindex++;
    }

    // Did we find any cutpoints?
    if (cpindex == 0) {
      m_CutPoints[index] = null;
    } else {
      double[] cp = new double[cpindex];
      for (int i = 0; i < cpindex; i++) {
	cp[i] = cutPoints[i];
      }
      m_CutPoints[index] = cp;
    }
  }

  /**
   * Optimizes the number of bins using leave-one-out cross-validation.
   *
   * @param index the attribute index
   */
  protected void findNumBins(int index) {

    double min = Double.MAX_VALUE, max = -Double.MAX_VALUE, binWidth = 0, 
      entropy, bestEntropy = Double.MAX_VALUE, currentVal;
    double[] distribution;
    int bestNumBins  = 1;
    Instance currentInstance;

    // Find minimum and maximum
    for (int i = 0; i < getInputFormat().numInstances(); i++) {
      currentInstance = getInputFormat().instance(i);
      if (!currentInstance.isMissing(index)) {
	currentVal = currentInstance.value(index);
	if (currentVal > max) {
	  max = currentVal;
	}
	if (currentVal < min) {
	  min = currentVal;
	}
      }
    }

    // Find best number of bins
    for (int i = 0; i < m_NumBins; i++) {
      distribution = new double[i + 1];
      binWidth = (max - min) / (i + 1);

      // Compute distribution
      for (int j = 0; j < getInputFormat().numInstances(); j++) {
	currentInstance = getInputFormat().instance(j);
	if (!currentInstance.isMissing(index)) {
	  for (int k = 0; k < i + 1; k++) {
	    if (currentInstance.value(index) <= 
		(min + (((double)k + 1) * binWidth))) {
	      distribution[k] += currentInstance.weight();
	      break;
	    }
	  }
	}
      }

      // Compute cross-validated entropy
      entropy = 0;
      for (int k = 0; k < i + 1; k++) {
	if (distribution[k] < 2) {
	  entropy = Double.MAX_VALUE;
	  break;
	}
	entropy -= distribution[k] * Math.log((distribution[k] - 1) / 
					      binWidth);
      }

      // Best entropy so far?
      if (entropy < bestEntropy) {
	bestEntropy = entropy;
	bestNumBins = i + 1;
      }
    }

    // Compute cut points
    double [] cutPoints = null;
    if ((bestNumBins > 1) && (binWidth > 0)) {
      cutPoints = new double [bestNumBins - 1];
      for(int i = 1; i < bestNumBins; i++) {
	cutPoints[i - 1] = min + binWidth * i;
      }
    }
    m_CutPoints[index] = cutPoints;
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
	  && (getInputFormat().attribute(i).isNumeric())
	  && (getInputFormat().classIndex() != i)) {
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
	  getInputFormat().attribute(i).isNumeric() &&
	  (getInputFormat().classIndex() != i)) {
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
    } 
    //RANKING BEGIN
    if(instance instanceof PreferenceDenseInstance){
    	PreferenceDenseInstance pdi = (PreferenceDenseInstance) instance;
    	inst = new PreferenceDenseInstance(instance.weight(), vals, pdi.getHashMap());
    }
    //RANKING END
    else {
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
    return RevisionUtils.extract("$Revision: 6567 $");
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
