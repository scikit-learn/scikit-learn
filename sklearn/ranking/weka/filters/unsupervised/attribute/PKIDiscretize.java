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
 *    PKIDiscretize.java
 *    Copyright (C) 2003 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.filters.unsupervised.attribute;

import weka.core.Instances;
import weka.core.Option;
import weka.core.RevisionUtils;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformationHandler;
import weka.core.Utils;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;

import java.util.Enumeration;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * Discretizes numeric attributes using equal frequency binning, where the number of bins is equal to the square root of the number of non-missing values.<br/>
 * <br/>
 * For more information, see:<br/>
 * <br/>
 * Ying Yang, Geoffrey I. Webb: Proportional k-Interval Discretization for Naive-Bayes Classifiers. In: 12th European Conference on Machine Learning, 564-575, 2001.
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;inproceedings{Yang2001,
 *    author = {Ying Yang and Geoffrey I. Webb},
 *    booktitle = {12th European Conference on Machine Learning},
 *    pages = {564-575},
 *    publisher = {Springer},
 *    series = {LNCS},
 *    title = {Proportional k-Interval Discretization for Naive-Bayes Classifiers},
 *    volume = {2167},
 *    year = {2001}
 * }
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -unset-class-temporarily
 *  Unsets the class index temporarily before the filter is
 *  applied to the data.
 *  (default: no)</pre>
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
 * @author Richard Kirkby (rkirkby@cs.waikato.ac.nz)
 * @version $Revision: 1.9 $
 */
public class PKIDiscretize 
  extends Discretize
  implements TechnicalInformationHandler {
  
  /** for serialization */
  static final long serialVersionUID = 6153101248977702675L;

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

    // alter child behaviour to do what we want
    m_FindNumBins = true;
    return super.setInputFormat(instanceInfo);
  }

  /**
   * Finds the number of bins to use and creates the cut points.
   *
   * @param index the attribute index
   */
  protected void findNumBins(int index) {

    Instances toFilter = getInputFormat();

    // Find number of instances for attribute where not missing
    int numOfInstances = toFilter.numInstances();
    for (int i = 0; i < toFilter.numInstances(); i++) {
      if (toFilter.instance(i).isMissing(index))
	numOfInstances--;
    }

    m_NumBins = (int)(Math.sqrt(numOfInstances));

    if (m_NumBins > 0) {
      calculateCutPointsByEqualFrequencyBinning(index);
    }
  }

  /**
   * Gets an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector result = new Vector();
    
    result.addElement(new Option(
	"\tUnsets the class index temporarily before the filter is\n"
	+ "\tapplied to the data.\n"
	+ "\t(default: no)",
	"unset-class-temporarily", 1, "-unset-class-temporarily"));
    
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

    setIgnoreClass(Utils.getFlag("unset-class-temporarily", options));
    setMakeBinary(Utils.getFlag('D', options));
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
  public String[] getOptions() {
    Vector        result;

    result = new Vector();

    if (getMakeBinary())
      result.add("-D");
    
    if (getInvertSelection())
      result.add("-V");
    
    if (!getAttributeIndices().equals("")) {
      result.add("-R");
      result.add(getAttributeIndices());
    }

    return (String[]) result.toArray(new String[result.size()]);
  }

  /**
   * Returns a string describing this filter
   *
   * @return a description of the filter suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {

    return "Discretizes numeric attributes using equal frequency binning,"
      + " where the number of bins is equal to the square root of the"
      + " number of non-missing values.\n\n"
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
    
    result = new TechnicalInformation(Type.INPROCEEDINGS);
    result.setValue(Field.AUTHOR, "Ying Yang and Geoffrey I. Webb");
    result.setValue(Field.TITLE, "Proportional k-Interval Discretization for Naive-Bayes Classifiers");
    result.setValue(Field.BOOKTITLE, "12th European Conference on Machine Learning");
    result.setValue(Field.YEAR, "2001");
    result.setValue(Field.PAGES, "564-575");
    result.setValue(Field.PUBLISHER, "Springer");
    result.setValue(Field.SERIES, "LNCS");
    result.setValue(Field.VOLUME, "2167");
    
    return result;
  }
  
  /**
   * Returns the tip text for this property
   *
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String findNumBinsTipText() {

    return "Ignored.";
  }

  /**
   * Get the value of FindNumBins.
   *
   * @return Value of FindNumBins.
   */
  public boolean getFindNumBins() {
    
    return false;
  }
  
  /**
   * Set the value of FindNumBins.
   *
   * @param newFindNumBins Value to assign to FindNumBins.
   */
  public void setFindNumBins(boolean newFindNumBins) {
    
  }
  
  /**
   * Returns the tip text for this property
   *
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String useEqualFrequencyTipText() {

    return "Always true.";
  }

  /**
   * Get the value of UseEqualFrequency.
   *
   * @return Value of UseEqualFrequency.
   */
  public boolean getUseEqualFrequency() {
    
    return true;
  }
  
  /**
   * Set the value of UseEqualFrequency.
   *
   * @param newUseEqualFrequency Value to assign to UseEqualFrequency.
   */
  public void setUseEqualFrequency(boolean newUseEqualFrequency) {
    
  }

  /**
   * Returns the tip text for this property
   *
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String binsTipText() {

    return "Ignored.";
  }

  /**
   * Ignored
   *
   * @return the number of bins.
   */
  public int getBins() {

    return 0;
  }

  /**
   * Ignored
   *
   * @param numBins the number of bins
   */
  public void setBins(int numBins) {

  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 1.9 $");
  }

  /**
   * Main method for testing this class.
   *
   * @param argv should contain arguments to the filter: use -h for help
   */
  public static void main(String [] argv) {
    runFilter(new PKIDiscretize(), argv);
  }
}
