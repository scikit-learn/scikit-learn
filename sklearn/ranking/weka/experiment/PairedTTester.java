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
 *    PairedTTester.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.experiment;

import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.Range;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.core.Utils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.Serializable;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Enumeration;
import java.util.Vector;

/**
 * Calculates T-Test statistics on data stored in a set of instances. <p/>
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -D &lt;index,index2-index4,...&gt;
 *  Specify list of columns that specify a unique
 *  dataset.
 *  First and last are valid indexes. (default none)</pre>
 * 
 * <pre> -R &lt;index&gt;
 *  Set the index of the column containing the run number</pre>
 * 
 * <pre> -F &lt;index&gt;
 *  Set the index of the column containing the fold number</pre>
 * 
 * <pre> -G &lt;index1,index2-index4,...&gt;
 *  Specify list of columns that specify a unique
 *  'result generator' (eg: classifier name and options).
 *  First and last are valid indexes. (default none)</pre>
 * 
 * <pre> -S &lt;significance level&gt;
 *  Set the significance level for comparisons (default 0.05)</pre>
 * 
 * <pre> -V
 *  Show standard deviations</pre>
 * 
 * <pre> -L
 *  Produce table comparisons in Latex table format</pre>
 * 
 * <pre> -csv
 *  Produce table comparisons in CSV table format</pre>
 * 
 * <pre> -html
 *  Produce table comparisons in HTML table format</pre>
 * 
 * <pre> -significance
 *  Produce table comparisons with only the significance values</pre>
 * 
 * <pre> -gnuplot
 *  Produce table comparisons output suitable for GNUPlot</pre>
 * 
 <!-- options-end -->
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 6432 $
 */
public class PairedTTester 
  implements OptionHandler, Tester, RevisionHandler {
  
  /** for serialization */
  static final long serialVersionUID = 8370014624008728610L;

  /** The set of instances we will analyse */
  protected Instances m_Instances;

  /** The index of the column containing the run number */
  protected int m_RunColumn = 0;

  /** The option setting for the run number column (-1 means last) */
  protected int m_RunColumnSet = -1;

  /** The option setting for the fold number column (-1 means none) */
  protected int m_FoldColumn = -1;

  /** The column to sort on (-1 means default sorting) */
  protected int m_SortColumn = -1;

  /** The sorting of the datasets (according to the sort column) */
  protected int[] m_SortOrder = null;

  /** The sorting of the columns (test base is always first) */
  protected int[] m_ColOrder = null;

  /** The significance level for comparisons */
  protected double m_SignificanceLevel = 0.05;

  /**
   * The range of columns that specify a unique "dataset"
   * (eg: scheme plus configuration)
   */
  protected Range m_DatasetKeyColumnsRange = new Range();

  /** An array containing the indexes of just the selected columns */ 
  protected int [] m_DatasetKeyColumns;

  /** The list of dataset specifiers */
  protected DatasetSpecifiers m_DatasetSpecifiers = 
    new DatasetSpecifiers();

  /**
   * The range of columns that specify a unique result set
   * (eg: scheme plus configuration)
   */
  protected Range m_ResultsetKeyColumnsRange = new Range();

  /** An array containing the indexes of just the selected columns */ 
  protected int [] m_ResultsetKeyColumns;

  /** An array containing the indexes of the datasets to display */
  protected int[] m_DisplayedResultsets = null;

  /** Stores a vector for each resultset holding all instances in each set */
  protected FastVector m_Resultsets = new FastVector();

  /** Indicates whether the instances have been partitioned */
  protected boolean m_ResultsetsValid;

  /** Indicates whether standard deviations should be displayed */
  protected boolean m_ShowStdDevs = false;
  
  /** the instance of the class to produce the output. */
  protected ResultMatrix m_ResultMatrix = new ResultMatrixPlainText();
  
  /** A list of unique "dataset" specifiers that have been observed */
  protected class DatasetSpecifiers
    implements RevisionHandler, Serializable {

    /** for serialization. */
    private static final long serialVersionUID = -9020938059902723401L;
    
    /** the specifiers that have been observed */
    FastVector m_Specifiers = new FastVector();

    /**
     * Removes all specifiers.
     */
    protected void removeAllSpecifiers() {

      m_Specifiers.removeAllElements();
    }

    /** 
     * Add an instance to the list of specifiers (if necessary)
     * 
     * @param inst	the instance to add
     */
    protected void add(Instance inst) {
      
      for (int i = 0; i < m_Specifiers.size(); i++) {
	Instance specifier = (Instance)m_Specifiers.elementAt(i);
	boolean found = true;
	for (int j = 0; j < m_DatasetKeyColumns.length; j++) {
	  if (inst.value(m_DatasetKeyColumns[j]) !=
	      specifier.value(m_DatasetKeyColumns[j])) {
	    found = false;
	  }
	}
	if (found) {
	  return;
	}
      }
      m_Specifiers.addElement(inst);
    }

    /**
     * Get the template at the given position.
     * 
     * @param i		the index
     * @return		the template
     */
    protected Instance specifier(int i) {

      return (Instance)m_Specifiers.elementAt(i);
    }

    /**
     * Gets the number of specifiers.
     * 
     * @return		the current number of specifiers
     */
    protected int numSpecifiers() {

      return m_Specifiers.size();
    }
    
    /**
     * Returns the revision string.
     * 
     * @return		the revision
     */
    public String getRevision() {
      return RevisionUtils.extract("$Revision: 6432 $");
    }
  }

  /** Utility class to store the instances pertaining to a dataset */
  protected class Dataset
    implements RevisionHandler, Serializable {

    /** for serialization. */
    private static final long serialVersionUID = -2801397601839433282L;

    /** the template */
    Instance m_Template;

    /** the dataset */
    FastVector m_Dataset;

    /**
     * Constructor
     * 
     * @param template	the template
     */
    public Dataset(Instance template) {

      m_Template = template;
      m_Dataset = new FastVector();
      add(template);
    }
    
    /**
     * Returns true if the two instances match on those attributes that have
     * been designated key columns (eg: scheme name and scheme options)
     *
     * @param first the first instance
     * @return true if first and second match on the currently set key columns
     */
    protected boolean matchesTemplate(Instance first) {
      
      for (int i = 0; i < m_DatasetKeyColumns.length; i++) {
	if (first.value(m_DatasetKeyColumns[i]) !=
	    m_Template.value(m_DatasetKeyColumns[i])) {
	  return false;
	}
      }
      return true;
    }

    /**
     * Adds the given instance to the dataset
     * 
     * @param inst	the instance to add
     */
    protected void add(Instance inst) {
      
      m_Dataset.addElement(inst);
    }

    /**
     * Returns a vector containing the instances in the dataset
     * 
     * @return 		the current contents
     */
    protected FastVector contents() {

      return m_Dataset;
    }

    /**
     * Sorts the instances in the dataset by the run number.
     *
     * @param runColumn a value of type 'int'
     */
    public void sort(int runColumn) {

      double [] runNums = new double [m_Dataset.size()];
      for (int j = 0; j < runNums.length; j++) {
	runNums[j] = ((Instance) m_Dataset.elementAt(j)).value(runColumn);
      }
      int [] index = Utils.stableSort(runNums);
      FastVector newDataset = new FastVector(runNums.length);
      for (int j = 0; j < index.length; j++) {
	newDataset.addElement(m_Dataset.elementAt(index[j]));
      }
      m_Dataset = newDataset;
    }
    
    /**
     * Returns the revision string.
     * 
     * @return		the revision
     */
    public String getRevision() {
      return RevisionUtils.extract("$Revision: 6432 $");
    }
  }
 
  /** Utility class to store the instances in a resultset */
  protected class Resultset
    implements RevisionHandler, Serializable {

    /** for serialization. */
    private static final long serialVersionUID = 1543786683821339978L;

    /** the template */
    Instance m_Template;
    
    /** the dataset */
    FastVector m_Datasets;

    /**
     * Constructir
     * 
     * @param template		the template
     */
    public Resultset(Instance template) {

      m_Template = template;
      m_Datasets = new FastVector();
      add(template);
    }
    
    /**
     * Returns true if the two instances match on those attributes that have
     * been designated key columns (eg: scheme name and scheme options)
     *
     * @param first the first instance
     * @return true if first and second match on the currently set key columns
     */
    protected boolean matchesTemplate(Instance first) {
      
      for (int i = 0; i < m_ResultsetKeyColumns.length; i++) {
	if (first.value(m_ResultsetKeyColumns[i]) !=
	    m_Template.value(m_ResultsetKeyColumns[i])) {
	  return false;
	}
      }
      return true;
    }

    /**
     * Returns a string descriptive of the resultset key column values
     * for this resultset
     *
     * @return a value of type 'String'
     */
    protected String templateString() {

      String result = "";
      String tempResult = "";
      for (int i = 0; i < m_ResultsetKeyColumns.length; i++) {
	tempResult = m_Template.toString(m_ResultsetKeyColumns[i]) + ' ';

	// compact the string
        tempResult = Utils.removeSubstring(tempResult, "weka.classifiers.");
        tempResult = Utils.removeSubstring(tempResult, "weka.filters.");
        tempResult = Utils.removeSubstring(tempResult, "weka.attributeSelection.");
	result += tempResult;
      }
      return result.trim();
    }
    
    /**
     * Returns a vector containing all instances belonging to one dataset.
     *
     * @param inst a template instance
     * @return a value of type 'FastVector'
     */
    public FastVector dataset(Instance inst) {

      for (int i = 0; i < m_Datasets.size(); i++) {
	if (((Dataset)m_Datasets.elementAt(i)).matchesTemplate(inst)) {
	  return ((Dataset)m_Datasets.elementAt(i)).contents();
	} 
      }
      return null;
    }
    
    /**
     * Adds an instance to this resultset
     *
     * @param newInst a value of type 'Instance'
     */
    public void add(Instance newInst) {
      
      for (int i = 0; i < m_Datasets.size(); i++) {
	if (((Dataset)m_Datasets.elementAt(i)).matchesTemplate(newInst)) {
	  ((Dataset)m_Datasets.elementAt(i)).add(newInst);
	  return;
	}
      }
      Dataset newDataset = new Dataset(newInst);
      m_Datasets.addElement(newDataset);
    }

    /**
     * Sorts the instances in each dataset by the run number.
     *
     * @param runColumn a value of type 'int'
     */
    public void sort(int runColumn) {

      for (int i = 0; i < m_Datasets.size(); i++) {
	((Dataset)m_Datasets.elementAt(i)).sort(runColumn);
      }
    }
    
    /**
     * Returns the revision string.
     * 
     * @return		the revision
     */
    public String getRevision() {
      return RevisionUtils.extract("$Revision: 6432 $");
    }
  } // Resultset


  /**
   * Returns a string descriptive of the key column values for
   * the "datasets
   *
   * @param template the template
   * @return a value of type 'String'
   */
  protected String templateString(Instance template) {
    
    String result = "";
    for (int i = 0; i < m_DatasetKeyColumns.length; i++) {
      result += template.toString(m_DatasetKeyColumns[i]) + ' ';
    }
    if (result.startsWith("weka.classifiers.")) {
      result = result.substring("weka.classifiers.".length());
    }
    return result.trim();
  }

  /**
   * Sets the matrix to use to produce the output.
   * @param matrix the instance to use to produce the output
   * @see ResultMatrix
   */
  public void setResultMatrix(ResultMatrix matrix) {
    m_ResultMatrix = matrix;
  }

  /**
   * Gets the instance that produces the output.
   * @return the instance to produce the output
   */
  public ResultMatrix getResultMatrix() {
    return m_ResultMatrix;
  }

  /**
   * Set whether standard deviations are displayed or not.
   * @param s true if standard deviations are to be displayed
   */
  public void setShowStdDevs(boolean s) {
    m_ShowStdDevs = s;
  }

  /**
   * Returns true if standard deviations have been requested.
   * @return true if standard deviations are to be displayed.
   */
  public boolean getShowStdDevs() {
    return m_ShowStdDevs;
  }
  
  /**
   * Separates the instances into resultsets and by dataset/run.
   *
   * @throws Exception if the TTest parameters have not been set.
   */
  protected void prepareData() throws Exception {

    if (m_Instances == null) {
      throw new Exception("No instances have been set");
    }
    if (m_RunColumnSet == -1) {
      m_RunColumn = m_Instances.numAttributes() - 1;
    } else {
      m_RunColumn = m_RunColumnSet;
    }

    if (m_ResultsetKeyColumnsRange == null) {
      throw new Exception("No result specifier columns have been set");
    }
    m_ResultsetKeyColumnsRange.setUpper(m_Instances.numAttributes() - 1);
    m_ResultsetKeyColumns = m_ResultsetKeyColumnsRange.getSelection();

    if (m_DatasetKeyColumnsRange == null) {
      throw new Exception("No dataset specifier columns have been set");
    }
    m_DatasetKeyColumnsRange.setUpper(m_Instances.numAttributes() - 1);
    m_DatasetKeyColumns = m_DatasetKeyColumnsRange.getSelection();
    
    //  Split the data up into result sets
    m_Resultsets.removeAllElements();  
    m_DatasetSpecifiers.removeAllSpecifiers();
    for (int i = 0; i < m_Instances.numInstances(); i++) {
      Instance current = m_Instances.instance(i);
      if (current.isMissing(m_RunColumn)) {
	throw new Exception("Instance has missing value in run "
			    + "column!\n" + current);
      } 
      for (int j = 0; j < m_ResultsetKeyColumns.length; j++) {
	if (current.isMissing(m_ResultsetKeyColumns[j])) {
	  throw new Exception("Instance has missing value in resultset key "
			      + "column " + (m_ResultsetKeyColumns[j] + 1)
			      + "!\n" + current);
	}
      }
      for (int j = 0; j < m_DatasetKeyColumns.length; j++) {
	if (current.isMissing(m_DatasetKeyColumns[j])) {
	  throw new Exception("Instance has missing value in dataset key "
			      + "column " + (m_DatasetKeyColumns[j] + 1)
			      + "!\n" + current);
	}
      }
      boolean found = false;
      for (int j = 0; j < m_Resultsets.size(); j++) {
	Resultset resultset = (Resultset) m_Resultsets.elementAt(j);
	if (resultset.matchesTemplate(current)) {
	  resultset.add(current);
	  found = true;
	  break;
	}
      }
      if (!found) {
	Resultset resultset = new Resultset(current);
	m_Resultsets.addElement(resultset);
      }

      m_DatasetSpecifiers.add(current);
    }

    // Tell each resultset to sort on the run column
    for (int j = 0; j < m_Resultsets.size(); j++) {
      Resultset resultset = (Resultset) m_Resultsets.elementAt(j);
      if (m_FoldColumn >= 0) {
        // sort on folds first in case they are out of order
        resultset.sort(m_FoldColumn);
      }
      resultset.sort(m_RunColumn);
    }

    m_ResultsetsValid = true;
  }

  /**
   * Gets the number of datasets in the resultsets
   *
   * @return the number of datasets in the resultsets
   */
  public int getNumDatasets() {

    if (!m_ResultsetsValid) {
      try {
	prepareData();
      } catch (Exception ex) {
	ex.printStackTrace();
	return 0;
      }
    }
    return m_DatasetSpecifiers.numSpecifiers();
  }

  /**
   * Gets the number of resultsets in the data.
   *
   * @return the number of resultsets in the data
   */
  public int getNumResultsets() {

    if (!m_ResultsetsValid) {
      try {
  prepareData();
      } catch (Exception ex) {
  ex.printStackTrace();
  return 0;
      }
    }
    return m_Resultsets.size();
  }

  /**
   * Gets a string descriptive of the specified resultset.
   *
   * @param index the index of the resultset
   * @return a descriptive string for the resultset
   */
  public String getResultsetName(int index) {

    if (!m_ResultsetsValid) {
      try {
	prepareData();
      } catch (Exception ex) {
	ex.printStackTrace();
	return null;
      }
    }
    return ((Resultset) m_Resultsets.elementAt(index)).templateString();
  }
  
  /**
   * Checks whether the resultset with the given index shall be displayed.
   * 
   * @param index the index of the resultset to check whether it shall be displayed 
   * @return whether the specified resultset is displayed 
   */
  public boolean displayResultset(int index) {
    boolean       result;
    int           i;
    
    result = true;

    if (m_DisplayedResultsets != null) {
      result = false;
      for (i = 0; i < m_DisplayedResultsets.length; i++) {
        if (m_DisplayedResultsets[i] == index) {
          result = true;
          break;
        }
      }
    }
      
    return result;
  }
  
  /**
   * Computes a paired t-test comparison for a specified dataset between
   * two resultsets.
   *
   * @param datasetSpecifier the dataset specifier
   * @param resultset1Index the index of the first resultset
   * @param resultset2Index the index of the second resultset
   * @param comparisonColumn the column containing values to compare
   * @return the results of the paired comparison
   * @throws Exception if an error occurs
   */
  public PairedStats calculateStatistics(Instance datasetSpecifier,
					 int resultset1Index,
					 int resultset2Index,
					 int comparisonColumn) throws Exception {

    if (m_Instances.attribute(comparisonColumn).type()
	!= Attribute.NUMERIC) {
      throw new Exception("Comparison column " + (comparisonColumn + 1)
			  + " ("
			  + m_Instances.attribute(comparisonColumn).name()
			  + ") is not numeric");
    }
    if (!m_ResultsetsValid) {
      prepareData();
    }

    Resultset resultset1 = (Resultset) m_Resultsets.elementAt(resultset1Index);
    Resultset resultset2 = (Resultset) m_Resultsets.elementAt(resultset2Index);
    FastVector dataset1 = resultset1.dataset(datasetSpecifier);
    FastVector dataset2 = resultset2.dataset(datasetSpecifier);
    String datasetName = templateString(datasetSpecifier);
    if (dataset1 == null) {
      throw new Exception("No results for dataset=" + datasetName
			 + " for resultset=" + resultset1.templateString());
    } else if (dataset2 == null) {
      throw new Exception("No results for dataset=" + datasetName
			 + " for resultset=" + resultset2.templateString());
    } else if (dataset1.size() != dataset2.size()) {
      throw new Exception("Results for dataset=" + datasetName
			  + " differ in size for resultset="
			  + resultset1.templateString()
			  + " and resultset="
			  + resultset2.templateString()
			  );
    }
    
    PairedStats pairedStats = new PairedStats(m_SignificanceLevel);

    for (int k = 0; k < dataset1.size(); k ++) {
      Instance current1 = (Instance) dataset1.elementAt(k);
      Instance current2 = (Instance) dataset2.elementAt(k);
      if (current1.isMissing(comparisonColumn)) {
	System.err.println("Instance has missing value in comparison "
			   + "column!\n" + current1);
	continue;
      }
      if (current2.isMissing(comparisonColumn)) {
	System.err.println("Instance has missing value in comparison "
			   + "column!\n" + current2);
	continue;
      }
      if (current1.value(m_RunColumn) != current2.value(m_RunColumn)) {
	System.err.println("Run numbers do not match!\n"
			    + current1 + current2);
      }
      if (m_FoldColumn != -1) {
	if (current1.value(m_FoldColumn) != current2.value(m_FoldColumn)) {
	  System.err.println("Fold numbers do not match!\n"
			     + current1 + current2);
	}
      }
      double value1 = current1.value(comparisonColumn);
      double value2 = current2.value(comparisonColumn);
      pairedStats.add(value1, value2);
    }
    pairedStats.calculateDerived();
    //System.err.println("Differences stats:\n" + pairedStats.differencesStats);
    return pairedStats;

  }
  
  /**
   * Creates a key that maps resultset numbers to their descriptions.
   *
   * @return a value of type 'String'
   */
  public String resultsetKey() {

    if (!m_ResultsetsValid) {
      try {
	prepareData();
      } catch (Exception ex) {
	ex.printStackTrace();
	return ex.getMessage();
      }
    }
    String result = "";
    for (int j = 0; j < getNumResultsets(); j++) {
      result += "(" + (j + 1) + ") " + getResultsetName(j) + '\n';
    }
    return result + '\n';
  }
  
  /**
   * Creates a "header" string describing the current resultsets.
   *
   * @param comparisonColumn a value of type 'int'
   * @return a value of type 'String'
   */
  public String header(int comparisonColumn) {

    if (!m_ResultsetsValid) {
      try {
	prepareData();
      } catch (Exception ex) {
	ex.printStackTrace();
	return ex.getMessage();
      }
    }
    
    initResultMatrix();
    m_ResultMatrix.addHeader("Tester", getClass().getName());
    m_ResultMatrix.addHeader("Analysing", m_Instances.attribute(comparisonColumn).name());
    m_ResultMatrix.addHeader("Datasets", Integer.toString(getNumDatasets()));
    m_ResultMatrix.addHeader("Resultsets", Integer.toString(getNumResultsets()));
    m_ResultMatrix.addHeader("Confidence", getSignificanceLevel() + " (two tailed)");
    m_ResultMatrix.addHeader("Sorted by", getSortColumnName());
    m_ResultMatrix.addHeader("Date", (new SimpleDateFormat()).format(new Date()));

    return m_ResultMatrix.toStringHeader() + "\n";
  }

  /**
   * Carries out a comparison between all resultsets, counting the number
   * of datsets where one resultset outperforms the other.
   *
   * @param comparisonColumn the index of the comparison column
   * @param nonSigWin for storing the non-significant wins
   * @return a 2d array where element [i][j] is the number of times resultset
   * j performed significantly better than resultset i.
   * @throws Exception if an error occurs
   */
  public int [][] multiResultsetWins(int comparisonColumn, int [][] nonSigWin)
    throws Exception {

    int numResultsets = getNumResultsets();
    int [][] win = new int [numResultsets][numResultsets];
    //    int [][] nonSigWin = new int [numResultsets][numResultsets];
    for (int i = 0; i < numResultsets; i++) {
      for (int j = i + 1; j < numResultsets; j++) {
	System.err.print("Comparing (" + (i + 1) + ") with ("
			 + (j + 1) + ")\r");
	System.err.flush();
	for (int k = 0; k < getNumDatasets(); k++) {
	  try {
	    PairedStats pairedStats = 
	      calculateStatistics(m_DatasetSpecifiers.specifier(k), i, j,
				  comparisonColumn);
	    if (pairedStats.differencesSignificance < 0) {
	      win[i][j]++;
	    } else if (pairedStats.differencesSignificance > 0) {
	      win[j][i]++;
	    }

	    if (pairedStats.differencesStats.mean < 0) {
	      nonSigWin[i][j]++;
	    } else if (pairedStats.differencesStats.mean > 0) {
	      nonSigWin[j][i]++;
	    }
	  } catch (Exception ex) {
	    //ex.printStackTrace();
	    System.err.println(ex.getMessage());
	  }
	}
      }
    }
    return win;
  }

  /**
   * clears the content and fills the column and row names according to the
   * given sorting
   */
  protected void initResultMatrix() {
    m_ResultMatrix.setSize(getNumResultsets(), getNumDatasets());
    m_ResultMatrix.setShowStdDev(m_ShowStdDevs);

    for (int i = 0; i < getNumDatasets(); i++)
      m_ResultMatrix.setRowName(i, 
          templateString(m_DatasetSpecifiers.specifier(i)));

    for (int j = 0; j < getNumResultsets(); j++) {
      m_ResultMatrix.setColName(j, getResultsetName(j));
      m_ResultMatrix.setColHidden(j, !displayResultset(j));
    }
  }
  
  /**
   * Carries out a comparison between all resultsets, counting the number
   * of datsets where one resultset outperforms the other. The results
   * are summarized in a table.
   *
   * @param comparisonColumn the index of the comparison column
   * @return the results in a string
   * @throws Exception if an error occurs
   */
  public String multiResultsetSummary(int comparisonColumn)
    throws Exception {
    
    int[][] nonSigWin = new int [getNumResultsets()][getNumResultsets()];
    int[][] win = multiResultsetWins(comparisonColumn, nonSigWin);
    
    initResultMatrix();    
    m_ResultMatrix.setSummary(nonSigWin, win);
    
    return m_ResultMatrix.toStringSummary();
  }

  /**
   * returns a ranking of the resultsets
   * 
   * @param comparisonColumn	the column to compare with
   * @return			the ranking
   * @throws Exception		if something goes wrong
   */
  public String multiResultsetRanking(int comparisonColumn)
    throws Exception {
    
    int[][] nonSigWin = new int [getNumResultsets()][getNumResultsets()];
    int[][] win       = multiResultsetWins(comparisonColumn, nonSigWin);
    
    initResultMatrix();    
    m_ResultMatrix.setRanking(win);

    return m_ResultMatrix.toStringRanking();
  }
				    
  /**
   * Creates a comparison table where a base resultset is compared to the
   * other resultsets. Results are presented for every dataset.
   *
   * @param baseResultset the index of the base resultset
   * @param comparisonColumn the index of the column to compare over
   * @return the comparison table string
   * @throws Exception if an error occurs
   */
  public String multiResultsetFull(int baseResultset,
				   int comparisonColumn) throws Exception {

    int maxWidthMean = 2;
    int maxWidthStdDev = 2;
    
    double[] sortValues = new double[getNumDatasets()];
      
    // determine max field width
    for (int i = 0; i < getNumDatasets(); i++) {
      sortValues[i] = Double.POSITIVE_INFINITY;  // sorts skipped cols to end
      
      for (int j = 0; j < getNumResultsets(); j++) {
        if (!displayResultset(j))
          continue;
	try {
	  PairedStats pairedStats = 
	    calculateStatistics(m_DatasetSpecifiers.specifier(i), 
				baseResultset, j, comparisonColumn);
          if (!Double.isInfinite(pairedStats.yStats.mean) &&
              !Double.isNaN(pairedStats.yStats.mean)) {
            double width = ((Math.log(Math.abs(pairedStats.yStats.mean)) / 
                             Math.log(10))+1);
            if (width > maxWidthMean) {
              maxWidthMean = (int)width;
            }
          }

          if (j == baseResultset) {
            if (getSortColumn() != -1)
              sortValues[i] = calculateStatistics(
                                m_DatasetSpecifiers.specifier(i), 
                                baseResultset, j, getSortColumn()).xStats.mean;
            else
              sortValues[i] = i;
          }
	  
	  if (m_ShowStdDevs &&
              !Double.isInfinite(pairedStats.yStats.stdDev) &&
              !Double.isNaN(pairedStats.yStats.stdDev)) {
	    double width = ((Math.log(Math.abs(pairedStats.yStats.stdDev)) / 
                             Math.log(10))+1);
	    if (width > maxWidthStdDev) {
	      maxWidthStdDev = (int)width;
	    }
	  }
	}  catch (Exception ex) {
	  //ex.printStackTrace();
          System.err.println(ex);
	}
      }
    }

    // sort rows according to sort column
    m_SortOrder = Utils.sort(sortValues);

    // determine column order
    m_ColOrder = new int[getNumResultsets()];
    m_ColOrder[0] = baseResultset;
    int index = 1;
    for (int i = 0; i < getNumResultsets(); i++) {
      if (i == baseResultset)
        continue;
      m_ColOrder[index] = i;
      index++;
    }

    // setup matrix
    initResultMatrix();    
    m_ResultMatrix.setRowOrder(m_SortOrder);
    m_ResultMatrix.setColOrder(m_ColOrder);
    m_ResultMatrix.setMeanWidth(maxWidthMean);
    m_ResultMatrix.setStdDevWidth(maxWidthStdDev);
    m_ResultMatrix.setSignificanceWidth(1);

    // make sure that test base is displayed, even though it might not be
    // selected
    for (int i = 0; i < m_ResultMatrix.getColCount(); i++) {
      if (    (i == baseResultset)
           && (m_ResultMatrix.getColHidden(i)) ) {
        m_ResultMatrix.setColHidden(i, false);
        System.err.println("Note: test base was hidden - set visible!");
      }
    }
    
    // the data
    for (int i = 0; i < getNumDatasets(); i++) {
      m_ResultMatrix.setRowName(i, 
          templateString(m_DatasetSpecifiers.specifier(i)));

      for (int j = 0; j < getNumResultsets(); j++) {
        try {
          // calc stats
          PairedStats pairedStats = 
            calculateStatistics(m_DatasetSpecifiers.specifier(i), 
                baseResultset, j, comparisonColumn);

          // count
          m_ResultMatrix.setCount(i, pairedStats.count);

          // mean
          m_ResultMatrix.setMean(j, i, pairedStats.yStats.mean);
          
          // std dev
          m_ResultMatrix.setStdDev(j, i, pairedStats.yStats.stdDev);

          // significance
          if (pairedStats.differencesSignificance < 0)
            m_ResultMatrix.setSignificance(j, i, ResultMatrix.SIGNIFICANCE_WIN);
          else if (pairedStats.differencesSignificance > 0)
            m_ResultMatrix.setSignificance(j, i, ResultMatrix.SIGNIFICANCE_LOSS);
          else
            m_ResultMatrix.setSignificance(j, i, ResultMatrix.SIGNIFICANCE_TIE);
        }
        catch (Exception e) {
          //e.printStackTrace();
          System.err.println(e);
        }
      }
    }

    // generate output
    StringBuffer result = new StringBuffer(1000);
    try {
      result.append(m_ResultMatrix.toStringMatrix());
    }
    catch (Exception e) {
      e.printStackTrace();
    }
    
    // append a key so that we can tell the difference between long
    // scheme+option names
    result.append("\n\n" + m_ResultMatrix.toStringKey());

    return result.toString();
  }

  /**
   * Lists options understood by this object.
   *
   * @return an enumeration of Options.
   */
  public Enumeration listOptions() {
    
    Vector newVector = new Vector();

    newVector.addElement(new Option(
             "\tSpecify list of columns that specify a unique\n"
	      + "\tdataset.\n"
	      + "\tFirst and last are valid indexes. (default none)",
              "D", 1, "-D <index,index2-index4,...>"));
    newVector.addElement(new Option(
	      "\tSet the index of the column containing the run number",
              "R", 1, "-R <index>"));
    newVector.addElement(new Option(
	      "\tSet the index of the column containing the fold number",
              "F", 1, "-F <index>"));
    newVector.addElement(new Option(
              "\tSpecify list of columns that specify a unique\n"
	      + "\t'result generator' (eg: classifier name and options).\n"
	      + "\tFirst and last are valid indexes. (default none)",
              "G", 1, "-G <index1,index2-index4,...>"));
    newVector.addElement(new Option(
	      "\tSet the significance level for comparisons (default 0.05)",
              "S", 1, "-S <significance level>"));
    newVector.addElement(new Option(
	      "\tShow standard deviations",
              "V", 0, "-V"));
    newVector.addElement(new Option(
	      "\tProduce table comparisons in Latex table format",
              "L", 0, "-L"));
    newVector.addElement(new Option(
         "\tProduce table comparisons in CSV table format",
         "csv", 0, "-csv"));
    newVector.addElement(new Option(
         "\tProduce table comparisons in HTML table format",
         "html", 0, "-html"));
    newVector.addElement(new Option(
         "\tProduce table comparisons with only the significance values",
         "significance", 0, "-significance"));
    newVector.addElement(new Option(
         "\tProduce table comparisons output suitable for GNUPlot",
         "gnuplot", 0, "-gnuplot"));

    return newVector.elements();
  }

  /**
   * Parses a given list of options. <p/>
   *
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -D &lt;index,index2-index4,...&gt;
   *  Specify list of columns that specify a unique
   *  dataset.
   *  First and last are valid indexes. (default none)</pre>
   * 
   * <pre> -R &lt;index&gt;
   *  Set the index of the column containing the run number</pre>
   * 
   * <pre> -F &lt;index&gt;
   *  Set the index of the column containing the fold number</pre>
   * 
   * <pre> -G &lt;index1,index2-index4,...&gt;
   *  Specify list of columns that specify a unique
   *  'result generator' (eg: classifier name and options).
   *  First and last are valid indexes. (default none)</pre>
   * 
   * <pre> -S &lt;significance level&gt;
   *  Set the significance level for comparisons (default 0.05)</pre>
   * 
   * <pre> -V
   *  Show standard deviations</pre>
   * 
   * <pre> -L
   *  Produce table comparisons in Latex table format</pre>
   * 
   * <pre> -csv
   *  Produce table comparisons in CSV table format</pre>
   * 
   * <pre> -html
   *  Produce table comparisons in HTML table format</pre>
   * 
   * <pre> -significance
   *  Produce table comparisons with only the significance values</pre>
   * 
   * <pre> -gnuplot
   *  Produce table comparisons output suitable for GNUPlot</pre>
   * 
   <!-- options-end -->
   *
   * @param options an array containing options to set.
   * @throws Exception if invalid options are given
   */
  public void setOptions(String[] options) throws Exception {

    setShowStdDevs(Utils.getFlag('V', options));
    if (Utils.getFlag('L', options))
      setResultMatrix(new ResultMatrixLatex());
    if (Utils.getFlag("csv", options))
      setResultMatrix(new ResultMatrixCSV());
    if (Utils.getFlag("html", options))
      setResultMatrix(new ResultMatrixHTML());
    if (Utils.getFlag("significance", options))
      setResultMatrix(new ResultMatrixSignificance());

    String datasetList = Utils.getOption('D', options);
    Range datasetRange = new Range();
    if (datasetList.length() != 0) {
      datasetRange.setRanges(datasetList);
    }
    setDatasetKeyColumns(datasetRange);

    String indexStr = Utils.getOption('R', options);
    if (indexStr.length() != 0) {
      if (indexStr.equals("first")) {
	setRunColumn(0);
      } else if (indexStr.equals("last")) {
	setRunColumn(-1);
      } else {
	setRunColumn(Integer.parseInt(indexStr) - 1);
      }    
    } else {
      setRunColumn(-1);
    }

    String foldStr = Utils.getOption('F', options);
    if (foldStr.length() != 0) {
      setFoldColumn(Integer.parseInt(foldStr) - 1);
    } else {
      setFoldColumn(-1);
    }

    String sigStr = Utils.getOption('S', options);
    if (sigStr.length() != 0) {
      setSignificanceLevel((new Double(sigStr)).doubleValue());
    } else {
      setSignificanceLevel(0.05);
    }
    
    String resultsetList = Utils.getOption('G', options);
    Range generatorRange = new Range();
    if (resultsetList.length() != 0) {
      generatorRange.setRanges(resultsetList);
    }
    setResultsetKeyColumns(generatorRange);
  }
  
  /**
   * Gets current settings of the PairedTTester.
   *
   * @return an array of strings containing current options.
   */
  public String[] getOptions() {

    String [] options = new String [11];
    int current = 0;

    if (!getResultsetKeyColumns().getRanges().equals("")) {
      options[current++] = "-G";
      options[current++] = getResultsetKeyColumns().getRanges();
    }
    if (!getDatasetKeyColumns().getRanges().equals("")) {
      options[current++] = "-D";
      options[current++] = getDatasetKeyColumns().getRanges();
    }
    options[current++] = "-R";
    options[current++] = "" + (getRunColumn() + 1);
    options[current++] = "-S";
    options[current++] = "" + getSignificanceLevel();
    
    if (getShowStdDevs()) {
      options[current++] = "-V";
    }

    if (getResultMatrix() instanceof ResultMatrixLatex)
      options[current++] = "-L";

    if (getResultMatrix() instanceof ResultMatrixCSV)
      options[current++] = "-csv";
   
    if (getResultMatrix() instanceof ResultMatrixHTML)
      options[current++] = "-html";
   
    if (getResultMatrix() instanceof ResultMatrixSignificance)
      options[current++] = "-significance";
   
    while (current < options.length) {
      options[current++] = "";
    }
    return options;
  }

  /**
   * Get the value of ResultsetKeyColumns.
   *
   * @return Value of ResultsetKeyColumns.
   */
  public Range getResultsetKeyColumns() {
    
    return m_ResultsetKeyColumnsRange;
  }
  
  /**
   * Set the value of ResultsetKeyColumns.
   *
   * @param newResultsetKeyColumns Value to assign to ResultsetKeyColumns.
   */
  public void setResultsetKeyColumns(Range newResultsetKeyColumns) {
    
    m_ResultsetKeyColumnsRange = newResultsetKeyColumns;
    m_ResultsetsValid = false;
  }
  
  /**
   * Gets the indices of the the datasets that are displayed (if <code>null</code>
   * then all are displayed). The base is always displayed.
   * 
   * @return the indices of the datasets to display
   */
  public int[] getDisplayedResultsets() {
    return m_DisplayedResultsets;
  }
  
  /**
   * Sets the indicies of the datasets to display (<code>null</code> means all).
   * The base is always displayed.
   * 
   * @param cols the indices of the datasets to display
   */
  public void setDisplayedResultsets(int[] cols) {
    m_DisplayedResultsets = cols;
  }
  
  /**
   * Get the value of SignificanceLevel.
   *
   * @return Value of SignificanceLevel.
   */
  public double getSignificanceLevel() {
    
    return m_SignificanceLevel;
  }
  
  /**
   * Set the value of SignificanceLevel.
   *
   * @param newSignificanceLevel Value to assign to SignificanceLevel.
   */
  public void setSignificanceLevel(double newSignificanceLevel) {
    
    m_SignificanceLevel = newSignificanceLevel;
  }

  /**
   * Get the value of DatasetKeyColumns.
   *
   * @return Value of DatasetKeyColumns.
   */
  public Range getDatasetKeyColumns() {
    
    return m_DatasetKeyColumnsRange;
  }
  
  /**
   * Set the value of DatasetKeyColumns.
   *
   * @param newDatasetKeyColumns Value to assign to DatasetKeyColumns.
   */
  public void setDatasetKeyColumns(Range newDatasetKeyColumns) {
    
    m_DatasetKeyColumnsRange = newDatasetKeyColumns;
    m_ResultsetsValid = false;
  }
  
  /**
   * Get the value of RunColumn.
   *
   * @return Value of RunColumn.
   */
  public int getRunColumn() {
    
    return m_RunColumnSet;
  }
  
  /**
   * Set the value of RunColumn.
   *
   * @param newRunColumn Value to assign to RunColumn.
   */
  public void setRunColumn(int newRunColumn) {
    
    m_RunColumnSet = newRunColumn;
    m_ResultsetsValid = false;
  }

  /**
   * Get the value of FoldColumn.
   *
   * @return Value of FoldColumn.
   */
  public int getFoldColumn() {
    
    return m_FoldColumn;
  }
  
  /**
   * Set the value of FoldColumn.
   *
   * @param newFoldColumn Value to assign to FoldColumn.
   */
  public void setFoldColumn(int newFoldColumn) {
    
    m_FoldColumn = newFoldColumn;
    m_ResultsetsValid = false;
  }

  /**
   * Returns the name of the column to sort on.
   *
   * @return the name of the column to sort on.
   */
  public String getSortColumnName() {
    if (getSortColumn() == -1)
      return "-";
    else
      return m_Instances.attribute(getSortColumn()).name();
  }

  /**
   * Returns the column to sort on, -1 means the default sorting.
   *
   * @return the column to sort on.
   */
  public int getSortColumn() {
    return m_SortColumn;
  }
  
  /**
   * Set the column to sort on, -1 means the default sorting.
   *
   * @param newSortColumn the new sort column.
   */
  public void setSortColumn(int newSortColumn) {
    if (newSortColumn >= -1)
      m_SortColumn = newSortColumn;
  }
  
  /**
   * Get the value of Instances.
   *
   * @return Value of Instances.
   */
  public Instances getInstances() {
    
    return m_Instances;
  }
  
  /**
   * Set the value of Instances.
   *
   * @param newInstances Value to assign to Instances.
   */
  public void setInstances(Instances newInstances) {
    
    m_Instances = newInstances;
    m_ResultsetsValid = false;
  }

  /**
   * retrieves all the settings from the given Tester
   *
   * @param tester      the Tester to get the settings from
   */
  public void assign(Tester tester) {
    setInstances(tester.getInstances());
    setResultMatrix(tester.getResultMatrix());
    setShowStdDevs(tester.getShowStdDevs());
    setResultsetKeyColumns(tester.getResultsetKeyColumns());
    setDisplayedResultsets(tester.getDisplayedResultsets());
    setSignificanceLevel(tester.getSignificanceLevel());
    setDatasetKeyColumns(tester.getDatasetKeyColumns());
    setRunColumn(tester.getRunColumn());
    setFoldColumn(tester.getFoldColumn());
    setSortColumn(tester.getSortColumn());
  }

  /**
   * returns a string that is displayed as tooltip on the "perform test"
   * button in the experimenter
   * 
   * @return	the tool tip
   */
  public String getToolTipText() {
    return "Performs test using t-test statistic";
  }

  /**
   * returns the name of the tester
   * 
   * @return	the display name
   */
  public String getDisplayName() {
    return "Paired T-Tester";
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6432 $");
  }
  
  /**
   * Test the class from the command line.
   *
   * @param args contains options for the instance ttests
   */
  public static void main(String args[]) {

    try {
      PairedTTester tt = new PairedTTester();
      String datasetName = Utils.getOption('t', args);
      String compareColStr = Utils.getOption('c', args);
      String baseColStr = Utils.getOption('b', args);
      boolean summaryOnly = Utils.getFlag('s', args);
      boolean rankingOnly = Utils.getFlag('r', args);
      try {
	if ((datasetName.length() == 0)
	    || (compareColStr.length() == 0)) {
	  throw new Exception("-t and -c options are required");
	}
	tt.setOptions(args);
	Utils.checkForRemainingOptions(args);
      } catch (Exception ex) {
	String result = "";
	Enumeration enu = tt.listOptions();
	while (enu.hasMoreElements()) {
	  Option option = (Option) enu.nextElement();
	  result += option.synopsis() + '\n'
	    + option.description() + '\n';
	}
	throw new Exception(
	      "Usage:\n\n"
	      + "-t <file>\n"
	      + "\tSet the dataset containing data to evaluate\n"
	      + "-b <index>\n"
	      + "\tSet the resultset to base comparisons against (optional)\n"
	      + "-c <index>\n"
	      + "\tSet the column to perform a comparison on\n"
	      + "-s\n"
	      + "\tSummarize wins over all resultset pairs\n\n"
	      + "-r\n"
	      + "\tGenerate a resultset ranking\n\n"
	      + result);
      }
      Instances data = new Instances(new BufferedReader(
				  new FileReader(datasetName)));
      tt.setInstances(data);
      //      tt.prepareData();
      int compareCol = Integer.parseInt(compareColStr) - 1;
      System.out.println(tt.header(compareCol));
      if (rankingOnly) {
	System.out.println(tt.multiResultsetRanking(compareCol));
      } else if (summaryOnly) {
	System.out.println(tt.multiResultsetSummary(compareCol));
      } else {
	System.out.println(tt.resultsetKey());
	if (baseColStr.length() == 0) {
	  for (int i = 0; i < tt.getNumResultsets(); i++) {
            if (!tt.displayResultset(i))
              continue;
	    System.out.println(tt.multiResultsetFull(i, compareCol));
	  }
	} else {
	  int baseCol = Integer.parseInt(baseColStr) - 1;
	  System.out.println(tt.multiResultsetFull(baseCol, compareCol));
	}
      }
    } catch(Exception e) {
      e.printStackTrace();
      System.err.println(e.getMessage());
    }
  }
}
