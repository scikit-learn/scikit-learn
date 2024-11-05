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
 * Tester.java
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.experiment;

import weka.core.Instance;
import weka.core.Instances;
import weka.core.Range;

import java.io.Serializable;

/**
 * Interface for different kinds of Testers in the Experimenter.
 * 
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 1.4 $
 */
public interface Tester
  extends Serializable {
  
  /**
   * returns the name of the testing algorithm
   */
  public String getDisplayName();

  /**
   * returns a string that is displayed as tooltip on the "perform test"
   * button in the experimenter
   */
  public String getToolTipText();

  /**
   * retrieves all the settings from the given Tester
   *
   * @param tester      the Tester to get the settings from
   */
  public void assign(Tester tester);

  /**
   * Sets the matrix to use to produce the output.
   * @param matrix the instance to use to produce the output
   * @see ResultMatrix
   */
  public void setResultMatrix(ResultMatrix matrix);

  /**
   * Gets the instance that produces the output.
   * @return the instance to produce the output
   */
  public ResultMatrix getResultMatrix();

  /**
   * Set whether standard deviations are displayed or not.
   * @param s true if standard deviations are to be displayed
   */
  public void setShowStdDevs(boolean s);

  /**
   * Returns true if standard deviations have been requested.
   * @return true if standard deviations are to be displayed.
   */
  public boolean getShowStdDevs();

  /**
   * Gets the number of datasets in the resultsets
   *
   * @return the number of datasets in the resultsets
   */
  public int getNumDatasets();

  /**
   * Gets the number of resultsets in the data.
   *
   * @return the number of resultsets in the data
   */
  public int getNumResultsets();

  /**
   * Gets a string descriptive of the specified resultset.
   *
   * @param index the index of the resultset
   * @return a descriptive string for the resultset
   */
  public String getResultsetName(int index);
  
  /**
   * Checks whether the resultset with the given index shall be displayed.
   * 
   * @param index the index of the resultset to check whether it shall be displayed 
   * @return whether the specified resultset is displayed 
   */
  public boolean displayResultset(int index);
  
  /**
   * Computes a paired t-test comparison for a specified dataset between
   * two resultsets.
   *
   * @param datasetSpecifier the dataset specifier
   * @param resultset1Index the index of the first resultset
   * @param resultset2Index the index of the second resultset
   * @param comparisonColumn the column containing values to compare
   * @return the results of the paired comparison
   * @exception Exception if an error occurs
   */
  public PairedStats calculateStatistics(Instance datasetSpecifier,
					 int resultset1Index,
					 int resultset2Index,
					 int comparisonColumn) throws Exception;
  
  /**
   * Creates a key that maps resultset numbers to their descriptions.
   *
   * @return a value of type 'String'
   */
  public String resultsetKey();
  
  /**
   * Creates a "header" string describing the current resultsets.
   *
   * @param comparisonColumn a value of type 'int'
   * @return a value of type 'String'
   */
  public String header(int comparisonColumn);

  /**
   * Carries out a comparison between all resultsets, counting the number
   * of datsets where one resultset outperforms the other.
   *
   * @param comparisonColumn the index of the comparison column
   * @return a 2d array where element [i][j] is the number of times resultset
   * j performed significantly better than resultset i.
   * @exception Exception if an error occurs
   */
  public int [][] multiResultsetWins(int comparisonColumn, int [][] nonSigWin)
    throws Exception;
  
  /**
   * Carries out a comparison between all resultsets, counting the number
   * of datsets where one resultset outperforms the other. The results
   * are summarized in a table.
   *
   * @param comparisonColumn the index of the comparison column
   * @return the results in a string
   * @exception Exception if an error occurs
   */
  public String multiResultsetSummary(int comparisonColumn)
    throws Exception;

  public String multiResultsetRanking(int comparisonColumn)
    throws Exception;
				    
  /**
   * Creates a comparison table where a base resultset is compared to the
   * other resultsets. Results are presented for every dataset.
   *
   * @param baseResultset the index of the base resultset
   * @param comparisonColumn the index of the column to compare over
   * @return the comparison table string
   * @exception Exception if an error occurs
   */
  public String multiResultsetFull(int baseResultset,
				   int comparisonColumn) throws Exception;

  /**
   * Get the value of ResultsetKeyColumns.
   *
   * @return Value of ResultsetKeyColumns.
   */
  public Range getResultsetKeyColumns();
  
  /**
   * Set the value of ResultsetKeyColumns.
   *
   * @param newResultsetKeyColumns Value to assign to ResultsetKeyColumns.
   */
  public void setResultsetKeyColumns(Range newResultsetKeyColumns);
  
  /**
   * Gets the indices of the the datasets that are displayed (if <code>null</code>
   * then all are displayed). The base is always displayed.
   * 
   * @return the indices of the datasets to display
   */
  public int[] getDisplayedResultsets();
  
  /**
   * Sets the indicies of the datasets to display (<code>null</code> means all).
   * The base is always displayed.
   * 
   * @param cols the indices of the datasets to display
   */
  public void setDisplayedResultsets(int[] cols);
  
  /**
   * Get the value of SignificanceLevel.
   *
   * @return Value of SignificanceLevel.
   */
  public double getSignificanceLevel();
  
  /**
   * Set the value of SignificanceLevel.
   *
   * @param newSignificanceLevel Value to assign to SignificanceLevel.
   */
  public void setSignificanceLevel(double newSignificanceLevel);

  /**
   * Get the value of DatasetKeyColumns.
   *
   * @return Value of DatasetKeyColumns.
   */
  public Range getDatasetKeyColumns();
  
  /**
   * Set the value of DatasetKeyColumns.
   *
   * @param newDatasetKeyColumns Value to assign to DatasetKeyColumns.
   */
  public void setDatasetKeyColumns(Range newDatasetKeyColumns);
  
  /**
   * Get the value of RunColumn.
   *
   * @return Value of RunColumn.
   */
  public int getRunColumn();
  
  /**
   * Set the value of RunColumn.
   *
   * @param newRunColumn Value to assign to RunColumn.
   */
  public void setRunColumn(int newRunColumn);

  /**
   * Get the value of FoldColumn.
   *
   * @return Value of FoldColumn.
   */
  public int getFoldColumn();
  
  /**
   * Set the value of FoldColumn.
   *
   * @param newFoldColumn Value to assign to FoldColumn.
   */
  public void setFoldColumn(int newFoldColumn);

  /**
   * Returns the name of the column to sort on.
   *
   * @return the name of the column to sort on.
   */
  public String getSortColumnName();

  /**
   * Returns the column to sort on, -1 means the default sorting.
   *
   * @return the column to sort on.
   */
  public int getSortColumn();
  
  /**
   * Set the column to sort on, -1 means the default sorting.
   *
   * @param newSortColumn the new sort column.
   */
  public void setSortColumn(int newSortColumn);
  
  /**
   * Get the value of Instances.
   *
   * @return Value of Instances.
   */
  public Instances getInstances();
  
  /**
   * Set the value of Instances.
   *
   * @param newInstances Value to assign to Instances.
   */
  public void setInstances(Instances newInstances);
}
