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
 * ResultMatrixGnuPlot.java
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.experiment;

import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.core.Version;

/**
 <!-- globalinfo-start -->
 * Generates output for a data and script file for GnuPlot.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -mean-prec &lt;int&gt;
 *  The number of decimals after the decimal point for the mean.
 *  (default: 2)</pre>
 * 
 * <pre> -stddev-prec &lt;int&gt;
 *  The number of decimals after the decimal point for the mean.
 *  (default: 2)</pre>
 * 
 * <pre> -col-name-width &lt;int&gt;
 *  The maximum width for the column names (0 = optimal).
 *  (default: 50)</pre>
 * 
 * <pre> -row-name-width &lt;int&gt;
 *  The maximum width for the row names (0 = optimal).
 *  (default: 50)</pre>
 * 
 * <pre> -mean-width &lt;int&gt;
 *  The width of the mean (0 = optimal).
 *  (default: 0)</pre>
 * 
 * <pre> -stddev-width &lt;int&gt;
 *  The width of the standard deviation (0 = optimal).
 *  (default: 0)</pre>
 * 
 * <pre> -sig-width &lt;int&gt;
 *  The width of the significance indicator (0 = optimal).
 *  (default: 0)</pre>
 * 
 * <pre> -count-width &lt;int&gt;
 *  The width of the counts (0 = optimal).
 *  (default: 0)</pre>
 * 
 * <pre> -show-stddev
 *  Whether to display the standard deviation column.
 *  (default: no)</pre>
 * 
 * <pre> -show-avg
 *  Whether to show the row with averages.
 *  (default: no)</pre>
 * 
 * <pre> -remove-filter
 *  Whether to remove the classname package prefixes from the
 *  filter names in datasets.
 *  (default: no)</pre>
 * 
 * <pre> -print-col-names
 *  Whether to output column names or just numbers representing them.
 *  (default: no)</pre>
 * 
 * <pre> -print-row-names
 *  Whether to output row names or just numbers representing them.
 *  (default: no)</pre>
 * 
 * <pre> -enum-col-names
 *  Whether to enumerate the column names (prefixing them with 
 *  '(x)', with 'x' being the index).
 *  (default: no)</pre>
 * 
 * <pre> -enum-row-names
 *  Whether to enumerate the row names (prefixing them with 
 *  '(x)', with 'x' being the index).
 *  (default: no)</pre>
 * 
 <!-- options-end -->
 *
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5346 $
 */
public class ResultMatrixGnuPlot
  extends ResultMatrix {

  /** for serialization. */
  private static final long serialVersionUID = -234648254944790097L;
  
  /**
   * initializes the matrix as 1x1 matrix.
   */
  public ResultMatrixGnuPlot() {
    this(1, 1);
  }

  /**
   * initializes the matrix with the given dimensions.
   * 
   * @param cols	the number of columns
   * @param rows	the number of rows
   */
  public ResultMatrixGnuPlot(int cols, int rows) {
    super(cols, rows);
  }

  /**
   * initializes the matrix with the values from the given matrix.
   * 
   * @param matrix      the matrix to get the values from
   */
  public ResultMatrixGnuPlot(ResultMatrix matrix) {
    super(matrix);
  }
  
  /**
   * Returns a string describing the matrix.
   * 
   * @return 		a description suitable for
   * 			displaying in the experimenter gui
   */
  public String globalInfo() {
    return "Generates output for a data and script file for GnuPlot.";
  }

  /**
   * returns the name of the output format.
   * 
   * @return		the display name
   */
  public String getDisplayName() {
    return "GNUPlot";
  }

  /**
   * removes the stored data but retains the dimensions of the matrix.
   */
  public void clear() {
    super.clear();
    LEFT_PARENTHESES = "";
    RIGHT_PARENTHESES = "";
  }

  /**
   * returns the default width for the row names.
   * 
   * @return		the width
   */
  public int getDefaultRowNameWidth() {
    return 50;
  }

  /**
   * returns the default width for the column names.
   * 
   * @return		the width
   */
  public int getDefaultColNameWidth() {
    return 50;
  }

  /**
   * returns the default of whether column names are prefixed with the index.
   * 
   * @return		true if the names are prefixed
   */
  public boolean getDefaultEnumerateColNames() {
    return false;
  }

  /**
   * returns the default of whether row names are prefixed with the index.
   * 
   * @return		true if the names are prefixed
   */
  public boolean getDefaultEnumerateRowNames() {
    return false;
  }
  
  /**
   * returns the header of the matrix as a string.
   * 
   * @return		the header
   * @see 		#m_HeaderKeys
   * @see 		#m_HeaderValues
   */
  public String toStringHeader() {
    return new ResultMatrixPlainText(this).toStringHeader();
  }

  /**
   * returns the matrix in CSV format.
   * 
   * @return		the matrix
   */
  public String toStringMatrix() {
    StringBuffer        result;
    String[][]          cells;
    int                 i;
    int                 n;
    String              line;
    String              title;
    String              generated;

    result = new StringBuffer();
    cells  = toArray();

    // generation comment
    generated = "# generated by WEKA " + Version.VERSION + "\n";

    // data
    result.append("\n");
    result.append("##################\n");
    result.append("# file: plot.dat #\n");
    result.append("##################\n");
    result.append(generated);
    result.append("# contains the data for the plot\n");
    // key for x-axis
    result.append("\n");
    result.append("# key for the x-axis\n");
    for (i = 1; i < cells.length - 1; i++)
      result.append("# " + i + " - " + cells[i][0] + "\n");
    // the data itself
    result.append("\n");
    result.append("# data for the plot\n");
    for (i = 1; i < cells.length - 1; i++) {
      result.append(Integer.toString(i));
      for (n = 1; n < cells[i].length; n++) {
        if (isSignificance(n))
          continue;
        result.append(" ");
        result.append(Utils.quote(cells[i][n]));
      }
      result.append("\n");
    }
    result.append("#######\n");
    result.append("# end #\n");
    result.append("#######\n");

    // script
    result.append("\n");
    result.append("##################\n");
    result.append("# file: plot.scr #\n");
    result.append("##################\n");
    result.append(generated);
    result.append("# script to plot the data\n");
    result.append("\n");
    result.append("# display it in a window:\n");
    result.append("set terminal x11\n");
    result.append("set output\n");
    result.append("\n");
    result.append("# to display all data rows:\n");
    result.append("set xrange [0:" + ((cells.length - 2) + 1) + "]\n");
    result.append("\n");
    result.append("# axis labels, e.g.:\n");
    result.append("#set xlabel \"Datasets\"\n");
    result.append("#set ylabel \"Accuracy in %\"\n");
    result.append("\n");
    result.append("# the plot commands\n");
    n = 1;
    i = 0;
    while (i < cells[0].length - 1) {
      i++;

      if (isSignificance(i))
        continue;

      n++;
      
      // plot
      if (i == 1)
        line = "plot";
      else
        line = "replot";
      line += " \"plot.dat\"";

      // title
      title = "title \"" + cells[0][i] + "\"";
      
      // columns
      line += " using 1:" + n;
      if (getShowStdDev()) {
        n++;
        i++;
        // errorbars
        line += ":" + n;
      }
      
      // options
      line += " with";
      if (getShowStdDev())
        line += " yerrorbars";
      else
        line += " lines";
      line += " " + title;
      
      result.append(line + "\n");
    }
    result.append("\n");
    result.append("# generate ps:\n");
    result.append("#set terminal postscript\n");
    result.append("#set output \"plot.ps\"\n");
    result.append("#replot\n");
    result.append("\n");
    result.append("# generate png:\n");
    result.append("#set terminal png size 800,600\n");
    result.append("#set output \"plot.png\"\n");
    result.append("#replot\n");
    result.append("\n");
    result.append("# wait for user to hit <Return>\n");
    result.append("pause -1\n");
    result.append("#######\n");
    result.append("# end #\n");
    result.append("#######\n");
    
    return result.toString();
  }

  /**
   * returns returns a key for all the col names, for better readability if
   * the names got cut off.
   * 
   * @return		the key
   */
  public String toStringKey() {
    return new ResultMatrixPlainText(this).toStringKey();
  }

  /**
   * returns the summary as string.
   * 
   * @return		the summary
   */
  public String toStringSummary() {
    return new ResultMatrixPlainText(this).toStringSummary();
  }

  /**
   * returns the ranking in a string representation.
   * 
   * @return		the ranking
   */
  public String toStringRanking() {
    return new ResultMatrixPlainText(this).toStringRanking();
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5346 $");
  }

  /**
   * for testing only.
   * 
   * @param args	ignored
   */
  public static void main(String[] args) {
    ResultMatrix        matrix;
    int                 i;
    int                 n;
    
    matrix = new ResultMatrixGnuPlot(3, 3);

    // set header
    matrix.addHeader("header1", "value1");
    matrix.addHeader("header2", "value2");
    matrix.addHeader("header2", "value3");
    
    // set values
    for (i = 0; i < matrix.getRowCount(); i++) {
      for (n = 0; n < matrix.getColCount(); n++) {
        matrix.setMean(n, i, (i+1)*n);
        matrix.setStdDev(n, i, ((double) (i+1)*n) / 100);
        if (i == n) {
          if (i % 2 == 1)
            matrix.setSignificance(n, i, SIGNIFICANCE_WIN);
          else
            matrix.setSignificance(n, i, SIGNIFICANCE_LOSS);
        }
      }
    }

    System.out.println("\n\n--> " + matrix.getDisplayName());
    
    System.out.println("\n1. complete\n");
    System.out.println(matrix.toStringHeader() + "\n");
    System.out.println(matrix.toStringMatrix() + "\n");
    System.out.println(matrix.toStringKey());
    
    System.out.println("\n2. complete with std deviations\n");
    matrix.setShowStdDev(true);
    System.out.println(matrix.toStringMatrix());
    
    System.out.println("\n3. cols numbered\n");
    matrix.setPrintColNames(false);
    System.out.println(matrix.toStringMatrix());
    
    System.out.println("\n4. second col missing\n");
    matrix.setColHidden(1, true);
    System.out.println(matrix.toStringMatrix());
    
    System.out.println("\n5. last row missing, rows numbered too\n");
    matrix.setRowHidden(2, true);
    matrix.setPrintRowNames(false);
    System.out.println(matrix.toStringMatrix());
    
    System.out.println("\n6. mean prec to 3\n");
    matrix.setMeanPrec(3);
    matrix.setPrintRowNames(false);
    System.out.println(matrix.toStringMatrix());
  }
}
