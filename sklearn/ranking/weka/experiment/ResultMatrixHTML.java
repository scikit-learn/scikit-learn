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
 * ResultMatrixHTML.java
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.experiment;

import weka.core.RevisionUtils;
import weka.core.Utils;

/**
 <!-- globalinfo-start -->
 * Generates the matrix output as HTML.
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
 *  (default: 0)</pre>
 * 
 * <pre> -row-name-width &lt;int&gt;
 *  The maximum width for the row names (0 = optimal).
 *  (default: 25)</pre>
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
public class ResultMatrixHTML
  extends ResultMatrix {

  /** for serialization. */
  private static final long serialVersionUID = 6672380422544799990L;

  /**
   * initializes the matrix as 1x1 matrix.
   */
  public ResultMatrixHTML() {
    this(1, 1);
  }

  /**
   * initializes the matrix with the given dimensions.
   * 
   * @param cols	the number of columns
   * @param rows	the number of rows
   */
  public ResultMatrixHTML(int cols, int rows) {
    super(cols, rows);
  }

  /**
   * initializes the matrix with the values from the given matrix.
   * 
   * @param matrix      the matrix to get the values from
   */
  public ResultMatrixHTML(ResultMatrix matrix) {
    super(matrix);
  }
  
  /**
   * Returns a string describing the matrix.
   * 
   * @return 		a description suitable for
   * 			displaying in the experimenter gui
   */
  public String globalInfo() {
    return "Generates the matrix output as HTML.";
  }

  /**
   * returns the name of the output format.
   * 
   * @return		the display name
   */
  public String getDisplayName() {
    return "HTML";
  }

  /**
   * returns the default width for the row names.
   * 
   * @return		the width
   */
  public int getDefaultRowNameWidth() {
    return 25;
  }

  /**
   * returns the default of whether column names or numbers instead are printed.
   * 
   * @return		true if names instead of numbers are printed
   */
  public boolean getDefaultPrintColNames() {
    return false;
  }

  /**
   * returns the default of whether column names are prefixed with the index.
   * 
   * @return		true if the names are prefixed
   */
  public boolean getDefaultEnumerateColNames() {
    return true;
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
   * returns the matrix in an HTML table.
   * 
   * @return		the matrix
   */
  public String toStringMatrix() {
    StringBuffer        result;
    String[][]          cells;
    int                 i;
    int                 n;
    int                 cols;

    result = new StringBuffer();
    cells  = toArray();

    result.append("<table border=\"1\" cellpadding=\"3\" cellspacing=\"0\">\n");
    
    // headings
    result.append("   <tr>");
    for (n = 0; n < cells[0].length; n++) {
      if (isRowName(n)) {
        result.append("<td><b>" + cells[0][n] + "</b></td>");
      }
      else if (isMean(n)) {
        if (n == 1)
          cols = 1;
        else
          cols = 2;
        if (getShowStdDev())
          cols++;
        result.append("<td align=\"center\" colspan=\"" + cols + "\">");
        result.append("<b>" + cells[0][n] + "</b>");
        result.append("</td>");
      }
    }
    result.append("</tr>\n");
      
    // data
    for (i = 1; i < cells.length; i++) {
      result.append("   <tr>");
      for (n = 0; n < cells[i].length; n++) {
        if (isRowName(n))
          result.append("<td>");
        else if (isMean(n) || isStdDev(n))
          result.append("<td align=\"right\">");
        else if (isSignificance(n))
          result.append("<td align=\"center\">");
        else
          result.append("<td>");
        
        // content
        if (cells[i][n].trim().equals(""))
          result.append("&nbsp;");
        else if (isStdDev(n))
          result.append("&plusmn;&nbsp;" + cells[i][n]);
        else
          result.append(cells[i][n]);
        
        result.append("</td>");
      }
      result.append("</tr>\n");
    }
    result.append("</table>\n");
    
    return result.toString();
  }

  /**
   * returns returns a key for all the col names, for better readability if
   * the names got cut off.
   * 
   * @return		the key
   */
  public String toStringKey() {
    String          result;
    int             i;

    result =   "<table border=\"1\" cellpadding=\"3\" cellspacing=\"0\">\n" 
             + "   <tr><td colspan=\"2\"><b>Key</b></td></tr>\n";
    for (i = 0; i < getColCount(); i++) {
      if (getColHidden(i))
        continue;

      result +=   "   <tr>"
                + "<td><b>(" + (i+1) + ")</b></td>"
                + "<td>" + removeFilterName(m_ColNames[i]) + "</td>" 
                + "</tr>\n";
    }

    result += "</table>\n";

    return result;
  }

  /**
   * returns the summary as string.
   * 
   * @return		the summary
   */
  public String toStringSummary() {
    String      result;
    String      titles;
    int         resultsetLength;
    int         i;
    int         j;
    String      content;

    if (m_NonSigWins == null)
      return "-summary data not set-";
    
    result = "<table border=\"1\" cellpadding=\"3\" cellspacing=\"0\">\n";
    titles = "   <tr>";
    resultsetLength = 1 + Math.max((int)(Math.log(getColCount())/Math.log(10)),
                                   (int)(Math.log(getRowCount())/Math.log(10)));

    for (i = 0; i < getColCount(); i++) {
      if (getColHidden(i))
        continue;
      titles += "<td align=\"center\"><b>" + getSummaryTitle(i) + "</b></td>";
    }
    result +=   titles 
              + "<td><b>(No. of datasets where [col] &gt;&gt; [row])</b></td></tr>\n";

    for (i = 0; i < getColCount(); i++) {
      if (getColHidden(i))
        continue;

      result += "   <tr>";

      for (j = 0; j < getColCount(); j++) {
        if (getColHidden(j))
          continue;

	if (j == i)
	  content = Utils.padLeft("-", resultsetLength * 2 + 3);
	else
	  content = Utils.padLeft("" + m_NonSigWins[i][j] 
                                  + " (" + m_Wins[i][j] + ")",
				  resultsetLength * 2 + 3);
        result += "<td>" + content.replaceAll(" ", "&nbsp;") + "</td>";
      }

      result += "<td><b>" + getSummaryTitle(i) + "</b> = " + removeFilterName(m_ColNames[i]) + "</td></tr>\n";
    }

    result += "</table>\n";

    return result;
  }

  /**
   * returns the ranking in a string representation.
   * 
   * @return		the ranking
   */
  public String toStringRanking() {
    String        result;
    int[]         ranking;
    int           i;
    int           curr;

    if (m_RankingWins == null)
      return "-ranking data not set-";

    result = "<table border=\"1\" cellpadding=\"3\" cellspacing=\"0\">\n";
    result +=  "   <tr>" 
             + "<td align=\"center\"><b>&gt;-&lt;</b></td>"
             + "<td align=\"center\"><b>&gt;</b></td>"
             + "<td align=\"center\"><b>&lt;</b></td>"
             + "<td><b>Resultset</b></td>"
             + "</tr>\n";

    ranking = Utils.sort(m_RankingDiff);

    for (i = getColCount() - 1; i >= 0; i--) {
      curr = ranking[i];

      if (getColHidden(curr))
        continue;

      result += "   <tr>"
        + "<td align=\"right\">" + m_RankingDiff[curr] + "</td>"
        + "<td align=\"right\">" + m_RankingWins[curr] + "</td>"
        + "<td align=\"right\">" + m_RankingLosses[curr] + "</td>"
        + "<td>" + removeFilterName(m_ColNames[curr]) + "</td>"
        + "<tr>\n";
    }

    result += "</table>\n";

    return result;
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
    
    matrix = new ResultMatrixHTML(3, 3);
    
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
