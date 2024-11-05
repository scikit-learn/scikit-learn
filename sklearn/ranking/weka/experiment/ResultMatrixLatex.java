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
 * ResultMatrixLatex.java
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.experiment;

import weka.core.RevisionUtils;
import weka.core.Utils;

/**
 <!-- globalinfo-start -->
 * Generates the matrix output in LaTeX-syntax.
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
 *  (default: 0)</pre>
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
public class ResultMatrixLatex
  extends ResultMatrix {

  /** for serialization. */
  private static final long serialVersionUID = 777690788447600978L;
  
  /**
   * initializes the matrix as 1x1 matrix.
   */
  public ResultMatrixLatex() {
    this(1, 1);
  }

  /**
   * initializes the matrix with the given dimensions.
   * 
   * @param cols	the number of columns
   * @param rows	the number of rows
   */
  public ResultMatrixLatex(int cols, int rows) {
    super(cols, rows);
  }

  /**
   * initializes the matrix with the values from the given matrix.
   * 
   * @param matrix      the matrix to get the values from
   */
  public ResultMatrixLatex(ResultMatrix matrix) {
    super(matrix);
  }
  
  /**
   * Returns a string describing the matrix.
   * 
   * @return 		a description suitable for
   * 			displaying in the experimenter gui
   */
  public String globalInfo() {
    return "Generates the matrix output in LaTeX-syntax.";
  }

  /**
   * returns the name of the output format.
   * 
   * @return		the display name
   */
  public String getDisplayName() {
    return "LaTeX";
  }

  /**
   * removes the stored data but retains the dimensions of the matrix.
   */
  public void clear() {
    super.clear();
    TIE_STRING  = " ";
    WIN_STRING  = "$\\circ$";
    LOSS_STRING = "$\\bullet$";
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
   * returns the matrix as latex table.
   * 
   * @return		the matrix
   */
  public String toStringMatrix() {
    StringBuffer    result;
    String[][]      cells;
    int             i;
    int             j;
    int             n;
    int             size;

    result  = new StringBuffer();
    cells   = toArray();

    result.append(  "\\begin{table}[thb]\n\\caption{\\label{labelname}"
                  + "Table Caption}\n");
    if (!getShowStdDev())
      result.append("\\footnotesize\n");
    else
      result.append("\\scriptsize\n");

    // output the column alignment characters
    // one for the dataset name and one for the comparison column
    if (!getShowStdDev()) {
      result.append(  "{\\centering \\begin{tabular}{"
                    + "l"                     // dataset
                    + ""                      // separator
                    + "r"                     // mean
                    );
    } else {
      // dataset, mean, std dev
      result.append(  "{\\centering \\begin{tabular}{" 
                    + "l"                     // dataset
                    + ""                      // separator
                    + "r"                     // mean
                    + "@{\\hspace{0cm}}"      // separator
                    + "c"                     // +/-
                    + "@{\\hspace{0cm}}"      // separator
                    + "r"                     // stddev
                    );
    }

    for (j = 1; j < getColCount(); j++) {
      if (getColHidden(j))
        continue;
      if (!getShowStdDev())
        result.append(  "r"                   // mean
                      + "@{\\hspace{0.1cm}}"  // separator
                      + "c"                   // significance
                      );
      else 
        result.append(  "r"                   // mean
                      + "@{\\hspace{0cm}}"    // separator
                      + "c"                   // +/-
                      + "@{\\hspace{0cm}}"    // separator
                      + "r"                   // stddev
                      + "@{\\hspace{0.1cm}}"  // separator
                      + "c"                   // significance
                      );
    }
    result.append("}\n\\\\\n\\hline\n");
    if (!getShowStdDev())
      result.append("Dataset & " + cells[0][1]);
    else
      result.append("Dataset & \\multicolumn{3}{c}{" + cells[0][1] + "}");

    // now do the column names (numbers)
    for (j = 2; j < cells[0].length; j++) {
      if (!isMean(j))
        continue;
      if (!getShowStdDev())
        result.append("& " + cells[0][j] + " & ");
      else
        result.append("& \\multicolumn{4}{c}{" + cells[0][j] + "} ");
    }
    result.append("\\\\\n\\hline\n");

    // change "_" to "-" in names
    for (i = 1; i < cells.length; i++)
      cells[i][0] = cells[i][0].replace('_', '-');

    // pad numbers
    for (n = 1; n < cells[0].length; n++) {
      size = getColSize(cells, n);
      for (i = 1; i < cells.length; i++)
        cells[i][n] = padString(cells[i][n], size, true);
    }

    // output data (w/o wins/ties/losses)
    for (i = 1; i < cells.length - 1; i++) {
      if (isAverage(i))
        result.append("\\hline\n");
      for (n = 0; n < cells[0].length; n++) {
        if (n == 0) {
          result.append(padString(cells[i][n], getRowNameWidth()));
        }
        else {
          if (getShowStdDev()) {
            if (isMean(n - 1)) {
              if (!cells[i][n].trim().equals(""))
                result.append(" & $\\pm$ & ");
              else
                result.append(" &       & ");
            }
            else
              result.append(" & ");
          }
          else {
            result.append(" & ");
          }
          result.append(cells[i][n]);
        }
      }
      
      result.append("\\\\\n");
    }

    result.append("\\hline\n\\multicolumn{" + cells[0].length + "}{c}{$\\circ$, $\\bullet$"
		  +" statistically significant improvement or degradation}"
		  +"\\\\\n\\end{tabular} ");
    if (!getShowStdDev())     
      result.append("\\footnotesize ");
    else
      result.append("\\scriptsize ");
    
    result.append("\\par}\n\\end{table}"
		  +"\n");
     
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

    result =   "\\begin{table}[thb]\n\\caption{\\label{labelname}"
             + "Table Caption (Key)}\n";
    result += "\\scriptsize\n";
    result += "{\\centering\n";
    result += "\\begin{tabular}{cl}\\\\\n";
    for (i = 0; i < getColCount(); i++) {
      if (getColHidden(i))
        continue;

      result +=   LEFT_PARENTHESES + (i+1) + RIGHT_PARENTHESES 
                + " & " + removeFilterName(m_ColNames[i]).replace('_', '-')
                                       .replaceAll("\\\\", "\\\\textbackslash") 
                + " \\\\\n";
    }
    result += "\\end{tabular}\n";
    result += "}\n";
    result += "\\end{table}\n";

    return result;
  }

  /**
   * returns the summary as string.
   * 
   * @return		the summary
   */
  public String toStringSummary() {
    int           resultsetLength;
    String        result;
    String        titles;
    int           i;
    int           j;

    if (m_NonSigWins == null)
      return "-summary data not set-";
    
    resultsetLength = 1 + Math.max((int)(Math.log(getColCount())/Math.log(10)),
                                   (int)(Math.log(getRowCount())/Math.log(10)));
    result = "";
    titles = "";

    result += "{\\centering\n";
    result += "\\begin{table}[thb]\n\\caption{\\label{labelname}"
                +"Table Caption}\n";
    result += "\\footnotesize\n";
    result += "\\begin{tabular}{l";

    for (i = 0; i < getColCount(); i++) {
      if (getColHidden(i))
        continue;

      titles += " &";
      result += "c";
      titles += ' ' + Utils.padLeft("" + getSummaryTitle(i),
				    resultsetLength * 2 + 3);
    }
    result += "}\\\\\n\\hline\n";
    result += titles + " \\\\\n\\hline\n";
    
    for (i = 0; i < getColCount(); i++) {
      if (getColHidden(i))
        continue;

      for (j = 0; j < getColCount(); j++) {
        if (getColHidden(j))
          continue;

	if (j == 0)
	  result +=  (char)((int)'a' + i % 26);

	if (j == i)
	  result += " & - ";
	else
	  result += "& " + m_NonSigWins[i][j] + " (" + m_Wins[i][j] + ") ";
      }
      result += "\\\\\n";
    }

    result += "\\hline\n\\end{tabular} \\footnotesize \\par\n\\end{table}}";

    return result;
  }

  /**
   * returns the ranking in a string representation.
   * 
   * @return 		the ranking
   */
  public String toStringRanking() {
    int           biggest;
    int           width;
    String        result;
    int[]         ranking;
    int           i;
    int           curr;

    if (m_RankingWins == null)
      return "-ranking data not set-";

    biggest = Math.max(m_RankingWins[Utils.maxIndex(m_RankingWins)],
                       m_RankingLosses[Utils.maxIndex(m_RankingLosses)]);
    width = Math.max(2 + (int)(Math.log(biggest) / Math.log(10)),
			 ">-<".length());
    result = "\\begin{table}[thb]\n\\caption{\\label{labelname}Table Caption"
      + "}\n\\footnotesize\n{\\centering \\begin{tabular}{rlll}\\\\\n\\hline\n";
    result +=   "Resultset & Wins$-$ & Wins & Losses \\\\\n& Losses & & "
              + "\\\\\n\\hline\n";

    ranking = Utils.sort(m_RankingDiff);
    for (i = getColCount() - 1; i >= 0; i--) {
      curr = ranking[i];
      
      if (getColHidden(curr))
        continue;

      result +=   "(" + (curr + 1) + ") & " 
                + Utils.padLeft("" + m_RankingDiff[curr], width) 
                + " & " + Utils.padLeft("" + m_RankingWins[curr], width)
                + " & " + Utils.padLeft("" + m_RankingLosses[curr], width)
                + "\\\\\n";
    }

    result += "\\hline\n\\end{tabular} \\footnotesize \\par}\n\\end{table}";

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
    
    matrix = new ResultMatrixLatex(3, 3);
    
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
