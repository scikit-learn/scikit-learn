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
 *    PairedCorrectedTTester.java
 *    Copyright (C) 2003 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.experiment;

import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformation.Type;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformationHandler;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Enumeration;

/**
 * Behaves the same as PairedTTester, only it uses the corrected
 * resampled t-test statistic.<p/>
 *
 * For more information see:<p/>
 *
 <!-- technical-plaintext-start -->
 * Claude Nadeau, Yoshua Bengio (2001). Inference for the Generalization Error. Machine Learning..
 <!-- technical-plaintext-end -->
 *
 * <p/>
 * 
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;article{Nadeau2001,
 *    author = {Claude Nadeau and Yoshua Bengio},
 *    journal = {Machine Learning},
 *    title = {Inference for the Generalization Error},
 *    year = {2001},
 *    PDF = {http://www.iro.umontreal.ca/\~lisa/bib/pub_subject/comparative/pointeurs/nadeau_MLJ1597.pdf}
 * }
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
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
 * @author Richard Kirkby (rkirkby@cs.waikato.ac.nz)
 * @version $Revision: 1.13 $
 */
public class PairedCorrectedTTester 
  extends PairedTTester
  implements TechnicalInformationHandler {
  
  /** for serialization */
  static final long serialVersionUID = -3105268939845653323L;

  /**
   * Returns an instance of a TechnicalInformation object, containing 
   * detailed information about the technical background of this class,
   * e.g., paper reference or book this class is based on.
   * 
   * @return the technical information about this class
   */
  public TechnicalInformation getTechnicalInformation() {
    TechnicalInformation 	result;
    
    result = new TechnicalInformation(Type.ARTICLE);
    result.setValue(Field.AUTHOR, "Claude Nadeau and Yoshua Bengio");
    result.setValue(Field.YEAR, "2001");
    result.setValue(Field.TITLE, "Inference for the Generalization Error");
    result.setValue(Field.JOURNAL, "Machine Learning");
    result.setValue(Field.PDF, "http://www.iro.umontreal.ca/~lisa/bib/pub_subject/comparative/pointeurs/nadeau_MLJ1597.pdf");

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

    // calculate the test/train ratio
    double testTrainRatio = 0.0;
    int trainSizeIndex = -1;
    int testSizeIndex = -1;
    // find the columns with the train/test sizes
    for (int i=0; i<m_Instances.numAttributes(); i++) {
      if (m_Instances.attribute(i).name().toLowerCase().equals("number_of_training_instances")) {
	trainSizeIndex = i;
      } else if (m_Instances.attribute(i).name().toLowerCase().equals("number_of_testing_instances")) {
	testSizeIndex = i;
      }
    }
    if (trainSizeIndex >= 0 && testSizeIndex >= 0) {
      double totalTrainSize = 0.0;
      double totalTestSize = 0.0;
      for (int k = 0; k < dataset1.size(); k ++) {
	Instance current = (Instance) dataset1.elementAt(k);
	totalTrainSize += current.value(trainSizeIndex);
	totalTestSize += current.value(testSizeIndex);
      }
      testTrainRatio = totalTestSize / totalTrainSize;
    }
    PairedStats pairedStats =
      new PairedStatsCorrected(m_SignificanceLevel, testTrainRatio);

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
    return pairedStats;
  }

  /**
   * Test the class from the command line.
   *
   * @param args contains options for the instance ttests
   */
  public static void main(String args[]) {
    
    try {
      PairedCorrectedTTester tt = new PairedCorrectedTTester();
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

  /**
   * returns the name of the tester
   * 
   * @return the display name
   */
  public String getDisplayName() {
    return "Paired T-Tester (corrected)";
  }

  /**
   * returns a string that is displayed as tooltip on the "perform test"
   * button in the experimenter
   * 
   * @return the string for the tool tip
   */
  public String getToolTipText() {
    return "Performs test using corrected resampled t-test statistic (Nadeau and Bengio)";
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 1.13 $");
  }
}
