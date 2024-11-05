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
 *    CfsSubsetEval.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package  weka.attributeSelection;

import weka.core.Capabilities;
import weka.core.ContingencyTables;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformationHandler;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;
import weka.filters.Filter;
import weka.filters.supervised.attribute.Discretize;

import java.util.BitSet;
import java.util.Enumeration;
import java.util.Vector;

/** 
 <!-- globalinfo-start -->
 * CfsSubsetEval :<br/>
 * <br/>
 * Evaluates the worth of a subset of attributes by considering the individual predictive ability of each feature along with the degree of redundancy between them.<br/>
 * <br/>
 * Subsets of features that are highly correlated with the class while having low intercorrelation are preferred.<br/>
 * <br/>
 * For more information see:<br/>
 * <br/>
 * M. A. Hall (1998). Correlation-based Feature Subset Selection for Machine Learning. Hamilton, New Zealand.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;phdthesis{Hall1998,
 *    address = {Hamilton, New Zealand},
 *    author = {M. A. Hall},
 *    school = {University of Waikato},
 *    title = {Correlation-based Feature Subset Selection for Machine Learning},
 *    year = {1998}
 * }
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -M
 *  Treat missing values as a separate value.</pre>
 * 
 * <pre> -L
 *  Don't include locally predictive attributes.</pre>
 * 
 <!-- options-end -->
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @version $Revision: 6132 $
 * @see Discretize
 */
public class CfsSubsetEval
  extends ASEvaluation
  implements SubsetEvaluator, 
             OptionHandler, 
             TechnicalInformationHandler {
  
  /** for serialization */
  static final long serialVersionUID = 747878400813276317L;

  /** The training instances */
  private Instances m_trainInstances;
  /** Discretise attributes when class in nominal */
  private Discretize m_disTransform;
  /** The class index */
  private int m_classIndex;
  /** Is the class numeric */
  private boolean m_isNumeric;
  /** Number of attributes in the training data */
  private int m_numAttribs;
  /** Number of instances in the training data */
  private int m_numInstances;
  /** Treat missing values as separate values */
  private boolean m_missingSeparate;
  /** Include locally predictive attributes */
  private boolean m_locallyPredictive;
  /** Holds the matrix of attribute correlations */
  //  private Matrix m_corr_matrix;
  private float [][] m_corr_matrix;
  /** Standard deviations of attributes (when using pearsons correlation) */
  private double[] m_std_devs;
  /** Threshold for admitting locally predictive features */
  private double m_c_Threshold;

  /**
   * Returns a string describing this attribute evaluator
   * @return a description of the evaluator suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return "CfsSubsetEval :\n\nEvaluates the worth of a subset of attributes "
      +"by considering the individual predictive ability of each feature "
      +"along with the degree of redundancy between them.\n\n"
      +"Subsets of features that are highly correlated with the class "
      +"while having low intercorrelation are preferred.\n\n"
      + "For more information see:\n\n"
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
    TechnicalInformation        result;
    
    result = new TechnicalInformation(Type.PHDTHESIS);
    result.setValue(Field.AUTHOR, "M. A. Hall");
    result.setValue(Field.YEAR, "1998");
    result.setValue(Field.TITLE, "Correlation-based Feature Subset Selection for Machine Learning");
    result.setValue(Field.SCHOOL, "University of Waikato");
    result.setValue(Field.ADDRESS, "Hamilton, New Zealand");
    
    return result;
  }

  /**
   * Constructor
   */
  public CfsSubsetEval () {
    resetOptions();
  }


  /**
   * Returns an enumeration describing the available options.
   * @return an enumeration of all the available options.
   *
   **/
  public Enumeration listOptions () {
    Vector newVector = new Vector(3);
    newVector.addElement(new Option("\tTreat missing values as a separate " 
                                    + "value.", "M", 0, "-M"));
    newVector.addElement(new Option("\tDon't include locally predictive attributes" 
                                    + ".", "L", 0, "-L"));
    return  newVector.elements();
  }


  /**
   * Parses and sets a given list of options. <p/>
   *
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -M
   *  Treat missing values as a separate value.</pre>
   * 
   * <pre> -L
   *  Don't include locally predictive attributes.</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   *
   **/
  public void setOptions (String[] options)
    throws Exception {

    resetOptions();
    setMissingSeparate(Utils.getFlag('M', options));
    setLocallyPredictive(!Utils.getFlag('L', options));
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String locallyPredictiveTipText() {
    return "Identify locally predictive attributes. Iteratively adds "
      +"attributes with the highest correlation with the class as long "
      +"as there is not already an attribute in the subset that has a "
      +"higher correlation with the attribute in question";
  }

  /**
   * Include locally predictive attributes
   *
   * @param b true or false
   */
  public void setLocallyPredictive (boolean b) {
    m_locallyPredictive = b;
  }


  /**
   * Return true if including locally predictive attributes
   *
   * @return true if locally predictive attributes are to be used
   */
  public boolean getLocallyPredictive () {
    return  m_locallyPredictive;
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String missingSeparateTipText() {
    return "Treat missing as a separate value. Otherwise, counts for missing "
      +"values are distributed across other values in proportion to their "
      +"frequency.";
  }

  /**
   * Treat missing as a separate value
   *
   * @param b true or false
   */
  public void setMissingSeparate (boolean b) {
    m_missingSeparate = b;
  }


  /**
   * Return true is missing is treated as a separate value
   *
   * @return true if missing is to be treated as a separate value
   */
  public boolean getMissingSeparate () {
    return  m_missingSeparate;
  }


  /**
   * Gets the current settings of CfsSubsetEval
   *
   * @return an array of strings suitable for passing to setOptions()
   */
  public String[] getOptions () {
    String[] options = new String[2];
    int current = 0;

    if (getMissingSeparate()) {
      options[current++] = "-M";
    }

    if (!getLocallyPredictive()) {
      options[current++] = "-L";
    }

    while (current < options.length) {
      options[current++] = "";
    }

    return  options;
  }

  /**
   * Returns the capabilities of this evaluator.
   *
   * @return            the capabilities of this evaluator
   * @see               Capabilities
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
    
    //RANKING BEGIN
    result.enable(Capability.PREFERENCE_ATTRIBUTE);
    result.enable(Capability.RANKING);
    //RANKING END
    
    return result;
  }

  /**
   * Generates a attribute evaluator. Has to initialize all fields of the 
   * evaluator that are not being set via options.
   *
   * CFS also discretises attributes (if necessary) and initializes
   * the correlation matrix.
   *
   * @param data set of instances serving as training data 
   * @throws Exception if the evaluator has not been 
   * generated successfully
   */
  public void buildEvaluator (Instances data)
    throws Exception {

    // can evaluator handle data?
    getCapabilities().testWithFail(data);

    m_trainInstances = new Instances(data);
    m_trainInstances.deleteWithMissingClass();
    m_classIndex = m_trainInstances.classIndex();
    m_numAttribs = m_trainInstances.numAttributes();
    m_numInstances = m_trainInstances.numInstances();
    m_isNumeric = m_trainInstances.attribute(m_classIndex).isNumeric();

    if (!m_isNumeric) {
      m_disTransform = new Discretize();
      m_disTransform.setUseBetterEncoding(true);
      m_disTransform.setInputFormat(m_trainInstances);
      m_trainInstances = Filter.useFilter(m_trainInstances, m_disTransform);
    }

    m_std_devs = new double[m_numAttribs];
    m_corr_matrix = new float [m_numAttribs][];
    for (int i = 0; i < m_numAttribs; i++) {
      m_corr_matrix[i] = new float [i+1];
    }

    for (int i = 0; i < m_corr_matrix.length; i++) {
      m_corr_matrix[i][i] = 1.0f;
      m_std_devs[i] = 1.0;
    }

    for (int i = 0; i < m_numAttribs; i++) {
      for (int j = 0; j < m_corr_matrix[i].length - 1; j++) {
        m_corr_matrix[i][j] = -999;
      }
    }
  }


  /**
   * evaluates a subset of attributes
   *
   * @param subset a bitset representing the attribute subset to be 
   * evaluated 
   * @return the merit
   * @throws Exception if the subset could not be evaluated
   */
  public double evaluateSubset (BitSet subset)
    throws Exception {
    double num = 0.0;
    double denom = 0.0;
    float corr;
    int larger, smaller;
    // do numerator
    for (int i = 0; i < m_numAttribs; i++) {
      if (i != m_classIndex) {
        if (subset.get(i)) {
          if (i > m_classIndex) {
            larger = i; smaller = m_classIndex;
          } else {
            smaller = i; larger = m_classIndex;
          }
          /*      int larger = (i > m_classIndex ? i : m_classIndex);
                  int smaller = (i > m_classIndex ? m_classIndex : i); */
          if (m_corr_matrix[larger][smaller] == -999) {
            corr = correlate(i, m_classIndex);
            m_corr_matrix[larger][smaller] = corr;
            num += (m_std_devs[i] * corr);
          }
          else {
            num += (m_std_devs[i] * m_corr_matrix[larger][smaller]);
          }
        }
      }
    }

    // do denominator
    for (int i = 0; i < m_numAttribs; i++) {
      if (i != m_classIndex) {
        if (subset.get(i)) {
          denom += (1.0 * m_std_devs[i] * m_std_devs[i]);

          for (int j = 0; j < m_corr_matrix[i].length - 1; j++) {
            if (subset.get(j)) {
              if (m_corr_matrix[i][j] == -999) {
                corr = correlate(i, j);
                m_corr_matrix[i][j] = corr;
                denom += (2.0 * m_std_devs[i] * m_std_devs[j] * corr);
              }
              else {
                denom += (2.0 * m_std_devs[i] * m_std_devs[j] * m_corr_matrix[i][j]);
              }
            }
          }
        }
      }
    }

    if (denom < 0.0) {
      denom *= -1.0;
    }

    if (denom == 0.0) {
      return  (0.0);
    }

    double merit = (num/Math.sqrt(denom));

    if (merit < 0.0) {
      merit *= -1.0;
    }

    return  merit;
  }


  private float correlate (int att1, int att2) {
    if (!m_isNumeric) {
      return  (float) symmUncertCorr(att1, att2);
    }

    boolean att1_is_num = (m_trainInstances.attribute(att1).isNumeric());
    boolean att2_is_num = (m_trainInstances.attribute(att2).isNumeric());

    if (att1_is_num && att2_is_num) {
      return  (float) num_num(att1, att2);
    }
    else {if (att2_is_num) {
      return  (float) num_nom2(att1, att2);
    }
    else {if (att1_is_num) {
      return  (float) num_nom2(att2, att1);
    }
    }
    }

    return (float) nom_nom(att1, att2);
  }


  private double symmUncertCorr (int att1, int att2) {
    int i, j, k, ii, jj;
    int ni, nj;
    double sum = 0.0;
    double sumi[], sumj[];
    double counts[][];
    Instance inst;
    double corr_measure;
    boolean flag = false;
    double temp = 0.0;

    if (att1 == m_classIndex || att2 == m_classIndex) {
      flag = true;
    }

    ni = m_trainInstances.attribute(att1).numValues() + 1;
    nj = m_trainInstances.attribute(att2).numValues() + 1;
    counts = new double[ni][nj];
    sumi = new double[ni];
    sumj = new double[nj];

    for (i = 0; i < ni; i++) {
      sumi[i] = 0.0;

      for (j = 0; j < nj; j++) {
        sumj[j] = 0.0;
        counts[i][j] = 0.0;
      }
    }

    // Fill the contingency table
    for (i = 0; i < m_numInstances; i++) {
      inst = m_trainInstances.instance(i);

      if (inst.isMissing(att1)) {
        ii = ni - 1;
      }
      else {
        ii = (int)inst.value(att1);
      }

      if (inst.isMissing(att2)) {
        jj = nj - 1;
      }
      else {
        jj = (int)inst.value(att2);
      }

      counts[ii][jj]++;
    }

    // get the row totals
    for (i = 0; i < ni; i++) {
      sumi[i] = 0.0;

      for (j = 0; j < nj; j++) {
        sumi[i] += counts[i][j];
        sum += counts[i][j];
      }
    }

    // get the column totals
    for (j = 0; j < nj; j++) {
      sumj[j] = 0.0;

      for (i = 0; i < ni; i++) {
        sumj[j] += counts[i][j];
      }
    }

    // distribute missing counts
    if (!m_missingSeparate && 
        (sumi[ni-1] < m_numInstances) && 
        (sumj[nj-1] < m_numInstances)) {
      double[] i_copy = new double[sumi.length];
      double[] j_copy = new double[sumj.length];
      double[][] counts_copy = new double[sumi.length][sumj.length];

      for (i = 0; i < ni; i++) {
        System.arraycopy(counts[i], 0, counts_copy[i], 0, sumj.length);
      }

      System.arraycopy(sumi, 0, i_copy, 0, sumi.length);
      System.arraycopy(sumj, 0, j_copy, 0, sumj.length);
      double total_missing = 
        (sumi[ni - 1] + sumj[nj - 1] - counts[ni - 1][nj - 1]);

      // do the missing i's
      if (sumi[ni - 1] > 0.0) {
        for (j = 0; j < nj - 1; j++) {
          if (counts[ni - 1][j] > 0.0) {
            for (i = 0; i < ni - 1; i++) {
              temp = ((i_copy[i]/(sum - i_copy[ni - 1]))*counts[ni - 1][j]);
              counts[i][j] += temp;
              sumi[i] += temp;
            }

            counts[ni - 1][j] = 0.0;
          }
        }
      }

      sumi[ni - 1] = 0.0;

      // do the missing j's
      if (sumj[nj - 1] > 0.0) {
        for (i = 0; i < ni - 1; i++) {
          if (counts[i][nj - 1] > 0.0) {
            for (j = 0; j < nj - 1; j++) {
              temp = ((j_copy[j]/(sum - j_copy[nj - 1]))*counts[i][nj - 1]);
              counts[i][j] += temp;
              sumj[j] += temp;
            }

            counts[i][nj - 1] = 0.0;
          }
        }
      }

      sumj[nj - 1] = 0.0;

      // do the both missing
      if (counts[ni - 1][nj - 1] > 0.0 && total_missing != sum) {
        for (i = 0; i < ni - 1; i++) {
          for (j = 0; j < nj - 1; j++) {
            temp = (counts_copy[i][j]/(sum - total_missing)) * 
              counts_copy[ni - 1][nj - 1];
            
            counts[i][j] += temp;
            sumi[i] += temp;
            sumj[j] += temp;
          }
        }

        counts[ni - 1][nj - 1] = 0.0;
      }
    }

    corr_measure = ContingencyTables.symmetricalUncertainty(counts);

    if (Utils.eq(corr_measure, 0.0)) {
      if (flag == true) {
        return  (0.0);
      }
      else {
        return  (1.0);
      }
    }
    else {
      return  (corr_measure);
    }
  }


  private double num_num (int att1, int att2) {
    int i;
    Instance inst;
    double r, diff1, diff2, num = 0.0, sx = 0.0, sy = 0.0;
    double mx = m_trainInstances.meanOrMode(m_trainInstances.attribute(att1));
    double my = m_trainInstances.meanOrMode(m_trainInstances.attribute(att2));

    for (i = 0; i < m_numInstances; i++) {
      inst = m_trainInstances.instance(i);
      diff1 = (inst.isMissing(att1))? 0.0 : (inst.value(att1) - mx);
      diff2 = (inst.isMissing(att2))? 0.0 : (inst.value(att2) - my);
      num += (diff1*diff2);
      sx += (diff1*diff1);
      sy += (diff2*diff2);
    }

    if (sx != 0.0) {
      if (m_std_devs[att1] == 1.0) {
        m_std_devs[att1] = Math.sqrt((sx/m_numInstances));
      }
    }

    if (sy != 0.0) {
      if (m_std_devs[att2] == 1.0) {
        m_std_devs[att2] = Math.sqrt((sy/m_numInstances));
      }
    }

    if ((sx*sy) > 0.0) {
      r = (num/(Math.sqrt(sx*sy)));
      return  ((r < 0.0)? -r : r);
    }
    else {
      if (att1 != m_classIndex && att2 != m_classIndex) {
        return  1.0;
      }
      else {
        return  0.0;
      }
    }
  }


  private double num_nom2 (int att1, int att2) {
    int i, ii, k;
    double temp;
    Instance inst;
    int mx = (int)m_trainInstances.
      meanOrMode(m_trainInstances.attribute(att1));
    double my = m_trainInstances.
      meanOrMode(m_trainInstances.attribute(att2));
    double stdv_num = 0.0;
    double diff1, diff2;
    double r = 0.0, rr;
    int nx = (!m_missingSeparate) 
      ? m_trainInstances.attribute(att1).numValues() 
      : m_trainInstances.attribute(att1).numValues() + 1;

    double[] prior_nom = new double[nx];
    double[] stdvs_nom = new double[nx];
    double[] covs = new double[nx];

    for (i = 0; i < nx; i++) {
      stdvs_nom[i] = covs[i] = prior_nom[i] = 0.0;
    }

    // calculate frequencies (and means) of the values of the nominal 
    // attribute
    for (i = 0; i < m_numInstances; i++) {
      inst = m_trainInstances.instance(i);

      if (inst.isMissing(att1)) {
        if (!m_missingSeparate) {
          ii = mx;
        }
        else {
          ii = nx - 1;
        }
      }
      else {
        ii = (int)inst.value(att1);
      }

      // increment freq for nominal
      prior_nom[ii]++;
    }

    for (k = 0; k < m_numInstances; k++) {
      inst = m_trainInstances.instance(k);
      // std dev of numeric attribute
      diff2 = (inst.isMissing(att2))? 0.0 : (inst.value(att2) - my);
      stdv_num += (diff2*diff2);

      // 
      for (i = 0; i < nx; i++) {
        if (inst.isMissing(att1)) {
          if (!m_missingSeparate) {
            temp = (i == mx)? 1.0 : 0.0;
          }
          else {
            temp = (i == (nx - 1))? 1.0 : 0.0;
          }
        }
        else {
          temp = (i == inst.value(att1))? 1.0 : 0.0;
        }

        diff1 = (temp - (prior_nom[i]/m_numInstances));
        stdvs_nom[i] += (diff1*diff1);
        covs[i] += (diff1*diff2);
      }
    }

    // calculate weighted correlation
    for (i = 0, temp = 0.0; i < nx; i++) {
      // calculate the weighted variance of the nominal
      temp += ((prior_nom[i]/m_numInstances)*(stdvs_nom[i]/m_numInstances));

      if ((stdvs_nom[i]*stdv_num) > 0.0) {
        //System.out.println("Stdv :"+stdvs_nom[i]);
        rr = (covs[i]/(Math.sqrt(stdvs_nom[i]*stdv_num)));

        if (rr < 0.0) {
          rr = -rr;
        }

        r += ((prior_nom[i]/m_numInstances)*rr);
      }
      /* if there is zero variance for the numeric att at a specific 
         level of the catergorical att then if neither is the class then 
         make this correlation at this level maximally bad i.e. 1.0. 
         If either is the class then maximally bad correlation is 0.0 */
      else {if (att1 != m_classIndex && att2 != m_classIndex) {
        r += ((prior_nom[i]/m_numInstances)*1.0);
      }
      }
    }

    // set the standard deviations for these attributes if necessary
    // if ((att1 != classIndex) && (att2 != classIndex)) // =============
    if (temp != 0.0) {
      if (m_std_devs[att1] == 1.0) {
        m_std_devs[att1] = Math.sqrt(temp);
      }
    }

    if (stdv_num != 0.0) {
      if (m_std_devs[att2] == 1.0) {
        m_std_devs[att2] = Math.sqrt((stdv_num/m_numInstances));
      }
    }

    if (r == 0.0) {
      if (att1 != m_classIndex && att2 != m_classIndex) {
        r = 1.0;
      }
    }

    return  r;
  }


  private double nom_nom (int att1, int att2) {
    int i, j, ii, jj, z;
    double temp1, temp2;
    Instance inst;
    int mx = (int)m_trainInstances.
      meanOrMode(m_trainInstances.attribute(att1));
    int my = (int)m_trainInstances.
      meanOrMode(m_trainInstances.attribute(att2));
    double diff1, diff2;
    double r = 0.0, rr;
    int nx = (!m_missingSeparate) 
      ? m_trainInstances.attribute(att1).numValues() 
      : m_trainInstances.attribute(att1).numValues() + 1;

    int ny = (!m_missingSeparate)
      ? m_trainInstances.attribute(att2).numValues() 
      : m_trainInstances.attribute(att2).numValues() + 1;

    double[][] prior_nom = new double[nx][ny];
    double[] sumx = new double[nx];
    double[] sumy = new double[ny];
    double[] stdvsx = new double[nx];
    double[] stdvsy = new double[ny];
    double[][] covs = new double[nx][ny];

    for (i = 0; i < nx; i++) {
      sumx[i] = stdvsx[i] = 0.0;
    }

    for (j = 0; j < ny; j++) {
      sumy[j] = stdvsy[j] = 0.0;
    }

    for (i = 0; i < nx; i++) {
      for (j = 0; j < ny; j++) {
        covs[i][j] = prior_nom[i][j] = 0.0;
      }
    }

    // calculate frequencies (and means) of the values of the nominal 
    // attribute
    for (i = 0; i < m_numInstances; i++) {
      inst = m_trainInstances.instance(i);

      if (inst.isMissing(att1)) {
        if (!m_missingSeparate) {
          ii = mx;
        }
        else {
          ii = nx - 1;
        }
      }
      else {
        ii = (int)inst.value(att1);
      }

      if (inst.isMissing(att2)) {
        if (!m_missingSeparate) {
          jj = my;
        }
        else {
          jj = ny - 1;
        }
      }
      else {
        jj = (int)inst.value(att2);
      }

      // increment freq for nominal
      prior_nom[ii][jj]++;
      sumx[ii]++;
      sumy[jj]++;
    }

    for (z = 0; z < m_numInstances; z++) {
      inst = m_trainInstances.instance(z);

      for (j = 0; j < ny; j++) {
        if (inst.isMissing(att2)) {
          if (!m_missingSeparate) {
            temp2 = (j == my)? 1.0 : 0.0;
          }
          else {
            temp2 = (j == (ny - 1))? 1.0 : 0.0;
          }
        }
        else {
          temp2 = (j == inst.value(att2))? 1.0 : 0.0;
        }

        diff2 = (temp2 - (sumy[j]/m_numInstances));
        stdvsy[j] += (diff2*diff2);
      }

      // 
      for (i = 0; i < nx; i++) {
        if (inst.isMissing(att1)) {
          if (!m_missingSeparate) {
            temp1 = (i == mx)? 1.0 : 0.0;
          }
          else {
            temp1 = (i == (nx - 1))? 1.0 : 0.0;
          }
        }
        else {
          temp1 = (i == inst.value(att1))? 1.0 : 0.0;
        }

        diff1 = (temp1 - (sumx[i]/m_numInstances));
        stdvsx[i] += (diff1*diff1);

        for (j = 0; j < ny; j++) {
          if (inst.isMissing(att2)) {
            if (!m_missingSeparate) {
              temp2 = (j == my)? 1.0 : 0.0;
            }
            else {
              temp2 = (j == (ny - 1))? 1.0 : 0.0;
            }
          }
          else {
            temp2 = (j == inst.value(att2))? 1.0 : 0.0;
          }

          diff2 = (temp2 - (sumy[j]/m_numInstances));
          covs[i][j] += (diff1*diff2);
        }
      }
    }

    // calculate weighted correlation
    for (i = 0; i < nx; i++) {
      for (j = 0; j < ny; j++) {
        if ((stdvsx[i]*stdvsy[j]) > 0.0) {
          //System.out.println("Stdv :"+stdvs_nom[i]);
          rr = (covs[i][j]/(Math.sqrt(stdvsx[i]*stdvsy[j])));

          if (rr < 0.0) {
            rr = -rr;
          }

          r += ((prior_nom[i][j]/m_numInstances)*rr);
        }
        // if there is zero variance for either of the categorical atts then if
        // neither is the class then make this
        // correlation at this level maximally bad i.e. 1.0. If either is 
        // the class then maximally bad correlation is 0.0
        else {if (att1 != m_classIndex && att2 != m_classIndex) {
          r += ((prior_nom[i][j]/m_numInstances)*1.0);
        }
        }
      }
    }

    // calculate weighted standard deviations for these attributes
    // (if necessary)
    for (i = 0, temp1 = 0.0; i < nx; i++) {
      temp1 += ((sumx[i]/m_numInstances)*(stdvsx[i]/m_numInstances));
    }

    if (temp1 != 0.0) {
      if (m_std_devs[att1] == 1.0) {
        m_std_devs[att1] = Math.sqrt(temp1);
      }
    }

    for (j = 0, temp2 = 0.0; j < ny; j++) {
      temp2 += ((sumy[j]/m_numInstances)*(stdvsy[j]/m_numInstances));
    }

    if (temp2 != 0.0) {
      if (m_std_devs[att2] == 1.0) {
        m_std_devs[att2] = Math.sqrt(temp2);
      }
    }

    if (r == 0.0) {
      if (att1 != m_classIndex && att2 != m_classIndex) {
        r = 1.0;
      }
    }

    return  r;
  }


  /**
   * returns a string describing CFS
   *
   * @return the description as a string
   */
  public String toString () {
    StringBuffer text = new StringBuffer();

    if (m_trainInstances == null) {
      text.append("CFS subset evaluator has not been built yet\n");
    }
    else {
      text.append("\tCFS Subset Evaluator\n");

      if (m_missingSeparate) {
        text.append("\tTreating missing values as a separate value\n");
      }

      if (m_locallyPredictive) {
        text.append("\tIncluding locally predictive attributes\n");
      }
    }

    return  text.toString();
  }


  private void addLocallyPredictive (BitSet best_group) {
    int i, j;
    boolean done = false;
    boolean ok = true;
    double temp_best = -1.0;
    float corr;
    j = 0;
    BitSet temp_group = (BitSet)best_group.clone();
    int larger, smaller;

    while (!done) {
      temp_best = -1.0;

      // find best not already in group
      for (i = 0; i < m_numAttribs; i++) {
        if (i > m_classIndex) {
          larger = i; smaller = m_classIndex;
        } else {
          smaller = i; larger = m_classIndex;
        }
        /*      int larger = (i > m_classIndex ? i : m_classIndex);
                int smaller = (i > m_classIndex ? m_classIndex : i); */
        if ((!temp_group.get(i)) && (i != m_classIndex)) {
          if (m_corr_matrix[larger][smaller] == -999) {
            corr = correlate(i, m_classIndex);
            m_corr_matrix[larger][smaller] = corr;
          }

          if (m_corr_matrix[larger][smaller]  > temp_best) {
            temp_best = m_corr_matrix[larger][smaller];
            j = i;
          }
        }
      }

      if (temp_best == -1.0) {
        done = true;
      }
      else {
        ok = true;
        temp_group.set(j);

        // check the best against correlations with others already
        // in group 
        for (i = 0; i < m_numAttribs; i++) {
          if (i > j) {
            larger = i; smaller = j;
          } else {
            larger = j; smaller = i;
          }
          /*  int larger = (i > j ? i : j);
              int smaller = (i > j ? j : i); */
          if (best_group.get(i)) {
            if (m_corr_matrix[larger][smaller] == -999) {
              corr = correlate(i, j);
              m_corr_matrix[larger][smaller] = corr;
            }

            if (m_corr_matrix[larger][smaller] > temp_best - m_c_Threshold) {
              ok = false;
              break;
            }
          }
        }

        // if ok then add to best_group
        if (ok) {
          best_group.set(j);
        }
      }
    }
  }


  /**
   * Calls locallyPredictive in order to include locally predictive
   * attributes (if requested).
   *
   * @param attributeSet the set of attributes found by the search
   * @return a possibly ranked list of postprocessed attributes
   * @throws Exception if postprocessing fails for some reason
   */
  public int[] postProcess (int[] attributeSet)
    throws Exception {
    int j = 0;

    if (!m_locallyPredictive) {
      //      m_trainInstances = new Instances(m_trainInstances,0);
      return  attributeSet;
    }

    BitSet bestGroup = new BitSet(m_numAttribs);

    for (int i = 0; i < attributeSet.length; i++) {
      bestGroup.set(attributeSet[i]);
    }

    addLocallyPredictive(bestGroup);

    // count how many are set
    for (int i = 0; i < m_numAttribs; i++) {
      if (bestGroup.get(i)) {
        j++;
      }
    }

    int[] newSet = new int[j];
    j = 0;

    for (int i = 0; i < m_numAttribs; i++) {
      if (bestGroup.get(i)) {
        newSet[j++] = i;
      }
    }

    //    m_trainInstances = new Instances(m_trainInstances,0);
    return  newSet;
  }


  protected void resetOptions () {
    m_trainInstances = null;
    m_missingSeparate = false;
    m_locallyPredictive = true;
    m_c_Threshold = 0.0;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6132 $");
  }

  /**
   * Main method for testing this class.
   *
   * @param args the options
   */
  public static void main (String[] args) {
    runEvaluator(new CfsSubsetEval(), args);
  }
}
