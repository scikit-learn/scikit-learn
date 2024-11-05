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
 *    Ranker.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.attributeSelection;

import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.Range;
import weka.core.RevisionUtils;
import weka.core.Utils;

import java.util.Enumeration;
import java.util.Vector;

/** 
 <!-- globalinfo-start -->
 * Ranker : <br/>
 * <br/>
 * Ranks attributes by their individual evaluations. Use in conjunction with attribute evaluators (ReliefF, GainRatio, Entropy etc).<br/>
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -P &lt;start set&gt;
 *  Specify a starting set of attributes.
 *  Eg. 1,3,5-7.
 *  Any starting attributes specified are
 *  ignored during the ranking.</pre>
 * 
 * <pre> -T &lt;threshold&gt;
 *  Specify a theshold by which attributes
 *  may be discarded from the ranking.</pre>
 * 
 * <pre> -N &lt;num to select&gt;
 *  Specify number of attributes to select</pre>
 * 
 <!-- options-end -->
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @version $Revision: 1.26 $
 */
public class Ranker 
  extends ASSearch 
  implements RankedOutputSearch, StartSetHandler, OptionHandler {
  
  /** for serialization */
  static final long serialVersionUID = -9086714848510751934L;

  /** Holds the starting set as an array of attributes */
  private int[] m_starting;

  /** Holds the start set for the search as a range */
  private Range m_startRange;

  /** Holds the ordered list of attributes */
  private int[] m_attributeList;

  /** Holds the list of attribute merit scores */
  private double[] m_attributeMerit;

  /** Data has class attribute---if unsupervised evaluator then no class */
  private boolean m_hasClass;

  /** Class index of the data if supervised evaluator */
  private int m_classIndex;

  /** The number of attribtes */
  private int m_numAttribs;

  /** 
   * A threshold by which to discard attributes---used by the
   * AttributeSelection module
   */
  private double m_threshold;

  /** The number of attributes to select. -1 indicates that all attributes
      are to be retained. Has precedence over m_threshold */
  private int m_numToSelect = -1;

  /** Used to compute the number to select */
  private int m_calculatedNumToSelect = -1;

  /**
   * Returns a string describing this search method
   * @return a description of the search suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return "Ranker : \n\nRanks attributes by their individual evaluations. "
      +"Use in conjunction with attribute evaluators (ReliefF, GainRatio, "
      +"Entropy etc).\n";
  }

  /**
   * Constructor
   */
  public Ranker () {
    resetOptions();
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String numToSelectTipText() {
    return "Specify the number of attributes to retain. The default value "
      +"(-1) indicates that all attributes are to be retained. Use either "
      +"this option or a threshold to reduce the attribute set.";
  }

  /**
   * Specify the number of attributes to select from the ranked list. -1
   * indicates that all attributes are to be retained.
   * @param n the number of attributes to retain
   */
  public void setNumToSelect(int n) {
    m_numToSelect = n;
  }

  /**
   * Gets the number of attributes to be retained.
   * @return the number of attributes to retain
   */
  public int getNumToSelect() {
    return m_numToSelect;
  }

  /**
   * Gets the calculated number to select. This might be computed
   * from a threshold, or if < 0 is set as the number to select then
   * it is set to the number of attributes in the (transformed) data.
   * @return the calculated number of attributes to select
   */
  public int getCalculatedNumToSelect() {
    if (m_numToSelect >= 0) {
      m_calculatedNumToSelect = m_numToSelect;
    }
    return m_calculatedNumToSelect;
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String thresholdTipText() {
    return "Set threshold by which attributes can be discarded. Default value "
      + "results in no attributes being discarded. Use either this option or "
      +"numToSelect to reduce the attribute set.";
  }

  /**
   * Set the threshold by which the AttributeSelection module can discard
   * attributes.
   * @param threshold the threshold.
   */
  public void setThreshold(double threshold) {
    m_threshold = threshold;
  }

  /**
   * Returns the threshold so that the AttributeSelection module can
   * discard attributes from the ranking.
   */
  public double getThreshold() {
    return m_threshold;
  }
  
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String generateRankingTipText() {
    return "A constant option. Ranker is only capable of generating "
      +" attribute rankings.";
  }

  /**
   * This is a dummy set method---Ranker is ONLY capable of producing
   * a ranked list of attributes for attribute evaluators.
   * @param doRank this parameter is N/A and is ignored
   */
  public void setGenerateRanking(boolean doRank) {
    
  }

  /**
   * This is a dummy method. Ranker can ONLY be used with attribute
   * evaluators and as such can only produce a ranked list of attributes
   * @return true all the time.
   */
  public boolean getGenerateRanking() {
    return true;
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String startSetTipText() {
    return "Specify a set of attributes to ignore. "
      +" When generating the ranking, Ranker will not evaluate the attributes "
      +" in this list. "
      +"This is specified as a comma " 
      +"seperated list off attribute indexes starting at 1. It can include "
      +"ranges. Eg. 1,2,5-9,17.";
  }

  /**
   * Sets a starting set of attributes for the search. It is the
   * search method's responsibility to report this start set (if any)
   * in its toString() method.
   * @param startSet a string containing a list of attributes (and or ranges),
   * eg. 1,2,6,10-15.
   * @throws Exception if start set can't be set.
   */
  public void setStartSet (String startSet) throws Exception {
    m_startRange.setRanges(startSet);
  }

  /**
   * Returns a list of attributes (and or attribute ranges) as a String
   * @return a list of attributes (and or attribute ranges)
   */
  public String getStartSet () {
    return m_startRange.getRanges();
  }

  /**
   * Returns an enumeration describing the available options.
   * @return an enumeration of all the available options.
   **/
  public Enumeration listOptions () {
    Vector newVector = new Vector(3);

    newVector
      .addElement(new Option("\tSpecify a starting set of attributes.\n" 
                             + "\tEg. 1,3,5-7.\n"
                             +"\tAny starting attributes specified are\n"
                             +"\tignored during the ranking."
                             ,"P",1
                             , "-P <start set>"));
    newVector
      .addElement(new Option("\tSpecify a theshold by which attributes\n" 
                             + "\tmay be discarded from the ranking.","T",1
                             , "-T <threshold>"));

    newVector
      .addElement(new Option("\tSpecify number of attributes to select" 
                             ,"N",1
                             , "-N <num to select>"));

    return newVector.elements();

  }
  
  /**
   * Parses a given list of options. <p/>
   *
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -P &lt;start set&gt;
   *  Specify a starting set of attributes.
   *  Eg. 1,3,5-7.
   *  Any starting attributes specified are
   *  ignored during the ranking.</pre>
   * 
   * <pre> -T &lt;threshold&gt;
   *  Specify a theshold by which attributes
   *  may be discarded from the ranking.</pre>
   * 
   * <pre> -N &lt;num to select&gt;
   *  Specify number of attributes to select</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions (String[] options)
    throws Exception {
    String optionString;
    resetOptions();

    optionString = Utils.getOption('P', options);
    if (optionString.length() != 0) {
      setStartSet(optionString);
    }

    optionString = Utils.getOption('T', options);
    if (optionString.length() != 0) {
      Double temp;
      temp = Double.valueOf(optionString);
      setThreshold(temp.doubleValue());
    }

    optionString = Utils.getOption('N', options);
    if (optionString.length() != 0) {
      setNumToSelect(Integer.parseInt(optionString));
    }
  }

  /**
   * Gets the current settings of ReliefFAttributeEval.
   *
   * @return an array of strings suitable for passing to setOptions()
   */
  public String[] getOptions () {
    String[] options = new String[6];
    int current = 0;

    if (!(getStartSet().equals(""))) {
      options[current++] = "-P";
      options[current++] = ""+startSetToString();
    }

    options[current++] = "-T";
    options[current++] = "" + getThreshold();

    options[current++] = "-N";
    options[current++] = ""+getNumToSelect();

    while (current < options.length) {
      options[current++] = "";
    }
    return  options;
  }

  /**
   * converts the array of starting attributes to a string. This is
   * used by getOptions to return the actual attributes specified
   * as the starting set. This is better than using m_startRanges.getRanges()
   * as the same start set can be specified in different ways from the
   * command line---eg 1,2,3 == 1-3. This is to ensure that stuff that
   * is stored in a database is comparable.
   * @return a comma seperated list of individual attribute numbers as a String
   */
  private String startSetToString() {
    StringBuffer FString = new StringBuffer();
    boolean didPrint;
    
    if (m_starting == null) {
      return getStartSet();
    }

    for (int i = 0; i < m_starting.length; i++) {
      didPrint = false;
      
      if ((m_hasClass == false) || 
          (m_hasClass == true && i != m_classIndex)) {
        FString.append((m_starting[i] + 1));
        didPrint = true;
      }
      
      if (i == (m_starting.length - 1)) {
        FString.append("");
      }
      else {
        if (didPrint) {
          FString.append(",");
        }
      }
    }
    
    return FString.toString();
  }

  /**
   * Kind of a dummy search algorithm. Calls a Attribute evaluator to
   * evaluate each attribute not included in the startSet and then sorts
   * them to produce a ranked list of attributes.
   *
   * @param ASEval the attribute evaluator to guide the search
   * @param data the training instances.
   * @return an array (not necessarily ordered) of selected attribute indexes
   * @throws Exception if the search can't be completed
   */
  public int[] search (ASEvaluation ASEval, Instances data)
    throws Exception {
    int i, j;

    if (!(ASEval instanceof AttributeEvaluator)) {
      throw  new Exception(ASEval.getClass().getName() 
                           + " is not a" 
                           + "Attribute evaluator!");
    }

    m_numAttribs = data.numAttributes();

    if (ASEval instanceof UnsupervisedAttributeEvaluator) {
      m_hasClass = false;
    }
    else {
      m_classIndex = data.classIndex();
      if (m_classIndex >= 0) {  
        m_hasClass = true;
      } else {
        m_hasClass = false;
      }
    }

    // get the transformed data and check to see if the transformer
    // preserves a class index
    if (ASEval instanceof AttributeTransformer) {
      data = ((AttributeTransformer)ASEval).transformedHeader();
      if (m_classIndex >= 0 && data.classIndex() >= 0) {
        m_classIndex = data.classIndex();
        m_hasClass = true;
      }
    }


    m_startRange.setUpper(m_numAttribs - 1);
    if (!(getStartSet().equals(""))) {
      m_starting = m_startRange.getSelection();
    }
    
    int sl=0;
    if (m_starting != null) {
      sl = m_starting.length;
    }
    if ((m_starting != null) && (m_hasClass == true)) {
      // see if the supplied list contains the class index
      boolean ok = false;
      for (i = 0; i < sl; i++) {
        if (m_starting[i] == m_classIndex) {
          ok = true;
          break;
        }
      }
      
      if (ok == false) {
        sl++;
      }
    }
    else {
      if (m_hasClass == true) {
        sl++;
      }
    }


    m_attributeList = new int[m_numAttribs - sl];
    m_attributeMerit = new double[m_numAttribs - sl];

    // add in those attributes not in the starting (omit list)
    for (i = 0, j = 0; i < m_numAttribs; i++) {
      if (!inStarting(i)) {
        m_attributeList[j++] = i;
      }
    }

    AttributeEvaluator ASEvaluator = (AttributeEvaluator)ASEval;

    for (i = 0; i < m_attributeList.length; i++) {
      m_attributeMerit[i] = ASEvaluator.evaluateAttribute(m_attributeList[i]);
    }

    double[][] tempRanked = rankedAttributes();
    int[] rankedAttributes = new int[m_attributeList.length];

    for (i = 0; i < m_attributeList.length; i++) {
      rankedAttributes[i] = (int)tempRanked[i][0];
    }

    return  rankedAttributes;
  }


  /**
   * Sorts the evaluated attribute list
   *
   * @return an array of sorted (highest eval to lowest) attribute indexes
   * @throws Exception of sorting can't be done.
   */
  public double[][] rankedAttributes ()
    throws Exception {
    int i, j;

    if (m_attributeList == null || m_attributeMerit == null) {
      throw  new Exception("Search must be performed before a ranked " 
                           + "attribute list can be obtained");
    }

    int[] ranked = Utils.sort(m_attributeMerit);
    // reverse the order of the ranked indexes
    double[][] bestToWorst = new double[ranked.length][2];

    for (i = ranked.length - 1, j = 0; i >= 0; i--) {
      bestToWorst[j++][0] = ranked[i];
    }

    // convert the indexes to attribute indexes
    for (i = 0; i < bestToWorst.length; i++) {
      int temp = ((int)bestToWorst[i][0]);
      bestToWorst[i][0] = m_attributeList[temp];
      bestToWorst[i][1] = m_attributeMerit[temp];
    }
    
    if (m_numToSelect > bestToWorst.length) {
      throw new Exception("More attributes requested than exist in the data");
    }

    if (m_numToSelect <= 0) {
      if (m_threshold == -Double.MAX_VALUE) {
        m_calculatedNumToSelect = bestToWorst.length;
      } else {
        determineNumToSelectFromThreshold(bestToWorst);
      }
    }
    /*    if (m_numToSelect > 0) {
      determineThreshFromNumToSelect(bestToWorst);
      } */

    return  bestToWorst;
  }

  private void determineNumToSelectFromThreshold(double [][] ranking) {
    int count = 0;
    for (int i = 0; i < ranking.length; i++) {
      if (ranking[i][1] > m_threshold) {
        count++;
      }
    }
    m_calculatedNumToSelect = count;
  }

  private void determineThreshFromNumToSelect(double [][] ranking) 
    throws Exception {
    if (m_numToSelect > ranking.length) {
      throw new Exception("More attributes requested than exist in the data");
    }

    if (m_numToSelect == ranking.length) {
      return;
    }

    m_threshold = (ranking[m_numToSelect-1][1] + 
                   ranking[m_numToSelect][1]) / 2.0;
  }

  /**
   * returns a description of the search as a String
   * @return a description of the search
   */
  public String toString () {
    StringBuffer BfString = new StringBuffer();
    BfString.append("\tAttribute ranking.\n");

    if (m_starting != null) {
      BfString.append("\tIgnored attributes: ");

      BfString.append(startSetToString());
      BfString.append("\n");
    }

    if (m_threshold != -Double.MAX_VALUE) {
      BfString.append("\tThreshold for discarding attributes: "
                      + Utils.doubleToString(m_threshold,8,4)+"\n");
    }

    return BfString.toString();
  }


  /**
   * Resets stuff to default values
   */
  protected void resetOptions () {
    m_starting = null;
    m_startRange = new Range();
    m_attributeList = null;
    m_attributeMerit = null;
    m_threshold = -Double.MAX_VALUE;
  }


  private boolean inStarting (int feat) {
    // omit the class from the evaluation
    if ((m_hasClass == true) && (feat == m_classIndex)) {
      return  true;
    }

    if (m_starting == null) {
      return  false;
    }

    for (int i = 0; i < m_starting.length; i++) {
      if (m_starting[i] == feat) {
        return  true;
      }
    }

    return  false;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 1.26 $");
  }
}
