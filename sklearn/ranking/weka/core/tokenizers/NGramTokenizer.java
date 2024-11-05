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
 * NGramTokenizer.java
 * Copyright (C) 2007 University of Waikato
 */

package weka.core.tokenizers;

import weka.core.Option;
import weka.core.RevisionUtils;
import weka.core.Utils;

import java.util.Enumeration;
import java.util.LinkedList;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * Splits a string into an n-gram with min and max grams.
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -delimiters &lt;value&gt;
 *  The delimiters to use
 *  (default ' \r\n\t.,;:'"()?!').</pre>
 * 
 * <pre> -max &lt;int&gt;
 *  The max size of the Ngram (default = 3).</pre>
 * 
 * <pre> -min &lt;int&gt;
 *  The min size of the Ngram (default = 1).</pre>
 * 
 <!-- options-end -->
 *
 * @author  Sebastian Germesin (sebastian.germesin@dfki.de)
 * @author  FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5953 $
 */
public class NGramTokenizer
  extends CharacterDelimitedTokenizer {

  /** for serialization */
  private static final long serialVersionUID = -2181896254171647219L;

  /** the maximum number of N */
  protected int m_NMax = 3;
  
  /** the minimum number of N */
  protected int m_NMin = 1;
  
  /** the current length of the N-grams */
  protected int m_N;
  
  /** the number of strings available */
  protected int m_MaxPosition;
  
  /** the current position for returning elements */
  protected int m_CurrentPosition;
  
  /** all the available grams */
  protected String[] m_SplitString;
  
  /**
   * Returns a string describing the stemmer
   * 
   * @return 		a description suitable for displaying in the 
   * 			explorer/experimenter gui
   */
  public String globalInfo() {
    return "Splits a string into an n-gram with min and max grams.";
  }
  
  /**
   * Returns an enumeration of all the available options..
   *
   * @return 		an enumeration of all available options.
   */
  public Enumeration listOptions() {
    Vector<Option>	result;
    Enumeration enm;
    
    result = new Vector<Option>();
    
    enm = super.listOptions();
    while (enm.hasMoreElements())
      result.addElement((Option)enm.nextElement());

    result.addElement(new Option(
	"\tThe max size of the Ngram (default = 3).",
	"max", 1, "-max <int>"));

    result.addElement(new Option(
	"\tThe min size of the Ngram (default = 1).",
	"min", 1, "-min <int>"));
    
    return result.elements();
  }
  
  /**
   * Gets the current option settings for the OptionHandler.
   *
   * @return 		the list of current option settings as an array of 
   * 			strings
   */
  public String[] getOptions() {
    Vector<String>	result;
    String[]		options;
    int			i;
    
    result = new Vector<String>();
    
    options = super.getOptions();
    for (i = 0; i < options.length; i++)
      result.add(options[i]);
    
    result.add("-max");
    result.add("" + getNGramMaxSize());

    result.add("-min");
    result.add("" + getNGramMinSize());

    return result.toArray(new String[result.size()]);
  }

  /**
   * Parses a given list of options. <p/>
   *
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -delimiters &lt;value&gt;
   *  The delimiters to use
   *  (default ' \r\n\t.,;:'"()?!').</pre>
   * 
   * <pre> -max &lt;int&gt;
   *  The max size of the Ngram (default = 3).</pre>
   * 
   * <pre> -min &lt;int&gt;
   *  The min size of the Ngram (default = 1).</pre>
   * 
   <!-- options-end -->
   *
   * @param options 	the list of options as an array of strings
   * @throws Exception 	if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String	value;
    
    super.setOptions(options);

    value = Utils.getOption("max", options);
    if (value.length() != 0)
      setNGramMaxSize(Integer.parseInt(value));
    else
      setNGramMaxSize(3);

    value = Utils.getOption("min", options);
    if (value.length() != 0)
      setNGramMinSize(Integer.parseInt(value));
    else
      setNGramMinSize(1);
  }
  
  /**
   * Gets the max N of the NGram.
   * 
   * @return 		the size (N) of the NGram.
   */
  public int getNGramMaxSize() {
    return m_NMax;
  }

  /**
   * Sets the max size of the Ngram.
   * 
   * @param value 	the size of the NGram.
   */
  public void setNGramMaxSize(int value) {
    if (value < 1)
      m_NMax = 1;
    else
      m_NMax = value;
  }

  /**
   * Returns the tip text for this property.
   *
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String NGramMaxSizeTipText() {
    return "The max N of the NGram.";
  }

  /**
   * Sets the min size of the Ngram.
   * 
   * @param value 	the size of the NGram.
   */
  public void setNGramMinSize(int value) {
    if (value < 1)
      m_NMin = 1;
    else
      m_NMin = value;
  }

  /**
   * Gets the min N of the NGram.
   * 
   * @return 		the size (N) of the NGram.
   */
  public int getNGramMinSize() {
    return m_NMin;
  }

  /**
   * Returns the tip text for this property.
   *
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String NGramMinSizeTipText() {
    return "The min N of the NGram.";
  }

  /**
   * returns true if there's more elements available
   * 
   * @return		true if there are more elements available
   */
  public boolean hasMoreElements() {
    return (m_CurrentPosition < m_MaxPosition && 
	m_N - 1 + m_CurrentPosition < m_MaxPosition && 
	m_N >= m_NMin);
  }
  
  /**
   * Returns N-grams and also (N-1)-grams and .... and 1-grams.
   * 
   * @return		the next element
   */
  public Object nextElement() {
    String retValue = "";
    
    for (int i = 0; i < m_N && i + m_CurrentPosition < m_MaxPosition; i++)
      retValue += " " + m_SplitString[m_CurrentPosition + i];
    
    m_CurrentPosition++;
    
    if (m_CurrentPosition + m_N - 1 == m_MaxPosition) {
      m_CurrentPosition = 0;
      m_N--;
    }

    return retValue.trim();
  }

  /** 
   * filters out empty strings in m_SplitString and
   * replaces m_SplitString with the cleaned version.
   * 
   * @see #m_SplitString
   */
  protected void filterOutEmptyStrings() {
    String[] newSplit;
    LinkedList<String> clean = new LinkedList<String>();

    for (int i = 0; i < m_SplitString.length; i++) {
      if (!m_SplitString[i].equals(""))
	clean.add(m_SplitString[i]);
    }

    newSplit = new String[clean.size()];
    for (int i = 0; i < clean.size(); i++) 
      newSplit[i] = clean.get(i);

    m_SplitString = newSplit;
  }
  
  /**
   * Sets the string to tokenize. Tokenization happens immediately.
   * 
   * @param s		the string to tokenize
   */
  public void tokenize(String s) {
    m_N           = m_NMax;
    m_SplitString = s.split("[" + getDelimiters() + "]");
    
    filterOutEmptyStrings();

    m_CurrentPosition = 0;
    m_MaxPosition     = m_SplitString.length;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5953 $");
  }

  /**
   * Runs the tokenizer with the given options and strings to tokenize.
   * The tokens are printed to stdout.
   * 
   * @param args	the commandline options and strings to tokenize
   */
  public static void main(String[] args) {
    runTokenizer(new NGramTokenizer(), args);
  }
}

