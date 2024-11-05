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
 * SimpleStringTokenizer.java
 * Copyright (C) 2007 University of Waikato, Hamilton, New Zealand
 */

package weka.core.tokenizers;

import weka.core.RevisionUtils;

import java.util.StringTokenizer;

/**
 <!-- globalinfo-start -->
 * A simple tokenizer that is using the java.util.StringTokenizer class to tokenize the strings.
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
 <!-- options-end -->
 *
 * @author  FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5953 $
 */
public class WordTokenizer
  extends CharacterDelimitedTokenizer {

  /** for serialization */
  private static final long serialVersionUID = -930893034037880773L;
  
  /** the actual tokenizer */
  protected transient StringTokenizer m_Tokenizer;
  
  /**
   * Returns a string describing the stemmer
   * 
   * @return 		a description suitable for displaying in the 
   * 			explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "A simple tokenizer that is using the java.util.StringTokenizer "
      + "class to tokenize the strings.";
  }

  /**
   * Tests if this enumeration contains more elements.
   * 
   * @return 		true if and only if this enumeration object contains 
   * 			at least one more element to provide; false otherwise.
   */
  public boolean hasMoreElements() {
    return m_Tokenizer.hasMoreElements();
  }

  /**
   * Returns the next element of this enumeration if this enumeration object 
   * has at least one more element to provide.
   * 
   * @return		the next element of this enumeration.
   */
  public Object nextElement() {
    return m_Tokenizer.nextElement();
  }
  
  /**
   * Sets the string to tokenize. Tokenization happens immediately.
   * 
   * @param s		the string to tokenize
   */
  public void tokenize(String s) {
    m_Tokenizer = new StringTokenizer(s, getDelimiters());
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
    runTokenizer(new WordTokenizer(), args);
  }
}

