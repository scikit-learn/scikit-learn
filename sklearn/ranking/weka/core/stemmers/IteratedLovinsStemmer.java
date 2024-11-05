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
 * IteratedLovinsStemmer.java
 * Copyright (C) 2001 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core.stemmers;

import weka.core.RevisionUtils;

/**
 <!-- globalinfo-start -->
 * An iterated version of the Lovins stemmer. It stems the word (in case it's longer than 2 characters) until it no further changes.<br/>
 * <br/>
 * For more information about the Lovins stemmer see:<br/>
 * <br/>
 * Julie Beth Lovins (1968). Development of a stemming algorithm. Mechanical Translation and Computational Linguistics. 11:22-31.
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;article{Lovins1968,
 *    author = {Julie Beth Lovins},
 *    journal = {Mechanical Translation and Computational Linguistics},
 *    pages = {22-31},
 *    title = {Development of a stemming algorithm},
 *    volume = {11},
 *    year = {1968}
 * }
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
 *
 * @author  Eibe Frank (eibe at cs dot waikato dot ac dot nz)
 * @version $Revision: 5953 $
 * @see     LovinsStemmer
 */
public class IteratedLovinsStemmer 
  extends LovinsStemmer {

  /** for serialization */
  static final long serialVersionUID = 960689687163788264L;
  
  /**
   * Returns a string describing the stemmer
   * @return a description suitable for
   *         displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "An iterated version of the Lovins stemmer. It stems the word (in "
      + "case it's longer than 2 characters) until it no further changes.\n\n"
      + "For more information about the Lovins stemmer see:\n\n"
      + getTechnicalInformation().toString();
  }

  /**
   * Iterated stemming of the given word.
   * Word is converted to lower case.
   * 
   * @param str 	the word to stem
   * @return 		the stemmed word
   */
  public String stem(String str) {

    if (str.length() <= 2) {
      return str;
    }
    String stemmed = super.stem(str);
    while (!stemmed.equals(str)) {
      str = stemmed;
      stemmed = super.stem(stemmed);
    }
    return stemmed;
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
   * Runs the stemmer with the given options
   *
   * @param args      the options
   */
  public static void main(String[] args) {
    try {
      Stemming.useStemmer(new IteratedLovinsStemmer(), args);
    }
    catch (Exception e) {
      e.printStackTrace();
    }
  }
}

