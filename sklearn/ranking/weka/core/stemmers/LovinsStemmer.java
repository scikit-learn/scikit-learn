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
 * LovinsStemmer.java
 * Copyright (C) 2001 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core.stemmers;

import weka.core.RevisionUtils;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformation.Type;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformationHandler;

import java.util.HashMap;

/**
 <!-- globalinfo-start -->
 * A stemmer based on the Lovins stemmer, described here:<br/>
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
 */
public class LovinsStemmer 
  implements Stemmer, TechnicalInformationHandler {

  /** for serialization */
  static final long serialVersionUID = -6113024782588197L;
  
  /** Enters C version compatibility mode if set to true (emulates
    features of the original C implementation that are inconsistent
    with the algorithm as described in Lovins's paper) */
  private static boolean m_CompMode = false;

  /** The hash tables containing the list of endings. */
  private static HashMap<String,String> m_l11 = null;
  private static HashMap<String,String> m_l10 = null;
  private static HashMap<String,String> m_l9 = null;
  private static HashMap<String,String> m_l8 = null;
  private static HashMap<String,String> m_l7 = null;
  private static HashMap<String,String> m_l6 = null;
  private static HashMap<String,String> m_l5 = null;
  private static HashMap<String,String> m_l4 = null;
  private static HashMap<String,String> m_l3 = null;
  private static HashMap<String,String> m_l2 = null;
  private static HashMap<String,String> m_l1 = null;

  static {

    m_l11 = new HashMap<String,String>();
    m_l11.put("alistically", "B");
    m_l11.put("arizability", "A");
    m_l11.put("izationally", "B");
    m_l10 = new HashMap<String,String>();
    m_l10.put("antialness", "A");
    m_l10.put("arisations", "A");
    m_l10.put("arizations", "A");
    m_l10.put("entialness", "A");
    m_l9 = new HashMap<String,String>();
    m_l9.put("allically", "C");
    m_l9.put("antaneous", "A");
    m_l9.put("antiality", "A");
    m_l9.put("arisation", "A");
    m_l9.put("arization", "A");
    m_l9.put("ationally", "B");
    m_l9.put("ativeness", "A");
    m_l9.put("eableness", "E");
    m_l9.put("entations", "A");
    m_l9.put("entiality", "A");
    m_l9.put("entialize", "A");
    m_l9.put("entiation", "A");
    m_l9.put("ionalness", "A");
    m_l9.put("istically", "A");
    m_l9.put("itousness", "A");
    m_l9.put("izability", "A");
    m_l9.put("izational", "A");
    m_l8 = new HashMap<String,String>();
    m_l8.put("ableness", "A");
    m_l8.put("arizable", "A");
    m_l8.put("entation", "A");
    m_l8.put("entially", "A");
    m_l8.put("eousness", "A");
    m_l8.put("ibleness", "A");
    m_l8.put("icalness", "A");
    m_l8.put("ionalism", "A");
    m_l8.put("ionality", "A");
    m_l8.put("ionalize", "A");
    m_l8.put("iousness", "A");
    m_l8.put("izations", "A");
    m_l8.put("lessness", "A");
    m_l7 = new HashMap<String,String>();
    m_l7.put("ability", "A");
    m_l7.put("aically", "A");
    m_l7.put("alistic", "B");
    m_l7.put("alities", "A");
    m_l7.put("ariness", "E");
    m_l7.put("aristic", "A");
    m_l7.put("arizing", "A");
    m_l7.put("ateness", "A");
    m_l7.put("atingly", "A");
    m_l7.put("ational", "B");
    m_l7.put("atively", "A");
    m_l7.put("ativism", "A");
    m_l7.put("elihood", "E");
    m_l7.put("encible", "A");
    m_l7.put("entally", "A");
    m_l7.put("entials", "A");
    m_l7.put("entiate", "A");
    m_l7.put("entness", "A");
    m_l7.put("fulness", "A");
    m_l7.put("ibility", "A");
    m_l7.put("icalism", "A");
    m_l7.put("icalist", "A");
    m_l7.put("icality", "A");
    m_l7.put("icalize", "A");
    m_l7.put("ication", "G");
    m_l7.put("icianry", "A");
    m_l7.put("ination", "A");
    m_l7.put("ingness", "A");
    m_l7.put("ionally", "A");
    m_l7.put("isation", "A");
    m_l7.put("ishness", "A");
    m_l7.put("istical", "A");
    m_l7.put("iteness", "A");
    m_l7.put("iveness", "A");
    m_l7.put("ivistic", "A");
    m_l7.put("ivities", "A");
    m_l7.put("ization", "F");
    m_l7.put("izement", "A");
    m_l7.put("oidally", "A");
    m_l7.put("ousness", "A");
    m_l6 = new HashMap<String,String>();
    m_l6.put("aceous", "A");
    m_l6.put("acious", "B");
    m_l6.put("action", "G");
    m_l6.put("alness", "A");
    m_l6.put("ancial", "A");
    m_l6.put("ancies", "A");
    m_l6.put("ancing", "B");
    m_l6.put("ariser", "A");
    m_l6.put("arized", "A");
    m_l6.put("arizer", "A");
    m_l6.put("atable", "A");
    m_l6.put("ations", "B");
    m_l6.put("atives", "A");
    m_l6.put("eature", "Z");
    m_l6.put("efully", "A");
    m_l6.put("encies", "A");
    m_l6.put("encing", "A");
    m_l6.put("ential", "A");
    m_l6.put("enting", "C");
    m_l6.put("entist", "A");
    m_l6.put("eously", "A");
    m_l6.put("ialist", "A");
    m_l6.put("iality", "A");
    m_l6.put("ialize", "A");
    m_l6.put("ically", "A");
    m_l6.put("icance", "A");
    m_l6.put("icians", "A");
    m_l6.put("icists", "A");
    m_l6.put("ifully", "A");
    m_l6.put("ionals", "A");
    m_l6.put("ionate", "D");
    m_l6.put("ioning", "A");
    m_l6.put("ionist", "A");
    m_l6.put("iously", "A");
    m_l6.put("istics", "A");
    m_l6.put("izable", "E");
    m_l6.put("lessly", "A");
    m_l6.put("nesses", "A");
    m_l6.put("oidism", "A");
    m_l5 = new HashMap<String,String>();
    m_l5.put("acies", "A");
    m_l5.put("acity", "A");
    m_l5.put("aging", "B");
    m_l5.put("aical", "A");
    if (!m_CompMode) {
      m_l5.put("alist", "A");
    }
    m_l5.put("alism", "B");
    m_l5.put("ality", "A");
    m_l5.put("alize", "A");
    m_l5.put("allic", "b");
    m_l5.put("anced", "B");
    m_l5.put("ances", "B");
    m_l5.put("antic", "C");
    m_l5.put("arial", "A");
    m_l5.put("aries", "A");
    m_l5.put("arily", "A");
    m_l5.put("arity", "B");
    m_l5.put("arize", "A");
    m_l5.put("aroid", "A");
    m_l5.put("ately", "A");
    m_l5.put("ating", "I");
    m_l5.put("ation", "B");
    m_l5.put("ative", "A");
    m_l5.put("ators", "A");
    m_l5.put("atory", "A");
    m_l5.put("ature", "E");
    m_l5.put("early", "Y");
    m_l5.put("ehood", "A");
    m_l5.put("eless", "A");
    if (!m_CompMode) {
      m_l5.put("elily", "A");
    } else {
      m_l5.put("elity", "A");
    }
    m_l5.put("ement", "A");
    m_l5.put("enced", "A");
    m_l5.put("ences", "A");
    m_l5.put("eness", "E");
    m_l5.put("ening", "E");
    m_l5.put("ental", "A");
    m_l5.put("ented", "C");
    m_l5.put("ently", "A");
    m_l5.put("fully", "A");
    m_l5.put("ially", "A");
    m_l5.put("icant", "A");
    m_l5.put("ician", "A");
    m_l5.put("icide", "A");
    m_l5.put("icism", "A");
    m_l5.put("icist", "A");
    m_l5.put("icity", "A");
    m_l5.put("idine", "I");
    m_l5.put("iedly", "A");
    m_l5.put("ihood", "A");
    m_l5.put("inate", "A");
    m_l5.put("iness", "A");
    m_l5.put("ingly", "B");
    m_l5.put("inism", "J");
    m_l5.put("inity", "c");
    m_l5.put("ional", "A");
    m_l5.put("ioned", "A");
    m_l5.put("ished", "A");
    m_l5.put("istic", "A");
    m_l5.put("ities", "A");
    m_l5.put("itous", "A");
    m_l5.put("ively", "A");
    m_l5.put("ivity", "A");
    m_l5.put("izers", "F");
    m_l5.put("izing", "F");
    m_l5.put("oidal", "A");
    m_l5.put("oides", "A");
    m_l5.put("otide", "A");
    m_l5.put("ously", "A");
    m_l4 = new HashMap<String,String>();
    m_l4.put("able", "A");
    m_l4.put("ably", "A");
    m_l4.put("ages", "B");
    m_l4.put("ally", "B");
    m_l4.put("ance", "B");
    m_l4.put("ancy", "B");
    m_l4.put("ants", "B");
    m_l4.put("aric", "A");
    m_l4.put("arly", "K");
    m_l4.put("ated", "I");
    m_l4.put("ates", "A");
    m_l4.put("atic", "B");
    m_l4.put("ator", "A");
    m_l4.put("ealy", "Y");
    m_l4.put("edly", "E");
    m_l4.put("eful", "A");
    m_l4.put("eity", "A");
    m_l4.put("ence", "A");
    m_l4.put("ency", "A");
    m_l4.put("ened", "E");
    m_l4.put("enly", "E");
    m_l4.put("eous", "A");
    m_l4.put("hood", "A");
    m_l4.put("ials", "A");
    m_l4.put("ians", "A");
    m_l4.put("ible", "A");
    m_l4.put("ibly", "A");
    m_l4.put("ical", "A");
    m_l4.put("ides", "L");
    m_l4.put("iers", "A");
    m_l4.put("iful", "A");
    m_l4.put("ines", "M");
    m_l4.put("ings", "N");
    m_l4.put("ions", "B");
    m_l4.put("ious", "A");
    m_l4.put("isms", "B");
    m_l4.put("ists", "A");
    m_l4.put("itic", "H");
    m_l4.put("ized", "F");
    m_l4.put("izer", "F");
    m_l4.put("less", "A");
    m_l4.put("lily", "A");
    m_l4.put("ness", "A");
    m_l4.put("ogen", "A");
    m_l4.put("ward", "A");
    m_l4.put("wise", "A");
    m_l4.put("ying", "B");
    m_l4.put("yish", "A");
    m_l3 = new HashMap<String,String>();
    m_l3.put("acy", "A");
    m_l3.put("age", "B");
    m_l3.put("aic", "A");
    m_l3.put("als", "b");
    m_l3.put("ant", "B");
    m_l3.put("ars", "O");
    m_l3.put("ary", "F");
    m_l3.put("ata", "A");
    m_l3.put("ate", "A");
    m_l3.put("eal", "Y");
    m_l3.put("ear", "Y");
    m_l3.put("ely", "E");
    m_l3.put("ene", "E");
    m_l3.put("ent", "C");
    m_l3.put("ery", "E");
    m_l3.put("ese", "A");
    m_l3.put("ful", "A");
    m_l3.put("ial", "A");
    m_l3.put("ian", "A");
    m_l3.put("ics", "A");
    m_l3.put("ide", "L");
    m_l3.put("ied", "A");
    m_l3.put("ier", "A");
    m_l3.put("ies", "P");
    m_l3.put("ily", "A");
    m_l3.put("ine", "M");
    m_l3.put("ing", "N");
    m_l3.put("ion", "Q");
    m_l3.put("ish", "C");
    m_l3.put("ism", "B");
    m_l3.put("ist", "A");
    m_l3.put("ite", "a");
    m_l3.put("ity", "A");
    m_l3.put("ium", "A");
    m_l3.put("ive", "A");
    m_l3.put("ize", "F");
    m_l3.put("oid", "A");
    m_l3.put("one", "R");
    m_l3.put("ous", "A");
    m_l2 = new HashMap<String,String>();
    m_l2.put("ae", "A"); 
    m_l2.put("al", "b");
    m_l2.put("ar", "X");
    m_l2.put("as", "B");
    m_l2.put("ed", "E");
    m_l2.put("en", "F");
    m_l2.put("es", "E");
    m_l2.put("ia", "A");
    m_l2.put("ic", "A");
    m_l2.put("is", "A");
    m_l2.put("ly", "B");
    m_l2.put("on", "S");
    m_l2.put("or", "T");
    m_l2.put("um", "U");
    m_l2.put("us", "V");
    m_l2.put("yl", "R");
    m_l2.put("s\'", "A");
    m_l2.put("\'s", "A");
    m_l1 = new HashMap<String,String>();
    m_l1.put("a", "A");
    m_l1.put("e", "A");
    m_l1.put("i", "A");
    m_l1.put("o", "A");
    m_l1.put("s", "W");
    m_l1.put("y", "B");	
  }

  /**
   * Returns a string describing the stemmer
   * @return a description suitable for
   *         displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "A stemmer based on the Lovins stemmer, described here:\n\n"
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
    TechnicalInformation 	result;
    
    result = new TechnicalInformation(Type.ARTICLE);
    result.setValue(Field.AUTHOR, "Julie Beth Lovins");
    result.setValue(Field.YEAR, "1968");
    result.setValue(Field.TITLE, "Development of a stemming algorithm");
    result.setValue(Field.JOURNAL, "Mechanical Translation and Computational Linguistics");
    result.setValue(Field.VOLUME, "11");
    result.setValue(Field.PAGES, "22-31");

    return result;
  }

  /**
   * Finds and removes ending from given word.
   * 
   * @param word	the word to work on
   * @return 		the processed word
   */
  private String removeEnding(String word) {

    int length = word.length();
    int el = 11;

    while (el > 0) {
      if (length - el > 1) {
        String ending = word.substring(length - el);
        String conditionCode = null;
        switch (el) {
          case 11: conditionCode = (String)m_l11.get(ending);
                   break;
          case 10: conditionCode = (String)m_l10.get(ending);
                   break; 
          case 9: conditionCode = (String)m_l9.get(ending);
                  break;
          case 8: conditionCode = (String)m_l8.get(ending);
                  break;   
          case 7: conditionCode = (String)m_l7.get(ending);
                  break;   
          case 6: conditionCode = (String)m_l6.get(ending);
                  break;   
          case 5: conditionCode = (String)m_l5.get(ending);
                  break;   
          case 4: conditionCode = (String)m_l4.get(ending);
                  break;   
          case 3: conditionCode = (String)m_l3.get(ending);
                  break;   
          case 2: conditionCode = (String)m_l2.get(ending);
                  break;   
          case 1: conditionCode = (String)m_l1.get(ending);
                  break;   
          default:
        }
        if (conditionCode != null) {
          switch (conditionCode.charAt(0)) {
            case 'A':
              return word.substring(0, length - el);
            case 'B':
              if (length - el > 2) {
                return word.substring(0, length - el);
              }
              break;
            case 'C':
              if (length - el > 3) {
                return word.substring(0, length - el);
              }
              break;
            case 'D':
              if (length - el > 4) {
                return word.substring(0, length - el);
              }
              break;
            case 'E':
              if (word.charAt(length - el - 1) != 'e') {
                return word.substring(0, length - el);
              }
              break;
            case 'F':
              if ((length - el > 2) &&
                  (word.charAt(length - el - 1) != 'e')) {
                return word.substring(0, length - el);
                  }
              break;
            case 'G':
              if ((length - el > 2) &&
                  (word.charAt(length - el - 1) == 'f')) {
                return word.substring(0, length - el);
                  }
              break;
            case 'H':
              if ((word.charAt(length - el - 1) == 't') ||
                  ((word.charAt(length - el - 1) == 'l') &&
                   (word.charAt(length - el - 2) == 'l'))) {
                return word.substring(0, length - el);
                   }
              break;
            case 'I':
              if ((word.charAt(length - el - 1) != 'o') &&
                  (word.charAt(length - el - 1) != 'e')) { 
                return word.substring(0, length - el);
                  }
              break;
            case 'J':
              if ((word.charAt(length - el - 1) != 'a') &&
                  (word.charAt(length - el - 1) != 'e')) { 
                return word.substring(0, length - el);
                  }
              break;
            case 'K':
              if ((length - el > 2) &&
                  ((word.charAt(length - el - 1) == 'l') ||
                   (word.charAt(length - el - 1) == 'i') ||
                   ((word.charAt(length - el - 1) == 'e') &&
                    (word.charAt(length - el - 3) == 'u')))) {
                return word.substring(0, length - el);
                    }
              break;
            case 'L':
              if ((word.charAt(length - el - 1) != 'u') &&
                  (word.charAt(length - el - 1) != 'x') &&
                  ((word.charAt(length - el - 1) != 's') ||
                   (word.charAt(length - el - 2) == 'o'))) {
                return word.substring(0, length - el);
                   }
              break;
            case 'M':
              if ((word.charAt(length - el - 1) != 'a') &&
                  (word.charAt(length - el - 1) != 'c') &&
                  (word.charAt(length - el - 1) != 'e') &&
                  (word.charAt(length - el - 1) != 'm')) {
                return word.substring(0, length - el);
                  }
              break;
            case 'N':
              if ((length - el > 3) || 
                  ((length - el == 3) &&
                   ((word.charAt(length - el - 3) != 's')))) {
                return word.substring(0, length - el);
                   }
              break;
            case 'O':
              if ((word.charAt(length - el - 1) == 'l') ||
                  (word.charAt(length - el - 1) == 'i')) {
                return word.substring(0, length - el);
                  }
              break;
            case 'P':
              if (word.charAt(length - el - 1) != 'c') {
                return word.substring(0, length - el);
              }
              break;
            case 'Q':
              if ((length - el > 2) &&
                  (word.charAt(length - el - 1) != 'l') &&
                  (word.charAt(length - el - 1) != 'n')) {
                return word.substring(0, length - el);
                  }
              break;
            case 'R':
              if ((word.charAt(length - el - 1) == 'n') ||
                  (word.charAt(length - el - 1) == 'r')) {
                return word.substring(0, length - el);
                  }
              break;
            case 'S':
              if (((word.charAt(length - el - 1) == 'r') &&
                    (word.charAt(length - el - 2) == 'd')) ||
                  ((word.charAt(length - el - 1) == 't') &&
                   (word.charAt(length - el - 2) != 't'))) {
                return word.substring(0, length - el);
                   }
              break;
            case 'T':
              if ((word.charAt(length - el - 1) == 's') ||
                  ((word.charAt(length - el - 1) == 't') &&
                   (word.charAt(length - el - 2) != 'o'))) {
                return word.substring(0, length - el);
                   }
              break;
            case 'U':
              if ((word.charAt(length - el - 1) == 'l') ||
                  (word.charAt(length - el - 1) == 'm') ||
                  (word.charAt(length - el - 1) == 'n') ||
                  (word.charAt(length - el - 1) == 'r')) {
                return word.substring(0, length - el);
                  }
              break;
            case 'V':
              if (word.charAt(length - el - 1) == 'c') {
                return word.substring(0, length - el);
              }
              break;
            case 'W':
              if ((word.charAt(length - el - 1) != 's') &&
                  (word.charAt(length - el - 1) != 'u')) {
                return word.substring(0, length - el);
                  }
              break;
            case 'X':
              if ((word.charAt(length - el - 1) == 'l') ||
                  (word.charAt(length - el - 1) == 'i') ||
                  ((length - el > 2) &&
                   (word.charAt(length - el - 1) == 'e') &&
                   (word.charAt(length - el - 3) == 'u'))) {
                return word.substring(0, length - el);
                   }
              break;
            case 'Y':
              if ((word.charAt(length - el - 1) == 'n') &&
                  (word.charAt(length - el - 2) == 'i')) {
                return word.substring(0, length - el);
                  }
              break;
            case 'Z':
              if (word.charAt(length - el - 1) != 'f') {
                return word.substring(0, length - el);
              }
              break;
            case 'a':
              if ((word.charAt(length - el - 1) == 'd') ||
                  (word.charAt(length - el - 1) == 'f') ||
                  (((word.charAt(length - el - 1) == 'h') &&
                    (word.charAt(length - el - 2) == 'p'))) ||
                  (((word.charAt(length - el - 1) == 'h') &&
                    (word.charAt(length - el - 2) == 't'))) ||
                  (word.charAt(length - el - 1) == 'l') ||
                  (((word.charAt(length - el - 1) == 'r') &&
                    (word.charAt(length - el - 2) == 'e'))) ||
                  (((word.charAt(length - el - 1) == 'r') &&
                    (word.charAt(length - el - 2) == 'o'))) ||
                  (((word.charAt(length - el - 1) == 's') &&
                    (word.charAt(length - el - 2) == 'e'))) ||
                  (word.charAt(length - el - 1) == 't')) {
                return word.substring(0, length - el);
                  }
              break;
            case 'b':
              if (m_CompMode) {
                if (((length - el == 3 ) &&
                      (!((word.charAt(length - el - 1) == 't') &&
                         (word.charAt(length - el - 2) == 'e') &&
                         (word.charAt(length - el - 3) == 'm')))) ||
                    ((length - el > 3) &&
                     (!((word.charAt(length - el - 1) == 't') &&
                        (word.charAt(length - el - 2) == 's') &&
                        (word.charAt(length - el - 3) == 'y') &&
                        (word.charAt(length - el - 4) == 'r'))))) {
                  return word.substring(0, length - el);
                        }
              } else {
                if ((length - el > 2) &&
                    (!((word.charAt(length - el - 1) == 't') &&
                       (word.charAt(length - el - 2) == 'e') &&
                       (word.charAt(length - el - 3) == 'm'))) &&
                    ((length - el < 4) ||
                     (!((word.charAt(length - el - 1) == 't') &&
                        (word.charAt(length - el - 2) == 's') &&
                        (word.charAt(length - el - 3) == 'y') &&
                        (word.charAt(length - el - 4) == 'r'))))) {
                  return word.substring(0, length - el);
                        }
              } 
              break;
            case 'c':
              if (word.charAt(length - el - 1) == 'l') {
                return word.substring(0, length - el);
              }
              break;
            default:
              throw new IllegalArgumentException("Fatal error.");
          }
        }
      }
      el--;
    }
    return word;
  }

  /**
   * Recodes ending of given word.
   * 
   * @param word	the word to work on
   * @return		the processed word
   */
  private String recodeEnding(String word) {

    int lastPos = word.length() - 1;

    // Rule 1
    if (word.endsWith("bb") ||
        word.endsWith("dd") ||
        word.endsWith("gg") ||
        word.endsWith("ll") ||
        word.endsWith("mm") ||
        word.endsWith("nn") ||
        word.endsWith("pp") ||
        word.endsWith("rr") ||
        word.endsWith("ss") ||
        word.endsWith("tt")) {
      word = word.substring(0, lastPos);
      lastPos--;
        }

    // Rule 2
    if (word.endsWith("iev")) {
      word = word.substring(0, lastPos - 2).concat("ief");
    }

    // Rule 3
    if (word.endsWith("uct")) {
      word = word.substring(0, lastPos - 2).concat("uc");
      lastPos--;
    }

    // Rule 4
    if (word.endsWith("umpt")) {
      word = word.substring(0, lastPos - 3).concat("um");
      lastPos -= 2;
    }

    // Rule 5
    if (word.endsWith("rpt")) {
      word = word.substring(0, lastPos - 2).concat("rb");
      lastPos--;
    }

    // Rule 6
    if (word.endsWith("urs")) {
      word = word.substring(0, lastPos - 2).concat("ur");
      lastPos--;
    }

    // Rule 7
    if (word.endsWith("istr")) {
      word = word.substring(0, lastPos - 3).concat("ister");
      lastPos++;
    }

    // Rule 7a
    if (word.endsWith("metr")) {
      word = word.substring(0, lastPos - 3).concat("meter");
      lastPos++;
    }

    // Rule 8
    if (word.endsWith("olv")) {
      word = word.substring(0, lastPos - 2).concat("olut");
      lastPos++;
    }

    // Rule 9
    if (word.endsWith("ul")) {
      if ((lastPos - 2 < 0) ||
          ((word.charAt(lastPos - 2) != 'a') &&
           (word.charAt(lastPos - 2) != 'i') &&
           (word.charAt(lastPos - 2) != 'o'))) {
        word = word.substring(0, lastPos - 1).concat("l");
        lastPos--;
           }
    }

    // Rule 10
    if (word.endsWith("bex")) {
      word = word.substring(0, lastPos - 2).concat("bic");
    }

    // Rule 11
    if (word.endsWith("dex")) {
      word = word.substring(0, lastPos - 2).concat("dic");
    }

    // Rule 12
    if (word.endsWith("pex")) {
      word = word.substring(0, lastPos - 2).concat("pic");
    }

    // Rule 13
    if (word.endsWith("tex")) {
      word = word.substring(0, lastPos - 2).concat("tic");
    }

    // Rule 14
    if (word.endsWith("ax")) {
      word = word.substring(0, lastPos - 1).concat("ac");
    }

    // Rule 15
    if (word.endsWith("ex")) {
      word = word.substring(0, lastPos - 1).concat("ec");
    }

    // Rule 16
    if (word.endsWith("ix")) {
      word = word.substring(0, lastPos - 1).concat("ic");
    }

    // Rule 17
    if (word.endsWith("lux")) {
      word = word.substring(0, lastPos - 2).concat("luc");
    }

    // Rule 18
    if (word.endsWith("uad")) {
      word = word.substring(0, lastPos - 2).concat("uas");
    }

    // Rule 19
    if (word.endsWith("vad")) {
      word = word.substring(0, lastPos - 2).concat("vas");
    }

    // Rule 20
    if (word.endsWith("cid")) {
      word = word.substring(0, lastPos - 2).concat("cis");
    }

    // Rule 21
    if (word.endsWith("lid")) {
      word = word.substring(0, lastPos - 2).concat("lis");
    }

    // Rule 22
    if (word.endsWith("erid")) {
      word = word.substring(0, lastPos - 3).concat("eris");
    }

    // Rule 23
    if (word.endsWith("pand")) {
      word = word.substring(0, lastPos - 3).concat("pans");
    }

    // Rule 24
    if (word.endsWith("end")) {
      if ((lastPos - 3 < 0) ||
          (word.charAt(lastPos - 3) != 's')) {
        word = word.substring(0, lastPos - 2).concat("ens");
          }
    }

    // Rule 25
    if (word.endsWith("ond")) {
      word = word.substring(0, lastPos - 2).concat("ons");
    }

    // Rule 26
    if (word.endsWith("lud")) {
      word = word.substring(0, lastPos - 2).concat("lus");
    }

    // Rule 27
    if (word.endsWith("rud")) {
      word = word.substring(0, lastPos - 2).concat("rus");
    }

    // Rule 28
    if (word.endsWith("her")) {
      if ((lastPos - 3 < 0) ||
          ((word.charAt(lastPos - 3) != 'p') &&
           (word.charAt(lastPos - 3) != 't'))) {
        word = word.substring(0, lastPos - 2).concat("hes");
           }
    }

    // Rule 29
    if (word.endsWith("mit")) {
      word = word.substring(0, lastPos - 2).concat("mis");
    }

    // Rule 30
    if (word.endsWith("end")) {
      if ((lastPos - 3 < 0) ||
          (word.charAt(lastPos - 3) != 'm')) {
        word = word.substring(0, lastPos - 2).concat("ens");
          }
    }

    // Rule 31
    if (word.endsWith("ert")) {
      word = word.substring(0, lastPos - 2).concat("ers");
    }

    // Rule 32
    if (word.endsWith("et")) {
      if ((lastPos - 2 < 0) ||
          (word.charAt(lastPos - 2) != 'n')) {
        word = word.substring(0, lastPos - 1).concat("es");
          }
    }

    // Rule 33
    if (word.endsWith("yt")) {
      word = word.substring(0, lastPos - 1).concat("ys");
    }

    // Rule 34
    if (word.endsWith("yz")) {
      word = word.substring(0, lastPos - 1).concat("ys");
    }

    return word;
  }

  /**
   * Returns the stemmed version of the given word.
   * Word is converted to lower case before stemming.
   * 
   * @param word 	a string consisting of a single word
   * @return 		the stemmed word
   */
  public String stem(String word) {

    if (word.length() > 2) {
      return recodeEnding(removeEnding(word.toLowerCase()));
    } else {
      return word.toLowerCase();
    }
  }

  /**
   * Stems everything in the given string. String
   * is converted to lower case before stemming.
   * 
   * @param str		the string to stem
   * @return 		the processed string
   */
  public String stemString(String str) {

    StringBuffer result = new StringBuffer();
    int start = -1;
    for (int j = 0; j < str.length(); j++) {
      char c = str.charAt(j);
      if (Character.isLetterOrDigit(c)) {
        if (start == -1) {
          start = j;
        }
      } else if (c == '\'') {
        if (start == -1) {
          result.append(c);
        }
      } else {
        if (start != -1) {
          result.append(stem(str.substring(start, j)));
          start = -1;
        }
        result.append(c);
      }
    }
    if (start != -1) {
      result.append(stem(str.substring(start, str.length())));
    }
    return result.toString();  
  }

  /**
   * returns a string representation of the stemmer
   * 
   * @return a string representation of the stemmer
   */
  public String toString() {
    return getClass().getName();
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
      Stemming.useStemmer(new LovinsStemmer(), args);
    }
    catch (Exception e) {
      e.printStackTrace();
    }
  }
}

