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
 *    M5Rules.java
 *    Copyright (C) 2001 University of Waikato, Hamilton, New Zealand
 */

package weka.classifiers.rules;

import weka.classifiers.trees.m5.M5Base;
import weka.core.RevisionUtils;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformationHandler;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;

/**
 <!-- globalinfo-start -->
 * Generates a decision list for regression problems using separate-and-conquer. In each iteration it builds a model tree using M5 and makes the "best" leaf into a rule.<br/>
 * <br/>
 * For more information see:<br/>
 * <br/>
 * Geoffrey Holmes, Mark Hall, Eibe Frank: Generating Rule Sets from Model Trees. In: Twelfth Australian Joint Conference on Artificial Intelligence, 1-12, 1999.<br/>
 * <br/>
 * Ross J. Quinlan: Learning with Continuous Classes. In: 5th Australian Joint Conference on Artificial Intelligence, Singapore, 343-348, 1992.<br/>
 * <br/>
 * Y. Wang, I. H. Witten: Induction of model trees for predicting continuous classes. In: Poster papers of the 9th European Conference on Machine Learning, 1997.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;inproceedings{Holmes1999,
 *    author = {Geoffrey Holmes and Mark Hall and Eibe Frank},
 *    booktitle = {Twelfth Australian Joint Conference on Artificial Intelligence},
 *    pages = {1-12},
 *    publisher = {Springer},
 *    title = {Generating Rule Sets from Model Trees},
 *    year = {1999}
 * }
 * 
 * &#64;inproceedings{Quinlan1992,
 *    address = {Singapore},
 *    author = {Ross J. Quinlan},
 *    booktitle = {5th Australian Joint Conference on Artificial Intelligence},
 *    pages = {343-348},
 *    publisher = {World Scientific},
 *    title = {Learning with Continuous Classes},
 *    year = {1992}
 * }
 * 
 * &#64;inproceedings{Wang1997,
 *    author = {Y. Wang and I. H. Witten},
 *    booktitle = {Poster papers of the 9th European Conference on Machine Learning},
 *    publisher = {Springer},
 *    title = {Induction of model trees for predicting continuous classes},
 *    year = {1997}
 * }
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -N
 *  Use unpruned tree/rules</pre>
 * 
 * <pre> -U
 *  Use unsmoothed predictions</pre>
 * 
 * <pre> -R
 *  Build regression tree/rule rather than a model tree/rule</pre>
 * 
 * <pre> -M &lt;minimum number of instances&gt;
 *  Set minimum number of instances per leaf
 *  (default 4)</pre>
 * 
 <!-- options-end -->
 *
 * @author <a href="mailto:mhall@cs.waikato.ac.nz">Mark Hall</a>
 * @version $Revision: 1.11 $
 */
public class M5Rules 
  extends M5Base
  implements TechnicalInformationHandler {
    
  /** for serialization */
  static final long serialVersionUID = -1746114858746563180L;
  
  /**
   * Returns a string describing classifier
   * @return a description suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {

    return "Generates a decision list for regression problems using " 
      + "separate-and-conquer. In each iteration it builds a "
      + "model tree using M5 and makes the \"best\" "
      + "leaf into a rule.\n\n"
      + "For more information see:\n\n"
      + getTechnicalInformation().toString();
  }

  /**
   * Constructor
   */
  public M5Rules() {
    super();
    setGenerateRules(true);
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
    
    result = new TechnicalInformation(Type.INPROCEEDINGS);
    result.setValue(Field.AUTHOR, "Geoffrey Holmes and Mark Hall and Eibe Frank");
    result.setValue(Field.TITLE, "Generating Rule Sets from Model Trees");
    result.setValue(Field.BOOKTITLE, "Twelfth Australian Joint Conference on Artificial Intelligence");
    result.setValue(Field.YEAR, "1999");
    result.setValue(Field.PAGES, "1-12");
    result.setValue(Field.PUBLISHER, "Springer");
    
    result.add(super.getTechnicalInformation());
    
    return result;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 1.11 $");
  }

  /**
   * Main method by which this class can be tested
   * 
   * @param args an array of options
   */
  public static void main(String[] args) {
    runClassifier(new M5Rules(), args);
  } 
}
