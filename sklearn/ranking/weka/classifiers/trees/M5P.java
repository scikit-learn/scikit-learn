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
 *    M5P.java
 *    Copyright (C) 2001 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.trees;

import weka.classifiers.trees.m5.M5Base;
import weka.classifiers.trees.m5.Rule;
import weka.core.Drawable;
import weka.core.Option;
import weka.core.RevisionUtils;
import weka.core.Utils;

import java.util.Enumeration;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * M5Base. Implements base routines for generating M5 Model trees and rules<br/>
 * The original algorithm M5 was invented by R. Quinlan and Yong Wang made improvements.<br/>
 * <br/>
 * For more information see:<br/>
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
 * <pre> -L
 *  Save instances at the nodes in
 *  the tree (for visualization purposes)</pre>
 * 
 <!-- options-end -->
 *
 * @author <a href="mailto:mhall@cs.waikato.ac.nz">Mark Hall</a>
 * @version $Revision: 1.10 $
 */
public class M5P 
  extends M5Base 
  implements Drawable {

  /** for serialization */
  static final long serialVersionUID = -6118439039768244417L;
  
  /**
   * Creates a new <code>M5P</code> instance.
   */
  public M5P() {
    super();
    setGenerateRules(false);
  }

  /**
   *  Returns the type of graph this classifier
   *  represents.
   *  @return Drawable.TREE
   */   
  public int graphType() {
      return Drawable.TREE;
  }

  /**
   * Return a dot style String describing the tree.
   *
   * @return a <code>String</code> value
   * @throws Exception if an error occurs
   */
  public String graph() throws Exception {
    StringBuffer text = new StringBuffer();
    
    text.append("digraph M5Tree {\n");
    Rule temp = (Rule)m_ruleSet.elementAt(0);
    temp.topOfTree().graph(text);
    text.append("}\n");
    return text.toString();
  }

  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String saveInstancesTipText() {
    return 
        "Whether to save instance data at each node in the tree for "
      + "visualization purposes.";
  }

  /**
   * Set whether to save instance data at each node in the
   * tree for visualization purposes
   *
   * @param save a <code>boolean</code> value
   */
  public void setSaveInstances(boolean save) {
    m_saveInstances = save;
  }

  /**
   * Get whether instance data is being save.
   *
   * @return a <code>boolean</code> value
   */
  public boolean getSaveInstances() {
    return m_saveInstances;
  }

  /**
   * Returns an enumeration describing the available options
   * 
   * @return an enumeration of all the available options
   */
  public Enumeration listOptions() {
    Enumeration superOpts = super.listOptions();
    
    Vector newVector = new Vector();
    while (superOpts.hasMoreElements()) {
      newVector.addElement((Option)superOpts.nextElement());
    }

    newVector.addElement(new Option("\tSave instances at the nodes in\n"
				    +"\tthe tree (for visualization purposes)",
				    "L", 0, "-L"));
    return newVector.elements();
  }

  /**
   * Parses a given list of options. <p/>
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
   * <pre> -L
   *  Save instances at the nodes in
   *  the tree (for visualization purposes)</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    setSaveInstances(Utils.getFlag('L', options));
    super.setOptions(options);
  }

  /**
   * Gets the current settings of the classifier.
   * 
   * @return an array of strings suitable for passing to setOptions
   */
  public String [] getOptions() {
    String[] superOpts = super.getOptions();
    String [] options = new String [superOpts.length+1];
    int current = superOpts.length;
    for (int i = 0; i < current; i++) {
      options[i] = superOpts[i];
    }
    
    if (getSaveInstances()) {
      options[current++] = "-L";
    }

    while (current < options.length) {
      options[current++] = "";
    }

    return options;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 1.10 $");
  }

  /**
   * Main method by which this class can be tested
   * 
   * @param args an array of options
   */
  public static void main(String[] args) {
    runClassifier(new M5P(), args);
  } 
}
