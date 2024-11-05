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
 * BallSplitter.java
 * Copyright (C) 2007 University of Waikato
 */

package weka.core.neighboursearch.balltrees;

import weka.core.EuclideanDistance;
import weka.core.Instances;
import weka.core.OptionHandler;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;

import java.io.Serializable;
import java.util.Enumeration;
import java.util.Vector;

/**
 * Abstract class for splitting a ball tree's BallNode.
 * 
 * @author Ashraf M. Kibriya (amk14[at-the-rate]cs[dot]waikato[dot]ac[dot]nz)
 * @version $Revision: 5953 $
 */
public abstract class BallSplitter
  implements Serializable, OptionHandler, RevisionHandler {
  
  /** The instance on which the tree is built. */
  protected Instances m_Instances;
  
  /** The distance function (metric) from which
   * the tree is (OR is to be) built. */
  protected EuclideanDistance m_DistanceFunction;
  
  /** 
   * The master index array that'll be reshuffled as nodes
   * are split (and the tree is constructed). 
   */
  protected int[] m_Instlist; 

  /**
   * default constructor.
   */
  public BallSplitter() {
  }
  
  /**
   * Creates a new instance of BallSplitter.
   * @param instList The master index array.
   * @param insts The instances on which the tree
   * is (or is to be) built.
   * @param e The Euclidean distance function to 
   * use for splitting.
   */
  public BallSplitter(int[] instList, Instances insts, EuclideanDistance e) { 
    m_Instlist = instList;
    m_Instances = insts;
    m_DistanceFunction = e;
  }

  /**
   * Checks whether if this ball splitter is 
   * correctly intialized or not (i.e. master index
   * array, instances, and distance function is 
   * supplied or not)
   * @throws Exception If the object is not correctly
   * initialized.
   */
  protected void correctlyInitialized() throws Exception {
    if(m_Instances==null)
      throw new Exception("No instances supplied.");
    else if(m_Instlist==null) 
      throw new Exception("No instance list supplied.");
    else if(m_DistanceFunction==null)
      throw new Exception("No Euclidean distance function supplied.");
    else if(m_Instances.numInstances() != m_Instlist.length)
      throw new Exception("The supplied instance list doesn't seem to match " +
                          "the supplied instances");
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return 		an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    return new Vector().elements();
  }

  /**
   * Parses a given list of options.
   * 
   * @param options 	the list of options as an array of strings
   * @throws Exception 	if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
  }

  /**
   * Gets the current settings of the object.
   *
   * @return 		an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    return new String[0];
  }
  
  /** 
   * Splits a node into two. 
   * @param node The node to split.
   * @param numNodesCreated The number of nodes that so far have been
   * created for the tree, so that the newly created nodes are 
   * assigned correct/meaningful node numbers/ids.
   * @throws Exception If there is some problem in splitting the
   * given node.
   */
  public abstract void splitNode(BallNode node, int numNodesCreated) 
      throws Exception;
  
  /**
   * Sets the training instances on which the tree is 
   * (or is to be) built. 
   * @param inst The training instances.
   */
  public void setInstances(Instances inst) {
    m_Instances = inst;
  }
  
  /** 
   * Sets the master index array containing indices of the
   * training instances. This array will be rearranged as 
   * the tree is built (or a node is split_), so that each 
   * node is assigned a portion in this array which 
   * contain the instances insides the node's region.
   * @param instList The master index array.
   */
  public void setInstanceList(int[] instList) {
    m_Instlist = instList;
  }
  
  /**
   * Sets the distance function used to (or to be used 
   * to) build the tree. 
   * @param func The distance function. 
   */
  public void setEuclideanDistanceFunction(EuclideanDistance func) {
    m_DistanceFunction = func;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5953 $");
  }
}
