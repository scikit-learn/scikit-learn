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
 * CoverTree.java
 * Copyright (C) 2006 Alina Beygelzimer and Sham Kakade and John Langford
 */

package weka.core.neighboursearch;

import weka.core.DistanceFunction;
import weka.core.EuclideanDistance;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformationHandler;
import weka.core.Utils;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;
import weka.core.converters.CSVLoader;
import weka.core.neighboursearch.covertrees.Stack;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.Serializable;
import java.util.Enumeration;
import java.util.List;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * Class implementing the CoverTree datastructure.<br/>
 * The class is very much a translation of the c source code made available by the authors.<br/>
 * <br/>
 * For more information and original source code see:<br/>
 * <br/>
 * Alina Beygelzimer, Sham Kakade, John Langford: Cover trees for nearest neighbor. In: ICML'06: Proceedings of the 23rd international conference on Machine learning, New York, NY, USA, 97-104, 2006.
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;inproceedings{Beygelzimer2006,
 *    address = {New York, NY, USA},
 *    author = {Alina Beygelzimer and Sham Kakade and John Langford},
 *    booktitle = {ICML'06: Proceedings of the 23rd international conference on Machine learning},
 *    pages = {97-104},
 *    publisher = {ACM Press},
 *    title = {Cover trees for nearest neighbor},
 *    year = {2006},
 *    location = {Pittsburgh, Pennsylvania},
 *    HTTP = {http://hunch.net/\~jl/projects/cover_tree/cover_tree.html}
 * }
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
 * 
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -B &lt;value&gt;
 *  Set base of the expansion constant
 *  (default = 1.3).</pre>
 * 
 <!-- options-end -->
 * 
 * @author Alina Beygelzimer (original C++ code)
 * @author Sham Kakade (original C++ code)
 * @author John Langford (original C++ code)
 * @author Ashraf M. Kibriya (amk14[at-the-rate]cs[dot]waikato[dot]ac[dot]nz) (Java port)
 * @version $Revision: 5953 $
 */
public class CoverTree
  extends NearestNeighbourSearch
  implements TechnicalInformationHandler {

  /** for serialization. */
  private static final long serialVersionUID = 7617412821497807586L;

  /**
   * class representing a node of the cover tree.
   * 
   * @author Ashraf M. Kibriya (amk14[at-the-rate]cs[dot]waikato[dot]ac[dot]nz)
   * @version $Revision: 5953 $
   */
  public class CoverTreeNode
    implements Serializable, RevisionHandler {
    
    /** for serialization. */
    private static final long serialVersionUID = 1808760031169036512L;
    
    /** ID for the node. */
    private int nodeid;
    
    /** Index of the instance represented by this node in the index array. */
    private Integer idx;
    
    /** The distance of the furthest descendant of the node. */
    private double max_dist; // The maximum distance to any grandchild.

    /** The distance to the nodes parent. */ 
    private double parent_dist; // The distance to the parent.

    /** The children of the node. */
    private Stack<CoverTreeNode> children;

    /** The number of children node has.  */
    private int num_children; // The number of children.

    /** The min i that makes base^i &lt;= max_dist. */
    private int scale; // Essentially, an upper bound on the distance to any child.

    /** Constructor for the class. */
    public CoverTreeNode() {
    }
    
    /**
     * Constructor.
     * @param i The index of the Instance this node is
     * associated with.
     * @param md The distance of the furthest descendant.
     * @param pd The distance of the node to its parent.
     * @param childs Children of the node in a stack.
     * @param numchilds The number of children of the 
     * node.
     * @param s The scale/level of the node in the tree.
     */
    public CoverTreeNode(Integer i, double md, double pd,
      Stack<CoverTreeNode> childs, int numchilds, int s) {
      idx = i;
      max_dist = md;
      parent_dist = pd;
      children = childs;
      num_children = numchilds;
      scale = s;
    }
    
    /** Returns the instance represented by the node.
     * @return The instance represented by the node.
     */
    public Instance p() {
      return m_Instances.instance(idx);
    }
    
    /** Returns whether if the node is a leaf or not.
     * @return true if the node is a leaf node. 
     */
    public boolean isALeaf() {
      return num_children==0;
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

  /**
   * Private class holding a point's distance to the current reference
   * point p.
   * 
   * @author Ashraf M. Kibriya (amk14[at-the-rate]cs[dot]waikato[dot]ac[dot]nz)
   * @version $Revision: 5953 $
   */
  private class DistanceNode
    implements RevisionHandler {
    
    /**
     * The last distance is to the current reference point
     * (potential current parent). The previous ones are
     * to reference points that were previously looked at
     * (all potential ancestors).      
     */
    Stack<Double> dist;
    
    /** The index of the instance represented by this node. */
    Integer idx;
    
    /**
     * Returns the instance represent by this DistanceNode.
     * @return The instance represented by this node. 
     */
    public Instance q() {
      return m_Instances.instance(idx);
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

  /** The euclidean distance function to use. */
  protected EuclideanDistance m_EuclideanDistance;
  { // to make sure we have only one object of EuclideanDistance
    if (m_DistanceFunction instanceof EuclideanDistance)
      m_EuclideanDistance = (EuclideanDistance) m_DistanceFunction;
    else
      m_DistanceFunction = m_EuclideanDistance = new EuclideanDistance();
  }

  /** The root node. */
  protected CoverTreeNode m_Root;

  /** 
   * Array holding the distances of the nearest neighbours. It is filled up
   *  both by nearestNeighbour() and kNearestNeighbours(). 
   */
  protected double [] m_DistanceList;

  /** Number of nodes in the tree. */
  protected int m_NumNodes, m_NumLeaves, m_MaxDepth;
  
  /** Tree Stats variables. */
  protected TreePerformanceStats m_TreeStats = null;

  /**
   * The base of our expansion constant. In other words the 2 in 2^i used
   * in covering tree and separation invariants of a cover tree. P.S.: In
   * paper it's suggested the separation invariant is relaxed in batch
   * construction.
   */
  protected double m_Base = 1.3;

  /**
   * if we have base 2 then this can be viewed as 1/ln(2), which can be
   * used later on to do il2*ln(d) instead of ln(d)/ln(2), to get log2(d),
   * in get_scale method.
   */
  protected double il2 = 1.0 / Math.log(m_Base);

  /**
   * default constructor.
   */
  public CoverTree() {
    super();
    if(getMeasurePerformance())
      m_Stats = m_TreeStats = new TreePerformanceStats();
  }

  /**
   * Returns a string describing this nearest neighbour search algorithm.
   * 
   * @return 		a description of the algorithm for displaying in the 
   * 			explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "Class implementing the CoverTree datastructure.\n"
      + "The class is very much a translation of the c source code made "
      + "available by the authors.\n\n"
      + "For more information and original source code see:\n\n"
      + getTechnicalInformation().toString();
  }

  /**
   * Returns an instance of a TechnicalInformation object, containing detailed
   * information about the technical background of this class, e.g., paper
   * reference or book this class is based on.
   * 
   * @return 		the technical information about this class
   */
  public TechnicalInformation getTechnicalInformation() {
    TechnicalInformation result;

    result = new TechnicalInformation(Type.INPROCEEDINGS);
    result.setValue(Field.AUTHOR, "Alina Beygelzimer and Sham Kakade and John Langford");
    result.setValue(Field.TITLE, "Cover trees for nearest neighbor");
    result.setValue(Field.BOOKTITLE, "ICML'06: Proceedings of the 23rd international conference on Machine learning");
    result.setValue(Field.PAGES, "97-104");
    result.setValue(Field.YEAR, "2006");
    result.setValue(Field.PUBLISHER, "ACM Press");
    result.setValue(Field.ADDRESS, "New York, NY, USA");
    result.setValue(Field.LOCATION, "Pittsburgh, Pennsylvania");
    result.setValue(Field.HTTP, "http://hunch.net/~jl/projects/cover_tree/cover_tree.html");

    return result;
  }
  
  /**
   * Returns an enumeration describing the available options.
   * 
   * @return 		an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector<Option> newVector = new Vector<Option>();

    newVector.addElement(new Option(
	"\tSet base of the expansion constant\n"
	+ "\t(default = 1.3).",
	"B", 1, "-B <value>"));
    
    return newVector.elements();
  }
  
  /**
   * Parses a given list of options. <p/>
   *
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -B &lt;value&gt;
   *  Set base of the expansion constant
   *  (default = 1.3).</pre>
   * 
   <!-- options-end -->
   * 
   * @param options 	the list of options as an array of strings
   * @throws Exception	if an option is not supported
   */
  public void setOptions(String[] options)
    throws Exception {    
    
    super.setOptions(options);
    
    String optionString = Utils.getOption('B', options);
    if (optionString.length() != 0)
      setBase(Double.parseDouble(optionString));
    else
      setBase(1.3);      
  }

  /**
   * Gets the current settings of KDtree.
   * 
   * @return 		an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    Vector<String>	result;
    String[]		options;
    int			i;
    
    result = new Vector<String>();
    
    options = super.getOptions();
    for (i = 0; i < options.length; i++)
      result.add(options[i]);
    
    result.add("-B");
    result.add("" + getBase());

    return result.toArray(new String[result.size()]);
  }

  /**
   * Returns the distance/value of a given scale/level. I.e. the value of
   * base^i (e.g. 2^i).
   * 
   * @param s 		the level/scale
   * @return 		base^s
   */
  protected double dist_of_scale(int s) {
    return Math.pow(m_Base, s);
  }

  /**
   * Finds the scale/level of a given value. I.e. the "i" in base^i.
   * 
   * @param d 		the value whose scale/level is to be determined.
   * @return 		the scale/level of the given value.
   */
  protected int get_scale(double d) {
    return (int) Math.ceil(il2 * Math.log(d));
  }

  /**
   * Creates a new internal node for a given Instance/point p.
   * @param idx The index of the instance the node represents.
   * @return Newly created CoverTreeNode. 
   */
  protected CoverTreeNode new_node(Integer idx) { // const point &p)
    CoverTreeNode new_node = new CoverTreeNode();
    new_node.idx = idx;
    return new_node;
  }

  /**
   * Creates a new leaf node for a given Instance/point p.
   * @param idx The index of the instance this leaf node 
   * represents.
   * @return Newly created leaf CoverTreeNode.
   */
  protected CoverTreeNode new_leaf(Integer idx) { // (const point &p)
    CoverTreeNode new_leaf = new CoverTreeNode(idx, 0.0, 0.0, null, 0, 100);
    return new_leaf;
  }

  /**
   * Returns the max distance of the reference point p in current node to
   * it's children nodes.
   * @param v The stack of DistanceNode objects.
   * @return Distance of the furthest child.
   */
  protected double max_set(Stack<DistanceNode> v) { // rename to
                                                        // maxChildDist
    double max = 0.0;
    for (int i = 0; i < v.length; i++) {
      DistanceNode n = v.element(i);
      if (max < n.dist.element(n.dist.length - 1).floatValue()) { // v[i].dist.last())
        max = n.dist.element(n.dist.length - 1).floatValue(); // v[i].dist.last();
      }
    }
    return max;
  }

  /**
   * Splits a given point_set into near and far based on the given
   * scale/level. All points with distance > base^max_scale would be moved
   * to far set. In other words, all those points that are not covered by the 
   * next child ball of a point p (ball made of the same point p but of 
   * smaller radius at the next lower level) are removed from the supplied
   * current point_set and put into far_set.  
   * 
   * @param point_set The supplied set from which all far points 
   * would be removed.
   * @param far_set The set in which all far points having distance
   * > base^max_scale would be put into. 
   * @param max_scale The given scale based on which the distances
   * of points are judged to be far or near.   
   */
  protected void split(Stack<DistanceNode> point_set,
      Stack<DistanceNode> far_set, int max_scale) {
    int new_index = 0;
    double fmax = dist_of_scale(max_scale);
    for (int i = 0; i < point_set.length; i++) {
      DistanceNode n = point_set.element(i);
      if (n.dist.element(n.dist.length - 1).doubleValue() <= fmax) {
        point_set.set(new_index++, point_set.element(i));
      } else
        far_set.push(point_set.element(i)); // point_set[i]);
    }
    List<DistanceNode> l = new java.util.LinkedList<DistanceNode>();
    for (int i = 0; i < new_index; i++)
      l.add(point_set.element(i));
    //removing all and adding only the near points
    point_set.clear();
    point_set.addAll(l); // point_set.index=new_index;
  }

  /**
   * Moves all the points in point_set covered by (the ball of) new_point 
   * into new_point_set, based on the given scale/level.
   * 
   * @param point_set The supplied set of instances from which
   * all points covered by new_point will be removed.
   * @param new_point_set The set in which all points covered by
   * new_point will be put into.
   * @param new_point The given new point.
   * @param max_scale The scale based on which distances are 
   * judged (radius of cover ball is calculated).
   */
  protected void dist_split(Stack<DistanceNode> point_set,
      Stack<DistanceNode> new_point_set, 
      DistanceNode new_point, int max_scale) {
    int new_index = 0;
    double fmax = dist_of_scale(max_scale);
    for (int i = 0; i < point_set.length; i++) {
      double new_d =  Math.sqrt(m_DistanceFunction.distance(new_point.q(), 
	  	       point_set.element(i).q(), fmax*fmax));
      if (new_d <= fmax) {
        point_set.element(i).dist.push(new_d);
        new_point_set.push(point_set.element(i));
      } else
        point_set.set(new_index++, point_set.element(i));
    }
    List<DistanceNode> l = new java.util.LinkedList<DistanceNode>();
    for (int i = 0; i < new_index; i++)
      l.add(point_set.element(i));
    point_set.clear();
    point_set.addAll(l);
  }

  /**
   * Creates a cover tree recursively using batch insert method. 
   * 
   * @param p The index of the instance from which to create the
   * first node. All other points will be inserted beneath this node
   * for p.
   * @param max_scale The current scale/level where the node is to be
   * created (Also determines the radius of the cover balls created at 
   * this level).
   * @param top_scale The max scale in the whole tree.
   * @param point_set The set of unprocessed points from which child nodes
   * need to be created.  
   * @param consumed_set The set of processed points from which child
   * nodes have already been created. This would be used to find the 
   * radius of the cover ball of p. 
   * @return the node of cover tree created with p.
   */
  protected CoverTreeNode batch_insert(Integer p, int max_scale, // current
                                                                 // scale/level
      int top_scale, // max scale/level for this dataset
      Stack<DistanceNode> point_set, // set of points that are nearer to p
                                        // [will also contain returned unused
                                        // points]
      Stack<DistanceNode> consumed_set) // to return the set of points that have
                                        // been used to calc. max_dist to a
                                        // descendent
      // Stack<Stack<DistanceNode>> stack) //may not be needed
      {
    if (point_set.length == 0) {
      CoverTreeNode leaf = new_leaf(p);
      leaf.nodeid = m_NumNodes;
      m_NumNodes++; // incrementing node count
      m_NumLeaves++; // incrementing leaves count
      return leaf;
    } else {
      double max_dist = max_set(point_set); // O(|point_set|) the max dist
      // in point_set to point "p".
      int next_scale = Math.min(max_scale - 1, get_scale(max_dist));
      if (next_scale == Integer.MIN_VALUE) { // We have points with distance
        // 0. if max_dist is 0.
        Stack<CoverTreeNode> children = new Stack<CoverTreeNode>();
        CoverTreeNode leaf = new_leaf(p);
        leaf.nodeid = m_NumNodes;
        children.push(leaf);
        m_NumLeaves++;
        m_NumNodes++; // incrementing node and leaf count
        while (point_set.length > 0) {
          DistanceNode tmpnode = point_set.pop();
          leaf = new_leaf(tmpnode.idx);
          leaf.nodeid = m_NumNodes;
          children.push(leaf);
          m_NumLeaves++;
          m_NumNodes++; // incrementing node and leaf count
          consumed_set.push(tmpnode);
        }
        CoverTreeNode n = new_node(p); // make a new node out of p and assign
        // it the children.
        n.nodeid = m_NumNodes;
        m_NumNodes++; // incrementing node count
        n.scale = 100; // A magic number meant to be larger than all scales.
        n.max_dist = 0; // since all points have distance 0 to p
        n.num_children = children.length;
        n.children = children;
        return n;
      } else {
        Stack<DistanceNode> far = new Stack<DistanceNode>();
        split(point_set, far, max_scale); // O(|point_set|)

        CoverTreeNode child = batch_insert(p, next_scale, top_scale, point_set,
            consumed_set);

        if (point_set.length == 0) { // not creating any node in this
          // recursive call
          // push(stack,point_set);
          point_set.replaceAllBy(far); // point_set=far;
          return child;
        } else {
          CoverTreeNode n = new_node(p);
          n.nodeid = m_NumNodes;
          m_NumNodes++; // incrementing node count
          Stack<CoverTreeNode> children = new Stack<CoverTreeNode>();
          children.push(child);

          while (point_set.length != 0) { // O(|point_set| * num_children)
            Stack<DistanceNode> new_point_set = new Stack<DistanceNode>();
            Stack<DistanceNode> new_consumed_set = new Stack<DistanceNode>();
            DistanceNode tmpnode = point_set.pop();
            double new_dist = tmpnode.dist.last();
            consumed_set.push(tmpnode);

            // putting points closer to new_point into new_point_set (and
            // removing them from point_set)
            dist_split(point_set, new_point_set, tmpnode, max_scale); // O(|point_saet|)
            // putting points closer to new_point into new_point_set (and
            // removing them from far)
            dist_split(far, new_point_set, tmpnode, max_scale); // O(|far|)

            CoverTreeNode new_child = batch_insert(tmpnode.idx, next_scale,
                top_scale, new_point_set, new_consumed_set);
            new_child.parent_dist = new_dist;

            children.push(new_child);

            // putting the unused points from new_point_set back into
            // point_set and far
            double fmax = dist_of_scale(max_scale);
            tmpnode = null;
            for (int i = 0; i < new_point_set.length; i++) { // O(|new_point_set|)
              tmpnode = new_point_set.element(i);
              tmpnode.dist.pop();
              if (tmpnode.dist.last() <= fmax)
                point_set.push(tmpnode);
              else
                far.push(tmpnode);
            }
            // putting the points consumed while recursing for new_point
            // into consumed_set
            tmpnode = null;
            for (int i = 0; i < new_consumed_set.length; i++) { // O(|new_point_set|)
              tmpnode = new_consumed_set.element(i);
              tmpnode.dist.pop();
              consumed_set.push(tmpnode);
            }
          }// end while(point_size.size!=0)
          point_set.replaceAllBy(far); // point_set=far;
          n.scale = top_scale - max_scale;
          n.max_dist = max_set(consumed_set);
          n.num_children = children.length;
          n.children = children;
          return n;
        }// end else if(pointset!=0)
      }// end else if(next_scale != -214....
    }// end else if(pointset!=0)
  }

  /** 
   * Builds the tree on the given set of instances.
   * P.S.: For internal use only. Outside classes 
   * should call setInstances(). 
   * @param insts The instances on which to build 
   * the cover tree.
   * @throws Exception If the supplied set of 
   * Instances is empty, or if there are missing
   * values. 
   */
  protected void buildCoverTree(Instances insts) throws Exception {
    if (insts.numInstances() == 0)
      throw new Exception(
	  "CoverTree: Empty set of instances. Cannot build tree.");
    checkMissing(insts);
    if (m_EuclideanDistance == null)
      m_DistanceFunction = m_EuclideanDistance = new EuclideanDistance(insts);
    else
      m_EuclideanDistance.setInstances(insts);
    
    Stack<DistanceNode> point_set = new Stack<DistanceNode>();
    Stack<DistanceNode> consumed_set = new Stack<DistanceNode>();

    Instance point_p = insts.instance(0); int p_idx = 0;
    double max_dist=-1, dist=0.0; Instance max_q=point_p;
    
    for (int i = 1; i < insts.numInstances(); i++) {
      DistanceNode temp = new DistanceNode();
      temp.dist = new Stack<Double>();
      dist = Math.sqrt(m_DistanceFunction.distance(point_p, insts.instance(i), Double.POSITIVE_INFINITY));
      if(dist > max_dist) {
        max_dist = dist; max_q = insts.instance(i);
      }
      temp.dist.push(dist);
      temp.idx = i;
      point_set.push(temp);
    }
    
      max_dist = max_set(point_set);
      m_Root = batch_insert(p_idx, get_scale(max_dist), get_scale(max_dist),
                            point_set, consumed_set);
  }

/*********************************NNSearch related stuff********************/

  /**
   * A class for a heap to store the nearest k neighbours to an instance. 
   * The heap also takes care of cases where multiple neighbours are the same 
   * distance away.
   * i.e. the minimum size of the heap is k.
   * 
   * @author Ashraf M. Kibriya (amk14[at-the-rate]cs[dot]waikato[dot]ac[dot]nz)
   * @version $Revision: 5953 $
   */
  protected class MyHeap
    implements RevisionHandler {
    
    /** the heap. */
    MyHeapElement m_heap[] = null;
    
    /**
     * constructor.
     * @param maxSize   the maximum size of the heap
     */
    public MyHeap(int maxSize) {
      if((maxSize%2)==0)
        maxSize++;
      
      m_heap = new MyHeapElement[maxSize+1];
      m_heap[0] = new MyHeapElement(-1);
    }
    
    /**
     * returns the size of the heap.
     * @return the size
     */
    public int size() {
      return m_heap[0].index;
    }
    
    /**
     * peeks at the first element.
     * @return the first element
     */
    public MyHeapElement peek() {
      return m_heap[1];
    }
    
    /**
     * returns the first element and removes it from the heap.
     * @return the first element
     * @throws Exception  if no elements in heap
     */
    public MyHeapElement get() throws Exception  {
      if(m_heap[0].index==0)
        throw new Exception("No elements present in the heap");
      MyHeapElement r = m_heap[1];
      m_heap[1] = m_heap[m_heap[0].index];
      m_heap[0].index--;
      downheap();
      return r;
    }
    
    /**
     * adds the distance value to the heap.
     * 
     * @param d the distance value 
     * @throws Exception  if the heap gets too large
     */
    public void put(double d) throws Exception {
      if((m_heap[0].index+1)>(m_heap.length-1))
        throw new Exception("the number of elements cannot exceed the "+
        "initially set maximum limit");
      m_heap[0].index++;
      m_heap[m_heap[0].index] = new MyHeapElement(d);
      upheap();
    }
    
    /**
     * Puts an element by substituting it in place of 
     * the top most element.
     * 
     * @param d The distance value.
     * @throws Exception If distance is smaller than that of the head
     *         element.
     */
    public void putBySubstitute(double d) throws Exception {
      MyHeapElement head = get();
      put(d);
      if(head.distance == m_heap[1].distance) {
        putKthNearest(head.distance);
      }
      else if(head.distance > m_heap[1].distance) {
        m_KthNearest = null;
        m_KthNearestSize = 0;
        initSize = 10;
      }
      else if(head.distance < m_heap[1].distance) {
        throw new Exception("The substituted element is greater than the "+
        "head element. put() should have been called "+
        "in place of putBySubstitute()");
      }
    }
    
    /** the kth nearest ones. */
    MyHeapElement m_KthNearest[] = null;
    
    /** The number of kth nearest elements. */
    int m_KthNearestSize = 0;
    
    /** the initial size of the heap. */
    int initSize=10;
    
    /**
     * returns the number of k nearest.
     * 
     * @return the number of k nearest
     * @see     #m_KthNearestSize
     */
    public int noOfKthNearest() {
      return m_KthNearestSize;
    }
    
    /**
     * Stores kth nearest elements (if there are 
     * more than one).
     * @param d the distance 
     */
    public void putKthNearest(double d) {
      if(m_KthNearest==null) {
        m_KthNearest = new MyHeapElement[initSize];
      }
      if(m_KthNearestSize>=m_KthNearest.length) {
        initSize += initSize;
        MyHeapElement temp[] = new MyHeapElement[initSize];
        System.arraycopy(m_KthNearest, 0, temp, 0, m_KthNearest.length);
        m_KthNearest = temp;
      }
      m_KthNearest[m_KthNearestSize++] = new MyHeapElement(d);
    }
    
    /**
     * returns the kth nearest element or null if none there.
     * 
     * @return      the kth nearest element
     */
    public MyHeapElement getKthNearest() {
      if(m_KthNearestSize==0)
        return null;
      m_KthNearestSize--;
      return m_KthNearest[m_KthNearestSize];
    }
    
    /** 
     * performs upheap operation for the heap 
     * to maintian its properties. 
     */
    protected void upheap() {
      int i = m_heap[0].index;
      MyHeapElement temp;
      while( i > 1  && m_heap[i].distance>m_heap[i/2].distance) {
        temp = m_heap[i];
        m_heap[i] = m_heap[i/2];
        i = i/2;
        m_heap[i] = temp; //this is i/2 done here to avoid another division.
      }
    }
    
    /** 
     * performs downheap operation for the heap 
     * to maintian its properties. 
     */
    protected void downheap() {
      int i = 1;
      MyHeapElement temp;
      while( ( (2*i) <= m_heap[0].index &&
      m_heap[i].distance < m_heap[2*i].distance )
      ||
      ( (2*i+1) <= m_heap[0].index &&
      m_heap[i].distance < m_heap[2*i+1].distance) ) {
        if((2*i+1)<=m_heap[0].index) {
          if(m_heap[2*i].distance>m_heap[2*i+1].distance) {
            temp = m_heap[i];
            m_heap[i] = m_heap[2*i];
            i = 2*i;
            m_heap[i] = temp;
          }
          else {
            temp = m_heap[i];
            m_heap[i] = m_heap[2*i+1];
            i = 2*i+1;
            m_heap[i] = temp;
          }
        }
        else {
          temp = m_heap[i];
          m_heap[i] = m_heap[2*i];
          i = 2*i;
          m_heap[i] = temp;
        }
      }
    }
    
    /**
     * returns the total size.
     * 
     * @return      the total size
     */
    public int totalSize() {
      return size()+noOfKthNearest();
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
  
  /**
   * A class for storing data about a neighboring instance.
   * 
   * @author Ashraf M. Kibriya (amk14[at-the-rate]cs[dot]waikato[dot]ac[dot]nz)
   * @version $Revision: 5953 $
   */
  protected class MyHeapElement
    implements RevisionHandler {
    
    /** the distance. */
    public double distance;
    
    /** 
     * The index of this element. Also used as 
     * the size of the heap in the first element.
     */
    int index = 0;
    
    /**
     * constructor.
     * 
     * @param d   the distance
     */
    public MyHeapElement(double d) {
      distance = d;
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
  
  /**
   * stores a CoverTreeNode and its distance to the current query node.
   * 
   * @author Ashraf M. Kibriya (amk14[at-the-rate]cs[dot]waikato[dot]ac[dot]nz)
   * @version $Revision: 5953 $
   */
  private class d_node
    implements RevisionHandler {
    
    /** The distance of the node's point to the query point. */
    double dist;
    
    /** The node. */
    CoverTreeNode n;
    
    /** 
     * Constructor.
     * @param d The distance of the node to the query.
     * @param node The node. 
     */
    public d_node(double d, CoverTreeNode node) {
      dist = d;
      n = node;
    }
    
    /**
     * Returns the revision string.
     * 
     * @return		the revision
     */
    public String getRevision() {
      return RevisionUtils.extract("$Revision: 5953 $");
    }
  };

  /** 
   * Initializes a heap with k values of the the given upper_bound.
   * 
   * @param heap The heap to put values into.
   * @param upper_bound The value to put into heap (the value with 
   * which it should be initialized).
   * @param k The number of times upper_bound should be put into
   * heap for initialization.
   * @throws Exception If there is some problem in initializing 
   * the heap (if k &gt; size of the heap).
   */
  protected void setter(MyHeap heap, double upper_bound, final int k) throws Exception {
    if(heap.size()>0)
      heap.m_heap[0].index=0;

    while(heap.size() < k) {
      heap.put(upper_bound);
    }
  }

  /** 
   * Replaces the current top/max value in the heap with the new one.
   * The new max value should be &lt;= the old one.
   * 
   * @param upper_bound The heap.
   * @param new_bound The new value that should replace the old top one.
   * @throws Exception if the new value is greater than the old value.
   */
  protected void update(MyHeap upper_bound, double new_bound) throws Exception {
    upper_bound.putBySubstitute(new_bound);
  }
  
  /**
   * Returns a cover set for a given level/scale.
   * A cover set for a level consists of nodes whose 
   * Instances/centres are which are inside the query
   * ball at that level. If no cover set exists for the
   * given level (if it is the first time it is going 
   * to be used), than a new one is created.  
   * 
   * @param idx The level/scale for which the cover set 
   * is required.
   * @param cover_sets The covers sets. Consists of stack 
   * of a stack of d_node objects. 
   * @return The cover set for the given level/scale.
   */
  protected Stack<d_node> getCoverSet(int idx, Stack<Stack<d_node>> cover_sets) {
    if (cover_sets.length <= idx) {
      int i = cover_sets.length - 1;
      while (i < idx) {
        i++;
        Stack<d_node> new_cover_set = new Stack<d_node>();
        cover_sets.push(new_cover_set);
      }
    }
    return cover_sets.element(idx);
  }
  
  /**
   * Copies the contents of one zero set to the other. This
   * is required if we are going to inspect child of some query node 
   * (if the queries are given in batch in the form of a cover tree).
   * Only those nodes are copied to the new zero set that are inside
   * the query ball of query_chi.
   * P.S.: A zero set is a set of all leaf nodes that are found
   * to be inside the query ball.  
   *   
   * @param query_chi The child node of our query node that we are 
   * going to inspect. 
   * @param new_upper_k New heap that will store the distances of the
   * k NNs for query_chi.
   * @param zero_set The zero set of query_chi's parent that needs
   * to be copied.
   * @param new_zero_set The new zero set of query_chi where old zero
   * sets need to be copied into.
   * @throws Exception If there is some problem.
   */
  protected void copy_zero_set(CoverTreeNode query_chi, MyHeap new_upper_k, 
      			Stack<d_node> zero_set, Stack<d_node> new_zero_set) throws Exception {
    new_zero_set.clear();
    d_node ele;
    for (int i = 0; i < zero_set.length; i++) {
      ele = zero_set.element(i);
      double upper_dist = new_upper_k.peek().distance + query_chi.max_dist;
      if (shell(ele.dist, query_chi.parent_dist, upper_dist)) {
        double d = Math.sqrt(m_DistanceFunction.distance(query_chi.p(), ele.n
            .p(), upper_dist * upper_dist));
        if (m_TreeStats != null)
          m_TreeStats.incrPointCount();
        if (d <= upper_dist) {
          if (d < new_upper_k.peek().distance)
            update(new_upper_k, d);
          d_node temp = new d_node(d, ele.n);
          new_zero_set.push(temp);
          if (m_TreeStats != null)
            m_TreeStats.incrLeafCount();
        }//end if(d<newupperbound)
      }//end if(shell(...
    }//end for
  }
  

  /**
   * Copies the contents of one set of cover sets to the other. It
   * is required if we are going to inspect child of some query node 
   * (if the queries are given in batch in the form of a cover tree).
   * For each level, only those nodes are copied to the new set 
   * which are inside the query ball of query_chi at that level.
   * 
   * @param query_chi The child node of our query node that we are 
   * going to inspect. 
   * @param new_upper_k New heap that will store the distances of the
   * k NNs for query_chi.
   * @param cover_sets The cover_sets of query_chi's parent, which
   * need to be copied to new_cover_sets.
   * @param new_cover_sets The new set of cover_sets that need to
   * contain contents of cover_sets. 
   * @param current_scale The scale/level we are inspecting in our 
   * cover tree.
   * @param max_scale The maximum level so far possible in our 
   * search (this is only updated as we descend and a deeper
   * child is found inside the query ball).   
   * @throws Exception If there is problem.
   */
  protected void copy_cover_sets(CoverTreeNode query_chi, MyHeap new_upper_k,
      		Stack<Stack<d_node>> cover_sets,
      		Stack<Stack<d_node>> new_cover_sets,
      		int current_scale, int max_scale) throws Exception {
    new_cover_sets.clear();
    for (; current_scale <= max_scale; current_scale++) {
      d_node ele;
      Stack<d_node> cover_set_currentscale = getCoverSet(current_scale,
          cover_sets);
      for (int i = 0; i < cover_set_currentscale.length; i++) { // ; ele != end;
                                                                // ele++) {
        ele = cover_set_currentscale.element(i);
        double upper_dist = new_upper_k.peek().distance + query_chi.max_dist
            + ele.n.max_dist;
        if (shell(ele.dist, query_chi.parent_dist, upper_dist)) {
          double d = Math.sqrt(m_DistanceFunction.distance(query_chi.p(), ele.n
              .p(), upper_dist * upper_dist));
          if (m_TreeStats != null)
            m_TreeStats.incrPointCount();
          if (d <= upper_dist) {
            if (d < new_upper_k.peek().distance)
              update(new_upper_k, d);
            d_node temp = new d_node(d, ele.n);
            new_cover_sets.element(current_scale).push(temp);
            if (m_TreeStats != null)
              m_TreeStats.incrIntNodeCount();
          }// end if(d<=..
        }// end if(shell(...
      }// end for(coverset_i)
    }// end for(scales)
  }
  

  /**
   * Prints the given cover sets and zero set.
   * 
   * @param cover_sets The cover sets to print.
   * @param zero_set The zero set to print.  
   * @param current_scale The scale/level to start printing
   * the cover sets from. 
   * @param max_scale The max scale/level to print the cover
   * sets upto. 
   */
  void print_cover_sets(Stack<Stack<d_node>> cover_sets,
      Stack<d_node> zero_set, int current_scale, int max_scale) {
    d_node ele;
    println("cover set = ");
    for (; current_scale <= max_scale; current_scale++) {
      println("" + current_scale);
      for (int i = 0; i < cover_sets.element(current_scale).length; i++) {
        ele = cover_sets.element(current_scale).element(i);
        CoverTreeNode n = ele.n;
        println(n.p());
      }
    }
    println("infinity");
    for (int i = 0; i < zero_set.length; i++) {
      ele = zero_set.element(i);
      CoverTreeNode n = ele.n;
      println(n.p());
    }
  }
  
  
  /**
   * Swap two nodes in a cover set.
   * 
   * @param a The index first node.
   * @param b The index of second node.
   * @param cover_set The cover set in which the two nodes are.
   */

  protected void SWAP(int a, int b, Stack<d_node>cover_set) {				
    d_node tmp = cover_set.element(a);
    cover_set.set(a, cover_set.element(b));
    cover_set.set(b, tmp);
  }
  
  
  
  /**
   * Returns the difference of two given nodes distance to 
   * the query. It is used in half-sorting a cover set. 
   *   
   * @param p1 The index of first node.
   * @param p2 The index of second node.
   * @param cover_set The cover set containing the two given
   * nodes.
   * @return dist_to_query_of_p1 - dist_to_query_of_p2
   */
  
  protected double compare(final int p1, final int p2, Stack<d_node> cover_set) {
    return cover_set.element(p1).dist - cover_set.element(p2).dist;
  }
  
  /**
   * Half-sorts a cover set, so that nodes nearer to the query
   * are at the front. 
   * @param cover_set The cover set to sort.
   */

  protected void halfsort(Stack<d_node> cover_set) {
    if(cover_set.length <= 1)
      return;
    int start=0;
    int hi = cover_set.length-1;
    int right = hi;
    int left;
    
    while (right > start) {
      int mid = start + ((hi - start) >> 1);

      boolean jumpover = false;
      if (compare(mid, start, cover_set) < 0.0)
        SWAP(mid, start, cover_set);
      if (compare(hi, mid, cover_set) < 0.0)
        SWAP(mid, hi, cover_set);
      else
        jumpover = true;
      if (!jumpover && compare(mid, start, cover_set) < 0.0)
        SWAP(mid, start, cover_set);
      jump_over:
      ;

      left = start + 1;
      right = hi - 1;

      do {
        while (compare(left, mid, cover_set) < 0.0)
          left++;

        while (compare(mid, right, cover_set) < 0.0)
          right--;

        if (left < right) {
          SWAP(left, right, cover_set);
          if (mid == left)
            mid = right;
          else if (mid == right)
            mid = left;
          left++;
          right--;
        } else if (left == right) {
          left++;
          right--;
          break;
        }
      } while (left <= right);
      hi = right;
    }
  }

  /**
   * Function to check if a child node can be inside a query ball, 
   * without calculating the child node's distance to the query.
   * This further avoids unnecessary distance calculation. 
   *  
   * @param parent_query_dist The distance of parent to the query
   * @param child_parent_dist The distance of child to the parent.
   * @param upper_bound The distance to the query of the best kth 
   * NN found so far.
   * @return true If child can be inside the query ball.
   */
  protected boolean shell(double parent_query_dist, double child_parent_dist, double upper_bound) {
    return parent_query_dist - child_parent_dist <= upper_bound;
  }
  
  /**
   * This functions adds nodes for inspection at the next level during NN 
   * search. The internal nodes are added to one of the cover sets (at 
   * the level of the child node which is added) and leaf nodes are
   * added to the zero set.  
   *  
   * An optimization to consider:
   * Make all distance evaluations occur in descend.
   * 
   * Instead of passing a cover_set, pass a stack of cover sets.  The
   * last element holds d_nodes with your distance.  The next lower
   * element holds a d_node with the distance to your query parent,
   * next = query grand parent, etc..
   * 
   * Compute distances in the presence of the tighter upper bound.
   * @param query The query (in shape of a cover tree node, as we 
   * are doing batch searching).
   * @param upper_k Heap containing distances of best k-NNs found so 
   * far.
   * @param current_scale The current scale/level being looked at in 
   * the tree.
   * @param max_scale The max scale/level that has so far been looked
   * at.
   * @param cover_sets The cover sets of tree nodes for each level of 
   * our trees for.
   * @param zero_set The set containing leaf nodes.
   * @return A new max_scale, if we descend to a deeper level.
   * @throws Exception If there is some problem (in updating the 
   * heap upper_k).
   */
  protected int descend(final CoverTreeNode query, MyHeap upper_k,
      int current_scale, int max_scale, // amk14comment: make sure this gets
                                        // passed by reference in Java
      Stack<Stack<d_node>> cover_sets, // amk14comment: contains children in
                                        // set Q in paper
      Stack<d_node> zero_set) // amk14comment: zeroset contains the children at
                              // the lowest level i.e. -infinity
      throws Exception {
    d_node parent;
    Stack<d_node> cover_set_currentscale = getCoverSet(current_scale,
        cover_sets);
    for (int i = 0; i < cover_set_currentscale.length; i++) {
      parent = cover_set_currentscale.element(i);
      CoverTreeNode par = parent.n;
      double upper_dist = upper_k.peek().distance + query.max_dist
          + query.max_dist; // *upper_bound + query->max_dist + query->max_dist;
      if (parent.dist <= upper_dist + par.max_dist) {
        CoverTreeNode chi;
        if (par == m_Root && par.num_children == 0) // if our tree consists of
                                                    // only one root(which is
                                                    // also leaf) node
          chi = par;
        else
          chi = par.children.element(0);
        if (parent.dist <= upper_dist + chi.max_dist) { // amk14comment: looking
                                                        // at child_0 (which is
                                                        // the parent itself)
          if (chi.num_children > 0) {
            if (max_scale < chi.scale) {
              max_scale = chi.scale;
            }
            d_node temp = new d_node(parent.dist, chi);
            getCoverSet(chi.scale, cover_sets).push(temp);
            if (m_TreeStats != null)
              m_TreeStats.incrIntNodeCount();
          } else if (parent.dist <= upper_dist) {
            d_node temp = new d_node(parent.dist, chi);
            zero_set.push(temp);
            if (m_TreeStats != null)
              m_TreeStats.incrLeafCount();
          }
        }
        for (int c = 1; c < par.num_children; c++) {
          chi = par.children.element(c);
          double upper_chi = upper_k.peek().distance + chi.max_dist
              + query.max_dist + query.max_dist; // *upper_bound + chi.max_dist
                                                  // + query.max_dist +
                                                  // query.max_dist;
          if (shell(parent.dist, chi.parent_dist, upper_chi)) { // amk14comment:parent_query_dist
                                                                // -
                                                                // child_parent_dist
                                                                // <= upper_chi - if child can be 
                                                                // inside the shrunk query ball 
            // NOT the same as above parent->dist <= upper_dist + chi->max_dist
            double d = Math.sqrt(m_DistanceFunction.distance(query.p(),
                chi.p(), upper_chi * upper_chi, m_TreeStats));
            if (m_TreeStats != null)
              m_TreeStats.incrPointCount();
            if (d <= upper_chi) { //if child is inside the shrunk query ball
              if (d < upper_k.peek().distance) // *upper_bound)
                update(upper_k, d);
              if (chi.num_children > 0) {
                if (max_scale < chi.scale) {
                  max_scale = chi.scale;
                }
                d_node temp = new d_node(d, chi);
                getCoverSet(chi.scale, cover_sets).push(temp);
                if (m_TreeStats != null)
                  m_TreeStats.incrIntNodeCount();
              } else if (d <= upper_chi - chi.max_dist) {
                d_node temp = new d_node(d, chi);
                zero_set.push(temp);
                if (m_TreeStats != null)
                  m_TreeStats.incrLeafCount();
              }
            }//end if(d<=upper_chi)
          }//end if(shell(parent.dist,...
        }//end for(child_1 to n)
      }//end if(parent.dist<=upper_dist..
    }//end for(covers_sets[current_scale][i])
    return max_scale;
  }
  
  /**
   * Does a brute force NN search on the nodes in the given zero set.
   * A zero set might have some nodes added to it that were not k-NNs,
   * so need to do a brute-force to pick only the k-NNs (without 
   * calculating distances, as each node in the zero set already had 
   * its distance calculated to the query, which is stored with the
   * node).
   *  
   * @param k The k in kNN.
   * @param query The query. 
   * @param zero_set The zero set on which the brute force NN search
   * is performed.
   * @param upper_k The heap storing distances of k-NNs found during
   * the search.
   * @param results The returned k-NNs.
   * @throws Exception If there is somem problem.
   */
  protected void brute_nearest(final int k, final CoverTreeNode query,
      Stack<d_node> zero_set, MyHeap upper_k, Stack<NeighborList> results)
      throws Exception {
    if (query.num_children > 0) {
      Stack<d_node> new_zero_set = new Stack<d_node>();
      CoverTreeNode query_chi = query.children.element(0);
      brute_nearest(k, query_chi, zero_set, upper_k, results);
      MyHeap new_upper_k = new MyHeap(k);

      for (int i = 1; i < query.children.length; i++) {
        query_chi = query.children.element(i);
        setter(new_upper_k, upper_k.peek().distance + query_chi.parent_dist, k);
        copy_zero_set(query_chi, new_upper_k, zero_set, new_zero_set);
        brute_nearest(k, query_chi, new_zero_set, new_upper_k, results);
      }
    } else {
      NeighborList temp = new NeighborList(k);
      d_node ele;
      for (int i = 0; i < zero_set.length; i++) {
        ele = zero_set.element(i);
        if (ele.dist <= upper_k.peek().distance) {
          temp.insertSorted(ele.dist, ele.n.p()); // temp.push(ele.n.p());
        }
      }
      results.push(temp);
    }
  }
  
  /**
   * Performs a recursive k-NN search for a given batch of queries provided in the
   * form of a cover tree. P.S.: This function should not be called from outside. 
   * Outside classes should use kNearestNeighbours() instead.
   *  
   * @param k The number of NNs to find.
   * @param query_node The node of the query tree to start the search from.
   * @param cover_sets The set of sets that contains internal
   * nodes that were found to be inside the query ball at previous scales/levels
   * (intially there would be just the root node at root level).
   * @param zero_set The set that'll contain the leaf nodes that are found to
   * be inside the query ball.
   * @param current_scale The level/scale to do the search from (this value
   * would be used to inspect the cover set in the provided set of cover sets).
   * @param max_scale The max scale/level that has so far been inspected.
   * @param upper_k The heap containing distances of the best k-NNs found so
   * far (initialized to Double.POSITIVE_INFINITY).
   * @param results The list of returned k-NNs.
   * @throws Exception If there is some problem during the search.
   */
  protected void internal_batch_nearest_neighbor(final int k, 
      					final CoverTreeNode query_node,
      					Stack<Stack<d_node>> cover_sets,
      					Stack<d_node> zero_set,
      					int current_scale,
      					int max_scale,
      					MyHeap upper_k,
      					Stack<NeighborList> results) throws Exception {
    if (current_scale > max_scale) { // All remaining points are in the zero set.
      brute_nearest(k, query_node, zero_set, upper_k, results);
    } else {
      // Our query_node has too much scale. Reduce.
      if (query_node.scale <= current_scale && query_node.scale != 100) { // amk14comment:if j>=i in paper
        CoverTreeNode query_chi;
        Stack<d_node> new_zero_set = new Stack<d_node>();
        Stack<Stack<d_node>> new_cover_sets = new Stack<Stack<d_node>>();
        MyHeap new_upper_k = new MyHeap(k);

        for (int i = 1; i < query_node.num_children; i++) { //processing child_1 and onwards
          query_chi = query_node.children.element(i);
          setter(new_upper_k, upper_k.peek().distance + query_chi.parent_dist, k);
          //copy the zero set that satisfy a certain bound to the new zero set
          copy_zero_set(query_chi, new_upper_k, zero_set, new_zero_set);
          //copy the coversets[current_scale] nodes that satisfy a certain
          //bound to the new_cover_sets[current_scale]
          copy_cover_sets(query_chi, new_upper_k, cover_sets, new_cover_sets,
              current_scale, max_scale);
          //search for the query_node child in the nodes nearer to it.
          internal_batch_nearest_neighbor(k, query_chi, new_cover_sets,
              new_zero_set, current_scale, max_scale, new_upper_k, results);
        }
        new_cover_sets = null;
        new_zero_set = null;
        new_upper_k = null;
        // now doing child_0 //which is the parent itself, that's why we don't
        // need new_zero_set or new_cover_sets
        internal_batch_nearest_neighbor(k, query_node.children.element(0),
            cover_sets, zero_set, current_scale, max_scale, upper_k, results);
      } else { // reduce cover set scale -- amk14comment: if j<i in paper
        Stack<d_node> cover_set_i = getCoverSet(current_scale, cover_sets);
        // println("sorting");
        halfsort(cover_set_i);
        max_scale = descend(query_node, upper_k, current_scale, max_scale,
            cover_sets, zero_set);
        cover_set_i.clear();
        current_scale++;
        internal_batch_nearest_neighbor(k, query_node, cover_sets, zero_set,
            current_scale, max_scale, upper_k, results);
      }
    }
  }
  
  /**
   * Performs k-NN search for a batch of queries provided in the form
   * of a cover tree. P.S.: Outside classes should call 
   * kNearestNeighbours().
   * 
   * @param k The number of k-NNs to find.
   * @param tree_root The root of the cover tree on which k-NN search
   * is to be performed.
   * @param query_root The root of the cover tree consisting of queries. 
   * @param results The list of returned k-NNs.
   * @throws Exception If there is some problem during the search.
   */
  protected void batch_nearest_neighbor(final int k, CoverTreeNode tree_root, CoverTreeNode query_root, 
      			      Stack<NeighborList> results) throws Exception {
    //amk14comment: These contain the covering nodes at each level    
    Stack<Stack<d_node>> cover_sets = new Stack<Stack<d_node>>(100);  
    //amk14comment: These contain the nodes thought to be nearest at the leaf level
    Stack<d_node> zero_set = new Stack<d_node>(); 
    MyHeap upper_k = new MyHeap(k);
    //probably not needed //amk14comment:initializes the array to MAXFLOAT
    setter(upper_k, Double.POSITIVE_INFINITY, k); 

    // amk14comment:distance from top query point to top node point
    double treeroot_to_query_dist = Math.sqrt(m_DistanceFunction.distance(
        query_root.p(), tree_root.p(), Double.POSITIVE_INFINITY));
    // amk14comment:probably stores the kth smallest distances encountered so
    // far
    update(upper_k, treeroot_to_query_dist);

    d_node temp = new d_node(treeroot_to_query_dist, tree_root);
    getCoverSet(0, cover_sets).push(temp);

    // incrementing counts for the root node
    if (m_TreeStats != null) {
      m_TreeStats.incrPointCount();
      if (tree_root.num_children > 0)
        m_TreeStats.incrIntNodeCount();
      else
        m_TreeStats.incrLeafCount();
    }

    internal_batch_nearest_neighbor(k, query_root, cover_sets, zero_set, 0, 0,
        upper_k, results);
  }
  
  /**
   * Performs k-NN serach for a single given query/test Instance.
   * 
   * @param target The query/test instance.
   * @param k Number of k-NNs to find.
   * @return List of k-NNs.
   * @throws Exception If there is some problem during the search
   * for k-NNs.
   */
  protected NeighborList findKNearest(final Instance target, final int k) throws Exception {
    Stack<d_node> cover_set_current = new Stack<d_node>(),
    	           cover_set_next,
    	           zero_set = new Stack<d_node>();
    CoverTreeNode parent, child; d_node par;
    MyHeap upper_k = new MyHeap(k);    
    double d = Math.sqrt(m_DistanceFunction.distance(m_Root.p(), target, Double.POSITIVE_INFINITY, m_TreeStats)),
           upper_bound;
    cover_set_current.push(new d_node(d, m_Root));    
    setter(upper_k, Double.POSITIVE_INFINITY, k);
    this.update(upper_k, d);
    //updating stats for the root node
    if(m_TreeStats!=null) {
      	if(m_Root.num_children > 0)
      	  m_TreeStats.incrIntNodeCount();
      	else
      	  m_TreeStats.incrLeafCount();
      	m_TreeStats.incrPointCount();
    }
    
    //if root is the only node
    if(m_Root.num_children==0) {
      NeighborList list = new NeighborList(k);
      list.insertSorted(d, m_Root.p());
      return list;
    }
    //else
    while(cover_set_current.length>0) {
      cover_set_next = new Stack<d_node>();
      for(int i=0; i<cover_set_current.length; i++) {
	par = cover_set_current.element(i);
	parent = par.n;
	for(int c=0; c<parent.num_children; c++) {
	  child = parent.children.element(c);
	  upper_bound = upper_k.peek().distance;
	  if(c==0)
	    d = par.dist;
	  else {
	    d = upper_bound + child.max_dist;
	    d = Math.sqrt(m_DistanceFunction.distance(child.p(), target, d*d, m_TreeStats));
	      if(m_TreeStats!=null)
		m_TreeStats.incrPointCount();
	  }
	  if(d <= (upper_bound + child.max_dist)) {
	    if(c>0 && d < upper_bound) {
	      update(upper_k, d);
	    }
	    if(child.num_children > 0) {
	      cover_set_next.push(new d_node(d, child));
	      if(m_TreeStats!=null)
		m_TreeStats.incrIntNodeCount();
	    }
	    else if (d <= upper_bound){
	      zero_set.push(new d_node(d, child));
	      if(m_TreeStats!=null)
		m_TreeStats.incrLeafCount();
	    }
	  }
	} //end for current_set children
      } //end for current_set elements
      cover_set_current = cover_set_next;
    } //end while(curret_set not empty)
    
    NeighborList list = new NeighborList(k);
    d_node tmpnode;
    upper_bound = upper_k.peek().distance;      
    for(int i=0; i<zero_set.length; i++) {
      tmpnode = zero_set.element(i);
      if(tmpnode.dist <= upper_bound)
	list.insertSorted(tmpnode.dist, tmpnode.n.p());
    }
    
    if(list.currentLength()<=0)
      throw new Exception("Error: No neighbour found. This cannot happen");
    
    return list;
  }
  
/*********************************NNSearch related stuff above.********************/  

  /**
   * Returns k-NNs of a given target instance, from among the previously
   * supplied training instances (supplied through setInstances method)
   * P.S.: May return more than k-NNs if more one instances have
   * the same distance to the target as the kth NN.
   * 
   * @param target The instance for which k-NNs are required.
   * @param k The number of k-NNs to find.
   * @return The k-NN instances of the given target instance. 
   * @throws Exception If there is some problem find the k-NNs.
   */
  public Instances kNearestNeighbours(Instance target, int k) throws Exception {
    if(m_Stats!=null)
      m_Stats.searchStart();
    CoverTree querytree = new CoverTree();
    Instances insts = new Instances(m_Instances, 0);
    insts.add(target);
    querytree.setInstances(insts);
    Stack<NeighborList> result = new Stack<NeighborList>();
    batch_nearest_neighbor(k, this.m_Root, querytree.m_Root, result);
    if(m_Stats!=null)
      m_Stats.searchFinish();

    insts = new Instances(m_Instances, 0);
    NeighborNode node = result.element(0).getFirst();
    m_DistanceList = new double[result.element(0).currentLength()];
    int i=0;
    while(node != null) {
      insts.add(node.m_Instance);
      m_DistanceList[i] = node.m_Distance;
      i++; node = node.m_Next;
    }
    return insts;
  }
  
  /**
   * Returns the NN instance of a given target instance, from among
   * the previously supplied training instances.
   * 
   * @param target The instance for which NN is required.
   * @throws Exception If there is some problem finding the nearest
   * neighbour.
   * @return The NN instance of the target instance.
   */
  public Instance nearestNeighbour(Instance target) throws Exception {
    return kNearestNeighbours(target, 1).instance(0);
  }

  /**
   * Returns the distances of the (k)-NN(s) found earlier
   * by kNearestNeighbours()/nearestNeighbour().
   * 
   * @throws Exception If the tree hasn't been built (by calling 
   * setInstances()), or none of kNearestNeighbours() or 
   * nearestNeighbour() has been called before. 
   * @return The distances (in the same order) of the k-NNs. 
   */
  public double[] getDistances() throws Exception {
    if(m_Instances==null || m_DistanceList==null)
      throw new Exception("The tree has not been supplied with a set of " +
	  		  "instances or getDistances() has been called " +
      			  "before calling kNearestNeighbours().");
    return m_DistanceList;
  }
  
  /**
   * Checks if there is any instance with missing values. Throws an
   * exception if there is, as KDTree does not handle missing values.
   * 
   * @param instances 	the instances to check
   * @throws Exception 	if missing values are encountered
   */
  protected void checkMissing(Instances instances) throws Exception {
    for (int i = 0; i < instances.numInstances(); i++) {
      Instance ins = instances.instance(i);
      for (int j = 0; j < ins.numValues(); j++) {
	if (ins.index(j) != ins.classIndex())
	  if (ins.isMissingSparse(j)) {
	    throw new Exception("ERROR: KDTree can not deal with missing "
		+ "values. Please run ReplaceMissingValues filter "
		+ "on the dataset before passing it on to the KDTree.");
	  }
      }
    }
  }

  /**
   * Builds the Cover Tree on the given set of instances.
   * 
   * @param instances The insts on which the Cover Tree is to be 
   * built. 
   * @throws Exception If some error occurs while 
   * building the Cover Tree
   */
  public void setInstances(Instances instances) throws Exception {
    super.setInstances(instances);
    buildCoverTree(instances);
  }

  /** 
   * Adds an instance to the cover tree. 
   * P.S.: The current version doesn't allow
   * addition of instances after batch construction.
   * 
   * @param ins The instance to add.
   * @throws Exception Alway throws this, as current 
   * implementation doesn't allow addition of instances 
   * after building.
   */
  public void update(Instance ins) throws Exception {
    throw new Exception("BottomUpConstruction method does not allow addition " +
    "of new Instances.");
  }

  /** 
   * Adds the given instance info. This implementation updates only the 
   * range datastructures of the EuclideanDistance. Nothing is 
   * required to be updated in the built Cover Tree.
   * 
   * @param ins 	The instance to add the information of. Usually this is
   * 			the test instance supplied to update the range of 
   * 			attributes in the distance function.
   */
  public void addInstanceInfo(Instance ins) {
    if(m_Instances!=null) {
      try {
      m_DistanceFunction.update(ins);
      } catch(Exception ex) { ex.printStackTrace(); }
    }
    else 
      if(m_Instances==null)
	      throw new IllegalStateException("No instances supplied yet. Cannot update without"+
	                          "supplying a set of instances first.");
  }
  
  /**
   * Sets the distance function to use for nearest neighbour search.
   * Currently only EuclideanDistance is supported.
   * 
   * @param df 		the distance function to use 
   * @throws Exception 	if not EuclideanDistance
   */
  public void setDistanceFunction(DistanceFunction df) throws Exception {
    if (!(df instanceof EuclideanDistance))
      throw new Exception("CoverTree currently only works with "
	  + "EuclideanDistanceFunction.");
    m_DistanceFunction = m_EuclideanDistance = (EuclideanDistance) df;
  }
  
  /**
   * Returns the tip text for this property.
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String baseTipText() {
    return "The base for the expansion constant.";
  }

  /**
   * Returns the base in use for expansion constant.
   * 
   * @return base 	currently in use.
   */
  public double getBase() {
    return m_Base;
  }
  
  /**
   * Sets the base to use for expansion constant.
   * The 2 in 2^i in the paper.
   * 
   * @param b 		the new base;
   */
  public void setBase(double b) {
    m_Base = b;
  }
  
  /**
   * Returns the size of the tree. 
   * (number of internal nodes + number of leaves)
   * 
   * @return 		the size of the tree
   */
  public double measureTreeSize() {
    return m_NumNodes;
  }
  
  /**
   * Returns the number of leaves.
   * 
   * @return 		the number of leaves
   */
  public double measureNumLeaves() {
    return m_NumLeaves;
  }
  
  /**
   * Returns the depth of the tree.
   * 
   * @return 		the number of rules
   */
  public double measureMaxDepth() {
    return m_MaxDepth;
  }
    
  /**
   * Returns an enumeration of the additional measure names.
   * 
   * @return 		an enumeration of the measure names
   */
  public Enumeration enumerateMeasures() {
    Vector<String> newVector = new Vector<String>();
    newVector.addElement("measureTreeSize");
    newVector.addElement("measureNumLeaves");
    newVector.addElement("measureMaxDepth");
    if(m_Stats!=null) {
      for(Enumeration e = m_Stats.enumerateMeasures(); e.hasMoreElements();) {
        newVector.addElement((String)e.nextElement());
      }
    }
    return newVector.elements();
  }
  
  /**
   * Returns the value of the named measure.
   * 
   * @param additionalMeasureName 	the name of the measure to query for 
   * 					its value
   * @return 				the value of the named measure
   * @throws IllegalArgumentException 	if the named measure is not supported
   */
  public double getMeasure(String additionalMeasureName) {
    if (additionalMeasureName.compareToIgnoreCase("measureMaxDepth") == 0) {
      return measureMaxDepth();
    } else if (additionalMeasureName.compareToIgnoreCase("measureTreeSize") == 0) {
      return measureTreeSize();
    } else if (additionalMeasureName.compareToIgnoreCase("measureNumLeaves") == 0) {
      return measureNumLeaves();
    } else if(m_Stats!=null) {
      return m_Stats.getMeasure(additionalMeasureName);
    } else {
      throw new IllegalArgumentException(additionalMeasureName 
			  + " not supported (KDTree)");
    }
  }
  
  /********Utility print functions.****** */
  /**
   * Prints a string to stdout. 
   * 
   * @param s The string to print.
   */
  protected static void print(String s) {
    System.out.print(s);
  }

  /** 
   * Prints a string to stdout followed by 
   * newline.
   * 
   * @param s The string to print. 
   */
  protected static void println(String s) {
    System.out.println(s);
  }

  /** 
   * Prints an object to stdout.
   * 
   * @param o The object to print. 
   */
  protected static void print(Object o) {
    System.out.print(o);
  }

  /** 
   * Prints an object to stdout followed by 
   * newline.
   * 
   * @param o The object to print.  
   */
  protected static void println(Object o) {
    System.out.println(o);
  }

  /** 
   * Prints the specified number of spaces.
   * 
   * @param s The number of space characters to print.  
   */
  protected static void print_space(int s) {
    for (int i = 0; i < s; i++)
      System.out.print(" ");
  }

  /**
   * Prints a cover tree starting from the given node.
   * 
   * @param depth The depth of top_node.
   * @param top_node The node to start printing from. 
   */
  protected static void print(int depth, CoverTreeNode top_node) {
    print_space(depth);
    println(top_node.p());
    if (top_node.num_children > 0) {
      print_space(depth);
      print("scale = " + top_node.scale + "\n");
      print_space(depth);
      print("num children = " + top_node.num_children + "\n");
      System.out.flush();
      for (int i = 0; i < top_node.num_children; i++)
        print(depth + 1, top_node.children.element(i)); // top_node.children[i]);
    }
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
   * Method for testing the class from command line. 
   * 
   * @param args The supplied command line arguments.
   */
  public static void main(String[] args) {
    if (args.length != 1) {
      System.err.println("Usage: CoverTree <ARFF file>");
      System.exit(-1);
    }
    try {
      Instances insts = null;
      if (args[0].endsWith(".csv")) {
        CSVLoader csv = new CSVLoader();
        csv.setFile(new File(args[0]));
        insts = csv.getDataSet();
      } else {
        insts = new Instances(new BufferedReader(new FileReader(args[0])));
      }

      CoverTree tree = new CoverTree();
      tree.setInstances(insts);
      print("Created data tree:\n");
      print(0, tree.m_Root);
      println("");
    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }
}
