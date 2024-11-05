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
 *    BestFirst.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.attributeSelection;

import weka.core.FastVector;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.Range;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.core.SelectedTag;
import weka.core.Tag;
import weka.core.Utils;

import java.io.Serializable;
import java.util.BitSet;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Vector;

/** 
 <!-- globalinfo-start -->
 * BestFirst:<br/>
 * <br/>
 * Searches the space of attribute subsets by greedy hillclimbing augmented with a backtracking facility. Setting the number of consecutive non-improving nodes allowed controls the level of backtracking done. Best first may start with the empty set of attributes and search forward, or start with the full set of attributes and search backward, or start at any point and search in both directions (by considering all possible single attribute additions and deletions at a given point).<br/>
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -P &lt;start set&gt;
 *  Specify a starting set of attributes.
 *  Eg. 1,3,5-7.</pre>
 * 
 * <pre> -D &lt;0 = backward | 1 = forward | 2 = bi-directional&gt;
 *  Direction of search. (default = 1).</pre>
 * 
 * <pre> -N &lt;num&gt;
 *  Number of non-improving nodes to
 *  consider before terminating search.</pre>
 * 
 * <pre> -S &lt;num&gt;
 *  Size of lookup cache for evaluated subsets.
 *  Expressed as a multiple of the number of
 *  attributes in the data set. (default = 1)</pre>
 * 
 <!-- options-end -->
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 *         Martin Guetlein (cashing merit of expanded nodes) 
 * @version $Revision: 1.29 $
 */
public class BestFirst 
  extends ASSearch 
  implements OptionHandler, StartSetHandler {
  
  /** for serialization */
  static final long serialVersionUID = 7841338689536821867L;

  // Inner classes
  /**
   * Class for a node in a linked list. Used in best first search.
   * @author Mark Hall (mhall@cs.waikato.ac.nz)
   **/
  public class Link2 
    implements Serializable, RevisionHandler {

    /** for serialization */
    static final long serialVersionUID = -8236598311516351420L;
    
    /*    BitSet group; */
    Object [] m_data;
    double m_merit;


    /** 
     * Constructor
     */
    public Link2 (Object [] data, double mer) {
      //      group = (BitSet)gr.clone();
      m_data = data;
      m_merit = mer;
    }


    /** Get a group */
    public Object [] getData () {
      return  m_data;
    }


    public String toString () {
      return  ("Node: " + m_data.toString() + "  " + m_merit);
    }
    
    /**
     * Returns the revision string.
     * 
     * @return		the revision
     */
    public String getRevision() {
      return RevisionUtils.extract("$Revision: 1.29 $");
    }
  }


  /**
   * Class for handling a linked list. Used in best first search.
   * Extends the Vector class.
   * @author Mark Hall (mhall@cs.waikato.ac.nz)
   **/
  public class LinkedList2
    extends FastVector {
    
    /** for serialization */
    static final long serialVersionUID = 3250538292330398929L;
    
    /** Max number of elements in the list */
    int m_MaxSize;

    // ================
    // Public methods
    // ================
    public LinkedList2 (int sz) {
      super();
      m_MaxSize = sz;
    }


    /**
     * removes an element (Link) at a specific index from the list.
     * @param index the index of the element to be removed.
     **/
    public void removeLinkAt (int index)
      throws Exception {
      
      if ((index >= 0) && (index < size())) {
        removeElementAt(index);
      }
      else {
        throw  new Exception("index out of range (removeLinkAt)");
      }
    }


    /**
     * returns the element (Link) at a specific index from the list.
     * @param index the index of the element to be returned.
     **/
    public Link2 getLinkAt (int index)
      throws Exception {
      
      if (size() == 0) {
        throw  new Exception("List is empty (getLinkAt)");
      }
      else {if ((index >= 0) && (index < size())) {
        return  ((Link2)(elementAt(index)));
      }
      else {
        throw  new Exception("index out of range (getLinkAt)");
      }
      }
    }


    /**
     * adds an element (Link) to the list.
     * @param data the attribute set specification
     * @param mer the "merit" of this attribute set
     **/
    public void addToList (Object [] data, double mer)
      throws Exception {
      Link2 newL = new Link2(data, mer);

      if (size() == 0) {
	addElement(newL);
      }
      else {if (mer > ((Link2)(firstElement())).m_merit) {
	if (size() == m_MaxSize) {
	  removeLinkAt(m_MaxSize - 1);
	}

	//----------
	insertElementAt(newL, 0);
      }
      else {
	int i = 0;
	int size = size();
	boolean done = false;

	//------------
	// don't insert if list contains max elements an this
	// is worst than the last
	if ((size == m_MaxSize) && (mer <= ((Link2)(lastElement())).m_merit)) {

	}
	//---------------
	else {
	  while ((!done) && (i < size)) {
	    if (mer > ((Link2)(elementAt(i))).m_merit) {
	      if (size == m_MaxSize) {
		removeLinkAt(m_MaxSize - 1);
	      }

	      // ---------------------
	      insertElementAt(newL, i);
	      done = true;
	    }
	    else {if (i == size - 1) {
	      addElement(newL);
	      done = true;
	    }
	    else {
	      i++;
	    }
	    }
	  }
	}
      }
      }
    }
    
    /**
     * Returns the revision string.
     * 
     * @return		the revision
     */
    public String getRevision() {
      return RevisionUtils.extract("$Revision: 1.29 $");
    }
  }

  // member variables
  /** maximum number of stale nodes before terminating search */
  protected int m_maxStale;

  /** 0 == backward search, 1 == forward search, 2 == bidirectional */
  protected int m_searchDirection;

  /** search direction: backward */
  protected static final int SELECTION_BACKWARD = 0;
  /** search direction: forward */
  protected static final int SELECTION_FORWARD = 1;
  /** search direction: bidirectional */
  protected static final int SELECTION_BIDIRECTIONAL = 2;
  /** search directions */
  public static final Tag [] TAGS_SELECTION = {
    new Tag(SELECTION_BACKWARD, "Backward"),
    new Tag(SELECTION_FORWARD, "Forward"),
    new Tag(SELECTION_BIDIRECTIONAL, "Bi-directional"),
  };

  /** holds an array of starting attributes */
  protected int[] m_starting;

  /** holds the start set for the search as a Range */
  protected Range m_startRange;

  /** does the data have a class */
  protected boolean m_hasClass;

  /** holds the class index */
  protected int m_classIndex;

  /** number of attributes in the data */
  protected int m_numAttribs;

  /** total number of subsets evaluated during a search */
  protected int m_totalEvals;

  /** for debugging */
  protected boolean m_debug;

  /** holds the merit of the best subset found */
  protected double m_bestMerit;

  /** holds the maximum size of the lookup cache for evaluated subsets */
  protected int m_cacheSize;
  
  /**
   * Returns a string describing this search method
   * @return a description of the search method suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return "BestFirst:\n\n"
      +"Searches the space of attribute subsets by greedy hillclimbing "
      +"augmented with a backtracking facility. Setting the number of "
      +"consecutive non-improving nodes allowed controls the level of "
      +"backtracking done. Best first may start with the empty set of "
      +"attributes and search forward, or start with the full set of "
      +"attributes and search backward, or start at any point and search "
      +"in both directions (by considering all possible single attribute "
      +"additions and deletions at a given point).\n";
  }

  /** 
   *Constructor 
   */
  public BestFirst () {
    resetOptions();
  }

  /**
   * Returns an enumeration describing the available options.
   * @return an enumeration of all the available options.
   *
   **/
  public Enumeration listOptions () {
    Vector newVector = new Vector(4);
    
    newVector.addElement(new Option("\tSpecify a starting set of attributes." 
				    + "\n\tEg. 1,3,5-7."
				    ,"P",1
				    , "-P <start set>"));
    newVector.addElement(new Option("\tDirection of search. (default = 1)."
				    , "D", 1
				    , "-D <0 = backward | 1 = forward " 
				    + "| 2 = bi-directional>"));
    newVector.addElement(new Option("\tNumber of non-improving nodes to" 
				    + "\n\tconsider before terminating search."
				    , "N", 1, "-N <num>"));
    newVector.addElement(new Option("\tSize of lookup cache for evaluated subsets."
				    +"\n\tExpressed as a multiple of the number of"
				    +"\n\tattributes in the data set. (default = 1)",
				    "S", 1, "-S <num>"));
				    
    return  newVector.elements();
  }


  /**
   * Parses a given list of options. <p/>
   *
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -P &lt;start set&gt;
   *  Specify a starting set of attributes.
   *  Eg. 1,3,5-7.</pre>
   * 
   * <pre> -D &lt;0 = backward | 1 = forward | 2 = bi-directional&gt;
   *  Direction of search. (default = 1).</pre>
   * 
   * <pre> -N &lt;num&gt;
   *  Number of non-improving nodes to
   *  consider before terminating search.</pre>
   * 
   * <pre> -S &lt;num&gt;
   *  Size of lookup cache for evaluated subsets.
   *  Expressed as a multiple of the number of
   *  attributes in the data set. (default = 1)</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   *
   **/
  public void setOptions (String[] options)
    throws Exception {
    String optionString;
    resetOptions();

    optionString = Utils.getOption('P', options);
    if (optionString.length() != 0) {
      setStartSet(optionString);
    }

    optionString = Utils.getOption('D', options);

    if (optionString.length() != 0) {
      setDirection(new SelectedTag(Integer.parseInt(optionString),
				   TAGS_SELECTION));
    } else {
      setDirection(new SelectedTag(SELECTION_FORWARD, TAGS_SELECTION));
    }

    optionString = Utils.getOption('N', options);

    if (optionString.length() != 0) {
      setSearchTermination(Integer.parseInt(optionString));
    }

    optionString = Utils.getOption('S', options);
    if (optionString.length() != 0) {
      setLookupCacheSize(Integer.parseInt(optionString));
    }

    m_debug = Utils.getFlag('Z', options);
  }

  /**
   * Set the maximum size of the evaluated subset cache (hashtable). This is
   * expressed as a multiplier for the number of attributes in the data set.
   * (default = 1).
   *
   * @param size the maximum size of the hashtable
   */
  public void setLookupCacheSize(int size) {
    if (size >= 0) {
      m_cacheSize = size;
    }
  }

  /**
   * Return the maximum size of the evaluated subset cache (expressed as a multiplier
   * for the number of attributes in a data set.
   *
   * @return the maximum size of the hashtable.
   */
  public int getLookupCacheSize() {
    return m_cacheSize;
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String lookupCacheSizeTipText() {
    return "Set the maximum size of the lookup cache of evaluated subsets. This is "
      +"expressed as a multiplier of the number of attributes in the data set. "
      +"(default = 1).";
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String startSetTipText() {
    return "Set the start point for the search. This is specified as a comma "
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
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String searchTerminationTipText() {
    return "Set the amount of backtracking. Specify the number of ";
  }

  /**
   * Set the numnber of non-improving nodes to consider before terminating
   * search.
   *
   * @param t the number of non-improving nodes
   * @throws Exception if t is less than 1
   */
  public void setSearchTermination (int t)
    throws Exception {
    if (t < 1) {
      throw  new Exception("Value of -N must be > 0.");
    }

    m_maxStale = t;
  }


  /**
   * Get the termination criterion (number of non-improving nodes).
   *
   * @return the number of non-improving nodes
   */
  public int getSearchTermination () {
    return  m_maxStale;
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String directionTipText() {
    return "Set the direction of the search.";
  }

  /**
   * Set the search direction
   *
   * @param d the direction of the search
   */
  public void setDirection (SelectedTag d) {
    
    if (d.getTags() == TAGS_SELECTION) {
      m_searchDirection = d.getSelectedTag().getID();
    }
  }


  /**
   * Get the search direction
   *
   * @return the direction of the search
   */
  public SelectedTag getDirection () {

    return new SelectedTag(m_searchDirection, TAGS_SELECTION);
  }


  /**
   * Gets the current settings of BestFirst.
   * @return an array of strings suitable for passing to setOptions()
   */
  public String[] getOptions () {
    String[] options = new String[6];
    int current = 0;

    if (!(getStartSet().equals(""))) {
      options[current++] = "-P";
      options[current++] = ""+startSetToString();
    }
    options[current++] = "-D";
    options[current++] = "" + m_searchDirection;
    options[current++] = "-N";
    options[current++] = "" + m_maxStale;

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
   * returns a description of the search as a String
   * @return a description of the search
   */
  public String toString () {
    StringBuffer BfString = new StringBuffer();
    BfString.append("\tBest first.\n\tStart set: ");

    if (m_starting == null) {
      BfString.append("no attributes\n");
    }
    else {
      BfString.append(startSetToString()+"\n");
    }

    BfString.append("\tSearch direction: ");

    if (m_searchDirection == SELECTION_BACKWARD) {
      BfString.append("backward\n");
    }
    else {if (m_searchDirection == SELECTION_FORWARD) {
      BfString.append("forward\n");
    }
    else {
      BfString.append("bi-directional\n");
    }
    }

    BfString.append("\tStale search after " 
		    + m_maxStale + " node expansions\n");
    BfString.append("\tTotal number of subsets evaluated: " 
		    + m_totalEvals + "\n");
    BfString.append("\tMerit of best subset found: "
		    +Utils.doubleToString(Math.abs(m_bestMerit),8,3)+"\n");
    return  BfString.toString();
  }


  protected void printGroup (BitSet tt, int numAttribs) {
    int i;

    for (i = 0; i < numAttribs; i++) {
      if (tt.get(i) == true) {
	System.out.print((i + 1) + " ");
      }
    }

    System.out.println();
  }


  /**
   * Searches the attribute subset space by best first search
   *
   * @param ASEval the attribute evaluator to guide the search
   * @param data the training instances.
   * @return an array (not necessarily ordered) of selected attribute indexes
   * @throws Exception if the search can't be completed
   */
  public int[] search (ASEvaluation ASEval, Instances data)
    throws Exception {
    m_totalEvals = 0;
    if (!(ASEval instanceof SubsetEvaluator)) {
      throw  new Exception(ASEval.getClass().getName() 
			   + " is not a " 
			   + "Subset evaluator!");
    }

    if (ASEval instanceof UnsupervisedSubsetEvaluator) {
      m_hasClass = false;
    } else {
      m_hasClass = true;
      m_classIndex = data.classIndex();
    }

    SubsetEvaluator ASEvaluator = (SubsetEvaluator)ASEval;
    m_numAttribs = data.numAttributes();
    int i, j;
    int best_size = 0;
    int size = 0;
    int done;
    int sd = m_searchDirection;
    BitSet best_group, temp_group;
    int stale;
    double best_merit;
    double merit;
    boolean z;
    boolean added;
    Link2 tl;
    Hashtable lookup = new Hashtable(m_cacheSize * m_numAttribs);
    int insertCount = 0;
    int cacheHits = 0;
    LinkedList2 bfList = new LinkedList2(m_maxStale);
    best_merit = -Double.MAX_VALUE;
    stale = 0;
    best_group = new BitSet(m_numAttribs);

    m_startRange.setUpper(m_numAttribs-1);
    if (!(getStartSet().equals(""))) {
      m_starting = m_startRange.getSelection();
    }
    // If a starting subset has been supplied, then initialise the bitset
    if (m_starting != null) {
      for (i = 0; i < m_starting.length; i++) {
	if ((m_starting[i]) != m_classIndex) {
	  best_group.set(m_starting[i]);
	}
      }

      best_size = m_starting.length;
      m_totalEvals++;
    } else {
      if (m_searchDirection == SELECTION_BACKWARD) {
	setStartSet("1-last");
	m_starting = new int[m_numAttribs];

	// init initial subset to all attributes
	for (i = 0, j = 0; i < m_numAttribs; i++) {
	  if (i != m_classIndex) {
	    best_group.set(i);
	    m_starting[j++] = i;
	  }
	}

	best_size = m_numAttribs - 1;
	m_totalEvals++;
      }
    }

    // evaluate the initial subset
    best_merit = ASEvaluator.evaluateSubset(best_group);
    // add the initial group to the list and the hash table
    Object [] best = new Object[1];
    best[0] = best_group.clone();
    bfList.addToList(best, best_merit);
    BitSet tt = (BitSet)best_group.clone();
    String hashC = tt.toString();
    lookup.put(hashC, new Double(best_merit));

    while (stale < m_maxStale) {
      added = false;

      if (m_searchDirection == SELECTION_BIDIRECTIONAL) {
	// bi-directional search
        done = 2;
        sd = SELECTION_FORWARD;
      } else {
	done = 1;
      }

      // finished search?
      if (bfList.size() == 0) {
	stale = m_maxStale;
	break;
      }

      // copy the attribute set at the head of the list
      tl = bfList.getLinkAt(0);
      temp_group = (BitSet)(tl.getData()[0]);
      temp_group = (BitSet)temp_group.clone();
      // remove the head of the list
      bfList.removeLinkAt(0);
      // count the number of bits set (attributes)
      int kk;

      for (kk = 0, size = 0; kk < m_numAttribs; kk++) {
	if (temp_group.get(kk)) {
	  size++;
	}
      }

      do {
	for (i = 0; i < m_numAttribs; i++) {
	  if (sd == SELECTION_FORWARD) {
	    z = ((i != m_classIndex) && (!temp_group.get(i)));
	  } else {
	    z = ((i != m_classIndex) && (temp_group.get(i)));
	  }
          
	  if (z) {
	    // set the bit (attribute to add/delete)
	    if (sd == SELECTION_FORWARD) {
	      temp_group.set(i);
	      size++;
	    } else {
	      temp_group.clear(i);
	      size--;
	    }

	    /* if this subset has been seen before, then it is already 
	       in the list (or has been fully expanded) */
	    tt = (BitSet)temp_group.clone();
	    hashC = tt.toString();
	    
	    if (lookup.containsKey(hashC) == false) {
	      merit = ASEvaluator.evaluateSubset(temp_group);
	      m_totalEvals++;
	      
	      // insert this one in the hashtable
	      if (insertCount > m_cacheSize * m_numAttribs) {
		lookup = new Hashtable(m_cacheSize * m_numAttribs);
		insertCount = 0;
	      }
	      hashC = tt.toString();
    	      lookup.put(hashC, new Double(merit));
    	      insertCount++;
	    } else {
	      merit = ((Double)lookup.get(hashC)).doubleValue();
	      cacheHits++;  
	    }
	    
	    // insert this one in the list
	    Object[] add = new Object[1];
	    add[0] = tt.clone();
	    bfList.addToList(add, merit);
	    
	    if (m_debug) {
	      System.out.print("Group: ");
	      printGroup(tt, m_numAttribs);
	      System.out.println("Merit: " + merit);
	    }

	    // is this better than the best?
	    if (sd == SELECTION_FORWARD) {
	      z = ((merit - best_merit) > 0.00001);
	    } else {
	      if (merit == best_merit) {
		z = (size < best_size);
	      } else {
		z = (merit >  best_merit);
	      } 
	    }

	    if (z) {
	      added = true;
	      stale = 0;
	      best_merit = merit;
	      //		best_size = (size + best_size);
	      best_size = size;
	      best_group = (BitSet)(temp_group.clone());
	    }

	    // unset this addition(deletion)
	    if (sd == SELECTION_FORWARD) {
	      temp_group.clear(i);
	      size--;
	    } else {
	      temp_group.set(i);
	      size++;
	    }
	  }
	}

	if (done == 2) {
	  sd = SELECTION_BACKWARD;
	}

	done--;
      } while (done > 0);

      /* if we haven't added a new attribute subset then full expansion 
	 of this node hasen't resulted in anything better */
      if (!added) {
	stale++;
      }
    }

    m_bestMerit = best_merit;
    return  attributeList(best_group);
  }


  /**
   * Reset options to default values
   */
  protected void resetOptions () {
    m_maxStale = 5;
    m_searchDirection = SELECTION_FORWARD;
    m_starting = null;
    m_startRange = new Range();
    m_classIndex = -1;
    m_totalEvals = 0;
    m_cacheSize = 1;
    m_debug = false;
  }


  /**
   * converts a BitSet into a list of attribute indexes 
   * @param group the BitSet to convert
   * @return an array of attribute indexes
   **/
  protected int[] attributeList (BitSet group) {
    int count = 0;

    // count how many were selected
    for (int i = 0; i < m_numAttribs; i++) {
      if (group.get(i)) {
	count++;
      }
    }

    int[] list = new int[count];
    count = 0;

    for (int i = 0; i < m_numAttribs; i++) {
      if (group.get(i)) {
	list[count++] = i;
      }
    }

    return  list;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 1.29 $");
  }
}
