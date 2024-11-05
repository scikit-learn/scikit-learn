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
 *    FPGrowth.java
 *    Copyright (C) 2009 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.associations;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Vector;

import weka.core.Attribute;
import weka.core.Capabilities;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.SelectedTag;
import weka.core.SparseInstance;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformationHandler;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;

/**
 <!-- globalinfo-start -->
 * Class implementing the FP-growth algorithm for finding large item sets without candidate generation. Iteratively reduces the minimum support until it finds the required number of rules with the given minimum metric. For more information see:<br/>
 * <br/>
 * J. Han, J.Pei, Y. Yin: Mining frequent patterns without candidate generation. In: Proceedings of the 2000 ACM-SIGMID International Conference on Management of Data, 1-12, 2000.
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;inproceedings{Han2000,
 *    author = {J. Han and J.Pei and Y. Yin},
 *    booktitle = {Proceedings of the 2000 ACM-SIGMID International Conference on Management of Data},
 *    pages = {1-12},
 *    title = {Mining frequent patterns without candidate generation},
 *    year = {2000}
 * }
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -P &lt;attribute index of positive value&gt;
 *  Set the index of the attribute value to consider as 'positive'
 *  for binary attributes in normal dense instances. Index 2 is always
 *  used for sparse instances. (default = 2)</pre>
 * 
 * <pre> -I &lt;max items&gt;
 *  The maximum number of items to include in large items sets (and rules). (default = -1, i.e. no limit.)</pre>
 * 
 * <pre> -N &lt;require number of rules&gt;
 *  The required number of rules. (default = 10)</pre>
 * 
 * <pre> -T &lt;0=confidence | 1=lift | 2=leverage | 3=Conviction&gt;
 *  The metric by which to rank rules. (default = confidence)</pre>
 * 
 * <pre> -C &lt;minimum metric score of a rule&gt;
 *  The minimum metric score of a rule. (default = 0.9)</pre>
 * 
 * <pre> -U &lt;upper bound for minimum support&gt;
 *  Upper bound for minimum support. (default = 1.0)</pre>
 * 
 * <pre> -M &lt;lower bound for minimum support&gt;
 *  The lower bound for the minimum support. (default = 0.1)</pre>
 * 
 * <pre> -D &lt;delta for minimum support&gt;
 *  The delta by which the minimum support is decreased in
 *  each iteration. (default = 0.05)</pre>
 * 
 * <pre> -S
 *  Find all rules that meet the lower bound on
 *  minimum support and the minimum metric constraint.
 *  Turning this mode on will disable the iterative support reduction
 *  procedure to find the specified number of rules.</pre>
 * 
 * <pre> -transactions &lt;comma separated list of attribute names&gt;
 *  Only consider transactions that contain these items (default = no restriction)</pre>
 * 
 * <pre> -rules &lt;comma separated list of attribute names&gt;
 *  Only print rules that contain these items. (default = no restriction)</pre>
 * 
 * <pre> -use-or
 *  Use OR instead of AND for must contain list(s). Use in conjunction
 *  with -transactions and/or -rules</pre>
 * 
 <!-- options-end -->
 *
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com)
 * @version $Revision: 6599 $
 */
public class FPGrowth extends AbstractAssociator 
  implements AssociationRulesProducer, OptionHandler, 
    TechnicalInformationHandler {
  
  /** For serialization */
  private static final long serialVersionUID = 3620717108603442911L;
  
  /**
   * Class for maintaining a frequent item set.
   */
  protected static class FrequentBinaryItemSet 
    implements Serializable, Cloneable {
    
    /** For serialization */
    private static final long serialVersionUID = -6543815873565829448L;

    /** The list of items in the item set */
    protected ArrayList<BinaryItem> m_items = new ArrayList<BinaryItem>();
    
    /** the support of this item set **/
    protected int m_support;
    
    /**
     * Constructor
     * 
     * @param items the items that make up the frequent item set.
     * @param support the support of this item set.
     */
    public FrequentBinaryItemSet(ArrayList<BinaryItem> items, int support) {
      m_items = items;
      m_support = support;
      Collections.sort(m_items);
    }
    
    /**
     * Add an item to this item set.
     * 
     * @param i the item to add.
     */
    public void addItem(BinaryItem i) {
      m_items.add(i);
      Collections.sort(m_items);
    }
    
    /**
     * Set the support for this item set.
     * 
     * @param support the support for this item set.
     */
    public void setSupport(int support) {
      m_support = support;
    }
    
    /**
     * Get the support of this item set.
     * 
     * @return the support of this item set.
     */
    public int getSupport() {
      return m_support;
    }
    
    /**
     * Get the items in this item set.
     * 
     * @return the items in this item set.
     */
    public Collection<BinaryItem> getItems() {
      return m_items;
    }
    
    /**
     * Get a particular item from this item set.
     * 
     * @param index the index of the item to get.
     * @return the item.
     */
    public BinaryItem getItem(int index) {
      return m_items.get(index);
    }
    
    /**
     * Get the number of items in this item set.
     * 
     * @return the number of items in this item set.
     */
    public int numberOfItems() {
      return m_items.size();
    }
    
    /**
     * Get a textual description of this item set.
     * 
     * @return a textual description of this item set.
     */
    public String toString() {
      StringBuffer buff = new StringBuffer();
      Iterator<BinaryItem> i = m_items.iterator();
      
      while (i.hasNext()) {
        buff.append(i.next().toString() + " ");        
      }
      buff.append(": " + m_support);
      return buff.toString();
    }
    
    /**
     * Make a copy of this item set.
     * 
     * @return a copy of this item set.
     */
    public Object clone() {
      ArrayList<BinaryItem> items = new ArrayList<BinaryItem>(m_items);
      return new FrequentBinaryItemSet(items, m_support);
    }
  }
  
  /**
   * Maintains a list of frequent item sets.
   */
  protected static class FrequentItemSets implements Serializable {
    
    /** For serialization */
    private static final long serialVersionUID = 4173606872363973588L;

    /** The list of frequent item sets */
    protected ArrayList<FrequentBinaryItemSet> m_sets = 
      new ArrayList<FrequentBinaryItemSet>();
    
    /** The total number of transactions in the data */
    protected int m_numberOfTransactions;
    
    /**
     * Constructor.
     * 
     * @param numTransactions the total number of transactions in the data.
     */
    public FrequentItemSets(int numTransactions) {
      m_numberOfTransactions = numTransactions;
    }
        
    /**
     * Get an item set.
     * 
     * @param index the index of the item set to get.
     * @return an item set.
     */
    public FrequentBinaryItemSet getItemSet(int index) {
      return m_sets.get(index);
    }
    
    /**
     * Get an iterator that can be used to access all the item sets.
     * 
     * @return an iterator.
     */
    public Iterator<FrequentBinaryItemSet> iterator() {
      return m_sets.iterator();
    }
    
    /**
     * Get the total number of transactions in the data that the item
     * sets were derived from.
     * 
     * @return the total number of transactions in the data.
     */
    public int getNumberOfTransactions() {
      return m_numberOfTransactions;
    }
    
    /**
     * Add an item set.
     * 
     * @param setToAdd the item set to add.
     */
    public void addItemSet(FrequentBinaryItemSet setToAdd) {
      m_sets.add(setToAdd);
    }
    
    /**
     * Sort the item sets according to the supplied comparator.
     * 
     * @param comp the comparator to use.
     */
    public void sort(Comparator<FrequentBinaryItemSet> comp) {
      Collections.sort(m_sets, comp);
    }
    
    /**
     * Get the number of item sets.
     * 
     * @return the number of item sets.
     */
    public int size() {
      return m_sets.size();
    }
    
    /**
     * Sort the item sets. Sorts by item set length. Ties are broken by comparing
     * the items in the two item sets.
     */
    public void sort() {
      Comparator<FrequentBinaryItemSet> compF = new Comparator<FrequentBinaryItemSet>() {
        public int compare(FrequentBinaryItemSet one, FrequentBinaryItemSet two) {
          Collection<BinaryItem> compOne = one.getItems();
          Collection<BinaryItem> compTwo = two.getItems();
          
//          if (one.getSupport() == two.getSupport()) {
            // if supports are equal then list shorter item sets before longer ones
            if (compOne.size() < compTwo.size()) {
              return -1;
            } else if (compOne.size() > compTwo.size()) {
              return 1;
            } else {
              // compare items
              Iterator<BinaryItem> twoIterator = compTwo.iterator();
              for (BinaryItem oneI : compOne) {
                BinaryItem twoI = twoIterator.next();
                int result = oneI.compareTo(twoI);
                if (result != 0) {
                  return result;
                }
              }
              return 0; // equal
            }
            
//            return 0;
    /*      } else if (one.getSupport() > two.getSupport()) {
            // reverse ordering (i.e. descending by support)
            return -1;
          } */
          
    //      return 1;
        }
      };
      
      sort(compF);
    }
    
    /**
     * Get a textual description of this list of item sets.
     * 
     * @param numSets the number of item sets to display.
     * @return a textual description of the item sets.
     */
    public String toString(int numSets) {
      if (m_sets.size() == 0) {
        return "No frequent items sets found!";
      }
      
      StringBuffer result = new StringBuffer();
      result.append("" + m_sets.size() + " frequent item sets found");
      if (numSets > 0) {
        result.append(" , displaying " + numSets);
      }
      result.append(":\n\n");
      
      int count = 0;
      for (FrequentBinaryItemSet i : m_sets) {
        if (numSets > 0 && count > numSets) {
          break;
        }
        result.append(i.toString() + "\n");
        count++;
      }
      
      return result.toString();
    }
  }
  
  /**
   * This class holds the counts for projected tree nodes
   * and header lists.
   */
  protected static class ShadowCounts implements Serializable {
    
    /** For serialization */
    private static final long serialVersionUID = 4435433714185969155L;
    
    /** Holds the counts at different recursion levels */
    private ArrayList<Integer> m_counts = new ArrayList<Integer>();
    
    /**
     * Get the count at the specified recursion depth.
     * 
     * @param recursionLevel the depth of the recursion.
     * @return the count.
     */
    public int getCount(int recursionLevel) {
      if (recursionLevel >= m_counts.size()) {
        return 0;
      } else {
        return m_counts.get(recursionLevel);
      }
    }
    
    /**
     * Increase the count at a given recursion level.
     * 
     * @param recursionLevel the level at which to increase the count.
     * @param incr the amount by which to increase the count.
     */
    public void increaseCount(int recursionLevel, int incr) {
      // basically treat the list like a stack where we
      // can add a new element, or increment the element
      // at the top
      
      if (recursionLevel == m_counts.size()) {
        // new element
        m_counts.add(incr);
      } else if (recursionLevel == m_counts.size() - 1) {
        // otherwise increment the top
        int n = m_counts.get(recursionLevel).intValue();
        m_counts.set(recursionLevel, (n + incr));
      }
    }
    
    /**
     * Remove the count at the given recursion level.
     * 
     * @param recursionLevel the level at which to remove the count.
     */
    public void removeCount(int recursionLevel) {
      if (recursionLevel < m_counts.size()) {
        m_counts.remove(recursionLevel);
      }
    }
  }
  
  /**
   * A node in the FP-tree.
   */
  protected static class FPTreeNode implements Serializable {
                
    /** For serialization */
    private static final long serialVersionUID = 4396315323673737660L;

    /** link to another sibling at this level in the tree */
    protected FPTreeNode m_levelSibling;
    
    /** link to the parent node */
    protected FPTreeNode m_parent;
    
    /** item at this node */
    protected BinaryItem m_item;
    
    /** ID (for graphing the tree) */
    protected int m_ID;
    
    /** the children of this node */
    protected Map<BinaryItem, FPTreeNode> m_children = 
      new HashMap<BinaryItem, FPTreeNode>();
    
    /** counts associated with projected versions of this node */
    protected ShadowCounts m_projectedCounts = new ShadowCounts();
    
    /**
     * Construct a new node with the given parent link and item.
     * 
     * @param parent a pointer to the parent of this node.
     * @param item the item at this node.
     */
    public FPTreeNode(FPTreeNode parent, BinaryItem item) {
      m_parent = parent;
      m_item = item;
    }
    
    /**
     * Insert an item set into the tree at this node. Removes the first item
     * from the supplied item set and makes a recursive call to insert the
     * remaining items.
     * 
     * @param itemSet the item set to insert.
     * @param headerTable the header table for the tree.
     * @param incr the amount by which to increase counts.
     */
    public void addItemSet(Collection<BinaryItem> itemSet, 
        Map<BinaryItem, FPTreeRoot.Header> headerTable, int incr) {
     
      Iterator<BinaryItem> i = itemSet.iterator();
      
      if (i.hasNext()) {
        BinaryItem first = i.next();
        
        FPTreeNode aChild;
        if (!m_children.containsKey(first)) {
          // not in the tree, so add it.
          aChild = new FPTreeNode(this, first);
          m_children.put(first, aChild);
          
          // update the header
          if (!headerTable.containsKey(first)) {
            headerTable.put(first, new FPTreeRoot.Header());
          }
          
          // append new node to header list
          headerTable.get(first).addToList(aChild);
        } else {
          // get the appropriate child node
          aChild = m_children.get(first);
        }
        
        // update counts in header table
        headerTable.get(first).getProjectedCounts().increaseCount(0, incr);
        
        // increase the child's count
        aChild.increaseProjectedCount(0, incr);
        
        // proceed recursively
        itemSet.remove(first);        
        aChild.addItemSet(itemSet, headerTable, incr);
      }
    }
    
    /**
     * Increase the projected count at the given recursion level at this
     * node
     * 
     * @param recursionLevel the recursion level to increase the node count
     * at.
     * @param incr the amount by which to increase the count.
     */
    public void increaseProjectedCount(int recursionLevel, int incr) {
      m_projectedCounts.increaseCount(recursionLevel, incr);
    }
    
    /**
     * Remove the projected count at the given recursion level for this
     * node.
     * 
     * @param recursionLevel the recursion level at which to remove the count.
     */
    public void removeProjectedCount(int recursionLevel) {
      m_projectedCounts.removeCount(recursionLevel);
    }
    
    /**
     * Get the projected count at the given recursion level for this node.
     * 
     * @param recursionLevel the recursion level at which to get the count.
     * @return the count.
     */
    public int getProjectedCount(int recursionLevel) {
      return m_projectedCounts.getCount(recursionLevel);
    }
    
    /**
     * Get the parent node.
     * 
     * @return the parent node.
     */
    public FPTreeNode getParent() {
      return m_parent;
    }
    
    /**
     * Get the item at this node.
     * 
     * @return the item at this node.
     */
    public BinaryItem getItem() {
      return m_item;
    }    
    
    /**
     * Return a textual description of this node for a given recursion
     * level.
     * 
     * @param recursionLevel the recursion depth to use.
     * @return a textual description of this node.
     */
    public String toString(int recursionLevel) {
      return toString("", recursionLevel);
    }

    /**
     * Return a textual description of this node for a given recursion
     * level.
     * 
     * @param prefix a prefix string to prepend.
     * @param recursionLevel the recursion level to use.
     * @return a textual description of this node. 
     */
    public String toString(String prefix, int recursionLevel) {
      StringBuffer buffer = new StringBuffer();
      buffer.append(prefix);
      buffer.append("|  ");
      buffer.append(m_item.toString());
      buffer.append(" (");
      buffer.append(m_projectedCounts.getCount(recursionLevel));
      buffer.append(")\n");
      
      for (FPTreeNode node : m_children.values()) {
        buffer.append(node.toString(prefix + "|  ", recursionLevel));
      }
      return buffer.toString();
    }
    
    protected int assignIDs(int lastID) {
      int currentLastID = lastID + 1;
      m_ID = currentLastID;
      if (m_children != null) {
        Collection<FPTreeNode> kids = m_children.values();
        for (FPTreeNode n : kids) {
          currentLastID = n.assignIDs(currentLastID);
        }
      }
      return currentLastID;
    }
    
    /**
     * Generate a dot graph description string for the tree.
     * 
     * @param text a StringBuffer to store the graph description
     * in.
     */
    public void graphFPTree(StringBuffer text) {
      if (m_children != null) {
        Collection<FPTreeNode> kids = m_children.values();
        for (FPTreeNode n : kids) {
          text.append("N" + n.m_ID);
          text.append(" [label=\"");
          text.append(n.getItem().toString() + " (" + n.getProjectedCount(0) + ")\\n");
          text.append("\"]\n");
          n.graphFPTree(text);
          text.append("N" + m_ID + "->" + "N" + n.m_ID + "\n");
        }
      }
    }
  }
  
  /**
   * Root of the FPTree
   */
  private static class FPTreeRoot extends FPTreeNode {
    
    /** For serialization */
    private static final long serialVersionUID = 632150939785333297L;

    /**
     * Stores a header entry for an FPTree 
     */
    protected static class Header implements Serializable {
      
      /** For serialization */
      private static final long serialVersionUID = -6583156284891368909L;
      
      /** The list of pointers into the tree structure */
      protected List<FPTreeNode> m_headerList = new LinkedList<FPTreeNode>();
      
      /** Projected header counts for this entry */
      protected ShadowCounts m_projectedHeaderCounts = new ShadowCounts();
      
      /**
       * Add a tree node into the list for this header entry.
       * 
       * @param toAdd the node to add.
       */
      public void addToList(FPTreeNode toAdd) {
        m_headerList.add(toAdd);
      }
      
      /**
       * Get the list of nodes for this header entry.
       * 
       * @return the list of nodes for this header entry.
       */
      public List<FPTreeNode> getHeaderList() {
        return m_headerList;
      }
      
      /**
       * Get the projected counts for this header entry.
       * 
       * @return the projected counts for this header entry.
       */
      public ShadowCounts getProjectedCounts() {
        return m_projectedHeaderCounts;
      }
    }
    
    /** Stores the header table as mapped Header entries */
    protected Map<BinaryItem, Header> m_headerTable = 
      new HashMap<BinaryItem, Header>();
    
    /**
     * Create a new FPTreeRoot.
     */
    public FPTreeRoot() {
      super(null, null);
    }
    
    /**
     * Insert an item set into the tree.
     * 
     * @param itemSet the item set to insert into the tree.
     * @param incr the increment by which to increase counters.
     */
    public void addItemSet(Collection<BinaryItem> itemSet, int incr) {
      super.addItemSet(itemSet, m_headerTable, incr);
    }
    
    /**
     * Get the header table for this tree.
     * 
     * @return the header table for this tree.
     */
    public Map<BinaryItem, Header> getHeaderTable() {
      return m_headerTable;
    }
    
    public boolean isEmpty(int recursionLevel) {
      for (FPTreeNode c : m_children.values()) {
        if (c.getProjectedCount(recursionLevel) > 0) {
          return false;
        }
      }
      return true;
    }
    
    /**
     * Get a textual description of the tree at a given recursion
     * (projection) level.
     * 
     * @param pad the string to use as a prefix for indenting nodes.
     * @param recursionLevel the recursion level (projection) to use.
     * @return the textual description of the tree.
     */
    public String toString(String pad, int recursionLevel) {
      StringBuffer result = new StringBuffer();
      result.append(pad);
      result.append("+ ROOT\n");

      for (FPTreeNode node : m_children.values()) {
        result.append(node.toString(pad + "|  ", recursionLevel));
      }
      return result.toString();
    }

    /**
     * Get a textual description of the header table for this tree.
     * 
     * @param recursionLevel the recursion level to use.
     * @return a textual description of the header table for this
     * tree at a given recursion level.
     */
    public String printHeaderTable(int recursionLevel) {
      StringBuffer buffer = new StringBuffer();
      for (BinaryItem item : m_headerTable.keySet()) {
        buffer.append(item.toString());
        buffer.append(" : ");
        buffer.append(m_headerTable.get(item).getProjectedCounts().getCount(recursionLevel));
        buffer.append("\n");
      }
      return buffer.toString();
    }
    
    public void graphHeaderTable(StringBuffer text, int maxID) {

      for (BinaryItem item : m_headerTable.keySet()) {
        Header h = m_headerTable.get(item);
        List<FPTreeNode> headerList = h.getHeaderList();
        if (headerList.size() > 1) {
          text.append("N" + maxID + " [label=\"" + headerList.get(0).getItem().toString() 
              + " (" + h.getProjectedCounts().getCount(0) + ")"
              + "\" shape=plaintext]\n");

          text.append("N" + maxID + "->" + "N" + headerList.get(1).m_ID + "\n");
          for (int i = 1; i < headerList.size() - 1; i++) {
            text.append("N" + headerList.get(i).m_ID + "->" + "N" + headerList.get(i+1).m_ID + "\n");
          }
          maxID++;
        }
      }
    }
  }
  
  private static void nextSubset(boolean[] subset) {
    for (int i = 0; i < subset.length; i++) {
      if (!subset[i]) {
        subset[i] = true;
        break;
      } else {
        subset[i] = false;
      }
    }
  }
  
  private static Collection<Item> getPremise(FrequentBinaryItemSet fis, 
      boolean[] subset) {
    boolean ok = false;
    for (int i = 0; i < subset.length; i++){
      if (!subset[i]) {
        ok = true;
        break;
      }
    }      
    
    if (!ok) {
      return null;
    }
    
    List<Item> premise = new ArrayList<Item>();
    ArrayList<Item> items = new ArrayList<Item>(fis.getItems());

    
    for (int i = 0; i < subset.length; i++) {
      if (subset[i]) {
        premise.add(items.get(i));
      }
    }
    return premise;
  }
  
  private static Collection<Item> getConsequence(FrequentBinaryItemSet fis,
      boolean[] subset) {
    List<Item> consequence = new ArrayList<Item>();
    ArrayList<Item> items = new ArrayList<Item>(fis.getItems());
    
    for (int i = 0; i < subset.length; i++) {
      if (!subset[i]) {
        consequence.add(items.get(i));
      }
    }
    return consequence;
  }
  
  
  /**
   * Generate all association rules, from the supplied frequet item sets,
   * that meet a given minimum metric threshold. Uses a brute force approach.
   * 
   * @param largeItemSets the set of frequent item sets
   * @param metricToUse the metric to use
   * @param metricThreshold the threshold value that a rule must meet
   * @param upperBoundMinSuppAsInstances the upper bound on the support
   * in order to accept the rule
   * @param lowerBoundMinSuppAsInstances the lower bound on the support
   * in order to accept the rule
   * @param totalTransactions the total number of transactions in the data
   * @return a list of association rules
   */
  public static List<AssociationRule> 
    generateRulesBruteForce(FrequentItemSets largeItemSets, 
        DefaultAssociationRule.METRIC_TYPE metricToUse, 
        double metricThreshold, int upperBoundMinSuppAsInstances,
        int lowerBoundMinSuppAsInstances, int totalTransactions) {
    
    List<AssociationRule> rules = new ArrayList<AssociationRule>();
    largeItemSets.sort();
    Map<Collection<BinaryItem>, Integer> frequencyLookup =
      new HashMap<Collection<BinaryItem>, Integer>();
    
    Iterator<FrequentBinaryItemSet> setI = largeItemSets.iterator();
    // process each large item set
    while (setI.hasNext()) {
      FrequentBinaryItemSet fis = setI.next();
      frequencyLookup.put(fis.getItems(), fis.getSupport());
      if (fis.getItems().size() > 1) {
        // generate all the possible subsets for the premise
        boolean[] subset = new boolean[fis.getItems().size()];
        Collection<Item> premise = null;
        Collection<Item> consequence = null;
        while ((premise = getPremise(fis, subset)) != null) {
          if (premise.size() > 0 && premise.size() < fis.getItems().size()) {
            consequence = getConsequence(fis, subset);
            int totalSupport = fis.getSupport();
            int supportPremise = frequencyLookup.get(premise).intValue();
            int supportConsequence = frequencyLookup.get(consequence).intValue();
            
            // a candidate rule
            DefaultAssociationRule candidate = 
              new DefaultAssociationRule(premise, consequence, metricToUse, supportPremise,
                  supportConsequence, totalSupport, totalTransactions);
            if (candidate.getPrimaryMetricValue() > metricThreshold &&
                candidate.getTotalSupport() >= lowerBoundMinSuppAsInstances &&
                candidate.getTotalSupport() <= upperBoundMinSuppAsInstances) {
              // accept this rule
              rules.add(candidate);
            }              
          }
          nextSubset(subset);
        }
      }
    }
    return rules;
  }
  
  public static List<AssociationRule> pruneRules(List<AssociationRule> rulesToPrune,
      ArrayList<Item> itemsToConsider, boolean useOr) {
    ArrayList<AssociationRule> result = new ArrayList<AssociationRule>();
    
    for (AssociationRule r : rulesToPrune) {
      if (r.containsItems(itemsToConsider, useOr)) {
        result.add(r);
      }
    }
    
    return result;
  }
  

  /** The number of rules to find */
  protected int m_numRulesToFind = 10;
  //protected double m_upperBoundMinSupport = 0.36;
  
  /** The upper bound on the minimum support */
  protected double m_upperBoundMinSupport = 1.0;
  
  /** The lower bound on minimum support */
  protected double m_lowerBoundMinSupport = 0.1;
  
  /** The amount by which to decrease the support in each iteration */
  protected double m_delta = 0.05;
  
  /** The number of instances in the data */
  protected int m_numInstances;
  
  /** 
   * When processing data off of disk report progress
   * this frequently (number of instances).
   */
  protected int m_offDiskReportingFrequency = 10000;
  
  /** 
   * If true, just all rules meeting the lower bound on the minimum
   * support will be found. The number of rules to find will be
   * ignored and the iterative reduction of support will not
   * be done. 
   */
  protected boolean m_findAllRulesForSupportLevel = false;
  
  //protected double m_lowerBoundMinSupport = 0.0;
  
  /** The index (1 based) of binary attributes to treat as the positive value */
  protected int m_positiveIndex = 2;
  
  protected DefaultAssociationRule.METRIC_TYPE m_metric = 
    DefaultAssociationRule.METRIC_TYPE.CONFIDENCE;
  
  protected double m_metricThreshold = 0.9;
  
  /** Holds the large item sets found */
  protected FrequentItemSets m_largeItemSets;
  
  /** Holds the rules */
  protected List<AssociationRule> m_rules;
  
  // maximum number of items in a large item set (zero means no limit)
  protected int m_maxItems = -1;
  
  /**
   *  If set, limit the transactions (instances) input to the
   *  algorithm to those that contain these items
   */
  protected String m_transactionsMustContain = "";
  
  /** Use OR rather than AND when considering must contain lists */
  protected boolean m_mustContainOR = false;
  
  /** If set, then only output rules containing these itmes */
  protected String m_rulesMustContain = "";
  
  /**
   * Returns default capabilities of the classifier.
   *
   * @return      the capabilities of this classifier
   */
  public Capabilities getCapabilities() {
    Capabilities result = super.getCapabilities();
    result.disableAll();

    // enable what we can handle
    
    // attributes
    result.enable(Capability.UNARY_ATTRIBUTES);
    result.enable(Capability.BINARY_ATTRIBUTES);
    result.enable(Capability.MISSING_VALUES);

    result.enable(Capability.NO_CLASS);
    
    return result;
  }
  
  /**
   * Returns a string describing this associator
   * 
   * @return a description of the evaluator suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return "Class implementing the FP-growth algorithm for finding" +
    		" large item sets without candidate generation. Iteratively" +
    		" reduces the minimum support until it finds the required" +
    		" number of rules with the given minimum metric." +
    		" For more information see:\n\n" +
    		getTechnicalInformation().toString();
  }
  
  /**
   * Returns an instance of a TechnicalInformation object, containing 
   * detailed information about the technical background of this class,
   * e.g., paper reference or book this class is based on.
   * 
   * @return the technical information about this class
   */
  public TechnicalInformation getTechnicalInformation() {
    TechnicalInformation        result;
    
    result = new TechnicalInformation(Type.INPROCEEDINGS);
    result.setValue(Field.AUTHOR, "J. Han and J.Pei and Y. Yin");
    result.setValue(Field.TITLE, "Mining frequent patterns without candidate generation");
    result.setValue(Field.BOOKTITLE, "Proceedings of the 2000 ACM-SIGMID International" +
    		" Conference on Management of Data");
    result.setValue(Field.YEAR, "2000");
    result.setValue(Field.PAGES, "1-12");
    
    return result;
  }
  
  private boolean passesMustContain(Instance inst, 
      boolean[] transactionsMustContainIndexes, 
      int numInTransactionsMustContainList) {
    
    boolean result = false;
    
    if (inst instanceof SparseInstance) {
      int containsCount = 0;
      for (int i = 0; i < inst.numValues(); i++) {
        int attIndex = inst.index(i);
        
        if (m_mustContainOR) {
          if (transactionsMustContainIndexes[attIndex]) {
            // break here since the operator is OR and this
            // instance contains at least one of the items
            return true;
          }
        } else {
          if (transactionsMustContainIndexes[attIndex]) {
            containsCount++;
          }
        }
      }
      
      if (!m_mustContainOR) {
        if (containsCount == numInTransactionsMustContainList) {
          return true;
        }
      }
    } else {
      int containsCount = 0;
      for (int i = 0; i < transactionsMustContainIndexes.length; i++) {
        if (transactionsMustContainIndexes[i]) {
          if ((int)inst.value(i) == m_positiveIndex - 1) {
            if (m_mustContainOR) {
              // break here since the operator is OR and
              // this instance contains at least one of the
              // requested items
              return true;
            } else {
              containsCount++;
            }
          }
        }
      }
      
      if (!m_mustContainOR) {
        if (containsCount == numInTransactionsMustContainList) {
          return true;
        }
      }
    }
    
    return result;
  }
  
  private void processSingleton(Instance current, 
      ArrayList<BinaryItem> singletons) throws Exception {
    
    if (current instanceof SparseInstance) {
      for (int j = 0; j < current.numValues(); j++) {
        int attIndex = current.index(j);
        singletons.get(attIndex).increaseFrequency();
      }
    } else {
      for (int j = 0; j < current.numAttributes(); j++) {
        if (!current.isMissing(j)) {
          if (current.attribute(j).numValues() == 1 
              || current.value(j) == m_positiveIndex - 1) {
            singletons.get(j).increaseFrequency();
          }
        }
      }
    }
  }
  
  /**
   * Get the singleton items in the data
   * 
   * @param source the source of the data (either Instances or
   * an ArffLoader).
   * @return a list of singleton item sets
   * @throws Exception if the singletons can't be found for some reason
   */
  protected ArrayList<BinaryItem> getSingletons(Object source) 
    throws Exception {
    ArrayList<BinaryItem> singletons = new ArrayList<BinaryItem>();
    Instances data = null;
    
    if (source instanceof Instances) {
      data = (Instances)source;
    } else if (source instanceof weka.core.converters.ArffLoader) {
      data = ((weka.core.converters.ArffLoader)source).getStructure();
    }
    
    for (int i = 0; i < data.numAttributes(); i++) {
      singletons.add(new BinaryItem(data.attribute(i), m_positiveIndex - 1));
    }
    
    if (source instanceof Instances) {
      // set the number of instances
      m_numInstances = data.numInstances();
      
      for (int i = 0; i < data.numInstances(); i++) {
        Instance current = data.instance(i);
        processSingleton(current, singletons);
      }
    } else if (source instanceof weka.core.converters.ArffLoader) {
      weka.core.converters.ArffLoader loader = (weka.core.converters.ArffLoader)source;
      Instance current = null;
      int count = 0;
      while ((current = loader.getNextInstance(data)) != null) {
        processSingleton(current, singletons);
        count++;
        if (count % m_offDiskReportingFrequency == 0) {
          System.err.println("Singletons: done " + count);
        }
      }
      
      // set the number of instances
      m_numInstances = count;
      
      loader.reset();
    }
    
    return singletons;
  }
  
  /**
   * Get the singleton items in the data
   * 
   * @param data the Instances to process
   * @return a list of singleton item sets
   * @throws Exception if the singletons can't be found for some reason
   */
  protected ArrayList<BinaryItem> getSingletons(Instances data) throws Exception {
    return getSingletons((Object)data);
    /*ArrayList<BinaryItem> singletons = new ArrayList<BinaryItem>();
    
    for (int i = 0; i < data.numAttributes(); i++) {
      singletons.add(new BinaryItem(data.attribute(i), m_positiveIndex - 1));
    }
    
    for (int i = 0; i < data.numInstances(); i++) {
      Instance current = data.instance(i);
      if (current instanceof SparseInstance) {
        for (int j = 0; j < current.numValues(); j++) {
          int attIndex = current.index(j);
          singletons.get(attIndex).increaseFrequency();
        }
      } else {
        for (int j = 0; j < data.numAttributes(); j++) {
          if (!current.isMissing(j)) {
            if (current.attribute(j).numValues() == 1 
                || current.value(j) == m_positiveIndex - 1) {
              singletons.get(j).increaseFrequency();
            }
          }
        }
      }
    }
    
    return singletons;*/
  }
  
  /*protected ArrayList<BinaryItem> getFrequent(ArrayList<BinaryItem> items, int minSupport) {
    ArrayList<BinaryItem> frequent = new ArrayList<BinaryItem>();
    for (BinaryItem b : items) {
      if (b.getFrequency() > minSupport) {
        frequent.add(b);
      }
    }
    
    // sort in descending order of support
    Collections.sort(frequent);
    return frequent;
  } */

  /**
   * Inserts a single instance into the FPTree.
   * 
   * @param current the instance to insert
   * @param singletons the singleton item sets
   * @param tree the tree to insert into
   * @param minSupport the minimum support threshold
   */
  private void insertInstance(Instance current, ArrayList<BinaryItem> singletons, 
      FPTreeRoot tree, int minSupport) {
    ArrayList<BinaryItem> transaction = new ArrayList<BinaryItem>();
    if (current instanceof SparseInstance) {
      for (int j = 0; j < current.numValues(); j++) {
        int attIndex = current.index(j);
        if (singletons.get(attIndex).getFrequency() >= minSupport) {
          transaction.add(singletons.get(attIndex));
        }
      }
      Collections.sort(transaction);
      tree.addItemSet(transaction, 1);
    } else {
      for (int j = 0; j < current.numAttributes(); j++) {
        if (!current.isMissing(j)) {
          if (current.attribute(j).numValues() == 1 
              || current.value(j) == m_positiveIndex - 1) {
            if (singletons.get(j).getFrequency() >= minSupport) {
              transaction.add(singletons.get(j));
            }
          }
        }
      }
      Collections.sort(transaction);
      tree.addItemSet(transaction, 1);
    }
  }
  
  /**
   * Construct the frequent pattern tree by inserting each transaction
   * in the data into the tree. Only those items from each transaction that
   * meet the minimum support threshold are inserted.
   * 
   * @param singletons the singleton item sets
   * @param data the Instances containing the transactions
   * @param minSupport the minimum support
   * @return the root of the tree
   */
  protected FPTreeRoot buildFPTree(ArrayList<BinaryItem> singletons,
      Object dataSource, int minSupport) throws Exception {
    
    FPTreeRoot tree = new FPTreeRoot();    
    Instances data = null;
    if (dataSource instanceof Instances) {
      data = (Instances)dataSource;
    } else if (dataSource instanceof weka.core.converters.ArffLoader) {
      data = ((weka.core.converters.ArffLoader)dataSource).getStructure();
    }
    
    if (dataSource instanceof Instances) {
      for (int i = 0; i < data.numInstances(); i++) {
        insertInstance(data.instance(i), singletons, tree, minSupport);
      }
    } else if (dataSource instanceof weka.core.converters.ArffLoader) {
      weka.core.converters.ArffLoader loader = 
        (weka.core.converters.ArffLoader)dataSource;
      Instance current = null;
      int count = 0;
      while ((current = loader.getNextInstance(data)) != null) {
        insertInstance(current, singletons, tree, minSupport);
        count++;
        if (count % m_offDiskReportingFrequency == 0) {
          System.err.println("build tree done: " + count);
        }
      }
    }
    
    return tree;
  }
  
  /**
   * Construct the frequent pattern tree by inserting each transaction
   * in the data into the tree. Only those items from each transaction that
   * meet the minimum support threshold are inserted.
   * 
   * @param singletons the singleton item sets
   * @param data the Instances containing the transactions
   * @param minSupport the minimum support
   * @return the root of the tree
   */
  /*protected FPTreeRoot buildFPTree(ArrayList<BinaryItem> singletons, 
      Instances data, int minSupport) {
    
    FPTreeRoot tree = new FPTreeRoot();
   
    for (int i = 0; i < data.numInstances(); i++) {
      Instance current = data.instance(i);
      ArrayList<BinaryItem> transaction = new ArrayList<BinaryItem>();
      if (current instanceof SparseInstance) {
        for (int j = 0; j < current.numValues(); j++) {
          int attIndex = current.index(j);
          if (singletons.get(attIndex).getFrequency() >= minSupport) {
            transaction.add(singletons.get(attIndex));
          }
        }
        Collections.sort(transaction);
        tree.addItemSet(transaction, 1);
      } else {
        for (int j = 0; j < data.numAttributes(); j++) {
          if (!current.isMissing(j)) {
            if (current.attribute(j).numValues() == 1 
                || current.value(j) == m_positiveIndex - 1) {
              if (singletons.get(j).getFrequency() >= minSupport) {
                transaction.add(singletons.get(j));
              }
            }
          }
        }
        Collections.sort(transaction);
        tree.addItemSet(transaction, 1);
      }
    }
    
    return tree;
  }*/
  
  /**
   * Find large item sets in the FP-tree.
   * 
   * @param tree the root of the tree to mine
   * @param largeItemSets holds the large item sets found
   * @param recursionLevel the recursion level for the current projected
   * counts
   * @param conditionalItems the current set of items that the current
   * (projected) tree is conditional on
   * @param minSupport the minimum acceptable support
   */
  protected void mineTree(FPTreeRoot tree, FrequentItemSets largeItemSets, 
      int recursionLevel, FrequentBinaryItemSet conditionalItems, int minSupport) {
    
    if (!tree.isEmpty(recursionLevel)) {
      if (m_maxItems > 0 && recursionLevel >= m_maxItems) {
        // don't mine any further
        return;
      }
      
      Map<BinaryItem, FPTreeRoot.Header> headerTable = tree.getHeaderTable();
      Set<BinaryItem> keys = headerTable.keySet();
//      System.err.println("Number of freq item sets collected " + largeItemSets.size());
      Iterator<BinaryItem> i = keys.iterator();
      while (i.hasNext()) {
        BinaryItem item = i.next();
        FPTreeRoot.Header itemHeader = headerTable.get(item);
        
        // check for minimum support at this level
        int support = itemHeader.getProjectedCounts().getCount(recursionLevel);
        if (support >= minSupport) {          
          // process header list at this recursion level
          for (FPTreeNode n : itemHeader.getHeaderList()) {
            // push count up path to root
            int currentCount = n.getProjectedCount(recursionLevel);
            if (currentCount > 0) {                            
              FPTreeNode temp = n.getParent();
              while (temp != tree) {
                // set/increase for the node
                temp.increaseProjectedCount(recursionLevel + 1, currentCount);

                // set/increase for the header table
                headerTable.get(temp.getItem()).
                getProjectedCounts().increaseCount(recursionLevel + 1, currentCount);
                
                temp = temp.getParent();
              }
            }
          }
          
          FrequentBinaryItemSet newConditional = 
            (FrequentBinaryItemSet) conditionalItems.clone();
          
          // this item gets added to the conditional items
          newConditional.addItem(item);
          newConditional.setSupport(support);
          
          // now add this conditional item set to the list of large item sets
          largeItemSets.addItemSet(newConditional);
          
          // now recursively process the new tree
          mineTree(tree, largeItemSets, recursionLevel + 1, newConditional,
              minSupport);
          
          // reverse the propagated counts
          for (FPTreeNode n : itemHeader.getHeaderList()) {
            FPTreeNode temp = n.getParent();
            while (temp != tree) {
              temp.removeProjectedCount(recursionLevel + 1);
              temp = temp.getParent();
            }
          }
          
          // reverse the propagated counts in the header list
          // at this recursion level
          for (FPTreeRoot.Header h : headerTable.values()) {
            h.getProjectedCounts().removeCount(recursionLevel + 1);
          }          
        }
      }
    }
  }
  
  /**
   * Construct a new FPGrowth object.
   */
  public FPGrowth() {
    resetOptions();
  }
  
  /**
   * Reset all options to their default values.
   */
  public void resetOptions() {
    m_delta = 0.05;
    m_metricThreshold = 0.9;
    m_numRulesToFind = 10;
    m_lowerBoundMinSupport = 0.1;
    m_upperBoundMinSupport = 1.0;
//    m_minSupport = -1;
    m_positiveIndex = 2;
    m_transactionsMustContain = "";
    m_rulesMustContain = "";
    m_mustContainOR = false;
  }
  
  /**
   * Tip text for this property suitable for displaying
   * in the GUI.
   * 
   * @return the tip text for this property.
   */
  public String positiveIndexTipText() {
    return "Set the index of binary valued attributes that is to be considered" +
    		" the positive index. Has no effect for sparse data (in this case" +
    		" the first index (i.e. non-zero values) is always treated as " +
    		" positive. Also has no effect for unary valued attributes (i.e." +
    		" when using the Weka Apriori-style format for market basket data," +
    		" which uses missing value \"?\" to indicate" +
    		" absence of an item.";
  }
  
  /**
   * Set the index of the attribute value to consider as positive
   * for binary attributes in normal dense instances. Index 1 is always
   * used for sparse instances.
   * 
   * @param index the index to use for positive values in binary attributes.
   */
  public void setPositiveIndex(int index) {
    m_positiveIndex = index;
  }
  
  /**
   * Get the index of the attribute value to consider as positive
   * for binary attributes in normal dense instances. Index 1 is always
   * used for sparse instances.
   * 
   * @return the index to use for positive values in binary attributes.
   */
  public int getPositiveIndex() {
    return m_positiveIndex;
  }
  
  /**
   * Set the desired number of rules to find.
   * 
   * @param numR the number of rules to find.
   */
  public void setNumRulesToFind(int numR) {
    m_numRulesToFind = numR;
  }
  
  /**
   * Get the number of rules to find.
   * 
   * @return the number of rules to find.
   */
  public int getNumRulesToFind() {
    return m_numRulesToFind;
  }
  
  /**
   * Tip text for this property suitable for displaying
   * in the GUI.
   * 
   * @return the tip text for this property.
   */
  public String numRulesToFindTipText() {
    return "The number of rules to output";
  }
  
  /**
   * Set the metric type to use.
   * 
   * @param d the metric type
   */
  public void setMetricType(SelectedTag d) {
    int ordinal =  d.getSelectedTag().getID();
    for (DefaultAssociationRule.METRIC_TYPE m : DefaultAssociationRule.METRIC_TYPE.values()) {
      if (m.ordinal() == ordinal) {
        m_metric = m;
        break;
      }
    }
  }
  
  /**
   * Set the maximum number of items to include in large items sets.
   * 
   * @param max the maxim number of items to include in large item sets.
   */
  public void setMaxNumberOfItems(int max) {
    m_maxItems = max;
  }
  
  /**
   * Gets the maximum number of items to be included in large item sets.
   * 
   * @return the maximum number of items to be included in large items sets.
   */
  public int getMaxNumberOfItems() {
    return m_maxItems;
  }
  
  /**
   * Tip text for this property suitable for displaying
   * in the GUI.
   * 
   * @return the tip text for this property.
   */
  public String maxNumberOfItemsTipText() {
    return "The maximum number of items to include in frequent item sets. -1 " +
    		"means no limit.";
  }
  
  /**
   * Get the metric type to use.
   * 
   * @return the metric type to use.
   */
  public SelectedTag getMetricType() {
    return new SelectedTag(m_metric.ordinal(), DefaultAssociationRule.TAGS_SELECTION);
  }
  
  /**
   * Tip text for this property suitable for displaying
   * in the GUI.
   * 
   * @return the tip text for this property.
   */
  public String metricTypeTipText() {
    return "Set the type of metric by which to rank rules. Confidence is "
    +"the proportion of the examples covered by the premise that are also "
    +"covered by the consequence(Class association rules can only be mined using confidence). Lift is confidence divided by the "
    +"proportion of all examples that are covered by the consequence. This "
    +"is a measure of the importance of the association that is independent "
    +"of support. Leverage is the proportion of additional examples covered "
    +"by both the premise and consequence above those expected if the "
    +"premise and consequence were independent of each other. The total "
    +"number of examples that this represents is presented in brackets "
    +"following the leverage. Conviction is "
    +"another measure of departure from independence.";
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String minMetricTipText() {
    return "Minimum metric score. Consider only rules with scores higher than "
      +"this value.";
  }

  /**
   * Get the value of minConfidence.
   *
   * @return Value of minConfidence.
   */
  public double getMinMetric() {
    
    return m_metricThreshold;
  }
  
  /**
   * Set the value of minConfidence.
   *
   * @param v  Value to assign to minConfidence.
   */
  public void setMinMetric(double v) {
    
    m_metricThreshold = v;
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String transactionsMustContainTipText() {
    return "Limit input to FPGrowth to those transactions (instances)" +
    		" that contain these items. Provide a comma separated" +
    		" list of attribute names.";
  }
  
  /**
   * Set the comma separated list of items that transactions
   * must contain in order to be considered for large
   * item sets and rules.
   * 
   * @param list a comma separated list of items (empty
   * string indicates no restriction on the transactions).
   */
  public void setTransactionsMustContain(String list) {
    m_transactionsMustContain = list;
  }
  
  /**
   * Gets the comma separated list of items that
   * transactions must contain in order to be considered
   * for large item sets and rules.
   * 
   * @return return the comma separated list of
   * items that transactions must contain.
   */
  public String getTransactionsMustContain() {
    return m_transactionsMustContain;
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String rulesMustContainTipText() {
    return "Only print rules that contain these items. Provide " +
    		"a comma separated list of attribute names.";
  }
  
  /**
   * Set the comma separated list of items that rules
   * must contain in order to be output. 
   * 
   * @param list a comma separated list of items (empty
   * string indicates no restriction on the rules).
   */
  public void setRulesMustContain(String list) {
    m_rulesMustContain = list;
  }
  
  /**
   * Get the comma separated list of items that
   * rules must contain in order to be output.
   * 
   * @return the comma separated list of items
   * that rules must contain in order to be output.
   */
  public String getRulesMustContain() {
    return m_rulesMustContain;
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String useORForMustContainListTipText() {
    return "Use OR instead of AND for transactions/rules must contain lists.";
  }
  
  /**
   * Set whether to use OR rather than AND when considering
   * must contain lists.
   * 
   * @param b true if OR should be used instead of AND when
   * considering transaction and rules must contain lists.
   */
  public void setUseORForMustContainList(boolean b) {
    m_mustContainOR = b;
  }
  
  /**
   * Gets whether OR is to be used rather than AND when
   * considering must contain lists.
   * 
   * @return true if OR is used instead of AND.
   */
  public boolean getUseORForMustContainList() {
    return m_mustContainOR;
  }
  
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying, in the explorer/experimenter gui
   */
  public String deltaTipText() {
    return "Iteratively decrease support by this factor. Reduces support "
      +"until min support is reached or required number of rules has been "
      +"generated.";
  }
    
  /**
   * Get the value of delta.
   *
   * @return Value of delta.
   */
  public double getDelta() {
    
    return m_delta;
  }
  
  /**
   * Set the value of delta.
   *
   * @param v  Value to assign to delta.
   */
  public void setDelta(double v) {
    
    m_delta = v;
  }
  
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String lowerBoundMinSupportTipText() {
    return "Lower bound for minimum support as a fraction or number of instances.";
  }

  /**
   * Get the value of lowerBoundMinSupport.
   *
   * @return Value of lowerBoundMinSupport.
   */
  public double getLowerBoundMinSupport() {
    
    return m_lowerBoundMinSupport;
  }
  
  /**
   * Set the value of lowerBoundMinSupport.
   *
   * @param v  Value to assign to lowerBoundMinSupport.
   */
  public void setLowerBoundMinSupport(double v) {
    
    m_lowerBoundMinSupport = v;
  }
  
  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String upperBoundMinSupportTipText() {
    return "Upper bound for minimum support as a fraction or number of instances. " +
      "Start iteratively decreasing " +
      "minimum support from this value.";
  }

  /**
   * Get the value of upperBoundMinSupport.
   *
   * @return Value of upperBoundMinSupport.
   */
  public double getUpperBoundMinSupport() {
    
    return m_upperBoundMinSupport;
  }
  
  /**
   * Set the value of upperBoundMinSupport.
   *
   * @param v  Value to assign to upperBoundMinSupport.
   */
  public void setUpperBoundMinSupport(double v) {
    
    m_upperBoundMinSupport = v;
  }
  
  /**
   * Tip text for this property suitable for displaying
   * in the GUI.
   * 
   * @return the tip text for this property.
   */
  public String findAllRulesForSupportLevelTipText() {
    return "Find all rules that meet " +
    "the lower bound on minimum support and the minimum metric constraint. " +
    "Turning this mode on will disable the iterative support reduction " +
    "procedure to find the specified number of rules.";
  }
  
  /**
   * If true then turn off the iterative support reduction method
   * of finding x rules that meet the minimum support and metric
   * thresholds and just return all the rules that meet the
   * lower bound on minimum support and the minimum metric.
   * 
   * @param s true if all rules meeting the lower bound on the support
   * and minimum metric thresholds are to be found.
   */
  public void setFindAllRulesForSupportLevel(boolean s) {
    m_findAllRulesForSupportLevel = s;
  }
  
  /**
   * Get whether all rules meeting the lower bound on min support
   * and the minimum metric threshold are to be found.
   * 
   * @return true if all rules meeting the lower bound on min
   * support and the min metric threshold are to be found.
   */
  public boolean getFindAllRulesForSupportLevel() {
    return m_findAllRulesForSupportLevel;
  }
  
  /**
   * Set how often to report some progress when the data is
   * being read incrementally off of the disk rather than
   * loaded into memory.
   * 
   * @param freq the frequency to print progress.
   */
  public void setOffDiskReportingFrequency(int freq) {
    m_offDiskReportingFrequency = freq;
  }
  
  /* public void setMinimumSupport(double minSupp) {
    m_minSupport = minSupp;
  }
  
  public double getMinimumSupport() {
    return m_minSupport;
  } */    
  
  /**
   * Gets the list of mined association rules.
   * 
   * @return the list of association rules discovered during mining.
   * Returns null if mining hasn't been performed yet.
   */
  public AssociationRules getAssociationRules() {
    List<AssociationRule> rulesToReturn = new ArrayList<AssociationRule>();
    
    int count = 0;
    for (AssociationRule r : m_rules) {
      rulesToReturn.add(r);
      count++;
      if (!m_findAllRulesForSupportLevel && count == m_numRulesToFind) {
        break;
      }
    }
    
    return new AssociationRules(rulesToReturn, this);
  }
  
  /**
   * Gets a list of the names of the metrics output for
   * each rule. This list should be the same (in terms of
   * the names and order thereof) as that produced by
   * AssociationRule.getMetricNamesForRule().
   * 
   * @return an array of the names of the metrics available
   * for each rule learned by this producer.
   */
  public String[] getRuleMetricNames() {
    String[] metricNames = new String[DefaultAssociationRule.TAGS_SELECTION.length];
    
    for (int i = 0; i < DefaultAssociationRule.TAGS_SELECTION.length; i++) {
      metricNames[i] = DefaultAssociationRule.TAGS_SELECTION[i].getReadable();
    }
    
    return metricNames;
  }
  
  /**
   * Returns true if this AssociationRulesProducer can actually
   * produce rules. Most implementing classes will always return
   * true from this method (obviously :-)). However, an implementing
   * class that actually acts as a wrapper around things that may
   * or may not implement AssociationRulesProducer will want to
   * return false if the thing they wrap can't produce rules.
   * 
   * @return true if this producer can produce rules in its current
   * configuration
   */
  public boolean canProduceRules() {
    return true;
  }
  
  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration<Option> listOptions() {
    Vector<Option> newVector = new Vector<Option>();
    
    String string00 = "\tSet the index of the attribute value to consider as 'positive'\n\t"
   + "for binary attributes in normal dense instances. Index 2 is always\n\t"
   + "used for sparse instances. (default = 2)";    
    String string0 = "\tThe maximum number of items to include " +
    		"in large items sets (and rules). (default " +
    		"= -1, i.e. no limit.)"; 
      
    String string1 = "\tThe required number of rules. (default = " 
      + m_numRulesToFind + ")";
    String string2 = "\tThe minimum metric score of a rule. (default" +
    		" = " + m_metricThreshold + ")";
    String string3 = "\tThe metric by which to rank rules. (default"
      + " = confidence)";
    String string4 = "\tThe lower bound for the minimum support as a fraction" +
    		" or number of instances. (default = "
      + m_lowerBoundMinSupport + ")";
    String string5 = "\tUpper bound for minimum support as a fraction or number of instances. "
      + "(default = 1.0)";
    String string6 = "\tThe delta by which the minimum support is decreased in\n"
     + "\teach iteration as a fraction or number of instances. (default = " + m_delta + ")";
    String string7 = "\tFind all rules that meet the lower bound on\n\t" +
    		"minimum support and the minimum metric constraint.\n\t" +
    		"Turning this mode on will disable the iterative support reduction\n\t" +
    		"procedure to find the specified number of rules.";
    String string8 = "\tOnly consider transactions that contain these items (default = no restriction)";
    String string9 = "\tOnly print rules that contain these items. (default = no restriction)";
    String string10 = "\tUse OR instead of AND for must contain list(s). Use in conjunction" +
    		"\n\twith -transactions and/or -rules";
    
    newVector.add(new Option(string00, "P", 1, "-P <attribute index of positive value>"));
    newVector.add(new Option(string0, "I", 1, "-I <max items>"));
    newVector.add(new Option(string1, "N", 1, "-N <require number of rules>"));
    newVector.add(new Option(string3, "T", 1, "-T <0=confidence | 1=lift | "
                                    + "2=leverage | 3=Conviction>"));
    newVector.add(new Option(string2, "C", 1, "-C <minimum metric score of a rule>"));
    newVector.add(new Option(string5, "U", 1, "-U <upper bound for minimum support>"));
    newVector.add(new Option(string4, "M", 1, "-M <lower bound for minimum support>"));
    newVector.add(new Option(string6, "D", 1, "-D <delta for minimum support>"));
    newVector.add(new Option(string7, "S", 0, "-S"));
    newVector.add(new Option(string8, "transactions", 1, "-transactions <comma separated " +
    		"list of attribute names>"));
    newVector.add(new Option(string9, "rules", 1, "-rules <comma separated list " +
    		"of attribute names>"));
    newVector.add(new Option(string10, "use-or", 0, "-use-or"));
    
    return newVector.elements();
  }
  
  /**
   * 
   * Parses a given list of options. <p/>
   * 
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -P &lt;attribute index of positive value&gt;
   *  Set the index of the attribute value to consider as 'positive'
   *  for binary attributes in normal dense instances. Index 2 is always
   *  used for sparse instances. (default = 2)</pre>
   * 
   * <pre> -I &lt;max items&gt;
   *  The maximum number of items to include in large items sets (and rules). (default = -1, i.e. no limit.)</pre>
   * 
   * <pre> -N &lt;require number of rules&gt;
   *  The required number of rules. (default = 10)</pre>
   * 
   * <pre> -T &lt;0=confidence | 1=lift | 2=leverage | 3=Conviction&gt;
   *  The metric by which to rank rules. (default = confidence)</pre>
   * 
   * <pre> -C &lt;minimum metric score of a rule&gt;
   *  The minimum metric score of a rule. (default = 0.9)</pre>
   * 
   * <pre> -U &lt;upper bound for minimum support&gt;
   *  Upper bound for minimum support. (default = 1.0)</pre>
   * 
   * <pre> -M &lt;lower bound for minimum support&gt;
   *  The lower bound for the minimum support. (default = 0.1)</pre>
   * 
   * <pre> -D &lt;delta for minimum support&gt;
   *  The delta by which the minimum support is decreased in
   *  each iteration. (default = 0.05)</pre>
   * 
   * <pre> -S
   *  Find all rules that meet the lower bound on
   *  minimum support and the minimum metric constraint.
   *  Turning this mode on will disable the iterative support reduction
   *  procedure to find the specified number of rules.</pre>
   * 
   * <pre> -transactions &lt;comma separated list of attribute names&gt;
   *  Only consider transactions that contain these items (default = no restriction)</pre>
   * 
   * <pre> -rules &lt;comma separated list of attribute names&gt;
   *  Only print rules that contain these items. (default = no restriction)</pre>
   * 
   * <pre> -use-or
   *  Use OR instead of AND for must contain list(s). Use in conjunction
   *  with -transactions and/or -rules</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    resetOptions();
    String positiveIndexString = Utils.getOption('P', options);
    String maxItemsString = Utils.getOption('I', options);
    String numRulesString = Utils.getOption('N', options);
    String minMetricString = Utils.getOption('C', options);
    String metricTypeString = Utils.getOption("T", options);
    String lowerBoundSupportString = Utils.getOption("M", options);
    String upperBoundSupportString = Utils.getOption("U", options);
    String deltaString = Utils.getOption("D", options);
    String transactionsString = Utils.getOption("transactions", options);
    String rulesString = Utils.getOption("rules", options);

    if (positiveIndexString.length() != 0) {
      setPositiveIndex(Integer.parseInt(positiveIndexString));
    }
    
    if (maxItemsString.length() != 0) {
      setMaxNumberOfItems(Integer.parseInt(maxItemsString));
    }
    
    if (metricTypeString.length() != 0) {
      setMetricType(new SelectedTag(Integer.parseInt(metricTypeString),
          DefaultAssociationRule.TAGS_SELECTION));
    }
    
    if (numRulesString.length() != 0) {
      setNumRulesToFind(Integer.parseInt(numRulesString));
    }
    
    if (minMetricString.length() != 0) {
      setMinMetric(Double.parseDouble(minMetricString));
    }
    
    if (deltaString.length() != 0) {
      setDelta(Double.parseDouble(deltaString));
    }
    
    if (lowerBoundSupportString.length() != 0) {
      setLowerBoundMinSupport(Double.parseDouble(lowerBoundSupportString));
    }
    
    if (upperBoundSupportString.length() != 0) {
      setUpperBoundMinSupport(Double.parseDouble(upperBoundSupportString));
    }
    
    if (transactionsString.length() != 0) {
      setTransactionsMustContain(transactionsString);
    }
    
    if (rulesString.length() > 0) {
      setRulesMustContain(rulesString);
    }
    
    setUseORForMustContainList(Utils.getFlag("use-or", options));
    
    setFindAllRulesForSupportLevel(Utils.getFlag('S', options));
  }
  
  /**
   * Gets the current settings of the classifier.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    ArrayList<String> options = new ArrayList<String>();
    
    options.add("-P"); options.add("" + getPositiveIndex());
    options.add("-I"); options.add("" + getMaxNumberOfItems());
    options.add("-N"); options.add("" + getNumRulesToFind());
    options.add("-T"); options.add("" + getMetricType().getSelectedTag().getID());
    options.add("-C"); options.add("" + getMinMetric());
    options.add("-D"); options.add("" + getDelta());
    options.add("-U"); options.add("" + getUpperBoundMinSupport());
    options.add("-M"); options.add("" + getLowerBoundMinSupport());
    if (getFindAllRulesForSupportLevel()) {
      options.add("-S");
    }
    
    if (getTransactionsMustContain().length() > 0) {
      options.add("-transactions"); options.add(getTransactionsMustContain());
    }
    
    if (getRulesMustContain().length() > 0) {
      options.add("-rules"); options.add(getRulesMustContain());
    }
    
    if (getUseORForMustContainList()) {
      options.add("-use-or");
    }
    
    return options.toArray(new String[1]);
  }
  
  private Instances parseTransactionsMustContain(Instances data) {
    String[] split = m_transactionsMustContain.trim().split(",");
    boolean[] transactionsMustContainIndexes = new boolean[data.numAttributes()];
    int numInTransactionsMustContainList = split.length;
    
    for (int i = 0; i < split.length; i++) {
      String attName = split[i].trim();
      Attribute att = data.attribute(attName);
      if (att == null) {
        System.err.println("[FPGrowth] : WARNING - can't find attribute " 
            + attName + " in the data.");
        numInTransactionsMustContainList--;
      } else {
        transactionsMustContainIndexes[att.index()] = true;
      }
    }
    
    if (numInTransactionsMustContainList == 0) {
      return data;
    } else {
      Instances newInsts = new Instances(data, 0);
      for (int i = 0; i < data.numInstances(); i++) {
        if (passesMustContain(data.instance(i), 
            transactionsMustContainIndexes, numInTransactionsMustContainList)) {
          newInsts.add(data.instance(i));
        }
      }
      newInsts.compactify();
      return newInsts;
    }
  }
  
  private ArrayList<Item> parseRulesMustContain(Instances data) {
    ArrayList<Item> result = new ArrayList<Item>();
    
    String[] split = m_rulesMustContain.trim().split(",");
    
    for (int i = 0; i < split.length; i++) {
      String attName = split[i].trim();
      Attribute att = data.attribute(attName);
      if (att == null) {
        System.err.println("[FPGrowth] : WARNING - can't find attribute " 
            + attName + " in the data.");
      } else {
        BinaryItem tempI = null;
        try {
          tempI = new BinaryItem(att, m_positiveIndex - 1);
        } catch (Exception e) {
          // this should never happen
          e.printStackTrace();
        }
        result.add(tempI);
      }
    }
    
    return result;
  }
  
  /**
   * Method that generates all large item sets with a minimum support, and from
   * these all association rules with a minimum metric (i.e. confidence, 
   * lift etc.).
   *
   * @param source the source of the data. May be an Instances object or
   * an ArffLoader. In the case of the latter, the two passes over the 
   * data that FPGrowth requires will be done off of disk (i.e. only one
   * instance will be in memory at any one time).
   * @throws Exception if rules can't be built successfully
   */
  private void buildAssociations(Object source) throws Exception {
    Instances data = null;
    Capabilities capabilities = getCapabilities();
    boolean arffLoader = false;
    
    if (source instanceof weka.core.converters.ArffLoader) {
      data = ((weka.core.converters.ArffLoader)source).getStructure();
      capabilities.setMinimumNumberInstances(0);
      arffLoader = true;
    } else {
      data = (Instances)source;
    }

    // can we handle the data?
    capabilities.testWithFail(data);
    
    // prune any instances that don't contain the requested items (if any)
    // can only do this if we are not reading the data incrementally
    if (m_transactionsMustContain.length() > 0 && (source instanceof Instances)) {
      data = parseTransactionsMustContain(data);
      getCapabilities().testWithFail(data);
    }
    
    ArrayList<Item> rulesMustContain = null;
    if (m_rulesMustContain.length() > 0) {
      rulesMustContain = parseRulesMustContain(data);
    }
    
    ArrayList<BinaryItem> singletons = getSingletons(source);
    
    int upperBoundMinSuppAsInstances = (m_upperBoundMinSupport > 1) 
      ? (int) m_upperBoundMinSupport
      : (int)Math.ceil(m_upperBoundMinSupport * m_numInstances);
      
    int lowerBoundMinSuppAsInstances = (m_lowerBoundMinSupport > 1)
      ? (int)m_lowerBoundMinSupport
      : (int)Math.ceil(m_lowerBoundMinSupport * m_numInstances);
      
    double upperBoundMinSuppAsFraction = (m_upperBoundMinSupport > 1)
      ? m_upperBoundMinSupport / m_numInstances
      : m_upperBoundMinSupport;
      
    double lowerBoundMinSuppAsFraction = (m_lowerBoundMinSupport > 1)
      ? m_lowerBoundMinSupport / m_numInstances
      : m_lowerBoundMinSupport;
      
    double deltaAsFraction = (m_delta > 1)
      ? m_delta / m_numInstances
      : m_delta;
      
    double currentSupport = upperBoundMinSuppAsFraction;      
             
    if (m_findAllRulesForSupportLevel) {
      currentSupport = lowerBoundMinSuppAsFraction;
    }
        
    do {
      if (arffLoader) {
        ((weka.core.converters.ArffLoader)source).reset();
      }
      
      int currentSupportAsInstances = (currentSupport > 1)
      ? (int)currentSupport
      : (int)Math.ceil(currentSupport * m_numInstances);
      
      // build the FPTree
      if (arffLoader) {
        System.err.println("Building FP-tree...");
      }
      FPTreeRoot tree = buildFPTree(singletons, source, currentSupportAsInstances);
      
      FrequentItemSets largeItemSets = new FrequentItemSets(m_numInstances);
      
      if (arffLoader) {
        System.err.println("Mining tree for min supp " + currentSupport);
      }

      // mine the tree
      FrequentBinaryItemSet conditionalItems = 
        new FrequentBinaryItemSet(new ArrayList<BinaryItem>(), 0);
      mineTree(tree, largeItemSets, 0, conditionalItems, currentSupportAsInstances);      

      m_largeItemSets = largeItemSets;
      
      if (arffLoader) {
        System.err.println("Number of large item sets: " + m_largeItemSets.size());
      }
      
      // save memory
      tree = null;

      m_rules = 
        generateRulesBruteForce(m_largeItemSets, m_metric, 
            m_metricThreshold, upperBoundMinSuppAsInstances, 
            lowerBoundMinSuppAsInstances, m_numInstances);
      
      if (arffLoader) {
        System.err.println("Number of rules found " + m_rules.size());
      }
      
      if (rulesMustContain != null && rulesMustContain.size() > 0) {
        m_rules = pruneRules(m_rules, rulesMustContain, 
            m_mustContainOR);
      }
      
      if (!m_findAllRulesForSupportLevel) {
        currentSupport -= deltaAsFraction;
        if (currentSupport < lowerBoundMinSuppAsFraction) {
          if (currentSupport + deltaAsFraction > lowerBoundMinSuppAsFraction) {
            // ensure that the lower bound does get evaluated
            currentSupport = lowerBoundMinSuppAsFraction;
          } else {
            break;
          }
        }
      } else {
        // just break out of the loop as we are just finding all rules
        // with a minimum support + metric
        break;
      }      
    } while (m_rules.size() < m_numRulesToFind);
    
    Collections.sort(m_rules);
  }

  /**
   * Method that generates all large item sets with a minimum support, and from
   * these all association rules with a minimum metric (i.e. confidence, 
   * lift etc.).
   *
   * @param data the instances to be used for generating the associations
   * @throws Exception if rules can't be built successfully
   */
  public void buildAssociations(Instances data) throws Exception {
    
    buildAssociations((Object)data);
    return;
  }
      
  /**
   * Output the association rules.
   * 
   * @return a string representation of the model.
   */
  public String toString() {
//    return m_largeItemSets.toString(m_numItemSetsToFind);
    if (m_rules == null) {
      return "FPGrowth hasn't been trained yet!";
    }

    StringBuffer result = new StringBuffer();
    int numRules = (m_rules.size() < m_numRulesToFind)
      ? m_rules.size()
      : m_numRulesToFind;
      
    if (m_rules.size() == 0) {
      return "No rules found!";
    } else {      
      result.append("FPGrowth found " + m_rules.size() + " rules");
      if (!m_findAllRulesForSupportLevel) {
        result.append(" (displaying top " + numRules + ")");
      }
      
      if (m_transactionsMustContain.length() > 0 || 
          m_rulesMustContain.length() > 0) {
        result.append("\n");
        if (m_transactionsMustContain.length() > 0) {
          result.append("\nUsing only transactions that contain: " 
              + m_transactionsMustContain);
        }
        if (m_rulesMustContain.length() > 0) {
          result.append("\nShowing only rules that contain: " 
              + m_rulesMustContain);
        }
      }
      
      result.append("\n\n");
    }

    int count = 0;
    for (AssociationRule r : m_rules) {
      result.append(Utils.doubleToString((double)count+1,
          (int)(Math.log(numRules)/Math.log(10)+1), 0) + ". ");
      result.append(r + "\n");
      count++;
      if (!m_findAllRulesForSupportLevel && count == m_numRulesToFind) {
        break;
      }
    }
    return result.toString();
  }
  
  /**
   * Assemble a dot graph representation of the FP-tree.
   * 
   * @param tree the root of the FP-tree
   * @return a graph representation as a String in dot format.
   */
  public String graph(FPTreeRoot tree) {
    //int maxID = tree.assignIDs(-1);
    
    
    StringBuffer text = new StringBuffer();
    text.append("digraph FPTree {\n");
    text.append("N0 [label=\"ROOT\"]\n");
    tree.graphFPTree(text);
    
//    tree.graphHeaderTable(text, maxID+1);
    text.append("}\n");
    
    return text.toString();
  }
    
  /**
   * Returns the revision string.
   * 
   * @return            the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6599 $");
  }
  
  /**
   * Main method.
   * 
   * @param args the commandline options
   */
  public static void main(String[] args) {
    try {
      String[] argsCopy = args.clone();
      if (Utils.getFlag('h', argsCopy) || Utils.getFlag("help", argsCopy)) {
        runAssociator(new FPGrowth(), args);
        System.out.println("-disk\n\tProcess data off of disk instead of loading\n\t" +
        		"into main memory. This is a command line only option.");
        return;
      }
        
      if (!Utils.getFlag("disk", args)) {
        runAssociator(new FPGrowth(), args);
      } else {
        String filename;
        filename = Utils.getOption('t', args);
        weka.core.converters.ArffLoader loader = null;
        if (filename.length() != 0) {
          loader = new weka.core.converters.ArffLoader();
          loader.setFile(new java.io.File(filename));
        } else {
          throw new Exception("No training file specified!");
        }
        FPGrowth fpGrowth = new FPGrowth();
        fpGrowth.setOptions(args);
        Utils.checkForRemainingOptions(args);
        fpGrowth.buildAssociations(loader);
        System.out.print(fpGrowth.toString());
      }
    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }
}

