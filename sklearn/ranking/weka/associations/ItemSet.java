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
 *    ItemSet.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.associations;

import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;

import java.io.Serializable;
import java.util.Enumeration;
import java.util.Hashtable;

/**
 * Class for storing a set of items. Item sets are stored in a lexicographic
 * order, which is determined by the header information of the set of instances
 * used for generating the set of items. All methods in this class assume that
 * item sets are stored in lexicographic order.
 * The class provides the general methods used for item sets in class - and 
 * standard association rule mining.
 *
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version $Revision: 6514 $
 */
public class ItemSet
  implements Serializable, RevisionHandler {

  /** for serialization */
  private static final long serialVersionUID = 2724000045282835791L;
  
  /** The items stored as an array of of ints. */
  protected int[] m_items;

  /** Counter for how many transactions contain this item set. */
  protected int m_counter;

  /** The total number of transactions */
  protected int m_totalTransactions;
  
  /** 
   * Treat zeros as missing (rather than a value in their
   * own right)
   */
  protected boolean m_treatZeroAsMissing = false;

  /**
   * Constructor
   * @param totalTrans the total number of transactions in the data
   */
  public ItemSet(int totalTrans) {
    m_totalTransactions = totalTrans;
  }

   /**
   * Constructor
   * @param totalTrans the total number of transactions in the data
   * @param array the attribute values encoded in an int array
   */
  public ItemSet(int totalTrans, int[] array){
       
       m_totalTransactions = totalTrans;
       m_items = array;
       m_counter =1;
   }
  
  /** Contsructor
   * @param array the item set represented as an int array
   */  
  public ItemSet(int[] array){
       
       m_items = array;
       m_counter = 0;
   }

  /**
   * Checks if an instance contains an item set.
   *
   * @param instance the instance to be tested
   * @return true if the given instance contains this item set
   */
  public boolean containedBy(Instance instance) {

    if (instance instanceof weka.core.SparseInstance && m_treatZeroAsMissing) {
      int numInstVals = instance.numValues();
      int numItemSetVals = m_items.length;

      for (int p1 = 0, p2 = 0; p1 < numInstVals || p2 < numItemSetVals; ) {
        int instIndex = Integer.MAX_VALUE;
        if (p1 < numInstVals) {
          instIndex = instance.index(p1);
        }
        int itemIndex = p2;

        if (m_items[itemIndex] > -1) {
          if (itemIndex != instIndex) {
            return false;
          } else {
            if (instance.isMissingSparse(p1)) {
              return false;
            }
            if (m_items[itemIndex] != (int)instance.valueSparse(p1)) {
              return false;
            }
          }

          p1++;
          p2++;
        } else {
          if (itemIndex < instIndex) {
            p2++;
          } else if (itemIndex == instIndex){
            p2++;
            p1++;
          }
        }      
      }
    } else {
      for (int i = 0; i < instance.numAttributes(); i++) 
        if (m_items[i] > -1) {
          if (instance.isMissing(i) || 
              (m_treatZeroAsMissing && (int)instance.value(i) == 0))
            return false;
          if (m_items[i] != (int)instance.value(i))
            return false;
        }
    }

    return true;
  }

  /** Deletes all item sets that don't have minimum support.
   * @return the reduced set of item sets
   * @param maxSupport the maximum support
   * @param itemSets the set of item sets to be pruned
   * @param minSupport the minimum number of transactions to be covered
   */
  public static FastVector deleteItemSets(FastVector itemSets, 
					  int minSupport,
					  int maxSupport) {

    FastVector newVector = new FastVector(itemSets.size());

    for (int i = 0; i < itemSets.size(); i++) {
      ItemSet current = (ItemSet)itemSets.elementAt(i);
      if ((current.m_counter >= minSupport) 
	  && (current.m_counter <= maxSupport))
	newVector.addElement(current);
    }
    return newVector;
  }

  /**
   * Tests if two item sets are equal.
   *
   * @param itemSet another item set
   * @return true if this item set contains the same items as the given one
   */
  public boolean equals(Object itemSet) {

    if ((itemSet == null) || !(itemSet.getClass().equals(this.getClass()))) {
      return false;
    }
    if (m_items.length != ((ItemSet)itemSet).m_items.length)
      return false;
    for (int i = 0; i < m_items.length; i++)
      if (m_items[i] != ((ItemSet)itemSet).m_items[i])
	return false;
    return true;
  }

  /**
   * Return a hashtable filled with the given item sets.
   *
   * @param itemSets the set of item sets to be used for filling the hash table
   * @param initialSize the initial size of the hashtable
   * @return the generated hashtable
   */
  public static Hashtable getHashtable(FastVector itemSets, int initialSize) {

    Hashtable hashtable = new Hashtable(initialSize);

    for (int i = 0; i < itemSets.size(); i++) {
      ItemSet current = (ItemSet)itemSets.elementAt(i);
      hashtable.put(current, new Integer(current.m_counter));
    }
    return hashtable;
  }

  /**
   * Produces a hash code for a item set.
   *
   * @return a hash code for a set of items
   */
  public int hashCode() {

    long result = 0;

    for (int i = m_items.length-1; i >= 0; i--)
      result += (i * m_items[i]);
    return (int)result;
  }

  /** Merges all item sets in the set of (k-1)-item sets
   * to create the (k)-item sets and updates the counters.
   * @return the generated (k)-item sets
   * @param totalTrans thetotal number of transactions
   * @param itemSets the set of (k-1)-item sets
   * @param size the value of (k-1)
   */
  public static FastVector mergeAllItemSets(FastVector itemSets, int size, 
					    int totalTrans) {

    FastVector newVector = new FastVector();
    ItemSet result;
    int numFound, k;

    for (int i = 0; i < itemSets.size(); i++) {
      ItemSet first = (ItemSet)itemSets.elementAt(i);
    out:
      for (int j = i+1; j < itemSets.size(); j++) {
	ItemSet second = (ItemSet)itemSets.elementAt(j);
	result = new ItemSet(totalTrans);
	result.m_items = new int[first.m_items.length];
        
	// Find and copy common prefix of size 'size'
	numFound = 0;
	k = 0;
	while (numFound < size) {
	  if (first.m_items[k] == second.m_items[k]) {
	    if (first.m_items[k] != -1) 
	      numFound++;
	    result.m_items[k] = first.m_items[k];
	  } else 
	    break out;
	  k++;
	}
	
	// Check difference
	while (k < first.m_items.length) {
	  if ((first.m_items[k] != -1) && (second.m_items[k] != -1))
	    break;
	  else {
	    if (first.m_items[k] != -1)
	      result.m_items[k] = first.m_items[k];
	    else
	      result.m_items[k] = second.m_items[k];
	  }
	  k++;
	}
	if (k == first.m_items.length) {
	  result.m_counter = 0;
         
	  newVector.addElement(result);
	}
      }
    }
    return newVector;
  }

  /**
   * Prunes a set of (k)-item sets using the given (k-1)-item sets.
   *
   * @param toPrune the set of (k)-item sets to be pruned
   * @param kMinusOne the (k-1)-item sets to be used for pruning
   * @return the pruned set of item sets
   */
  public static FastVector pruneItemSets(FastVector toPrune, Hashtable kMinusOne) {

    FastVector newVector = new FastVector(toPrune.size());
    int help, j;

    for (int i = 0; i < toPrune.size(); i++) {
      ItemSet current = (ItemSet)toPrune.elementAt(i);
      for (j = 0; j < current.m_items.length; j++)
	if (current.m_items[j] != -1) {
	  help = current.m_items[j];
	  current.m_items[j] = -1;
	  if (kMinusOne.get(current) == null) {
	    current.m_items[j] = help;
	    break;
	  } else{ 
	    current.m_items[j] = help;
          }
	}
      if (j == current.m_items.length) 
	newVector.addElement(current);
    }
    return newVector;
  }

  /**
   * Prunes a set of rules.
   *
   * @param rules a two-dimensional array of lists of item sets. The first list
   * of item sets contains the premises, the second one the consequences.
   * @param minConfidence the minimum confidence the rules have to have
   */
  public static void pruneRules(FastVector[] rules, double minConfidence) {

    FastVector newPremises = new FastVector(rules[0].size()),
      newConsequences = new FastVector(rules[1].size()),
      newConf = new FastVector(rules[2].size());
    
    FastVector newLift = null, newLev = null, newConv = null;
    if (rules.length > 3) {
      newLift = new FastVector(rules[3].size());
      newLev = new FastVector(rules[4].size());
      newConv = new FastVector(rules[5].size());
    }

    for (int i = 0; i < rules[0].size(); i++) 
      if (!(((Double)rules[2].elementAt(i)).doubleValue() <
	    minConfidence)) {
	newPremises.addElement(rules[0].elementAt(i));
	newConsequences.addElement(rules[1].elementAt(i));
	newConf.addElement(rules[2].elementAt(i));
	
	if (rules.length > 3) {
	  newLift.addElement(rules[3].elementAt(i));
	  newLev.addElement(rules[4].elementAt(i));
	  newConv.addElement(rules[5].elementAt(i));
	}
      }
    rules[0] = newPremises;
    rules[1] = newConsequences;
    rules[2] = newConf;
    
    if (rules.length > 3) {
      rules[3] = newLift;
      rules[4] = newLev;
      rules[5] = newConv;
    }
  }

  /**
   * Converts the header info of the given set of instances into a set 
   * of item sets (singletons). The ordering of values in the header file 
   * determines the lexicographic order.
   *
   * @param instances the set of instances whose header info is to be used
   * @return a set of item sets, each containing a single item
   * @exception Exception if singletons can't be generated successfully
   */
  public static FastVector singletons(Instances instances) throws Exception {

    FastVector setOfItemSets = new FastVector();
    ItemSet current;

    for (int i = 0; i < instances.numAttributes(); i++) {
      if (instances.attribute(i).isNumeric())
	throw new Exception("Can't handle numeric attributes!");
      for (int j = 0; j < instances.attribute(i).numValues(); j++) {
	current = new ItemSet(instances.numInstances());
	current.m_items = new int[instances.numAttributes()];
	for (int k = 0; k < instances.numAttributes(); k++)
	  current.m_items[k] = -1;
	current.m_items[i] = j;
        
	setOfItemSets.addElement(current);
      }
    }
    return setOfItemSets;
  }

  /**
   * Outputs the support for an item set.
   *
   * @return the support
   */
  public int support() {

    return m_counter;
  }

  /**
   * Returns the contents of an item set as a string.
   *
   * @param instances contains the relevant header information
   * @return string describing the item set
   */
  public String toString(Instances instances) {

    StringBuffer text = new StringBuffer();

    for (int i = 0; i < instances.numAttributes(); i++)
      if (m_items[i] != -1) {
	text.append(instances.attribute(i).name()+'=');
	text.append(instances.attribute(i).value(m_items[i])+' ');
      }
    text.append(m_counter);
    return text.toString();
  }

  /**
   * Updates counter of item set with respect to given transaction.
   *
   * @param instance the instance to be used for ubdating the counter
   */
  public void upDateCounter(Instance instance) {

    if (containedBy(instance))
      m_counter++;
  }

  /**
   * Updates counters for a set of item sets and a set of instances.
   *
   * @param itemSets the set of item sets which are to be updated
   * @param instances the instances to be used for updating the counters
   */
  public static void upDateCounters(FastVector itemSets, Instances instances) {

    for (int i = 0; i < instances.numInstances(); i++) {
      Enumeration enu = itemSets.elements();
      while (enu.hasMoreElements()) 
	((ItemSet)enu.nextElement()).upDateCounter(instances.instance(i));
    }
  }

  /** Gets the counter
   * @return the counter
   */  
  public int counter(){
      
      return m_counter;
  }
  
  /** Gest the item set as an int array
   * @return int array represneting an item set
   */  
  public int[] items(){
      
      return m_items;
  }
  
  /** Gest the index of the value of the specified attribute
   * @param k the attribute index
   * @return the index of the attribute value
   */  
  public int itemAt(int k){
      
      return m_items[k];
  }
  
  /** Sets the counter
   * @param count the counter
   */  
  public void setCounter(int count){
      
      m_counter = count;
  }
  
  /** Sets an item sets
   * @param items an int array representing an item set
   */  
  public void setItem(int[] items){
      
      m_items = items;
  }
  
  /** Sets the index of an attribute value
   * @param value the inex of the attribute value
   * @param k the index of the attribute
   */  
  public void setItemAt(int value, int k){
      
      m_items[k] = value;
  }
  
  /**
   * Sets whether zeros (i.e. the first value of a nominal attribute)
   * should be treated as missing values.
   * 
   * @param z true if zeros should be treated as missing values.
   */
  public void setTreatZeroAsMissing(boolean z) {
    m_treatZeroAsMissing = z;
  }
  
  /**
   * Gets whether zeros (i.e. the first value of a nominal attribute)
   * is to be treated int he same way as missing values.
   * 
   * @return true if zeros are to be treated like missing values.
   */
  public boolean getTreatZeroAsMissing() {
    return m_treatZeroAsMissing;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6514 $");
  }
}
