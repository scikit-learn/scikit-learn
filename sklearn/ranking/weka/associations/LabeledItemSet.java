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
 * LabeledItemSet.java
 * Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
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
 * Class for storing a set of items together with a  class label. Item sets are stored in a lexicographic
 * order, which is determined by the header information of the set of instances
 * used for generating the set of items. All methods in this class assume that
 * item sets are stored in lexicographic order.
 * The class provides the methods used for item sets in class association rule mining.
 * Because every item set knows its class label the training set can be splitted up virtually.
 *
 * @author Stefan Mutter (mutter@cs.waikato.ac.nz)
 * @version $Revision: 1.5 $
 */

public class LabeledItemSet
  extends ItemSet
  implements Serializable, RevisionHandler {

  /** for serialization */
  private static final long serialVersionUID = 4158771925518299903L;

    /** The class label. */
    protected int m_classLabel;

     /** The support of the rule. */
    protected int m_ruleSupCounter;    
  
    /**
     * Constructor
     *
     * @param totalTrans the total number of transactions
     * @param classLabel the class lebel
     */    
    public LabeledItemSet(int totalTrans, int classLabel){
	
	super(totalTrans);
	m_classLabel = classLabel;
    }

    

    /**
     * Deletes all item sets that don't have minimum support and have more than maximum support
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
      LabeledItemSet current = (LabeledItemSet)itemSets.elementAt(i);
      if ((current.m_ruleSupCounter >= minSupport) 
	  && (current.m_ruleSupCounter <= maxSupport))
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
  public final boolean equals(Object itemSet) {

    if (!(this.equalCondset(itemSet)))
      return false;
    if(m_classLabel != ((LabeledItemSet)itemSet).m_classLabel)
	return false;
    
    return true;
  }

  /**
   * Compares two item sets
   * @param itemSet an item set
   * @return true if the the item sets are equal, false otherwise
   */  
    public final boolean equalCondset(Object itemSet) {

    if ((itemSet == null) || !(itemSet.getClass().equals(this.getClass()))) {
      return false;
    }
    if (m_items.length != ((ItemSet)itemSet).items().length)
      return false;
    for (int i = 0; i < m_items.length; i++)
      if (m_items[i] != ((ItemSet)itemSet).itemAt(i))
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
      LabeledItemSet current = (LabeledItemSet)itemSets.elementAt(i);
      hashtable.put(current, new Integer(current.m_classLabel));
    }
   
    return hashtable;
  }



    /**
     * Merges all item sets in the set of (k-1)-item sets
     * to create the (k)-item sets and updates the counters.
     * @return the generated (k)-item sets
     * @param totalTrans the total number of transactions
     * @param itemSets the set of (k-1)-item sets
     * @param size the value of (k-1)
     */
  public static FastVector mergeAllItemSets(FastVector itemSets, int size, 
					    int totalTrans) {
      
      FastVector newVector = new FastVector();
      LabeledItemSet result;
      int numFound, k;
      
      for (int i = 0; i < itemSets.size(); i++) {
	  LabeledItemSet first = (LabeledItemSet)itemSets.elementAt(i);
      out:
	  for (int j = i+1; j < itemSets.size(); j++) {
	      LabeledItemSet second = (LabeledItemSet)itemSets.elementAt(j);
	      while(first.m_classLabel != second.m_classLabel){
		  j++;
		  if(j == itemSets.size())
		      break out;
		  second = (LabeledItemSet)itemSets.elementAt(j);
	      }
	      result = new LabeledItemSet(totalTrans,first.m_classLabel);
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
		  result.m_ruleSupCounter = 0;
		  result.m_counter = 0;
		  newVector.addElement(result);
	      }
	  }
      }
      
    return newVector;
  }
    
  /**
   * Splits the class attribute away. Depending on the invert flag, the instances without class attribute or only the class attribute of all instances is returned
   * @param instances the instances
   * @param invert flag; if true only the class attribute remains, otherweise the class attribute is the only attribute that is deleted.
   * @throws Exception exception if instances cannot be splitted
   * @return Instances without the class attribute or instances with only the class attribute
   */  
    public static Instances divide(Instances instances, boolean invert) throws Exception{
    
	Instances newInstances = new Instances(instances);
	if(instances.classIndex() < 0)
            throw new Exception("For class association rule mining a class attribute has to be specified.");
	if(invert){
		for(int i=0;i<newInstances.numAttributes();i++){
		    if(i!=newInstances.classIndex()){
			newInstances.deleteAttributeAt(i);
			i--;
		    }
		}
	    return newInstances;
	}
	else{
	    newInstances.setClassIndex(-1);
            newInstances.deleteAttributeAt(instances.classIndex());
	    return newInstances;
	}
    }

   
    
     /**
      * Converts the header info of the given set of instances into a set
      * of item sets (singletons). The ordering of values in the header file
      * determines the lexicographic order. Each item set knows its class label.
      * @return a set of item sets, each containing a single item
      * @param instancesNoClass instances without the class attribute
      * @param classes the values of the class attribute sorted according to instances
      * @exception Exception if singletons can't be generated successfully
      */
    public static FastVector singletons(Instances instancesNoClass, Instances classes) throws Exception {

    FastVector cSet, setOfItemSets = new FastVector();
    LabeledItemSet current;
    
    
    //make singletons
    for (int i = 0; i < instancesNoClass.numAttributes(); i++) {
      if (instancesNoClass.attribute(i).isNumeric())
	  throw new Exception("Can't handle numeric attributes!");
      for (int j = 0; j < instancesNoClass.attribute(i).numValues(); j++){
	  for(int k =0; k < (classes.attribute(0)).numValues(); k++){
	      current = new LabeledItemSet(instancesNoClass.numInstances(),k);
	      current.m_items = new int[instancesNoClass.numAttributes()];
	      for (int l = 0; l < instancesNoClass.numAttributes(); l++)
		  current.m_items[l] = -1;
	      current.m_items[i] = j;
	      setOfItemSets.addElement(current);
	  }      
      }
    }
    return setOfItemSets;
  }
    
 
    
  
    /**
   * Prunes a set of (k)-item sets using the given (k-1)-item sets.
   *
   * @param toPrune the set of (k)-item sets to be pruned
   * @param kMinusOne the (k-1)-item sets to be used for pruning
   * @return the pruned set of item sets
   */
  public static FastVector pruneItemSets(FastVector toPrune, Hashtable kMinusOne){

    FastVector newVector = new FastVector(toPrune.size());
    int help, j;

   
    for (int i = 0; i < toPrune.size(); i++) {
      LabeledItemSet current = (LabeledItemSet)toPrune.elementAt(i);
      
      for (j = 0; j < current.m_items.length; j++){
	      if (current.m_items[j] != -1) {
		  help = current.m_items[j];
		  current.m_items[j] = -1;
		  if(kMinusOne.get(current) != null && (current.m_classLabel == (((Integer)kMinusOne.get(current)).intValue())))
		      current.m_items[j] = help;
                  else{
		      current.m_items[j] = help;
		      break;
		  } 
	  }
      }
      if (j == current.m_items.length) 
	newVector.addElement(current);
    }
    return newVector;
  }

    
  /**
   * Outputs the support for an item set.
   *
   * @return the support
   */
  public final int support() {

    return m_ruleSupCounter;
  }

  
  
     /**
      * Updates counter of item set with respect to given transaction.
      * @param instanceNoClass instances without the class attribute
      * @param instanceClass the values of the class attribute sorted according to instances
      */
  public final void upDateCounter(Instance instanceNoClass, Instance instanceClass) {

     
      if (containedBy(instanceNoClass)){
	  m_counter++;
	  if(this.m_classLabel == instanceClass.value(0))
	      m_ruleSupCounter++;
      }
  }

  /**
   * Updates counter of a specific item set
   * @param itemSets an item sets
   * @param instancesNoClass instances without the class attribute
   * @param instancesClass the values of the class attribute sorted according to instances
   */  
   public static void upDateCounters(FastVector itemSets, Instances instancesNoClass, Instances instancesClass){

    for (int i = 0; i < instancesNoClass.numInstances(); i++) {
        Enumeration enu = itemSets.elements();
	while (enu.hasMoreElements())
            ((LabeledItemSet)enu.nextElement()).upDateCounter(instancesNoClass.instance(i),instancesClass.instance(i));
    }

  }

   /**
    * Generates rules out of item sets
    * @param minConfidence the minimum confidence
    * @param noPrune flag indicating whether the rules are pruned accoridng to the minimum confidence value
    * @return a set of rules
    */   
    public final FastVector[] generateRules(double minConfidence, boolean noPrune) {

	FastVector premises = new FastVector(),consequences = new FastVector(),
	    conf = new FastVector();
	FastVector[] rules = new FastVector[3];
	ItemSet premise, consequence;
	
    // Generate all rules with class in the consequence. 
	premise = new ItemSet(m_totalTransactions);
	consequence = new ItemSet(m_totalTransactions);
	int[] premiseItems = new int[m_items.length];
	int[] consequenceItems = new int[1];
	System.arraycopy(m_items, 0, premiseItems, 0, m_items.length);
        consequence.setItem(consequenceItems);
        premise.setItem(premiseItems);
	consequence.setItemAt(m_classLabel,0);
	consequence.setCounter(this.m_ruleSupCounter);
	premise.setCounter(this.m_counter);
	premises.addElement(premise);
	consequences.addElement(consequence);
	conf.addElement(new Double((double)this.m_ruleSupCounter/(double)this.m_counter));
	    
	rules[0] = premises;
	rules[1] = consequences;
	rules[2] = conf;
	if(!noPrune)
	    pruneRules(rules, minConfidence);
	
    
	return rules;
    }
    
    /**
     * Returns the revision string.
     * 
     * @return		the revision
     */
    public String getRevision() {
      return RevisionUtils.extract("$Revision: 1.5 $");
    }
}
