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

/**
 *    KStarCache.java
 *    Copyright (C) 1995 University of Waikato
 *    Java port to Weka by Abdelaziz Mahoui (am14@cs.waikato.ac.nz).
 *
 */


package weka.classifiers.lazy.kstar;

import weka.core.RevisionHandler;
import weka.core.RevisionUtils;

import java.io.Serializable;

/**
 * A class representing the caching system used to keep track of each attribute
 * value and its corresponding scale factor or stop parameter.
 *
 * @author Len Trigg (len@reeltwo.com)
 * @author Abdelaziz Mahoui (am14@cs.waikato.ac.nz)
 * @version $Revision: 1.11 $
 */
public class KStarCache
  implements Serializable, RevisionHandler {

  /** for serialization */
  private static final long serialVersionUID = -7693632394267140678L;
  
  /**
   * cache table
   */
  CacheTable m_Cache = new CacheTable();
  
  /**
   * Stores the specified values in the cahce table for easy retrieval.
   *
   * @param key attribute value used key to lookup the cache table.
   * @param value cache parameter: attribute scale/stop parameter.
   * @param pmiss cache parameter: transformation probability to 
   * attribute with missing value.
   */
  public void store(double key, double value, double pmiss) {
    if ( !m_Cache.containsKey(key) ) {
      m_Cache.insert(key, value, pmiss);
    }
  }
  
  /**
   * Checks if the specified key maps with an entry in the cache table
   *
   * @param key the key to map with an entry in the hashtable.
   */
  public boolean containsKey(double key) {
    if ( m_Cache.containsKey(key) ) {
      return true;
    }
    return false;
  }
  
  /**
   * Returns the values in the cache mapped by the specified key
   *
   * @param key the key used to retrieve the table entry.
   */
  public TableEntry getCacheValues( double key ) {
    if ( m_Cache.containsKey(key) ) {
      return m_Cache.getEntry(key);
    }
    return null;
  }

  /**
   * A custom hashtable class to support the caching system.
   *
   */
  public class CacheTable
    implements Serializable, RevisionHandler {

    /** for serialization */
    private static final long serialVersionUID = -8086106452588253423L;

    /** The hash table data. */
    private TableEntry [] m_Table;

    /** The total number of entries in the hash table. */
    private int m_Count;

    /** Rehashes the table when count exceeds this threshold. */
    private int m_Threshold;

    /** The load factor for the hashtable. */
    private float m_LoadFactor;

    /** The default size of the hashtable */
    private final int DEFAULT_TABLE_SIZE = 101;

    /** The default load factor for the hashtable */
    private final float DEFAULT_LOAD_FACTOR = 0.75f;

    //    private final float DEFAULT_LOAD_FACTOR = 0.5f;
    /** Accuracy value for equality */
    private final double EPSILON = 1.0E-5;
    
    /**
     * Constructs a new hashtable with a default capacity and load factor.
     */
    public CacheTable(int size, float loadFactor) {
      m_Table = new TableEntry[size];
      m_LoadFactor = loadFactor;
      m_Threshold = (int)(size * loadFactor);
      m_Count = 0;
    }
    
    /**
     * Constructs a new hashtable with a default capacity and load factor.
     */
    public CacheTable() {
      this(101, 0.75f);
    }
    
    /**
     * Tests if the specified double is a key in this hashtable.
     */
    public boolean containsKey(double key) {
      TableEntry [] table = m_Table;
      int hash = hashCode(key);
      int index = (hash & 0x7FFFFFFF) % table.length;
      for (TableEntry e = table[index] ; e != null ; e = e.next) {
	if ((e.hash == hash) && (Math.abs(e.key - key) < EPSILON)) {
	  return true;
	}
      }
      return false;
    }
    
    /**
     * Inserts a new entry in the hashtable using the specified key. 
     * If the key already exist in the hashtable, do nothing.
     */
    public void insert(double key, double value, double pmiss) {
      // Makes sure the key is not already in the hashtable.
      TableEntry e, ne;
      TableEntry [] table = m_Table;
      int hash = hashCode(key);
      int index = (hash & 0x7FFFFFFF) % table.length;
      // start looking along the chain
      for (e = table[index] ; e != null ; e = e.next) {
	if ((e.hash == hash) && (Math.abs(e.key - key) < EPSILON)) {
	  return;
	}
      }
      // At this point, key is not in table.
      // Creates a new entry.
      ne = new TableEntry( hash, key, value, pmiss, table[index] );
      // Put entry at the head of the chain.
      table[index] = ne;
      m_Count++;
      // Rehash the table if the threshold is exceeded
      if (m_Count >= m_Threshold) {
	rehash();
      }
    }
    
    /**
     * Returns the table entry to which the specified key is mapped in 
     * this hashtable.
     * @return a table entry.
     */
    public TableEntry getEntry(double key) {
      TableEntry [] table = m_Table;
      int hash = hashCode(key);
      int index = (hash & 0x7FFFFFFF) % table.length;
      for (TableEntry e = table[index] ; e != null ; e = e.next) {
	if ((e.hash == hash) && (Math.abs(e.key - key) < EPSILON)) {
	  return e;
	}
      }
      return null;
    }
    
    /**
     * Returns the number of keys in this hashtable.
     * @return the number of keys in this hashtable.
     */
    public int size() {
      return m_Count;
    }
    
    /**
     * Tests if this hashtable maps no keys to values.
     * @return true if this hastable maps no keys to values.
     */
    public boolean isEmpty() {
      return m_Count == 0;
    }
    
    /**
     * Clears this hashtable so that it contains no keys.
     */
    public void clear() {
      TableEntry table[] = m_Table;
      for (int index = table.length; --index >= 0; ) {
	table[index] = null;
      }
      m_Count = 0;
    }
    
    /**
     * Rehashes the contents of the hashtable into a hashtable with a 
     * larger capacity. This method is called automatically when the 
     * number of keys in the hashtable exceeds this hashtable's capacity 
     * and load factor. 
     */
    private void rehash() {
      int oldCapacity = m_Table.length;
      TableEntry [] oldTable = m_Table;    
      int newCapacity = oldCapacity * 2 + 1;
      TableEntry [] newTable = new TableEntry[newCapacity];
      m_Threshold = (int)(newCapacity * m_LoadFactor);
      m_Table = newTable;
      TableEntry e, old;
      for (int i = oldCapacity ; i-- > 0 ;) {
	for (old = oldTable[i] ; old != null ; ) {
	  e = old;
	  old = old.next;
	  int index = (e.hash & 0x7FFFFFFF) % newCapacity;
	  e.next = newTable[index];
	  newTable[index] = e;
	}
      }
    }
    
    /**
     * Returns the hash code of the specified double.
     * @return the hash code of the specified double.
     */
    private int hashCode(double key) {
      long bits = Double.doubleToLongBits(key);
      return (int)(bits ^ (bits >> 32));
    }
    
    /**
     * Returns the revision string.
     * 
     * @return		the revision
     */
    public String getRevision() {
      return RevisionUtils.extract("$Revision: 1.11 $");
    }
  } // CacheTable
  
  /**
   * Hashtable collision list.
   */
  public class TableEntry
    implements Serializable, RevisionHandler {

    /** for serialization */
    private static final long serialVersionUID = 4057602386766259138L;

    /** attribute value hash code */
    public int hash;

    /** attribute value */
    public double key;

    /** scale factor or stop parameter */
    public double value;

    /** transformation probability to missing value */
    public double pmiss;

    /** next table entry (separate chaining) */
    public TableEntry next = null;

    /** Constructor */
    public TableEntry(int hash, double key, double value, 
		      double pmiss, TableEntry next) {
      this.hash  = hash;
      this.key   = key;
      this.value = value;
      this.pmiss = pmiss;
      this.next  = next;
    }
    
    /**
     * Returns the revision string.
     * 
     * @return		the revision
     */
    public String getRevision() {
      return RevisionUtils.extract("$Revision: 1.11 $");
    }
  }  // TableEntry
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 1.11 $");
  }
  
} // Cache
