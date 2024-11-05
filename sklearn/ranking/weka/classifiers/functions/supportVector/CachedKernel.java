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
 * CachedKernel.java
 * 
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 */

package weka.classifiers.functions.supportVector;

import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.Utils;

import java.util.Enumeration;
import java.util.Vector;

/**
 * Base class for RBFKernel and PolyKernel that implements a simple LRU.
 * (least-recently-used) cache if the cache size is set to a value > 0.
 * Otherwise it uses a full cache.
 * 
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @author Shane Legg (shane@intelligenesis.net) (sparse vector code)
 * @author Stuart Inglis (stuart@reeltwo.com) (sparse vector code)
 * @author J. Lindgren (jtlindgr{at}cs.helsinki.fi) (RBF kernel)
 * @author Steven Hugg (hugg@fasterlight.com) (refactored, LRU cache)
 * @author Bernhard Pfahringer (bernhard@cs.waikato.ac.nz) (full cache)
 * @version $Revision: 4549 $
 */
public abstract class CachedKernel 
  extends Kernel {

  /** for serialization */
  private static final long serialVersionUID = 702810182699015136L;
    
  /** Counts the number of kernel evaluations. */
  protected int m_kernelEvals;

  /** Counts the number of kernel cache hits. */
  protected int m_cacheHits;

  /** The size of the cache (a prime number) */
  protected int m_cacheSize = 250007;

  /** Kernel cache */
  protected double[] m_storage;
  protected long[] m_keys;

  /** The kernel matrix if full cache is used (i.e. size is set to 0) */
  protected double[][] m_kernelMatrix;

  /** The number of instance in the dataset */
  protected int m_numInsts;

  /** number of cache slots in an entry */
  protected int m_cacheSlots = 4;

  /**
   * default constructor - does nothing.
   */
  public CachedKernel() {
    super();
  }
  
  /**
   * Initializes the kernel cache. The actual size of the cache in bytes is
   * (64 * cacheSize).
   * 
   * @param data	the data to use
   * @param cacheSize	the cache size
   * @throws Exception	if something goes wrong
   */
  protected CachedKernel(Instances data, int cacheSize) throws Exception {
    super();
    
    setCacheSize(cacheSize);
    
    buildKernel(data);
  }
  
  /**
   * Returns an enumeration describing the available options.
   *
   * @return 		an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector		result;
    Enumeration		en;
    
    result = new Vector();

    en = super.listOptions();
    while (en.hasMoreElements())
      result.addElement(en.nextElement());

    result.addElement(new Option(
	"\tThe size of the cache (a prime number), 0 for full cache and \n"
	+ "\t-1 to turn it off.\n"
	+ "\t(default: 250007)",
	"C", 1, "-C <num>"));

    return result.elements();
  }

  /**
   * Parses a given list of options. <p/>
   * 
   * @param options 	the list of options as an array of strings
   * @throws Exception 	if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String	tmpStr;
    
    tmpStr = Utils.getOption('C', options);
    if (tmpStr.length() != 0)
      setCacheSize(Integer.parseInt(tmpStr));
    else
      setCacheSize(250007);
    
    super.setOptions(options);
  }

  /**
   * Gets the current settings of the Kernel.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    int       i;
    Vector    result;
    String[]  options;

    result = new Vector();
    options = super.getOptions();
    for (i = 0; i < options.length; i++)
      result.add(options[i]);

    result.add("-C");
    result.add("" + getCacheSize());

    return (String[]) result.toArray(new String[result.size()]);	  
  }

  /**
   * This method is overridden in subclasses to implement specific kernels.
   * 
   * @param id1   	the index of instance 1
   * @param id2		the index of instance 2
   * @param inst1	the instance 1 object
   * @return 		the dot product
   * @throws Exception 	if something goes wrong
   */
  protected abstract double evaluate(int id1, int id2, Instance inst1)
    throws Exception;

  /**
   * Implements the abstract function of Kernel using the cache. This method
   * uses the evaluate() method to do the actual dot product.
   *
   * @param id1 	the index of the first instance in the dataset
   * @param id2 	the index of the second instance in the dataset
   * @param inst1 	the instance corresponding to id1 (used if id1 == -1)
   * @return 		the result of the kernel function
   * @throws Exception 	if something goes wrong
   */
  public double eval(int id1, int id2, Instance inst1) throws Exception {
		
    double result = 0;
    long key = -1;
    int location = -1;

    // we can only cache if we know the indexes and caching is not 
    // disbled (m_cacheSize == -1)
    if ( (id1 >= 0) && (m_cacheSize != -1) ) {

      // Use full cache?
      if (m_cacheSize == 0) {
	if (m_kernelMatrix == null) {
	  m_kernelMatrix = new double[m_data.numInstances()][];
	  for(int i = 0; i < m_data.numInstances(); i++) {
	    m_kernelMatrix[i] = new double[i + 1];
	    for(int j = 0; j <= i; j++) {
	      m_kernelEvals++;
	      m_kernelMatrix[i][j] = evaluate(i, j, m_data.instance(i));
	    }
	  }
	} 
	m_cacheHits++;
	result = (id1 > id2) ? m_kernelMatrix[id1][id2] : m_kernelMatrix[id2][id1];
	return result;
      }

      // Use LRU cache
      if (id1 > id2) {
	key = (id1 + ((long) id2 * m_numInsts));
      } else {
	key = (id2 + ((long) id1 * m_numInsts));
      }
      location = (int) (key % m_cacheSize) * m_cacheSlots;
      int loc = location;
      for (int i = 0; i < m_cacheSlots; i++) {
	long thiskey = m_keys[loc];
	if (thiskey == 0)
	  break; // empty slot, so break out of loop early
	if (thiskey == (key + 1)) {
	  m_cacheHits++;
	  // move entry to front of cache (LRU) by swapping
	  // only if it's not already at the front of cache
	  if (i > 0) {
	    double tmps = m_storage[loc];
	    m_storage[loc] = m_storage[location];
	    m_keys[loc] = m_keys[location];
	    m_storage[location] = tmps;
	    m_keys[location] = thiskey;
	    return tmps;
	  } else
	    return m_storage[loc];
	}
	loc++;
      }
    }

    result = evaluate(id1, id2, inst1);

    m_kernelEvals++;

    // store result in cache
    if ( (key != -1) && (m_cacheSize != -1) ) {
      // move all cache slots forward one array index
      // to make room for the new entry
      System.arraycopy(m_keys, location, m_keys, location + 1,
		       m_cacheSlots - 1);
      System.arraycopy(m_storage, location, m_storage, location + 1,
		       m_cacheSlots - 1);
      m_storage[location] = result;
      m_keys[location] = (key + 1);
    }
    return result;
  }

  /**
   * Returns the number of time Eval has been called.
   * 
   * @return 		the number of kernel evaluation.
   */
  public int numEvals() {
    return m_kernelEvals;
  }

  /**
   * Returns the number of cache hits on dot products.
   * 
   * @return 		the number of cache hits.
   */
  public int numCacheHits() {
    return m_cacheHits;
  }

  /**
   * Frees the cache used by the kernel.
   */
  public void clean() {
    m_storage = null;
    m_keys = null;
    m_kernelMatrix = null;
  }

  /**
   * Calculates a dot product between two instances
   * 
   * @param inst1	the first instance
   * @param inst2	the second instance
   * @return 		the dot product of the two instances.
   * @throws Exception	if an error occurs
   */
  protected final double dotProd(Instance inst1, Instance inst2)
    throws Exception {

    double result = 0;

    // we can do a fast dot product
    int n1 = inst1.numValues();
    int n2 = inst2.numValues();
    int classIndex = m_data.classIndex();
    for (int p1 = 0, p2 = 0; p1 < n1 && p2 < n2;) {
      int ind1 = inst1.index(p1);
      int ind2 = inst2.index(p2);
      if (ind1 == ind2) {
	if (ind1 != classIndex) {
	  result += inst1.valueSparse(p1) * inst2.valueSparse(p2);
	}
	p1++;
	p2++;
      } else if (ind1 > ind2) {
	p2++;
      } else {
	p1++;
      }
    }
    return (result);
  }

  /**
   * Sets the size of the cache to use (a prime number)
   * 
   * @param value	the size of the cache
   */
  public void setCacheSize(int value) {
    if (value >= -1) {
      m_cacheSize = value;
      clean();
    }
    else {
      System.out.println(
	  "Cache size cannot be smaller than -1 (provided: " + value + ")!");
    }
  }
  
  /**
   * Gets the size of the cache
   * 
   * @return 		the cache size
   */
  public int getCacheSize() {
    return m_cacheSize;
  }

  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String cacheSizeTipText() {
    return "The size of the cache (a prime number), 0 for full cache and -1 to turn it off.";
  }

  /**
   * initializes variables etc.
   * 
   * @param data	the data to use
   */
  protected void initVars(Instances data) {
    super.initVars(data);
    
    m_kernelEvals = 0;
    m_cacheHits   = 0;
    m_numInsts    = m_data.numInstances();

    if (getCacheSize() > 0) {
      // Use LRU cache
      m_storage = new double[m_cacheSize * m_cacheSlots];
      m_keys    = new long[m_cacheSize * m_cacheSlots];
    } 
    else {
      m_storage      = null;
      m_keys         = null;
      m_kernelMatrix = null;
    }
  }
  
  /**
   * builds the kernel with the given data. Initializes the kernel cache. 
   * The actual size of the cache in bytes is (64 * cacheSize).
   * 
   * @param data	the data to base the kernel on
   * @throws Exception	if something goes wrong
   */
  public void buildKernel(Instances data) throws Exception {
    // does kernel handle the data?
    if (!getChecksTurnedOff())
      getCapabilities().testWithFail(data);

    initVars(data);
  }
}
