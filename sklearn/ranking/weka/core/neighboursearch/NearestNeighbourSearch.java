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
 *    NearestNeighbourSearch.java
 *    Copyright (C) 1999-2007 University of Waikato
 */

package weka.core.neighboursearch;

import weka.core.AdditionalMeasureProducer;
import weka.core.DistanceFunction;
import weka.core.EuclideanDistance;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.core.Utils;

import java.io.Serializable;
import java.util.Enumeration;
import java.util.Vector;

/**
 * Abstract class for nearest neighbour search. All algorithms (classes) that
 * do nearest neighbour search should extend this class.
 *
 * @author Ashraf M. Kibriya (amk14[at-the-rate]cs[dot]waikato[dot]ac[dot]nz)
 * @version $Revision: 5953 $
 */
public abstract class NearestNeighbourSearch
  implements Serializable, OptionHandler, AdditionalMeasureProducer,
             RevisionHandler {

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
     * 
     * @param maxSize		the maximum size of the heap
     */
    public MyHeap(int maxSize) {
      if((maxSize%2)==0)
        maxSize++;
      
      m_heap = new MyHeapElement[maxSize+1];
      m_heap[0] = new MyHeapElement(0, 0);
    }
    
    /**
     * returns the size of the heap.
     * 
     * @return			the size
     */
    public int size() {
      return m_heap[0].index;
    }
    
    /**
     * peeks at the first element.
     * 
     * @return			the first element
     */
    public MyHeapElement peek() {
      return m_heap[1];
    }
    
    /**
     * returns the first element and removes it from the heap.
     * 
     * @return			the first element
     * @throws Exception	if no elements in heap
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
     * adds the value to the heap.
     * 
     * @param i			the index
     * @param d			the distance
     * @throws Exception	if the heap gets too large
     */
    public void put(int i, double d) throws Exception {
      if((m_heap[0].index+1)>(m_heap.length-1))
        throw new Exception("the number of elements cannot exceed the "+
        "initially set maximum limit");
      m_heap[0].index++;
      m_heap[m_heap[0].index] = new MyHeapElement(i, d);
      upheap();
    }
    
    /**
     * Puts an element by substituting it in place of 
     * the top most element.
     * 
     * @param i			the index
     * @param d			the distance
     * @throws Exception	if distance is smaller than that of the head
     * 				element
     */
    public void putBySubstitute(int i, double d) throws Exception {
      MyHeapElement head = get();
      put(i, d);
      //      System.out.println("previous: "+head.distance+" current: "+m_heap[1].distance);
      if(head.distance == m_heap[1].distance) { //Utils.eq(head.distance, m_heap[1].distance)) {
        putKthNearest(head.index, head.distance);
      }
      else if(head.distance > m_heap[1].distance) { //Utils.gr(head.distance, m_heap[1].distance)) {
        m_KthNearest = null;
        m_KthNearestSize = 0;
        initSize = 10;
      }
      else if(head.distance < m_heap[1].distance) {
        throw new Exception("The substituted element is smaller than the "+
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
     * @see			#m_KthNearestSize
     */
    public int noOfKthNearest() {
      return m_KthNearestSize;
    }
    
    /**
     * Stores kth nearest elements (if there are 
     * more than one).
     * @param i			the index
     * @param d			the distance
     */
    public void putKthNearest(int i,  double d) {
      if(m_KthNearest==null) {
        m_KthNearest = new MyHeapElement[initSize];
      }
      if(m_KthNearestSize>=m_KthNearest.length) {
        initSize += initSize;
        MyHeapElement temp[] = new MyHeapElement[initSize];
        System.arraycopy(m_KthNearest, 0, temp, 0, m_KthNearest.length);
        m_KthNearest = temp;
      }
      m_KthNearest[m_KthNearestSize++] = new MyHeapElement(i, d);
    }
    
    /**
     * returns the kth nearest element or null if none there.
     * 
     * @return			the kth nearest element
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
     * @return			the total size
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
    
    /** the index of this element. */
    public int index;
    
    /** the distance of this element. */
    public double distance;
    
    /**
     * constructor.
     * 
     * @param i		the index
     * @param d		the distance
     */
    public MyHeapElement(int i, double d) {
      distance = d;
      index = i;
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
   */ //better to change this into a heap element
  protected class NeighborNode
    implements RevisionHandler {

    /** The neighbor instance. */
    public Instance m_Instance;

    /** The distance from the current instance to this neighbor. */
    public double m_Distance;

    /** A link to the next neighbor instance. */
    public NeighborNode m_Next;
    
    /**
     * Create a new neighbor node.
     *
     * @param distance 		the distance to the neighbor
     * @param instance 		the neighbor instance
     * @param next 		the next neighbor node
     */
    public NeighborNode(double distance, Instance instance, NeighborNode next) {
      m_Distance = distance;
      m_Instance = instance;
      m_Next = next;
    }

    /**
     * Create a new neighbor node that doesn't link to any other nodes.
     *
     * @param distance 		the distance to the neighbor
     * @param instance 		the neighbor instance
     */
    public NeighborNode(double distance, Instance instance) {

      this(distance, instance, null);
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
   * A class for a linked list to store the nearest k neighbours
   * to an instance. We use a list so that we can take care of
   * cases where multiple neighbours are the same distance away.
   * i.e. the minimum length of the list is k.
   *
   * @author Ashraf M. Kibriya (amk14[at-the-rate]cs[dot]waikato[dot]ac[dot]nz)
   * @version $Revision: 5953 $
   */ //better to change this into a heap
  protected class NeighborList
    implements RevisionHandler {

    /** The first node in the list. */
    protected NeighborNode m_First;

    /** The last node in the list. */
    protected NeighborNode m_Last;

    /** The number of nodes to attempt to maintain in the list. */
    protected int m_Length = 1;
        
    /**
     * Creates the neighborlist with a desired length.
     *
     * @param length 		the length of list to attempt to maintain
     */
    public NeighborList(int length) {
      m_Length = length;
    }

    /**
     * Gets whether the list is empty.
     *
     * @return 			true if list is empty
     */
    public boolean isEmpty() {
      return (m_First == null);
    }

    /**
     * Gets the current length of the list.
     *
     * @return 			the current length of the list
     */
    public int currentLength() {
      int i = 0;
      NeighborNode current = m_First;
      while (current != null) {
        i++;
        current = current.m_Next;
      }
      return i;
    }

    /**
     * Inserts an instance neighbor into the list, maintaining the list
     * sorted by distance.
     *
     * @param distance 		the distance to the instance
     * @param instance 		the neighboring instance
     */
    public void insertSorted(double distance, Instance instance) {
      
      if (isEmpty()) {
        m_First = m_Last = new NeighborNode(distance, instance);
      } else {
        NeighborNode current = m_First;
        if (distance < m_First.m_Distance) {// Insert at head
          m_First = new NeighborNode(distance, instance, m_First);
        } else { // Insert further down the list
          for( ;(current.m_Next != null) &&
          (current.m_Next.m_Distance < distance);
          current = current.m_Next);
          current.m_Next = new NeighborNode(distance, instance,
          current.m_Next);
          if (current.equals(m_Last)) {
            m_Last = current.m_Next;
          }
        }
        
        // Trip down the list until we've got k list elements (or more if the
        // distance to the last elements is the same).
        int valcount = 0;
        for(current = m_First; current.m_Next != null;
        current = current.m_Next) {
          valcount++;
          if ((valcount >= m_Length) && (current.m_Distance !=
          current.m_Next.m_Distance)) {
            m_Last = current;
            current.m_Next = null;
            break;
          }
        }
      }
    }

    /**
     * Prunes the list to contain the k nearest neighbors. If there are
     * multiple neighbors at the k'th distance, all will be kept.
     *
     * @param k 		the number of neighbors to keep in the list.
     */
    public void pruneToK(int k) {
      
      if (isEmpty()) {
        return;
      }
      if (k < 1) {
        k = 1;
      }
      int currentK = 0;
      double currentDist = m_First.m_Distance;
      NeighborNode current = m_First;
      for(; current.m_Next != null; current = current.m_Next) {
        currentK++;
        currentDist = current.m_Distance;
        if ((currentK >= k) && (currentDist != current.m_Next.m_Distance)) {
          m_Last = current;
          current.m_Next = null;
          break;
        }
      }
    }
    
    /**
     * Prints out the contents of the neighborlist.
     */
    public void printList() {
      
      if (isEmpty()) {
        System.out.println("Empty list");
      } else {
        NeighborNode current = m_First;
        while (current != null) {
          System.out.println("Node: instance " + current.m_Instance
          + ", distance " + current.m_Distance);
          current = current.m_Next;
        }
        System.out.println();
      }
    }
    
    /**
     * returns the first element in the list.
     * 
     * @return			the first element
     */
    public NeighborNode getFirst() {
      return m_First;
    }
    
    /**
     * returns the last element in the list.
     * 
     * @return			the last element
     */
    public NeighborNode getLast() {
      return m_Last;
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
  
  /** The neighbourhood of instances to find neighbours in. */
  protected Instances m_Instances;
  
  /** The number of neighbours to find. */
  protected int m_kNN;

  /** the distance function used. */
  protected DistanceFunction m_DistanceFunction = new EuclideanDistance();

  /** Performance statistics. */
  protected PerformanceStats m_Stats = null;
  
  /** Should we measure Performance. */
  protected boolean m_MeasurePerformance = false;
  
  /**
   * Constructor.
   */
  public NearestNeighbourSearch() {
    if(m_MeasurePerformance)
      m_Stats = new PerformanceStats();
  }
  
  /**
   * Constructor. 
   * 
   * @param insts 	The set of instances that constitute the neighbourhood.
   */
  public NearestNeighbourSearch(Instances insts) {
    this();
    m_Instances = insts;
  }
  
  /**
   * Returns a string describing this nearest neighbour search algorithm.
   * 
   * @return 		a description of the algorithm for displaying in the 
   * 			explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "Abstract class for nearest neighbour search. All algorithms (classes) that "
      + "do nearest neighbour search should extend this class.";
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return 		an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector<Option> newVector = new Vector<Option>();

    newVector.add(new Option(
	"\tDistance function to use.\n"
	+ "\t(default: weka.core.EuclideanDistance)",
	"A", 1,"-A <classname and options>"));
    
    newVector.add(new Option(
	"\tCalculate performance statistics.",
	"P", 0,"-P"));
    
    return newVector.elements();
  }
  
  /**
   * Parses a given list of options. Valid options are:
   *
   <!-- options-start -->
   <!-- options-end -->
   *
   * @param options 	the list of options as an array of strings
   * @throws Exception 	if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String nnSearchClass = Utils.getOption('A', options);
    if(nnSearchClass.length() != 0) {
      String nnSearchClassSpec[] = Utils.splitOptions(nnSearchClass);
      if(nnSearchClassSpec.length == 0) { 
        throw new Exception("Invalid DistanceFunction specification string."); 
      }
      String className = nnSearchClassSpec[0];
      nnSearchClassSpec[0] = "";

      setDistanceFunction( (DistanceFunction)
                            Utils.forName( DistanceFunction.class, 
                                           className, nnSearchClassSpec) );
    }
    else {
      setDistanceFunction(new EuclideanDistance());
    }
    
    setMeasurePerformance(Utils.getFlag('P',options));
  }

  /**
   * Gets the current settings.
   *
   * @return 		an array of strings suitable for passing to setOptions()
   */
  public String [] getOptions() {
    Vector<String>	result;
    
    result = new Vector<String>();

    result.add("-A");
    result.add((m_DistanceFunction.getClass().getName() + " " +
                   Utils.joinOptions(m_DistanceFunction.getOptions())).trim());
    
    if(getMeasurePerformance())
      result.add("-P");
    
    return result.toArray(new String[result.size()]);
  }

  /** 
   * Returns the tip text for this property.
   * 
   * @return 		tip text for this property suitable for
   *         		displaying in the explorer/experimenter gui
   */
  public String distanceFunctionTipText() {
    return "The distance function to use for finding neighbours " +
           "(default: weka.core.EuclideanDistance). ";
  }
  
  /**
   * returns the distance function currently in use.
   * 
   * @return		the distance function
   */
  public DistanceFunction getDistanceFunction() {
    return m_DistanceFunction;
  }
  
  /**
   * sets the distance function to use for nearest neighbour search.
   * 
   * @param df		the new distance function to use
   * @throws Exception	if instances cannot be processed
   */
  public void setDistanceFunction(DistanceFunction df) throws Exception {
    m_DistanceFunction = df;
  }

  /** 
   * Returns the tip text for this property.
   * 
   * @return 		tip text for this property suitable for
   *         		displaying in the explorer/experimenter gui
   */
  public String measurePerformanceTipText() {
    return "Whether to calculate performance statistics " +
           "for the NN search or not";
  }
  
  /**
   * Gets whether performance statistics are being calculated or not.
   * 
   * @return		true if the measure performance is calculated
   */
  public boolean getMeasurePerformance() {
    return m_MeasurePerformance;
  }
  
  /**
   * Sets whether to calculate the performance statistics or not.
   * 
   * @param measurePerformance	if true then the performance is calculated
   */
  public void setMeasurePerformance(boolean measurePerformance) {
    m_MeasurePerformance = measurePerformance;
    if(m_MeasurePerformance) {
      if(m_Stats==null)
        m_Stats = new PerformanceStats();
    }
    else
      m_Stats = null;
  }
    
  /** 
   * Returns the nearest instance in the current neighbourhood to the supplied
   * instance.
   * 
   * @param target 	The instance to find the nearest neighbour for.
   * @return		the nearest neighbor
   * @throws Exception 	if the nearest neighbour could not be found.
   */
  public abstract Instance nearestNeighbour(Instance target) throws Exception;
  
  /**
   * Returns k nearest instances in the current neighbourhood to the supplied
   * instance.
   *  
   * @param target 	The instance to find the k nearest neighbours for.
   * @param k		The number of nearest neighbours to find.
   * @return		the k nearest neighbors
   * @throws Exception 	if the neighbours could not be found.
   */
  public abstract Instances kNearestNeighbours(Instance target, int k) throws Exception;
 
  /**
   * Returns the distances of the k nearest neighbours. The kNearestNeighbours
   * or nearestNeighbour needs to be called first for this to work.
   *
   * @return		the distances
   * @throws Exception 	if called before calling kNearestNeighbours
   *            	or nearestNeighbours.
   */
  public abstract double[] getDistances() throws Exception;
  
  /**
   * Updates the NearNeighbourSearch algorithm for the new added instance.
   * P.S.: The method assumes the instance has already been added to the 
   * m_Instances object by the caller.
   * 
   * @param ins		the instance to add
   * @throws Exception	if updating fails
   */
  public abstract void update(Instance ins) throws Exception;

  /** 
   * Adds information from the given instance without modifying the 
   * datastructure a lot.
   * 
   * @param ins		the instance to add the information from
   */
  public void addInstanceInfo(Instance ins) {
  }
  
  /**
   * Sets the instances.
   * 
   * @param insts	the instances to use
   * @throws Exception	if setting fails
   */
  public void setInstances(Instances insts) throws Exception {
    m_Instances = insts;
  }
  
  /** 
   * returns the instances currently set.
   * 
   * @return		the current instances
   */
  public Instances getInstances() {
    return m_Instances;
  }

  /** 
   * Gets the class object that contains the performance statistics of
   * the search method. 
   * 
   * @return		the performance statistics
   */
  public PerformanceStats getPerformanceStats() {
    return m_Stats;
  }

  /**
   * Returns an enumeration of the additional measure names.
   * 
   * @return 		an enumeration of the measure names
   */
  public Enumeration enumerateMeasures() {
    Vector<String> newVector; 
    if(m_Stats == null) {
      newVector = new Vector<String>(0);
    }
    else {
      newVector = new Vector<String>();
      Enumeration en = m_Stats.enumerateMeasures();
      while(en.hasMoreElements())
        newVector.add((String)en.nextElement());
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
    if(m_Stats==null)
      throw new IllegalArgumentException(additionalMeasureName 
        + " not supported (NearestNeighbourSearch)");
    else
      return m_Stats.getMeasure(additionalMeasureName);
  }

  /** 
   * sorts the two given arrays.
   * 
   * @param arrayToSort 	The array sorting should be based on.
   * @param linkedArray		The array that should have the same ordering as 
   * 				arrayToSort.
   */
  public static void combSort11(double arrayToSort[], int linkedArray[]) {
    int switches, j, top, gap;
    double hold1; int hold2;
    gap = arrayToSort.length;
    do {
      gap=(int)(gap/1.3);
      switch(gap) {
        case 0:
          gap = 1;
          break;
        case 9:
        case 10:
          gap=11;
          break;
        default:
          break;
      }
      switches=0;
      top = arrayToSort.length-gap;
      for(int i=0; i<top; i++) {
        j=i+gap;
        if(arrayToSort[i] > arrayToSort[j]) {
          hold1=arrayToSort[i];
          hold2=linkedArray[i];
          arrayToSort[i]=arrayToSort[j];
          linkedArray[i]=linkedArray[j];
          arrayToSort[j]=hold1;
          linkedArray[j]=hold2;
          switches++;
        }//endif
      }//endfor
    } while(switches>0 || gap>1);
  }
   
  /**
   * Partitions the instances around a pivot. Used by quicksort and
   * kthSmallestValue.
   *
   * @param arrayToSort 	the array of doubles to be sorted
   * @param linkedArray		the linked array
   * @param l 			the first index of the subset 
   * @param r 			the last index of the subset 
   * @return 			the index of the middle element
   */
  protected static int partition(double[] arrayToSort, double[] linkedArray, int l, int r) {
    double pivot = arrayToSort[(l + r) / 2];
    double help;

    while (l < r) {
      while ((arrayToSort[l] < pivot) && (l < r)) {
        l++;
      }
      while ((arrayToSort[r] > pivot) && (l < r)) {
        r--;
      }
      if (l < r) {
        help = arrayToSort[l];
        arrayToSort[l] = arrayToSort[r];
        arrayToSort[r] = help;
        help = linkedArray[l];
        linkedArray[l] = linkedArray[r];
        linkedArray[r] = help;
        l++;
        r--;
      }
    }
    if ((l == r) && (arrayToSort[r] > pivot)) {
      r--;
    } 

    return r;
  }
  
  /**
   * performs quicksort.
   * 
   * @param arrayToSort		the array to sort
   * @param linkedArray		the linked array
   * @param left 		the first index of the subset 
   * @param right		the last index of the subset 
   */
  public static void quickSort(double[] arrayToSort, double[] linkedArray, int left, int right) {
    if (left < right) {
      int middle = partition(arrayToSort, linkedArray, left, right);
      quickSort(arrayToSort, linkedArray, left, middle);
      quickSort(arrayToSort, linkedArray, middle + 1, right);
    }
  }
}
