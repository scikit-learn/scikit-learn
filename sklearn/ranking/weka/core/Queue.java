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
 *    Queue.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 *    Modified March-May 2004 by Mark Utting to add JML specs
 *    (this was done as the example solution of an assignment for a
 *     software engineering course, so the specifications are more precise
 *     and complete than one would normally do).
 *    Passed a static analysis using ESC/Java-2.0a6 with no warnings.
 */

package weka.core;

import java.io.Serializable;

/** 
 * Class representing a FIFO queue.
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 5953 $
 */
public class Queue
  extends Object
  implements Serializable, RevisionHandler {

  /** for serialization */
  private static final long serialVersionUID = -1141282001146389780L;

  /**
   * Represents one node in the queue.
   */
  protected class QueueNode
    implements Serializable, RevisionHandler {

    /** for serialization */
    private static final long serialVersionUID = -5119358279412097455L;

    /** The next node in the queue */
    protected /*@ spec_public @*/ QueueNode m_Next;

    /** The nodes contents
     */
    protected /*@ non_null spec_public @*/ Object m_Contents;

    /** 
     * Creates a queue node with the given contents 
     */
    //@ requires contents != null;
    //@ assignable m_Contents, m_Next;
    //@ ensures m_Contents == contents;
    //@ ensures m_Next == null;
    public QueueNode(Object contents) {
      m_Contents = contents;
      next(null);
    }

    /**
     * Sets the next node in the queue, and returns it.
     */
    //@ requires next != this ;
    //@ assignable m_Next;
    //@ ensures m_Next==next && \result==next;
    public QueueNode next(QueueNode next) {
      return m_Next = next;
    } //@ nowarn Invariant; // Because it stupidly checks the Queue invariant!

    /**
     * Gets the next node in the queue. 
     */
    //@ ensures \result == m_Next;
    public /*@ pure @*/ QueueNode next() {
      return m_Next;
    }

    /**
     * Sets the contents of the node.
     */
    //@ requires contents != null;
    //@ assignable m_Contents;
    //@ ensures  m_Contents == contents && \result == contents;
    public Object contents(Object contents) {
      return m_Contents = contents;
    }

    /**
     * Returns the contents in the node.
     */
      //@ ensures \result == m_Contents;
    public /*@ pure @*/ Object contents() {
      return m_Contents;
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

  /** Store a reference to the head of the queue */
  protected /*@ spec_public @*/ QueueNode m_Head = null;

  /** Store a reference to the tail of the queue */
  protected /*@ spec_public @*/ QueueNode m_Tail = null;

  /** Store the c m_Tail.m_Nexturrent number of elements in the queue */
  protected /*@ spec_public @*/ int m_Size = 0;

  //@ public invariant m_Head == null <==> m_Tail == null;
  //@public invariant m_Tail != null ==> m_Tail.m_Next == null;
  //@ public invariant m_Size >= 0;
  //@ public invariant m_Size == 0 <==> m_Head == null;
  //@ public invariant m_Size == 1 <==> m_Head != null && m_Head == m_Tail;
  //@ public invariant m_Size > 1 ==> m_Head != m_Tail;
  //@ public invariant m_Size > 1 <== m_Head != m_Tail;



  /**
   * Removes all objects from the queue m_Tail.m_Next.
   */
  //@ assignable m_Size, m_Head, m_Tail;
  //@ ensures m_Size == 0;
  //@ ensures m_Head == null;
  //@ ensures m_Tail == null;
  public final synchronized void removeAllElements() {
    m_Size = 0;
    m_Head = null;
    m_Tail = null;
  }

  /**
   * Appends an object to the back of the queue.
   *
   * @param item the object to be appended
   * @return the object appended
   */
  //@ requires item != null;
  //@ assignable m_Head, m_Tail, m_Tail.m_Next, m_Head.m_Next, m_Size;
  //@ ensures m_Head != null;
  //@ ensures m_Tail != \old(m_Tail);
  //@ ensures m_Size == \old(m_Size) + 1;
  //@ ensures \old(m_Size) == 0 ==> m_Head == m_Tail; 
  //@ ensures \old(m_Size) != 0 ==> m_Head == \old(m_Head);
  //@ ensures m_Tail.contents() == \old(item);
  //@ ensures \result == item;
  public synchronized Object push(Object item) {
    QueueNode newNode = new QueueNode(item);
    
    if (m_Head == null) {
      m_Head = m_Tail = newNode;
    } else {
      m_Tail = m_Tail.next(newNode);
    }
    m_Size++;
    return item;
  }

  /**
   * Pops an object from the front of the queue.
   *
   * @return the object at the front of the queue
   * @exception RuntimeException if the queue is empty
   */
  //@ assignable m_Head, m_Tail, m_Size;
  //@ ensures m_Size == \old(m_Size) - 1;
  //@ ensures m_Head == \old(m_Head.m_Next);
  //@ ensures m_Head != null ==> m_Tail == \old(m_Tail);
  //@ ensures \result == \old(m_Head.m_Contents);
  //@ signals (RuntimeException) \old(m_Head) == null;
  public synchronized Object pop() 
      throws RuntimeException   // REDUNDANT, BUT ESCJAVA REQUIRES THIS
  {
    if (m_Head == null) {
	throw new RuntimeException("Queue is empty");
    }
    Object retval = m_Head.contents();
    m_Size--;
    m_Head = m_Head.next();
    // Here we need to either tell ESC/Java some facts about
    // the contents of the list after popping off the head,
    // or turn off the 'invariant' warnings.
    //
    //@ assume m_Size == 0 <==> m_Head == null;
    //@ assume m_Size == 1 <==> m_Head == m_Tail;
    if (m_Head == null) {
      m_Tail = null;
    }
    return retval;
  }

  /**
   * Gets object from the front of the queue.
   *
   * @return the object at the front of the queue
   * @exception RuntimeException if the queue is empty
   */
  //@ ensures \result == \old(m_Head.m_Contents);
  //@ signals (RuntimeException) \old(m_Head) == null;
  public /*@ pure @*/ synchronized Object peek() 
    throws RuntimeException
  { 
    if (m_Head == null) {
      throw new RuntimeException("Queue is empty");
    }
    return m_Head.contents();
  }

  /**
   * Checks if queue is empty.
   * 
   * @return true if queue is empty
   */
  //@ ensures \result <==> m_Head == null;
  public /*@ pure @*/ boolean empty() {
    return m_Head == null;
  }

  /**
   * Gets queue's size.
   *
   * @return size of queue
   */
  //@ ensures \result == m_Size;
  public /*@ pure @*/ int size() {
    return m_Size;
  }

  /**
   * Produces textual description of queue.
   *
   * @return textual description of queue
   */
  //@ also
  //@ ensures \result != null;
  //@ ensures (* \result == textual description of the queue *);
  public  /*@ pure @*/ String toString() {

    String retval = "Queue Contents "+m_Size+" elements\n";
    QueueNode current = m_Head;
    if (current == null) {
      return retval + "Empty\n";
    } else {
      while (current != null) {
        retval += current.contents().toString()+"\n"; //@nowarn Modifies;
	current = current.next();
      }
    }
    return retval;
  } //@ nowarn Post;
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5953 $");
  }

  /**
   * Main method for testing this class.
   *
   * @param argv a set of strings that are pushed on a test queue
   */
  //@ requires argv.length >= 0;
  //@ requires argv != null;
  //@ requires (\forall int i; 0 <= i && i < argv.length; argv[i] != null);
  public static void main(String [] argv) {

    try {
      Queue queue = new Queue();
      for(int i = 0; i < argv.length; i++) {
	queue.push(argv[i]);
      }
      System.out.println("After pushing command line arguments");
      System.out.println(queue.toString());
      while (!queue.empty()) {
	System.out.println("Pop: " + queue.pop().toString());
      }
      // try one more pop, to make sure we get an exception
      try 
	{
	  queue.pop();
	  System.out.println("ERROR: pop did not throw exception!");
	}
      catch (RuntimeException ex)
        {
	  System.out.println("Pop on empty queue correctly gave exception.");
	}
    } catch (Exception ex) {
      System.out.println(ex.getMessage());
    }
  }
}
