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
 *    InstanceJoiner.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.streams;

import weka.core.Instance;
import weka.core.Instances;

import java.io.Serializable;
import java.util.Vector;

/** 
 * A bean that joins two streams of instances into one.
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 4997 $
 */
public class InstanceJoiner
  implements Serializable, InstanceProducer, SerialInstanceListener {

  /** for serialization */
  private static final long serialVersionUID = -6529972700291329656L;

  /** The listeners */
  private Vector listeners;

  /** Debugging mode */
  private boolean b_Debug;

  /** The input format for instances */
  protected Instances m_InputFormat;

  /** The current output instance */
  private Instance m_OutputInstance;


  /** Whether the first input batch has finished */
  private boolean b_FirstInputFinished;
  private boolean b_SecondInputFinished;


  /** Setup the initial states of the member variables */
  public InstanceJoiner() {
    
    listeners = new Vector();
    m_InputFormat = null;
    m_OutputInstance = null;
    b_Debug = false;
    b_FirstInputFinished = false;
    b_SecondInputFinished = false;
  }


  /**
   * Sets the format of the input instances. If the filter is able to determine
   * the output format before seeing any input instances, it does so here. This
   * default implementation assumes the output format is determined when 
   * batchFinished() is called.
   *
   * @param instanceInfo an Instances object containing the input instance
   * structure (any instances contained in the object are ignored - only the
   * structure is required).
   * @return true if the outputFormat may be collected immediately
   */
  public boolean inputFormat(Instances instanceInfo) {
    
    m_InputFormat = new Instances(instanceInfo,0);
    notifyInstanceProduced(new InstanceEvent(this, InstanceEvent.FORMAT_AVAILABLE));
    b_FirstInputFinished = false;
    b_SecondInputFinished = false;
    return true;
  }

  /**
   * Gets the format of the output instances. This should only be called
   * after input() or batchFinished() has returned true.
   *
   * @return an Instances object containing the output instance
   * structure only.
   * @throws Exception if no input structure has been defined (or the output
   * format hasn't been determined yet)
   */
  public Instances outputFormat() throws Exception {
    
    if (m_InputFormat == null) {
      throw new Exception("No output format defined.");
    }
    return new Instances(m_InputFormat,0);
  }
  
  public boolean input(Instance instance) throws Exception {
    
    if (m_InputFormat == null) {
      throw new Exception("No input instance format defined");
    }
    if (instance != null) {
      m_OutputInstance = (Instance)instance.copy();
      notifyInstanceProduced(new InstanceEvent(this,
				InstanceEvent.INSTANCE_AVAILABLE));
      return true;
    }
    return false;
  }

  /**
   * Signify that this batch of input to the filter is finished. If the filter
   * requires all instances prior to filtering, output() may now be called
   * to retrieve the filtered instances. Any subsequent instances filtered
   * should be filtered based on setting obtained from the first batch
   * (unless the inputFormat has been re-assigned or new options have been
   * set). This default implementation assumes all instance processing occurs
   * during inputFormat() and input().
   *
   * @throws Exception if no input structure has been defined
   */
  public void batchFinished() throws Exception {
    
    if (m_InputFormat == null) {
      throw new Exception("No input instance format defined");
    }
    notifyInstanceProduced(new InstanceEvent(this,
					     InstanceEvent.BATCH_FINISHED));
  }


  /**
   * Output an instance after filtering but do not remove from the output
   * queue.
   *
   * @return the instance that has most recently been filtered (or null if
   * the queue is empty).
   * @throws Exception if no input structure has been defined
   */
  public Instance outputPeek() throws Exception {
    
    if (m_InputFormat == null) {
      throw new Exception("No output instance format defined");
    }
    if (m_OutputInstance == null) {
      return null;
    }
    return (Instance)m_OutputInstance.copy();
  }


  public void setDebug(boolean debug) {
    
    b_Debug = debug;
  }
  
  public boolean getDebug() {
    
    return b_Debug;
  }

  public synchronized void addInstanceListener(InstanceListener ipl) {
    
    listeners.addElement(ipl);
  }
  
  public synchronized void removeInstanceListener(InstanceListener ipl) {
    
    listeners.removeElement(ipl);
  }
  
  protected void notifyInstanceProduced(InstanceEvent e) {
    
    if (listeners.size() > 0) {
      if (b_Debug) {
	System.err.println(this.getClass().getName()
			   + "::notifyInstanceProduced()");
      }
      Vector l;
      synchronized (this) {
	l = (Vector)listeners.clone();
      }
      for(int i = 0; i < l.size(); i++) {
	((InstanceListener)l.elementAt(i)).instanceProduced(e);
      }
      // If there are any listeners, and the event is an INSTANCE_AVAILABLE,
      // they should have retrieved the instance with outputPeek();
      try {
	if (e.getID() == InstanceEvent.INSTANCE_AVAILABLE) {
	  m_OutputInstance = null;
	}
      } catch (Exception ex) {
	System.err.println("Problem: notifyInstanceProduced() was\n"
			   + "called with INSTANCE_AVAILABLE, but output()\n"
			   + "threw an exception: " + ex.getMessage());
      }
    }
  }

  public void instanceProduced(InstanceEvent e) {
    
    Object source = e.getSource();
    if (source instanceof InstanceProducer) { 
      try {
	InstanceProducer a = (InstanceProducer) source;
	switch (e.getID()) {
	case InstanceEvent.FORMAT_AVAILABLE:
	  if (b_Debug) {
	    System.err.println(this.getClass().getName()
			+ "::firstInstanceProduced() - Format available");
	  }
	  inputFormat(a.outputFormat());
	  break;
	case InstanceEvent.INSTANCE_AVAILABLE:
	  if (b_Debug) {
	    System.err.println(this.getClass().getName()
			+ "::firstInstanceProduced() - Instance available");
	  }
	  input(a.outputPeek());
	  break;
	case InstanceEvent.BATCH_FINISHED:
	  if (b_Debug) {
	    System.err.println(this.getClass().getName()
			+ "::firstInstanceProduced() - End of instance batch");
	  }
	  batchFinished();
	  b_FirstInputFinished = true;
	  break;
	default:
	  System.err.println(this.getClass().getName()
	       + "::firstInstanceProduced() - unknown event type");
	  break;
	}
      } catch (Exception ex) {
	System.err.println(ex.getMessage());
      }
    } else {
      System.err.println(this.getClass().getName()
	     + "::firstInstanceProduced() - Unknown source object type");
    }
  }

  public void secondInstanceProduced(InstanceEvent e) {
    
    Object source = e.getSource();
    if (source instanceof InstanceProducer) { 
      try {
	if (!b_FirstInputFinished) {
	  throw new Exception(this.getClass().getName()
	  + "::secondInstanceProduced() - Input received from"
	  + " second stream before first stream finished");
	}
	InstanceProducer a = (InstanceProducer) source;
	switch (e.getID()) {
	case InstanceEvent.FORMAT_AVAILABLE:
	  if (b_Debug) {
	    System.err.println(this.getClass().getName()
	    + "::secondInstanceProduced() - Format available");
	  }
	  // Check the formats are compatible
	  if (!(a.outputFormat()).equalHeaders(outputFormat())) {
	    throw new Exception(this.getClass().getName()
	    + "::secondInstanceProduced() - incompatible instance streams\n" + (a.outputFormat()).equalHeadersMsg(outputFormat()));
	  }
	  break;
	case InstanceEvent.INSTANCE_AVAILABLE:
	  if (b_Debug) {
	    System.err.println(this.getClass().getName()
	    + "::secondInstanceProduced() - Instance available");
	  }
	  input(a.outputPeek());
	  break;
	case InstanceEvent.BATCH_FINISHED:
	  if (b_Debug) {
	    System.err.println(this.getClass().getName()
	    + "::secondInstanceProduced() - End of instance batch");
	  }
	  batchFinished();
	  break;
	default:
	  System.err.println(this.getClass().getName()
		+ "::secondInstanceProduced() - unknown event type");
	  break;
	}
      } catch (Exception ex) {
	System.err.println(ex.getMessage());
      }
    } else {
      System.err.println(this.getClass().getName()
	  + "::secondInstanceProduced() - Unknown source object type");
    }
  }
}








