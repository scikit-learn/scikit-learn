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
 *    Filter.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

import weka.core.Instance;
import weka.core.Instances;
import weka.core.OptionHandler;
import weka.core.Utils;
import weka.filters.AllFilter;
import weka.filters.StreamableFilter;
import weka.filters.SupervisedFilter;
import weka.gui.Logger;

import java.awt.BorderLayout;
import java.beans.EventSetDescriptor;
import java.io.Serializable;
import java.util.Enumeration;
import java.util.EventObject;
import java.util.Hashtable;
import java.util.Vector;

import javax.swing.JPanel;

/**
 * A wrapper bean for Weka filters
 *
 * @author <a href="mailto:mhall@cs.waikato.ac.nz">Mark Hall</a>
 * @version $Revision: 6458 $
 */
public class Filter
  extends JPanel
  implements BeanCommon, Visible, WekaWrapper,
	     Serializable, UserRequestAcceptor,
	     TrainingSetListener, TestSetListener,
	     TrainingSetProducer, TestSetProducer,
	     DataSource, DataSourceListener, 
	     InstanceListener, EventConstraints,
	     ConfigurationProducer {

  /** for serialization */
  private static final long serialVersionUID = 8249759470189439321L;

  protected BeanVisual m_visual = 
    new BeanVisual("Filter",
		   BeanVisual.ICON_PATH+"DefaultFilter.gif",
		   BeanVisual.ICON_PATH+"DefaultFilter_animated.gif");

  private static int IDLE = 0;
  private static int FILTERING_TRAINING = 1;
  private static int FILTERING_TEST = 2;
  private int m_state = IDLE;

  protected Thread m_filterThread = null;

  private transient Instances m_trainingSet;
  private transient Instances m_testingSet;

  /**
   * Global info for the wrapped filter (if it exists).
   */
  protected String m_globalInfo;

  /**
   * Objects talking to us
   */
  private Hashtable m_listenees = new Hashtable();

  /**
   * Objects listening for training set events
   */
  private Vector m_trainingListeners = new Vector();

  /**
   * Objects listening for test set events
   */
  private Vector m_testListeners = new Vector();

  /**
   * Objects listening for instance events
   */
  private Vector m_instanceListeners = new Vector();

  /**
   * Objects listening for data set events
   */
  private Vector m_dataListeners = new Vector();

  /**
   * The filter to use.
   */
  private weka.filters.Filter m_Filter = new AllFilter();

  /**
   * Instance event object for passing on filtered instance streams
   */
  private InstanceEvent m_ie = new InstanceEvent(this);

  /**
   * Logging.
   */
  private transient Logger m_log = null;
  
  /**
   * Counts incoming streamed instances.
   */
  private transient int m_instanceCount;
  
  /**
   * Global info (if it exists) for the wrapped filter
   *
   * @return the global info
   */
  public String globalInfo() {
    return m_globalInfo;
  }

  public Filter() {
    setLayout(new BorderLayout());
    add(m_visual, BorderLayout.CENTER);
    setFilter(m_Filter);
  }

  /**
   * Set a custom (descriptive) name for this bean
   * 
   * @param name the name to use
   */
  public void setCustomName(String name) {
    m_visual.setText(name);
  }

  /**
   * Get the custom (descriptive) name for this bean (if one has been set)
   * 
   * @return the custom name (or the default name)
   */
  public String getCustomName() {
    return m_visual.getText();
  }

  /**
   * Set the filter to be wrapped by this bean
   *
   * @param c a <code>weka.filters.Filter</code> value
   */
  public void setFilter(weka.filters.Filter c) {
    boolean loadImages = true;
    if (c.getClass().getName().
	compareTo(m_Filter.getClass().getName()) == 0) {
      loadImages = false;
    }
    m_Filter = c;
    String filterName = c.getClass().toString();
    filterName = filterName.substring(filterName.
				      indexOf('.')+1, 
				      filterName.length());
    if (loadImages) {
      if (m_Filter instanceof Visible) {
        m_visual = ((Visible) m_Filter).getVisual();
      } else {
        if (!m_visual.loadIcons(BeanVisual.ICON_PATH+filterName+".gif",
                                BeanVisual.ICON_PATH+filterName+"_animated.gif")) {
          useDefaultVisual();
        }
      }
    }
    m_visual.setText(filterName.substring(filterName.lastIndexOf('.')+1,
					  filterName.length()));

    if (m_Filter instanceof LogWriter && m_log != null) {
      ((LogWriter) m_Filter).setLog(m_log);
    }

    if (!(m_Filter instanceof StreamableFilter) &&
	(m_listenees.containsKey("instance"))) {
      if (m_log != null) {
	m_log.logMessage("[Filter] " + 
	    statusMessagePrefix() + " WARNING : "
	    + m_Filter.getClass().getName()
	    + " is not an incremental filter");
	m_log.statusMessage(statusMessagePrefix()
	    + "WARNING: Not an incremental filter.");
      }
    }
    
    // get global info
    m_globalInfo = KnowledgeFlowApp.getGlobalInfo(m_Filter);
  }

  public weka.filters.Filter getFilter() {
    return m_Filter;
  }

  /**
   * Set the filter to be wrapped by this bean
   *
   * @param algorithm a weka.filters.Filter
   * @exception IllegalArgumentException if an error occurs
   */
  public void setWrappedAlgorithm(Object algorithm) {
    
    if (!(algorithm instanceof weka.filters.Filter)) { 
      throw new IllegalArgumentException(algorithm.getClass()+" : incorrect "
					 +"type of algorithm (Filter)");
    }
    setFilter((weka.filters.Filter)algorithm);
  }

  /**
   * Get the filter wrapped by this bean
   *
   * @return an <code>Object</code> value
   */
  public Object getWrappedAlgorithm() {
    return getFilter();
  }

  /**
   * Accept a training set
   *
   * @param e a <code>TrainingSetEvent</code> value
   */
  public void acceptTrainingSet(TrainingSetEvent e) {
    processTrainingOrDataSourceEvents(e);
  }

  private boolean m_structurePassedOn = false;
  /**
   * Accept an instance for processing by StreamableFilters only
   *
   * @param e an <code>InstanceEvent</code> value
   */
  public void acceptInstance(InstanceEvent e) {
    // to do!
    if (m_filterThread != null) {
      String messg = "[Filter] " + statusMessagePrefix() 
        + " is currently batch processing!";
      if (m_log != null) {
	m_log.logMessage(messg);
	m_log.statusMessage(statusMessagePrefix()
	    + "WARNING: Filter is currently batch processing.");
      } else {
	System.err.println(messg);
      }
      return;
    }
    if (!(m_Filter instanceof StreamableFilter)) {
      stop(); // stop all processing
      if (m_log != null) {
	m_log.logMessage("[Filter] " + statusMessagePrefix() 
	    + " ERROR : "+m_Filter.getClass().getName()
	    +"can't process streamed instances; can't continue");
	m_log.statusMessage(statusMessagePrefix()
	    + "ERROR: Can't process streamed instances; can't continue.");
      }
      return;
    }
    if (e.getStatus() == InstanceEvent.FORMAT_AVAILABLE) {
      try {
        m_instanceCount = 0;
        //notifyInstanceListeners(e);
	//	Instances dataset = e.getInstance().dataset();
	Instances dataset = e.getStructure();
	if (m_Filter instanceof SupervisedFilter) {
	  // defualt to last column if no class is set
	  if (dataset.classIndex() < 0) {
	    dataset.setClassIndex(dataset.numAttributes()-1);
	  }
	}
	// initialize filter
	m_Filter.setInputFormat(dataset);
	// attempt to determine post-filtering
	// structure. If successful this can be passed on to instance
	// listeners as a new FORMAT_AVAILABLE event.
	m_structurePassedOn = false;
	try {
	  if (m_Filter.isOutputFormatDefined()) {
//	    System.err.println("Filter - passing on output format...");
	    //	    System.err.println(m_Filter.getOutputFormat());
	    m_ie.setStructure(m_Filter.getOutputFormat());
	    notifyInstanceListeners(m_ie);
	    m_structurePassedOn = true;
	  }
	} catch (Exception ex) {
	  stop(); // stop all processing
	  if (m_log != null) {
	    m_log.logMessage("[Filter] " + statusMessagePrefix() 
                + " Error in obtaining post-filter structure. " 
                + ex.getMessage());
	    m_log.statusMessage(statusMessagePrefix()
	        +"ERROR (See log for details).");
	  } else {
	    System.err.println("[Filter] " + statusMessagePrefix() 
	        + " Error in obtaining post-filter structure");
	  }
	}
      } catch (Exception ex) {
	ex.printStackTrace();
      }
      return;
    }
   
    if (e.getStatus() == InstanceEvent.BATCH_FINISHED) {
      // get the last instance (if available)
      try {
        if (m_log != null) {
          m_log.statusMessage(statusMessagePrefix() 
              + "Stream finished.");
        }
        if (m_Filter.input(e.getInstance())) {
          Instance filteredInstance = m_Filter.output();
          if (filteredInstance != null) {
            if (!m_structurePassedOn) {
              // pass on the new structure first
              m_ie.setStructure(new Instances(filteredInstance.dataset(), 0));
              notifyInstanceListeners(m_ie);
              m_structurePassedOn = true;
            }

            m_ie.setInstance(filteredInstance);
            
            // if there are instances pending for output don't want to send
            // a batch finisehd at this point...
            //System.err.println("Filter - in batch finisehd...");
            if (m_Filter.batchFinished() && m_Filter.numPendingOutput() > 0) {
              m_ie.setStatus(InstanceEvent.INSTANCE_AVAILABLE);
            } else {
              m_ie.setStatus(e.getStatus());
            }
            notifyInstanceListeners(m_ie);
          }
        }
        if (m_log != null) {
          m_log.statusMessage(statusMessagePrefix() + "Done.");
        }
      } catch (Exception ex) {
        stop(); // stop all processing
        if (m_log != null) {
          m_log.logMessage("[Filter] " 
              + statusMessagePrefix() + ex.getMessage());
          m_log.statusMessage(statusMessagePrefix()
              + "ERROR (See log for details).");
        }
        ex.printStackTrace();
      }
      
      // check for any pending instances that we might need to pass on
      try {
        if (m_Filter.batchFinished() && m_Filter.numPendingOutput() > 0) {
          if (m_log != null) {
            m_log.statusMessage(statusMessagePrefix() 
                + "Passing on pending instances...");
          }
          Instance filteredInstance = m_Filter.output();
          if (filteredInstance != null) {
            if (!m_structurePassedOn) {
              // pass on the new structure first
              m_ie.setStructure(new Instances(filteredInstance.dataset(), 0));
              notifyInstanceListeners(m_ie);
              m_structurePassedOn = true;
            }

            m_ie.setInstance(filteredInstance);
            
            m_ie.setStatus(InstanceEvent.INSTANCE_AVAILABLE);
            notifyInstanceListeners(m_ie);
          }
          while (m_Filter.numPendingOutput() > 0) {
            filteredInstance = m_Filter.output();
            m_ie.setInstance(filteredInstance);
//            System.err.println("Filter - sending pending...");
            if (m_Filter.numPendingOutput() == 0) {
              m_ie.setStatus(InstanceEvent.BATCH_FINISHED);
            } else {
              m_ie.setStatus(InstanceEvent.INSTANCE_AVAILABLE);
            }
            notifyInstanceListeners(m_ie);
          }
          if (m_log != null) {
            m_log.statusMessage(statusMessagePrefix() + "Finished.");
          }
        }
      } catch (Exception ex) {
        stop(); // stop all processing
        if (m_log != null) {
          m_log.logMessage("[Filter] " 
              + statusMessagePrefix() 
              + ex.toString());
          m_log.statusMessage(statusMessagePrefix()
              + "ERROR (See log for details.");
        }
        ex.printStackTrace();
      }
    } else {
      // pass instance through the filter
      try {
        if (!m_Filter.input(e.getInstance())) {
//          System.err.println("Filter - inputing instance into filter...");
          /* if (m_log != null) {
            m_log.logMessage("ERROR : filter not ready to output instance");
          } */
          
          // quietly return. Filter might be able to output some instances
          // once the batch is finished.
          return;
        }

        // collect output instance.
        Instance filteredInstance = m_Filter.output();
        if (filteredInstance == null) {
          return;
        }
        m_instanceCount++;
        
        if (!m_structurePassedOn) {
          // pass on the new structure first
          m_ie.setStructure(new Instances(filteredInstance.dataset(), 0));
          notifyInstanceListeners(m_ie);
          m_structurePassedOn = true;
        }

        m_ie.setInstance(filteredInstance);
        m_ie.setStatus(e.getStatus());
        
        if (m_log != null && (m_instanceCount % 10000 == 0)) {
          m_log.statusMessage(statusMessagePrefix()
              + "Received " + m_instanceCount + " instances.");
        }
        notifyInstanceListeners(m_ie);
      } catch (Exception ex) {
        stop(); // stop all processing
        if (m_log != null) {
          m_log.logMessage("[Filter] " + statusMessagePrefix() 
              + ex.toString());
          m_log.statusMessage(statusMessagePrefix()
              + "ERROR (See log for details).");
        }
        ex.printStackTrace();
      }
    }
  }

  private void processTrainingOrDataSourceEvents(final EventObject e) {
    boolean structureOnly = false;
    if (e instanceof DataSetEvent) {
      structureOnly = ((DataSetEvent)e).isStructureOnly();
      if(structureOnly){
       notifyDataOrTrainingListeners(e);  
      }
    }
    if (e instanceof TrainingSetEvent) {
      structureOnly = ((TrainingSetEvent)e).isStructureOnly();
      if(structureOnly){
       notifyDataOrTrainingListeners(e);  
      }
    }
    if (structureOnly && !(m_Filter instanceof StreamableFilter)) {
      return; // nothing can be done
    }
    
    if (m_filterThread == null) {
      try {
	if (m_state == IDLE) {
	  synchronized (this) {
	    m_state = FILTERING_TRAINING;
	  }
	  m_trainingSet = (e instanceof TrainingSetEvent) 
	    ? ((TrainingSetEvent)e).getTrainingSet()
	    : ((DataSetEvent)e).getDataSet();

//	  final String oldText = m_visual.getText();
	  m_filterThread = new Thread() {
	      public void run() {
		try {
		  if (m_trainingSet != null) {
		    m_visual.setAnimated();
//		    m_visual.setText("Filtering training data...");
		    if (m_log != null) {
		      m_log.statusMessage(statusMessagePrefix() 
		          + "Filtering training data ("
		          + m_trainingSet.relationName() + ")");
		    }
		    m_Filter.setInputFormat(m_trainingSet);
		    Instances filteredData = 
		      weka.filters.Filter.useFilter(m_trainingSet, m_Filter);
//		    m_visual.setText(oldText);
		    m_visual.setStatic();
		    EventObject ne;
		    if (e instanceof TrainingSetEvent) {
		      ne = new TrainingSetEvent(weka.gui.beans.Filter.this, 
						filteredData);
		      ((TrainingSetEvent)ne).m_setNumber =
			((TrainingSetEvent)e).m_setNumber;
		      ((TrainingSetEvent)ne).m_maxSetNumber = 
			((TrainingSetEvent)e).m_maxSetNumber;
		    } else {
		      ne = new DataSetEvent(weka.gui.beans.Filter.this,
					    filteredData);
		    }

		    notifyDataOrTrainingListeners(ne);
		  }
		} catch (Exception ex) {
		  ex.printStackTrace();
                  if (m_log != null) {
                    m_log.logMessage("[Filter] " + statusMessagePrefix() 
                        + ex.getMessage());
                    m_log.statusMessage(statusMessagePrefix()
                        + "ERROR (See log for details).");
//                    m_log.statusMessage("Problem filtering: see log for details.");
                  }
                  Filter.this.stop(); // stop all processing
		} finally {
//		  m_visual.setText(oldText);
		  m_visual.setStatic();
		  m_state = IDLE;
		  if (isInterrupted()) {
		    m_trainingSet = null;
		    if (m_log != null) {
		      m_log.logMessage("[Filter] " + statusMessagePrefix()
                                       + " training set interrupted!");
		      m_log.statusMessage(statusMessagePrefix()
		          + "INTERRUPTED");
		    }		    
		  } else {
		    if (m_log != null) {
		      m_log.statusMessage(statusMessagePrefix() + "Finished.");
		    }
		  }
		  block(false);
		}
	      }
	    };
	  m_filterThread.setPriority(Thread.MIN_PRIORITY);
	  m_filterThread.start();
	  block(true);
	  m_filterThread = null;
	  m_state = IDLE;
	}
      } catch (Exception ex) {
	ex.printStackTrace();
      }
    }
  }

  /**
   * Accept a test set
   *
   * @param e a <code>TestSetEvent</code> value
   */
  public void acceptTestSet(final TestSetEvent e) {
      if(e.isStructureOnly())
            notifyTestListeners(e);
      if (m_trainingSet != null && 
	m_trainingSet.equalHeaders(e.getTestSet()) && 
	m_filterThread == null) {
      try {
	if (m_state == IDLE) {
	  m_state = FILTERING_TEST;
	}
	m_testingSet = e.getTestSet();
        //	final String oldText = m_visual.getText();
	m_filterThread = new Thread() {
	    public void run() {
	      try {
		if (m_testingSet != null) {
		  m_visual.setAnimated();
                  //		  m_visual.setText("Filtering test data...");
		  if (m_log != null) {
		    m_log.statusMessage(statusMessagePrefix() 
		        + "Filtering test data ("
			+ m_testingSet.relationName() + ")");
		  }
		  Instances filteredTest = 
		    weka.filters.Filter.useFilter(m_testingSet, m_Filter);
                  //		  m_visual.setText(oldText);
		  m_visual.setStatic();
		  TestSetEvent ne =
		    new TestSetEvent(weka.gui.beans.Filter.this,
				     filteredTest);
		  ne.m_setNumber = e.m_setNumber;
		  ne.m_maxSetNumber = e.m_maxSetNumber;
		  notifyTestListeners(ne);
		}
	      } catch (Exception ex) {
		ex.printStackTrace();
                if (m_log != null) {
                  m_log.logMessage("[Filter] " + statusMessagePrefix() 
                      + ex.getMessage());
                  m_log.statusMessage(statusMessagePrefix() 
                      + "ERROR (See log for details).");
                }
                Filter.this.stop();
	      } finally {
                //		m_visual.setText(oldText);
		m_visual.setStatic();
		m_state = IDLE;
		if (isInterrupted()) {
		  m_trainingSet = null;
		  if (m_log != null) {
		      m_log.logMessage("[Filter] " + statusMessagePrefix()
                                       + " test set interrupted!");
		      m_log.statusMessage(statusMessagePrefix()
		          + "INTERRUPTED");
//		    m_log.statusMessage("OK");
		  }
		} else {
		  if (m_log != null) {
		    m_log.statusMessage(statusMessagePrefix() + "Finished.");
		  }
		}
		block(false);
	      }
	    }
	  };
	m_filterThread.setPriority(Thread.MIN_PRIORITY);
	m_filterThread.start();
	block(true);
	m_filterThread = null;
	m_state = IDLE;
      } catch (Exception ex) {
	ex.printStackTrace();
      }
    }
  }

  /**
   * Accept a data set
   *
   * @param e a <code>DataSetEvent</code> value
   */
  public void acceptDataSet(DataSetEvent e) {
    processTrainingOrDataSourceEvents(e);
  }

  /**
   * Set the visual appearance of this bean
   *
   * @param newVisual a <code>BeanVisual</code> value
   */
  public void setVisual(BeanVisual newVisual) {
    m_visual = newVisual;
  }

  /**
   * Get the visual appearance of this bean
   *
   * @return a <code>BeanVisual</code> value
   */
  public BeanVisual getVisual() {
    return m_visual;
  }

  /**
   * Use the default visual appearance
   */
  public void useDefaultVisual() {
    m_visual.loadIcons(BeanVisual.ICON_PATH+"DefaultFilter.gif",
		       BeanVisual.ICON_PATH+"DefaultFilter_animated.gif");
  }

  /**
   * Add a training set listener
   *
   * @param tsl a <code>TrainingSetListener</code> value
   */
  public synchronized void addTrainingSetListener(TrainingSetListener tsl) {
    m_trainingListeners.addElement(tsl);
  }
  
  /**
   * Remove a training set listener
   *
   * @param tsl a <code>TrainingSetListener</code> value
   */
  public synchronized void removeTrainingSetListener(TrainingSetListener tsl) {
     m_trainingListeners.removeElement(tsl);
  }

  /**
   * Add a test set listener
   *
   * @param tsl a <code>TestSetListener</code> value
   */
  public synchronized void addTestSetListener(TestSetListener tsl) {
    m_testListeners.addElement(tsl);
  }
  
  /**
   * Remove a test set listener
   *
   * @param tsl a <code>TestSetListener</code> value
   */
  public synchronized void removeTestSetListener(TestSetListener tsl) {
    m_testListeners.removeElement(tsl);
  }

  /**
   * Add a data source listener
   *
   * @param dsl a <code>DataSourceListener</code> value
   */
  public synchronized void addDataSourceListener(DataSourceListener dsl) {
    m_dataListeners.addElement(dsl);
  }

  /**
   * Remove a data source listener
   *
   * @param dsl a <code>DataSourceListener</code> value
   */
  public synchronized void removeDataSourceListener(DataSourceListener dsl) {
    m_dataListeners.remove(dsl);
  }

  /**
   * Add an instance listener
   *
   * @param tsl an <code>InstanceListener</code> value
   */
  public synchronized void addInstanceListener(InstanceListener tsl) {
    m_instanceListeners.addElement(tsl);
  }

  /**
   * Remove an instance listener
   *
   * @param tsl an <code>InstanceListener</code> value
   */
  public synchronized void removeInstanceListener(InstanceListener tsl) {
    m_instanceListeners.removeElement(tsl);
  }
  
  /**
   * We don't have to keep track of configuration listeners (see the
   * documentation for ConfigurationListener/ConfigurationEvent).
   * 
   * @param cl a ConfigurationListener.
   */
  public synchronized void addConfigurationListener(ConfigurationListener cl) {
    
  }
  
  /**
   * We don't have to keep track of configuration listeners (see the
   * documentation for ConfigurationListener/ConfigurationEvent).
   * 
   * @param cl a ConfigurationListener.
   */
  public synchronized void removeConfigurationListener(ConfigurationListener cl) {
    
  }

  private void notifyDataOrTrainingListeners(EventObject ce) {
    Vector l;
    synchronized (this) {
      l = (ce instanceof TrainingSetEvent)
	? (Vector)m_trainingListeners.clone()
	: (Vector)m_dataListeners.clone();
    }
    if (l.size() > 0) {
      for(int i = 0; i < l.size(); i++) {
	if (ce instanceof TrainingSetEvent) {
	  ((TrainingSetListener)l.elementAt(i)).
	    acceptTrainingSet((TrainingSetEvent)ce);
	} else {
	  ((DataSourceListener)l.elementAt(i)).acceptDataSet((DataSetEvent)ce);
	}
      }
    }
  }

  private void notifyTestListeners(TestSetEvent ce) {
    Vector l;
    synchronized (this) {
      l = (Vector)m_testListeners.clone();
    }
    if (l.size() > 0) {
      for(int i = 0; i < l.size(); i++) {
	((TestSetListener)l.elementAt(i)).acceptTestSet(ce);
      }
    }
  }

  protected void notifyInstanceListeners(InstanceEvent tse) {
    Vector l;
    synchronized (this) {
      l = (Vector)m_instanceListeners.clone();
    }
    if (l.size() > 0) {
      for(int i = 0; i < l.size(); i++) {
	//	System.err.println("Notifying instance listeners "
	//			   +"(Filter)");
	((InstanceListener)l.elementAt(i)).acceptInstance(tse);
      }
    }
  }
  
  /**
   * Returns true if, at this time, 
   * the object will accept a connection with respect to the supplied
   * event name
   *
   * @param eventName the event
   * @return true if the object will accept a connection
   */
  public boolean connectionAllowed(String eventName) {

    if (m_listenees.containsKey(eventName)) {
      return false;
    }

    /* reject a test event if we don't have a training or data set event
    if (eventName.compareTo("testSet") == 0) {
      if (!m_listenees.containsKey("trainingSet") &&
	  !m_listenees.containsKey("dataSet")) {
	return false;
      }
      } */
    
    // will need to reject train/test listener if we have a
    // data source listener and vis versa
    if (m_listenees.containsKey("dataSet") &&
	(eventName.compareTo("trainingSet") == 0 ||
	 eventName.compareTo("testSet") == 0 ||
	eventName.compareTo("instance") == 0)) {
      return false;
    }

    if ((m_listenees.containsKey("trainingSet") ||
	 m_listenees.containsKey("testSet")) &&
	(eventName.compareTo("dataSet") == 0 || 
	eventName.compareTo("instance") == 0)) {
      return false;
    }

    if (m_listenees.containsKey("instance") &&
	(eventName.compareTo("trainingSet") == 0 ||
	 eventName.compareTo("testSet") == 0 ||
	 eventName.compareTo("dataSet") == 0)) {
      return false;
    }

    // reject an instance event connection if our filter isn't
    // streamable
    if (eventName.compareTo("instance") == 0 &&
	!(m_Filter instanceof StreamableFilter)) {
      return false;
    }
    return true;
  }

  /**
   * Returns true if, at this time, 
   * the object will accept a connection according to the supplied
   * EventSetDescriptor
   *
   * @param esd the EventSetDescriptor
   * @return true if the object will accept a connection
   */
  public boolean connectionAllowed(EventSetDescriptor esd) {
    return connectionAllowed(esd.getName());
  }

  /**
   * Notify this object that it has been registered as a listener with
   * a source with respect to the supplied event name
   *
   * @param eventName
   * @param source the source with which this object has been registered as
   * a listener
   */
  public synchronized void connectionNotification(String eventName,
						  Object source) {
    if (connectionAllowed(eventName)) {
      m_listenees.put(eventName, source);
      if (m_Filter instanceof ConnectionNotificationConsumer) {
        ((ConnectionNotificationConsumer) m_Filter).
          connectionNotification(eventName, source);
      }
    }
  }

  /**
   * Notify this object that it has been deregistered as a listener with
   * a source with respect to the supplied event name
   *
   * @param eventName the event
   * @param source the source with which this object has been registered as
   * a listener
   */
  public synchronized void disconnectionNotification(String eventName,
						     Object source) {
    if (m_Filter instanceof ConnectionNotificationConsumer) {
      ((ConnectionNotificationConsumer) m_Filter).
        disconnectionNotification(eventName, source);
    }
    m_listenees.remove(eventName);
  }

  /**
   * Function used to stop code that calls acceptTrainingSet, acceptTestSet,
   * or acceptDataSet. This is 
   * needed as filtering is performed inside a separate
   * thread of execution.
   *
   * @param tf a <code>boolean</code> value
   */
  private synchronized void block(boolean tf) {
    if (tf) {
      try {
	// only block if thread is still doing something useful!
	if (m_filterThread.isAlive() && m_state != IDLE) {
	  wait();
	}
      } catch (InterruptedException ex) {
      }
    } else {
      notifyAll();
    }
  }

  /**
   * Stop all action if possible
   */
  public void stop() {
    // tell all listenees (upstream beans) to stop
    Enumeration en = m_listenees.keys();
    while (en.hasMoreElements()) {
      Object tempO = m_listenees.get(en.nextElement());
      if (tempO instanceof BeanCommon) {
	((BeanCommon)tempO).stop();
      }
    }
    
    // stop the filter thread
    if (m_filterThread != null) {
      m_filterThread.interrupt();
      m_filterThread.stop();
      m_filterThread = null;
      m_visual.setStatic();
    }
  }
  
  /**
   * Returns true if. at this time, the bean is busy with some
   * (i.e. perhaps a worker thread is performing some calculation).
   * 
   * @return true if the bean is busy.
   */
  public boolean isBusy() {
    return (m_filterThread != null);
  }
  
  /**
   * Set a logger
   *
   * @param logger a <code>Logger</code> value
   */
  public void setLog(Logger logger) {
    m_log = logger;

    if (m_Filter != null && m_Filter instanceof LogWriter) {
      ((LogWriter) m_Filter).setLog(m_log);
    }
  }

  /**
   * Return an enumeration of user requests
   *
   * @return an <code>Enumeration</code> value
   */
  public Enumeration enumerateRequests() {
    Vector newVector = new Vector(0);
    if (m_filterThread != null) {
      newVector.addElement("Stop");
    }
    return newVector.elements();
  }

  /**
   * Perform the named request
   *
   * @param request a <code>String</code> value
   * @exception IllegalArgumentException if an error occurs
   */
  public void performRequest(String request) {
    if (request.compareTo("Stop") == 0) {
      stop();
    } else {
      throw new IllegalArgumentException(request
					 + " not supported (Filter)");
    }
  }

  /**
   * Returns true, if at the current time, the named event could
   * be generated. Assumes that supplied event names are names of
   * events that could be generated by this bean.
   *
   * @param eventName the name of the event in question
   * @return true if the named event could be generated at this point in
   * time
   */
  public boolean eventGeneratable(String eventName) {
    
    if (eventName.equals("configuration") && m_Filter != null) {
      return true;
    }
    
    // can't generate the named even if we are not receiving it as an
    // input!
    if (!m_listenees.containsKey(eventName)) {
      return false;
    }
    Object source = m_listenees.get(eventName);
    if (source instanceof EventConstraints) {
      if (!((EventConstraints)source).eventGeneratable(eventName)) {
	return false;
      }
    }
    if (eventName.compareTo("instance") == 0) {
      if (!(m_Filter instanceof StreamableFilter)) {
	return false;
      }
    }
    return true;
  }
  
  private String statusMessagePrefix() {
    return getCustomName() + "$" + hashCode() + "|"
    + ((m_Filter instanceof OptionHandler &&
        Utils.joinOptions(((OptionHandler)m_Filter).getOptions()).length() > 0) 
        ? Utils.joinOptions(((OptionHandler)m_Filter).getOptions()) + "|"
            : "");
  }
}
