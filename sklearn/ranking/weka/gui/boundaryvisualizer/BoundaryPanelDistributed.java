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
 *   BoundaryPanelDistrubuted.java
 *   Copyright (C) 2003 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.boundaryvisualizer;

import weka.classifiers.Classifier;
import weka.classifiers.AbstractClassifier;
import weka.core.FastVector;
import weka.core.Instances;
import weka.core.Utils;
import weka.experiment.Compute;
import weka.experiment.RemoteExperimentEvent;
import weka.experiment.RemoteExperimentListener;
import weka.experiment.TaskStatusInfo;

import java.awt.BorderLayout;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.ObjectInputStream;
import java.rmi.Naming;
import java.util.Vector;

/**
 * This class extends BoundaryPanel with code for distributing the
 * processing necessary to create a visualization among a list of
 * remote machines. Specifically, a visualization is broken down and
 * processed row by row using the available remote computers.
 *
 * @author <a href="mailto:mhall@cs.waikato.ac.nz">Mark Hall</a>
 * @version $Revision: 5928 $
 * @since 1.0
 * @see BoundaryPanel
 */
public class BoundaryPanelDistributed
  extends BoundaryPanel {

  /** for serialization */
  private static final long serialVersionUID = -1743284397893937776L;

  /** a list of RemoteExperimentListeners */
  protected Vector m_listeners = new Vector();

  /** Holds the names of machines with remoteEngine servers running */
  protected Vector m_remoteHosts = new Vector();
  
  /** The queue of available hosts */
  private weka.core.Queue m_remoteHostsQueue = new weka.core.Queue();

  /** The status of each of the remote hosts */
  private int [] m_remoteHostsStatus;

  /** The number of times tasks have failed on each remote host */
  private int [] m_remoteHostFailureCounts;

  protected static final int AVAILABLE=0;
  protected static final int IN_USE=1;
  protected static final int CONNECTION_FAILED=2;
  protected static final int SOME_OTHER_FAILURE=3;

  protected static final int MAX_FAILURES=3;

  /** Set to true if MAX_FAILURES exceeded on all hosts or connections fail 
      on all hosts or user aborts plotting */
  private boolean m_plottingAborted = false;

  /** The number of hosts removed due to exceeding max failures */
  private int m_removedHosts;

  /** The count of failed sub-tasks */
  private int m_failedCount;

  /** The count of successfully completed sub-tasks */
  private int m_finishedCount;

  /** The queue of sub-tasks waiting to be processed */
  private weka.core.Queue m_subExpQueue = new weka.core.Queue();

  /** number of seconds between polling server */
  private int m_minTaskPollTime = 1000;

  private int [] m_hostPollingTime;

  /**
   * Creates a new <code>BoundaryPanelDistributed</code> instance.
   *
   * @param panelWidth width of the display
   * @param panelHeight height of the display
   */
  public BoundaryPanelDistributed(int panelWidth, int panelHeight) {
    super(panelWidth, panelHeight);
  }

  /**
   * Set a list of host names of machines to distribute processing to
   *
   * @param remHosts a Vector of host names (Strings)
   */
  public void setRemoteHosts(Vector remHosts) {
    m_remoteHosts = remHosts;
  }

  /**
   * Add an object to the list of those interested in recieving update
   * information from the RemoteExperiment
   * @param r a listener
   */
  public void addRemoteExperimentListener(RemoteExperimentListener r) {
    m_listeners.addElement(r);
  }

  protected void initialize() {
    super.initialize();

    m_plottingAborted = false;
    m_finishedCount = 0;
    m_failedCount = 0;

    // initialize all remote hosts to available
    m_remoteHostsStatus = new int [m_remoteHosts.size()];    
    m_remoteHostFailureCounts = new int [m_remoteHosts.size()];

    m_remoteHostsQueue = new weka.core.Queue();

    if (m_remoteHosts.size() == 0) {
      System.err.println("No hosts specified!");
      System.exit(1);
    }

    // prime the hosts queue
    m_hostPollingTime = new int [m_remoteHosts.size()];
    for (int i=0;i<m_remoteHosts.size();i++) {
      m_remoteHostsQueue.push(new Integer(i));
      m_hostPollingTime[i] = m_minTaskPollTime;
    }

    // set up sub taskss (just holds the row numbers to be processed
    m_subExpQueue = new weka.core.Queue();
    for (int i = 0; i < m_panelHeight; i++) {
      m_subExpQueue.push(new Integer(i));
    }
    
    try {
      // need to build classifier and data generator
      m_classifier.buildClassifier(m_trainingData);
    } catch (Exception ex) {
      ex.printStackTrace();
      System.exit(1);
    }
    
    boolean [] attsToWeightOn;
    // build DataGenerator
    attsToWeightOn = new boolean[m_trainingData.numAttributes()];
    attsToWeightOn[m_xAttribute] = true;
    attsToWeightOn[m_yAttribute] = true;
    
    m_dataGenerator.setWeightingDimensions(attsToWeightOn);    
    try {
      m_dataGenerator.buildGenerator(m_trainingData);
    } catch (Exception ex) {
      ex.printStackTrace();
      System.exit(1);
    }
  }

  /**
   * Start processing
   *
   * @exception Exception if an error occurs
   */
  public void start() throws Exception {
    // done in the sub task
    /*     m_numOfSamplesPerGenerator = 
	   (int)Math.pow(m_samplesBase, m_trainingData.numAttributes()-3); */

    m_stopReplotting = true;
    if (m_trainingData == null) {
      throw new Exception("No training data set (BoundaryPanel)");
    }
    if (m_classifier == null) {
      throw new Exception("No classifier set (BoundaryPanel)");
    }
    if (m_dataGenerator == null) {
      throw new Exception("No data generator set (BoundaryPanel)");
    }
    if (m_trainingData.attribute(m_xAttribute).isNominal() || m_trainingData.attribute(m_xAttribute).isRanking() || 
	m_trainingData.attribute(m_yAttribute).isNominal() || m_trainingData.attribute(m_yAttribute).isRanking()) {
      throw new Exception("Visualization dimensions must be numeric "
			  +"(BoundaryPanel)");
    }
    
    computeMinMaxAtts();
    initialize();

    // launch tasks on all available hosts
    int totalHosts = m_remoteHostsQueue.size();
    for (int i = 0; i < totalHosts; i++) {
      availableHost(-1);
      Thread.sleep(70);
    }
  }

  /**
   * Push a host back onto the list of available hosts and launch a waiting
   * Task (if any).
   *
   * @param hostNum the number of the host to return to the queue. -1
   * if no host to return.
   */
  protected synchronized void availableHost(int hostNum) {
    if (hostNum >= 0) { 
      if (m_remoteHostFailureCounts[hostNum] < MAX_FAILURES) {
	m_remoteHostsQueue.push(new Integer(hostNum));
      } else {
	notifyListeners(false,true,false,"Max failures exceeded for host "
			+((String)m_remoteHosts.elementAt(hostNum))
			+". Removed from host list.");
	m_removedHosts++;
      }
    }

    // check for all sub exp complete or all hosts failed or failed count
    // exceeded
    if (m_failedCount == (MAX_FAILURES * m_remoteHosts.size())) {
      m_plottingAborted = true;
      notifyListeners(false,true,true,"Plotting aborted! Max failures "
		      +"exceeded on all remote hosts.");
      return;
    }

    /*    System.err.println("--------------");
    System.err.println("exp q :"+m_subExpQueue.size());
    System.err.println("host list size "+m_remoteHosts.size());
    System.err.println("actual host list size "+m_remoteHostsQueue.size());
    System.err.println("removed hosts "+m_removedHosts); */
    if (m_subExpQueue.size() == 0 && 
	(m_remoteHosts.size() == 
	 (m_remoteHostsQueue.size() + m_removedHosts))) {
      if (m_plotTrainingData) {
	plotTrainingData();
      }
      notifyListeners(false,true,true,"Plotting completed successfully.");

      return;
    }


    if (checkForAllFailedHosts()) {
      return;
    }

    if (m_plottingAborted && 
	(m_remoteHostsQueue.size() + m_removedHosts) == 
	m_remoteHosts.size()) {
      notifyListeners(false,true,true,"Plotting aborted. All remote tasks "
		      +"finished.");
    }

    if (!m_subExpQueue.empty() && !m_plottingAborted) {
      if (!m_remoteHostsQueue.empty()) {
	int availHost, waitingTask;
	try {
	  availHost = ((Integer)m_remoteHostsQueue.pop()).intValue();
	  waitingTask = ((Integer)m_subExpQueue.pop()).intValue();
	  launchNext(waitingTask, availHost);
	} catch (Exception ex) {
	  ex.printStackTrace();
	}
      }
    }    
  }

  /**
   * Inform all listeners of progress
   * @param status true if this is a status type of message
   * @param log true if this is a log type of message
   * @param finished true if the remote task has finished
   * @param message the message.
   */
  private synchronized void notifyListeners(boolean status, 
					    boolean log, 
					    boolean finished,
					    String message) {
    if (m_listeners.size() > 0) {
      for (int i=0;i<m_listeners.size();i++) {
	RemoteExperimentListener r = 
	  (RemoteExperimentListener)(m_listeners.elementAt(i));
	r.remoteExperimentStatus(new RemoteExperimentEvent(status,
							   log,
							   finished,
							   message));
      }
    } else {
      System.err.println(message);
    }
  }

  /**
   * Check to see if we have failed to connect to all hosts
   */
  private boolean checkForAllFailedHosts() {
    boolean allbad = true;
    for (int i = 0; i < m_remoteHostsStatus.length; i++) {
      if (m_remoteHostsStatus[i] != CONNECTION_FAILED) {
	allbad = false;
	break;
      }
    }
    if (allbad) {
      m_plottingAborted = true;
      notifyListeners(false,true,true,"Plotting aborted! All connections "
		      +"to remote hosts failed.");
    }
    return allbad;
  }

  /**
   * Increment the number of successfully completed sub experiments
   */
  protected synchronized void incrementFinished() {
    m_finishedCount++;
  }

  /**
   * Increment the overall number of failures and the number of failures for
   * a particular host
   * @param hostNum the index of the host to increment failure count
   */
  protected synchronized void incrementFailed(int hostNum) {
    m_failedCount++;
    m_remoteHostFailureCounts[hostNum]++;
  }

  /**
   * Push an experiment back on the queue of waiting experiments
   * @param expNum the index of the experiment to push onto the queue
   */
  protected synchronized void waitingTask(int expNum) {
    m_subExpQueue.push(new Integer(expNum));
  }

  protected void launchNext(final int wtask, final int ah) {
    Thread subTaskThread;
    subTaskThread = new Thread() {
	public void run() {
	  m_remoteHostsStatus[ah] = IN_USE;
	  //	  m_subExpComplete[wtask] = TaskStatusInfo.PROCESSING;
	  RemoteBoundaryVisualizerSubTask vSubTask = 
	    new RemoteBoundaryVisualizerSubTask();
	  vSubTask.setXAttribute(m_xAttribute);
	  vSubTask.setYAttribute(m_yAttribute);
	  vSubTask.setRowNumber(wtask);
	  vSubTask.setPanelWidth(m_panelWidth);
	  vSubTask.setPanelHeight(m_panelHeight);
	  vSubTask.setPixHeight(m_pixHeight);
	  vSubTask.setPixWidth(m_pixWidth);
	  vSubTask.setClassifier(m_classifier);
	  vSubTask.setDataGenerator(m_dataGenerator);
	  vSubTask.setInstances(m_trainingData);
	  vSubTask.setMinMaxX(m_minX, m_maxX);
	  vSubTask.setMinMaxY(m_minY, m_maxY);
	  vSubTask.setNumSamplesPerRegion(m_numOfSamplesPerRegion);
	  vSubTask.setGeneratorSamplesBase(m_samplesBase);
	  try {
	    String name = "//"
	      +((String)m_remoteHosts.elementAt(ah))
	      +"/RemoteEngine";
	    Compute comp = (Compute) Naming.lookup(name);
	    // assess the status of the sub-exp
	    notifyListeners(false,true,false,"Starting row "
			    +wtask
			    +" on host "
			    +((String)m_remoteHosts.elementAt(ah)));
	    Object subTaskId = comp.executeTask(vSubTask);
	    boolean finished = false;
	    TaskStatusInfo is = null;
	    long startTime = System.currentTimeMillis();
	    while (!finished) {
	      try {
		Thread.sleep(Math.max(m_minTaskPollTime, 
				      m_hostPollingTime[ah]));
		
		TaskStatusInfo cs = (TaskStatusInfo)comp.
		  checkStatus(subTaskId);
		if (cs.getExecutionStatus() == TaskStatusInfo.FINISHED) {
		  // push host back onto queue and try launching any waiting 
		  // sub-experiments
		  long runTime = System.currentTimeMillis() - startTime;
		  runTime /= 4;
		  if (runTime < 1000) {
		    runTime = 1000;
		  }
		  m_hostPollingTime[ah] = (int)runTime;

		  // Extract the row from the result
		  RemoteResult rr =  (RemoteResult)cs.getTaskResult();
		  double [][] probs = rr.getProbabilities();
		  
		  for (int i = 0; i < m_panelWidth; i++) {
		    m_probabilityCache[wtask][i] = probs[i];
		    if (i < m_panelWidth-1) {
		      plotPoint(i, wtask, probs[i], false);
		    } else {
		      plotPoint(i, wtask, probs[i], true);
		    }
		  }
		  notifyListeners(false, true, false,  cs.getStatusMessage());
		  m_remoteHostsStatus[ah] = AVAILABLE;
		  incrementFinished();
		  availableHost(ah);
		  finished = true;
		} else if (cs.getExecutionStatus() == 
			   TaskStatusInfo.FAILED) {
		  // a non connection related error---possibly host doesn't have
		  // access to data sets or security policy is not set up
		  // correctly or classifier(s) failed for some reason
		  notifyListeners(false, true, false,  
				  cs.getStatusMessage());
		  m_remoteHostsStatus[ah] = SOME_OTHER_FAILURE;
		  //		  m_subExpComplete[wexp] = TaskStatusInfo.FAILED;
		  notifyListeners(false,true,false,"Row "+wtask
				  +" "+cs.getStatusMessage()
				  +". Scheduling for execution on another host.");
		  incrementFailed(ah);
		  // push experiment back onto queue
		  waitingTask(wtask);	
		  // push host back onto queue and try launching any waiting 
		  // Tasks. Host is pushed back on the queue as the
		  // failure may be temporary.
		  availableHost(ah);
		  finished = true;
		} else {
		  if (is == null) {
		    is = cs;
		    notifyListeners(false, true, false, cs.getStatusMessage());
		  } else {
		    RemoteResult rr = (RemoteResult)cs.getTaskResult();
		    if (rr != null) {
		      int percentComplete = rr.getPercentCompleted();
		      String timeRemaining = "";
		      if (percentComplete > 0 && percentComplete < 100) {
			double timeSoFar = (double)System.currentTimeMillis() -
			  (double)startTime;
			double timeToGo = 
			  ((100.0 - percentComplete) 
			   / (double)percentComplete) * timeSoFar;
			if (timeToGo < m_hostPollingTime[ah]) {
			  m_hostPollingTime[ah] = (int)timeToGo;
			}
			String units = "seconds";
			timeToGo /= 1000.0;
			if (timeToGo > 60) {
			  units = "minutes";
			  timeToGo /= 60.0;
			}
			if (timeToGo > 60) {
			  units = "hours";
			  timeToGo /= 60.0;
			}
			timeRemaining = " (approx. time remaining "
			  +Utils.doubleToString(timeToGo, 1)+" "+units+")";
		      }
		      if (percentComplete < 25 
			  /*&& minTaskPollTime < 30000*/) {		
			if (percentComplete > 0) {
			  m_hostPollingTime[ah] = 
			    (int)((25.0 / (double)percentComplete) * 
				  m_hostPollingTime[ah]);
			} else {
			  m_hostPollingTime[ah] *= 2;
			}
			if (m_hostPollingTime[ah] > 60000) {
			  m_hostPollingTime[ah] = 60000;
			}
		      }
		      notifyListeners(false, true, false,
				      "Row "+wtask+" "+percentComplete
				      +"% complete"+timeRemaining+".");
		    } else {
		      notifyListeners(false, true, false,
				      "Row "+wtask+" queued on "
				      +((String)m_remoteHosts.
					elementAt(ah)));
		      if (m_hostPollingTime[ah] < 60000) {
			m_hostPollingTime[ah] *= 2;
		      }
		    }

		    is = cs;
		  }
		}
	      } catch (InterruptedException ie) {
		ie.printStackTrace();
	      }
	    }
	  } catch (Exception ce) {
	    m_remoteHostsStatus[ah] = CONNECTION_FAILED;
	    m_removedHosts++;
	    System.err.println(ce);
	    ce.printStackTrace();
	    notifyListeners(false,true,false,"Connection to "
			    +((String)m_remoteHosts.elementAt(ah))
			    +" failed. Scheduling row "
			    +wtask
			    +" for execution on another host.");
	    checkForAllFailedHosts();
	    waitingTask(wtask);
	  } finally {
	    if (isInterrupted()) {
	      System.err.println("Sub exp Interupted!");
	    }
	  }
	}
      };
    subTaskThread.setPriority(Thread.MIN_PRIORITY);
    subTaskThread.start();
  }

  /**
   * Main method for testing this class
   *
   * @param args a <code>String[]</code> value
   */
  public static void main (String [] args) {
    try {
      if (args.length < 8) {
	System.err.println("Usage : BoundaryPanelDistributed <dataset> "
			   +"<class col> <xAtt> <yAtt> "
			   +"<base> <# loc/pixel> <kernel bandwidth> "
			   +"<display width> "
			   +"<display height> <classifier "
			   +"[classifier options]>");
	System.exit(1);
      }
      
      Vector hostNames = new Vector();
      // try loading hosts file
      try {
	BufferedReader br = new BufferedReader(new FileReader("hosts.vis"));
	String hostName = br.readLine();
	while (hostName != null) {
	  System.out.println("Adding host "+hostName);
	  hostNames.add(hostName);
	  hostName = br.readLine();
	}
	br.close();
      } catch (Exception ex) {
	System.err.println("No hosts.vis file - create this file in "
			   +"the current directory with one host name "
			   +"per line, or use BoundaryPanel instead.");
	System.exit(1);
      }

      final javax.swing.JFrame jf = 
	new javax.swing.JFrame("Weka classification boundary visualizer");
      jf.getContentPane().setLayout(new BorderLayout());

      System.err.println("Loading instances from : "+args[0]);
      java.io.Reader r = new java.io.BufferedReader(
			 new java.io.FileReader(args[0]));
      final Instances i = new Instances(r);
      i.setClassIndex(Integer.parseInt(args[1]));

      //      bv.setClassifier(new Logistic());
      final int xatt = Integer.parseInt(args[2]);
      final int yatt = Integer.parseInt(args[3]);
      int base = Integer.parseInt(args[4]);
      int loc = Integer.parseInt(args[5]);

      int bandWidth = Integer.parseInt(args[6]);
      int panelWidth = Integer.parseInt(args[7]);
      int panelHeight = Integer.parseInt(args[8]);

      final String classifierName = args[9];
      final BoundaryPanelDistributed bv = 
	new BoundaryPanelDistributed(panelWidth,panelHeight);
      bv.addRemoteExperimentListener(new RemoteExperimentListener() {
	  public void remoteExperimentStatus(RemoteExperimentEvent e) {
	    if (e.m_experimentFinished) {
	      String classifierNameNew = 
		classifierName.substring(classifierName.lastIndexOf('.')+1, 
					 classifierName.length());
	      bv.saveImage(classifierNameNew+"_"+i.relationName()
			   +"_X"+xatt+"_Y"+yatt+".jpg");
	    } else {
	      System.err.println(e.m_messageString);
	    }
	  }
	});
      bv.setRemoteHosts(hostNames);

      jf.getContentPane().add(bv, BorderLayout.CENTER);
      jf.setSize(bv.getMinimumSize());
      //      jf.setSize(200,200);
      jf.addWindowListener(new java.awt.event.WindowAdapter() {
	  public void windowClosing(java.awt.event.WindowEvent e) {
	    jf.dispose();
	    System.exit(0);
	  }
	});

      jf.pack();
      jf.setVisible(true);
      //      bv.initialize();
      bv.repaint();
      

      String [] argsR = null;
      if (args.length > 10) {
	argsR = new String [args.length-10];
	for (int j = 10; j < args.length; j++) {
	  argsR[j-10] = args[j];
	}
      }
      Classifier c = AbstractClassifier.forName(args[9], argsR);
      KDDataGenerator dataGen = new KDDataGenerator();
      dataGen.setKernelBandwidth(bandWidth);
      bv.setDataGenerator(dataGen);
      bv.setNumSamplesPerRegion(loc);
      bv.setGeneratorSamplesBase(base);
      bv.setClassifier(c);
      bv.setTrainingData(i);
      bv.setXAttribute(xatt);
      bv.setYAttribute(yatt);

      try {
	// try and load a color map if one exists
	FileInputStream fis = new FileInputStream("colors.ser");
	ObjectInputStream ois = new ObjectInputStream(fis);
	FastVector colors = (FastVector)ois.readObject();
	bv.setColors(colors);	
      } catch (Exception ex) {
	System.err.println("No color map file");
      }
      bv.start();
    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }
}
