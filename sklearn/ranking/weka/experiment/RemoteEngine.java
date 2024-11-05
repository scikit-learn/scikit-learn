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
 *    RemoteEngine.java
 *    Copyright (C) 2000 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.experiment;

import weka.core.Queue;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.core.Utils;

import java.net.InetAddress;
import java.net.URL;
import java.net.URLClassLoader;
import java.rmi.Naming;
import java.rmi.RMISecurityManager;
import java.rmi.RemoteException;
import java.rmi.server.UnicastRemoteObject;
import java.util.Enumeration;
import java.util.Hashtable;

/**
 * A general purpose server for executing Task objects sent via RMI.
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @version $Revision: 1.12 $
 */
public class RemoteEngine
  extends UnicastRemoteObject
  implements Compute, RevisionHandler {

  /** for serialization */
  private static final long serialVersionUID = -1021538162895448259L;

  /** The name of the host that this engine is started on */
  private String m_HostName = "local";

  /** A queue of waiting tasks */
  private Queue m_TaskQueue = new Queue();

  /** A queue of corresponding ID's for tasks */
  private Queue m_TaskIdQueue = new Queue();

  /** A hashtable of experiment status */
  private Hashtable m_TaskStatus = new Hashtable();

  /** Is there a task running */
  private boolean m_TaskRunning = false;
  
  /** Clean up interval (in ms) */
  protected static long CLEANUPTIMEOUT = 3600000;

  /**
   * Constructor
   * @param hostName name of the host
   * @exception RemoteException if something goes wrong
   */
  public RemoteEngine(String hostName) throws RemoteException {
    super();
    m_HostName = hostName;

    /* launch a clean-up thread. Will purge any failed or finished 
       tasks still in the TaskStatus hashtable after an hour */
       
    Thread cleanUpThread;
    cleanUpThread = new Thread() {
	public void run() {
	  while (true) {
	    try {
	      // sleep for a while
	      Thread.sleep(CLEANUPTIMEOUT);
	    } catch (InterruptedException ie) {}

	    if (m_TaskStatus.size() > 0) {
	      purge();
	    } else {
	      System.err.println("RemoteEngine : purge - no tasks to check.");
	    }
	  }
	}
      };
    cleanUpThread.setPriority(Thread.MIN_PRIORITY);
    cleanUpThread.setDaemon(true);
    cleanUpThread.start();
  }
  
  /**
   * Takes a task object and queues it for execution
   * @param t the Task object to execute
   * @return an identifier for the Task that can be used when querying
   * Task status
   */
  public synchronized Object executeTask(Task t) throws RemoteException {
    String taskId = ""+System.currentTimeMillis()+":";
    taskId += t.hashCode();
    addTaskToQueue(t, taskId);

    return taskId;
    //    return t.execute();
  }

  /**
   * Returns status information on a particular task
   *
   * @param taskId the ID of the task to check
   * @return a <code>TaskStatusInfo</code> encapsulating task status info
   * @exception Exception if an error occurs
   */
  public Object checkStatus(Object taskId) throws Exception {
    
    TaskStatusInfo inf = (TaskStatusInfo)m_TaskStatus.get(taskId);

    if (inf == null) {
      throw new Exception("RemoteEngine ("+m_HostName+") : Task not found.");
    }
    
    TaskStatusInfo result = new TaskStatusInfo();
    result.setExecutionStatus(inf.getExecutionStatus());
    result.setStatusMessage(inf.getStatusMessage());
    result.setTaskResult(inf.getTaskResult());

    if (inf.getExecutionStatus() == TaskStatusInfo.FINISHED ||
	inf.getExecutionStatus() == TaskStatusInfo.FAILED) {
      System.err.println("Finished/failed Task id : " 
			 + taskId + " checked by client. Removing.");
      inf.setTaskResult(null);
      inf = null;
      m_TaskStatus.remove(taskId);
    }
    inf = null;
    return result;
  }

  /**
   * Adds a new task to the queue.
   *
   * @param t a <code>Task</code> value to be added
   * @param taskId the id of the task to be added
   */
  private synchronized void addTaskToQueue(Task t, String taskId) {
    TaskStatusInfo newTask = t.getTaskStatus();
    if (newTask == null) {
      newTask = new TaskStatusInfo();
    }
    m_TaskQueue.push(t);
    m_TaskIdQueue.push(taskId);
    newTask.setStatusMessage("RemoteEngine ("
			     +m_HostName
			     +") : task " + taskId + " queued at postion: "
			     +m_TaskQueue.size());
    // add task status to HashTable
    m_TaskStatus.put(taskId, newTask);
    System.err.println("Task id : " + taskId + " Queued.");
    if (m_TaskRunning == false) {
      startTask();
    }
  }

  /**
   * Checks to see if there are any waiting tasks, and if no task is
   * currently running starts a waiting task.
   */
  private void startTask() {

    if (m_TaskRunning == false && m_TaskQueue.size() > 0) {
      Thread activeTaskThread;
      activeTaskThread = new Thread() {
	  public void run() {
	    m_TaskRunning = true;
	    Task currentTask = (Task)m_TaskQueue.pop();
	    String taskId = (String)m_TaskIdQueue.pop();
	    TaskStatusInfo tsi = (TaskStatusInfo)m_TaskStatus.get(taskId);
	    tsi.setExecutionStatus(TaskStatusInfo.PROCESSING);
	    tsi.setStatusMessage("RemoteEngine ("
				 +m_HostName
				 +") : task " + taskId + " running...");
	    try {
	      System.err.println("Launching task id : "
				 + taskId + "...");
	      currentTask.execute();
	      TaskStatusInfo runStatus = currentTask.getTaskStatus();
	      tsi.setExecutionStatus(runStatus.getExecutionStatus());
	      tsi.setStatusMessage("RemoteExperiment ("
				   +m_HostName+") "
				   +runStatus.getStatusMessage());
	      tsi.setTaskResult(runStatus.getTaskResult());
	    } catch (Error er) {
              // Object initialization can raise Error, which are not subclass of Exception
	      tsi.setExecutionStatus(TaskStatusInfo.FAILED);
              if (er.getCause() instanceof java.security.AccessControlException) {
                tsi.setStatusMessage("RemoteEngine ("
                                     +m_HostName
                                     +") : security error, check remote policy file.");
                System.err.println("Task id " + taskId + " Failed! Check remote policy file");
              }
              else {
                tsi.setStatusMessage("RemoteEngine ("
                                     +m_HostName
                                     +") : unknown initialization error.");
                System.err.println("Task id " + taskId + " Unknown initialization error");
              }
	    } catch (Exception ex) {
	      tsi.setExecutionStatus(TaskStatusInfo.FAILED);
              if (ex instanceof java.io.FileNotFoundException) {
                tsi.setStatusMessage("RemoteEngine ("
                                     +m_HostName
                                     +") : " + ex.getMessage());
                System.err.println("Task id " + taskId + " Failed, " + ex.getMessage());
              }
              else {
                tsi.setStatusMessage("RemoteEngine ("
                                     +m_HostName
                                     +") : task " + taskId + " failed.");
                System.err.println("Task id " + taskId + " Failed!");
              }
	    } finally {
	      if (m_TaskStatus.size() == 0) {
		purgeClasses();
	      }
	      m_TaskRunning = false;
	      // start any waiting tasks
	      startTask();
	    }
	  }
	};
      activeTaskThread.setPriority(Thread.MIN_PRIORITY);
      activeTaskThread.start();
    }
  }

  /**
   * Attempts to purge class types from the virtual machine. May take some
   * time as it relies on garbage collection
   */
  private void purgeClasses() {
    try {
      // see if we can purge classes
      ClassLoader prevCl = 
	Thread.currentThread().getContextClassLoader();
      ClassLoader urlCl = 
	URLClassLoader.newInstance(new URL[] {new URL("file:.")}, prevCl);
      Thread.currentThread().setContextClassLoader(urlCl);
    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }
  
  /**
   * Checks the hash table for failed/finished tasks. Any that have been
   * around for an @seeCLEANUPTIMEOUT or more are removed. Clients are expected to check
   * on the status of their remote tasks. Checking on the status of a
   * finished/failed task will remove it from the hash table, therefore
   * any failed/finished tasks left lying around for more than an hour
   * suggest that their client has died..
   *
   */
  private void purge() {
    Enumeration keys = m_TaskStatus.keys();
    long currentTime = System.currentTimeMillis();
    System.err.println("RemoteEngine purge. Current time : " + currentTime);
    while (keys.hasMoreElements()) {
      String taskId = (String)keys.nextElement();
      System.err.print("Examining task id : " + taskId + "... ");
      String timeString = taskId.substring(0, taskId.indexOf(':'));
      long ts = Long.valueOf(timeString).longValue();
      if (currentTime - ts > CLEANUPTIMEOUT) {
	TaskStatusInfo tsi = (TaskStatusInfo)m_TaskStatus.get(taskId);
	if ((tsi != null) 
	    && (tsi.getExecutionStatus() == TaskStatusInfo.FINISHED ||
                tsi.getExecutionStatus() == TaskStatusInfo.FAILED)) {
	  System.err.println("\nTask id : " 
			     + taskId + " has gone stale. Removing.");
	  m_TaskStatus.remove(taskId);
	  tsi.setTaskResult(null);
	  tsi = null;
	}
      } else {
	System.err.println("ok.");
      }
    }
    if (m_TaskStatus.size() == 0) {
      purgeClasses();
    }
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 1.12 $");
  }

  /**
   * Main method. Gets address of the local host, creates a remote engine
   * object and binds it in the RMI registry. If there is no RMI registry,
   * then it tries to create one with default port 1099.
   *
   * @param args 
   */
  public static void main(String[] args) {
    if (System.getSecurityManager() == null) {
      System.setSecurityManager(new RMISecurityManager());
    }
    
    int port = 1099;
    InetAddress localhost = null;
    try {
      localhost = InetAddress.getLocalHost();
      System.err.println("Host name : "+localhost.getHostName());
    } catch (Exception ex) {
      ex.printStackTrace();
    }
    String name;
    if (localhost != null) {
      name = localhost.getHostName();
    } else {
      name = "localhost";
    }
    
    // get optional port
    try {
      String portOption = Utils.getOption("p", args);
      if (!portOption.equals("")) 
        port = Integer.parseInt(portOption);
    } catch (Exception ex) {
      System.err.println("Usage : -p <port>");
    }

    if (port != 1099) {
      name = name + ":" + port;
    }
    name = "//"+name+"/RemoteEngine";
    
    try {
      Compute engine = new RemoteEngine(name);
      
      try {      
        Naming.rebind(name, engine);
        System.out.println("RemoteEngine bound in RMI registry");
      } catch (RemoteException ex) {
        // try to bootstrap a new registry
        System.err.println("Attempting to start RMI registry on port " + port + "...");
        java.rmi.registry.LocateRegistry.createRegistry(port);
        Naming.bind(name, engine);
        System.out.println("RemoteEngine bound in RMI registry");
      }
      
    } catch (Exception e) {
      System.err.println("RemoteEngine exception: " + 
			 e.getMessage());
      e.printStackTrace();
    }
  }
}
