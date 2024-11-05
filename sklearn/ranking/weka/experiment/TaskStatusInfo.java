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
 *    TaskStatusInfo.java
 *    Copyright (C) 2001 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.experiment;

import weka.core.RevisionHandler;
import weka.core.RevisionUtils;

import java.io.Serializable;

/**
 * A class holding information for tasks being executed
 * on RemoteEngines. Also holds an object encapsulating any returnable result
 * produced by the task (Note: result object must be serializable). Task
 * objects execute methods return instances of this class. RemoteEngines also
 * use this class for storing progress information for tasks that they
 * execute.
 *
 * @author <a href="mailto:mhall@cs.waikato.ac.nz">Mark Hall</a>
 * @version $Revision: 1.5 $
 */
public class TaskStatusInfo
  implements Serializable, RevisionHandler {

  /** for serialization */
  private static final long serialVersionUID = -6129343303703560015L;
  
  public static final int TO_BE_RUN = 0;
  public static final int PROCESSING=1;
  public static final int FAILED=2;
  public static final int FINISHED=3;

  /**
   * Holds current execution status.
   */
  private int m_ExecutionStatus = TO_BE_RUN;

  /**
   * Holds current status message.
   */
  private String m_StatusMessage = "New Task";

  /**
   * Holds task result. Set to null for no returnable result.
   */
  private Object m_TaskResult = null;

  /**
   * Set the execution status of this Task.
   *
   * @param newStatus the new execution status code
   */
  public void setExecutionStatus(int newStatus) {
    m_ExecutionStatus = newStatus;
  }

  /**
   * Get the execution status of this Task.
   * @return the execution status
   */
  public int getExecutionStatus() {
    return m_ExecutionStatus;
  }

  /**
   * Set the status message.
   *
   * @param newMessage the new status message
   */
  public void setStatusMessage(String newMessage) {
    m_StatusMessage = newMessage;
  }

  /**
   * Get the status message.
   *
   * @return the status message
   */
  public String getStatusMessage() {
    return m_StatusMessage;
  }

  /**
   * Set the returnable result for this task..
   *
   * @param taskResult the new returnable result for the task. null if no
   * result is returnable.
   */
  public void setTaskResult(Object taskResult) {
    m_TaskResult = taskResult;
  }

  /**
   * Get the returnable result of this task.
   *
   * @return an object encapsulating the result of executing the task. May
   * be null if the task has no returnable result (eg. a remote experiment
   * task that sends its results to a data base).
   */
  public Object getTaskResult() {
    return m_TaskResult;
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
