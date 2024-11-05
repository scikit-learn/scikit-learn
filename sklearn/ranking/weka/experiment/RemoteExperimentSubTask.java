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
 *    RemoteExperimentSubTask.java
 *    Copyright (C) 2000 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.experiment;

import weka.core.RevisionHandler;
import weka.core.RevisionUtils;

import java.io.File;

/**
 * Class to encapsulate an experiment as a task that can be executed on
 * a remote host.
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @version $Revision: 1.10 $
 */
public class RemoteExperimentSubTask
  implements Task, RevisionHandler {

  /* Info on the task */
  private TaskStatusInfo m_result = new TaskStatusInfo();
  
  /* The (sub) experiment to execute */
  private Experiment m_experiment;
  
  public RemoteExperimentSubTask() {
    m_result.setStatusMessage("Not running.");
    m_result.setExecutionStatus(TaskStatusInfo.TO_BE_RUN);
  }

  /**
   * Set the experiment for this sub task
   * @param task the experiment
   */
  public void setExperiment(Experiment task) {
    m_experiment = task;
  }
  
  /**
   * Get the experiment for this sub task
   * @return this sub task's experiment
   */
  public Experiment getExperiment() {
    return m_experiment;
  }
  
  /**
   * Run the experiment
   */
  public void execute() {
    //      FastVector result = new FastVector();
    m_result = new TaskStatusInfo();
    m_result.setStatusMessage("Running...");
    String goodResult = "(sub)experiment completed successfully";
    String subTaskType;
    if (m_experiment.getRunLower() != m_experiment.getRunUpper()) {
      subTaskType = "(dataset "
	+ ((File)m_experiment.getDatasets().elementAt(0)).getName();
    } else {
      subTaskType = "(exp run # "+
	m_experiment.getRunLower();
    }
    try {	
      System.err.println("Initializing " + subTaskType + ")...");
      m_experiment.initialize();
      System.err.println("Iterating " + subTaskType + ")...");
      // Do not invoke runExperiment(): every exception will be lost
      while (m_experiment.hasMoreIterations()) {
        m_experiment.nextIteration();
      }
      System.err.println("Postprocessing " + subTaskType + ")...");
      m_experiment.postProcess();
    } catch (Exception ex) {
      ex.printStackTrace();
      String badResult =  "(sub)experiment " + subTaskType 
	+ ") failed : "+ex.toString();
      m_result.setExecutionStatus(TaskStatusInfo.FAILED);
      //	m_result.addElement(new Integer(RemoteExperiment.FAILED));
      //	m_result.addElement(badResult);
      m_result.setStatusMessage(badResult);
      m_result.setTaskResult("Failed");
      //      return m_result;
      return;
    }            
    //      m_result.addElement(new Integer(RemoteExperiment.FINISHED));
    //      m_result.addElement(goodResult);
    m_result.setExecutionStatus(TaskStatusInfo.FINISHED);
    m_result.setStatusMessage(goodResult+" "+subTaskType+").");
    m_result.setTaskResult("No errors");
    //    return m_result;
  }

  public TaskStatusInfo getTaskStatus() {
    return m_result;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 1.10 $");
  }
}






