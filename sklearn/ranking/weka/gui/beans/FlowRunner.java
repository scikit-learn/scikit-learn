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
 *    FlowRunner.java
 *    Copyright (C) 2008 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.Set;
import java.util.TreeMap;
import java.util.Vector;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStream;

import weka.core.Environment;
import weka.core.EnvironmentHandler;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.gui.Logger;
import weka.gui.beans.xml.*;

/**
 * Small utility class for executing KnowledgeFlow
 * flows outside of the KnowledgeFlow application
 *
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}org)
 * @version $Revision: 5928 $
 */
public class FlowRunner implements RevisionHandler {

  /** The potential flow(s) to execute */
  protected Vector m_beans;

  protected int m_runningCount = 0;

  protected transient Logger m_log = null;
  
  protected transient Environment m_env;
  
  /** run each Startable bean sequentially? (default in parallel) */
  protected boolean m_startSequentially = false;
  
  public static class SimpleLogger implements weka.gui.Logger {
    SimpleDateFormat m_DateFormat = 
      new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
    
    public void logMessage(String lm) {
      System.out.println(m_DateFormat.format(new Date()) + ": " + lm);
    }
    
    public void statusMessage(String lm) {
      System.out.println(m_DateFormat.format(new Date()) + ": " + lm);  
    }
  }

  /**
   * Constructor
   */
  public FlowRunner() {
    // make sure that properties and plugins are loaded
    KnowledgeFlowApp.loadProperties();
  }

  public void setLog(Logger log) {
    m_log = log;
  }
  
  protected void runSequentially(TreeMap<Integer, Startable> startables) {
    Set<Integer> s = startables.keySet();
    for (Integer i : s) {
      try {
        Startable startPoint = startables.get(i);
        startPoint.start();
        Thread.sleep(200);
        waitUntilFinished();
      } catch (Exception ex) {
        ex.printStackTrace();
        if (m_log != null) {
          m_log.logMessage(ex.getMessage());
          m_log.logMessage("Aborting...");
        } else {
          System.err.println(ex.getMessage());
          System.err.println("Aborting...");
        }
        break;
      }
    }
  }

  protected synchronized void launchThread(final Startable s, final int flowNum) {
    Thread t = new Thread() {
        private int m_num = flowNum;
        public void run() {
          try {
            s.start();
          } catch (Exception ex) {
            ex.printStackTrace();
            if (m_log != null) {
              m_log.logMessage(ex.getMessage());
            } else {
              System.err.println(ex.getMessage());
            }
          } finally {
            /*
            if (m_log != null) { 
              m_log.logMessage("[FlowRunner] flow " + m_num + " finished.");
            } else {
              System.out.println("[FlowRunner] Flow " + m_num + " finished.");
            }
            */
            decreaseCount();
          }
        }
      };
    m_runningCount++;
    t.setPriority(Thread.MIN_PRIORITY);
    t.start();
  }

  protected synchronized void decreaseCount() {
    m_runningCount--;
  }

  public synchronized void stopAllFlows() {
    for (int i = 0; i < m_beans.size(); i++) {
      BeanInstance temp = (BeanInstance)m_beans.elementAt(i);
      if (temp.getBean() instanceof BeanCommon) {
        // try to stop any execution
        ((BeanCommon)temp.getBean()).stop();
      }
    }
  }

  /**
   * Waits until all flows have finished executing before returning
   *
   */
  public void waitUntilFinished() {
    try {
      while (m_runningCount > 0) {
        Thread.sleep(200);
      }
      
      // now poll beans to see if there are any that are still busy
      // (i.e. any multi-threaded ones that queue data instead of blocking)
      while (true) {
        boolean busy = false;
        for (int i = 0; i < m_beans.size(); i++) {
          BeanInstance temp = (BeanInstance)m_beans.elementAt(i);
          if (temp.getBean() instanceof BeanCommon) {
            if (((BeanCommon)temp.getBean()).isBusy()) {
              busy = true;
              break; // for
            }            
          }
        }
        if (busy) {
          Thread.sleep(3000);
        } else {
          break; // while
        }
      }
    } catch (Exception ex) {
      if (m_log != null) {
        m_log.logMessage("[FlowRunner] Attempting to stop all flows...");
      } else {
        System.err.println("[FlowRunner] Attempting to stop all flows...");
      }
      stopAllFlows();
      //      ex.printStackTrace();
    }
  }

  /**
   * Load a serialized KnowledgeFlow (either binary or xml)
   *
   * @param fileName the name of the file to load from
   * @throws Exception if something goes wrong
   */
  public void load(String fileName) throws Exception {
    if (!fileName.endsWith(".kf") && !fileName.endsWith(".kfml")) {
      throw new Exception("Can only load and run binary or xml serialized KnowledgeFlows "
                          + "(*.kf | *.kfml)");
    }
    
    if (fileName.endsWith(".kf")) {
      loadBinary(fileName);
    } else if (fileName.endsWith(".kfml")) {
      loadXML(fileName);
    }
  }

  /**
   * Load a binary serialized KnowledgeFlow
   *
   * @param fileName the name of the file to load from
   * @throws Exception if something goes wrong
   */
  public void loadBinary(String fileName) throws Exception {
    if (!fileName.endsWith(".kf")) {
      throw new Exception("File must be a binary flow (*.kf)");
    }

    InputStream is = new FileInputStream(fileName);
    ObjectInputStream ois = new ObjectInputStream(is);
    m_beans = (Vector)ois.readObject();
    
    // don't need the graphical connections
    ois.close();
    
    if (m_env != null) {
      String parentDir = (new File(fileName)).getParent();
      if (parentDir == null) {
        parentDir = "./";
      }
      m_env.addVariable("Internal.knowledgeflow.directory", 
          parentDir);
    }
  }

  /**
   * Load an XML serialized KnowledgeFlow
   *
   * @param fileName the name of the file to load from
   * @throws Exception if something goes wrong
   */
  public void loadXML(String fileName) throws Exception {
    if (!fileName.endsWith(".kfml")) {
      throw new Exception("File must be an XML flow (*.kfml)");
    }

    XMLBeans xml = new XMLBeans(null, null);
    Vector v = (Vector) xml.read(new File(fileName));
    m_beans = (Vector) v.get(XMLBeans.INDEX_BEANINSTANCES);

    if (m_env != null) {
      String parentDir = (new File(fileName)).getParent();
      if (parentDir == null) {
        parentDir = "./";
      }
      m_env.addVariable("Internal.knowledgeflow.directory", 
          parentDir);
    } else {
      System.err.println("++++++++++++ Environment variables null!!...");
    }
  }
  
  /**
   * Get the vector holding the flow(s)
   *
   * @return the Vector holding the flow(s)
   */
  public Vector getFlows() {
    return m_beans;
  }

  /**
   * Set the vector holding the flows(s) to run
   *
   * @param beans the Vector holding the flows to run
   */
  public void setFlows(Vector beans) {
    m_beans = beans;
  }
  
  /**
   * Set the environment variables to use. NOTE: this needs
   * to be called BEFORE a load method is invoked to ensure
   * that the ${Internal.knowledgeflow.directory} variable get
   * set in the supplied Environment object.
   * 
   * @param env the environment variables to use.
   */
  public void setEnvironment(Environment env) {
    m_env = env;
  }
  
  /**
   * Get the environment variables that are in use.
   * 
   * @return the environment variables that are in ues.
   */
  public Environment getEnvironment() {
    return m_env;
  }
  
  /**
   * Set whether to launch Startable beans one after the other
   * or all in parallel.
   * 
   * @param s true if Startable beans are to be launched sequentially
   */
  public void setStartSequentially(boolean s) {
    m_startSequentially = s;
  }
  
  /**
   * Gets whether Startable beans will be launched sequentially
   * or all in parallel.
   * 
   * @return true if Startable beans will be launched sequentially
   */
  public boolean getStartSequentially() {
    return m_startSequentially;
  }

  /**
   * Launch all loaded KnowledgeFlow
   *
   * @throws Exception if something goes wrong during execution
   */
  public void run() throws Exception {
    if (m_beans == null) {
      throw new Exception("Don't seem to have any beans I can execute.");
    }
    
    // register the log (if set) with the beans
    for (int i = 0; i < m_beans.size(); i++) {
      BeanInstance tempB = (BeanInstance)m_beans.elementAt(i);
      if (m_log != null) {
        if (tempB.getBean() instanceof BeanCommon) {
          ((BeanCommon)tempB.getBean()).setLog(m_log);
        }
      }
        
      if (tempB.getBean() instanceof EnvironmentHandler) {
        ((EnvironmentHandler)tempB.getBean()).setEnvironment(m_env);
      }
    }
    
    int numFlows = 1;

    if (m_log != null) {
      if (m_startSequentially) {
        m_log.logMessage("[FlowRunner] launching flow start points sequentially...");
      } else {
        m_log.logMessage("[FlowRunner] launching flow start points in parallel...");
      }
    }
    TreeMap<Integer, Startable> startables = new TreeMap<Integer, Startable>();
    // look for a Startable bean...
    for (int i = 0; i < m_beans.size(); i++) {
      BeanInstance tempB = (BeanInstance)m_beans.elementAt(i);
      if (tempB.getBean() instanceof Startable) {
        Startable s = (Startable)tempB.getBean();
        // start that sucker (if it's happy to be started)...
        if (!m_startSequentially) {
          if (s.getStartMessage().charAt(0) != '$') {
            if (m_log != null) {
              m_log.logMessage("[FlowRunner] Launching flow "+numFlows+"...");
            } else {
              System.out.println("[FlowRunner] Launching flow "+numFlows+"...");
            }
            launchThread(s, numFlows);
            numFlows++;
          } else {
            String beanName = s.getClass().getName();
            if (s instanceof BeanCommon) {
              String customName = ((BeanCommon)s).getCustomName();
              beanName = customName;
            }
            if (m_log != null) {
              m_log.logMessage("[FlowRunner] WARNING: Can't start " + beanName + " at this time.");
            } else {
              System.out.println("[FlowRunner] WARNING: Can't start " + beanName + " at this time.");
            }
          }
        } else {
          boolean ok = false;
          Integer position = null;
          String beanName = s.getClass().getName();
          if (s instanceof BeanCommon) {
            String customName = ((BeanCommon)s).getCustomName();
            beanName = customName;
            // see if we have a parseable integer at the start of the name
            if (customName.indexOf(':') > 0) {
              String startPos = customName.substring(0, customName.indexOf(':'));
              try {
                position = new Integer(startPos);
                ok = true;
              } catch (NumberFormatException n) {
              }
            }            
          }
          
          if (!ok) {
            if (startables.size() == 0) {
              position = new Integer(0);
            } else {
              int newPos = startables.lastKey().intValue();
              newPos++;
              position = new Integer(newPos);
            }
          }
          
          if (s.getStartMessage().charAt(0) != '$') {
            if (m_log != null) {
              m_log.logMessage("[FlowRunner] adding start point " + beanName
                  + " to the execution list (position " + position + ")");
            } else {
              System.out.println("[FlowRunner] adding start point " + beanName
                  + " to the execution list (position " + position + ")");
            }
            startables.put(position, s);
          } else {
            if (m_log != null) {
              m_log.logMessage("[FlowRunner] WARNING: Can't start " + beanName + " at this time.");
            } else {
              System.out.println("[FlowRunner] WARNING: Can't start " + beanName + " at this time.");
            }
          }
        }
      }
    }
    
    if (m_startSequentially) {
      runSequentially(startables);
    }
  }

  /**
   * Main method for testing this class. <p>
   * <br>Usage:<br><br>
   * <pre>Usage:\n\nFlowRunner <serialized kf file></pre>
   *
   * @param args command line arguments
   */
  public static void main(String[] args) {
    weka.core.logging.Logger.log(weka.core.logging.Logger.Level.INFO, "Logging started");
    if (args.length < 1) {
      System.err.println("Usage:\n\nFlowRunner <serialized kf file> [-s]\n\n" 
          + "\tUse -s to launch start points sequentially (default launches "
          + "in parallel).");
    } else {
      try {
        FlowRunner fr = new FlowRunner();
        FlowRunner.SimpleLogger sl = new FlowRunner.SimpleLogger();
        String fileName = args[0];
        
        if (args.length == 2 && args[1].equals("-s")) {
          fr.setStartSequentially(true);
        }
        
        // start with the system-wide vars
        Environment env = Environment.getSystemWide();

        fr.setLog(sl);
        fr.setEnvironment(env);
        
        fr.load(fileName);
        fr.run();
        fr.waitUntilFinished();
        System.out.println("Finished all flows.");
        System.exit(1);
      } catch (Exception ex) {
        ex.printStackTrace();
        System.err.println(ex.getMessage());
      }
    }                         
  }

  public String getRevision() {
    return "$Revision: 5928 $";
  }
}
