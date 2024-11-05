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
 *    InstanceLoader.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.streams;

import weka.core.Instance;
import weka.core.Instances;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.Reader;
import java.util.Vector;

import javax.swing.JButton;
import javax.swing.JPanel;
import javax.swing.JTextField;

/** 
 * A bean that produces a stream of instances from a file.
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 1.5 $
 */
public class InstanceLoader
  extends JPanel 
  implements ActionListener, InstanceProducer {

  /** for serialization */
  private static final long serialVersionUID = -8725567310271862492L;
  
  private Vector m_Listeners;
  private Thread m_LoaderThread;
  private Instance m_OutputInstance;
  private Instances m_OutputInstances;
  private boolean m_Debug;
  private JButton m_StartBut;
  private JTextField m_FileNameTex;

  private class LoadThread extends Thread {
    
    private InstanceProducer m_IP;
    public LoadThread(InstanceProducer ip) {
      
      m_IP = ip;
    }

    public void run() {
      
      try {
	m_StartBut.setText("Stop");
	m_StartBut.setBackground(Color.red);
	if (m_Debug) {
	  System.err.println("InstanceLoader::LoadThread::run()");
	}
	// Load the instances one at a time and pass them on to the listener
	Reader input = new BufferedReader(
		       new FileReader(m_FileNameTex.getText()));
	m_OutputInstances = new Instances(input, 1);
	if (m_Debug) {
	  System.err.println("InstanceLoader::LoadThread::run()"
			     + " - Instances opened from: "
			     + m_FileNameTex.getText());
	}
	InstanceEvent ie = new InstanceEvent(m_IP,
					     InstanceEvent.FORMAT_AVAILABLE);
	notifyInstanceProduced(ie);
	while (m_OutputInstances.readInstance(input)) {
	  if (m_LoaderThread != this) {
	    return;
	  }
	  if (m_Debug) {
	    System.err.println("InstanceLoader::LoadThread::run()"
			       + " - read instance");
	  }
	  // put the instance into a queue?
	  m_OutputInstance = m_OutputInstances.instance(0);
	  m_OutputInstances.delete(0);
	  ie = new InstanceEvent(m_IP, InstanceEvent.INSTANCE_AVAILABLE);
	  notifyInstanceProduced(ie);
	}
	ie = new InstanceEvent(m_IP, InstanceEvent.BATCH_FINISHED);
	notifyInstanceProduced(ie);
      } catch (Exception ex) {
	System.err.println(ex.getMessage());
      } finally {
	m_LoaderThread = null;
	m_StartBut.setText("Start");
	m_StartBut.setBackground(Color.green);
      }
    }
  }

  public InstanceLoader() {
    setLayout(new BorderLayout());
    m_StartBut = new JButton("Start");
    m_StartBut.setBackground(Color.green);
    add("West",m_StartBut);
    m_StartBut.addActionListener(this);
    m_FileNameTex = new JTextField("/home/trigg/datasets/UCI/iris.arff");
    add("Center",m_FileNameTex);
    m_Listeners = new Vector();
    //    setSize(60,40);
  }

  public void setDebug(boolean debug) {
    
    m_Debug = debug;
  }
  
  public boolean getDebug() {
    
    return m_Debug;
  }

  public void setArffFile(String newArffFile) {
    
    m_FileNameTex.setText(newArffFile);
  }
  
  public String getArffFile() {
    return m_FileNameTex.getText();
  }

  public synchronized void addInstanceListener(InstanceListener ipl) {
    
    m_Listeners.addElement(ipl);
  }
  
  public synchronized void removeInstanceListener(InstanceListener ipl) {
    
    m_Listeners.removeElement(ipl);
  }
  
  protected void notifyInstanceProduced(InstanceEvent e) {
    
    if (m_Debug) {
      System.err.println("InstanceLoader::notifyInstanceProduced()");
    }
    Vector l;
    synchronized (this) {
      l = (Vector)m_Listeners.clone();
    }
    if (l.size() > 0) {
      for(int i = 0; i < l.size(); i++) {
	((InstanceListener)l.elementAt(i)).instanceProduced(e);
      }
      if (e.getID() == InstanceEvent.INSTANCE_AVAILABLE) {
	m_OutputInstance = null;
      }
    }
  }

  public Instances outputFormat() throws Exception {
    
    if (m_OutputInstances == null) {
      throw new Exception("No output format defined.");
    }
    return new Instances(m_OutputInstances,0);
  }
  
  public Instance outputPeek() throws Exception {
    
    if ((m_OutputInstances == null)
	|| (m_OutputInstance == null)) {
      return null;
    }
    return (Instance)m_OutputInstance.copy();
  }

  public void actionPerformed(ActionEvent e) {
    
    Object source = e.getSource();

    if (source == m_StartBut) {
      // load the arff file and send the instances out to the listener
      if (m_LoaderThread == null) {
	m_LoaderThread = new LoadThread(this);
	m_LoaderThread.setPriority(Thread.MIN_PRIORITY);
	m_LoaderThread.start();
      } else {
	m_LoaderThread = null;
      }
    }
  }
}
