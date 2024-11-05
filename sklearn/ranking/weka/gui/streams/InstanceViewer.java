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
 *    InstanceViewer.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.streams;

import weka.core.Instance;
import weka.core.Instances;

import java.awt.BorderLayout;

import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;

/**
 * This is a very simple instance viewer - just displays the dataset as
 * text output as it would be written to a file. A more complex viewer
 * might be more spreadsheet-like
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 1.5 $
 */
public class InstanceViewer
  extends JPanel
  implements InstanceListener {

  /** for serialization */
  private static final long serialVersionUID = -4925729441294121772L;
  
  private JTextArea m_OutputTex;
  private boolean m_Debug;
  private boolean m_Clear;
  private String m_UpdateString;

  private void updateOutput() {
    
    m_OutputTex.append(m_UpdateString);
    m_UpdateString = "";
  }

  private void clearOutput() {
    
    m_UpdateString = "";
    m_OutputTex.setText("");
  }

  public void inputFormat(Instances instanceInfo) {
    
    if (m_Debug) {
      System.err.println("InstanceViewer::inputFormat()\n"
			 + instanceInfo.toString());
    }
    if (m_Clear) {
      clearOutput();
    }
    m_UpdateString += instanceInfo.toString();
    updateOutput();
  }

  public void input(Instance instance) throws Exception {
    
    if (m_Debug) {
      System.err.println("InstanceViewer::input(" + instance +")");
    }
    m_UpdateString += instance.toString() + "\n";
    updateOutput();
  }
  
  public void batchFinished() {
    
    updateOutput();
    if (m_Debug) {
      System.err.println("InstanceViewer::batchFinished()");
    }
  }

  public InstanceViewer() {
    
    setLayout(new BorderLayout());
    m_UpdateString = "";
    setClearEachDataset(true);
    m_OutputTex = new JTextArea(10,20);
    m_OutputTex.setEditable(false);
    add("Center", new JScrollPane(m_OutputTex));
  }

  public void setClearEachDataset(boolean clear) {
    
    m_Clear = clear;
  }
  
  public boolean getClearEachDataset() {
    
    return m_Clear;
  }
  
  public void setDebug(boolean debug) {
    
    m_Debug = debug;
  }
  
  public boolean getDebug() {
    
    return m_Debug;
  }

  public void instanceProduced(InstanceEvent e) {
    
    Object source = e.getSource();
    if (source instanceof InstanceProducer) { 
      try {
	InstanceProducer a = (InstanceProducer) source;
	switch (e.getID()) {
	case InstanceEvent.FORMAT_AVAILABLE:
	  inputFormat(a.outputFormat());
	  break;
	case InstanceEvent.INSTANCE_AVAILABLE:
	  input(a.outputPeek());
	  break;
	case InstanceEvent.BATCH_FINISHED:
	  batchFinished();
	  break;
	default:
	  System.err.println("InstanceViewer::instanceProduced()"
			     + " - unknown event type");
	  break;
	}
      } catch (Exception ex) {
	System.err.println(ex.getMessage());
      }
    } else {
      System.err.println("InstanceViewer::instanceProduced()"
			 + " - Unknown source object type");
    }
  }
}
