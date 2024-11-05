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
 *    InstanceCounter.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.streams;

import weka.core.Instance;
import weka.core.Instances;

import java.awt.Color;

import javax.swing.JLabel;
import javax.swing.JPanel;

/** 
 * A bean that counts instances streamed to it.
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 1.5 $
 */
public class InstanceCounter
  extends JPanel
  implements InstanceListener {

  /** for serialization */
  private static final long serialVersionUID = -6084967152645935934L;
  
  private JLabel m_Count_Lab;
  private int m_Count;
  private boolean m_Debug;

  public void input(Instance instance) throws Exception {
    
    if (m_Debug) {
      System.err.println("InstanceCounter::input(" + instance +")");
    }
    m_Count++;
    m_Count_Lab.setText(""+m_Count+" instances");
    repaint();
  }
  
  public void inputFormat(Instances instanceInfo) {
    
    if (m_Debug) {
      System.err.println("InstanceCounter::inputFormat()");
    }
    Instances inputInstances = new Instances(instanceInfo,0);
    m_Count = 0;
    m_Count_Lab.setText(""+m_Count+" instances");
  }

  public void setDebug(boolean debug) {
    
    m_Debug = debug;
  }
  
  public boolean getDebug() {
    
    return m_Debug;
  }

  public InstanceCounter() {
    
    m_Count = 0;
    m_Count_Lab = new JLabel("no instances");
    add(m_Count_Lab);
    //    setSize(60,40);
    setBackground(Color.lightGray);
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
	  if (m_Debug)
	    System.err.println("InstanceCounter::instanceProduced() - End of instance batch");
	  break;
	default:
	  System.err.println("InstanceCounter::instanceProduced() - unknown event type");
	  break;
	}
      } catch (Exception ex) {
	System.err.println(ex.getMessage());
      }
    } else {
      System.err.println("InstanceCounter::instanceProduced() - Unknown source object type");
    }
  }
}








