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
 *    ResultHistoryPanel.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.gui;

import weka.gui.visualize.PrintableComponent;

import java.awt.BorderLayout;
import java.awt.Font;
import java.awt.Point;
import java.awt.event.InputEvent;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.Serializable;
import java.util.Hashtable;

import javax.swing.BorderFactory;
import javax.swing.DefaultListModel;
import javax.swing.JFrame;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JViewport;
import javax.swing.ListSelectionModel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.text.JTextComponent;

/** 
 * A component that accepts named stringbuffers and displays the name in a list
 * box. When a name is right-clicked, a frame is popped up that contains
 * the string held by the stringbuffer. Optionally a text component may be
 * provided that will have it's text set to the named result text on a
 * left-click.
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 1.26 $
 */
public class ResultHistoryPanel
  extends JPanel {
  
  /** for serialization */
  static final long serialVersionUID = 4297069440135326829L;
  
  /** An optional component for single-click display */
  protected JTextComponent m_SingleText;

  /** The named result being viewed in the single-click display */
  protected String m_SingleName;
  
  /** The list model */
  protected DefaultListModel m_Model = new DefaultListModel();

  /** The list component */
  protected JList m_List = new JList(m_Model);
  
  /** A Hashtable mapping names to result buffers */
  protected Hashtable m_Results = new Hashtable();

  /** A Hashtable mapping names to output text components */
  protected Hashtable m_FramedOutput = new Hashtable();

  /** A hashtable mapping names to arbitrary objects */
  protected Hashtable m_Objs = new Hashtable();

  /** Let the result history list handle right clicks in the default
      manner---ie, pop up a window displaying the buffer */
  protected boolean m_HandleRightClicks = true;

  /** for printing the output to files */
  protected PrintableComponent m_Printer = null;
  
  /**
   * Extension of MouseAdapter that implements Serializable.
   */
  public static class RMouseAdapter 
    extends MouseAdapter implements Serializable {
    
    /** for serialization */
    static final long serialVersionUID = -8991922650552358669L;    
  }
 
  
  /**
   * Extension of KeyAdapter that implements Serializable.
   */
  public static class RKeyAdapter 
    extends KeyAdapter implements Serializable {
    
    /** for serialization */
    static final long serialVersionUID = -8675332541861828079L;
  }

  /**
   * Create the result history object
   *
   * @param text the optional text component for single-click display
   */
  public ResultHistoryPanel(JTextComponent text) {
    m_SingleText = text;
    if (text != null) {
      m_Printer = new PrintableComponent(m_SingleText);
    }
    m_List.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
    m_List.addMouseListener(new RMouseAdapter() {
      private static final long serialVersionUID = -9015397020486290479L;
      
      public void mouseClicked(MouseEvent e) {
	if ((e.getModifiers() & InputEvent.BUTTON1_MASK)
	    == InputEvent.BUTTON1_MASK) {
            if (    ((e.getModifiers() & InputEvent.SHIFT_DOWN_MASK) == 0)
                 && ((e.getModifiers() & InputEvent.CTRL_DOWN_MASK) == 0) ) {
              int index = m_List.locationToIndex(e.getPoint());
              if ((index != -1) && (m_SingleText != null)) {
                setSingle((String)m_Model.elementAt(index));
            }
          }
	} else {
	  // if there are stored objects then assume that the storer
	  // will handle popping up the text in a seperate frame
	  if (m_HandleRightClicks) {
	    int index = m_List.locationToIndex(e.getPoint());
	    if (index != -1) {
	      String name = (String)m_Model.elementAt(index);
	      openFrame(name);
	    }
	  }
	}
      }
    });

    m_List.addKeyListener(new RKeyAdapter() {
      private static final long serialVersionUID = 7910681776999302344L;
      
      public void keyReleased(KeyEvent e) {
        if (e.getKeyCode() == KeyEvent.VK_DELETE) {
          int selected = m_List.getSelectedIndex();
          if (selected != -1) {
            removeResult((String)m_Model.elementAt(selected));
          }
        }
      }
    });
    m_List.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
      public void valueChanged(ListSelectionEvent e) {
	if (!e.getValueIsAdjusting()) {
	  ListSelectionModel lm = (ListSelectionModel) e.getSource();
	  for (int i = e.getFirstIndex(); i <= e.getLastIndex(); i++) {
	    if (lm.isSelectedIndex(i)) {
	      //m_AttSummaryPanel.setAttribute(i);
	       if ((i != -1) && (m_SingleText != null)) {
		 setSingle((String)m_Model.elementAt(i));
	       }
	      break;
	    }
	  }
	}
      }
    });


    setLayout(new BorderLayout());
    //    setBorder(BorderFactory.createTitledBorder("Result history"));
    final JScrollPane js = new JScrollPane(m_List);
    js.getViewport().addChangeListener(new ChangeListener() {
      private int lastHeight;
      public void stateChanged(ChangeEvent e) {
	JViewport vp = (JViewport)e.getSource();
	int h = vp.getViewSize().height; 
	if (h != lastHeight) { // i.e. an addition not just a user scrolling
	  lastHeight = h;
	  int x = h - vp.getExtentSize().height;
	  vp.setViewPosition(new Point(0, x));
	}
      }
    });
    add(js, BorderLayout.CENTER);
  }

  /**
   * Adds a new result to the result list.
   *
   * @param name the name to associate with the result
   * @param result the StringBuffer that contains the result text
   */
  public void addResult(String name, StringBuffer result) {
    
    m_Model.addElement(name);
    m_Results.put(name, result);
  }

  /**
   * Removes one of the result buffers from the history. Any windows currently
   * displaying the contents of the buffer are not affected.
   *
   * @param name the name of the buffer to remove.
   */
  public void removeResult(String name) {

    StringBuffer buff = (StringBuffer) m_Results.get(name);
    if (buff != null) {
      m_Results.remove(name);
      m_Model.removeElement(name);
      m_Objs.remove(name);
      System.gc();
    } 
  }

  /**
   * Removes all of the result buffers from the history. Any windows currently
   * displaying the contents of the buffer are not affected.
   */
  public void clearResults() {
    m_Results.clear();
    m_Model.clear();
    m_Objs.clear();
    System.gc();
  }

  /**
   * Adds an object to the results list
   * @param name the name to associate with the object
   * @param o the object
   */
  public void addObject(String name, Object o) {
    m_Objs.put(name, o);
  }

  /**
   * Get the named object from the list
   * @param name the name of the item to retrieve the stored object
   * for
   * @return the object or null if there is no object at this index
   */
  public Object getNamedObject(String name) {
    Object v = null;
    v = m_Objs.get(name);
    return v;
  }

  /**
   * Gets the object associated with the currently
   * selected item in the list.
   * @return the object or null if there is no
   * object corresponding to the current selection in
   * the list
   */
  public Object getSelectedObject() {
    Object v = null;
    int index = m_List.getSelectedIndex();
    if (index != -1) {
      String name = (String)(m_Model.elementAt(index));
      v = m_Objs.get(name);
    }
    
    return v;
  }

  /**
   * Gets the named buffer
   * @return the buffer or null if there are no items in
   * the list
   */
  public StringBuffer getNamedBuffer(String name) {
    StringBuffer b = null;
    b = (StringBuffer)(m_Results.get(name));
    return b;
  }

  /**
   * Gets the buffer associated with the currently
   * selected item in the list.
   * @return the buffer or null if there are no items in
   * the list
   */
  public StringBuffer getSelectedBuffer() {
    StringBuffer b = null;
    int index = m_List.getSelectedIndex();
    if (index != -1) {
      String name = (String)(m_Model.elementAt(index));
      b = (StringBuffer)(m_Results.get(name));
    }
    return b;
  }

  /**
   * Get the name of the currently selected item in the list
   * @return the name of the currently selected item or null if no
   * item selected
   */
  public String getSelectedName() {
    int index = m_List.getSelectedIndex();
    if (index != -1) {
      return (String)(m_Model.elementAt(index));
    }
    return null;
  }

  /**
   * Gets the name of theitem in the list at the specified index
   * @return the name of item or null if there is no item at that index
   */
  public String getNameAtIndex(int index) {
    if (index != -1) {
      return (String)(m_Model.elementAt(index));
    }
    return null;
  }

  /**
   * Sets the single-click display to view the named result.
   *
   * @param name the name of the result to display.
   */
  public void setSingle(String name) {

    StringBuffer buff = (StringBuffer) m_Results.get(name);
    if (buff != null) {
      m_SingleName = name;
      m_SingleText.setText(buff.toString());
      m_List.setSelectedValue(name, true);
    }
  }
  
  /**
   * Opens the named result in a separate frame.
   *
   * @param name the name of the result to open.
   */
  public void openFrame(String name) {

    StringBuffer buff = (StringBuffer) m_Results.get(name);
    JTextComponent currentText = (JTextComponent) m_FramedOutput.get(name);
    if ((buff != null) && (currentText == null)) {
      // Open the frame.
      JTextArea ta = new JTextArea();
      ta.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
      ta.setFont(new Font("Monospaced", Font.PLAIN, 12));
      ta.setEditable(false);
      ta.setText(buff.toString());
      m_FramedOutput.put(name, ta);
      final JFrame jf = new JFrame(name);
      jf.addWindowListener(new WindowAdapter() {
	public void windowClosing(WindowEvent e) {
	  m_FramedOutput.remove(jf.getTitle());
	  jf.dispose();
	}
      });
      jf.getContentPane().setLayout(new BorderLayout());
      jf.getContentPane().add(new JScrollPane(ta), BorderLayout.CENTER);
      jf.pack();
      jf.setSize(450, 350);
      jf.setVisible(true);
    }
  }

  /**
   * Tells any component currently displaying the named result that the
   * contents of the result text in the StringBuffer have been updated.
   *
   * @param name the name of the result that has been updated.
   */
  public void updateResult(String name) {

    StringBuffer buff = (StringBuffer) m_Results.get(name);
    if (buff == null) {
      return;
    }
    if (m_SingleName == name) {
      m_SingleText.setText(buff.toString());
    }
    JTextComponent currentText = (JTextComponent) m_FramedOutput.get(name);
    if (currentText != null) {
      currentText.setText(buff.toString());
    }
  }

  /**
   * Gets the selection model used by the results list.
   *
   * @return a value of type 'ListSelectionModel'
   */
  public ListSelectionModel getSelectionModel() {
    
    return m_List.getSelectionModel();
  }

  /**
   * Gets the JList used by the results list
   * @return the JList
   */
  public JList getList() {
    return m_List;
  }

  /**
   * Set whether the result history list should handle right clicks
   * or whether the parent object will handle them.
   * @param tf false if parent object will handle right clicks
   */
  public void setHandleRightClicks(boolean tf) {
    m_HandleRightClicks = tf;
  }


  /**
   * Tests out the result history from the command line.
   *
   * @param args ignored
   */
  public static void main(String [] args) {

    try {
      final javax.swing.JFrame jf =
	new javax.swing.JFrame("Weka Explorer: Classifier");
      jf.getContentPane().setLayout(new BorderLayout());
      final ResultHistoryPanel jd = new ResultHistoryPanel(null);
      jd.addResult("blah", new StringBuffer("Nothing to see here"));
      jd.addResult("blah1", new StringBuffer("Nothing to see here1"));
      jd.addResult("blah2", new StringBuffer("Nothing to see here2"));
      jd.addResult("blah3", new StringBuffer("Nothing to see here3"));
      jf.getContentPane().add(jd, BorderLayout.CENTER);
      jf.addWindowListener(new java.awt.event.WindowAdapter() {
	public void windowClosing(java.awt.event.WindowEvent e) {
	  jf.dispose();
	  System.exit(0);
	}
      });
      jf.pack();
      jf.setVisible(true);
    } catch (Exception ex) {
      ex.printStackTrace();
      System.err.println(ex.getMessage());
    }
  }
}
