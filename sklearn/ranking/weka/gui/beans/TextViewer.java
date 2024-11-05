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
 *    TextViewer.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

import weka.gui.ResultHistoryPanel;
import weka.gui.SaveBuffer;
import weka.gui.Logger;

import java.awt.BorderLayout;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.beans.EventSetDescriptor;
import java.beans.PropertyChangeListener;
import java.beans.VetoableChangeListener;
import java.beans.beancontext.BeanContext;
import java.beans.beancontext.BeanContextChild;
import java.beans.beancontext.BeanContextChildSupport;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Enumeration;
import java.util.Vector;

import javax.swing.BorderFactory;
import javax.swing.JFrame;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;

/**
 * Bean that collects and displays pieces of text
 *
 * @author <a href="mailto:mhall@cs.waikato.ac.nz">Mark Hall</a>
 * @version $Revision: 4806 $
 */
public class TextViewer 
  extends JPanel
  implements TextListener, DataSourceListener, 
	     TrainingSetListener, TestSetListener,
	     Visible, UserRequestAcceptor, 
	     BeanContextChild,
             BeanCommon,
             EventConstraints {

  /** for serialization */
  private static final long serialVersionUID = 104838186352536832L;

  protected BeanVisual m_visual;


  private transient JFrame m_resultsFrame = null;

  /**
   * Output area for a piece of text
   */
  private transient JTextArea m_outText = null;// = new JTextArea(20, 80);

  /**
   * List of text revieved so far
   */
  protected transient ResultHistoryPanel m_history;

  /**
   * True if this bean's appearance is the design mode appearance
   */
  protected boolean m_design;
  
  /**
   * BeanContex that this bean might be contained within
   */
  protected transient BeanContext m_beanContext = null;

  /**
   * BeanContextChild support
   */
  protected BeanContextChildSupport m_bcSupport = 
    new BeanContextChildSupport(this);

  /**
   * Objects listening for text events
   */
  private Vector m_textListeners = new Vector();

  private transient Logger m_log = null;

  
  public TextViewer() {
    /*    setUpResultHistory();
    setLayout(new BorderLayout());
    add(m_visual, BorderLayout.CENTER); */
    java.awt.GraphicsEnvironment ge = 
      java.awt.GraphicsEnvironment.getLocalGraphicsEnvironment(); 
    if (!ge.isHeadless()) {
      appearanceFinal();
    }
  }

  protected void appearanceDesign() {
    setUpResultHistory();
    removeAll();
    m_visual =  new BeanVisual("TextViewer", 
			       BeanVisual.ICON_PATH+"DefaultText.gif",
			       BeanVisual.ICON_PATH+"DefaultText_animated.gif");
    setLayout(new BorderLayout());
    add(m_visual, BorderLayout.CENTER);
  }

  protected void appearanceFinal() {
    removeAll();
    setLayout(new BorderLayout());
    setUpFinal();
  }

  protected void setUpFinal() {
    setUpResultHistory();
    JPanel holder = new JPanel();
    holder.setLayout(new BorderLayout());
    JScrollPane js = new JScrollPane(m_outText);
    js.setBorder(BorderFactory.createTitledBorder("Text"));
    holder.add(js, BorderLayout.CENTER);
    holder.add(m_history, BorderLayout.WEST);

    add(holder, BorderLayout.CENTER);
  }

  /**
   * Global info for this bean
   *
   * @return a <code>String</code> value
   */
  public String globalInfo() {
    return "General purpose text display.";
  }

  private void setUpResultHistory() {
    java.awt.GraphicsEnvironment ge = 
      java.awt.GraphicsEnvironment.getLocalGraphicsEnvironment(); 
    if(!ge.isHeadless()) {
      if (m_outText == null) {
        m_outText = new JTextArea(20, 80);
        m_history = new ResultHistoryPanel(m_outText);
      }
      m_outText.setEditable(false);
      m_outText.setFont(new Font("Monospaced", Font.PLAIN, 12));
      m_outText.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
      m_history.setBorder(BorderFactory.createTitledBorder("Result list"));
      m_history.setHandleRightClicks(false);
      m_history.getList().addMouseListener(new MouseAdapter() {
          public void mouseClicked(MouseEvent e) {
            if (((e.getModifiers() & InputEvent.BUTTON1_MASK)
                 != InputEvent.BUTTON1_MASK) || e.isAltDown()) {
              int index = m_history.getList().locationToIndex(e.getPoint());
              if (index != -1) {
                String name = m_history.getNameAtIndex(index);
                visualize(name, e.getX(), e.getY());
              } else {
                visualize(null, e.getX(), e.getY());
              }
            }
          }
        });
    }
  }

  /**
   * Handles constructing a popup menu with visualization options.
   * @param name the name of the result history list entry clicked on by
   * the user
   * @param x the x coordinate for popping up the menu
   * @param y the y coordinate for popping up the menu
   */
  protected void visualize(String name, int x, int y) {
    final JPanel panel = this;
    final String selectedName = name;
    JPopupMenu resultListMenu = new JPopupMenu();
    
    JMenuItem visMainBuffer = new JMenuItem("View in main window");
    if (selectedName != null) {
      visMainBuffer.addActionListener(new ActionListener() {
	  public void actionPerformed(ActionEvent e) {
	    m_history.setSingle(selectedName);
	  }
	});
    } else {
      visMainBuffer.setEnabled(false);
    }
    resultListMenu.add(visMainBuffer);
    
    JMenuItem visSepBuffer = new JMenuItem("View in separate window");
    if (selectedName != null) {
      visSepBuffer.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e) {
	  m_history.openFrame(selectedName);
	}
      });
    } else {
      visSepBuffer.setEnabled(false);
    }
    resultListMenu.add(visSepBuffer);
    
    JMenuItem saveOutput = new JMenuItem("Save result buffer");
    if (selectedName != null) {
      saveOutput.addActionListener(new ActionListener() {
	  public void actionPerformed(ActionEvent e) {
	    SaveBuffer m_SaveOut = new SaveBuffer(null, panel);
	    StringBuffer sb = m_history.getNamedBuffer(selectedName);
	    if (sb != null) {
	      m_SaveOut.save(sb);
	    }
	  }
	});
    } else {
      saveOutput.setEnabled(false);
    }
    resultListMenu.add(saveOutput);
    
    JMenuItem deleteOutput = new JMenuItem("Delete result buffer");
    if (selectedName != null) {
      deleteOutput.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e) {
	  m_history.removeResult(selectedName);
	}
      });
    } else {
      deleteOutput.setEnabled(false);
    }
    resultListMenu.add(deleteOutput);

    resultListMenu.show(m_history.getList(), x, y);
  }

  /**
   * Accept a data set for displaying as text
   *
   * @param e a <code>DataSetEvent</code> value
   */
  public synchronized void acceptDataSet(DataSetEvent e) {
    TextEvent nt = new TextEvent(e.getSource(), 
				 e.getDataSet().toString(),
				 e.getDataSet().relationName());
    acceptText(nt);
  }

  /**
   * Accept a training set for displaying as text
   *
   * @param e a <code>TrainingSetEvent</code> value
   */
  public synchronized void acceptTrainingSet(TrainingSetEvent e) {
    TextEvent nt = new TextEvent(e.getSource(), 
				 e.getTrainingSet().toString(),
				 e.getTrainingSet().relationName());
    acceptText(nt);
  }

  /**
   * Accept a test set for displaying as text
   *
   * @param e a <code>TestSetEvent</code> value
   */
  public synchronized void acceptTestSet(TestSetEvent e) {
    TextEvent nt = new TextEvent(e.getSource(), 
				 e.getTestSet().toString(),
				 e.getTestSet().relationName());
    acceptText(nt);
  }

  /**
   * Accept some text
   *
   * @param e a <code>TextEvent</code> value
   */
  public synchronized void acceptText(TextEvent e) {
    if (m_outText == null) {
      setUpResultHistory();
    }
    StringBuffer result = new StringBuffer();
    result.append(e.getText());
    //    m_resultsString.append(e.getText());
    //    m_outText.setText(m_resultsString.toString());
    String name = (new SimpleDateFormat("HH:mm:ss - "))
      .format(new Date());
    name += e.getTextTitle();
    //    System.err.println(name);
    if (name.length() > 30) {
      name = name.substring(0, 30);
    }

    if (m_outText != null) {
      // see if there is an entry with this name already in the list -
      // could happen if two items with the same name arrive at the same second
      int mod = 2;
      String nameOrig = new String(name);
      while (m_history.getNamedBuffer(name) != null) {
        name = new String(nameOrig+""+mod);
        mod++;
      }
      m_history.addResult(name, result);
      m_history.setSingle(name);
    }

    // pass on the event to any listeners
    notifyTextListeners(e);
  }

  /**
   * Describe <code>setVisual</code> method here.
   *
   * @param newVisual a <code>BeanVisual</code> value
   */
  public void setVisual(BeanVisual newVisual) {
    m_visual = newVisual;
  }

  /**
   * Get the visual appearance of this bean
   */
  public BeanVisual getVisual() {
    return m_visual;
  }
  
  /**
   * Use the default visual appearance for this bean
   */
  public void useDefaultVisual() {
    m_visual.loadIcons(BeanVisual.ICON_PATH+"DefaultText.gif",
		       BeanVisual.ICON_PATH+"DefaultText_animated.gif");
  }

  /**
   * Popup a component to display the selected text
   */
  public void showResults() {
    if (m_resultsFrame == null) {
      if (m_outText == null) {
	setUpResultHistory();
      }
      m_resultsFrame = new JFrame("Text Viewer");
      m_resultsFrame.getContentPane().setLayout(new BorderLayout());
      final JScrollPane js = new JScrollPane(m_outText);
      js.setBorder(BorderFactory.createTitledBorder("Text"));
      m_resultsFrame.getContentPane().add(js, BorderLayout.CENTER);
      m_resultsFrame.getContentPane().add(m_history, BorderLayout.WEST);
      m_resultsFrame.addWindowListener(new java.awt.event.WindowAdapter() {
	  public void windowClosing(java.awt.event.WindowEvent e) {
	    m_resultsFrame.dispose();
	    m_resultsFrame = null;
	  }
	});
      m_resultsFrame.pack();
      m_resultsFrame.setVisible(true);
    } else {
      m_resultsFrame.toFront();
    }
  }

  /**
   * Get a list of user requests
   *
   * @return an <code>Enumeration</code> value
   */
  public Enumeration enumerateRequests() {
    Vector newVector = new Vector(0);

    newVector.addElement("Show results");

    newVector.addElement("?Clear results");
    return newVector.elements();
  }

  /**
   * Perform the named request
   *
   * @param request a <code>String</code> value
   * @exception IllegalArgumentException if an error occurs
   */
  public void performRequest(String request) {
    if (request.compareTo("Show results") == 0) {
      showResults();
    } else if (request.compareTo("Clear results") == 0) {
      m_outText.setText("");
      m_history.clearResults();
    } else {
      throw new 
	IllegalArgumentException(request
		    + " not supported (TextViewer)");
    }
  }

  /**
   * Add a property change listener to this bean
   *
   * @param name the name of the property of interest
   * @param pcl a <code>PropertyChangeListener</code> value
   */
  public void addPropertyChangeListener(String name,
					PropertyChangeListener pcl) {
    m_bcSupport.addPropertyChangeListener(name, pcl);
  }

  /**
   * Remove a property change listener from this bean
   *
   * @param name the name of the property of interest
   * @param pcl a <code>PropertyChangeListener</code> value
   */
  public void removePropertyChangeListener(String name,
					   PropertyChangeListener pcl) {
    m_bcSupport.removePropertyChangeListener(name, pcl);
  }

  /**
   * Add a vetoable change listener to this bean
   *
   * @param name the name of the property of interest
   * @param vcl a <code>VetoableChangeListener</code> value
   */
  public void addVetoableChangeListener(String name,
				       VetoableChangeListener vcl) {
    m_bcSupport.addVetoableChangeListener(name, vcl);
  }
  
  /**
   * Remove a vetoable change listener from this bean
   *
   * @param name the name of the property of interest
   * @param vcl a <code>VetoableChangeListener</code> value
   */
  public void removeVetoableChangeListener(String name,
					   VetoableChangeListener vcl) {
    m_bcSupport.removeVetoableChangeListener(name, vcl);
  }

  /**
   * Set a bean context for this bean
   *
   * @param bc a <code>BeanContext</code> value
   */
  public void setBeanContext(BeanContext bc) {
    m_beanContext = bc;
    m_design = m_beanContext.isDesignTime();
    if (m_design) {
      appearanceDesign();
    } else {
      java.awt.GraphicsEnvironment ge = 
        java.awt.GraphicsEnvironment.getLocalGraphicsEnvironment(); 
      if (!ge.isHeadless()){
        appearanceFinal();
      }
    }
  }

  /**
   * Notify all text listeners of a text event
   *
   * @param ge a <code>TextEvent</code> value
   */
  private void notifyTextListeners(TextEvent ge) {
    Vector l;
    synchronized (this) {
      l = (Vector)m_textListeners.clone();
    }
    if (l.size() > 0) {
      for(int i = 0; i < l.size(); i++) {
	((TextListener)l.elementAt(i)).acceptText(ge);
      }
    }
  }

  /**
   * Return the bean context (if any) that this bean is embedded in
   *
   * @return a <code>BeanContext</code> value
   */
  public BeanContext getBeanContext() {
    return m_beanContext;
  }

  /**
   * Stop any processing that the bean might be doing.
   */
  public void stop() {
  }
  
  /**
   * Returns true if. at this time, the bean is busy with some
   * (i.e. perhaps a worker thread is performing some calculation).
   * 
   * @return true if the bean is busy.
   */
  public boolean isBusy() {
    return false;
  }

  /**
   * Set a logger
   *
   * @param logger a <code>Logger</code> value
   */
  public void setLog(Logger logger) {
    m_log = logger;
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
   * Returns true if, at this time, 
   * the object will accept a connection via the supplied
   * EventSetDescriptor
   *
   * @param esd the EventSetDescriptor
   * @return true if the object will accept a connection
   */
  public boolean connectionAllowed(EventSetDescriptor esd) {
    return connectionAllowed(esd.getName());
  }

  /**
   * Returns true if, at this time, 
   * the object will accept a connection via the named event
   *
   * @param eventName the name of the event
   * @return true if the object will accept a connection
   */
  public boolean connectionAllowed(String eventName) {
    return true;
  }

  /**
   * Notify this object that it has been registered as a listener with
   * a source for recieving events described by the named event
   * This object is responsible for recording this fact.
   *
   * @param eventName the event
   * @param source the source with which this object has been registered as
   * a listener
   */
  public void connectionNotification(String eventName, Object source) {
  }

  /**
   * Notify this object that it has been deregistered as a listener with
   * a source for named event. This object is responsible
   * for recording this fact.
   *
   * @param eventName the event
   * @param source the source with which this object has been registered as
   * a listener
   */
  public void disconnectionNotification(String eventName, Object source) {
  }

  /**
   * Returns true, if at the current time, the named event could
   * be generated. Assumes that the supplied event name is
   * an event that could be generated by this bean
   *
   * @param eventName the name of the event in question
   * @return true if the named event could be generated at this point in
   * time
   */
  public boolean eventGeneratable(String eventName) {
    if (eventName.equals("text")) {
      return true;
    }
    return false;
  }

  /**
   * Add a text listener
   *
   * @param cl a <code>TextListener</code> value
   */
  public synchronized void addTextListener(TextListener cl) {
    m_textListeners.addElement(cl);
  }

  /**
   * Remove a text listener
   *
   * @param cl a <code>TextListener</code> value
   */
  public synchronized void removeTextListener(TextListener cl) {
    m_textListeners.remove(cl);
  }

  public static void main(String [] args) {
    try {
      final javax.swing.JFrame jf = new javax.swing.JFrame();
      jf.getContentPane().setLayout(new java.awt.BorderLayout());

      final TextViewer tv = new TextViewer();

      tv.acceptText(new TextEvent(tv, "Here is some test text from the main "
				  +"method of this class.", "The Title"));
      jf.getContentPane().add(tv, java.awt.BorderLayout.CENTER);
      jf.addWindowListener(new java.awt.event.WindowAdapter() {
        public void windowClosing(java.awt.event.WindowEvent e) {
          jf.dispose();
          System.exit(0);
        }
      });
      jf.setSize(800,600);
      jf.setVisible(true);
    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }
}
