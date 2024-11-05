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
 *    LogPanel.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui;

import java.awt.BorderLayout;
import java.awt.Point;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.text.SimpleDateFormat;
import java.util.Date;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JViewport;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

/** 
 * This panel allows log and status messages to be posted. Log messages
 * appear in a scrollable text area, and status messages appear as one-line
 * transient messages.
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 4789 $
 */
public class LogPanel
  extends JPanel
  implements Logger, TaskLogger {

  /** for serialization */
  private static final long serialVersionUID = -4072464549112439484L;

  /** Displays the current status */
  protected JLabel m_StatusLab = new JLabel("OK");
  
  /** Displays the log messages */
  protected JTextArea m_LogText = new JTextArea(4, 20);

  /** The button for viewing the log */
  protected JButton m_logButton = new JButton("Log");

  /** An indicator for whether text has been output yet */
  protected boolean m_First = true;

  /** The panel for monitoring the number of running tasks (if supplied)*/
  protected WekaTaskMonitor m_TaskMonitor=null;
  
  /**
   * Creates the log panel with no task monitor and
   * the log always visible.
   */
  public LogPanel() {

    this(null, false, false, true);
  }

  /**
   * Creates the log panel with a task monitor,
   * where the log is hidden.
   *
   * @param tm the task monitor, or null for none
   */
  public LogPanel(WekaTaskMonitor tm) {

    this(tm, true, false, true);
  }

  /**
   * Creates the log panel, possibly with task monitor,
   * where the log is optionally hidden.
   *
   * @param tm the task monitor, or null for none
   * @param logHidden true if the log should be hidden and
   *                  acessible via a button, or false if the
   *                  log should always be visible.
   */
  public LogPanel(WekaTaskMonitor tm, boolean logHidden) {
    this(tm, logHidden, false, true);
  }

  /**
   * Creates the log panel, possibly with task monitor,
   * where the either the log is optionally hidden or the status
   * (having both hidden is not allowed).
   * 
   *
   * @param tm the task monitor, or null for none
   * @param logHidden true if the log should be hidden and
   *                  acessible via a button, or false if the
   *                  log should always be visible.
   * @param statusHidden true if the status bar should be hidden (i.e.
   * @param titledBorder true if the log should have a title
   * you only want the log part).
   */
  public LogPanel(WekaTaskMonitor tm, boolean logHidden, 
      boolean statusHidden, boolean titledBorder) {

    m_TaskMonitor = tm;
    m_LogText.setEditable(false);
    m_LogText.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
    m_StatusLab.setBorder(BorderFactory.createCompoundBorder(
			  BorderFactory.createTitledBorder("Status"),
			  BorderFactory.createEmptyBorder(0, 5, 5, 5)));

    // create scrolling log
    final JScrollPane js = new JScrollPane(m_LogText);
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

    if (logHidden) {

      // create log window
      final JFrame jf = new JFrame("Log");
      jf.addWindowListener(new WindowAdapter() {
	  public void windowClosing(WindowEvent e) {
	    jf.setVisible(false);
	  }
	});
      jf.getContentPane().setLayout(new BorderLayout());
      jf.getContentPane().add(js, BorderLayout.CENTER);
      jf.pack();
      jf.setSize(450, 350);
      
      // display log window on request
      m_logButton.addActionListener(new ActionListener() {
	  public void actionPerformed(ActionEvent e) {
	    jf.setVisible(true);
	  }
	});
      
      // do layout
      setLayout(new BorderLayout());
      JPanel logButPanel = new JPanel();
      logButPanel.setLayout(new BorderLayout());
      logButPanel.setBorder(BorderFactory.createEmptyBorder(10, 5, 10, 5));
      logButPanel.add(m_logButton, BorderLayout.CENTER);
      JPanel p1 = new JPanel();
      p1.setLayout(new BorderLayout());
      p1.add(m_StatusLab, BorderLayout.CENTER);
      p1.add(logButPanel, BorderLayout.EAST);
      
      if (tm == null) {
	add(p1, BorderLayout.SOUTH);
      } else {
	JPanel p2 = new JPanel();
	p2.setLayout(new BorderLayout());
	p2.add(p1, BorderLayout.CENTER);
	p2.add((java.awt.Component)m_TaskMonitor, BorderLayout.EAST);
	add(p2, BorderLayout.SOUTH);
      }
    } else {
      // log always visible
      
      JPanel p1 = new JPanel();
      if (titledBorder) {
        p1.setBorder(BorderFactory.createTitledBorder("Log"));
      }
      p1.setLayout(new BorderLayout());
      p1.add(js, BorderLayout.CENTER);
      setLayout(new BorderLayout());
      add(p1, BorderLayout.CENTER);

      if (tm == null) {
        if (!statusHidden) {
          add(m_StatusLab, BorderLayout.SOUTH);
        }
      } else {
        if (!statusHidden) {
          JPanel p2 = new JPanel();
          p2.setLayout(new BorderLayout());
          p2.add(m_StatusLab,BorderLayout.CENTER);
          p2.add((java.awt.Component)m_TaskMonitor, BorderLayout.EAST);
          add(p2, BorderLayout.SOUTH);
        }
      }
    }
    addPopup();
  }

  /**
   * adds thousand's-separators to the number
   * @param l       the number to print
   * @return        the number as string with separators
   */
  private String printLong(long l) {
    String        result;
    String        str;
    int           i;
    int           count;

    str    = Long.toString(l);
    result = "";
    count  = 0;

    for (i = str.length() - 1; i >= 0; i--) {
      count++;
      result = str.charAt(i) + result;
      if ( (count == 3) && (i > 0) ) {
        result = "," + result;
        count = 0;
      }
    }
    
    return result;
  }

  /**
   * Add a popup menu for displaying the amount of free memory
   * and running the garbage collector
   */
  private void addPopup() {
    addMouseListener(new MouseAdapter() {
	public void mouseClicked(MouseEvent e) {
	  if (((e.getModifiers() & InputEvent.BUTTON1_MASK)
	       != InputEvent.BUTTON1_MASK) || e.isAltDown()) {
	    JPopupMenu gcMenu = new JPopupMenu();
	    JMenuItem availMem = new JMenuItem("Memory information");
	    availMem.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent ee) {
		  System.gc();
		  Runtime currR = Runtime.getRuntime();
		  long freeM = currR.freeMemory();
		  long totalM = currR.totalMemory();
		  long maxM = currR.maxMemory();
		  logMessage("Memory (free/total/max.) in bytes: " + printLong(freeM) + " / " + printLong(totalM) + " / " + printLong(maxM));
		  statusMessage("Memory (free/total/max.) in bytes: " + printLong(freeM) + " / " + printLong(totalM) + " / " + printLong(maxM));
		}
	      });
	    gcMenu.add(availMem);
	    JMenuItem runGC = new JMenuItem("Run garbage collector");
	    runGC.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent ee) {
		  statusMessage("Running garbage collector");
		  System.gc();
		  statusMessage("OK");
		}
	      });
	    gcMenu.add(runGC);
	    gcMenu.show(LogPanel.this, e.getX(), e.getY());
	  }
	}
      });
  }

  /**
   * Record the starting of a new task
   */
  public void taskStarted() {
    if (m_TaskMonitor != null) {
      m_TaskMonitor.taskStarted();
    }
  }

  /**
   * Record a task ending
   */
  public void taskFinished() {
    if (m_TaskMonitor != null) {
      m_TaskMonitor.taskFinished();
    }
  }
    
  /**
   * Gets a string containing current date and time.
   *
   * @return a string containing the date and time.
   */
  protected static String getTimestamp() {

    return (new SimpleDateFormat("HH:mm:ss:")).format(new Date());
  }

  /**
   * Sends the supplied message to the log area. The current timestamp will
   * be prepended.
   *
   * @param message a value of type 'String'
   */
  public synchronized void logMessage(String message) {

    if (m_First) {
      m_First = false;
    } else {
      m_LogText.append("\n");
    }
    m_LogText.append(LogPanel.getTimestamp() + ' ' + message);
    weka.core.logging.Logger.log(weka.core.logging.Logger.Level.INFO, message);
  }

  /**
   * Sends the supplied message to the status line.
   *
   * @param message the status message
   */
  public synchronized void statusMessage(String message) {

    m_StatusLab.setText(message);
  }

  
  /**
   * Tests out the log panel from the command line.
   *
   * @param args ignored
   */
  public static void main(String [] args) {

    try {
      final javax.swing.JFrame jf = new javax.swing.JFrame("Log Panel");
      jf.getContentPane().setLayout(new BorderLayout());
      final LogPanel lp = new LogPanel();
      jf.getContentPane().add(lp, BorderLayout.CENTER);
      jf.addWindowListener(new java.awt.event.WindowAdapter() {
	public void windowClosing(java.awt.event.WindowEvent e) {
	  jf.dispose();
	  System.exit(0);
	}
      });
      jf.pack();
      jf.setVisible(true);
      lp.logMessage("Welcome to the generic log panel!");
      lp.statusMessage("Hi there");
      lp.logMessage("Funky chickens");
      
    } catch (Exception ex) {
      ex.printStackTrace();
      System.err.println(ex.getMessage());
    }
  }
}
