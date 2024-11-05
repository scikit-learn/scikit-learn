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
 *    LogPanel
 *    Copyright (C) 2008 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.EventQueue;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.HashMap;
import java.util.Iterator;

import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.SwingUtilities;
import javax.swing.Timer;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableCellRenderer;

import weka.gui.Logger;

/**
 * Class for displaying a status area (made up of a variable number of
 * lines) and a log area.
 * 
 * @author mhall (mhall{[at]}pentaho{[dot]}com)
 * @version $Revision: 5563 $
 */
public class LogPanel extends JPanel implements Logger {
  
  /**
   * Holds the index (line number) in the JTable of each component
   * being tracked. 
   */
  private HashMap<String,Integer> m_tableIndexes = 
    new HashMap<String, Integer>();
  
  /**
   * Holds the timers associated with each component being tracked.
   */
  private HashMap<String, Timer> m_timers =
    new HashMap<String, Timer>();
  
  /**
   * The table model for the JTable used in the status area
   */
  private final DefaultTableModel m_tableModel;
  
  /**
   * The table for the status area
   */
  private JTable m_table;
  
  /**
   * Tabbed pane to hold both the status and the log
   */
  private JTabbedPane m_tabs = new JTabbedPane();
      
  /**
   * The log panel to delegate log messages to. 
   */
  private weka.gui.LogPanel m_logPanel = 
    new weka.gui.LogPanel(null, false, true, false);
    
  public LogPanel() {
    
    String[] columnNames = {"Component", "Parameters", "Time", "Status"};
    m_tableModel = new DefaultTableModel(columnNames, 0);
    
    // JTable with error/warning indication for rows.
    m_table = new JTable() {
      public Class getColumnClass(int column) {
        return getValueAt(0, column).getClass();
      }

      public Component prepareRenderer(TableCellRenderer renderer, 
          int row, int column) {
        Component c = super.prepareRenderer(renderer, row, column);
        if (!c.getBackground().equals(getSelectionBackground()))
        {
          String type = (String)getModel().getValueAt(row, 3);
          Color backgroundIndicator = null;
          if (type.startsWith("ERROR")) {
            backgroundIndicator = Color.RED;
          } else if (type.startsWith("WARNING")) {
            backgroundIndicator = Color.YELLOW;
          } else if (type.startsWith("INTERRUPTED")) {
            backgroundIndicator = Color.MAGENTA;
          }
          c.setBackground(backgroundIndicator);
        }
        return c;
      }
    };
    
    m_table.setModel(m_tableModel);
    m_table.getColumnModel().getColumn(0).setPreferredWidth(100);
    m_table.getColumnModel().getColumn(1).setPreferredWidth(150);
    m_table.getColumnModel().getColumn(2).setPreferredWidth(2);
    m_table.getColumnModel().getColumn(3).setPreferredWidth(500);
    m_table.setShowVerticalLines(true);
    
    JPanel statusPan = new JPanel();
    statusPan.setLayout(new BorderLayout());
    statusPan.add(new JScrollPane(m_table), BorderLayout.CENTER);
    m_tabs.addTab("Status", statusPan);
    m_tabs.addTab("Log", m_logPanel);
    
    setLayout(new BorderLayout());
    add(m_tabs, BorderLayout.CENTER);
    
  }
  
  /**
   * Clear the status area.
   */
  public void clearStatus() {
    // stop any running timers
    Iterator<Timer> i = m_timers.values().iterator();
    while (i.hasNext()) {
      i.next().stop();
    }
    
    // clear the map entries
    m_timers.clear();
    m_tableIndexes.clear();
    
    // clear the rows from the table
    while (m_tableModel.getRowCount() > 0) {
      m_tableModel.removeRow(0);
    }
  }
  
  /**
   * The JTable used for the status messages (in case clients
   * want to attach listeners etc.)
   * 
   * @return the JTable used for the status messages.
   */
  public JTable getStatusTable() {
    return m_table;
  }
  
  /**
   * Sends the supplied message to the log area. These message will typically
   * have the current timestamp prepended, and be viewable as a history.
   *
   * @param message the log message
   */
  public synchronized void logMessage(String message) {
    // delegate to the weka.gui.LogPanel
    m_logPanel.logMessage(message);
  }
  
  /**
   * Sends the supplied message to the status area. These messages are
   * typically one-line status messages to inform the user of progress
   * during processing (i.e. it doesn't matter if the user doesn't happen
   * to look at each message). These messages have the following format:
   * 
   * <Component name (needs to be unique)>|<Parameter string (optional)|<Status message>
   *
   * @param message the status message.
   */
  public synchronized void statusMessage(String message) {

    boolean hasDelimiters = (message.indexOf('|') > 0);
    String stepName = "";
    String stepHash = "";
    String stepParameters = "";
    String stepStatus = "";
    
    if (!hasDelimiters) {
      stepName = "Unknown";
      stepHash = "Unknown";
      stepStatus = message;
    } else {
      // Extract the fields of the status message
      stepHash = message.substring(0, message.indexOf('|'));
      message = message.substring(message.indexOf('|') + 1,
          message.length());
      // See if there is a unique object ID in the stepHash string
      if (stepHash.indexOf('$') > 0) {
        // Extract the step name
        stepName = stepHash.substring(0, stepHash.indexOf('$'));
      } else {
        stepName = stepHash;
      }
      
      // See if there are any step parameters to extract
      if (message.indexOf('|') > 0) {
        stepParameters = message.substring(0, message.indexOf('|'));
        stepStatus = message.substring(message.indexOf('|') + 1, 
            message.length());
      } else {
        // set the status message to the remainder
        stepStatus = message;
      }
    }
    
    // Now see if this step is in the hashmap
    if (m_tableIndexes.containsKey(stepHash)) {
      // Get the row number and update the table model...
      final Integer rowNum = m_tableIndexes.get(stepHash);
      if (stepStatus.equalsIgnoreCase("remove") ||
          stepStatus.equalsIgnoreCase("remove.")) {
        
        //m_tableModel.fireTableDataChanged();
        m_tableIndexes.remove(stepHash);
        m_timers.get(stepHash).stop();
        m_timers.remove(stepHash);
        
        // now need to decrement all the row indexes of
        // any rows greater than this one
        Iterator<String> i = m_tableIndexes.keySet().iterator();
        while (i.hasNext()) {
          String nextKey = i.next();
          int index = m_tableIndexes.get(nextKey).intValue();
          if (index > rowNum.intValue()) {
            index--;
            //System.err.println("*** " + nextKey + " needs decrementing to " + index);
            m_tableIndexes.put(nextKey, index);
//            System.err.println("new index " + m_rows.get(nextKey).intValue());
          }
        }
        
        // Remove the entry...
        if (!SwingUtilities.isEventDispatchThread()) {
          try {
            SwingUtilities.invokeLater(new Runnable() {
              public void run() {
                m_tableModel.removeRow(rowNum);
              }
            });
          } catch (Exception ex) {
            ex.printStackTrace();
          }
        } else {
          m_tableModel.removeRow(rowNum);
        }
      } else {
        final String stepNameCopy = stepName;
        final String stepStatusCopy = stepStatus;
        final String stepParametersCopy = stepParameters;

        if (!SwingUtilities.isEventDispatchThread()) {
          try {
            SwingUtilities.invokeLater(new Runnable() {
              public void run() {
                // ERROR overrides INTERRUPTED
                if (!(stepStatusCopy.startsWith("INTERRUPTED") &&
                    ((String)m_tableModel.getValueAt(rowNum.intValue(), 3)).startsWith("ERROR"))) {
                  m_tableModel.setValueAt(stepNameCopy, rowNum.intValue(), 0);
                  m_tableModel.setValueAt(stepParametersCopy, rowNum.intValue(), 1);
                  m_tableModel.setValueAt(m_table.getValueAt(rowNum.intValue(), 2), rowNum.intValue(), 2);
                  m_tableModel.setValueAt(stepStatusCopy, rowNum.intValue(), 3);
                }
              }
            });
          } catch (Exception ex) {
            ex.printStackTrace();
          }
        } else {
          if (!(stepStatusCopy.startsWith("INTERRUPTED") &&
              ((String)m_tableModel.getValueAt(rowNum.intValue(), 3)).startsWith("ERROR"))) {
            m_tableModel.setValueAt(stepNameCopy, rowNum.intValue(), 0);
            m_tableModel.setValueAt(stepParametersCopy, rowNum.intValue(), 1);
            m_tableModel.setValueAt(m_table.getValueAt(rowNum.intValue(), 2), rowNum.intValue(), 2);
            m_tableModel.setValueAt(stepStatusCopy, rowNum.intValue(), 3);
          }
        }
        if (stepStatus.startsWith("ERROR") ||
            stepStatus.startsWith("INTERRUPTED") ||
            stepStatus.equalsIgnoreCase("finished") ||
            stepStatus.equalsIgnoreCase("finished.") ||
            stepStatus.equalsIgnoreCase("done") ||
            stepStatus.equalsIgnoreCase("done.")) {
          // stop the timer.
          m_timers.get(stepHash).stop();
        } else if (!m_timers.get(stepHash).isRunning()) {
          // need to create a new one in order to reset the
          // elapsed time.
          installTimer(stepHash);
        }
      //  m_tableModel.fireTableCellUpdated(rowNum.intValue(), 3);
      }
    } else if (!stepStatus.equalsIgnoreCase("Remove") &&
        !stepStatus.equalsIgnoreCase("Remove.")) {
      // Add this one to the hash map
      int numKeys = m_tableIndexes.keySet().size();
      m_tableIndexes.put(stepHash, numKeys);
      
      // Now add a row to the table model
      final Object[] newRow = new Object[4];
      newRow[0] = stepName;
      newRow[1] = stepParameters;
      newRow[2] = "-";
      newRow[3] = stepStatus;
      final String stepHashCopy = stepHash;
      try {
        if (!SwingUtilities.isEventDispatchThread()) {
          SwingUtilities.invokeLater(new Runnable() {
            public void run() {
              m_tableModel.addRow(newRow);
              //m_tableModel.fireTableDataChanged();
            }
          });
        } else {
          m_tableModel.addRow(newRow);
        }
        
        installTimer(stepHashCopy);
      } catch (Exception ex) {
        ex.printStackTrace();
      }
    }
  }
  
  private void installTimer(final String stepHash) {
    final long startTime = System.currentTimeMillis();
    Timer newTimer = new Timer(1000, new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        synchronized(LogPanel.this) {
          if (m_tableIndexes.containsKey(stepHash)) {
            final Integer rn = m_tableIndexes.get(stepHash);
            long elapsed = System.currentTimeMillis() - startTime;
            long seconds = elapsed / 1000;
            long minutes = seconds / 60;
            final long hours = minutes / 60;
            seconds = seconds - (minutes * 60);
            minutes = minutes - (hours * 60);
            final long seconds2 = seconds;
            final long minutes2 = minutes;
            if (!SwingUtilities.isEventDispatchThread()) {
              try {
              SwingUtilities.invokeLater(new Runnable() {
                public void run() {
                  m_tableModel.
                    setValueAt("" + hours + ":" + minutes2 + ":" + seconds2, rn.intValue(), 2);
                }
              });
              } catch (Exception ex) {
                ex.printStackTrace();
              }
            } else {
              m_tableModel.
                setValueAt("" + hours + ":" + minutes2 + ":" + seconds2, rn.intValue(), 2);
            }
          }
        }
      }
    });
    m_timers.put(stepHash, newTimer);
    newTimer.start();
  }

  /**
   * Main method to test this class.
   * 
   * @param args any arguments (unused)
   */
  public static void main(String[] args) {
    try {
      final javax.swing.JFrame jf = new javax.swing.JFrame("Status/Log Panel");
      
      jf.getContentPane().setLayout(new BorderLayout());
      final LogPanel lp = new LogPanel();
      jf.getContentPane().add(lp, BorderLayout.CENTER);
      
      jf.getContentPane().add(lp, BorderLayout.CENTER);
      jf.addWindowListener(new java.awt.event.WindowAdapter() {
        public void windowClosing(java.awt.event.WindowEvent e) {
          jf.dispose();
          System.exit(0);
        }
      });
      jf.pack();
      jf.setVisible(true);
      lp.statusMessage("Step 1|Some options here|A status message");
      lp.statusMessage("Step 2$hashkey|Status message: no options");
      Thread.sleep(3000);
      lp.statusMessage("Step 2$hashkey|Funky Chickens!!!");
      Thread.sleep(3000);
      lp.statusMessage("Step 1|Some options here|finished");
      //lp.statusMessage("Step 1|Some options here|back again!");
      Thread.sleep(3000);
      lp.statusMessage("Step 2$hashkey|ERROR! More Funky Chickens!!!");
      Thread.sleep(3000);
      lp.statusMessage("Step 2$hashkey|WARNING - now a warning...");
      Thread.sleep(3000);
      lp.statusMessage("Step 2$hashkey|Back to normal.");
      Thread.sleep(3000);
      lp.statusMessage("Step 2$hashkey|INTERRUPTED.");
      
    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }
}
