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
 * ArffViewer.java
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.arffviewer;

import weka.core.Memory;
import weka.gui.ComponentHelper;
import weka.gui.LookAndFeel;

import java.awt.BorderLayout;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;

import javax.swing.JFrame;
import javax.swing.JOptionPane;

/**
 * A little tool for viewing ARFF files.
 *
 *
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 4719 $ 
 */

public class ArffViewer 
  extends JFrame
  implements WindowListener {
  
  /** for serialization */
  static final long serialVersionUID = -7455845566922685175L;

  /** the main panel */
  private ArffViewerMainPanel m_MainPanel;
  
  /** for monitoring the Memory consumption */
  private static Memory m_Memory = new Memory(true);
  
  /** the viewer if started from command line */
  private static ArffViewer m_Viewer;

  /** whether the files were already loaded */
  private static boolean m_FilesLoaded;

  /** the command line arguments */
  private static String[] m_Args;
  
  /**
   * initializes the object
   */
  public ArffViewer() {
    super("ARFF-Viewer");
    createFrame();
  }
  
  /**
   * creates all the components in the frame
   */
  protected void createFrame() {
    // basic setup
    setIconImage(ComponentHelper.getImage("weka_icon.gif"));
    setSize(ArffViewerMainPanel.WIDTH, ArffViewerMainPanel.HEIGHT);
    setCenteredLocation();
    setDefaultCloseOperation(DISPOSE_ON_CLOSE);
    // remove the listener - otherwise we get the strange behavior that one
    // frame receives a window-event for every single open frame!
    removeWindowListener(this);
    // add listener anew
    addWindowListener(this);
    
    getContentPane().setLayout(new BorderLayout());
    
    m_MainPanel = new ArffViewerMainPanel(this);
    m_MainPanel.setConfirmExit(false);
    getContentPane().add(m_MainPanel, BorderLayout.CENTER);
    
    setJMenuBar(m_MainPanel.getMenu());
  }
  
  /**
   * returns the left coordinate if the frame would be centered
   * 
   * @return 		the left coordinate
   */
  protected int getCenteredLeft() {
    int            width;
    int            x;
    
    width  = getBounds().width;
    x      = (getGraphicsConfiguration().getBounds().width  - width) / 2;
    
    if (x < 0) 
      x = 0;
    
    return x;
  }
  
  /**
   * returns the top coordinate if the frame would be centered
   * 
   * @return		the top coordinate
   */
  protected int getCenteredTop() {
    int            height;
    int            y;
    
    height = getBounds().height;
    y      = (getGraphicsConfiguration().getBounds().height - height) / 2;
    
    if (y < 0) 
      y = 0;
    
    return y;
  }
  
  /**
   * positions the window at the center of the screen
   */
  public void setCenteredLocation() { 
    setLocation(getCenteredLeft(), getCenteredTop());
  }
  
  /**
   * whether to present a MessageBox on Exit or not
   * @param confirm           whether a MessageBox pops up or not to confirm
   *                          exit
   */
  public void setConfirmExit(boolean confirm) {
    m_MainPanel.setConfirmExit(confirm);
  }
  
  /**
   * returns the setting of whether to display a confirm messagebox or not
   * on exit
   * @return                  whether a messagebox is displayed or not
   */
  public boolean getConfirmExit() {
    return m_MainPanel.getConfirmExit();
  }

  /**
   * whether to do a System.exit(0) on close
   * 
   * @param value	enables/disables the System.exit(0)
   */
  public void setExitOnClose(boolean value) {
    m_MainPanel.setExitOnClose(value);
  }

  /**
   * returns TRUE if a System.exit(0) is done on a close
   * 
   * @return		true if System.exit(0) is done
   */
  public boolean getExitOnClose() {
    return m_MainPanel.getExitOnClose();
  }
  
  /**
   * returns the main panel
   * 
   * @return		the main panel
   */
  public ArffViewerMainPanel getMainPanel() {
    return m_MainPanel;
  }
  
  /**
   * validates and repaints the frame
   */
  public void refresh() {
    validate();
    repaint();
  }

  /**
   * invoked when a window is activated
   * 
   * @param e		the window event
   */
  public void windowActivated(WindowEvent e) {
  }
  
  /**
   * invoked when a window is closed
   * 
   * @param e		the window event
   */
  public void windowClosed(WindowEvent e) {
  }
  
  /**
   * invoked when a window is in the process of closing
   * 
   * @param e		the window event
   */
  public void windowClosing(WindowEvent e) {
    int         button;
    
    while (getMainPanel().getTabbedPane().getTabCount() > 0)
      getMainPanel().closeFile(false);
    
    if (getConfirmExit()) {
      button = ComponentHelper.showMessageBox(
          this,
          "Quit - " + getTitle(),
          "Do you really want to quit?",
          JOptionPane.YES_NO_OPTION,
          JOptionPane.QUESTION_MESSAGE);
      if (button == JOptionPane.YES_OPTION)
        dispose();
    } 
    else {
      dispose();
    }

    if (getExitOnClose())
      System.exit(0);
  }
  
  /**
   * invoked when a window is deactivated
   * 
   * @param e		the window event
   */
  public void windowDeactivated(WindowEvent e) {
  }
  
  /**
   * invoked when a window is deiconified
   * 
   * @param e		the window event
   */
  public void windowDeiconified(WindowEvent e) {
  }
  
  /**
   * invoked when a window is iconified
   * 
   * @param e		the window event
   */
  public void windowIconified(WindowEvent e) {
  }
  
  /**
   * invoked when a window is has been opened
   * 
   * @param e		the window event
   */
  public void windowOpened(WindowEvent e) {
  }
  
  /**
   * returns only the classname
   * 
   * @return 		the classname
   */
  public String toString() {
    return this.getClass().getName();
  }
  
  /**
   * shows the frame and it tries to load all the arff files that were
   * provided as arguments.
   * 
   * @param args	the commandline parameters
   * @throws Exception	if something goes wrong
   */
  public static void main(String[] args) throws Exception {
    weka.core.logging.Logger.log(weka.core.logging.Logger.Level.INFO, "Logging started");
    LookAndFeel.setLookAndFeel();
    
    try {
      // uncomment to disable the memory management:
      //m_Memory.setEnabled(false);

      m_Viewer      = new ArffViewer();
      m_Viewer.setExitOnClose(true);
      m_Viewer.setVisible(true);
      m_FilesLoaded = false;
      m_Args        = args;

      Thread memMonitor = new Thread() {
        public void run() {
          while(true) {
            try {
              if ( (m_Args.length > 0) && (!m_FilesLoaded) ) {
                for (int i = 0; i < m_Args.length; i++) {
                  System.out.println("Loading " + (i+1) + "/" 
                      + m_Args.length +  ": '" + m_Args[i] + "'...");
                  m_Viewer.getMainPanel().loadFile(m_Args[i]);
                }
                m_Viewer.getMainPanel().getTabbedPane().setSelectedIndex(0);
                System.out.println("Finished!");
                m_FilesLoaded = true;
              }

              //System.out.println("before sleeping");
              Thread.sleep(4000);
              
              System.gc();
              
              if (m_Memory.isOutOfMemory()) {
                // clean up
                m_Viewer.dispose();
                m_Viewer = null;
                System.gc();

                // stop threads
                m_Memory.stopThreads();

                // display error
                System.err.println("\ndisplayed message:");
                m_Memory.showOutOfMemory();
                System.err.println("\nrestarting...");

                // restart GUI
                System.gc();
                m_Viewer = new ArffViewer();
                m_Viewer.setExitOnClose(true);
                m_Viewer.setVisible(true);
                // Note: no re-loading of datasets, otherwise we could end up
                //       in an endless loop!
              }
            }
            catch(InterruptedException ex) { 
              ex.printStackTrace(); 
            }
          }
        }
      };

      memMonitor.setPriority(Thread.NORM_PRIORITY);
      memMonitor.start();    
    } 
    catch (Exception ex) {
      ex.printStackTrace();
      System.err.println(ex.getMessage());
    }
  }
}
