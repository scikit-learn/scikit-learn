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
 *    Explorer.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.explorer;

import weka.core.Capabilities;
import weka.core.Copyright;
import weka.core.Instances;
import weka.core.Memory;
import weka.core.converters.AbstractFileLoader;
import weka.core.converters.ConverterUtils;
import weka.gui.LogPanel;
import weka.gui.Logger;
import weka.gui.LookAndFeel;
import weka.gui.WekaTaskMonitor;

import java.awt.BorderLayout;
import java.awt.Image;
import java.awt.Toolkit;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.EventListener;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Vector;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.event.ChangeEvent;

/** 
 * The main class for the Weka explorer. Lets the user create,
 * open, save, configure, datasets, and perform ML analysis.
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 6778 $
 */
public class Explorer
  extends JPanel {

  /** for serialization */
  private static final long serialVersionUID = -7674003708867909578L;

  /**
   * Interface for classes that listen for filter changes.
   * 
   * @author FracPete (fracpete at waikato dot ac dot nz)
   * @version $Revision: 6778 $
   */
  public static interface CapabilitiesFilterChangeListener 
    extends EventListener {
    
    /**
     * method gets called in case of a change event
     * 
     * @param e		the associated change event
     */
    public void capabilitiesFilterChanged(CapabilitiesFilterChangeEvent e);
  }

  /**
   * This event can be fired in case the capabilities filter got changed 
   * 
   * @author FracPete (fracpete at waikato dot ac dot nz)
   * @version $Revision: 6778 $
   */
  public static class CapabilitiesFilterChangeEvent
    extends ChangeEvent {

    /** for serialization */
    private static final long serialVersionUID = 1194260517270385559L;
    
    /** the capabilities filter */
    protected Capabilities m_Filter;
    
    /**
     * Constructs a GOECapabilitiesFilterChangeEvent object.
     * 
     * @param source	the Object that is the source of the event
     * @param filter	the responsible capabilities filter
     */
    public CapabilitiesFilterChangeEvent(Object source, Capabilities filter) {
      super(source);
      m_Filter = filter;
    }
    
    /**
     * returns the associated Capabilities filter
     * 
     * @return		the filter
     */
    public Capabilities getFilter() {
      return m_Filter;
    }
  }

  /**
   * A common interface for panels to be displayed in the Explorer
   * 
   * @author FracPete (fracpete at waikato dot ac dot nz)
   * @version $Revision: 6778 $
   */
  public static interface ExplorerPanel {

    /**
     * Sets the Explorer to use as parent frame (used for sending notifications
     * about changes in the data)
     * 
     * @param parent	the parent frame
     */
    public void setExplorer(Explorer parent);
    
    /**
     * returns the parent Explorer frame
     * 
     * @return		the parent
     */
    public Explorer getExplorer();
    
    /**
     * Tells the panel to use a new set of instances.
     *
     * @param inst a set of Instances
     */
    public void setInstances(Instances inst);
    
    /**
     * Returns the title for the tab in the Explorer
     * 
     * @return the title of this tab
     */
    public String getTabTitle();
    
    /**
     * Returns the tooltip for the tab in the Explorer
     * 
     * @return the tooltip of this tab
     */
    public String getTabTitleToolTip();
  }

  /**
   * A common interface for panels in the explorer that can handle logs
   * 
   * @author FracPete (fracpete at waikato dot ac dot nz)
   * @version $Revision: 6778 $
   */
  public static interface LogHandler {
    
    /**
     * Sets the Logger to receive informational messages
     *
     * @param newLog the Logger that will now get info messages
     */
    public void setLog(Logger newLog);
  }

  /** The panel for preprocessing instances */
  protected PreprocessPanel m_PreprocessPanel = new PreprocessPanel();
  
  /** Contains all the additional panels apart from the pre-processing panel */
  protected Vector<ExplorerPanel> m_Panels = new Vector<ExplorerPanel>();
  
  /** The tabbed pane that controls which sub-pane we are working with */
  protected JTabbedPane m_TabbedPane = new JTabbedPane();
  
  /** The panel for log and status messages */
  protected LogPanel m_LogPanel = new LogPanel(new WekaTaskMonitor());
  
  /** the listeners that listen to filter changes */
  protected HashSet<CapabilitiesFilterChangeListener> m_CapabilitiesFilterChangeListeners = new HashSet<CapabilitiesFilterChangeListener>();
  
  /**
   * Creates the experiment environment gui with no initial experiment
   */
  public Explorer() {
    
    String date = (new SimpleDateFormat("EEEE, d MMMM yyyy")).format(new Date());
    m_LogPanel.logMessage("Weka Explorer");
    m_LogPanel.logMessage("(c) " + Copyright.getFromYear() + "-" + Copyright.getToYear() 
	+ " " + Copyright.getOwner() + ", " + Copyright.getAddress());
    m_LogPanel.logMessage("web: " + Copyright.getURL());
    m_LogPanel.logMessage("Started on " + date);
    m_LogPanel.statusMessage("Welcome to the Weka Explorer");

    // intialize pre-processpanel
    m_PreprocessPanel.setLog(m_LogPanel);
    m_TabbedPane.addTab(
	m_PreprocessPanel.getTabTitle(),
	null,
	m_PreprocessPanel,
	m_PreprocessPanel.getTabTitleToolTip());
    
    // initialize additional panels
    String[] tabs = ExplorerDefaults.getTabs();
    Hashtable<String, HashSet> tabOptions = new Hashtable<String, HashSet>();
    for (int i = 0; i < tabs.length; i++) {
      try {
	// determine classname and additional options
	String[] optionsStr = tabs[i].split(":");
	String classname = optionsStr[0];
	HashSet options = new HashSet();
	tabOptions.put(classname, options);
	for (int n = 1; n < optionsStr.length; n++)
	  options.add(optionsStr[n]);
	  
	// setup panel
	ExplorerPanel panel = (ExplorerPanel) Class.forName(classname).newInstance();
	panel.setExplorer(this);
	m_Panels.add(panel);
	if (panel instanceof LogHandler)
	  ((LogHandler) panel).setLog(m_LogPanel);
	m_TabbedPane.addTab(
	    panel.getTabTitle(), null, (JPanel) panel, panel.getTabTitleToolTip());
      }
      catch (Exception e) {
	e.printStackTrace();
      }
    }

    // setup tabbed pane
    m_TabbedPane.setSelectedIndex(0);
    for (int i = 0; i < m_Panels.size(); i++) {
      HashSet options = tabOptions.get(m_Panels.get(i).getClass().getName());
      m_TabbedPane.setEnabledAt(i + 1, options.contains("standalone"));
    }

    // setup notification for dataset changes
    m_PreprocessPanel.addPropertyChangeListener(new PropertyChangeListener() {
      public void propertyChange(PropertyChangeEvent e) {
	for (int i = 0; i < m_Panels.size(); i++) {
	   m_Panels.get(i).setInstances(m_PreprocessPanel.getInstances());
	   m_TabbedPane.setEnabledAt(i + 1, true);
	}
      }
    });

    // add listeners for changes in the capabilities
    m_PreprocessPanel.setExplorer(this);
    addCapabilitiesFilterListener(m_PreprocessPanel);
    for (int i = 0; i < m_Panels.size(); i++) {
      if (m_Panels.get(i) instanceof CapabilitiesFilterChangeListener)
	addCapabilitiesFilterListener((CapabilitiesFilterChangeListener) m_Panels.get(i));
    }

    // add components to layout
    setLayout(new BorderLayout());
    add(m_TabbedPane, BorderLayout.CENTER);
    add(m_LogPanel, BorderLayout.SOUTH);
  }
  
  /**
   * returns all the panels, apart from the PreprocessPanel
   * 
   * @return		the currently displayed panels w/o PreprocessPanel
   */
  public Vector<ExplorerPanel> getPanels() {
    return m_Panels;
  }
  
  /**
   * returns the instance of the PreprocessPanel being used in this instance
   * of the Explorer
   * 
   * @return		the panel
   */
  public PreprocessPanel getPreprocessPanel() {
    return m_PreprocessPanel;
  }
  
  /**
   * returns the tabbed pane of the Explorer
   * 
   * @return		the tabbed pane
   */
  public JTabbedPane getTabbedPane() {
    return m_TabbedPane;
  }
  
  /**
   * adds the listener to the list of objects that listen for changes of the
   * CapabilitiesFilter
   * 
   * @param l		the listener to add
   * @see		#m_CapabilitiesFilterChangeListeners
   */
  public void addCapabilitiesFilterListener(CapabilitiesFilterChangeListener l) {
    m_CapabilitiesFilterChangeListeners.add(l);
  }

  /**
   * Removes the specified listener from the set of listeners if it is present.
   * 
   * @param l		the listener to remove
   * @return		true if the listener was registered
   */
  public boolean removeCapabilitiesFilterListener(CapabilitiesFilterChangeListener l) {
    return m_CapabilitiesFilterChangeListeners.remove(l);
  }
  
  /**
   * notifies all the listeners of a change
   * 
   * @param filter	the affected filter
   */
  public void notifyCapabilitiesFilterListener(Capabilities filter) {
    for (CapabilitiesFilterChangeListener l: m_CapabilitiesFilterChangeListeners) {
      if (l == this)
	continue;
      l.capabilitiesFilterChanged(new CapabilitiesFilterChangeEvent(this, filter));
    }
  }
  
  /** variable for the Explorer class which would be set to null by the memory 
      monitoring thread to free up some memory if we running out of memory
   */
  private static Explorer m_explorer;

  /** for monitoring the Memory consumption */
  private static Memory m_Memory = new Memory(true);

  /**
   * Tests out the explorer environment.
   *
   * @param args ignored.
   */
  public static void main(String [] args) {

    weka.core.logging.Logger.log(weka.core.logging.Logger.Level.INFO, "Logging started");
    
    LookAndFeel.setLookAndFeel();
    // make sure that packages are loaded and the GenericPropertiesCreator
    // executes to populate the lists correctly
    weka.gui.GenericObjectEditor.determineClasses();
    
    try {
      // uncomment to disable the memory management:
      //m_Memory.setEnabled(false);

      m_explorer = new Explorer();
      final JFrame jf = new JFrame("Weka Explorer");
      jf.getContentPane().setLayout(new BorderLayout());
      jf.getContentPane().add(m_explorer, BorderLayout.CENTER);
      jf.addWindowListener(new WindowAdapter() {
        public void windowClosing(WindowEvent e) {
          jf.dispose();
          System.exit(0);
        }
      });
      jf.pack();
      jf.setSize(800, 600);
      jf.setVisible(true);
      Image icon = Toolkit.getDefaultToolkit().
      getImage(ClassLoader.getSystemResource("weka/gui/weka_icon.gif"));
      jf.setIconImage(icon);

      if (args.length == 1) {
        System.err.println("Loading instances from " + args[0]);
        AbstractFileLoader loader = ConverterUtils.getLoaderForFile(args[0]);
	loader.setFile(new File(args[0]));
        m_explorer.m_PreprocessPanel.setInstancesFromFile(loader);
      }

      Thread memMonitor = new Thread() {
        public void run() {
          while(true) {
            try {
              //System.out.println("Before sleeping.");
              this.sleep(4000);

              System.gc();

              if (m_Memory.isOutOfMemory()) {
                // clean up
                jf.dispose();
                m_explorer = null;
                System.gc();

                // stop threads
                m_Memory.stopThreads();

                // display error
                System.err.println("\ndisplayed message:");
                m_Memory.showOutOfMemory();
                System.err.println("\nexiting");
                System.exit(-1);
              }

            } catch(InterruptedException ex) { ex.printStackTrace(); }
          }
        }
      };

      memMonitor.setPriority(Thread.MAX_PRIORITY);
      memMonitor.start();
    } catch (Exception ex) {
      ex.printStackTrace();
      System.err.println(ex.getMessage());
    }
  }
}
