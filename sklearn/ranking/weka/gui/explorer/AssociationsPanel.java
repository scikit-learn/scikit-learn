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
 *    AssociationsPanel.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.explorer;

import weka.associations.AssociationRule;
import weka.associations.AssociationRules;
import weka.associations.Associator;
import weka.core.Attribute;
import weka.core.Capabilities;
import weka.core.CapabilitiesHandler;
import weka.core.Drawable;
import weka.core.Instances;
import weka.core.OptionHandler;
import weka.core.Utils;
import weka.core.Version;
import weka.core.labelranking.PreferenceAttribute;
import weka.gui.GenericObjectEditor;
import weka.gui.Logger;
import weka.gui.PropertyPanel;
import weka.gui.ResultHistoryPanel;
import weka.gui.SaveBuffer;
import weka.gui.SysErrLog;
import weka.gui.TaskLogger;
import weka.gui.explorer.Explorer.CapabilitiesFilterChangeEvent;
import weka.gui.explorer.Explorer.CapabilitiesFilterChangeListener;
import weka.gui.explorer.Explorer.ExplorerPanel;
import weka.gui.explorer.Explorer.LogHandler;
import weka.gui.treevisualizer.PlaceNode2;
import weka.gui.treevisualizer.TreeVisualizer;
import weka.gui.visualize.plugins.TreeVisualizePlugin;
import weka.gui.visualize.plugins.AssociationRuleVisualizePlugin;

import java.awt.BorderLayout;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Point;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.List;
import java.util.Vector;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JViewport;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

/** 
 * This panel allows the user to select, configure, and run a scheme
 * that learns associations.
 *
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version $Revision: 6499 $
 */
public class AssociationsPanel 
  extends JPanel
  implements CapabilitiesFilterChangeListener, ExplorerPanel, LogHandler {
  
  /** for serialization */
  static final long serialVersionUID = -6867871711865476971L;

  /** the parent frame */
  protected Explorer m_Explorer = null;

  /** Lets the user configure the associator */
  protected GenericObjectEditor m_AssociatorEditor =
    new GenericObjectEditor();

  /** The panel showing the current associator selection */
  protected PropertyPanel m_CEPanel = new PropertyPanel(m_AssociatorEditor);
  
  /** The output area for associations */
  protected JTextArea m_OutText = new JTextArea(20, 40);

  /** The destination for log/status messages */
  protected Logger m_Log = new SysErrLog();

  /** The buffer saving object for saving output */
  protected SaveBuffer m_SaveOut = new SaveBuffer(m_Log, this);

  /** A panel controlling results viewing */
  protected ResultHistoryPanel m_History = new ResultHistoryPanel(m_OutText);

  /** Click to start running the associator */
  protected JButton m_StartBut = new JButton("Start");

  /** Click to stop a running associator */
  protected JButton m_StopBut = new JButton("Stop");
  
  /** 
   * Whether to store any graph or xml rules output in
   * the history list
   */
  protected JCheckBox m_storeOutput = 
    new JCheckBox("Store output for visualization");
  
  /** The main set of instances we're playing with */
  protected Instances m_Instances;

  /** The user-supplied test set (if any) */
  protected Instances m_TestInstances;
  
  /** A thread that associator runs in */
  protected Thread m_RunThread;

  /* Register the property editors we need */
  static {
     GenericObjectEditor.registerEditors();
  }
  
  /**
   * Creates the associator panel
   */
  public AssociationsPanel() {

    // Connect / configure the components
    m_OutText.setEditable(false);
    m_OutText.setFont(new Font("Monospaced", Font.PLAIN, 12));
    m_OutText.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
    m_OutText.addMouseListener(new MouseAdapter() {
      public void mouseClicked(MouseEvent e) {
	if ((e.getModifiers() & InputEvent.BUTTON1_MASK)
	    != InputEvent.BUTTON1_MASK) {
	  m_OutText.selectAll();
	}
      }
    });
    m_History.setBorder(BorderFactory.createTitledBorder("Result list (right-click for options)"));
    m_History.setHandleRightClicks(false);
    // see if we can popup a menu for the selected result
    m_History.getList().addMouseListener(new MouseAdapter() {
	public void mouseClicked(MouseEvent e) {
	  if (((e.getModifiers() & InputEvent.BUTTON1_MASK)
	       != InputEvent.BUTTON1_MASK) || e.isAltDown()) {
	    int index = m_History.getList().locationToIndex(e.getPoint());
	    if (index != -1) {
	      String name = m_History.getNameAtIndex(index);
	      historyRightClickPopup(name, e.getX(), e.getY());
	    } else {
	      historyRightClickPopup(null, e.getX(), e.getY());
	    }
	  }
	}
      });

    m_AssociatorEditor.setClassType(Associator.class);
    m_AssociatorEditor.setValue(ExplorerDefaults.getAssociator());
    m_AssociatorEditor.addPropertyChangeListener(new PropertyChangeListener() {
      public void propertyChange(PropertyChangeEvent e) {
        m_StartBut.setEnabled(true);
        // Check capabilities
        Capabilities currentFilter = m_AssociatorEditor.getCapabilitiesFilter();
        Associator associator = (Associator) m_AssociatorEditor.getValue();
        Capabilities currentSchemeCapabilities =  null;
        if (associator != null && currentFilter != null && 
            (associator instanceof CapabilitiesHandler)) {
          currentSchemeCapabilities = ((CapabilitiesHandler)associator).getCapabilities();
          
          if (!currentSchemeCapabilities.supportsMaybe(currentFilter) &&
              !currentSchemeCapabilities.supports(currentFilter)) {
            m_StartBut.setEnabled(false);
          }
        }
	repaint();
      }
    });

    m_StartBut.setToolTipText("Starts the associator");
    m_StopBut.setToolTipText("Stops the associator");
    m_StartBut.setEnabled(false);
    m_StopBut.setEnabled(false);
    m_StartBut.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
	startAssociator();
      }
    });
    m_StopBut.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
	stopAssociator();
      }
    });
    
    // check for any visualization plugins so that we
    // can add a checkbox for storing graphs or rules
    boolean showStoreOutput = 
      (GenericObjectEditor.
        getClassnames(AssociationRuleVisualizePlugin.class.getName()).size() > 0 ||
       GenericObjectEditor.
        getClassnames(TreeVisualizePlugin.class.getName()).size() > 0);

    // Layout the GUI
    JPanel p1 = new JPanel();
    p1.setBorder(BorderFactory.createCompoundBorder(
		 BorderFactory.createTitledBorder("Associator"),
		 BorderFactory.createEmptyBorder(0, 5, 5, 5)
		 ));
    p1.setLayout(new BorderLayout());
    p1.add(m_CEPanel, BorderLayout.NORTH);

    JPanel buttons = new JPanel();
    buttons.setLayout(new BorderLayout());
    JPanel buttonsP = new JPanel();
    buttonsP.setLayout(new GridLayout(1,2));
    
    JPanel ssButs = new JPanel();
    ssButs.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
    ssButs.setLayout(new GridLayout(1, 2, 5, 5));
    ssButs.add(m_StartBut);
    ssButs.add(m_StopBut);
    buttonsP.add(ssButs);
    buttons.add(buttonsP, BorderLayout.SOUTH);
    if (showStoreOutput) {
      buttons.add(m_storeOutput, BorderLayout.NORTH);
    }

    JPanel p3 = new JPanel();
    p3.setBorder(BorderFactory.createTitledBorder("Associator output"));
    p3.setLayout(new BorderLayout());
    final JScrollPane js = new JScrollPane(m_OutText);
    p3.add(js, BorderLayout.CENTER);
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
    
    GridBagLayout gbL = new GridBagLayout();
    GridBagConstraints gbC = new GridBagConstraints();
    JPanel mondo = new JPanel();
    gbL = new GridBagLayout();
    mondo.setLayout(gbL);
    gbC = new GridBagConstraints();
    gbC.anchor = GridBagConstraints.NORTH;
    gbC.fill = GridBagConstraints.HORIZONTAL;
    gbC.gridy = 1;     gbC.gridx = 0;
    gbL.setConstraints(buttons, gbC);
    mondo.add(buttons);
    gbC = new GridBagConstraints();
    gbC.fill = GridBagConstraints.BOTH;
    gbC.gridy = 2;     gbC.gridx = 0; gbC.weightx = 0;
    gbL.setConstraints(m_History, gbC);
    mondo.add(m_History);
    gbC = new GridBagConstraints();
    gbC.fill = GridBagConstraints.BOTH;
    gbC.gridy = 0;     gbC.gridx = 1;
    gbC.gridheight = 3;
    gbC.weightx = 100; gbC.weighty = 100;
    gbL.setConstraints(p3, gbC);
    mondo.add(p3);

    setLayout(new BorderLayout());
    add(p1, BorderLayout.NORTH);
    add(mondo, BorderLayout.CENTER);
  }
  
  /**
   * Sets the Logger to receive informational messages
   *
   * @param newLog the Logger that will now get info messages
   */
  public void setLog(Logger newLog) {

    m_Log = newLog;
  }
  
  /**
   * Tells the panel to use a new set of instances.
   *
   * @param inst a set of Instances
   */
  public void setInstances(Instances inst) {
    
    m_Instances = inst;
    String [] attribNames = new String [m_Instances.numAttributes()];
    for (int i = 0; i < attribNames.length; i++) {
      String type = "";
      switch (m_Instances.attribute(i).type()) {
      case Attribute.NOMINAL:
	type = "(Nom) ";
	break;
      case Attribute.NUMERIC:
	type = "(Num) ";
	break;
      case Attribute.STRING:
	type = "(Str) ";
	break;
      case Attribute.DATE:
	type = "(Dat) ";
	break;
      case Attribute.RELATIONAL:
	type = "(Rel) ";
	break;
      case PreferenceAttribute.RANKING:
    type = "(Rnk) ";
    break;
      default:
	type = "(???) ";
      }
      attribNames[i] = type + m_Instances.attribute(i).name();
    }
    m_StartBut.setEnabled(m_RunThread == null);
    m_StopBut.setEnabled(m_RunThread != null);
  }
  
  /**
   * Starts running the currently configured associator with the current
   * settings. This is run in a separate thread, and will only start if
   * there is no associator already running. The associator output is sent
   * to the results history panel.
   */
  protected void startAssociator() {

    if (m_RunThread == null) {
      m_StartBut.setEnabled(false);
      m_StopBut.setEnabled(true);
      m_RunThread = new Thread() {
	public void run() {
	  // Copy the current state of things
	  m_Log.statusMessage("Setting up...");
	  Instances inst = new Instances(m_Instances);
	  String grph = null;
	  //String xmlRules = null;
	  AssociationRules rulesList = null;
	  Associator associator = (Associator) m_AssociatorEditor.getValue();
	  StringBuffer outBuff = new StringBuffer();
	  String name = (new SimpleDateFormat("HH:mm:ss - "))
	  .format(new Date());
	  String cname = associator.getClass().getName();
	  if (cname.startsWith("weka.associations.")) {
	    name += cname.substring("weka.associations.".length());
	  } else {
	    name += cname;
	  }
          String cmd = m_AssociatorEditor.getValue().getClass().getName();
          if (m_AssociatorEditor.getValue() instanceof OptionHandler)
            cmd += " " + Utils.joinOptions(((OptionHandler) m_AssociatorEditor.getValue()).getOptions());
	  try {

	    // Output some header information
	    m_Log.logMessage("Started " + cname);
	    m_Log.logMessage("Command: " + cmd);
	    if (m_Log instanceof TaskLogger) {
	      ((TaskLogger)m_Log).taskStarted();
	    }
	    outBuff.append("=== Run information ===\n\n");
	    outBuff.append("Scheme:       " + cname);
	    if (associator instanceof OptionHandler) {
	      String [] o = ((OptionHandler) associator).getOptions();
	      outBuff.append(" " + Utils.joinOptions(o));
	    }
	    outBuff.append("\n");
	    outBuff.append("Relation:     " + inst.relationName() + '\n');
	    outBuff.append("Instances:    " + inst.numInstances() + '\n');
	    outBuff.append("Attributes:   " + inst.numAttributes() + '\n');
	    if (inst.numAttributes() < 100) {
	      for (int i = 0; i < inst.numAttributes(); i++) {
		outBuff.append("              " + inst.attribute(i).name()
			       + '\n');
	      }
	    } else {
	      outBuff.append("              [list of attributes omitted]\n");
	    }
	    m_History.addResult(name, outBuff);
	    m_History.setSingle(name);
	    
	    // Build the model and output it.
	    m_Log.statusMessage("Building model on training data...");
	    associator.buildAssociations(inst);
	    outBuff.append("=== Associator model (full training set) ===\n\n");
	    outBuff.append(associator.toString() + '\n');
	    m_History.updateResult(name);
	    if (m_storeOutput.isSelected()) {
	      if (associator instanceof Drawable) {
	        grph = null;
	        try {
	          grph = ((Drawable)associator).graph();
	        } catch (Exception ex) {	        
	        }
	      }

	      if (associator instanceof weka.associations.AssociationRulesProducer) {
	        // xmlRules = null;
	        rulesList = null;
	        try {
	          // xmlRules = ((weka.associations.XMLRulesProducer)associator).xmlRules();
	          rulesList = 
	            ((weka.associations.AssociationRulesProducer)associator).getAssociationRules();
	        } catch (Exception ex) {}
	      }
	    }
	    m_Log.logMessage("Finished " + cname);
	    m_Log.statusMessage("OK");
	  } catch (Exception ex) {
	    m_Log.logMessage(ex.getMessage());
	    m_Log.statusMessage("See error log");
	  } finally {
	    if (grph != null || rulesList != null) {
	      Vector<Object> visVect = new Vector<Object>();

	      if (grph != null) {
	        visVect.add(grph);
	      }

	      if (rulesList != null) {
	        visVect.add(rulesList);
	      }
	      m_History.addObject(name, visVect);
	    }
	    if (isInterrupted()) {
	      m_Log.logMessage("Interrupted " + cname);
	      m_Log.statusMessage("See error log");
	    }
	    m_RunThread = null;
	    m_StartBut.setEnabled(true);
	    m_StopBut.setEnabled(false);
	    if (m_Log instanceof TaskLogger) {
	      ((TaskLogger)m_Log).taskFinished();
	    }
	  }
	}
      };
      m_RunThread.setPriority(Thread.MIN_PRIORITY);
      m_RunThread.start();
    }
  }
  
  /**
   * Stops the currently running Associator (if any).
   */
  protected void stopAssociator() {

    if (m_RunThread != null) {
      m_RunThread.interrupt();
      
      // This is deprecated (and theoretically the interrupt should do).
      m_RunThread.stop();
      
    }
  }

  /**
   * Save the currently selected associator output to a file.
   * @param name the name of the buffer to save
   */
  protected void saveBuffer(String name) {
    StringBuffer sb = m_History.getNamedBuffer(name);
    if (sb != null) {
      if (m_SaveOut.save(sb)) {
	m_Log.logMessage("Save successful.");
      }
    }
  }
  
  /**
   * Pops up a TreeVisualizer for the associator from the currently
   * selected item in the results list
   * @param dottyString the description of the tree in dotty format
   * @param treeName the title to assign to the display
   */
  protected void visualizeTree(String dottyString, String treeName) {
    final javax.swing.JFrame jf = 
      new javax.swing.JFrame("Weka Classifier Tree Visualizer: "+treeName);
    jf.setSize(500,400);
    jf.getContentPane().setLayout(new BorderLayout());
    TreeVisualizer tv = new TreeVisualizer(null,
                                           dottyString,
                                           new PlaceNode2());
    jf.getContentPane().add(tv, BorderLayout.CENTER);
    jf.addWindowListener(new java.awt.event.WindowAdapter() {
        public void windowClosing(java.awt.event.WindowEvent e) {
          jf.dispose();
        }
      });
    
    jf.setVisible(true);
    tv.fitToScreen();
  }
    
  /**
   * Handles constructing a popup menu with visualization options.
   * @param name the name of the result history list entry clicked on by
   * the user
   * @param x the x coordinate for popping up the menu
   * @param y the y coordinate for popping up the menu
   */
  protected void historyRightClickPopup(String name, int x, int y) {
    final String selectedName = name;
    JPopupMenu resultListMenu = new JPopupMenu();
    
    JMenuItem visMainBuffer = new JMenuItem("View in main window");
    if (selectedName != null) {
      visMainBuffer.addActionListener(new ActionListener() {
	  public void actionPerformed(ActionEvent e) {
	    m_History.setSingle(selectedName);
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
	  m_History.openFrame(selectedName);
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
	    saveBuffer(selectedName);
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
	  m_History.removeResult(selectedName);
	}
      });
    } else {
      deleteOutput.setEnabled(false);
    }
    resultListMenu.add(deleteOutput);
    
//    String grph = null;
    Vector<Object> visVect = null;
    if (selectedName != null) {
      //grph = (String)m_History.getNamedObject(selectedName);
      visVect = (Vector)m_History.getNamedObject(selectedName);      
    }
    
/*    final String fgrph = grph;
    JMenuItem visTree = new JMenuItem("Visualize tree");
    if (grph != null) {
      visTree.addActionListener(new ActionListener() {
        public void actionPerformed(ActionEvent e) {
          String title = selectedName;
          visualizeTree(fgrph, title);
        }
      });
      resultListMenu.add(visTree);
    } */
    
    // plugins
    JMenu visPlugins = new JMenu("Plugins");
    boolean availablePlugins = false;
    
    // tree plugins
    if (visVect != null) {
      for (Object o : visVect) {
        if (o instanceof AssociationRules) {
          Vector pluginsVector = 
            GenericObjectEditor.getClassnames(AssociationRuleVisualizePlugin.class.getName());
          for (int i = 0; i < pluginsVector.size(); i++) {
            String className = (String) (pluginsVector.elementAt(i));
            try {
              AssociationRuleVisualizePlugin plugin = 
                (AssociationRuleVisualizePlugin)Class.forName(className).newInstance();
              if (plugin == null) {
                continue;
              }
              availablePlugins = true;
              JMenuItem pluginMenuItem = plugin.getVisualizeMenuItem((AssociationRules)o, selectedName);
              if (pluginMenuItem != null) {
                visPlugins.add(pluginMenuItem);
              }
            } catch (Exception ex) {
              //ex.printStackTrace();
            }
          }
        } else if (o instanceof String) {
          Vector pluginsVector = 
            GenericObjectEditor.getClassnames(TreeVisualizePlugin.class.getName());
          for (int i = 0; i < pluginsVector.size(); i++) {
            String className = (String) (pluginsVector.elementAt(i));
            try {
              TreeVisualizePlugin plugin = (TreeVisualizePlugin) Class.forName(className).newInstance();
              if (plugin == null)
                continue;
              availablePlugins = true;
              JMenuItem pluginMenuItem = plugin.getVisualizeMenuItem((String)o, selectedName);
              // Version version = new Version();
              if (pluginMenuItem != null) {
                /*  if (version.compareTo(plugin.getMinVersion()) < 0)
              pluginMenuItem.setText(pluginMenuItem.getText() + " (weka outdated)");
            if (version.compareTo(plugin.getMaxVersion()) >= 0)
              pluginMenuItem.setText(pluginMenuItem.getText() + " (plugin outdated)");*/
                visPlugins.add(pluginMenuItem);
              }
            }
            catch (Exception e) {
              //e.printStackTrace();
            }
          }
        }
      }
    }
    
    if (availablePlugins) {
      resultListMenu.add(visPlugins);
    }

    resultListMenu.show(m_History.getList(), x, y);
  }
  
  /**
   * updates the capabilities filter of the GOE
   * 
   * @param filter	the new filter to use
   */
  protected void updateCapabilitiesFilter(Capabilities filter) {
    Instances           tempInst;
    Capabilities        filterClass;
    
    if (filter == null) {
      m_AssociatorEditor.setCapabilitiesFilter(new Capabilities(null));
      return;
    }
    
    if (!ExplorerDefaults.getInitGenericObjectEditorFilter())
      tempInst = new Instances(m_Instances, 0);
    else
      tempInst = new Instances(m_Instances);
    tempInst.setClassIndex(-1);
    
    try {
      filterClass = Capabilities.forInstances(tempInst);
    }
    catch (Exception e) {
      filterClass = new Capabilities(null);
    }
    
    m_AssociatorEditor.setCapabilitiesFilter(filterClass);
    
    m_StartBut.setEnabled(true);
    // Check capabilities
    Capabilities currentFilter = m_AssociatorEditor.getCapabilitiesFilter();
    Associator associator = (Associator) m_AssociatorEditor.getValue();
    Capabilities currentSchemeCapabilities =  null;
    if (associator != null && currentFilter != null && 
        (associator instanceof CapabilitiesHandler)) {
      currentSchemeCapabilities = ((CapabilitiesHandler)associator).getCapabilities();
      
      if (!currentSchemeCapabilities.supportsMaybe(currentFilter) &&
          !currentSchemeCapabilities.supports(currentFilter)) {
        m_StartBut.setEnabled(false);
      }
    }
  }
  
  /**
   * method gets called in case of a change event
   * 
   * @param e		the associated change event
   */
  public void capabilitiesFilterChanged(CapabilitiesFilterChangeEvent e) {
    if (e.getFilter() == null)
      updateCapabilitiesFilter(null);
    else
      updateCapabilitiesFilter((Capabilities) e.getFilter().clone());
  }

  /**
   * Sets the Explorer to use as parent frame (used for sending notifications
   * about changes in the data)
   * 
   * @param parent	the parent frame
   */
  public void setExplorer(Explorer parent) {
    m_Explorer = parent;
  }
  
  /**
   * returns the parent Explorer frame
   * 
   * @return		the parent
   */
  public Explorer getExplorer() {
    return m_Explorer;
  }
  
  /**
   * Returns the title for the tab in the Explorer
   * 
   * @return 		the title of this tab
   */
  public String getTabTitle() {
    return "Associate";
  }
  
  /**
   * Returns the tooltip for the tab in the Explorer
   * 
   * @return 		the tooltip of this tab
   */
  public String getTabTitleToolTip() {
    return "Discover association rules";
  }

  /**
   * Tests out the Associator panel from the command line.
   *
   * @param args may optionally contain the name of a dataset to load.
   */
  public static void main(String [] args) {

    try {
      final javax.swing.JFrame jf =
	new javax.swing.JFrame("Weka Explorer: Associator");
      jf.getContentPane().setLayout(new BorderLayout());
      final AssociationsPanel sp = new AssociationsPanel();
      jf.getContentPane().add(sp, BorderLayout.CENTER);
      weka.gui.LogPanel lp = new weka.gui.LogPanel();
      sp.setLog(lp);
      jf.getContentPane().add(lp, BorderLayout.SOUTH);
      jf.addWindowListener(new java.awt.event.WindowAdapter() {
	public void windowClosing(java.awt.event.WindowEvent e) {
	  jf.dispose();
	  System.exit(0);
	}
      });
      jf.pack();
      jf.setVisible(true);
      if (args.length == 1) {
	System.err.println("Loading instances from " + args[0]);
	java.io.Reader r = new java.io.BufferedReader(
			   new java.io.FileReader(args[0]));
	Instances i = new Instances(r);
	sp.setInstances(i);
      }
    } catch (Exception ex) {
      ex.printStackTrace();
      System.err.println(ex.getMessage());
    }
  }
}
