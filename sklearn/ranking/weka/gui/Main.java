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
 * Main.java
 * Copyright (C) 2006 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui;

import weka.classifiers.bayes.net.GUI;
import weka.classifiers.evaluation.ThresholdCurve;
import weka.core.Copyright;
import weka.core.Instances;
import weka.core.Memory;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.SelectedTag;
import weka.core.SystemInfo;
import weka.core.Tag;
import weka.core.Utils;
import weka.core.Version;
import weka.core.scripting.Groovy;
import weka.core.scripting.Jython;
import weka.gui.arffviewer.ArffViewerMainPanel;
import weka.gui.beans.KnowledgeFlowApp;
import weka.gui.beans.StartUpListener;
import weka.gui.boundaryvisualizer.BoundaryVisualizer;
import weka.gui.experiment.Experimenter;
import weka.gui.explorer.Explorer;
import weka.gui.graphvisualizer.GraphVisualizer;
import weka.gui.scripting.GroovyPanel;
import weka.gui.scripting.JythonPanel;
import weka.gui.sql.SqlViewer;
import weka.gui.treevisualizer.Node;
import weka.gui.treevisualizer.NodePlace;
import weka.gui.treevisualizer.PlaceNode2;
import weka.gui.treevisualizer.TreeBuild;
import weka.gui.treevisualizer.TreeVisualizer;
import weka.gui.visualize.PlotData2D;
import weka.gui.visualize.ThresholdVisualizePanel;
import weka.gui.visualize.VisualizePanel;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.GridLayout;
import java.awt.Image;
import java.awt.LayoutManager;
import java.awt.Point;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.Reader;
import java.util.Collections;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Vector;

import javax.swing.BorderFactory;
import javax.swing.ImageIcon;
import javax.swing.JDesktopPane;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JInternalFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JTable;
import javax.swing.SwingConstants;
import javax.swing.WindowConstants;
import javax.swing.event.InternalFrameAdapter;
import javax.swing.event.InternalFrameEvent;

/**
 * Menu-based GUI for Weka, replacement for the GUIChooser.
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -gui &lt;MDI|SDI&gt;
 *  Determines the layout of the GUI:
 *  MDI = MDI Layout
 *  SDI = SDI Layout
 *  (default: MDI)</pre>
 * 
 <!-- options-end -->
 *
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5837 $
 */
public class Main
  extends JFrame
  implements OptionHandler {
  
  /** for serialization. */
  private static final long serialVersionUID = 1453813254824253849L;
  
  /**
   * DesktopPane with background image.
   * 
   * @author  fracpete (fracpete at waikato dot ac dot nz)
   * @version $Revision: 5837 $
   */
  public static class BackgroundDesktopPane
    extends JDesktopPane {
    
    /** for serialization. */
    private static final long serialVersionUID = 2046713123452402745L;
    
    /** the actual background image. */
    protected Image m_Background;
    
    /**
     * intializes the desktop pane.
     * 
     * @param image	the image to use as background
     */
    public BackgroundDesktopPane(String image) {
      super();
      
      try {
	m_Background = Toolkit.getDefaultToolkit().getImage(ClassLoader.getSystemResource(image));
      }
      catch (Exception e) {
	e.printStackTrace();
      }
    }
    
    /**
     * draws the background image.
     * 
     * @param g		the graphics context
     */
    public void paintComponent(Graphics g) {
      super.paintComponent(g);
     
      if (m_Background != null) {
	g.setColor(Color.WHITE);
	g.clearRect(0, 0, getWidth(), getHeight());
	
	int width  = m_Background.getWidth(null);
	int height = m_Background.getHeight(null);
	int x = (getWidth() - width) / 2;
	int y = (getHeight() - height) / 2;
	g.drawImage(m_Background, x, y, width, height, this);
      }
    }
  }
  
  /**
   * Specialized JFrame class.
   * 
   * @author  fracpete (fracpete at waikato dot ac dot nz)
   * @version $Revision: 5837 $
   */
  public static class ChildFrameSDI 
    extends JFrame {
    
    /** for serialization. */
    private static final long serialVersionUID = 8588293938686425618L;
    
    /** the parent frame. */
    protected Main m_Parent;
    
    /**
     * constructs a new internal frame that knows about its parent.
     * 
     * @param parent	the parent frame
     * @param title	the title of the frame
     */
    public ChildFrameSDI(Main parent, String title) {
      super(title);
      
      m_Parent = parent;

      addWindowListener(new WindowAdapter() {
	public void windowActivated(WindowEvent e) {
	  // update title of parent
	  if (getParentFrame() != null)
	    getParentFrame().createTitle(getTitle());
	}
      });
      
      // add to parent
      if (getParentFrame() != null) {
	getParentFrame().addChildFrame(this);
	setIconImage(getParentFrame().getIconImage());
      }
    }
    
    /**
     * returns the parent frame, can be null.
     * 
     * @return		the parent frame
     */
    public Main getParentFrame() {
      return m_Parent;
    }
    
    /**
     * de-registers the child frame with the parent first.
     */
    public void dispose() {
      if (getParentFrame() != null) {
	getParentFrame().removeChildFrame(this);
	getParentFrame().createTitle("");
      }
      
      super.dispose();
    }
  }
  
  /**
   * Specialized JInternalFrame class.
   * 
   * @author  fracpete (fracpete at waikato dot ac dot nz)
   * @version $Revision: 5837 $
   */
  public static class ChildFrameMDI
    extends JInternalFrame {
    
    /** for serialization. */
    private static final long serialVersionUID = 3772573515346899959L;
    
    /** the parent frame. */
    protected Main m_Parent;
    
    /**
     * constructs a new internal frame that knows about its parent.
     * 
     * @param parent	the parent frame
     * @param title	the title of the frame
     */
    public ChildFrameMDI(Main parent, String title) {
      super(title, true, true, true, true);
      
      m_Parent = parent;

      addInternalFrameListener(new InternalFrameAdapter() {
	public void internalFrameActivated(InternalFrameEvent e) {
	  // update title of parent
	  if (getParentFrame() != null)
	    getParentFrame().createTitle(getTitle());
	}
      });
      
      // add to parent
      if (getParentFrame() != null) {
	getParentFrame().addChildFrame(this);
	getParentFrame().jDesktopPane.add(this);
      }
    }
    
    /**
     * returns the parent frame, can be null.
     * 
     * @return		the parent frame
     */
    public Main getParentFrame() {
      return m_Parent;
    }
    
    /**
     * de-registers the child frame with the parent first.
     */
    public void dispose() {
      if (getParentFrame() != null) {
	getParentFrame().removeChildFrame(this);
	getParentFrame().createTitle("");
      }
      
      super.dispose();
    }
  }

  /** displays the GUI as MDI. */
  public final static int GUI_MDI = 0;
  /** displays the GUI as SDI. */
  public final static int GUI_SDI = 1;
  /** GUI tags. */
  public static final Tag[] TAGS_GUI = {
    new Tag(GUI_MDI, "MDI", "MDI Layout"),
    new Tag(GUI_SDI, "SDI", "SDI Layout")
  };
  
  /** the frame itself. */
  protected Main m_Self;
  
  /** the type of GUI to display. */
  protected int m_GUIType = GUI_MDI;
  
  /** variable for the Main class which would be set to null by the memory
   *  monitoring thread to free up some memory if we running out of memory. */
  protected static Main m_MainCommandline;
  
  /** singleton instance of the GUI. */
  protected static Main m_MainSingleton;

  /** list of things to be notified when the startup process of
   *  the KnowledgeFlow is complete. */
  protected static Vector m_StartupListeners = new Vector();

  /** for monitoring the Memory consumption. */
  protected static Memory m_Memory = new Memory(true);
  
  /** contains the child frames (title &lt;-&gt; object). */
  protected HashSet<Container> m_ChildFrames = new HashSet<Container>();

  /** The frame of the LogWindow. */
  protected static LogWindow m_LogWindow = new LogWindow();

  /** filechooser for the TreeVisualizer. */
  protected JFileChooser m_FileChooserTreeVisualizer = new JFileChooser(new File(System.getProperty("user.dir")));

  /** filechooser for the GraphVisualizer. */
  protected JFileChooser m_FileChooserGraphVisualizer = new JFileChooser(new File(System.getProperty("user.dir")));

  /** filechooser for Plots. */
  protected JFileChooser m_FileChooserPlot = new JFileChooser(new File(System.getProperty("user.dir")));

  /** filechooser for ROC curves. */
  protected JFileChooser m_FileChooserROC = new JFileChooser(new File(System.getProperty("user.dir")));
  
  // GUI components
  private JMenu jMenuHelp;
  private JMenu jMenuVisualization;
  private JMenu jMenuTools;
  private JDesktopPane jDesktopPane;
  private JMenu jMenuApplications;
  private JMenuItem jMenuItemHelpSystemInfo;
  private JMenuItem jMenuItemHelpAbout;
  private JMenuItem jMenuItemHelpHomepage;
  private JMenuItem jMenuItemHelpWekaWiki;
  private JMenuItem jMenuItemHelpWekaDoc;
  private JMenuItem jMenuItemHelpSourceforge;
  private JMenuItem jMenuItemVisualizationBoundaryVisualizer;
  private JMenuItem jMenuItemVisualizationGraphVisualizer;
  private JMenuItem jMenuItemVisualizationTreeVisualizer;
  private JMenuItem jMenuItemVisualizationROC;
  private JMenuItem jMenuItemVisualizationPlot;
  private JMenuItem jMenuItemToolsEnsembleLibrary;
  private JMenuItem jMenuItemToolsSqlViewer;
  private JMenuItem jMenuItemToolsGroovyConsole;
  private JMenuItem jMenuItemToolsJythonConsole;
  private JMenuItem jMenuItemToolsArffViewer;
  private JMenuItem jMenuItemApplicationsSimpleCLI;
  private JMenuItem jMenuItemApplicationsKnowledgeFlow;
  private JMenuItem jMenuItemApplicationsExperimenter;
  private JMenuItem jMenuItemApplicationsExplorer;
  private JMenuItem jMenuItemProgramExit;
  private JMenuItem jMenuItemProgramLogWindow;
  private JMenuItem jMenuItemProgramMemoryUsage;
  private JMenuItem jMenuItemProgramPreferences;  // TODO: see below
  private JMenu jMenuProgram;
  private JMenu jMenuExtensions;
  private JMenu jMenuWindows;
  private JMenuBar jMenuBar;
  
  /**
   * default constructor.
   */
  public Main() {
    super();
  }
  
  /**
   * creates a frame (depending on m_GUIType) and returns it.
   * 
   * @param parent		the parent of the generated frame
   * @param title		the title of the frame
   * @param c			the component to place, can be null
   * @param layout		the layout to use, e.g., BorderLayout
   * @param layoutConstraints	the layout constraints, e.g., BorderLayout.CENTER
   * @param width		the width of the frame, ignored if -1
   * @param height		the height of the frame, ignored if -1
   * @param menu		an optional menu
   * @param listener		if true a default listener is added
   * @param visible		if true then the frame is made visible immediately
   * @return			the generated frame
   * @see 			#m_GUIType
   */
  protected Container createFrame(
      Main parent, String title, Component c, LayoutManager layout, 
      Object layoutConstraints, int width, int height, JMenuBar menu,
      boolean listener, boolean visible) {

    Container result = null;
    
    if (m_GUIType == GUI_MDI) {
      final ChildFrameMDI frame = new ChildFrameMDI(parent, title);
      
      // layout
      frame.setLayout(layout);
      if (c != null)
	frame.getContentPane().add(c, layoutConstraints);
      
      // menu
      frame.setJMenuBar(menu);
      
      // size
      frame.pack();
      if ((width > -1) && (height > -1))
	frame.setSize(width, height);
      frame.validate();

      // listener?
      if (listener) {
	frame.addInternalFrameListener(new InternalFrameAdapter() {
	  public void internalFrameClosing(InternalFrameEvent e) {
	    frame.dispose();
	  }
	});
      }
      
      // display frame
      if (visible) {
	frame.setVisible(true);
	try {
	  frame.setSelected(true);
	}
	catch (Exception e) {
	  e.printStackTrace();
	}
      }
      
      result = frame;
    }
    else if (m_GUIType == GUI_SDI) {
      final ChildFrameSDI frame = new ChildFrameSDI(parent, title);
      
      // layout
      frame.setLayout(layout);
      if (c != null)
	frame.getContentPane().add(c, layoutConstraints);
      
      // menu
      frame.setJMenuBar(menu);
      
      // size
      frame.pack();
      if ((width > -1) && (height > -1))
	frame.setSize(width, height);
      frame.validate();

      // location
      int screenHeight = getGraphicsConfiguration().getBounds().height;
      int screenWidth  = getGraphicsConfiguration().getBounds().width;
      frame.setLocation(
	  (screenWidth - frame.getBounds().width) / 2,
	  (screenHeight - frame.getBounds().height) / 2);
      
      // listener?
      if (listener) {
	frame.addWindowListener(new WindowAdapter() {
	  public void windowClosing(WindowEvent e) {
	    frame.dispose();
	  }
	});
      }
      
      // display frame
      if (visible)
	frame.setVisible(true);

      result = frame;
    }
    
    return result;
  }
  
  /**
   * insert the menu item in a sorted fashion.
   * 
   * @param menu	the menu to add the item to
   * @param menuitem	the menu item to add
   */
  protected void insertMenuItem(JMenu menu, JMenuItem menuitem) {
    insertMenuItem(menu, menuitem, 0);
  }
  
  /**
   * insert the menu item in a sorted fashion.
   * 
   * @param menu	the menu to add the item to
   * @param menuitem	the menu item to add
   * @param startIndex	the index in the menu to start with (0-based)
   */
  protected void insertMenuItem(JMenu menu, JMenuItem menuitem, int startIndex) {
    boolean	inserted;
    int		i;
    JMenuItem	current;
    String	currentStr;
    String	newStr;
    
    inserted = false;
    newStr   = menuitem.getText().toLowerCase();
    
    // try to find a spot inbetween
    for (i = startIndex; i < menu.getMenuComponentCount(); i++) {
      if (!(menu.getMenuComponent(i) instanceof JMenuItem))
	continue;
      
      current    = (JMenuItem) menu.getMenuComponent(i);
      currentStr = current.getText().toLowerCase();
      if (currentStr.compareTo(newStr) > 0) {
	inserted = true;
	menu.insert(menuitem, i);
	break;
      }
    }
    
    // add it at the end if not yet inserted
    if (!inserted)
      menu.add(menuitem);
  }
  
  /**
   * initializes the GUI.
   */
  protected void initGUI() {
    m_Self = this;
    
    try {
      // main window
      createTitle("");
      this.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
      this.setIconImage(new ImageIcon(getClass().getClassLoader().getResource("weka/gui/weka_icon.gif")).getImage());

      // bits and pieces
      m_FileChooserGraphVisualizer.addChoosableFileFilter(
	  new ExtensionFileFilter(".bif", "BIF Files (*.bif)"));
      m_FileChooserGraphVisualizer.addChoosableFileFilter(
	  new ExtensionFileFilter(".xml", "XML Files (*.xml)"));

      m_FileChooserPlot.addChoosableFileFilter(
	  new ExtensionFileFilter(
	      Instances.FILE_EXTENSION,
	      "ARFF Files (*" + Instances.FILE_EXTENSION + ")"));
      m_FileChooserPlot.setMultiSelectionEnabled(true);

      m_FileChooserROC.addChoosableFileFilter(
	  new ExtensionFileFilter(
	      Instances.FILE_EXTENSION,
	      "ARFF Files (*" + Instances.FILE_EXTENSION + ")"));

      // Desktop
      if (m_GUIType == GUI_MDI) {
	jDesktopPane = new BackgroundDesktopPane("weka/gui/images/weka_background.gif");
	jDesktopPane.setDragMode(JDesktopPane.OUTLINE_DRAG_MODE);
	setContentPane(jDesktopPane);
      }
      else {
	jDesktopPane = null;
      }

      // Menu
      jMenuBar = new JMenuBar();
      setJMenuBar(jMenuBar);

      // Program
      jMenuProgram = new JMenu();
      jMenuBar.add(jMenuProgram);
      jMenuProgram.setText("Program");
      jMenuProgram.setMnemonic('P');

      // Program/Preferences
      // TODO: read all properties from all props file and display them
      /*
        jMenuItemProgramPreferences = new JMenuItem();
        jMenuProgram.add(jMenuItemProgramPreferences);
        jMenuItemProgramPreferences.setText("Preferences");
        jMenuItemProgramPreferences.setMnemonic('P');
        jMenuItemProgramPreferences.addActionListener(new ActionListener() {
  	  public void actionPerformed(ActionEvent evt) {
	    System.out.println("jMenuItemProgramPreferences.actionPerformed, event="+evt);
	    //TODO add your code for jMenuItemProgramPreferences.actionPerformed
  	  }
        });
       */

      // Program/LogWindow
      jMenuItemProgramLogWindow = new JMenuItem();
      jMenuProgram.add(jMenuItemProgramLogWindow);
      jMenuItemProgramLogWindow.setText("LogWindow");
      jMenuItemProgramLogWindow.setMnemonic('L');
      jMenuItemProgramLogWindow.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent evt) {
	  m_LogWindow.setVisible(true);
	}
      });
      
      jMenuItemProgramMemoryUsage = new JMenuItem();
      jMenuProgram.add(jMenuItemProgramMemoryUsage);
      jMenuItemProgramMemoryUsage.setText("Memory usage");
      jMenuItemProgramMemoryUsage.setMnemonic('M');
      jMenuItemProgramMemoryUsage.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent evt) {
	  String title = jMenuItemProgramMemoryUsage.getText();
	  if (!containsWindow(title)) {
	    final MemoryUsagePanel panel = new MemoryUsagePanel();
	    Container c = createFrame(
		m_Self, title, panel, new BorderLayout(), 
		BorderLayout.CENTER, 400, 50, null, true, true);
	    
	    // optimize size
	    Dimension size = c.getPreferredSize();
	    c.setSize(new Dimension((int) size.getWidth(), (int) size.getHeight()));

	    // stop threads
	    if (m_GUIType == GUI_MDI) {
	      final ChildFrameMDI frame = (ChildFrameMDI) c;
	      Point l = panel.getFrameLocation();
	      if ((l.x != -1) && (l.y != -1))
		frame.setLocation(l);
	      frame.addInternalFrameListener(new InternalFrameAdapter() {
		public void internalFrameClosing(InternalFrameEvent e) {
		  panel.stopMonitoring();
		}
	      });
	    }
	    else {
	      final ChildFrameSDI frame = (ChildFrameSDI) c;
	      Point l = panel.getFrameLocation();
	      if ((l.x != -1) && (l.y != -1))
		frame.setLocation(l);
	      frame.addWindowListener(new WindowAdapter() {
		public void windowClosing(WindowEvent e) {
		  panel.stopMonitoring();
		}
	      });
	    }
	  }
	  else {
	    showWindow(getWindow(title));
	  }
	}
      });

      jMenuProgram.add(new JSeparator());

      // Program/Exit
      jMenuItemProgramExit = new JMenuItem();
      jMenuProgram.add(jMenuItemProgramExit);
      jMenuItemProgramExit.setText("Exit");
      jMenuItemProgramExit.setMnemonic('E');
      jMenuItemProgramExit.addActionListener(new ActionListener() {	
	public void actionPerformed(ActionEvent evt) {
	  // close all children
	  Iterator iter = getWindowList();
	  Vector<Container> list = new Vector<Container>();
	  while (iter.hasNext())
	    list.add((Container) iter.next());
	  for (int i = 0; i < list.size(); i++) {
	    Container c = list.get(i);
	    if (c instanceof ChildFrameMDI)
	      ((ChildFrameMDI) c).dispose();
	    else if (c instanceof ChildFrameSDI)
	      ((ChildFrameSDI) c).dispose();
	  }
	  // close logwindow
	  m_LogWindow.dispose();
	  // close main window
	  m_Self.dispose();
	  // make sure we stop
	  System.exit(0);
	}
      });

      // Applications
      jMenuApplications = new JMenu();
      jMenuBar.add(jMenuApplications);
      jMenuApplications.setText("Applications");
      jMenuApplications.setMnemonic('A');

      // Applications/Explorer
      jMenuItemApplicationsExplorer = new JMenuItem();
      jMenuApplications.add(jMenuItemApplicationsExplorer);
      jMenuItemApplicationsExplorer.setText("Explorer");
      jMenuItemApplicationsExplorer.setMnemonic('E');
      jMenuItemApplicationsExplorer.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent evt) {
	  String title = jMenuItemApplicationsExplorer.getText();
	  if (!containsWindow(title)) {
	    createFrame(
		m_Self, title, new Explorer(), new BorderLayout(), 
		BorderLayout.CENTER, 800, 600, null, true, true);
	  }
	  else {
	    showWindow(getWindow(title));
	  }
	}
      });

      // Applications/Experimenter
      jMenuItemApplicationsExperimenter = new JMenuItem();
      jMenuApplications.add(jMenuItemApplicationsExperimenter);
      jMenuItemApplicationsExperimenter.setText("Experimenter");
      jMenuItemApplicationsExperimenter.setMnemonic('X');
      jMenuItemApplicationsExperimenter.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent evt) {
	  String title = jMenuItemApplicationsExperimenter.getText();
	  if (!containsWindow(title)) {
	    createFrame(
		m_Self, title, new Experimenter(false), new BorderLayout(), 
		BorderLayout.CENTER, 800, 600, null, true, true);
	  }
	  else {
	    showWindow(getWindow(title));
	  }
	}
      });

      // Applications/KnowledgeFlow
      jMenuItemApplicationsKnowledgeFlow = new JMenuItem();
      jMenuApplications.add(jMenuItemApplicationsKnowledgeFlow);
      jMenuItemApplicationsKnowledgeFlow.setText("KnowledgeFlow");
      jMenuItemApplicationsKnowledgeFlow.setMnemonic('K');
      jMenuItemApplicationsKnowledgeFlow.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent evt) {
	  String title = jMenuItemApplicationsKnowledgeFlow.getText();
	  if (!containsWindow(title)) {
	    KnowledgeFlowApp.createSingleton(new String[0]);
	    createFrame(
		m_Self, title, KnowledgeFlowApp.getSingleton(), new BorderLayout(), 
		BorderLayout.CENTER, 900, 600, null, true, true);
	  }
	  else {
	    showWindow(getWindow(title));
	  }
	}
      });

      // Applications/SimpleCLI
      jMenuItemApplicationsSimpleCLI = new JMenuItem();
      jMenuApplications.add(jMenuItemApplicationsSimpleCLI);
      jMenuItemApplicationsSimpleCLI.setText("SimpleCLI");
      jMenuItemApplicationsSimpleCLI.setMnemonic('S');
      jMenuItemApplicationsSimpleCLI.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent evt) {
	  String title = jMenuItemApplicationsSimpleCLI.getText();
	  if (!containsWindow(title)) {
	    try {
	      createFrame(
		  m_Self, title, new SimpleCLIPanel(), new BorderLayout(), 
		  BorderLayout.CENTER, 600, 500, null, true, true);
	    }
	    catch (Exception e) {
	      e.printStackTrace();
	      JOptionPane.showMessageDialog(
		  m_Self, "Error instantiating SimpleCLI:\n" + e.getMessage());
	      return;
	    }
	  }
	  else {
	    showWindow(getWindow(title));
	  }
	}
      });

      // Tools
      jMenuTools = new JMenu();
      jMenuBar.add(jMenuTools);
      jMenuTools.setText("Tools");
      jMenuTools.setMnemonic('T');

      // Tools/ArffViewer
      jMenuItemToolsArffViewer = new JMenuItem();
      jMenuTools.add(jMenuItemToolsArffViewer);
      jMenuItemToolsArffViewer.setText("ArffViewer");
      jMenuItemToolsArffViewer.setMnemonic('A');
      jMenuItemToolsArffViewer.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent evt) {
	  String title = jMenuItemToolsArffViewer.getText();
	  if (!containsWindow(title)) {
	    ArffViewerMainPanel panel = new ArffViewerMainPanel(null);
	    panel.setConfirmExit(false);
	    Container frame = createFrame(
		m_Self, title, panel, new BorderLayout(), 
		BorderLayout.CENTER, 800, 600, panel.getMenu(), true, true);
	    panel.setParent(frame);
	  }
	  else {
	    showWindow(getWindow(title));
	  }
	}
      });

      // Tools/SqlViewer
      jMenuItemToolsSqlViewer = new JMenuItem();
      jMenuTools.add(jMenuItemToolsSqlViewer);
      jMenuItemToolsSqlViewer.setText("SqlViewer");
      jMenuItemToolsSqlViewer.setMnemonic('S');
      jMenuItemToolsSqlViewer.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent evt) {
	  String title = jMenuItemToolsSqlViewer.getText();
	  if (!containsWindow(title)) {
	    final SqlViewer sql = new SqlViewer(null);
	    final Container frame = createFrame(
		m_Self, title, sql, new BorderLayout(), 
		BorderLayout.CENTER, -1, -1, null, false, true);

	    // custom listener
	    if (frame instanceof ChildFrameMDI) {
	      ((ChildFrameMDI) frame).addInternalFrameListener(new InternalFrameAdapter() {
		public void internalFrameClosing(InternalFrameEvent e) {
		  sql.saveSize();
		  ((ChildFrameMDI) frame).dispose();
		}
	      });
	    }
	    else if (frame instanceof ChildFrameSDI) {
	      ((ChildFrameSDI) frame).addWindowListener(new WindowAdapter() {
		public void windowClosing(WindowEvent e) {
		  sql.saveSize();
		  ((ChildFrameSDI) frame).dispose();
		}
	      });
	    }
	  }
	  else {
	    showWindow(getWindow(title));
	  }
	}
      });
      
      // Tools/Bayes net editor
      // Tools/Bayes net editor
      final JMenuItem jMenuItemBayesNet = new JMenuItem();
      jMenuTools.add(jMenuItemBayesNet);
      jMenuItemBayesNet.setText("Bayes net editor");
      jMenuItemBayesNet.setMnemonic('N');

      jMenuItemBayesNet.addActionListener(new ActionListener() {
        public void actionPerformed(ActionEvent e) {
          String title = jMenuItemBayesNet.getText();
          
          if (!containsWindow(title)) {
            final GUI bayesNetGUI = new GUI();
            final Container frame = createFrame(
                m_Self, title, bayesNetGUI, new BorderLayout(), 
                BorderLayout.CENTER, 800, 600, bayesNetGUI.getMenuBar(), false, true);
          }
          else {
            showWindow(getWindow(title));
          }          
        }       
      });

      // Tools/Groovy console
      if (Groovy.isPresent()) {
	jMenuItemToolsGroovyConsole = new JMenuItem();
	jMenuTools.add(jMenuItemToolsGroovyConsole);
	jMenuItemToolsGroovyConsole.setText("Groovy console");
	jMenuItemToolsGroovyConsole.setMnemonic('G');
	jMenuItemToolsGroovyConsole.addActionListener(new ActionListener() {
	  public void actionPerformed(ActionEvent evt) {
	    String title = jMenuItemToolsGroovyConsole.getText();
	    if (!containsWindow(title)) {
	      final GroovyPanel panel = new GroovyPanel();
	      final Container frame = createFrame(
		  m_Self, title, panel, new BorderLayout(), 
		  BorderLayout.CENTER, 800, 600, panel.getMenuBar(), false, true);

	      // custom listener
	      if (frame instanceof ChildFrameMDI) {
		((ChildFrameMDI) frame).addInternalFrameListener(new InternalFrameAdapter() {
		  public void internalFrameClosing(InternalFrameEvent e) {
		    ((ChildFrameMDI) frame).dispose();
		  }
		});
	      }
	      else if (frame instanceof ChildFrameSDI) {
		((ChildFrameSDI) frame).addWindowListener(new WindowAdapter() {
		  public void windowClosing(WindowEvent e) {
		    ((ChildFrameSDI) frame).dispose();
		  }
		});
	      }
	    }
	    else {
	      showWindow(getWindow(title));
	    }
	  }
	});
      }

      // Tools/Jython console
      if (Jython.isPresent()) {
	jMenuItemToolsJythonConsole = new JMenuItem();
	jMenuTools.add(jMenuItemToolsJythonConsole);
	jMenuItemToolsJythonConsole.setText("Jython console");
	jMenuItemToolsJythonConsole.setMnemonic('J');
	jMenuItemToolsJythonConsole.addActionListener(new ActionListener() {
	  public void actionPerformed(ActionEvent evt) {
	    String title = jMenuItemToolsJythonConsole.getText();
	    if (!containsWindow(title)) {
	      final JythonPanel panel = new JythonPanel();
	      final Container frame = createFrame(
		  m_Self, title, panel, new BorderLayout(), 
		  BorderLayout.CENTER, 800, 600, panel.getMenuBar(), false, true);

	      // custom listener
	      if (frame instanceof ChildFrameMDI) {
		((ChildFrameMDI) frame).addInternalFrameListener(new InternalFrameAdapter() {
		  public void internalFrameClosing(InternalFrameEvent e) {
		    ((ChildFrameMDI) frame).dispose();
		  }
		});
	      }
	      else if (frame instanceof ChildFrameSDI) {
		((ChildFrameSDI) frame).addWindowListener(new WindowAdapter() {
		  public void windowClosing(WindowEvent e) {
		    ((ChildFrameSDI) frame).dispose();
		  }
		});
	      }
	    }
	    else {
	      showWindow(getWindow(title));
	    }
	  }
      });
      }

      // Tools/EnsembleLibrary
      /* currently disabled due to bugs... FracPete
      jMenuItemToolsEnsembleLibrary = new JMenuItem();
      jMenuTools.add(jMenuItemToolsEnsembleLibrary);
      jMenuItemToolsEnsembleLibrary.setText("EnsembleLibrary");
      jMenuItemToolsEnsembleLibrary.setMnemonic('E');
      jMenuItemToolsEnsembleLibrary.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent evt) {
	  String title = jMenuItemToolsEnsembleLibrary.getText();
	  if (!containsWindow(title)) {
	    EnsembleLibrary value = new EnsembleLibrary();
	    EnsembleLibraryEditor libraryEditor = new EnsembleLibraryEditor();
	    libraryEditor.setValue(value);
	    createFrame(
		m_Self, title, libraryEditor.getCustomEditor(), new BorderLayout(), 
		BorderLayout.CENTER, 800, 600, null, true, true);
	  }
	  else {
	    showWindow(getWindow(title));
	  }
	}
      });
      */

      // Visualization
      jMenuVisualization = new JMenu();
      jMenuBar.add(jMenuVisualization);
      jMenuVisualization.setText("Visualization");
      jMenuVisualization.setMnemonic('V');

      // Visualization/Plot
      jMenuItemVisualizationPlot = new JMenuItem();
      jMenuVisualization.add(jMenuItemVisualizationPlot);
      jMenuItemVisualizationPlot.setText("Plot");
      jMenuItemVisualizationPlot.setMnemonic('P');
      jMenuItemVisualizationPlot.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent evt) {
	  // choose file
	  int retVal = m_FileChooserPlot.showOpenDialog(m_Self);
	  if (retVal != JFileChooser.APPROVE_OPTION)
	    return;

	  // build plot
	  VisualizePanel panel = new VisualizePanel();
	  String filenames = "";
	  File[] files = m_FileChooserPlot.getSelectedFiles();
	  for (int j = 0; j < files.length; j++) {
	    String filename = files[j].getAbsolutePath();
	    if (j > 0)
	      filenames += ", ";
	    filenames += filename;
	    System.err.println("Loading instances from " + filename);
	    try {
	      Reader r = new java.io.BufferedReader(new FileReader(filename));
	      Instances i = new Instances(r);
	      i.setClassIndex(i.numAttributes()-1);
	      PlotData2D pd1 = new PlotData2D(i);

	      if (j == 0) {
		pd1.setPlotName("Master plot");
		panel.setMasterPlot(pd1);
	      } else {
		pd1.setPlotName("Plot "+(j+1));
		pd1.m_useCustomColour = true;
		pd1.m_customColour = (j % 2 == 0) ? Color.red : Color.blue; 
		panel.addPlot(pd1);
	      }
	    }
	    catch (Exception e) {
	      e.printStackTrace();
	      JOptionPane.showMessageDialog(
		  m_Self, "Error loading file '" + files[j] + "':\n" + e.getMessage());
	      return;
	    }
	  }

	  // create frame
	  createFrame(
	      m_Self, jMenuItemVisualizationPlot.getText() + " - " + filenames, 
	      panel, new BorderLayout(), 
	      BorderLayout.CENTER, 800, 600, null, true, true);
	}
      });

      // Visualization/ROC
      // based on this Wiki article:
      // http://weka.sourceforge.net/wiki/index.php/Visualizing_ROC_curve
      jMenuItemVisualizationROC = new JMenuItem();
      jMenuVisualization.add(jMenuItemVisualizationROC);
      jMenuItemVisualizationROC.setText("ROC");
      jMenuItemVisualizationROC.setMnemonic('R');
      jMenuItemVisualizationROC.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent evt) {
	  // choose file
	  int retVal = m_FileChooserROC.showOpenDialog(m_Self);
	  if (retVal != JFileChooser.APPROVE_OPTION)
	    return;

	  // create plot
	  String filename  = m_FileChooserROC.getSelectedFile().getAbsolutePath();
	  Instances result = null;
	  try {
	    result = new Instances(new BufferedReader(new FileReader(filename)));
	  }
	  catch (Exception e) {
	    e.printStackTrace();
	    JOptionPane.showMessageDialog(
		m_Self, "Error loading file '" + filename + "':\n" + e.getMessage());
	    return;
	  }
	  result.setClassIndex(result.numAttributes() - 1);
	  ThresholdVisualizePanel vmc = new ThresholdVisualizePanel();
	  vmc.setROCString("(Area under ROC = " + 
	      Utils.doubleToString(ThresholdCurve.getROCArea(result), 4) + ")");
	  vmc.setName(result.relationName());
	  PlotData2D tempd = new PlotData2D(result);
	  tempd.setPlotName(result.relationName());
	  tempd.addInstanceNumberAttribute();
	  try {
	    vmc.addPlot(tempd);
	  }
	  catch (Exception e) {
	    e.printStackTrace();
	    JOptionPane.showMessageDialog(
		m_Self, "Error adding plot:\n" + e.getMessage());
	    return;
	  }

	  createFrame(
	      m_Self, jMenuItemVisualizationROC.getText() + " - " + filename, 
	      vmc, new BorderLayout(), 
	      BorderLayout.CENTER, 800, 600, null, true, true);
	}
      });

      // Visualization/TreeVisualizer
      jMenuItemVisualizationTreeVisualizer = new JMenuItem();
      jMenuVisualization.add(jMenuItemVisualizationTreeVisualizer);
      jMenuItemVisualizationTreeVisualizer.setText("TreeVisualizer");
      jMenuItemVisualizationTreeVisualizer.setMnemonic('T');
      jMenuItemVisualizationTreeVisualizer.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent evt) {
	  // choose file
	  int retVal = m_FileChooserTreeVisualizer.showOpenDialog(m_Self);
	  if (retVal != JFileChooser.APPROVE_OPTION)
	    return;

	  // build tree
	  String filename = m_FileChooserTreeVisualizer.getSelectedFile().getAbsolutePath();
	  TreeBuild builder = new TreeBuild();
	  Node top = null;
	  NodePlace arrange = new PlaceNode2();
	  try {
	    top = builder.create(new FileReader(filename));
	  }
	  catch (Exception e) {
	    e.printStackTrace();
	    JOptionPane.showMessageDialog(
		m_Self, "Error loading file '" + filename + "':\n" + e.getMessage());
	    return;
	  }

	  // create frame
	  createFrame(
	      m_Self, jMenuItemVisualizationTreeVisualizer.getText() + " - " + filename, 
	      new TreeVisualizer(null, top, arrange), new BorderLayout(), 
	      BorderLayout.CENTER, 800, 600, null, true, true);
	}
      });

      // Visualization/GraphVisualizer
      jMenuItemVisualizationGraphVisualizer = new JMenuItem();
      jMenuVisualization.add(jMenuItemVisualizationGraphVisualizer);
      jMenuItemVisualizationGraphVisualizer.setText("GraphVisualizer");
      jMenuItemVisualizationGraphVisualizer.setMnemonic('G');
      jMenuItemVisualizationGraphVisualizer.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent evt) {
	  // choose file
	  int retVal = m_FileChooserGraphVisualizer.showOpenDialog(m_Self);
	  if (retVal != JFileChooser.APPROVE_OPTION)
	    return;

	  // build graph
	  String filename = m_FileChooserGraphVisualizer.getSelectedFile().getAbsolutePath();
	  GraphVisualizer panel = new GraphVisualizer();
	  try{
	    if (    filename.toLowerCase().endsWith(".xml") 
		|| filename.toLowerCase().endsWith(".bif") ) {
	      panel.readBIF(new FileInputStream(filename));
	    }
	    else {
	      panel.readDOT(new FileReader(filename));
	    }
	  }
	  catch (Exception e) {
	    e.printStackTrace();
	    JOptionPane.showMessageDialog(
		m_Self, "Error loading file '" + filename + "':\n" + e.getMessage());
	    return;
	  }

	  // create frame
	  createFrame(
	      m_Self, jMenuItemVisualizationGraphVisualizer.getText() + " - " + filename, 
	      panel, new BorderLayout(), 
	      BorderLayout.CENTER, 800, 600, null, true, true);
	}
      });

      // Visualization/BoundaryVisualizer
      jMenuItemVisualizationBoundaryVisualizer = new JMenuItem();
      jMenuVisualization.add(jMenuItemVisualizationBoundaryVisualizer);
      jMenuItemVisualizationBoundaryVisualizer.setText("BoundaryVisualizer");
      jMenuItemVisualizationBoundaryVisualizer.setMnemonic('B');
      jMenuItemVisualizationBoundaryVisualizer.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent evt) {
	  String title = jMenuItemVisualizationBoundaryVisualizer.getText();
	  if (!containsWindow(title)) {
	    createFrame(
		m_Self, title, new BoundaryVisualizer(), new BorderLayout(), 
		BorderLayout.CENTER, 800, 600, null, true, true);
	    // dont' do a System.exit after last window got closed!
	    BoundaryVisualizer.setExitIfNoWindowsOpen(false);
	  }
	  else {
	    showWindow(getWindow(title));
	  }
	}
      });

      // Extensions
      jMenuExtensions = new JMenu("Extensions");
      jMenuExtensions.setMnemonic(java.awt.event.KeyEvent.VK_E);
      jMenuBar.add(jMenuExtensions);
      jMenuExtensions.setVisible(false);

      String extensions = GenericObjectEditor.EDITOR_PROPERTIES.getProperty(
	  MainMenuExtension.class.getName(), "");

      if (extensions.length() > 0) {
	jMenuExtensions.setVisible(true);
	String[] classnames = GenericObjectEditor.EDITOR_PROPERTIES.getProperty(
	    MainMenuExtension.class.getName(), "").split(",");
	Hashtable<String,JMenu> submenus = new Hashtable<String,JMenu>();

	// add all extensions
	for (int i = 0; i < classnames.length; i++) {
	  String classname = classnames[i];
	  try {
	    MainMenuExtension ext = (MainMenuExtension) Class.forName(classname).newInstance();

	    // menuitem in a submenu?
	    JMenu submenu = null;
	    if (ext.getSubmenuTitle() != null) {
	      submenu = submenus.get(ext.getSubmenuTitle());
	      if (submenu == null) {
		submenu = new JMenu(ext.getSubmenuTitle());
		submenus.put(ext.getSubmenuTitle(), submenu);
		insertMenuItem(jMenuExtensions, submenu);
	      }
	    }

	    // create menu item
	    JMenuItem menuitem = new JMenuItem();
	    menuitem.setText(ext.getMenuTitle());
	    // does the extension need a frame or does it have its own ActionListener?
	    ActionListener listener = ext.getActionListener(m_Self);
	    if (listener != null) {
	      menuitem.addActionListener(listener);
	    }
	    else {
	      final JMenuItem finalMenuitem = menuitem;
	      final MainMenuExtension finalExt = ext;
	      menuitem.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent e) {
		  Component frame = createFrame(
		      m_Self, finalMenuitem.getText(), 
		      null, null, null, -1, -1, null, false, false);
		  finalExt.fillFrame(frame);
		  frame.setVisible(true);
		}
	      });
	    }

	    // sorted insert of menu item
	    if (submenu != null)
	      insertMenuItem(submenu, menuitem);
	    else
	      insertMenuItem(jMenuExtensions, menuitem);
	  }
	  catch (Exception e) {
	    e.printStackTrace();
	  }
	}
      }

      // Windows
      jMenuWindows = new JMenu("Windows");
      jMenuWindows.setMnemonic(java.awt.event.KeyEvent.VK_W);
      jMenuBar.add(jMenuWindows);
      jMenuWindows.setVisible(false);  // initially, there are no windows open

      // Help
      jMenuHelp = new JMenu();
      jMenuBar.add(jMenuHelp);
      jMenuHelp.setText("Help");
      jMenuHelp.setMnemonic('H');

      // Help/Homepage
      jMenuItemHelpHomepage = new JMenuItem();
      jMenuHelp.add(jMenuItemHelpHomepage);
      jMenuItemHelpHomepage.setText("Weka homepage");
      jMenuItemHelpHomepage.setMnemonic('H');
      jMenuItemHelpHomepage.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent evt) {
	  BrowserHelper.openURL(m_Self, "http://www.cs.waikato.ac.nz/~ml/weka/");
	}
      });

      jMenuHelp.add(new JSeparator());

/*      // Help/WekaDoc
      jMenuItemHelpWekaDoc = new JMenuItem();
      jMenuHelp.add(jMenuItemHelpWekaDoc);
      jMenuItemHelpWekaDoc.setText("Online documentation");
      jMenuItemHelpWekaDoc.setMnemonic('D');
      jMenuItemHelpWekaDoc.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent evt) {
	  BrowserHelper.openURL(m_Self, "http://weka.sourceforge.net/wekadoc/");
	}
      }); */

      // Help/WekaWiki
      jMenuItemHelpWekaWiki = new JMenuItem();
      jMenuHelp.add(jMenuItemHelpWekaWiki);
      jMenuItemHelpWekaWiki.setText("HOWTOs, code snippets, etc.");
      jMenuItemHelpWekaWiki.setMnemonic('W');
      jMenuItemHelpWekaWiki.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent evt) {
	  BrowserHelper.openURL(m_Self, "http://weka.wikispaces.com/");
	}
      });

      // Help/Sourceforge
      jMenuItemHelpSourceforge = new JMenuItem();
      jMenuHelp.add(jMenuItemHelpSourceforge);
      jMenuItemHelpSourceforge.setText("Weka on SourceForge");
      jMenuItemHelpSourceforge.setMnemonic('F');
      jMenuItemHelpSourceforge.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent evt) {
	  BrowserHelper.openURL(m_Self, "http://sourceforge.net/projects/weka/");
	}
      });

      jMenuHelp.add(new JSeparator());

      // Help/SystemInfo
      jMenuItemHelpSystemInfo = new JMenuItem();
      jMenuHelp.add(jMenuItemHelpSystemInfo);
      jMenuItemHelpSystemInfo.setText("SystemInfo");
      jMenuItemHelpHomepage.setMnemonic('S');
      jMenuItemHelpSystemInfo.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent evt) {
	  String title = jMenuItemHelpSystemInfo.getText();
	  if (!containsWindow(title)) {
	    // get info
	    Hashtable info = new SystemInfo().getSystemInfo();

	    // sort names
	    Vector names = new Vector();
	    Enumeration enm = info.keys();
	    while (enm.hasMoreElements())
	      names.add(enm.nextElement());
	    Collections.sort(names);

	    // generate table
	    String[][] data = new String[info.size()][2];
	    for (int i = 0; i < names.size(); i++) {
	      data[i][0] = names.get(i).toString();
	      data[i][1] = info.get(data[i][0]).toString();
	    }
	    String[] titles = new String[]{"Key", "Value"};
	    JTable table = new JTable(data, titles);

	    createFrame(
		m_Self, title, new JScrollPane(table), new BorderLayout(), 
		BorderLayout.CENTER, 800, 600, null, true, true);
	  }
	  else {
	    showWindow(getWindow(title));
	  }
	}
      });

      jMenuHelp.add(new JSeparator());

      // Help/About
      jMenuItemHelpAbout = new JMenuItem();
      jMenuHelp.add(jMenuItemHelpAbout);
      jMenuItemHelpAbout.setText("About");
      jMenuItemHelpAbout.setMnemonic('A');
      jMenuItemHelpAbout.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent evt) {
	  String title = jMenuItemHelpAbout.getText();
	  if (!containsWindow(title)) {
	    JPanel wekaPan = new JPanel();
	    wekaPan.setToolTipText("Weka, a native bird of New Zealand");
	    ImageIcon wii = new ImageIcon(Toolkit.getDefaultToolkit().getImage(ClassLoader.getSystemResource("weka/gui/weka3.gif")));
	    JLabel wekaLab = new JLabel(wii);
	    wekaPan.add(wekaLab);
	    Container frame = createFrame(
		m_Self, title, wekaPan, new BorderLayout(), 
		BorderLayout.CENTER, -1, -1, null, true, true);

	    JPanel titlePan = new JPanel();
	    titlePan.setLayout(new GridLayout(8,1));
	    titlePan.setBorder(BorderFactory.createEmptyBorder(10, 5, 10, 5));
	    titlePan.add(new JLabel("Waikato Environment for", SwingConstants.CENTER));
	    titlePan.add(new JLabel("Knowledge Analysis", SwingConstants.CENTER));
	    titlePan.add(new JLabel(""));
	    titlePan.add(new JLabel("Version " + Version.VERSION, SwingConstants.CENTER));
	    titlePan.add(new JLabel(""));
	    titlePan.add(new JLabel("(c) " + Copyright.getFromYear() + " - " + Copyright.getToYear(), SwingConstants.CENTER));
	    titlePan.add(new JLabel(Copyright.getOwner(), SwingConstants.CENTER));
	    titlePan.add(new JLabel(Copyright.getAddress(), SwingConstants.CENTER));

	    if (frame instanceof ChildFrameMDI) {
	      ((ChildFrameMDI) frame).getContentPane().add(titlePan, BorderLayout.NORTH);
	      ((ChildFrameMDI) frame).pack();
	    }
	    else if (frame instanceof ChildFrameSDI) {
	      ((ChildFrameSDI) frame).getContentPane().add(titlePan, BorderLayout.NORTH);
	      ((ChildFrameSDI) frame).pack();
	    }
	  }
	  else {
	    showWindow(getWindow(title));
	  }
	}
      });

      // size + position
      int screenHeight = getGraphicsConfiguration().getBounds().height;
      int screenWidth  = getGraphicsConfiguration().getBounds().width;
      if (m_GUIType == GUI_MDI) {
	int newHeight = (int) (((double) screenHeight) * 0.75);
	int newWidth  = (int) (((double) screenWidth)  * 0.75);
	setSize(
	    1000 > newWidth  ? newWidth  : 1000,
		800  > newHeight ? newHeight : 800);
	setLocation(
	    (screenWidth - getBounds().width) / 2,
	    (screenHeight - getBounds().height) / 2);
      }
      else if (m_GUIType == GUI_SDI) {
	pack();
	setSize(screenWidth, getHeight());
	setLocation(0, 0);
      }
    } 
    catch (Exception e) {
      e.printStackTrace();
    }
  }
  
  /**
   * creates and displays the title.
   * 
   * @param title 	the additional part of the title
   */
  protected void createTitle(String title) {
    String	newTitle;
    
    newTitle = "Weka " + new Version();
    if (title.length() != 0)
      newTitle += " - " + title;
    
    setTitle(newTitle);
  }
  
  /**
   * adds the given child frame to the list of frames.
   * 
   * @param c 		the child frame to add
   */
  public void addChildFrame(Container c) {
    m_ChildFrames.add(c);
    windowListChanged();
  }
  
  /**
   * tries to remove the child frame, it returns true if it could do such.
   * 
   * @param c 		the child frame to remove
   * @return 		true if the child frame could be removed
   */
  public boolean removeChildFrame(Container c) {
    boolean result = m_ChildFrames.remove(c);
    windowListChanged();
    return result;
  }
  
  /**
   * brings child frame to the top.
   * 
   * @param c 		the frame to activate
   * @return 		true if frame was activated
   */
  public boolean showWindow(Container c) {
    boolean        	result;
    ChildFrameMDI	mdiFrame;
    ChildFrameSDI	sdiFrame;
    
    if (c != null) {
      try {
	if (c instanceof ChildFrameMDI) {
	  mdiFrame = (ChildFrameMDI) c;
	  mdiFrame.setIcon(false);
	  mdiFrame.toFront();
	  createTitle(mdiFrame.getTitle());
	}
	else if (c instanceof ChildFrameSDI) {
	  sdiFrame = (ChildFrameSDI) c;
	  sdiFrame.setExtendedState(JFrame.NORMAL);
	  sdiFrame.toFront();
	  createTitle(sdiFrame.getTitle());
	}
      }
      catch (Exception e) {
	e.printStackTrace();
      }
      result = true;
    }
    else {
      result = false;
    }
    
    return result;
  }
  
  /**
   * brings the first frame to the top that is of the specified
   * window class.
   *  
   * @param windowClass	the class to display the first child for
   * @return		true, if a child was found and brought to front
   */
  public boolean showWindow(Class windowClass) {
    return showWindow(getWindow(windowClass));
  }
  
  /**
   * returns all currently open frames.
   * 
   * @return 		an iterator over all currently open frame
   */
  public Iterator getWindowList() {
    return m_ChildFrames.iterator();
  }

  /**
   * returns the first instance of the given window class, null if none can be 
   * found.
   * 
   * @param windowClass	the class to retrieve the first instance for
   * @return		null, if no instance can be found
   */
  public Container getWindow(Class windowClass) {
    Container	result;
    Iterator	iter;
    Container	current;
    
    result = null;
    iter   = getWindowList();
    while (iter.hasNext()) {
      current = (Container) iter.next();
      if (current.getClass() == windowClass) {
        result = current;
        break;
      }
    }
    
    return result;
  }

  /**
   * returns the first window with the given title, null if none can be 
   * found.
   * 
   * @param title	the title to look for
   * @return		null, if no instance can be found
   */
  public Container getWindow(String title) {
    Container	result;
    Iterator	iter;
    Container	current;
    boolean	found;
    
    result = null;
    iter   = getWindowList();
    while (iter.hasNext()) {
      current = (Container) iter.next();
      found   = false;
      
      if (current instanceof ChildFrameMDI)
	found = ((ChildFrameMDI) current).getTitle().equals(title);
      else if (current instanceof ChildFrameSDI)
	found = ((ChildFrameSDI) current).getTitle().equals(title);
	
      if (found) {
        result = current;
        break;
      }
    }
    
    return result;
  }
  
  /**
   * checks, whether an instance of the given window class is already in
   * the Window list.
   * 
   * @param windowClass	the class to check for an instance in the current
   * 			window list
   * @return		true if the class is already listed in the Window list
   */
  public boolean containsWindow(Class windowClass) {
    return (getWindow(windowClass) != null);
  }
  
  /**
   * checks, whether a window with the given title is already in
   * the Window list.
   * 
   * @param title	the title to check for in the current window list
   * @return		true if a window with the given title is already 
   * 			listed in the Window list
   */
  public boolean containsWindow(String title) {
    return (getWindow(title) != null);
  }
  
  /**
   * minimizes all windows.
   */
  public void minimizeWindows() {
    Iterator	iter;
    Container	frame;
    
    iter = getWindowList();
    while (iter.hasNext()) {
      frame = (Container) iter.next();
      try {
	if (frame instanceof ChildFrameMDI)
	  ((ChildFrameMDI) frame).setIcon(true);
	else if (frame instanceof ChildFrameSDI)
	  ((ChildFrameSDI) frame).setExtendedState(JFrame.ICONIFIED);
      }
      catch (Exception e) {
	e.printStackTrace();
      }
    }
  }
  
  /**
   * restores all windows.
   */
  public void restoreWindows() {
    Iterator	iter;
    Container	frame;
    
    iter = getWindowList();
    while (iter.hasNext()) {
      frame = (Container) iter.next();
      try {
	if (frame instanceof ChildFrameMDI)
	  ((ChildFrameMDI) frame).setIcon(false);
	else if (frame instanceof ChildFrameSDI)
	  ((ChildFrameSDI) frame).setExtendedState(JFrame.NORMAL);
    }
      catch (Exception e) {
	e.printStackTrace();
      }
    }
  }
  
  /**
   * is called when window list changed somehow (add or remove).
   */
  public void windowListChanged() {
    createWindowMenu();
  }
  
  /**
   * creates the menu of currently open windows.
   */
  protected synchronized void createWindowMenu() {
    Iterator          iter;
    JMenuItem         menuItem;
    int	              startIndex;
    
    // remove all existing entries
    jMenuWindows.removeAll();
    
    // minimize + restore + separator
    menuItem = new JMenuItem("Minimize");
    menuItem.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent evt) {
        minimizeWindows();
      }
    });
    jMenuWindows.add(menuItem);
    
    menuItem = new JMenuItem("Restore");
    menuItem.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent evt) {
        restoreWindows();
      }
    });
    jMenuWindows.add(menuItem);
    
    jMenuWindows.addSeparator();
    
    // windows
    startIndex = jMenuWindows.getMenuComponentCount() - 1;
    iter = getWindowList();
    jMenuWindows.setVisible(iter.hasNext());
    while (iter.hasNext()) {
      Container frame = (Container) iter.next();
      if (frame instanceof ChildFrameMDI)
	menuItem = new JMenuItem(((ChildFrameMDI) frame).getTitle());
      else if (frame instanceof ChildFrameSDI)
	menuItem = new JMenuItem(((ChildFrameSDI) frame).getTitle());
      insertMenuItem(jMenuWindows, menuItem, startIndex);
      menuItem.setActionCommand(Integer.toString(frame.hashCode()));
      menuItem.addActionListener(new ActionListener() {
        public void actionPerformed(ActionEvent evt) {
          Container frame = null;
          Iterator iter = getWindowList();
          while (iter.hasNext()) {
            frame = (Container) iter.next();
            String hashFrame = Integer.toString(frame.hashCode());
            if (hashFrame.equals(evt.getActionCommand())) {
              showWindow(frame);
              break;
            }
          }
          showWindow(frame);
        }
      });
    }
  }
  
  /**
   * Shows or hides this component depending on the value of parameter b.
   * 
   * @param b		if true, shows this component; otherwise, hides this 
   * 			component
   */
  public void setVisible(boolean b) {
    super.setVisible(b);
    
    if (b)
      paint(this.getGraphics());
  }
  
  /**
   * Create the singleton instance of the Main GUI.
   * 
   * @param args 	commandline options
   */
  public static void createSingleton(String[] args) {
    if (m_MainSingleton == null)
      m_MainSingleton = new Main();
    
    // set options
    try {
      m_MainSingleton.setOptions(args);
    }
    catch (Exception e) {
      e.printStackTrace();
    }

    // notify listeners (if any)
    for (int i = 0; i < m_StartupListeners.size(); i++)
      ((StartUpListener) m_StartupListeners.elementAt(i)).startUpComplete();
  }

  /**
   * Return the singleton instance of the Main GUI.
   *
   * @return the singleton instance
   */
  public static Main getSingleton() {
    return m_MainSingleton;
  }

  /**
   * Add a listener to be notified when startup is complete.
   * 
   * @param s 		a listener to add
   */
  public static void addStartupListener(StartUpListener s) {
    m_StartupListeners.add(s);
  }

  /**
   * Gets an enumeration describing the available options.
   *
   * @return 		an enumeration of all the available options.
   */
  public Enumeration listOptions(){
    Vector        	result;
    String		desc;
    SelectedTag		tag;
    int			i;

    result = new Vector();

    desc  = "";
    for (i = 0; i < TAGS_GUI.length; i++) {
      tag = new SelectedTag(TAGS_GUI[i].getID(), TAGS_GUI);
      desc  +=   "\t" + tag.getSelectedTag().getIDStr() 
      	       + " = " + tag.getSelectedTag().getReadable()
      	       + "\n";
    }
    result.addElement(new Option(
	"\tDetermines the layout of the GUI:\n"
	+ desc
	+ "\t(default: " + new SelectedTag(GUI_MDI, TAGS_GUI) + ")",
	"gui", 1, "-gui " + Tag.toOptionList(TAGS_GUI)));

    return result.elements();
  }
  
  /**
   * returns the options of the current setup.
   *
   * @return		the current options
   */
  public String[] getOptions(){
    Vector<String>    	result;

    result = new Vector();

    result.add("-gui");
    result.add("" + getGUIType());

    return result.toArray(new String[result.size()]);	  
  }

  /**
   * Parses the options for this object. <p/>
   *
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -gui &lt;MDI|SDI&gt;
   *  Determines the layout of the GUI:
   *  MDI = MDI Layout
   *  SDI = SDI Layout
   *  (default: MDI)</pre>
   * 
   <!-- options-end -->
   *
   * @param options	the options to use
   * @throws Exception	if setting of options fails
   */
  public void setOptions(String[] options) throws Exception {
    String	tmpStr;

    tmpStr = Utils.getOption("gui", options);
    if (tmpStr.length() != 0)
      setGUIType(new SelectedTag(tmpStr, TAGS_GUI));
    else
      setGUIType(new SelectedTag(GUI_MDI, TAGS_GUI));
  }

  /**
   * Sets the type of GUI to use. 
   *
   * @param value 	.the GUI type
   */
  public void setGUIType(SelectedTag value) {
    if (value.getTags() == TAGS_GUI) {
      m_GUIType = value.getSelectedTag().getID();
      initGUI();
    }
  }

  /**
   * Gets the currently set type of GUI to display. 
   * 
   * @return 		the current GUI Type.
   */
  public SelectedTag getGUIType() {
    return new SelectedTag(m_GUIType, TAGS_GUI);
  }
  
  /**
   * starts the application.
   * 
   * @param args	the commandline arguments - ignored
   */
  public static void main(String[] args) {
    weka.core.logging.Logger.log(weka.core.logging.Logger.Level.INFO, "Logging started");
    
    LookAndFeel.setLookAndFeel();
    
    try {
      // uncomment the following line to disable the memory management:
      //m_Memory.setEnabled(false);

      // help?
      if (Utils.getFlag('h', args)) {
	System.out.println();
	System.out.println("Help requested.");
	System.out.println();
	System.out.println("General options:");
	System.out.println();
	System.out.println("-h");
	System.out.println("\tprints this help screen");
	System.out.println();

	Enumeration enu = new Main().listOptions();
	while (enu.hasMoreElements()) {
	  Option option = (Option) enu.nextElement();
	  System.out.println(option.synopsis());
	  System.out.println(option.description());
	}

	System.out.println();
	System.exit(0);
      }
      
      // setup splash screen
      Main.addStartupListener(new weka.gui.beans.StartUpListener() {
        public void startUpComplete() {
          m_MainCommandline = Main.getSingleton();
          m_MainCommandline.setVisible(true);
        }
      });
      Main.addStartupListener(new StartUpListener() {
        public void startUpComplete() {
          SplashWindow.disposeSplash();
        }
      });
      SplashWindow.splash(ClassLoader.getSystemResource("weka/gui/images/weka_splash.gif"));

      // start GUI
      final String[] options = (String[]) args.clone();
      Thread nt = new Thread() {
	public void run() {
	  weka.gui.SplashWindow.invokeMethod(
	      Main.class.getName(), "createSingleton", options);
	}
      };
      nt.start();
      
      Thread memMonitor = new Thread() {
	public void run() {
	  while(true) {
	    try {
	      Thread.sleep(4000);
	      System.gc();
	      
	      if (m_Memory.isOutOfMemory()) {
		// clean up
		m_MainCommandline = null;
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
    }
    catch (Exception ex) {
      ex.printStackTrace();
      System.err.println(ex.getMessage());
    }
  }
}
