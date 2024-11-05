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
 *    TreeVisualizer.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.treevisualizer;

import weka.core.Instances;
import weka.core.Utils;
import weka.gui.visualize.PrintablePanel;
import weka.gui.visualize.VisualizePanel;
import weka.gui.visualize.VisualizeUtils;

import java.awt.Color;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.io.FileReader;
import java.io.IOException;
import java.io.StringReader;
import java.util.Properties;

import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPopupMenu;
import javax.swing.JRadioButton;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.JTextField;
import javax.swing.Timer;

/**
 * Class for displaying a Node structure in Swing. <p>
 *
 * To work this class simply create an instance of it.<p>
 *
 * Assign it to a window or other such object.<p>
 *
 * Resize it to the desired size.<p>
 *
 *
 * When using the Displayer hold the left mouse button to drag the 
 * tree around. <p>
 *
 * Click the left mouse button with ctrl to shrink the size of the tree 
 * by half. <p>
 *
 * Click and drag with the left mouse button and shift to draw a box,
 * when the left mouse button is released the contents of the box 
 * will be magnified 
 * to fill the screen. <p> <p>
 *
 * Click the right mouse button to bring up a menu. <p>
 * Most options are self explanatory.<p>
 *
 * Select Auto Scale to set the tree to it's optimal display size.
 *
 * @author Malcolm Ware (mfw4@cs.waikato.ac.nz)
 * @version $Revision: 4960 $
 */
public class TreeVisualizer
  extends PrintablePanel
  implements MouseMotionListener, MouseListener, ActionListener, ItemListener {

  /** for serialization */
  private static final long serialVersionUID = -8668637962504080749L;

  /** the props file. */
  public final static String PROPERTIES_FILE = "weka/gui/treevisualizer/TreeVisualizer.props";
  
  /** The placement algorithm for the Node structure. */
  private NodePlace m_placer;  

  /** The top Node. */
  private Node m_topNode;
  
  /** The postion of the view relative to the tree. */
  private Dimension m_viewPos;        

  /** The size of the tree in pixels. */
  private Dimension m_viewSize;      
  
  /** The font used to display the tree. */
  private Font m_currentFont;        

  /** The size information for the current font. */
  private FontMetrics m_fontSize;    

  /** The number of Nodes in the tree. */
  private int m_numNodes;

  /** The number of levels in the tree. */
  private int m_numLevels;     

  /** An array with the Nodes sorted into it and display information 
   * about the Nodes. */
  private NodeInfo[] m_nodes;

  /** An array with the Edges sorted into it and display information 
   * about the Edges. */
  private EdgeInfo[] m_edges;
                     
  /** A timer to keep the frame rate constant. */
  private Timer m_frameLimiter;      
                         
  /** Describes the action the user is performing. */
  private int m_mouseState;           

  /** A variable used to tag the start pos of a user action. */
  private Dimension m_oldMousePos;

  /** A variable used to tag the most current point of a user action. */
  private Dimension m_newMousePos;

  /** A variable used to determine for the clicked method if any other 
   * mouse state has already taken place. */
  private boolean m_clickAvailable;

  /** A variable used to remember the desired view pos. */
  private Dimension m_nViewPos;     

  /** A variable used to remember the desired tree size. */
  private Dimension m_nViewSize;       

  /** The number of frames left to calculate. */
  private int m_scaling;         

  /** A right (or middle) click popup menu. */
  private JPopupMenu m_winMenu;

  /** An option on the win_menu */
  private JMenuItem m_topN;

  /** An option on the win_menu*/
  private JMenuItem m_fitToScreen;

  /** An option on the win_menu */
  private JMenuItem m_autoScale;

  /** A sub group on the win_menu */
  private JMenu m_selectFont;

  /** A grouping for the font choices */
  private ButtonGroup m_selectFontGroup;

  /** A font choice. */
  private JRadioButtonMenuItem m_size24;

  /** A font choice. */
  private JRadioButtonMenuItem m_size22;

  /** A font choice. */
  private JRadioButtonMenuItem m_size20;

  /** A font choice. */
  private JRadioButtonMenuItem m_size18;

  /** A font choice. */
  private JRadioButtonMenuItem m_size16;

  /** A font choice. */
  private JRadioButtonMenuItem m_size14;

  /** A font choice. */
  private JRadioButtonMenuItem m_size12;

  /** A font choice. */
  private JRadioButtonMenuItem m_size10;

  /** A font choice. */
  private JRadioButtonMenuItem m_size8;

  /** A font choice. */
  private JRadioButtonMenuItem m_size6;

  /** A font choice. */
  private JRadioButtonMenuItem m_size4;

  /** A font choice. */
  private JRadioButtonMenuItem m_size2;

  /** A font choice. */
  private JRadioButtonMenuItem m_size1;

  /** An option on the win menu. */
  private JMenuItem m_accept;

  /** A right or middle click popup menu for nodes. */
  private JPopupMenu m_nodeMenu;

  /** A visualize choice for the node, may not be available. */
  private JMenuItem m_visualise;

  /** 
   * An add children to Node choice, This is only available if the tree
   * display has a treedisplay listerner added to it.
   */
  private JMenuItem m_addChildren;

  /** Similar to add children but now it removes children. */
  private JMenuItem m_remChildren;

  /** Use this to have J48 classify this node. */
  private JMenuItem m_classifyChild;
  
  /** Use this to dump the instances from this node to the vis panel. */
  private JMenuItem m_sendInstances;

  /** The subscript for the currently selected node (this is an internal 
   * thing, so the user is unaware of this). */
  private int m_focusNode;

  /**
   * The Node the user is currently focused on , this is similar to 
   * focus node except that it is used by other 
   * classes rather than this one.
   */
  private int m_highlightNode;
  
  /* A pointer to this tree's classifier if a classifier is using it. */
  //private UserClassifier classer;
  private TreeDisplayListener m_listener;

  private JTextField m_searchString;
  private JDialog m_searchWin;
  private JRadioButton m_caseSen;

  /** the font color. */
  protected Color m_FontColor = null;

  /** the background color. */
  protected Color m_BackgroundColor = null;

  /** the node color. */
  protected Color m_NodeColor = null;

  /** the line color. */
  protected Color m_LineColor = null;

  /** the color of the zoombox. */
  protected Color m_ZoomBoxColor = null;

  /** the XOR color of the zoombox. */
  protected Color m_ZoomBoxXORColor = null;
  
  /** whether to show the border or not. */
  protected boolean m_ShowBorder = true;
  
  ///////////////////

  //this is the event fireing stuff


  /**
   * Constructs Displayer to display a tree provided in a dot format.
   * Uses the NodePlacer to place the Nodes.
   * @param tdl listener 
   * @param dot string containing the dot representation of the tree to
   * display
   * @param p the algorithm to be used to position the nodes.
   */
  public TreeVisualizer(TreeDisplayListener tdl, String dot, NodePlace p) {
    super();

    initialize();
    
    //generate the node structure in here
    if (m_ShowBorder)
      setBorder(BorderFactory.createTitledBorder("Tree View")); 
    m_listener = tdl;

    TreeBuild builder = new TreeBuild();
    
    Node n = null;
    NodePlace arrange = new PlaceNode2();
    n = builder.create(new StringReader(dot));
    //    System.out.println(n.getCount(n, 0));
    //if the size needs to be automatically alocated I will do it here
    m_highlightNode = 5;
    m_topNode = n;
    m_placer = p;
    m_placer.place(m_topNode);
    m_viewPos = new Dimension(0, 0);    //will be adjusted 
    m_viewSize = new Dimension(800, 600);   //I allocate this now so that
    //the tree will be visible
    //when the panel is enlarged

    m_nViewPos = new Dimension(0, 0);            
    m_nViewSize = new Dimension(800, 600);          
                                      
    m_scaling = 0;
    
    m_numNodes = m_topNode.getCount(m_topNode,0);   //note the second 
    //argument must be a zero, this is a 
    //recursive function

    m_numLevels = m_topNode.getHeight(m_topNode,0);
  
    m_nodes = new NodeInfo[m_numNodes];
    m_edges = new EdgeInfo[m_numNodes-1];
    

    arrayFill(m_topNode, m_nodes, m_edges);
    
    changeFontSize(12);

    m_mouseState = 0;
    m_oldMousePos = new Dimension(0, 0);
    m_newMousePos = new Dimension(0, 0);
    m_frameLimiter = new Timer(120, this);



    m_winMenu = new JPopupMenu();
    m_topN = new JMenuItem("Center on Top Node");           //note to change 
    //language change this line
    m_topN.setActionCommand("Center on Top Node");          //but not this one,
    //same for all menu items
    m_fitToScreen = new JMenuItem("Fit to Screen");
    m_fitToScreen.setActionCommand("Fit to Screen");
    //unhide = new JMenuItem("Unhide all Nodes");
    m_selectFont = new JMenu("Select Font");
    m_selectFont.setActionCommand("Select Font");
    m_autoScale = new JMenuItem("Auto Scale");
    m_autoScale.setActionCommand("Auto Scale");
    m_selectFontGroup = new ButtonGroup();
    
    m_accept = new JMenuItem("Accept The Tree");
    m_accept.setActionCommand("Accept The Tree");
    
    m_winMenu.add(m_topN);
    m_winMenu.addSeparator();
    m_winMenu.add(m_fitToScreen);
    m_winMenu.add(m_autoScale);
    //m_winMenu.addSeparator();
    //m_winMenu.add(unhide);
    m_winMenu.addSeparator();
    m_winMenu.add(m_selectFont);

    if (m_listener != null) {
      m_winMenu.addSeparator();
      m_winMenu.add(m_accept);
    }
    
    m_topN.addActionListener(this);
    m_fitToScreen.addActionListener(this);
    //unhide.addActionListener(this);
    m_autoScale.addActionListener(this);
    m_accept.addActionListener(this);
        
    m_size24 = new JRadioButtonMenuItem("Size 24",false);//,select_font_group);
    m_size22 = new JRadioButtonMenuItem("Size 22",false);//,select_font_group);
    m_size20 = new JRadioButtonMenuItem("Size 20",false);//,select_font_group);
    m_size18 = new JRadioButtonMenuItem("Size 18",false);//,select_font_group);
    m_size16 = new JRadioButtonMenuItem("Size 16",false);//,select_font_group);
    m_size14 = new JRadioButtonMenuItem("Size 14",false);//,select_font_group);
    m_size12 = new JRadioButtonMenuItem("Size 12",true);//,select_font_group);
    m_size10 = new JRadioButtonMenuItem("Size 10",false);//,select_font_group);
    m_size8 = new JRadioButtonMenuItem("Size 8",false);//,select_font_group);
    m_size6 = new JRadioButtonMenuItem("Size 6",false);//,select_font_group);
    m_size4 = new JRadioButtonMenuItem("Size 4",false);//,select_font_group);
    m_size2 = new JRadioButtonMenuItem("Size 2",false);//,select_font_group);
    m_size1 = new JRadioButtonMenuItem("Size 1",false);//,select_font_group);

    m_size24.setActionCommand("Size 24");//,select_font_group);
    m_size22.setActionCommand("Size 22");//,select_font_group);
    m_size20.setActionCommand("Size 20");//,select_font_group);
    m_size18.setActionCommand("Size 18");//,select_font_group);
    m_size16.setActionCommand("Size 16");//,select_font_group);
    m_size14.setActionCommand("Size 14");//,select_font_group);
    m_size12.setActionCommand("Size 12");//,select_font_group);
    m_size10.setActionCommand("Size 10");//,select_font_group);
    m_size8.setActionCommand("Size 8");//,select_font_group);
    m_size6.setActionCommand("Size 6");//,select_font_group);
    m_size4.setActionCommand("Size 4");//,select_font_group);
    m_size2.setActionCommand("Size 2");//,select_font_group);
    m_size1.setActionCommand("Size 1");//,select_font_group);
    
    
    m_selectFontGroup.add(m_size24);
    m_selectFontGroup.add(m_size22);
    m_selectFontGroup.add(m_size20);
    m_selectFontGroup.add(m_size18);
    m_selectFontGroup.add(m_size16);
    m_selectFontGroup.add(m_size14);
    m_selectFontGroup.add(m_size12);
    m_selectFontGroup.add(m_size10);
    m_selectFontGroup.add(m_size8);
    m_selectFontGroup.add(m_size6);
    m_selectFontGroup.add(m_size4);
    m_selectFontGroup.add(m_size2);
    m_selectFontGroup.add(m_size1);

    
    m_selectFont.add(m_size24);
    m_selectFont.add(m_size22);
    m_selectFont.add(m_size20);
    m_selectFont.add(m_size18);
    m_selectFont.add(m_size16);
    m_selectFont.add(m_size14);
    m_selectFont.add(m_size12);
    m_selectFont.add(m_size10);
    m_selectFont.add(m_size8);
    m_selectFont.add(m_size6);
    m_selectFont.add(m_size4);
    m_selectFont.add(m_size2);
    m_selectFont.add(m_size1);


    m_size24.addItemListener(this);
    m_size22.addItemListener(this);
    m_size20.addItemListener(this);
    m_size18.addItemListener(this);
    m_size16.addItemListener(this);
    m_size14.addItemListener(this);
    m_size12.addItemListener(this);
    m_size10.addItemListener(this);
    m_size8.addItemListener(this);
    m_size6.addItemListener(this);
    m_size4.addItemListener(this);
    m_size2.addItemListener(this);
    m_size1.addItemListener(this);

    /*
      search_string = new JTextField(22);
      search_win = new JDialog();
      case_sen = new JRadioButton("Case Sensitive");



      search_win.getContentPane().setLayout(null);
      search_win.setSize(300, 200);
 
      search_win.getContentPane().add(search_string);
      search_win.getContentPane().add(case_sen);

      search_string.setLocation(50, 70);
      case_sen.setLocation(50, 120);
      case_sen.setSize(100, 24); 
      search_string.setSize(100, 24);
      //search_string.setVisible(true);
      //case_sen.setVisible(true);

      //search_win.setVisible(true);
    */

    m_nodeMenu = new JPopupMenu();
    /* A visualize choice for the node, may not be available. */
    m_visualise = new JMenuItem("Visualize The Node");
    m_visualise.setActionCommand("Visualize The Node");
    m_visualise.addActionListener(this);
    m_nodeMenu.add(m_visualise);
   
    if (m_listener != null) {
      m_remChildren = new JMenuItem("Remove Child Nodes");
      m_remChildren.setActionCommand("Remove Child Nodes");
      m_remChildren.addActionListener(this);
      m_nodeMenu.add(m_remChildren);
      
      
      m_classifyChild = new JMenuItem("Use Classifier...");
      m_classifyChild.setActionCommand("classify_child");
      m_classifyChild.addActionListener(this);
      m_nodeMenu.add(m_classifyChild);
      
      /*m_sendInstances = new JMenuItem("Add Instances To Viewer");
      m_sendInstances.setActionCommand("send_instances");
      m_sendInstances.addActionListener(this);
      m_nodeMenu.add(m_sendInstances); */
      
    }
    
    m_focusNode = -1;
    m_highlightNode = -1;
    
    addMouseMotionListener(this);
    addMouseListener(this);
    //repaint();
    //frame_limiter.setInitialDelay();
    m_frameLimiter.setRepeats(false);
    m_frameLimiter.start();
  }
  
  /**
   * Constructs Displayer with the specified Node as the top 
   * of the tree, and uses the NodePlacer to place the Nodes.
   * @param tdl listener.
   * @param n the top Node of the tree to be displayed.
   * @param p the algorithm to be used to position the nodes.
   */  
  public TreeVisualizer(TreeDisplayListener tdl, Node n, NodePlace p) {
    super();

    initialize();
    
    //if the size needs to be automatically alocated I will do it here
    if (m_ShowBorder)
      setBorder(BorderFactory.createTitledBorder("Tree View")); 
    m_listener = tdl;
    m_topNode = n;
    m_placer = p;
    m_placer.place(m_topNode);
    m_viewPos = new Dimension(0, 0);    //will be adjusted 
    m_viewSize = new Dimension(800, 600);   //I allocate this now so that
    //the tree will be visible
    //when the panel is enlarged

    m_nViewPos = new Dimension(0, 0);            
    m_nViewSize = new Dimension(800, 600);          
                                      
    m_scaling = 0;
    
    m_numNodes = m_topNode.getCount(m_topNode,0);   //note the second argument 
    //must be a zero, this is a 
    //recursive function

    m_numLevels = m_topNode.getHeight(m_topNode,0);
  
    m_nodes = new NodeInfo[m_numNodes];
    m_edges = new EdgeInfo[m_numNodes-1];

    arrayFill(m_topNode, m_nodes, m_edges);
    
    changeFontSize(12);

    m_mouseState = 0;
    m_oldMousePos = new Dimension(0, 0);
    m_newMousePos = new Dimension(0, 0);
    m_frameLimiter = new Timer(120, this);





    m_winMenu = new JPopupMenu();
    m_topN = new JMenuItem("Center on Top Node");           //note to change 
    //language change this line
    m_topN.setActionCommand("Center on Top Node");          //but not this 
    //one, same for all menu items
    m_fitToScreen = new JMenuItem("Fit to Screen");
    m_fitToScreen.setActionCommand("Fit to Screen");
    //unhide = new JMenuItem("Unhide all Nodes");
    m_selectFont = new JMenu("Select Font");
    m_selectFont.setActionCommand("Select Font");
    m_autoScale = new JMenuItem("Auto Scale");
    m_autoScale.setActionCommand("Auto Scale");
    m_selectFontGroup = new ButtonGroup();
    
    m_accept = new JMenuItem("Accept The Tree");
    m_accept.setActionCommand("Accept The Tree");
    
    m_winMenu.add(m_topN);
    m_winMenu.addSeparator();
    m_winMenu.add(m_fitToScreen);
    m_winMenu.add(m_autoScale);
    m_winMenu.addSeparator();
    //m_winMenu.add(unhide);
    m_winMenu.addSeparator();
    m_winMenu.add(m_selectFont);
    m_winMenu.addSeparator();

    if (m_listener != null) {
      m_winMenu.add(m_accept);
    }
    
    m_topN.addActionListener(this);
    m_fitToScreen.addActionListener(this);
    //unhide.addActionListener(this);
    m_autoScale.addActionListener(this);
    m_accept.addActionListener(this);
        
    m_size24 = new JRadioButtonMenuItem("Size 24",false);//,select_font_group);
    m_size22 = new JRadioButtonMenuItem("Size 22",false);//,select_font_group);
    m_size20 = new JRadioButtonMenuItem("Size 20",false);//,select_font_group);
    m_size18 = new JRadioButtonMenuItem("Size 18",false);//,select_font_group);
    m_size16 = new JRadioButtonMenuItem("Size 16",false);//,select_font_group);
    m_size14 = new JRadioButtonMenuItem("Size 14",false);//,select_font_group);
    m_size12 = new JRadioButtonMenuItem("Size 12",true);//,select_font_group);
    m_size10 = new JRadioButtonMenuItem("Size 10",false);//,select_font_group);
    m_size8 = new JRadioButtonMenuItem("Size 8",false);//,select_font_group);
    m_size6 = new JRadioButtonMenuItem("Size 6",false);//,select_font_group);
    m_size4 = new JRadioButtonMenuItem("Size 4",false);//,select_font_group);
    m_size2 = new JRadioButtonMenuItem("Size 2",false);//,select_font_group);
    m_size1 = new JRadioButtonMenuItem("Size 1",false);//,select_font_group);

    m_size24.setActionCommand("Size 24");//,select_font_group);
    m_size22.setActionCommand("Size 22");//,select_font_group);
    m_size20.setActionCommand("Size 20");//,select_font_group);
    m_size18.setActionCommand("Size 18");//,select_font_group);
    m_size16.setActionCommand("Size 16");//,select_font_group);
    m_size14.setActionCommand("Size 14");//,select_font_group);
    m_size12.setActionCommand("Size 12");//,select_font_group);
    m_size10.setActionCommand("Size 10");//,select_font_group);
    m_size8.setActionCommand("Size 8");//,select_font_group);
    m_size6.setActionCommand("Size 6");//,select_font_group);
    m_size4.setActionCommand("Size 4");//,select_font_group);
    m_size2.setActionCommand("Size 2");//,select_font_group);
    m_size1.setActionCommand("Size 1");//,select_font_group);



    
    
    m_selectFontGroup.add(m_size24);
    m_selectFontGroup.add(m_size22);
    m_selectFontGroup.add(m_size20);
    m_selectFontGroup.add(m_size18);
    m_selectFontGroup.add(m_size16);
    m_selectFontGroup.add(m_size14);
    m_selectFontGroup.add(m_size12);
    m_selectFontGroup.add(m_size10);
    m_selectFontGroup.add(m_size8);
    m_selectFontGroup.add(m_size6);
    m_selectFontGroup.add(m_size4);
    m_selectFontGroup.add(m_size2);
    m_selectFontGroup.add(m_size1);



    
    m_selectFont.add(m_size24);
    m_selectFont.add(m_size22);
    m_selectFont.add(m_size20);
    m_selectFont.add(m_size18);
    m_selectFont.add(m_size16);
    m_selectFont.add(m_size14);
    m_selectFont.add(m_size12);
    m_selectFont.add(m_size10);
    m_selectFont.add(m_size8);
    m_selectFont.add(m_size6);
    m_selectFont.add(m_size4);
    m_selectFont.add(m_size2);
    m_selectFont.add(m_size1);


    m_size24.addItemListener(this);
    m_size22.addItemListener(this);
    m_size20.addItemListener(this);
    m_size18.addItemListener(this);
    m_size16.addItemListener(this);
    m_size14.addItemListener(this);
    m_size12.addItemListener(this);
    m_size10.addItemListener(this);
    m_size8.addItemListener(this);
    m_size6.addItemListener(this);
    m_size4.addItemListener(this);
    m_size2.addItemListener(this);
    m_size1.addItemListener(this);




    /*
      search_string = new JTextField(22);
      search_win = new JDialog();
      case_sen = new JRadioButton("Case Sensitive");



      search_win.getContentPane().setLayout(null);
      search_win.setSize(300, 200);
 
      search_win.getContentPane().add(search_string);
      search_win.getContentPane().add(case_sen);

      search_string.setLocation(50, 70);
      case_sen.setLocation(50, 120);
      case_sen.setSize(100, 24); 
      search_string.setSize(100, 24);
      //search_string.setVisible(true);
      //case_sen.setVisible(true);

      search_win.setVisible(true);
    */


    m_nodeMenu = new JPopupMenu();
    /* A visualize choice for the node, may not be available. */
    m_visualise = new JMenuItem("Visualize The Node");
    m_visualise.setActionCommand("Visualize The Node");
    m_visualise.addActionListener(this);
    m_nodeMenu.add(m_visualise);

    if (m_listener != null) {
      m_remChildren = new JMenuItem("Remove Child Nodes");
      m_remChildren.setActionCommand("Remove Child Nodes");
      m_remChildren.addActionListener(this);
      m_nodeMenu.add(m_remChildren);
      
      m_classifyChild = new JMenuItem("Use Classifier...");
      m_classifyChild.setActionCommand("classify_child");
      m_classifyChild.addActionListener(this);
      m_nodeMenu.add(m_classifyChild);
      
      m_sendInstances = new JMenuItem("Add Instances To Viewer");
      m_sendInstances.setActionCommand("send_instances");
      m_sendInstances.addActionListener(this);
      m_nodeMenu.add(m_sendInstances);
      
      
    }
  
    m_focusNode = -1;
    m_highlightNode = -1;
    

    addMouseMotionListener(this);
    addMouseListener(this);

  
    //repaint();

    //frame_limiter.setInitialDelay();
    m_frameLimiter.setRepeats(false);
    m_frameLimiter.start();
  }

  /**
   * Processes the color string. Returns null if empty.
   * 
   * @param colorStr	the string to process
   * @return		the processed color or null
   */
  protected Color getColor(String colorStr) {
    Color	result;
    
    result = null;
    
    if ((colorStr != null) && (colorStr.length() > 0))
      result = VisualizeUtils.processColour(colorStr, result);
    
    return result;
  }
  
  /**
   * Performs some initialization.
   */
  protected void initialize() {
    Properties	props;
    
    try {
      props = Utils.readProperties(PROPERTIES_FILE);
    }
    catch (Exception e) {
      e.printStackTrace();
      props = new Properties();
    }
    
    m_FontColor       = getColor(props.getProperty("FontColor", ""));
    m_BackgroundColor = getColor(props.getProperty("BackgroundColor", ""));
    m_NodeColor       = getColor(props.getProperty("NodeColor", ""));
    m_LineColor       = getColor(props.getProperty("LineColor", ""));
    m_ZoomBoxColor    = getColor(props.getProperty("ZoomBoxColor", ""));
    m_ZoomBoxXORColor = getColor(props.getProperty("ZoomBoxXORColor", ""));
    m_ShowBorder      = Boolean.parseBoolean(props.getProperty("ShowBorder", "true"));
  }

  /**
   * Fits the tree to the current screen size. Call this after
   * window has been created to get the entrire tree to be in view
   * upon launch.
   */
  public void fitToScreen() {

    getScreenFit(m_viewPos, m_viewSize);
    repaint();
  }

  /**
   * Calculates the dimensions needed to fit the entire tree into view.
   */
  private void getScreenFit(Dimension np, Dimension ns) {

    int leftmost = 1000000, rightmost = -1000000;
    int leftCenter = 1000000, rightCenter = -1000000, rightNode = 0;
    int highest = -1000000, highTop = -1000000;
    for (int noa = 0; noa < m_numNodes; noa++) {
      calcScreenCoords(noa);
      if (m_nodes[noa].m_center - m_nodes[noa].m_side < leftmost) {
	leftmost = m_nodes[noa].m_center - m_nodes[noa].m_side;
      }
      if (m_nodes[noa].m_center < leftCenter) {
	leftCenter = m_nodes[noa].m_center;
      }
      
      if (m_nodes[noa].m_center + m_nodes[noa].m_side > rightmost) {
	rightmost = m_nodes[noa].m_center + m_nodes[noa].m_side;	  
      }
      if (m_nodes[noa].m_center > rightCenter) {
	rightCenter = m_nodes[noa].m_center;
	rightNode = noa;
      }
	if (m_nodes[noa].m_top + m_nodes[noa].m_height > highest) {
	  highest = m_nodes[noa].m_top + m_nodes[noa].m_height;
	}
	if (m_nodes[noa].m_top > highTop) {
	  highTop = m_nodes[noa].m_top;
	}
    }
    
    ns.width = getWidth();
    ns.width -= leftCenter - leftmost + rightmost - rightCenter + 30;
    ns.height = getHeight() - highest + highTop - 40;
    
    if (m_nodes[rightNode].m_node.getCenter() != 0 
	&& leftCenter != rightCenter) {
	ns.width /= m_nodes[rightNode].m_node.getCenter();
    }
    if (ns.width < 10)
      {
	ns.width = 10;
      }
    if (ns.height < 10)
      {
	ns.height = 10;
      }
    
    np.width = (leftCenter - leftmost + rightmost - rightCenter) / 2 + 15;
    np.height = (highest - highTop) / 2 + 20;
  }

  /**
   * Performs the action associated with the ActionEvent.
   *
   * @param e the action event.
   */
  public void actionPerformed(ActionEvent e) {
    
    //JMenuItem m = (JMenuItem)e.getSource();
    
    if (e.getActionCommand() == null) {
      if (m_scaling == 0) {
	repaint();
      }
      else {
	animateScaling(m_nViewPos, m_nViewSize, m_scaling);
      }
    }
    else if (e.getActionCommand().equals("Fit to Screen")) {
      
      Dimension np = new Dimension();
      Dimension ns = new Dimension();

      getScreenFit(np, ns);

      animateScaling(np, ns, 10);
      
    }
    else if (e.getActionCommand().equals("Center on Top Node")) {
      
      int tpx = (int)(m_topNode.getCenter() * m_viewSize.width);   //calculate
      //the top nodes postion but don't adjust for where 
      int tpy = (int)(m_topNode.getTop() * m_viewSize.height);     //view is
      
      
      
      Dimension np = new Dimension(getSize().width / 2 - tpx, 
				   getSize().width / 6 - tpy);
      
      animateScaling(np, m_viewSize, 10);
      
    }
    else if (e.getActionCommand().equals("Auto Scale")) {
      autoScale();  //this will figure the best scale value 
      //keep the focus on the middle of the screen and call animate
    }
    else if (e.getActionCommand().equals("Visualize The Node")) {
      //send the node data to the visualizer 
      if (m_focusNode >= 0) {
	Instances inst;
	if ((inst = m_nodes[m_focusNode].m_node.getInstances()) != null) {
	  VisualizePanel pan = new VisualizePanel();
	  pan.setInstances(inst);
	  JFrame nf = new JFrame();
	  nf.setSize(400, 300);
	  nf.getContentPane().add(pan);
	  nf.setVisible(true);
	}
	else {
	  JOptionPane.showMessageDialog(this, "Sorry, there is no " + 
					"available Instances data for " +
					"this Node.", "Sorry!",
					JOptionPane.WARNING_MESSAGE); 
	}
      }
      else {
	JOptionPane.showMessageDialog(this, "Error, there is no " + 
				      "selected Node to perform " +
				      "this operation on.", "Error!",
				      JOptionPane.ERROR_MESSAGE); 
      }
    }
    else if (e.getActionCommand().equals("Create Child Nodes")) {
      if (m_focusNode >= 0) {
	if (m_listener != null) {
	  //then send message to the listener
	  m_listener.userCommand(new TreeDisplayEvent
	    (TreeDisplayEvent.ADD_CHILDREN, 
	     m_nodes[m_focusNode].m_node.getRefer()));
	}
	else {
	  JOptionPane.showMessageDialog(this, "Sorry, there is no " + 
					"available Decision Tree to " +
					"perform this operation on.",
					"Sorry!", 
					JOptionPane.WARNING_MESSAGE);
	}
      }
      else {
	JOptionPane.showMessageDialog(this, "Error, there is no " +
				      "selected Node to perform this " +
				      "operation on.", "Error!",
				      JOptionPane.ERROR_MESSAGE);
      }
    }
    else if (e.getActionCommand().equals("Remove Child Nodes")) {
      if (m_focusNode >= 0) {
	if (m_listener != null) {
	  //then send message to the listener
	  m_listener.userCommand(new 
	    TreeDisplayEvent(TreeDisplayEvent.REMOVE_CHILDREN, 
			     m_nodes[m_focusNode].m_node.getRefer()));
	}
	else {
	  JOptionPane.showMessageDialog(this, "Sorry, there is no " + 
					"available Decsion Tree to " +
					"perform this operation on.",
					"Sorry!", 
					JOptionPane.WARNING_MESSAGE);
	}
      }
      else {
	JOptionPane.showMessageDialog(this, "Error, there is no " +
				      "selected Node to perform this " +
				      "operation on.", "Error!",
				      JOptionPane.ERROR_MESSAGE);
      }
    }
    else if (e.getActionCommand().equals("classify_child")) {
      if (m_focusNode >= 0) {
	if (m_listener != null) {
	  //then send message to the listener
	  m_listener.userCommand(new TreeDisplayEvent
	    (TreeDisplayEvent.CLASSIFY_CHILD, 
	     m_nodes[m_focusNode].m_node.getRefer()));
	}
	else {
	  JOptionPane.showMessageDialog(this, "Sorry, there is no " + 
					"available Decsion Tree to " +
					"perform this operation on.",
					"Sorry!", 
					JOptionPane.WARNING_MESSAGE);
	}
      }
      else {
	JOptionPane.showMessageDialog(this, "Error, there is no " +
				      "selected Node to perform this " +
				      "operation on.", "Error!",
				      JOptionPane.ERROR_MESSAGE);
      }
    }
    else if (e.getActionCommand().equals("send_instances")) {
      if (m_focusNode >= 0) {
	if (m_listener != null) {
	  //then send message to the listener
	  m_listener.userCommand(new TreeDisplayEvent
	    (TreeDisplayEvent.SEND_INSTANCES, 
	     m_nodes[m_focusNode].m_node.getRefer()));
	}
	else {
	  JOptionPane.showMessageDialog(this, "Sorry, there is no " + 
					"available Decsion Tree to " +
					"perform this operation on.",
					"Sorry!", 
					JOptionPane.WARNING_MESSAGE);
	}
      }
      else {
	JOptionPane.showMessageDialog(this, "Error, there is no " +
				      "selected Node to perform this " +
				      "operation on.", "Error!",
				      JOptionPane.ERROR_MESSAGE);
      }
    }
    else if (e.getActionCommand().equals("Accept The Tree")) {
      if (m_listener != null) {
	//then send message to the listener saying that the tree is done
	m_listener.userCommand(new TreeDisplayEvent(TreeDisplayEvent.ACCEPT,
						  null));
      }
      else {
	JOptionPane.showMessageDialog(this, "Sorry, there is no " +
				      "available Decision Tree to " +
				      "perform this operation on.",
				      "Sorry!", 
				      JOptionPane.WARNING_MESSAGE);
      }
    }
  }

  /**
   * Performs the action associated with the ItemEvent.
   *
   * @param e the item event.
   */
  public void itemStateChanged(ItemEvent e)
  {
    JRadioButtonMenuItem c = (JRadioButtonMenuItem)e.getSource();
    if (c.getActionCommand().equals("Size 24")) {
      changeFontSize(24);
    }
    else if (c.getActionCommand().equals("Size 22")) {
      changeFontSize(22);
    }
    else if (c.getActionCommand().equals("Size 20")) {
      changeFontSize(20);
    }
    else if (c.getActionCommand().equals("Size 18")) {
      changeFontSize(18);
    } 
    else if (c.getActionCommand().equals("Size 16")) {
      changeFontSize(16);
    }
    else if (c.getActionCommand().equals("Size 14")) {
      changeFontSize(14);
    }
    else if (c.getActionCommand().equals("Size 12")) {
      changeFontSize(12);
    }
    else if (c.getActionCommand().equals("Size 10")) {
      changeFontSize(10);
    }
    else if (c.getActionCommand().equals("Size 8")) {
      changeFontSize(8);
    }
    else if (c.getActionCommand().equals("Size 6")) {
      changeFontSize(6);
    }
    else if (c.getActionCommand().equals("Size 4")) {
      changeFontSize(4);
    }
    else if (c.getActionCommand().equals("Size 2")) {
      changeFontSize(2);
    }
    else if (c.getActionCommand().equals("Size 1")) {
      changeFontSize(1);
    }
    else if (c.getActionCommand().equals("Hide Descendants")) {
      //focus_node.setCVisible(!c.isSelected());
      //no longer used...
    }
  }

  /**
   * Does nothing.
   * @param e the mouse event.
   */
  public void mouseClicked(MouseEvent e) {
    //if the mouse was left clicked on 
    //the node then 
    if (m_clickAvailable) {
      //determine if the click was on a node or not
      int s = -1;
      
      for (int noa = 0; noa < m_numNodes;noa++) {
	if (m_nodes[noa].m_quad == 18) {
	  //then is on the screen
	  calcScreenCoords(noa);
	  if (e.getX() <= m_nodes[noa].m_center + m_nodes[noa].m_side 
	      && e.getX() 
	      >= m_nodes[noa].m_center - m_nodes[noa].m_side &&
	      e.getY() >= m_nodes[noa].m_top && e.getY() 
	      <= m_nodes[noa].m_top + m_nodes[noa].m_height) {
	    //then it is this node that the mouse was clicked on
	    s = noa;
	  }
	  m_nodes[noa].m_top = 32000;
	}
      }
      m_focusNode = s;
      
      if (m_focusNode != -1) {
	if (m_listener != null) {
	  //then set this to be the selected node for editing
	  actionPerformed(new ActionEvent(this, 32000, "Create Child Nodes"));
	  
	}
	else {
	  //then open a visualize to display this nodes instances if possible
	  actionPerformed(new ActionEvent(this, 32000, "Visualize The Node"));
	}
      }
    }
  }
  
  /**
   * Determines what action the user wants to perform.
   *
   * @param e the mouse event.
   */  
  public void mousePressed(MouseEvent e) {
    m_frameLimiter.setRepeats(true);
    if ((e.getModifiers() & e.BUTTON1_MASK) != 0 && !e.isAltDown() && 
	m_mouseState == 0 
	&& m_scaling == 0) {
      //then the left mouse button has been pressed
      //check for modifiers
      
      if (((e.getModifiers() & e.CTRL_MASK) != 0) && ((e.getModifiers() & e.SHIFT_MASK) == 0)) {
	//then is in zoom out mode
	m_mouseState = 2;
      }
      else if (((e.getModifiers() & e.SHIFT_MASK) != 0) && ((e.getModifiers() & e.CTRL_MASK) == 0)) {
	//then is in zoom mode
	//note if both are pressed default action is to zoom out
	m_oldMousePos.width = e.getX();
	m_oldMousePos.height = e.getY();
	m_newMousePos.width = e.getX();
	m_newMousePos.height = e.getY();
	m_mouseState = 3;
	
	Graphics g = getGraphics();
	if (m_ZoomBoxColor == null)
	  g.setColor(Color.black);
	else
	  g.setColor(m_ZoomBoxColor);
	if (m_ZoomBoxXORColor == null)
	  g.setXORMode(Color.white);
	else
	  g.setXORMode(m_ZoomBoxXORColor);
	g.drawRect(m_oldMousePos.width, m_oldMousePos.height,
		   m_newMousePos.width - m_oldMousePos.width, 
		   m_newMousePos.height - m_oldMousePos.height);
	g.dispose();
      }
      else {
	//no modifiers drag area around
	m_oldMousePos.width = e.getX();
	m_oldMousePos.height = e.getY();
	m_newMousePos.width = e.getX();
	m_newMousePos.height = e.getY();
	m_mouseState = 1;
	m_frameLimiter.start();
      }
      
    }
    // pop up save dialog explicitly (is somehow overridden...)
    else if ( (e.getButton() == MouseEvent.BUTTON1) && e.isAltDown() && e.isShiftDown() && !e.isControlDown() ) {
      saveComponent();
    }
    else if (m_mouseState == 0 && m_scaling == 0) {
      //either middle or right mouse button pushed
      //determine menu to use
      
    }
  }
  
  /**
   * Performs the final stages of what the user wants to perform.
   *
   * @param e the mouse event.
   */
  public void mouseReleased(MouseEvent e) {
    if (m_mouseState == 1) {
      //this is used by mouseClicked to determine if it is alright to do 
      //something
      m_clickAvailable = true;
	//note that a standard click with the left mouse is pretty much the 
      //only safe input left to be assigned anything.
    }
    else {
      m_clickAvailable = false;
    }
    if (m_mouseState == 2 && mouseInBounds(e)) {
      //then zoom out;
      m_mouseState = 0;
      Dimension ns = new Dimension(m_viewSize.width / 2, m_viewSize.height 
				   / 2);
      if (ns.width < 10) {
	ns.width = 10;
      }
      if (ns.height < 10) {
	ns.height = 10;
      }
      
      Dimension d = getSize();
      Dimension np = new Dimension((int)(d.width / 2 
					 - ((double)d.width / 2 
					    - m_viewPos.width) / 2),
				   (int)(d.height / 2 
					 - ((double)d.height / 2
					    - m_viewPos.height) / 2));
      
      animateScaling(np, ns, 10);
      
      //view_pos.width += view_size.width / 2;
      //view_pos.height += view_size.height / 2;
      
    }
    else if (m_mouseState == 3) {
      //then zoom in
      m_mouseState = 0;
      Graphics g = getGraphics();
      if (m_ZoomBoxColor == null)
	g.setColor(Color.black);
      else
	g.setColor(m_ZoomBoxColor);
      if (m_ZoomBoxXORColor == null)
	g.setXORMode(Color.white);
      else
	g.setXORMode(m_ZoomBoxXORColor);
      g.drawRect(m_oldMousePos.width, m_oldMousePos.height, 
		 m_newMousePos.width - m_oldMousePos.width, 
		 m_newMousePos.height - m_oldMousePos.height);
      g.dispose();
      
      
      
      int cw = m_newMousePos.width - m_oldMousePos.width;
      int ch = m_newMousePos.height - m_oldMousePos.height;
      if (cw >= 1 && ch >= 1) {
	if (mouseInBounds(e) && 
	    (getSize().width / cw) <= 6 &&
	    (getSize().height / ch) <= 6) {
	  
	  //now calculate new position and size
	  Dimension ns = new Dimension();
	  Dimension np = new Dimension();
	  double nvsw = getSize().width / (double)(cw);
	  double nvsh = getSize().height / (double)(ch);
	  np.width = (int)((m_oldMousePos.width - m_viewPos.width) * -nvsw);
	  np.height = (int)((m_oldMousePos.height - m_viewPos.height) * -nvsh);
	  ns.width = (int)(m_viewSize.width * nvsw);
	  ns.height = (int)(m_viewSize.height * nvsh);
	  
	  animateScaling(np, ns, 10);
	  
	  
	}
      }
    }
    else if (m_mouseState == 0 && m_scaling == 0) {
      //menu
      m_mouseState = 0;
      setFont(new Font("A Name", 0, 12));
      //determine if the click was on a node or not
      int s = -1;
      
      for (int noa = 0; noa < m_numNodes;noa++) {
	if (m_nodes[noa].m_quad == 18) {
	  //then is on the screen
	  calcScreenCoords(noa);
	  if (e.getX() <= m_nodes[noa].m_center + m_nodes[noa].m_side 
	      && e.getX() 
	      >= m_nodes[noa].m_center - m_nodes[noa].m_side &&
	      e.getY() >= m_nodes[noa].m_top && e.getY() 
	      <= m_nodes[noa].m_top + m_nodes[noa].m_height) {
	    //then it is this node that the mouse was clicked on
	    s = noa;
	  }
	  m_nodes[noa].m_top = 32000;
	}
      }
      if (s == -1) {
	//the mouse wasn't clicked on a node
	m_winMenu.show(this,e.getX(),e.getY());
      }
      else {
	//the mouse was clicked on a node
	m_focusNode = s;
	m_nodeMenu.show(this, e.getX(), e.getY());
	
      }
      setFont(m_currentFont);
    }
    else if (m_mouseState == 1) {
      //dragging
      m_mouseState = 0;
      m_frameLimiter.stop();
      repaint();
    }
    
  }
  
  /**
   * Checks to see if the coordinates of the mouse lie on this JPanel.
   *
   * @param e the mouse event.
   * @return true if the mouse lies on this JPanel. 
   */
  private boolean mouseInBounds(MouseEvent e) {
    //this returns true if the mouse is currently over the canvas otherwise 
    //false
    
    if (e.getX() < 0 || e.getY() < 0 || e.getX() > getSize().width 
	|| e.getY() > getSize().height) {
      return false;
    }
    return true;
  }
  
  /**
   * Performs intermediate updates to what the user wishes to do.
   *
   * @param e the mouse event.
   */
  public void mouseDragged(MouseEvent e) {
    //use mouse state to determine what to do to the view of the tree
    
    if (m_mouseState == 1) {
      //then dragging view
      m_oldMousePos.width = m_newMousePos.width;
      m_oldMousePos.height = m_newMousePos.height;
      m_newMousePos.width = e.getX();
      m_newMousePos.height = e.getY();
      m_viewPos.width += m_newMousePos.width - m_oldMousePos.width;
      m_viewPos.height += m_newMousePos.height - m_oldMousePos.height;
      
      
    }
    else if (m_mouseState == 3) {
      //then zoom box being created
      //redraw the zoom box
      Graphics g = getGraphics();
      if (m_ZoomBoxColor == null)
	g.setColor(Color.black);
      else
	g.setColor(m_ZoomBoxColor);
      if (m_ZoomBoxXORColor == null)
	g.setXORMode(Color.white);
      else
	g.setXORMode(m_ZoomBoxXORColor);
      g.drawRect(m_oldMousePos.width, m_oldMousePos.height,
		 m_newMousePos.width - m_oldMousePos.width, 
		 m_newMousePos.height - m_oldMousePos.height);
      
      m_newMousePos.width = e.getX();
      m_newMousePos.height = e.getY();
      
      g.drawRect(m_oldMousePos.width, m_oldMousePos.height,
		 m_newMousePos.width - m_oldMousePos.width, 
		 m_newMousePos.height - m_oldMousePos.height);
      g.dispose();
    }
    
    
  }
  
  /**
   * Does nothing.
   *
   * @param e the mouse event.
   */
  public void mouseMoved(MouseEvent e) {
  }

  /**
   * Does nothing.
   *
   * @param e the mouse event.
   */
  public void mouseEntered(MouseEvent e) {
  }
  
  /**
   * Does nothing.
   * 
   * @param e the mouse event.
   */
  public void mouseExited(MouseEvent e) {
  }

  /**
   * Set the highlight for the node with the given id
   * @param id the id of the node to set the highlight for
   */
  public void setHighlight(String id) {
    //set the highlight for the node with the given id
    
    for (int noa = 0; noa < m_numNodes; noa++) {
      if (id.equals(m_nodes[noa].m_node.getRefer())) {
	//then highlight this node
	m_highlightNode = noa;
      }
    }
    //System.out.println("ahuh " + highlight_node + " " + 
    //nodes[0].node.getRefer());
    repaint();
    
  }
  
  /**
   * Updates the screen contents.
   *
   * @param g the drawing surface.
   */
  public void paintComponent(Graphics g) {
    Color oldBackground = ((Graphics2D) g).getBackground();
    if (m_BackgroundColor != null)
      ((Graphics2D) g).setBackground(m_BackgroundColor);
    g.clearRect(0, 0, getSize().width, getSize().height);
    ((Graphics2D) g).setBackground(oldBackground);
    g.setClip(3, 7, getWidth() - 6, getHeight() - 10);
    painter(g);
    g.setClip(0, 0, getWidth(), getHeight());
    
  }

  /**
   * Draws the tree to the graphics context
   *
   * @param g the drawing surface.
   */
  private void painter(Graphics g) {
    //I have moved what would normally be in the paintComponent 
    //function to here 
    //for now so that if I do in fact need to do double 
    //buffering or the like it will be easier

    //this will go through the table of edges and draw the edge if it deems the
    //two nodes attached to it could cause it to cut the screen or be on it.
    
    //in the process flagging all nodes so that they can quickly be put to the
    //screen if they lie on it

    //I do it in this order because in some circumstances I have seen a line
    //cut through a node , to make things look better the line will
    //be drawn under the node


    //converting the screen edges to the node scale so that they 
    //can be positioned relative to the screen
    //note I give a buffer around the edges of the screen.
    
    //when seeing
    //if a node is on screen I only bother to check the nodes top centre 
    //if it has large enough size it may still fall onto the screen
    double left_clip = (double)(-m_viewPos.width - 50) / m_viewSize.width;
    double right_clip = (double)(getSize().width - m_viewPos.width + 50) / 
      m_viewSize.width;
    double top_clip = (double)(-m_viewPos.height - 50) / m_viewSize.height;
    double bottom_clip = (double)(getSize().height - m_viewPos.height + 50) / 
      m_viewSize.height;
  

    
    //  12 10  9           //the quadrants
    //  20 18 17
    //  36 34 33


    //first the edges must be rendered

    Edge e;
    Node r,s;
    double ncent,ntop;    

    int row = 0, col = 0, pq, cq;
    for (int noa = 0 ; noa < m_numNodes ; noa++) {
      r = m_nodes[noa].m_node;
      if (m_nodes[noa].m_change) {
	//then recalc row component of quadrant
	ntop = r.getTop();
	if (ntop < top_clip) {
	  row = 8;
	}
	else if (ntop > bottom_clip) {
	  row = 32;
	}
	else {
	  row = 16;
	}
      }
      
      //calc the column the node falls in for the quadrant
      ncent = r.getCenter();
      if (ncent < left_clip) {
	col = 4;
      }
      else if (ncent > right_clip) {
	col = 1;
      }
      else {
	col = 2;
      }
      
      m_nodes[noa].m_quad = row | col;
      
      if (m_nodes[noa].m_parent >= 0) {
	//this will draw the edge if it should be drawn 
	//It will do this by eliminating all edges that definitely won't enter
	//the screen and then draw the rest
	
	pq = m_nodes[m_edges[m_nodes[noa].m_parent].m_parent].m_quad;
	cq = m_nodes[noa].m_quad;
	
	//note that this will need to be altered if more than 1 parent exists
	if ((cq & 8) == 8) {
	  //then child exists above screen
	}
	else if ((pq & 32) == 32) {
	  //then parent exists below screen
	}
	else if ((cq & 4) == 4 && (pq & 4) == 4) {
	  //then both child and parent exist to the left of the screen
	}
	else if ((cq & 1) == 1 && (pq & 1) == 1) {
	  //then both child and parent exist to the right of the screen
	}
	else {
	  //then draw the line
	  drawLine(m_nodes[noa].m_parent, g);
	}
      }
      
      //now draw the nodes
    }
    
    for (int noa = 0 ;noa < m_numNodes; noa++) {
      if (m_nodes[noa].m_quad == 18) {
	//then the node is on the screen , draw it
	drawNode(noa, g);
      }
    }
    
    if (m_highlightNode >= 0 && m_highlightNode < m_numNodes) {
      //then draw outline
      if (m_nodes[m_highlightNode].m_quad == 18) {
	Color acol;
	if (m_NodeColor == null)
	  acol = m_nodes[m_highlightNode].m_node.getColor();
	else
	  acol = m_NodeColor;
	g.setColor(new Color((acol.getRed() + 125) % 256, 
			     (acol.getGreen() + 125) % 256, 
			     (acol.getBlue() + 125) % 256));
	//g.setXORMode(Color.white);
	if (m_nodes[m_highlightNode].m_node.getShape() == 1) {
	  g.drawRect(m_nodes[m_highlightNode].m_center 
		     - m_nodes[m_highlightNode].m_side,
		     m_nodes[m_highlightNode].m_top, 
		     m_nodes[m_highlightNode].m_width, 
		     m_nodes[m_highlightNode].m_height);
	  
	  g.drawRect(m_nodes[m_highlightNode].m_center 
		     - m_nodes[m_highlightNode].m_side 
		     + 1, m_nodes[m_highlightNode].m_top + 1,
		     m_nodes[m_highlightNode].m_width - 2, 
		     m_nodes[m_highlightNode].m_height - 2);
	}
	else if (m_nodes[m_highlightNode].m_node.getShape() == 2) {
	  g.drawOval(m_nodes[m_highlightNode].m_center 
		     - m_nodes[m_highlightNode].m_side,
		     m_nodes[m_highlightNode].m_top, 
		     m_nodes[m_highlightNode].m_width, 
		     m_nodes[m_highlightNode].m_height);
	  
	  g.drawOval(m_nodes[m_highlightNode].m_center 
		     - m_nodes[m_highlightNode].m_side 
		     + 1, m_nodes[m_highlightNode].m_top + 1,
		     m_nodes[m_highlightNode].m_width - 2, 
		     m_nodes[m_highlightNode].m_height - 2);
	}
      }
    }
    
    for (int noa = 0;noa < m_numNodes;noa++) {
      //this resets the coords so that next time a refresh occurs 
      //they don't accidentally get used
      //I will use 32000 to signify that they are invalid, even if this 
      //coordinate occurs it doesn't
      //matter as it is only for the sake of the caching
      
      m_nodes[noa].m_top = 32000;
    }
  }
  
  /**
   * Determines the attributes of the node and draws it.
   *
   * @param n A subscript identifying the node in <i>nodes</i> array
   * @param g The drawing surface
   */
  private void drawNode(int n, Graphics g) {
    //this will draw a node and then print text on it
    
    if (m_NodeColor == null)
      g.setColor(m_nodes[n].m_node.getColor());
    else
      g.setColor(m_NodeColor);
    g.setPaintMode();
    calcScreenCoords(n);
    int x = m_nodes[n].m_center - m_nodes[n].m_side;
    int y = m_nodes[n].m_top;
    if (m_nodes[n].m_node.getShape() == 1) {
      g.fill3DRect(x, y, m_nodes[n].m_width, m_nodes[n].m_height, true);
      drawText(x, y, n, false, g);
      
    }
    else if (m_nodes[n].m_node.getShape() == 2) {
      
      g.fillOval(x, y, m_nodes[n].m_width, m_nodes[n].m_height);
      drawText(x, y + (int)(m_nodes[n].m_height * .15), n, false, g);
    }
  }
  
  /**
   * Determines the attributes of the edge and draws it.
   *
   * @param e A subscript identifying the edge in <i>edges</i> array.
   * @param g The drawing surface.
   */
  private void drawLine(int e, Graphics g) {
    //this will draw a line taking in the edge number and then getting 
    //the nodes subscript for the parent and child entries
    
    //this will draw a line that has been broken in the middle 
    //for the edge text to be displayed 
    //if applicable
    
    //first convert both parent and child node coords to screen coords
    int p = m_edges[e].m_parent;
    int c = m_edges[e].m_child;
    calcScreenCoords(c);
    calcScreenCoords(p);
    
    if (m_LineColor == null)
      g.setColor(Color.black);
    else
      g.setColor(m_LineColor);
    g.setPaintMode();
    
    if (m_currentFont.getSize() < 2) {
      //text to small to bother cutting the edge
      g.drawLine(m_nodes[p].m_center, m_nodes[p].m_top + m_nodes[p].m_height, 
		 m_nodes[c].m_center, m_nodes[c].m_top); 
      
    }
    else {
      //find where to cut the edge to insert text
      int e_width = m_nodes[c].m_center - m_nodes[p].m_center;
      int e_height = m_nodes[c].m_top - (m_nodes[p].m_top 
					 + m_nodes[p].m_height);
      int e_width2 = e_width / 2;
      int e_height2 = e_height / 2;
      int e_centerx = m_nodes[p].m_center + e_width2;
      int e_centery = m_nodes[p].m_top + m_nodes[p].m_height + e_height2;
      int e_offset = m_edges[e].m_tb;
      
      int tmp = (int)(((double)e_width / e_height) * 
		      (e_height2 - e_offset)) + m_nodes[p].m_center;
      //System.out.println(edges[e].m_height);
      
      //draw text now
      
      drawText(e_centerx - m_edges[e].m_side, e_centery - e_offset, e, true
	       , g);
      
      
      if (tmp > (e_centerx - m_edges[e].m_side) && tmp 
	  < (e_centerx + m_edges[e].m_side)) {
	//then cut line on top and bottom of text
	g.drawLine(m_nodes[p].m_center, m_nodes[p].m_top + m_nodes[p].m_height
		   , tmp, e_centery - e_offset);  //first segment
	g.drawLine(e_centerx * 2 - tmp, e_centery + e_offset, 
		   m_nodes[c].m_center, m_nodes[c].m_top);    //second segment
      }
      else {
	e_offset = m_edges[e].m_side;
	if (e_width < 0) {
	  e_offset *= -1;   //adjusting for direction which could otherwise 
	  //screw up the calculation
	}
	tmp = (int)(((double)e_height / e_width) * (e_width2 - e_offset)) + 
	  m_nodes[p].m_top + m_nodes[p].m_height;
	
	g.drawLine(m_nodes[p].m_center, m_nodes[p].m_top + m_nodes[p].m_height
		   , e_centerx - e_offset, tmp);   //first segment
	g.drawLine(e_centerx + e_offset, e_centery * 2 - tmp, 
		   m_nodes[c].m_center, m_nodes[c].m_top);  //second segment
	
      }
    }
    //System.out.println("here" + nodes[p].center);
  }

  /**
   * Draws the text for either an Edge or a Node.
   *
   * @param x1 the left side of the text area.
   * @param y1 the top of the text area.
   * @param s A subscript identifying either a Node or Edge.
   * @param e_or_n Distinguishes whether it is a node or edge.
   * @param g The drawing surface.
   */
  private void drawText(int x1, int y1, int s, boolean e_or_n, Graphics g) {
    //this function will take in the rectangle that the text should be 
    //drawn in as well as the subscript
    //for either the edge or node and a boolean variable to tell which
    
    // backup color
    Color oldColor = g.getColor();

    g.setPaintMode();
    if (m_FontColor == null)
      g.setColor(Color.black);
    else
      g.setColor(m_FontColor);
    String st;
    if (e_or_n) {
      //then paint for edge
      Edge e = m_edges[s].m_edge;
      for (int noa = 0;(st = e.getLine(noa)) != null; noa++) {
	g.drawString(st, (m_edges[s].m_width - m_fontSize.stringWidth(st)) / 2 
		     + x1,
		     y1 + (noa + 1) * m_fontSize.getHeight()); 
      }
    }
    else {
      //then paint for node
      Node e = m_nodes[s].m_node;
      for (int noa = 0;(st = e.getLine(noa)) != null; noa++) {
	g.drawString(st, (m_nodes[s].m_width - m_fontSize.stringWidth(st)) / 2 
		     + x1,
		     y1 + (noa + 1) * m_fontSize.getHeight()); 
      }
    }
    
    // restore color
    g.setColor(oldColor);
  }
  
  /**
   * Converts the internal coordinates of the node found from <i>n</i>
   * and converts them to the actual screen coordinates.
   *
   * @param n A subscript identifying the Node.
   */
  private void calcScreenCoords(int n) {
    //this converts the coordinate system the Node uses into screen coordinates
    // System.out.println(n + " " + view_pos.height + " " + 
    //nodes[n].node.getCenter());
    if (m_nodes[n].m_top == 32000) {
      m_nodes[n].m_top = ((int)(m_nodes[n].m_node.getTop() 
				* m_viewSize.height)) + m_viewPos.height;
      m_nodes[n].m_center = ((int)(m_nodes[n].m_node.getCenter() 
				   * m_viewSize.width)) + m_viewPos.width;
    }
  }

  /**
   * This Calculates the minimum size of the tree which will prevent any text
   * overlapping and make it readable, and then set the size of the tree to 
   * this.
   */
  private void autoScale() {
    //this function will determine the smallest scale value that keeps the text
    //from overlapping
    //it will leave the view centered
    
    int dist;
    Node ln,rn;
    Dimension temp = new Dimension(10, 10);
    
    if (m_numNodes <= 1) {
      return;
    }
    
    //calc height needed by first node
    dist = (m_nodes[0].m_height + 40) * m_numLevels;
    if (dist > temp.height) {
      temp.height = dist;
    }
    
    for (int noa = 0;noa < m_numNodes - 1;noa++) {
      calcScreenCoords(noa);  
      calcScreenCoords(noa+1);
      if (m_nodes[noa+1].m_change) {
	//then on a new level so don't check width this time round
      }
      else {
	
	dist = m_nodes[noa+1].m_center - m_nodes[noa].m_center; 
	//the distance between the node centers, along horiz
	if (dist <= 0) {
	  dist = 1;
	}
	dist = ((6 + m_nodes[noa].m_side + m_nodes[noa+1].m_side) 
		* m_viewSize.width) / dist; //calc optimal size for width
	
	if (dist > temp.width) {
	  
	  temp.width = dist;
	}
      }
      //now calc.. minimun hieght needed by nodes
      
      dist = (m_nodes[noa+1].m_height + 40) * m_numLevels;
      if (dist > temp.height) {
	
	temp.height = dist;
      }
    }
    
    int y1, y2, xa, xb;
    
    y1 = m_nodes[m_edges[0].m_parent].m_top;
    y2 = m_nodes[m_edges[0].m_child].m_top;
    
    dist = y2 - y1;
    if (dist <= 0) {
      dist = 1;
    }
    dist = ((60 + m_edges[0].m_height + m_nodes[m_edges[0].m_parent].m_height) 
	    * m_viewSize.height) / dist;
    if (dist > temp.height) {
      
      temp.height = dist;
    }
    
    for (int noa = 0;noa < m_numNodes - 2; noa++) {
      //check the edges now
      if (m_nodes[m_edges[noa+1].m_child].m_change) {
	//then edge is on a different level , so skip this one
      }
      else {
	//calc the width requirements of this pair of edges
	
	xa = m_nodes[m_edges[noa].m_child].m_center 
	  - m_nodes[m_edges[noa].m_parent].m_center;
	xa /= 2;
	xa += m_nodes[m_edges[noa].m_parent].m_center;
	
	xb = m_nodes[m_edges[noa+1].m_child].m_center - 
	  m_nodes[m_edges[noa+1].m_parent].m_center;
	xb /= 2;
	xb += m_nodes[m_edges[noa+1].m_parent].m_center;
	
	dist = xb - xa;
	if (dist <= 0) {
	  dist = 1;
	}
	dist = ((12 + m_edges[noa].m_side + m_edges[noa+1].m_side) 
		* m_viewSize.width) 
	  / dist;
	if (dist > temp.width) {
	  
	  temp.width = dist;
	}
      }
      //now calc height need by the edges
      y1 = m_nodes[m_edges[noa+1].m_parent].m_top;
      y2 = m_nodes[m_edges[noa+1].m_child].m_top;
      
      dist = y2 - y1;
      if (dist <= 0) {
	
	dist = 1;
      }
      dist = ((60 + m_edges[noa+1].m_height 
	       + m_nodes[m_edges[noa+1].m_parent].m_height) 
	      * m_viewSize.height) / dist;
      
      if (dist > temp.height) {
	
	temp.height = dist;
      }
    }

    Dimension e = getSize();
    
    Dimension np = new Dimension();
    np.width = (int)(e.width / 2 -  (((double)e.width / 2) - m_viewPos.width) /
		     ((double)m_viewSize.width) * (double)temp.width);
    np.height = (int)(e.height / 2 -  (((double)e.height / 2) - 
				       m_viewPos.height) /       
		      ((double)m_viewSize.height) * (double)temp.height);
    //animate_scaling(c_size,c_pos,25);
    
    for (int noa = 0;noa < m_numNodes;noa++) {
      //this resets the coords so that next time a refresh occurs they don't 
      //accidentally get used
      //I will use 32000 to signify that they are invalid, even if this 
      //coordinate occurs it doesn't
      //matter as it is only for the sake of the caching
      
      m_nodes[noa].m_top = 32000;
      
    }
    animateScaling(np, temp, 10);
  }

  /**
   * This will increment the size and position of the tree towards the 
   * desired size and position
   * a little (depending on the value of <i>frames</i>) everytime it is called.
   *
   * @param n_pos The final position of the tree wanted.
   * @param n_size The final size of the tree wanted.
   * @param frames The number of frames that shall occur before the final 
   * size and pos is reached.
   */
  private void animateScaling(Dimension n_pos,Dimension n_size,int frames) {
    //this function will take new size and position coords , and incrementally
    //scale the view to these
    //since I will be tying it in with the framelimiter I will simply call 
    //this function and increment it once
    //I will have to use a global variable since I am doing it proportionally

    if (frames == 0) {
      System.out.println("the timer didn't end in time");
      m_scaling = 0;
    }
    else {
      if (m_scaling == 0) {
	//new animate session
	//start timer and set scaling
	m_frameLimiter.start();
	m_nViewPos.width = n_pos.width;
	m_nViewPos.height = n_pos.height;
	m_nViewSize.width = n_size.width;
	m_nViewSize.height = n_size.height;
	
	m_scaling = frames;
      }
      
      int s_w = (n_size.width - m_viewSize.width) / frames;
      int s_h = (n_size.height - m_viewSize.height) / frames;
      int p_w = (n_pos.width - m_viewPos.width) / frames;
      int p_h = (n_pos.height - m_viewPos.height) / frames;
      
      m_viewSize.width += s_w;
      m_viewSize.height += s_h;
      
      m_viewPos.width += p_w;
      m_viewPos.height += p_h;
      
      repaint();
      
      m_scaling--;
      if (m_scaling == 0) {
	//all done 
	m_frameLimiter.stop();
      }
    }
  }
  
  /**
   * This will change the font size for displaying the tree to the one 
   * specified.
   *
   * @param s The new pointsize of the font.
   */
  private void changeFontSize(int s) {
    //this will set up the new font that shall be used
    //it will also recalculate the size of the nodes as these will change as 
    //a result of 
    //the new font size
    setFont(m_currentFont = new Font("A Name", 0, s));             

    m_fontSize = getFontMetrics(getFont());
    
    Dimension d;

    for (int noa = 0; noa < m_numNodes; noa++) {
      //this will set the size info for each node and edge
      
      d = m_nodes[noa].m_node.stringSize(m_fontSize);
      
      if (m_nodes[noa].m_node.getShape() == 1) {
	m_nodes[noa].m_height = d.height + 10;
	m_nodes[noa].m_width = d.width + 8;
	m_nodes[noa].m_side = m_nodes[noa].m_width / 2;
      }
      else if (m_nodes[noa].m_node.getShape() == 2) {
	m_nodes[noa].m_height = (int)((d.height + 2) * 1.6);
	m_nodes[noa].m_width = (int)((d.width + 2) * 1.6);
	m_nodes[noa].m_side = m_nodes[noa].m_width / 2;
      }
      
      if (noa < m_numNodes - 1) {
	//this will do the same for edges
	
	d = m_edges[noa].m_edge.stringSize(m_fontSize);
	
	m_edges[noa].m_height =  d.height + 8;
	m_edges[noa].m_width = d.width + 8;
	m_edges[noa].m_side = m_edges[noa].m_width / 2;
	m_edges[noa].m_tb = m_edges[noa].m_height / 2;
      }
    }
  }

  /**
   * This will fill two arrays with the Nodes and Edges from the tree
   * into a particular order.
   *
   * @param t The top Node of the tree.
   * @param l An array that has already been allocated, to be filled.
   * @param k An array that has already been allocated, to be filled.
   */
  private void arrayFill(Node t, NodeInfo[] l, EdgeInfo[] k) {
    
    //this will take the top node and the array to fill
    //it will go through the tree structure and and fill the array with the 
    //nodes 
    //from top to bottom left to right
    
    //note I do not believe this function will be able to deal with multiple 
    //parents

    if (t == null || l == null) {
      System.exit(1);      //this is just a preliminary safety check 
      //(i shouldn' need it)
    }
    
    Edge e;
    Node r,s;
    l[0] = new NodeInfo();
    l[0].m_node = t;
    l[0].m_parent = -1;
    l[0].m_change = true;
    
    int floater;       //this will point at a node that has previously been 
    //put in the list 
    //all of the children that this node has shall be put in the list , 
    //once this is done the floater shall point at the next node in the list
    //this will allow the nodes to be put into order from closest to top node
    //to furtherest from top node

    int free_space = 1; //the next empty array position

    double height = t.getTop(); //this will be used to determine if the node 
    //has a 
    //new height compared to the
    //previous one

    for (floater = 0;floater < free_space;floater++) {
      r = l[floater].m_node;
      for (int noa = 0;(e = r.getChild(noa)) != null;noa++) {
	//this loop pulls out each child of r
	
	//e points to a child edge, getTarget will return that edges child node
	s = e.getTarget();
	l[free_space] = new NodeInfo();
	l[free_space].m_node = s;
	l[free_space].m_parent = free_space - 1;
	
	k[free_space - 1] = new EdgeInfo();
	k[free_space - 1].m_edge = e;
	k[free_space - 1].m_parent = floater;
	k[free_space - 1].m_child = free_space;     //note although it's child 
	//will always have a subscript
	//of 1 more , I may not nessecarily have access to that
	//and it will need the subscr.. for multiple parents
	
	//determine if level of node has changed from previous one
	if (height != s.getTop()) {
	  l[free_space].m_change = true;
	  height = s.getTop();
	}
	else {
	  l[free_space].m_change = false;
	}
	free_space++;
      }
    }
  }

  /**
   * Internal Class for containing display information about a Node. 
   */
  private class NodeInfo {
    //this class contains a pointer to the node itself along with extra 
    //information
    //about the node used by the Displayer class

    /** The y pos of the node on screen. */
    int m_top = 32000;           //the main node coords calculated out

    /** The x pos of the node on screen. */
    int m_center;        // these coords will probably change each refresh 

    //and are the positioning coords
    //which the rest of the offsets use
    
    /** The offset to get to the left or right of the node. */
    int m_side;          //these are the screen offset for the dimensions of 

    //the node relative to the nodes 
    //internal top and center values (after they have been converted to 
    //screen coords
    /** The width of the node. */
    int m_width;

    /** The height of the node. */
    int m_height;

    /** True if the node is at the start (left) of a new level (not sibling 
     * group). */
    boolean m_change;    //this is quickly used to identify whether the node 
    //has chenged height from the
    //previous one to help speed up the calculation of what row it lies in

    /** The subscript number of the Nodes parent. */
    int m_parent;     //this is the index of the nodes parent edge in an array

    /** The rough position of the node relative to the screen. */
    int m_quad;       //what of nine quadrants is it in

    /*
      12 10  9
      20 18 17          //18 being the screen
      36 34 33          //this arrangement uses 6 bits, each bit represents a 
      row or column
    */

    /** The Node itself. */
    Node m_node;
  }

  /**
   * Internal Class for containing display information about an Edge. 
   */
  private class EdgeInfo {
    //this class contains a pointer to the edge along with all the other
    //extra info about the edge
    
    /** The parent subscript (for a Node). */
    int m_parent;            //array indexs for its two connections

    /** The child subscript (for a Node). */
    int m_child;
    

    /** The distance from the center of the text to either side. */
    int m_side;            //these are used to describe the dimensions of the 
    //text

    /** The distance from the center of the text to top or bottom. */
    int m_tb;              //tb stands for top , bottom, this is simply the 
    //distance from the middle to top bottom

    /** The width of the text. */
    int m_width;

    /** The height of the text. */
    int m_height;

    /** The Edge itself. */
    Edge m_edge;
  }

  /**
   * Main method for testing this class.
   * @param args first argument should be the name of a file that contains
   * a tree discription in dot format.
   */
  public static void main(String[] args) {
    try {
      weka.core.logging.Logger.log(weka.core.logging.Logger.Level.INFO, "Logging started");
      //put in the random data generator right here
      // this call with import java.lang gives me between 0 and 1 Math.random
      TreeBuild builder = new TreeBuild();
      Node top = null;
      NodePlace arrange = new PlaceNode2();
      //top = builder.create(new StringReader("digraph atree { top [label=\"the top\"] a [label=\"the first node\"] b [label=\"the second nodes\"] c [label=\"comes off of first\"] top->a top->b b->c }"));
      top = builder.create(new FileReader(args[0]));

      int num = top.getCount(top,0);
      //System.out.println("counter counted " + num + " nodes");
      //System.out.println("there are " + num + " nodes");
      TreeVisualizer a = new TreeVisualizer(null, top, arrange);
      a.setSize(800 ,600);
      //a.setTree(top);
      JFrame f;
      f = new JFrame();
      //a.addMouseMotionListener(a);
      //a.addMouseListener(a);
      //f.add(a);
      Container contentPane = f.getContentPane();
      contentPane.add(a);
      f.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
      f.setSize(800,600);
      f.setVisible(true);
      //f.
      //find_prop(top);
      //a.setTree(top,arrange);//,(num + 1000), num / 2 + 1000);
    }
    catch(IOException e) {
      // ignored
    }
  }
}
