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
 *    GraphVisualizer.java
 *    Copyright (C) 2003 University of Waikato, Hamilton, New Zealand
 *
 */
package weka.gui.graphvisualizer;

import weka.core.FastVector;
import weka.gui.ExtensionFileFilter;
import weka.gui.visualize.PrintablePanel;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Frame;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.LayoutManager;
import java.awt.Rectangle;
import java.awt.RenderingHints;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionAdapter;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.Reader;

import javax.swing.BorderFactory;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.JToolBar;
import javax.swing.table.AbstractTableModel;

/**
 * This class displays the graph we want to visualize. It should
 * be sufficient to use only this class in weka.gui.graphvisulizer
 * package to visualize a graph. The description of a graph should
 * be provided as a string argument using readBIF or readDOT method
 * in either XMLBIF03 or DOT format. Alternatively, an InputStream
 * in XMLBIF03 can also be provided to another variation of readBIF.
 * It would be necessary in case input is in DOT format to call the
 * layoutGraph() method to display the graph correctly after the call
 * to readDOT. It is also necessary to do so if readBIF is called and
 * the graph description doesn't have x y positions for nodes.
 * <p> The graph's data is held in two FastVectors, nodes are stored as
 * objects of GraphNode class and edges as objects of GraphEdge class.
 * <p> The graph is displayed by positioning and drawing each node
 * according to its x y position and then drawing all the edges coming
 * out of it give by its edges[][] array, the arrow heads are ofcourse
 * marked in the opposite(ie original direction) or both directions if
 * the edge is reversed or is in both directions. The graph is centered
 * if it is smaller than it's display area. The edges are drawn from the
 * bottom of the current node to the top of the node given by edges[][]
 * array in GraphNode class, to avoid edges crossing over other nodes.
 * This might need to be changed if another layout engine is added or
 * the current Hierarchical engine is updated to avoid such crossings
 * over nodes.
 *
 * @author Ashraf M. Kibriya (amk14@cs.waikato.ac.nz)
 * @version $Revision: 4723 $
 */
public class GraphVisualizer
  extends JPanel
  implements GraphConstants, LayoutCompleteEventListener {

  /** for serialization */
  private static final long serialVersionUID = -2038911085935515624L;
  
  /** Vector containing nodes */
  protected FastVector m_nodes=new FastVector();
  /** Vector containing edges */
  protected FastVector m_edges=new FastVector();
  /** The current LayoutEngine  */
  protected LayoutEngine m_le;
  /** Panel actually displaying the graph */
  protected GraphPanel m_gp;
  /** String containing graph's name */
  protected String graphID;
  
  /**
   * Save button to save the current graph in DOT or XMLBIF format.
   * The graph should be layed out again to get the original form
   * if reloaded from command line, as the formats do not allow
   * saving specific information for a properly layed out graph.
   */
  protected JButton m_jBtSave;
  
  /** path for icons */
  private final String ICONPATH = "weka/gui/graphvisualizer/icons/";
  
  private FontMetrics fm = this.getFontMetrics( this.getFont() );
  private double scale = 1;   //current zoom
  private int nodeHeight = 2*fm.getHeight(), nodeWidth = 24;
  private int paddedNodeWidth = 24+8;
  /** TextField for node's width */
  private final JTextField jTfNodeWidth = new JTextField(3);
  /** TextField for nodes height */
  private final JTextField jTfNodeHeight = new JTextField(3);
  /** Button for laying out the graph again, necessary after changing node's
   * size or some other property of the layout engine
   */
  private final JButton jBtLayout;
  /** used for setting appropriate node size */
  private int maxStringWidth=0;
  /** used when using zoomIn and zoomOut buttons */
  private int [] zoomPercents = { 10, 25, 50, 75, 100, 125, 150, 175, 200, 225,
  250, 275, 300, 350, 400, 450, 500, 550, 600, 650, 700, 800, 900, 999 };
  /** this contains the m_gp GraphPanel */
  JScrollPane m_js;
  
  /**
   * Constructor<br>
   * Sets up the gui and initializes all the other previously
   * uninitialized variables.
   */
  public GraphVisualizer() {
    m_gp = new GraphPanel();
    m_js = new JScrollPane(m_gp);
    
    //creating a new layout engine and adding this class as its listener
    // to receive layoutComplete events
    m_le=new HierarchicalBCEngine(m_nodes, m_edges, 
                                  paddedNodeWidth, nodeHeight);
    m_le.addLayoutCompleteEventListener(this);
    
    m_jBtSave = new JButton();
    java.net.URL tempURL = ClassLoader.getSystemResource(ICONPATH+"save.gif");
    if(tempURL!=null)
      m_jBtSave.setIcon(new ImageIcon(tempURL) );
    else
      System.err.println(ICONPATH+
      "save.gif not found for weka.gui.graphvisualizer.Graph");
    m_jBtSave.setToolTipText("Save Graph");
    m_jBtSave.addActionListener( new ActionListener() {
      public void actionPerformed(ActionEvent ae) {
        JFileChooser fc = new JFileChooser(System.getProperty("user.dir"));
        ExtensionFileFilter ef1 = new ExtensionFileFilter(".dot", "DOT files");
        ExtensionFileFilter ef2 = new ExtensionFileFilter(".xml",
        "XML BIF files");
        fc.addChoosableFileFilter(ef1);
        fc.addChoosableFileFilter(ef2);
        fc.setDialogTitle("Save Graph As");
        int rval = fc.showSaveDialog(GraphVisualizer.this);
        
        if (rval == JFileChooser.APPROVE_OPTION) {
          //System.out.println("Saving to file \""+
          //                   f.getAbsoluteFile().toString()+"\"");
          if(fc.getFileFilter()==ef2) {
            String filename = fc.getSelectedFile().toString();
            if(!filename.endsWith(".xml"))
              filename = filename.concat(".xml");
            BIFParser.writeXMLBIF03(filename, graphID, m_nodes, m_edges);
          }
          else {
            String filename = fc.getSelectedFile().toString();
            if(!filename.endsWith(".dot"))
              filename = filename.concat(".dot");
            DotParser.writeDOT(filename, graphID, m_nodes, m_edges);
          }
        }
      }
    });
    
    final JButton jBtZoomIn = new JButton();
    tempURL = ClassLoader.getSystemResource(ICONPATH+"zoomin.gif");
    if(tempURL!=null)
      jBtZoomIn.setIcon(new ImageIcon(tempURL) );
    else
      System.err.println(ICONPATH+
      "zoomin.gif not found for weka.gui.graphvisualizer.Graph");
    jBtZoomIn.setToolTipText("Zoom In");
    
    final JButton jBtZoomOut = new JButton();
    tempURL = ClassLoader.getSystemResource(ICONPATH+"zoomout.gif");
    if(tempURL!=null)
      jBtZoomOut.setIcon(new ImageIcon(tempURL) );
    else
      System.err.println(ICONPATH+
      "zoomout.gif not found for weka.gui.graphvisualizer.Graph");
    jBtZoomOut.setToolTipText("Zoom Out");
    
    final JTextField jTfZoom = new JTextField("100%");
    jTfZoom.setMinimumSize( jTfZoom.getPreferredSize() );
    jTfZoom.setHorizontalAlignment(JTextField.CENTER);
    jTfZoom.setToolTipText("Zoom");
    
    jTfZoom.addActionListener( new ActionListener() {
      public void actionPerformed(ActionEvent ae) {
        JTextField jt = (JTextField)ae.getSource();
        try {
          int i=-1;
          i = jt.getText().indexOf('%');
          if(i==-1)
            i = Integer.parseInt(jt.getText());
          else
            i = Integer.parseInt(jt.getText().substring(0,i));
          
          if(i<=999)
            scale = i/100D;
          
          jt.setText((int)(scale*100)+"%");
          
          if(scale>0.1){
            if(!jBtZoomOut.isEnabled())
              jBtZoomOut.setEnabled(true);
          }
          else
            jBtZoomOut.setEnabled(false);
          if(scale<9.99) {
            if(!jBtZoomIn.isEnabled())
              jBtZoomIn.setEnabled(true);
          }
          else
            jBtZoomIn.setEnabled(false);
          
          setAppropriateSize();
          //m_gp.clearBuffer();
          m_gp.repaint();
          m_gp.invalidate();
          m_js.revalidate();
        } catch(NumberFormatException ne) {
          JOptionPane.showMessageDialog(GraphVisualizer.this.getParent(),
          "Invalid integer entered for zoom.",
          "Error",
          JOptionPane.ERROR_MESSAGE);
          jt.setText((scale*100)+"%");
        }
      }
    });
    
    
    jBtZoomIn.addActionListener( new ActionListener() {
      public void actionPerformed(ActionEvent ae) {
        int i=0, s = (int)(scale*100);
        if(s<300)
          i = s/25;
        else if(s<700)
          i = 6 +  s/50;
        else
          i = 13 +s/100;
        
        if(s>=999) {
          JButton b = (JButton)ae.getSource();
          b.setEnabled(false);
          return;
        }
        else if(s>=10){
          if(i>=22) {
            JButton b = (JButton)ae.getSource();
            b.setEnabled(false);
          }
          if(s==10 && !jBtZoomOut.isEnabled())
            jBtZoomOut.setEnabled(true);
          //System.out.println("i: "+i+"Zoom is: "+zoomPercents[i+1]);
          jTfZoom.setText(zoomPercents[i+1]+"%");
          scale = zoomPercents[i+1]/100D;
        }
        else {
          if(!jBtZoomOut.isEnabled())
            jBtZoomOut.setEnabled(true);
          //System.out.println("i: "+i+"Zoom is: "+zoomPercents[0]);
          jTfZoom.setText(zoomPercents[0]+"%");
          scale = zoomPercents[0]/100D;
        }
        setAppropriateSize();
        m_gp.repaint();
        m_gp.invalidate();
        m_js.revalidate();
      }
    });
    
    
    jBtZoomOut.addActionListener( new ActionListener() {
      public void actionPerformed(ActionEvent ae) {
        int i=0, s = (int)(scale*100);
        if(s<300)
          i = (int) Math.ceil(s/25D);
        else if(s<700)
          i = 6 +  (int) Math.ceil(s/50D);
        else
          i = 13 + (int) Math.ceil(s/100D);
        
        if(s<=10) {
          JButton b = (JButton)ae.getSource();
          b.setEnabled(false);
        }
        else if(s<999) {
          if(i<=1) {
            JButton b = (JButton)ae.getSource();
            b.setEnabled(false);
          }
          //System.out.println("i: "+i+"Zoom is: "+zoomPercents[i-1]);
          jTfZoom.setText(zoomPercents[i-1]+"%");
          scale = zoomPercents[i-1]/100D;
        }
        else{
          if(!jBtZoomIn.isEnabled())
            jBtZoomIn.setEnabled(true);
          //System.out.println("i: "+i+"Zoom is: "+zoomPercents[22]);
          jTfZoom.setText(zoomPercents[22]+"%");
          scale = zoomPercents[22]/100D;
        }
        setAppropriateSize();
        m_gp.repaint();
        m_gp.invalidate();
        m_js.revalidate();
      }
    });
    
    
    //This button pops out the extra controls
    JButton jBtExtraControls = new JButton();
    tempURL = ClassLoader.getSystemResource(ICONPATH+"extra.gif");
    if(tempURL!=null)
      jBtExtraControls.setIcon(new ImageIcon(tempURL) );
    else
      System.err.println(ICONPATH+
      "extra.gif not found for weka.gui.graphvisualizer.Graph");
    jBtExtraControls.setToolTipText("Show/Hide extra controls");
    
    
    final JCheckBox jCbCustomNodeSize = new JCheckBox("Custom Node Size");
    final JLabel jLbNodeWidth = new JLabel("Width");
    final JLabel jLbNodeHeight = new JLabel("Height");
    
    jTfNodeWidth.setHorizontalAlignment(JTextField.CENTER);
    jTfNodeWidth.setText(""+nodeWidth);
    jTfNodeHeight.setHorizontalAlignment(JTextField.CENTER);
    jTfNodeHeight.setText(""+nodeHeight);
    jLbNodeWidth.setEnabled(false);
    jTfNodeWidth.setEnabled(false);
    jLbNodeHeight.setEnabled(false);
    jTfNodeHeight.setEnabled(false);
    
    jCbCustomNodeSize.addActionListener( new ActionListener() {
      public void actionPerformed(ActionEvent ae) {
        if( ((JCheckBox)ae.getSource()).isSelected() ) {
          jLbNodeWidth.setEnabled(true);
          jTfNodeWidth.setEnabled(true);
          jLbNodeHeight.setEnabled(true);
          jTfNodeHeight.setEnabled(true);
        }
        else {
          jLbNodeWidth.setEnabled(false);
          jTfNodeWidth.setEnabled(false);
          jLbNodeHeight.setEnabled(false);
          jTfNodeHeight.setEnabled(false);
          setAppropriateNodeSize();
        }
      }
    });
    
    
    jBtLayout  = new JButton("Layout Graph");
    jBtLayout.addActionListener( new ActionListener() {
      public void actionPerformed(ActionEvent ae) {
        int tmpW, tmpH;
        
        if(jCbCustomNodeSize.isSelected()) {
          try{ tmpW = Integer.parseInt(jTfNodeWidth.getText()); }
          catch(NumberFormatException ne) {
            JOptionPane.showMessageDialog(GraphVisualizer.this.getParent(),
            "Invalid integer entered for node width.",
            "Error",
            JOptionPane.ERROR_MESSAGE);
            tmpW = nodeWidth;
            jTfNodeWidth.setText(""+nodeWidth);
            
          }
          try{ tmpH = Integer.parseInt(jTfNodeHeight.getText()); }
          catch(NumberFormatException ne) {
            JOptionPane.showMessageDialog(GraphVisualizer.this.getParent(),
            "Invalid integer entered for node height.",
            "Error",
            JOptionPane.ERROR_MESSAGE);
            tmpH = nodeHeight;
            jTfNodeWidth.setText(""+nodeHeight);
          }
          
          if(tmpW!=nodeWidth || tmpH!=nodeHeight) {
            nodeWidth = tmpW; paddedNodeWidth = nodeWidth+8; nodeHeight = tmpH;
          }
        }
        JButton bt = (JButton)ae.getSource();
        bt.setEnabled(false);
        m_le.setNodeSize(paddedNodeWidth, nodeHeight);
        m_le.layoutGraph();
      }
    });
    
    
    GridBagConstraints gbc = new GridBagConstraints();
    
    final JPanel p = new JPanel(new GridBagLayout());
    gbc.gridwidth = gbc.REMAINDER;
    gbc.anchor = gbc.NORTHWEST;
    gbc.fill = gbc.NONE;
    p.add( m_le.getControlPanel(), gbc);
    gbc.gridwidth = 1;
    gbc.insets = new Insets(8,0,0,0);
    gbc.anchor = gbc.NORTHWEST;
    gbc.gridwidth = gbc.REMAINDER;
    
    p.add( jCbCustomNodeSize, gbc );
    gbc.insets = new Insets(0,0,0,0);
    gbc.gridwidth = gbc.REMAINDER;
    Container c = new Container();
    c.setLayout( new GridBagLayout() );
    gbc.gridwidth = gbc.RELATIVE;
    c.add(jLbNodeWidth, gbc);
    gbc.gridwidth = gbc.REMAINDER;
    c.add(jTfNodeWidth, gbc);
    gbc.gridwidth = gbc.RELATIVE;
    c.add(jLbNodeHeight, gbc);
    gbc.gridwidth = gbc.REMAINDER;
    c.add(jTfNodeHeight, gbc);
    gbc.fill = gbc.HORIZONTAL;
    p.add( c, gbc );
    
    gbc.anchor = gbc.NORTHWEST;
    gbc.insets = new Insets(8,0,0,0);
    gbc.fill = gbc.HORIZONTAL;
    p.add( jBtLayout, gbc );
    gbc.fill = gbc.NONE;
    p.setBorder(BorderFactory.createCompoundBorder(
    BorderFactory.createTitledBorder("ExtraControls"),
    BorderFactory.createEmptyBorder(4,4,4,4)
    ) );
    p.setPreferredSize( new Dimension(0, 0) );
    
    final JToolBar jTbTools = new JToolBar();
    jTbTools.setFloatable(false);
    jTbTools.setLayout( new GridBagLayout() );
    gbc.anchor = gbc.NORTHWEST;
    gbc.gridwidth = gbc.REMAINDER;
    gbc.insets = new Insets(0,0,0,0);
    jTbTools.add(p,gbc);
    gbc.gridwidth = 1;
    jTbTools.add(m_jBtSave, gbc);
    jTbTools.addSeparator(new Dimension(2,2));
    jTbTools.add(jBtZoomIn, gbc);
    
    gbc.fill = gbc.VERTICAL;
    gbc.weighty = 1;
    JPanel p2 = new JPanel(new BorderLayout());
    p2.setPreferredSize( jTfZoom.getPreferredSize() );
    p2.setMinimumSize( jTfZoom.getPreferredSize() );
    p2.add(jTfZoom, BorderLayout.CENTER);
    jTbTools.add(p2, gbc);
    gbc.weighty =0;
    gbc.fill = gbc.NONE;
    
    jTbTools.add(jBtZoomOut, gbc);
    jTbTools.addSeparator(new Dimension(2,2));
    jTbTools.add(jBtExtraControls, gbc);
    jTbTools.addSeparator(new Dimension(4,2));
    gbc.weightx = 1;
    gbc.fill = gbc.BOTH;
    jTbTools.add(m_le.getProgressBar(), gbc);
    
    jBtExtraControls.addActionListener( new ActionListener() {
      public void actionPerformed(ActionEvent ae) {
        Dimension d = p.getPreferredSize();
        if(d.width==0 || d.height==0) {
          LayoutManager lm = p.getLayout();
          Dimension d2 = lm.preferredLayoutSize(p);
          p.setPreferredSize(d2); jTbTools.revalidate();
          /*
          // this piece of code adds in an animation
          // for popping out the extra controls panel
          Thread th = new Thread() {
            int h = 0, w = 0;
            LayoutManager lm = p.getLayout();
            Dimension d2 = lm.preferredLayoutSize(p);
           
            int tow = (int)d2.getWidth(), toh = (int)d2.getHeight();
            //toh = (int)d2.getHeight();
            //tow = (int)d2.getWidth();
           
            public void run() {
              while(h<toh || w<tow) {
                if((h+10)<toh)
                  h += 10;
                else if(h<toh)
                  h = toh;
                if((w+10)<tow)
                  w += 10;
                else if(w<tow)
                  w = tow;
                p.setPreferredSize(new Dimension(w, h));
                //p.invalidate();
                jTbTools.revalidate();
                //paint(Temp4.this.getGraphics());
                try {this.sleep(30);}
                catch(InterruptedException ie) {ie.printStackTrace(); break;}
              }
              p.setPreferredSize(new Dimension(tow,toh)); jTbTools.revalidate();
            }
          };
          th.start();
           */
        }
        else {
          p.setPreferredSize( new Dimension(0,0) );
          jTbTools.revalidate();
          /*
          Thread th = new Thread() {
            int h = p.getHeight(), w = p.getWidth();
            LayoutManager lm = p.getLayout();
            int tow = 0, toh = 0;
           
            public void run() {
              while(h>toh || w>tow) {
                if((h-10)>toh)
                  h -= 10;
                else if(h>toh)
                  h = toh;
                if((w-10)>tow)
                  w -= 10;
                else if(w>tow)
                  w = tow;
           
                p.setPreferredSize(new Dimension(w, h));
                //p.invalidate();
                jTbTools.revalidate();
                //paint(Temp4.this.getGraphics());
                try {this.sleep(30);}
                catch(InterruptedException ie) {ie.printStackTrace(); break;}
              }
              p.setPreferredSize(new Dimension(tow,toh)); jTbTools.revalidate();
            }
          };
          th.start();
           */
        }
      }
    });
    this.setLayout( new BorderLayout() );
    this.add(jTbTools, BorderLayout.NORTH);
    this.add(m_js, BorderLayout.CENTER);
  }
  
  
  /**
   * This method sets the node size that is appropriate
   * considering the maximum label size that is present.
   * It is used internally when custom node size checkbox
   * is unchecked.
   */
  protected void setAppropriateNodeSize() {
    int strWidth;
    if(maxStringWidth==0)
      for(int i=0; i<m_nodes.size(); i++) {
        strWidth = fm.stringWidth(((GraphNode)m_nodes.elementAt(i)).lbl);
        if(strWidth>maxStringWidth)
          maxStringWidth=strWidth;
      }
    nodeWidth = maxStringWidth+4;
    paddedNodeWidth = nodeWidth+8;
    jTfNodeWidth.setText(""+nodeWidth);
    
    nodeHeight = 2*fm.getHeight();
    jTfNodeHeight.setText(""+nodeHeight);
  }
  
  /**
   * Sets the preferred size for m_gp GraphPanel to the
   * minimum size that is neccessary to display the graph.
   */
  protected void setAppropriateSize() {
    int maxX=0, maxY=0;
    
    m_gp.setScale(scale, scale);
    
    for(int i=0; i<m_nodes.size(); i++) {
      GraphNode n = (GraphNode)m_nodes.elementAt(i);
      if(maxX<n.x)
        maxX=n.x;
      if(maxY<n.y)
        maxY=n.y;
    }
    //System.out.println("Scale: "+scale+" paddedWidth: "+paddedNodeWidth+
    //                   " nodeHeight: "+nodeHeight+"\nmaxX: "+maxX+" maxY: "+
    //                   maxY+" final: "+(int)((maxX+paddedNodeWidth+2)*scale)+
    //                   ","+(int)((maxY+nodeHeight+2)*scale) );
    m_gp.setPreferredSize(new Dimension((int)((maxX+paddedNodeWidth+2)*scale),
    (int)((maxY+nodeHeight+2)*scale)));
    //System.out.println("Size set to "+this.getPreferredSize());
  }
  
  
  /**
   * This method is an implementation for LayoutCompleteEventListener
   * class. It sets the size appropriate for m_gp GraphPanel and
   * and revalidates it's container JScrollPane once a
   * LayoutCompleteEvent is received from the LayoutEngine.
   */
  public void layoutCompleted(LayoutCompleteEvent le) {
    setAppropriateSize();
    //m_gp.clearBuffer();
    m_gp.invalidate();
    m_js.revalidate();
    m_gp.repaint();
    jBtLayout.setEnabled(true);
  }
  
  
  /**
   * This method lays out the graph by calling the
   * LayoutEngine's layoutGraph() method. This method
   * should be called to display the graph nicely, unless
   * the input XMLBIF03 already contains some layout
   * information (ie the x,y positions of nodes.
   */
  public void layoutGraph() {
    if(m_le!=null)
      m_le.layoutGraph();
    
  }
  
  /*********************************************************
   *
   *  BIF reader<br>
   *  Reads a graph description in XMLBIF03 from a string
   *
   *********************************************************
   */
  public void readBIF(String instring) throws BIFFormatException {
    BIFParser bp = new BIFParser(instring, m_nodes, m_edges);
    try {
      graphID = bp.parse();
    } catch(BIFFormatException bf) {
      System.out.println("BIF format error");
      bf.printStackTrace();
    }
    catch(Exception ex) { ex.printStackTrace(); return; }
    
    setAppropriateNodeSize();
    if(m_le!=null) {
      m_le.setNodeSize(paddedNodeWidth, nodeHeight);
    }
  } //end readBIF1
  
  /**
   *
   *  BIF reader<br>
   *  Reads a graph description in XMLBIF03 from an InputStrem
   *
   *
   */
  public void readBIF(InputStream instream) throws BIFFormatException {
    BIFParser bp = new BIFParser(instream, m_nodes, m_edges);
    try {
      graphID = bp.parse();
    } catch(BIFFormatException bf) {
      System.out.println("BIF format error");
      bf.printStackTrace();
    }
    catch(Exception ex) { ex.printStackTrace(); return; }
    
    setAppropriateNodeSize();
    if(m_le!=null) {
      m_le.setNodeSize(paddedNodeWidth, nodeHeight);
    }
    setAppropriateSize();
  } //end readBIF2
  
  
  /*********************************************************
   *
   *  Dot reader<br>
   *  Reads a graph description in DOT format from a string
   *
   *********************************************************
   */
  public void readDOT(Reader input) {
    DotParser dp = new DotParser(input, m_nodes, m_edges);
    graphID = dp.parse();
    
    setAppropriateNodeSize();
    if(m_le!=null) {
      m_le.setNodeSize(paddedNodeWidth, nodeHeight);
      jBtLayout.setEnabled(false);
      layoutGraph();
    }
  }
  
  /**
   * The panel which contains the actual graph.
   */
  private class GraphPanel
    extends PrintablePanel {

    /** for serialization */
    private static final long serialVersionUID = -3562813603236753173L;

    public GraphPanel() {
      super();
      this.addMouseListener( new GraphVisualizerMouseListener() );
      this.addMouseMotionListener( new GraphVisualizerMouseMotionListener() );
      this.setToolTipText("");
    }
    
    public String getToolTipText(MouseEvent me) {
      int x, y, nx, ny;
      Rectangle r;
      GraphNode n;
      Dimension d = m_gp.getPreferredSize();
      //System.out.println("Preferred Size: "+this.getPreferredSize()+
      //                   " Actual Size: "+this.getSize());
      x=y=nx=ny=0;
      
      if(d.width < m_gp.getWidth())
        nx = (int)((nx + m_gp.getWidth()/2 - d.width/2)/scale);
      if(d.height < m_gp.getHeight())
        ny = (int)((ny + m_gp.getHeight()/2 - d.height/2)/scale);
      
      r = new Rectangle(0, 0, 
                       (int)(paddedNodeWidth*scale), (int)(nodeHeight*scale));
      x += me.getX(); y += me.getY();
      
      int i;
      for(i=0; i<m_nodes.size(); i++) {
        n = (GraphNode) m_nodes.elementAt(i);
        if(n.nodeType!=NORMAL)
          return null;
        r.x = (int)((nx+n.x)*scale); r.y = (int)((ny+n.y)*scale);
        if(r.contains(x,y)) {
          if(n.probs==null)
            return n.lbl;
          else
            return n.lbl+" (click to view the probability dist. table)";
        }
      }
      return null;
    }
    
    
    public void paintComponent(Graphics gr) {
      Graphics2D g = (Graphics2D)gr;
      RenderingHints rh = new RenderingHints(RenderingHints.KEY_ANTIALIASING,
      RenderingHints.VALUE_ANTIALIAS_ON);
      rh.put(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_SPEED);
      g.setRenderingHints(rh);
      g.scale(scale, scale);
      Rectangle r = g.getClipBounds();
      g.clearRect(r.x,r.y,r.width,r.height);
      //g.setColor(this.getBackground());
      //g.fillRect(0, 0, width+5, height+5);
      int x=0, y=0;
      Dimension d = this.getPreferredSize();
      //System.out.println("Preferred Size: "+this.getPreferredSize()+
      //                   " Actual Size: "+this.getSize());
      
      //initializing x & y to display the graph in the middle
      //if the display area is larger than the graph
      if(d.width < this.getWidth())
        x = (int)((x + this.getWidth()/2 - d.width/2)/scale);
      if(d.height < this.getHeight())
        y = (int)((y + this.getHeight()/2 - d.height/2)/scale);
      
      for(int index=0; index<m_nodes.size(); index++) {
        GraphNode n = (GraphNode) m_nodes.elementAt(index);
        if( n.nodeType==NORMAL) {
          g.setColor( this.getBackground().darker().darker() );
          g.fillOval(x+n.x+paddedNodeWidth-nodeWidth-
                          (paddedNodeWidth-nodeWidth)/2,
                     y+n.y,
                     nodeWidth, nodeHeight);
          
          g.setColor(Color.white);
          //g.setColor(Color.black);
          //System.out.println("drawing "+
          //                   ((GraphNode)m_nodes.elementAt(index)).ID+
          //                   " at "+" x: "+ (x+n.x+paddedNodeWidth/2-
          //      fm.stringWidth( ((GraphNode)m_nodes.elementAt(index)).ID )/2)+
          //		       " y: "+(y+n.y+nodeHeight/2+fm.getHeight()/2-2) );
          
          
          //Draw the node's label if it can fit inside the node's current
          // width otherwise display its ID or otherwise just display its
          // idx in the FastVector (to distinguish it from others)
          // if any can fit in node's current width
          if(fm.stringWidth(n.lbl)<=nodeWidth)
            g.drawString( n.lbl,
            x+n.x+paddedNodeWidth/2
            -fm.stringWidth( n.lbl )/2,
            y+n.y+nodeHeight/2+fm.getHeight()/2-2 );
          else if(fm.stringWidth(n.ID)<=nodeWidth)
            g.drawString( n.ID,
            x+n.x+paddedNodeWidth/2
            -fm.stringWidth( n.ID )/2,
            y+n.y+nodeHeight/2+fm.getHeight()/2-2 );
          else if(fm.stringWidth( Integer.toString(index) )<=nodeWidth)
            g.drawString( Integer.toString(index),
            x+n.x+paddedNodeWidth/2
            -fm.stringWidth( Integer.toString(index) )/2,
            y+n.y+nodeHeight/2+fm.getHeight()/2-2 );
          
          g.setColor(Color.black);
        }
        else {
          //g.draw( new java.awt.geom.QuadCurve2D.Double(n.x+paddedNodeWidth/2, 
          //                                             n.y,
          //				n.x+paddedNodeWidth-nodeSize
          //                                   -(paddedNodeWidth-nodeSize)/2,
          //                                  n.y+nodeHeight/2,
          //                           n.x+paddedNodeWidth/2, n.y+nodeHeight) );
          g.drawLine(x+n.x+paddedNodeWidth/2, y+n.y, 
                     x+n.x+paddedNodeWidth/2, y+n.y+nodeHeight);
          
        }
        
        GraphNode n2;
        int x1, y1, x2, y2;
        //System.out.println("Drawing edges of "+n.lbl);
        
        //Drawing all the edges coming out from the node,
        //including reversed and double ones
        if(n.edges!=null)
          for(int k=0; k<n.edges.length; k++) {
            if(n.edges[k][1]>0) {
              n2 = (GraphNode) m_nodes.elementAt(n.edges[k][0]); //m_nodes.elementAt(k);
              //System.out.println("  -->to "+n2.lbl);
              x1=n.x+paddedNodeWidth/2; y1=n.y+nodeHeight;
              x2=n2.x+paddedNodeWidth/2; y2=n2.y;
              g.drawLine(x+x1, y+y1, x+x2, y+y2);
              if(n.edges[k][1]==DIRECTED) {
                if(n2.nodeType==n2.NORMAL)
                  drawArrow(g, x+x1, y+y1, x+x2, y+y2);
              }
              else if(n.edges[k][1]==REVERSED) {
                if(n.nodeType==NORMAL)
                  drawArrow(g, x+x2, y+y2, x+x1, y+y1);
              }
              else if(n.edges[k][1]==DOUBLE) {
                if(n.nodeType==NORMAL)
                  drawArrow(g, x+x2, y+y2, x+x1, y+y1);
                if(n2.nodeType==NORMAL)
                  drawArrow(g, x+x1, y+y1, x+x2, y+y2);
              }
            }
          }
      }
    }
    
    /**
     * This method draws an arrow on a line from (x1,y1)
     * to (x2,y2). The arrow head is seated on (x2,y2) and
     * is in the direction of the line.
     * If the arrow is needed to be drawn in the opposite
     * direction then simply swap the order of (x1, y1)
     * and (x2, y2) when calling this function.
     */
    protected void drawArrow(Graphics g, int x1, int y1, int x2, int y2) {
      
      if(x1==x2) {
        if(y1<y2) {
          g.drawLine(x2, y2, x2+4, y2-8);
          g.drawLine(x2, y2, x2-4, y2-8);
        }
        else {
          g.drawLine(x2, y2, x2+4, y2+8);
          g.drawLine(x2, y2, x2-4, y2+8);
        }
      }
      else {
        //theta=line's angle from base, beta=angle of arrow's side from line
        double hyp=0, base=0, perp=0, theta, beta;
        int x3=0, y3=0;
        
        if(x2<x1) {
          base = x1-x2; hyp = Math.sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) );
          theta = Math.acos( base/hyp );
        }
        else { //x1>x2 as we already checked x1==x2 before
          base = x1-x2; hyp = Math.sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) );
          theta = Math.acos( base/hyp );
        }
        beta = 30*Math.PI/180;
        //System.out.println("Original base "+base+" perp "+perp+" hyp "+hyp+
        //                   "\ntheta "+theta+" beta "+beta);
        
        hyp = 8;
        base = Math.cos(theta-beta)*hyp;
        perp = Math.sin(theta-beta)*hyp;
        
        x3 = (int)(x2+base);
        if(y1<y2)
          y3 = (int)(y2-perp);
        else
          y3 = (int)(y2+perp);
        
        //System.out.println("Drawing 1 from "+x2+","+y2+" to "+x3+","+y3+
        //                   " x1,y1 is "+x1+","+y1+" base "+base+
        //		     " perp "+perp+" cos(theta-beta) "+
        //                   Math.cos(theta-beta));
        g.drawLine(x2, y2, x3, y3);
        
        base = Math.cos(theta+beta)*hyp;
        perp = Math.sin(theta+beta)*hyp;
        
        x3 = (int)(x2+base);
        if(y1<y2)
          y3 = (int)(y2-perp);
        else
          y3 = (int)(y2+perp);
        //System.out.println("Drawing 2 from "+x2+","+y2+" to "+x3+","+y3+
        //                   " x1,y1 is "+x1+","+y1+" base "+base+
        //		     " perp "+perp);
        g.drawLine(x2, y2, x3, y3);
      }
    }
    
    /**
     * This method highlights a given node and all its children
     * and the edges coming out of it.
     */
    public void highLight(GraphNode n) {
      Graphics2D g = (Graphics2D) this.getGraphics();
      RenderingHints rh = new RenderingHints(RenderingHints.KEY_ANTIALIASING,
      RenderingHints.VALUE_ANTIALIAS_ON);
      rh.put(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_SPEED);
      g.setRenderingHints(rh);
      g.setPaintMode();
      g.scale(scale, scale);
      int x=0, y=0;
      Dimension d = this.getPreferredSize();
      //System.out.println("Preferred Size: "+this.getPreferredSize()+
      //                   " Actual Size: "+this.getSize());
      
      //initializing x & y to display the graph in the middle
      //if the display area is larger than the graph
      if(d.width < this.getWidth())
        x = (int)((x + this.getWidth()/2 - d.width/2)/scale);
      if(d.height < this.getHeight())
        y = (int)((y + this.getHeight()/2 - d.height/2)/scale);
      
      //if the node is of type NORMAL only then highlight
      if(n.nodeType==NORMAL) {
        
        g.setXORMode(Color.green); //g.setColor(Color.green);
        
        g.fillOval(x+n.x+paddedNodeWidth-nodeWidth-
        (paddedNodeWidth-nodeWidth)/2,
        y+n.y, nodeWidth, nodeHeight);
        g.setXORMode(Color.red);
        
        //Draw the node's label if it can fit inside the node's current
        // width otherwise display its ID or otherwise just display its
        // idx in the FastVector (to distinguish it from others)
        // if any can fit in node's current width
        if(fm.stringWidth(n.lbl)<=nodeWidth)
          g.drawString( n.lbl,
          x+n.x+paddedNodeWidth/2
          -fm.stringWidth( n.lbl )/2,
          y+n.y+nodeHeight/2+fm.getHeight()/2-2 );
        else if(fm.stringWidth(n.ID)<=nodeWidth)
          g.drawString( n.ID,
          x+n.x+paddedNodeWidth/2
          -fm.stringWidth( n.ID )/2,
          y+n.y+nodeHeight/2+fm.getHeight()/2-2 );
        else if( fm.stringWidth( Integer.toString(m_nodes.indexOf(n)) ) <=
        nodeWidth )
          g.drawString( Integer.toString(m_nodes.indexOf(n)),
          x+n.x+paddedNodeWidth/2
          -fm.stringWidth( Integer.toString(m_nodes.indexOf(n)) )/2,
          y+n.y+nodeHeight/2+fm.getHeight()/2-2 );
        
        g.setXORMode(Color.green);
        
        
        GraphNode n2;
        int x1, y1, x2, y2;
        //System.out.println("Drawing edges of "+n.lbl);
        if(n.edges!=null)
          //Drawing all the edges from and upward ones coming to the node
          for(int k=0; k<n.edges.length; k++) {
            if(n.edges[k][1]==DIRECTED || n.edges[k][1]==DOUBLE) {
              n2 = (GraphNode) m_nodes.elementAt(n.edges[k][0]); //m_nodes.elementAt(k);
              //System.out.println("  -->to "+n2.lbl);
              x1=n.x+paddedNodeWidth/2; y1=n.y+nodeHeight;
              x2=n2.x+paddedNodeWidth/2; y2=n2.y;
              g.drawLine(x+x1, y+y1, x+x2, y+y2);
              if(n.edges[k][1]==DIRECTED) {
                if(n2.nodeType==n2.NORMAL) //!n2.dummy)
                  drawArrow(g, x+x1, y+y1, x+x2, y+y2);
              }
              else if(n.edges[k][1]==DOUBLE) {
                if(n.nodeType==NORMAL) //!n.dummy)
                  drawArrow(g, x+x2, y+y2, x+x1, y+y1);
                if(n2.nodeType==NORMAL) //!n2.dummy)
                  drawArrow(g, x+x1, y+y1, x+x2, y+y2);
              }
              if(n2.nodeType==NORMAL)
                g.fillOval(x+n2.x+paddedNodeWidth-nodeWidth-
                (paddedNodeWidth-nodeWidth)/2,
                y+n2.y, nodeWidth, nodeHeight);
              
              //If n2 is not of NORMAL type
              // then carry on drawing all the edges and add all the
              // dummy nodes encountered in a Vector until no
              // more dummy nodes are found and all the child nodes(node n2)
              // are of type normal
              java.util.Vector t = new java.util.Vector();
              while(n2.nodeType!=NORMAL || t.size()>0) { //n2.dummy==true) {
                //System.out.println("in while processing "+n2.ID);
                if(t.size()>0)
                { n2 = (GraphNode)t.elementAt(0);
                  t.removeElementAt(0); }
                if(n2.nodeType!=NORMAL) {
                  g.drawLine(x+n2.x+paddedNodeWidth/2, y+n2.y,
                  x+n2.x+paddedNodeWidth/2, y+n2.y+nodeHeight);
                  x1=n2.x+paddedNodeWidth/2; y1=n2.y+nodeHeight;
                  //System.out.println("Drawing from "+n2.lbl);
                  for(int m=0; m<n2.edges.length; m++) {
                    //System.out.println(" to "+n2.lbl+", "+
                    //                   graphMatrix[tmpIndex][m]);
                    if(n2.edges[m][1]>0) {
                      GraphNode n3 =
                      (GraphNode) m_nodes.elementAt(n2.edges[m][0]); //m_nodes.elementAt(m);
                      g.drawLine(x+x1, y+y1, x+n3.x+paddedNodeWidth/2, y+n3.y);
                      
                      if(n3.nodeType==NORMAL){ //!n2.dummy)
                        g.fillOval(x+n3.x+paddedNodeWidth-nodeWidth-
                        (paddedNodeWidth-nodeWidth)/2,
                        y+n3.y, nodeWidth, nodeHeight);
                        drawArrow(g, x+x1, y+y1,
                        x+n3.x+paddedNodeWidth/2, y+n3.y);
                      }
                      //if(n3.nodeType!=n3.NORMAL)
                      t.addElement(n3);
                      //break;
                    }
                  }
                }
              }
            }
            else if(n.edges[k][1]==-REVERSED || n.edges[k][1]==-DOUBLE) {
              //Drawing all the reversed and double edges which are going
              //upwards in the drawing.
              n2 = (GraphNode) m_nodes.elementAt(n.edges[k][0]); //m_nodes.elementAt(k);
              //System.out.println("  -->to "+n2.lbl);
              x1=n.x+paddedNodeWidth/2; y1=n.y;
              x2=n2.x+paddedNodeWidth/2; y2=n2.y+nodeHeight;
              g.drawLine(x+x1, y+y1, x+x2, y+y2);
              
              if(n.edges[k][1]==-DOUBLE) {
                drawArrow(g, x+x2, y+y2, x+x1, y+y1);
                if(n2.nodeType!=SINGULAR_DUMMY) //!n2.dummy)
                  drawArrow(g, x+x1, y+y1, x+x2, y+y2);
              }
              
              int tmpIndex=k;
              while(n2.nodeType!=NORMAL) { //n2.dummy==true) {
                g.drawLine(x+n2.x+paddedNodeWidth/2,
                y+n2.y+nodeHeight, x+n2.x+paddedNodeWidth/2, y+n2.y);
                x1=n2.x+paddedNodeWidth/2; y1=n2.y;
                for(int m=0; m<n2.edges.length; m++) {
                  if(n2.edges[m][1]<0) {
                    n2 = (GraphNode) m_nodes.elementAt(n2.edges[m][0]); //m_nodes.elementAt(m);
                    g.drawLine(x+x1, y+y1,
                    x+n2.x+paddedNodeWidth/2, y+n2.y+nodeHeight);
                    tmpIndex=m;
                    if(n2.nodeType!=SINGULAR_DUMMY) //!n2.dummy)
                      drawArrow(g, x+x1, y+y1,
                      x+n2.x+paddedNodeWidth/2, y+n2.y+nodeHeight);
                    break;
                  }
                }
              }
            }
          }
      }
    }
  }
  
  
  /**
   * Table Model for the Table that shows the probability
   * distribution for a node
   */
  private class GraphVisualizerTableModel
    extends AbstractTableModel {

    /** for serialization */
    private static final long serialVersionUID = -4789813491347366596L;
    
    final String[] columnNames;
    final double[][] data;
    
    
    public GraphVisualizerTableModel(double[][] d, String[] c) {
      data = d;
      columnNames = c;
    }
    
    public int getColumnCount() {
      return columnNames.length;
    }
    
    public int getRowCount() {
      return data.length;
    }
    
    public String getColumnName(int col) {
      return columnNames[col];
    }
    
    public Object getValueAt(int row, int col) {
      return new Double(data[row][col]);
    }
    
   /*
    * JTable uses this method to determine the default renderer/
    * editor for each cell.
    */
    public Class getColumnClass(int c) {
      return getValueAt(0, c).getClass();
    }
    
    /*
     * Implemented this to make sure the table is uneditable.
     */
    public boolean isCellEditable(int row, int col) {
      return false;
    }
  }
  
  
  
  /**
   * Listener class for processing mouseClicked
   */
  private class GraphVisualizerMouseListener extends MouseAdapter {
    int x, y, nx, ny; Rectangle r;
    
    /**
     * If the mouse is clicked on a node then this method
     * displays a dialog box with the probability distribution
     * table for that node IF it exists
     */
    public void mouseClicked(MouseEvent me) {
      GraphNode n;
      Dimension d = m_gp.getPreferredSize();
      //System.out.println("Preferred Size: "+this.getPreferredSize()+
      //                   " Actual Size: "+this.getSize());
      x=y=nx=ny=0;
      
      if(d.width < m_gp.getWidth())
        nx = (int)((nx + m_gp.getWidth()/2 - d.width/2)/scale);
      if(d.height < m_gp.getHeight())
        ny = (int)((ny + m_gp.getHeight()/2 - d.height/2)/scale);
      
      r=new Rectangle(0, 0, 
                     (int)(paddedNodeWidth*scale), (int)(nodeHeight*scale));
      x += me.getX(); y += me.getY();
      
      int i;
      for(i=0; i<m_nodes.size(); i++) {
        n = (GraphNode) m_nodes.elementAt(i);
        r.x = (int)((nx+n.x)*scale); r.y = (int)((ny+n.y)*scale);
        if(r.contains(x,y)) {
          if(n.probs==null)
            return;
          
          int noOfPrntsOutcomes = 1;
          if(n.prnts!=null) {
            for(int j=0; j<n.prnts.length; j++) {
              GraphNode n2 = (GraphNode)m_nodes.elementAt(n.prnts[j]);
              noOfPrntsOutcomes *= n2.outcomes.length;
            }
            if(noOfPrntsOutcomes>511) {
              System.err.println("Too many outcomes of parents ("+noOfPrntsOutcomes+
                                 ") can't display probabilities");
              return;
            }
          }
          
          GraphVisualizerTableModel tm = 
                             new GraphVisualizerTableModel(n.probs, n.outcomes);
          
          JTable jTblProbs = new JTable(tm); //JTable(probabilities, (Object[])n.outcomes);
          
          JScrollPane js = new JScrollPane(jTblProbs);
          
          if(n.prnts!=null) {
            GridBagConstraints gbc = new GridBagConstraints();
            JPanel jPlRowHeader = new JPanel( new GridBagLayout() );
            
            //indices of the parent nodes in the Vector
            int [] idx = new int[n.prnts.length]; 
            //max length of values of each parent
            int [] lengths = new int[n.prnts.length];  
            
            //System.out.println("n.probs.length "+n.probs.length+
            //                   " should be "+noOfPrntsOutcomes);
            //System.out.println("n.probs[0].length "+n.probs[0].length+
            //                   " should be "+n.outcomes.length);
            //System.out.println("probabilities are: ");
            //for(int j=0; j<probabilities.length; j++) {
            //    for(int k=0; k<probabilities[j].length; k++)
            //	   System.out.print(probabilities[j][k]+" ");
            //     System.out.println("");
            //}
            
            //Adding labels for rows
            gbc.anchor = gbc.NORTHWEST;
            gbc.fill = gbc.HORIZONTAL;
            gbc.insets = new Insets(0,1,0,0);
            int addNum=0, temp=0;
            boolean dark=false;
            while(true){
              GraphNode n2;
              gbc.gridwidth = 1;
              for(int k=0; k<n.prnts.length; k++) {
                n2 = (GraphNode)m_nodes.elementAt(n.prnts[k]);
                JLabel lb = new JLabel(n2.outcomes[idx[k]]);
                lb.setFont( new Font("Dialog", Font.PLAIN, 12) );
                lb.setOpaque( true );
                lb.setBorder( BorderFactory.createEmptyBorder( 1,2,1,1 ) );
                lb.setHorizontalAlignment( JLabel.CENTER );
                if(dark) {
                  lb.setBackground( lb.getBackground().darker() );
                  lb.setForeground( Color.white );
                }
                else
                  lb.setForeground( Color.black );
                
                temp = lb.getPreferredSize().width;
                //System.out.println("Preferred width "+temp+
                //                   " for "+n2.outcomes[idx[k]]);
                lb.setPreferredSize(
                                 new Dimension(temp, jTblProbs.getRowHeight())
                                 );
                if(lengths[k]<temp)
                  lengths[k] = temp;
                temp=0;
                
                if(k==n.prnts.length-1) {
                  gbc.gridwidth = gbc.REMAINDER;
                  dark = (dark==true) ?  false:true;
                }
                jPlRowHeader.add(lb, gbc);
                addNum++;
              }
              
              for(int k=n.prnts.length-1; k>=0; k--) {
                n2 = (GraphNode) m_nodes.elementAt(n.prnts[k]);
                if(idx[k]==n2.outcomes.length-1 && k!=0) {
                  idx[k]=0;
                  continue;
                }
                else {
                  idx[k]++;
                  break;
                }
              }
              
              n2 = (GraphNode) m_nodes.elementAt(n.prnts[0]);
              if(idx[0]==n2.outcomes.length) {
                JLabel lb= (JLabel) jPlRowHeader.getComponent(addNum-1);
                jPlRowHeader.remove(addNum-1);
                lb.setPreferredSize( new Dimension(lb.getPreferredSize().width, 
                                                   jTblProbs.getRowHeight()) );
                gbc.gridwidth = gbc.REMAINDER;
                gbc.weighty = 1;
                jPlRowHeader.add(lb, gbc);
                gbc.weighty=0;
                break;
              }
            }
            
            
            gbc.gridwidth = 1;
            //The following panel contains the names of the parents
            //and is displayed above the row names to identify
            //which value belongs to which parent
            JPanel jPlRowNames = new JPanel(new GridBagLayout());
            for(int j=0; j<n.prnts.length; j++) {
              JLabel lb2;
              JLabel lb1 = 
                   new JLabel( ((GraphNode)m_nodes.elementAt(n.prnts[j])).lbl );
              lb1.setBorder( BorderFactory.createEmptyBorder( 1,2,1,1 ) );
              Dimension tempd = lb1.getPreferredSize();
              //System.out.println("lengths[j]: "+lengths[j]+
              //                   " tempd.width: "+tempd.width);
              if(tempd.width<lengths[j]) {
                lb1.setPreferredSize( new Dimension(lengths[j], tempd.height) );
                lb1.setHorizontalAlignment( JLabel.CENTER );
                lb1.setMinimumSize( new Dimension(lengths[j], tempd.height) );
              }
              else if(tempd.width>lengths[j]) {
                lb2 = (JLabel) jPlRowHeader.getComponent(j);
                lb2.setPreferredSize( new Dimension(tempd.width, 
                                               lb2.getPreferredSize().height) );
              }
              jPlRowNames.add(lb1, gbc);
              //System.out.println("After adding "+lb1.getPreferredSize());
            }
            js.setRowHeaderView(jPlRowHeader);
            js.setCorner( JScrollPane.UPPER_LEFT_CORNER, jPlRowNames );
          }
          
          
          JDialog jd = 
                new JDialog((Frame)GraphVisualizer.this.getTopLevelAncestor(),
                            "Probability Distribution Table For "+n.lbl, true);
          jd.setSize(500, 400);
          jd.setLocation(GraphVisualizer.this.getLocation().x+
                            GraphVisualizer.this.getWidth()/2-250,
                         GraphVisualizer.this.getLocation().y+
                            GraphVisualizer.this.getHeight()/2-200 );
          
          jd.getContentPane().setLayout( new BorderLayout() );
          jd.getContentPane().add(js, BorderLayout.CENTER);
          jd.setVisible(true);
          
          return;
        }
      }
    }
    
  }
  
  
  /**
   * private class for handling mouseMoved events
   * to highlight nodes if the the mouse is moved on
   * one
   */
  private class GraphVisualizerMouseMotionListener extends MouseMotionAdapter {
    int x, y, nx, ny; Rectangle r;
    GraphNode lastNode;
    
    public void mouseMoved(MouseEvent me) {
      GraphNode n;
      Dimension d = m_gp.getPreferredSize();
      //System.out.println("Preferred Size: "+this.getPreferredSize()+
      //                   " Actual Size: "+this.getSize());
      x=y=nx=ny=0;
      
      if(d.width < m_gp.getWidth())
        nx = (int)((nx + m_gp.getWidth()/2 - d.width/2)/scale);
      if(d.height < m_gp.getHeight())
        ny = (int)((ny + m_gp.getHeight()/2 - d.height/2)/scale);
      
      r=new Rectangle(0, 0, 
                     (int)(paddedNodeWidth*scale), (int)(nodeHeight*scale));
      x += me.getX(); y += me.getY();
      
      int i;
      for(i=0; i<m_nodes.size(); i++) {
        n = (GraphNode) m_nodes.elementAt(i);
        r.x = (int)((nx+n.x)*scale); r.y = (int)((ny+n.y)*scale);
        if(r.contains(x,y)) {
          if(n!=lastNode) {
            m_gp.highLight(n);
            if(lastNode!=null)
              m_gp.highLight(lastNode);
            lastNode = n; //lastIndex = i;
          }
          break;
        }
      }
      if(i==m_nodes.size()  && lastNode!=null) {
        m_gp.repaint();
        //m_gp.highLight(lastNode);
        lastNode=null;
      }
    }
  }
  
  /**
   * Main method to load a text file with the
   * description of a graph from the command
   * line
   */
  public static void main(String [] args) {
    weka.core.logging.Logger.log(weka.core.logging.Logger.Level.INFO, "Logging started");
    JFrame jf = new JFrame("Graph Visualizer");
    GraphVisualizer g = new GraphVisualizer();
    
    try{
      if(args[0].endsWith(".xml")) {
        //StringBuffer sb = new StringBuffer();
        //FileReader infile = new FileReader(args[0]);
        //int i;
        //while( (i=infile.read())!=-1) {
        //    sb.append((char)i);
        //}
        //System.out.println(sb.toString());
        //g.readBIF(sb.toString() );
        g.readBIF( new FileInputStream(args[0]) );
      }
      else {
        //BufferedReader infile=new BufferedReader();
        g.readDOT(new FileReader(args[0])); //infile);
      }
    }
    catch(IOException ex) { ex.printStackTrace(); }
    catch(BIFFormatException bf) { bf.printStackTrace(); System.exit(-1); }
    
    jf.getContentPane().add(g);
    //RepaintManager.currentManager(jf.getRootPane()).setDoubleBufferingEnabled(false);
    jf.setDefaultCloseOperation( jf.EXIT_ON_CLOSE );
    jf.setSize(800,600);
    //jf.pack();
    jf.setVisible(true);
  }
}
