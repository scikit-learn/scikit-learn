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
 *    AttributePanel.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.visualize;

import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instances;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;

import javax.swing.JPanel;
import javax.swing.JScrollPane;

/**
 * This panel displays one dimensional views of the attributes in a
 * dataset. Colouring is done on the basis of a column in the dataset or
 * an auxiliary array (useful for colouring cluster predictions).
 * 
 * @author Malcolm Ware (mfw4@cs.waikato.ac.nz)
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @version $Revision: 5086 $
 */
public class AttributePanel
  extends JScrollPane {

  /** for serialization */
  private static final long serialVersionUID = 3533330317806757814L;
  
  /** The instances to be plotted */
  protected Instances m_plotInstances=null;
    
  /** Holds the min and max values of the colouring attributes */
  protected double m_maxC;
  protected double m_minC;
  protected int m_cIndex;
  protected int m_xIndex;
  protected int m_yIndex;

  /** The colour map to use for colouring points */
  protected FastVector m_colorList;

   /** default colours for colouring discrete class */
    protected Color [] m_DefaultColors = {Color.blue,
					Color.red,
					Color.green,
					Color.cyan,
					Color.pink,
					new Color(255, 0, 255),
					Color.orange,
					new Color(255, 0, 0),
					new Color(0, 255, 0),
					Color.white};
    
  /**
   *  If set, it allows this panel to avoid setting a color in
   * the color list that is equal to the background color
   */
  protected Color m_backgroundColor = null; 

  /** The list of things listening to this panel */
  protected FastVector m_Listeners = new FastVector();

  /** Holds the random height for each instance. */
  protected int[] m_heights;
  //protected Color[] colors_array;

  /** The container window for the attribute bars, and also where the
   * X,Y or B get printed.
   */ 
  protected JPanel m_span=null;

  /** The default colour to use for the background of the bars if
      a colour is not defined in Visualize.props */
  protected Color m_barColour=Color.black;
    
  /** inner inner class used for plotting the points 
   * into a bar for a particular attribute. 
   */
  protected class AttributeSpacing
    extends JPanel {

    /** for serialization */
    private static final long serialVersionUID = 7220615894321679898L;

    /** The min and max values for this attribute. */
    protected double m_maxVal;
    protected double m_minVal;

    /** The attribute itself. */
    protected Attribute m_attrib;
      
    /** The index for this attribute. */
    protected int m_attribIndex;
      
    /** The x position of each point. */
    protected int[] m_cached;
    //note for m_cached, if you wanted to speed up the drawing algorithm
    // and save memory, the system could be setup to drop out any
    // of the instances not being drawn (you would need to find a new way
    //of matching the height however).

    /** A temporary array used to strike any instances that would be 
     * drawn redundantly.
     */
    protected boolean[][] m_pointDrawn;
      
    /** Used to determine if the positions need to be recalculated. */
    protected int m_oldWidth=-9000;

    /** The container window for the attribute bars, and also where the
     * X,Y or B get printed.
     */
      
    /**
     * This constructs the bar with the specified attribute and
     * sets its index to be used for selecting by the mouse.
     * @param a The attribute this bar represents.
     * @param aind The index of this bar.
     */
    public AttributeSpacing(Attribute a, int aind) {
      m_attrib = a;
      m_attribIndex = aind;
      this.setBackground(m_barColour);
      this.setPreferredSize(new Dimension(0, 20));
      this.setMinimumSize(new Dimension(0, 20));
      m_cached = new int[m_plotInstances.numInstances()];
	
      //this will only get allocated if m_plotInstances != null
      //this is used to determine the min and max values for plotting
      double min=Double.POSITIVE_INFINITY;
      double max=Double.NEGATIVE_INFINITY;
      double value;
      if (m_plotInstances.attribute(m_attribIndex).isNominal() || m_plotInstances.attribute(m_attribIndex).isRanking()) {
	m_minVal = 0;
	m_maxVal = m_plotInstances.attribute(m_attribIndex).numValues()-1;
      } else {
	for (int i=0;i<m_plotInstances.numInstances();i++) {
	  if (!m_plotInstances.instance(i).isMissing(m_attribIndex)) {
	    value = m_plotInstances.instance(i).value(m_attribIndex);
	    if (value < min) {
	      min = value;
	    }
	    if (value > max) {
	      max = value;
	    }
	  }
	}
	m_minVal = min; m_maxVal = max;
	if (min == max) {
	  m_maxVal += 0.05;
	  m_minVal -= 0.05;
	}
      }
	
      this.addMouseListener(new MouseAdapter() {
	  public void mouseClicked(MouseEvent e) {
	    if ((e.getModifiers() & e.BUTTON1_MASK) == e.BUTTON1_MASK) {
	      setX(m_attribIndex);
	      if (m_Listeners.size() > 0) {
		for (int i=0;i<m_Listeners.size();i++) {
		  AttributePanelListener l = 
		    (AttributePanelListener)(m_Listeners.elementAt(i));
		  l.attributeSelectionChange(new AttributePanelEvent(true,
					     false, m_attribIndex));
		}
	      }
	    }
	    else {
	      //put it on the y axis
	      setY(m_attribIndex);
	      if (m_Listeners.size() > 0) {
		for (int i=0;i<m_Listeners.size();i++) {
		  AttributePanelListener l = 
		    (AttributePanelListener)(m_Listeners.elementAt(i));
		  l.attributeSelectionChange(new AttributePanelEvent(false,
					     true, m_attribIndex));
		}
	      }
	    }
	  }
	});
    }
      
    /**
     * Convert an raw x value to Panel x coordinate.
     * @param val the raw x value
     * @return an x value for plotting in the panel.
     */
    private double convertToPanel(double val) {
      double temp = (val - m_minVal)/(m_maxVal - m_minVal);
      double temp2 = temp * (this.getWidth() - 10);
	
      return temp2 + 4; 
    }
      
    /**
     * paints all the visible instances to the panel , and recalculates
     * their position if need be.
     * @param gx The graphics context.
     */
    public void paintComponent(Graphics gx) {
      super.paintComponent(gx);
      int xp, yp, h;
      h = this.getWidth();
      if (m_plotInstances != null 
	  && m_plotInstances.numAttributes() > 0
	  && m_plotInstances.numInstances() > 0) {

	if (m_oldWidth != h) {
	  m_pointDrawn = new boolean[h][20];
	  for (int noa = 0; noa < m_plotInstances.numInstances(); noa++) {
	    if (!m_plotInstances.instance(noa).isMissing(m_attribIndex)
		&& !m_plotInstances.instance(noa).isMissing(m_cIndex)) {
	      m_cached[noa] = (int)convertToPanel(m_plotInstances.
						  instance(noa).
						  value(m_attribIndex));
		
	      if (m_pointDrawn[m_cached[noa] % h][m_heights[noa]]) {
		m_cached[noa] = -9000;
	      }
	      else {
		m_pointDrawn[m_cached[noa]%h][m_heights[noa]] = true;
	      }
		
	    }
	    else {
	      m_cached[noa] = -9000; //this value will not happen 
	      //so it is safe
	    }
	  }
	  m_oldWidth = h;
	}
	  
	if (m_plotInstances.attribute(m_cIndex).isNominal() || m_plotInstances.attribute(m_cIndex).isRanking()) {
	  for (int noa = 0; noa < m_plotInstances.numInstances(); noa++) {
	      
	    if (m_cached[noa] != -9000) {
	      xp = m_cached[noa];
	      yp = m_heights[noa];
	      if (m_plotInstances.attribute(m_attribIndex).
		  isNominal() || m_plotInstances.attribute(m_attribIndex).
		  isRanking()) {
		xp += (int)(Math.random() * 5) - 2;
	      }
	      int ci = (int)m_plotInstances.instance(noa).value(m_cIndex);

	      gx.setColor((Color)m_colorList.elementAt
			  (ci % m_colorList.size()));
	      gx.drawRect(xp, yp, 1, 1);
	    }
	  }
	}
	else {
	  double r;
	  for (int noa = 0; noa < m_plotInstances.numInstances(); noa++) {
	    if (m_cached[noa] != -9000) {		  
		
	      r = (m_plotInstances.instance(noa).value(m_cIndex) 
		   - m_minC) / (m_maxC - m_minC);

	      r = (r * 240) + 15;

	      gx.setColor(new Color((int)r,150,(int)(255-r)));
		
	      xp = m_cached[noa];
	      yp = m_heights[noa];
	      if (m_plotInstances.attribute(m_attribIndex).
		  isNominal() || m_plotInstances.attribute(m_attribIndex).
		  isRanking()) {
		xp += (int)(Math.random() * 5) - 2;
	      }
	      gx.drawRect(xp, yp, 1, 1);
	    }
	  }
	}
      } 
    }
  }   
    
  /**
   * Set the properties for the AttributePanel
   */
  private void setProperties() {
    if (VisualizeUtils.VISUALIZE_PROPERTIES != null) {
      String thisClass = this.getClass().getName();
      String barKey = thisClass+".barColour";
      
      String barC = VisualizeUtils.VISUALIZE_PROPERTIES.
	      getProperty(barKey);
      if (barC == null) {
	/*
	System.err.println("Warning: no configuration property found in "
			   +VisualizeUtils.PROPERTY_FILE
			   +" for "+barKey);
	*/
      } else {
	//System.err.println("Setting attribute bar colour to: "+barC);
	m_barColour = VisualizeUtils.processColour(barC, m_barColour);
      }
    }
  }
  
  public AttributePanel() {
    this(null);
  }
 
  /**
   * This constructs an attributePanel.
   */
  public AttributePanel(Color background) {
    m_backgroundColor = background;
    
    setProperties();
    this.setBackground(Color.blue);
    setVerticalScrollBarPolicy(VERTICAL_SCROLLBAR_ALWAYS);
    m_colorList = new FastVector(10);

    for (int noa = m_colorList.size(); noa < 10; noa++) {
      Color pc = m_DefaultColors[noa % 10];
      int ija =  noa / 10;
      ija *= 2; 
      for (int j=0;j<ija;j++) {
	pc = pc.darker();
      }
      
      m_colorList.addElement(pc);
    }
  }

  /**
   * Add a listener to the list of things listening to this panel
   * @param a the listener to notify when attribute bars are clicked on
   */
  public void addAttributePanelListener(AttributePanelListener a) {
    m_Listeners.addElement(a);
  }
  
  /**
   * Set the index of the attribute by which to colour the data. Updates
   * the number of entries in the colour list if there are more values
   * for this new attribute than previous ones.
   * @param c the index of the attribute to colour on
   * @param h maximum value of this attribute
   * @param l minimum value of this attribute
   */
  public void setCindex(int c, double h, double l) {
    m_cIndex = c;
    m_maxC = h;
    m_minC = l;
    
    if (m_span != null) {
      if (m_plotInstances.numAttributes() > 0 &&
	  m_cIndex < m_plotInstances.numAttributes()) {
	if (m_plotInstances.attribute(m_cIndex).isNominal() || m_plotInstances.attribute(m_cIndex).isRanking()) {
	  if (m_plotInstances.attribute(m_cIndex).numValues() > 
	    m_colorList.size()) {
	    extendColourMap();
	  }
	}
      }
      this.repaint();
    }
  }

  /**
   * Set the index of the attribute by which to colour the data. Updates
   * the number of entries in the colour list if there are more values
   * for this new attribute than previous ones.
   * @param c the index of the attribute to colour on
   */
  public void setCindex(int c) {
    m_cIndex = c;
    /*    m_maxC = h;
	  m_minC = l; */

    if (m_span != null) {
      if (m_cIndex < m_plotInstances.numAttributes() && 
	  m_plotInstances.attribute(m_cIndex).isNumeric()) {
	double min=Double.POSITIVE_INFINITY;
	double max=Double.NEGATIVE_INFINITY;
	double value;

	for (int i=0;i<m_plotInstances.numInstances();i++) {
	  if (!m_plotInstances.instance(i).isMissing(m_cIndex)) {
	    value = m_plotInstances.instance(i).value(m_cIndex);
	    if (value < min) {
	      min = value;
	    }
	    if (value > max) {
	      max = value;
	    }
	  }
	}
    
	m_minC = min; m_maxC = max;
      } else {
	if (m_plotInstances.attribute(m_cIndex).numValues() > 
	    m_colorList.size()) {
	  extendColourMap();
	}
      }
    
      this.repaint();
    }
  }

  /**
   * Adds more colours to the colour list
   */
  private void extendColourMap() {
    if (m_plotInstances.attribute(m_cIndex).isNominal() || m_plotInstances.attribute(m_cIndex).isRanking()) {
      for (int i = m_colorList.size(); 
	   i < m_plotInstances.attribute(m_cIndex).numValues();
	   i++) {
	Color pc = m_DefaultColors[i % 10];
	int ija =  i / 10;
	ija *= 2; 
	for (int j=0;j<ija;j++) {
	  pc = pc.brighter();
	}
	
	if (m_backgroundColor != null) {
	  pc = Plot2D.checkAgainstBackground(pc, m_backgroundColor);
	}

	m_colorList.addElement(pc);
      }
    }
  }

  /**
   * Sets a list of colours to use for colouring data points
   * @param cols a list of java.awt.Color
   */
  public void setColours(FastVector cols) {
    m_colorList = cols;
  }
  
  protected void setDefaultColourList(Color[] list) {
    m_DefaultColors = list;
  }

  /** 
   * This sets the instances to be drawn into the attribute panel
   * @param ins The instances.
   */
  public void setInstances(Instances ins) throws Exception {
    if (ins.numAttributes() > 512) {
      throw new Exception("Can't display more than 512 attributes!");
    }

    if (m_span == null) {
      m_span = new JPanel() {
	  private static final long serialVersionUID = 7107576557995451922L;
	  
	  public void paintComponent(Graphics gx) {
	    super.paintComponent(gx);
	    gx.setColor(Color.red);
	    if (m_yIndex != m_xIndex) {
	      gx.drawString("X", 5, m_xIndex * 20 + 16);
	      gx.drawString("Y", 5, m_yIndex * 20 + 16);
	    }
	    else {
	      gx.drawString("B", 5, m_xIndex * 20 + 16);
	    }
	  }
	};
    }

    m_span.removeAll();
    m_plotInstances = ins;
    if (ins.numInstances() > 0 && ins.numAttributes() > 0) {
      JPanel padder = new JPanel();
      JPanel padd2 = new JPanel();
      
      /*    if (m_splitListener != null) {
	    m_plotInstances.randomize(new Random());
	    } */

      m_heights = new int[ins.numInstances()];

      m_cIndex = ins.numAttributes() - 1;
      for (int noa = 0; noa < ins.numInstances(); noa++) {
	m_heights[noa] = (int)(Math.random() * 19);
      }
      m_span.setPreferredSize(new Dimension(m_span.getPreferredSize().width, 
					    (m_cIndex + 1) * 20));
      m_span.setMaximumSize(new Dimension(m_span.getMaximumSize().width, 
					  (m_cIndex + 1) * 20));
      AttributeSpacing tmp;
      
      GridBagLayout gb = new GridBagLayout();
      GridBagLayout gb2 = new GridBagLayout();
      GridBagConstraints constraints = new GridBagConstraints();
      


      padder.setLayout(gb);
      m_span.setLayout(gb2);
      constraints.anchor = GridBagConstraints.CENTER;
      constraints.gridx=0;constraints.gridy=0;constraints.weightx=5;
      constraints.fill = GridBagConstraints.HORIZONTAL;
      constraints.gridwidth=1;constraints.gridheight=1;
      constraints.insets = new Insets(0, 0, 0, 0);
      padder.add(m_span, constraints);
      constraints.gridx=0;constraints.gridy=1;constraints.weightx=5;
      constraints.fill = GridBagConstraints.BOTH;
      constraints.gridwidth=1;constraints.gridheight=1;constraints.weighty=5;
      constraints.insets = new Insets(0, 0, 0, 0);
      padder.add(padd2, constraints);
      constraints.weighty=0;
      setViewportView(padder);
      //getViewport().setLayout(null);
      //m_span.setMinimumSize(new Dimension(100, (m_cIndex + 1) * 24));
      //m_span.setSize(100, (m_cIndex + 1) * 24);
      constraints.anchor = GridBagConstraints.CENTER;
      constraints.gridx=0;constraints.gridy=0;constraints.weightx=5;
      constraints.fill = GridBagConstraints.HORIZONTAL;
      constraints.gridwidth=1;constraints.gridheight=1;constraints.weighty=5;
      constraints.insets = new Insets(2,20,2,4);

      for (int noa = 0; noa < ins.numAttributes(); noa++) {
	tmp = new AttributeSpacing(ins.attribute(noa), noa);
	 
	constraints.gridy = noa;
	m_span.add(tmp, constraints);
      }
    }
  }
    
  /**
   * shows which bar is the current x attribute.
   * @param x The attributes index.
   */
  public void setX(int x) {
    if (m_span != null) {
      m_xIndex = x;
      m_span.repaint();
    }
  }
    
  /**
   * shows which bar is the current y attribute.
   * @param y The attributes index.
   */
  public void setY(int y) {
    if (m_span != null) {
      m_yIndex = y;
      m_span.repaint();
    }
  }

  /**
   * Main method for testing this class.
   * @param args first argument should be an arff file. Second argument
   * can be an optional class col
   */
  public static void main(String [] args) {
    try {
      if (args.length < 1) {
	System.err.println("Usage : weka.gui.visualize.AttributePanel "
			   +"<dataset> [class col]");
	System.exit(1);
      }
      final javax.swing.JFrame jf = 
	new javax.swing.JFrame("Weka Explorer: Attribute");
      jf.setSize(100,100);
      jf.getContentPane().setLayout(new BorderLayout());
      final AttributePanel p2 = new AttributePanel();
      p2.addAttributePanelListener(new AttributePanelListener() {
	  public void attributeSelectionChange(AttributePanelEvent e) {
	    if (e.m_xChange) {
	      System.err.println("X index changed to : "+e.m_indexVal);
	    } else {
	      System.err.println("Y index changed to : "+e.m_indexVal);
	    }
	  }
	});
      jf.getContentPane().add(p2, BorderLayout.CENTER);
      jf.addWindowListener(new java.awt.event.WindowAdapter() {
	  public void windowClosing(java.awt.event.WindowEvent e) {
	    jf.dispose();
	    System.exit(0);
	  }
	});
      if (args.length >= 1) {
	System.err.println("Loading instances from " + args[0]);
	java.io.Reader r = new java.io.BufferedReader(
			   new java.io.FileReader(args[0]));
	Instances i = new Instances(r);
	i.setClassIndex(i.numAttributes()-1);
	p2.setInstances(i);
      }
      if (args.length > 1) {
	p2.setCindex((Integer.parseInt(args[1]))-1);
      } else {
	p2.setCindex(0);
      }
      jf.setVisible(true);
    } catch (Exception ex) {
       ex.printStackTrace();
       System.err.println(ex.getMessage());
     }
  }
}
