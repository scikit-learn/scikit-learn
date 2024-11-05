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
 *    ClassPanel.java
 *    Copyright (C) 2000 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.visualize;

import weka.core.FastVector;
import weka.core.Instances;
import weka.core.Utils;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;

import javax.swing.JColorChooser;
import javax.swing.JLabel;
import javax.swing.JPanel;

/**
 * This panel displays coloured labels for nominal attributes and a spectrum
 * for numeric attributes. It can also be told to colour on the basis
 * of an array of doubles (this can be useful for displaying coloured labels
 * that correspond to a clusterers predictions).
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @author Malcolm Ware (mfw4@cs.waikato.ac.nz)
 * @version $Revision: 5086 $
 */
public class ClassPanel
  extends JPanel {

  /** for serialization */
  private static final long serialVersionUID = -7969401840501661430L;
    
  /** True when the panel has been enabled (ie after 
      setNumeric or setNominal has been called */
  private boolean m_isEnabled = false;

  /** True if the colouring attribute is numeric */
  private boolean m_isNumeric = false;
    
  /** The height of the spectrum for numeric class */
  private final int m_spectrumHeight = 5;

  /**  The maximum value for the colouring attribute */
  private double m_maxC;
    
  /** The minimum value for the colouring attribute */
  private double m_minC;

  /** The size of the ticks */
  private final int m_tickSize = 5;

  /** Font metrics */
  private FontMetrics m_labelMetrics = null;

  /** The font used in labeling */
  private Font m_labelFont = null;

  /** The amount of space to leave either side of the legend */ 
  private int m_HorizontalPad=0;

  /** The precision with which to display real values */
  private int m_precisionC;

  /** Field width for numeric values */
  private int m_fieldWidthC;

  /** The old width. */
  private int m_oldWidth = -9000;
    
  /** Instances being plotted */
  private Instances m_Instances=null;

  /** Index of the colouring attribute */
  private int m_cIndex;

  /** the list of colours to use for colouring nominal attribute labels */
  private FastVector m_colorList;

  /** An optional list of Components that use the colour list
      maintained by this class. If the user changes a colour
      using the colour chooser, then these components need to
      be repainted in order to display the change */
  private FastVector m_Repainters = new FastVector();

  /** An optional list of listeners who want to know when a colour
      changes. Listeners are notified via an ActionEvent */
  private FastVector m_ColourChangeListeners = new FastVector();

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
   *  if set, it allows this panel to steer away from setting up
   * a color in the color list that is equal to the background color
   */
  protected Color m_backgroundColor = null;

  /** Inner Inner class used to create labels for nominal attributes
   * so that there color can be changed.
   */
  private class NomLabel
    extends JLabel {

    /** for serialization */
    private static final long serialVersionUID = -4686613106474820655L;

    private int m_index = 0;

    /** 
     * Creates a label with its name and class index value.
     * @param name The name of the nominal class value.
     * @param id The index value for that class value.
     */
    public NomLabel(String name, int id) {
      super(name);
      m_index = id;

      this.addMouseListener(new MouseAdapter() {
	  public void mouseClicked(MouseEvent e) {
	      
	    if ((e.getModifiers() & e.BUTTON1_MASK) == e.BUTTON1_MASK) {
	      Color tmp = JColorChooser.showDialog
		(ClassPanel.this, "Select new Color", 
		 (Color)m_colorList.elementAt(m_index));
		
	      if (tmp != null) {
		m_colorList.setElementAt(tmp, m_index);
		m_oldWidth = -9000;
		ClassPanel.this.repaint();
		if (m_Repainters.size() > 0) {
		  for (int i=0;i<m_Repainters.size();i++) {
		    ((Component)(m_Repainters.elementAt(i))).repaint();
		  }
		}
		
		if (m_ColourChangeListeners.size() > 0) {
		  for (int i = 0; i < m_ColourChangeListeners.size(); i++) {
		    ((ActionListener)(m_ColourChangeListeners.elementAt(i))).
		      actionPerformed(new ActionEvent(this, 0, ""));
		  }
		}
	      }
	    }
	  }
	});
    }
  }
  
  public ClassPanel() {
    this(null);
  }

  public ClassPanel(Color background) {
    m_backgroundColor = background;
    
    /** Set up some default colours */
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
   * Adds a component that will need to be repainted if the user
   * changes the colour of a label.
   * @param c the component to be repainted
   */
  public void addRepaintNotify(Component c) {
    m_Repainters.addElement(c);
  }

  /**
   * Add an action listener that will be notified if the user changes the
   * colour of a label
   *
   * @param a an <code>ActionListener</code> value
   */
  public void addActionListener(ActionListener a) {
    m_ColourChangeListeners.addElement(a);
  }

  /**
   * Set up fonts and font metrics
   * @param gx the graphics context
   */
  private void setFonts(Graphics gx) {
    if (m_labelMetrics == null) {
      m_labelFont = new Font("Monospaced", Font.PLAIN, 12);
      m_labelMetrics = gx.getFontMetrics(m_labelFont);
      int hf = m_labelMetrics.getAscent();
      if (this.getHeight() < (3*hf)) {
	m_labelFont = new Font("Monospaced",Font.PLAIN,11);
	m_labelMetrics = gx.getFontMetrics(m_labelFont);
      }
    }
    gx.setFont(m_labelFont);
  }
    
  /**
   * Enables the panel
   * @param e true to enable the panel
   */
  public void setOn(boolean e) {
    m_isEnabled = e;
  }

  /**
   * Set the instances.
   * @param insts the instances
   */
  public void setInstances(Instances insts) {
    m_Instances = insts;
  }

  /**
   * Set the index of the attribute to display coloured labels for
   * @param cIndex the index of the attribute to display coloured labels for
   */
  public void setCindex(int cIndex) {
    if (m_Instances.numAttributes() > 0) {
      m_cIndex = cIndex;
      if (m_Instances.attribute(m_cIndex).isNumeric()) {
	setNumeric();
      } else {
	if (m_Instances.attribute(m_cIndex).numValues() > m_colorList.size()) {
	  extendColourMap();
	}
	setNominal();
      }
    }
  }

  /**
   * Extends the list of colours if a new attribute with more values than
   * the previous one is chosen
   */
  private void extendColourMap() {
    if (m_Instances.attribute(m_cIndex).isNominal() || m_Instances.attribute(m_cIndex).isRanking()) {
      for (int i = m_colorList.size(); 
	   i < m_Instances.attribute(m_cIndex).numValues();
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
    
  protected void setDefaultColourList(Color[] list) {
    m_DefaultColors = list;
  }

  /**
   * Set a list of colours to use for colouring labels
   * @param cols a list containing java.awt.Colors
   */
  public void setColours(FastVector cols) {
    m_colorList = cols;
  }
    
  /**
   * Sets the legend to be for a nominal variable
   */
  protected void setNominal() {
    m_isNumeric = false;
    m_HorizontalPad = 0;
    setOn(true);
    m_oldWidth = -9000;
     
    this.repaint();
  }

  /**
   * Sets the legend to be for a numeric variable
   */
  protected void setNumeric() {
    m_isNumeric = true;
    /*      m_maxC = mxC;
	    m_minC = mnC; */

    double min=Double.POSITIVE_INFINITY;
    double max=Double.NEGATIVE_INFINITY;
    double value;

    for (int i=0;i<m_Instances.numInstances();i++) {
      if (!m_Instances.instance(i).isMissing(m_cIndex)) {
	value = m_Instances.instance(i).value(m_cIndex);
	if (value < min) {
	  min = value;
	}
	if (value > max) {
	  max = value;
	}
      }
    }
     
    // handle case where all values are missing
    if (min == Double.POSITIVE_INFINITY) min = max = 0.0;

    m_minC = min; m_maxC = max;

    int whole = (int)Math.abs(m_maxC);
    double decimal = Math.abs(m_maxC) - whole;
    int nondecimal;
    nondecimal = (whole > 0) 
      ? (int)(Math.log(whole) / Math.log(10))
      : 1;
    
    m_precisionC = (decimal > 0) 
      ? (int)Math.abs(((Math.log(Math.abs(m_maxC)) / 
				      Math.log(10))))+2
      : 1;
    if (m_precisionC > VisualizeUtils.MAX_PRECISION) {
      m_precisionC = 1;
    }

    String maxStringC = Utils.doubleToString(m_maxC,
					     nondecimal+1+m_precisionC
					     ,m_precisionC);
    if (m_labelMetrics != null) {
      m_HorizontalPad = m_labelMetrics.stringWidth(maxStringC);
    }

    whole = (int)Math.abs(m_minC);
    decimal = Math.abs(m_minC) - whole;
    nondecimal = (whole > 0) 
      ? (int)(Math.log(whole) / Math.log(10))
      : 1;
    
     m_precisionC = (decimal > 0) 
       ? (int)Math.abs(((Math.log(Math.abs(m_minC)) / 
				      Math.log(10))))+2
      : 1;
     if (m_precisionC > VisualizeUtils.MAX_PRECISION) {
       m_precisionC = 1;
     }
    
     maxStringC = Utils.doubleToString(m_minC,
				       nondecimal+1+m_precisionC
				       ,m_precisionC);
     if (m_labelMetrics != null) {
       if (m_labelMetrics.stringWidth(maxStringC) > m_HorizontalPad) {
	 m_HorizontalPad = m_labelMetrics.stringWidth(maxStringC);
       }
     }

    setOn(true);
    this.repaint();
  }
    
  /**
   * Renders the legend for a nominal colouring attribute
   * @param gx the graphics context
   */
  protected void paintNominal(Graphics gx) {
    setFonts(gx);



    int numClasses;

    numClasses = m_Instances.attribute(m_cIndex).numValues();

    int maxLabelLen = 0;
    int idx=0;
    int legendHeight;
    int w = this.getWidth();
    int hf = m_labelMetrics.getAscent();


    for (int i=0;i<numClasses;i++) {
      if (m_Instances.attribute(m_cIndex).value(i).length() > 
	  maxLabelLen) {
	maxLabelLen = m_Instances.
	  attribute(m_cIndex).value(i).length();
	idx = i;
      }
    }
      
    maxLabelLen = m_labelMetrics.stringWidth(m_Instances.
					     attribute(m_cIndex).value(idx));
    

    if (((w-(2*m_HorizontalPad))/(maxLabelLen+5)) >= numClasses) {
      legendHeight = 1;
    } else {
      legendHeight = 2;
    }
	
    int x = m_HorizontalPad;
    int y = 1 + hf;

    // do the first row
    int ci, mp;
    Color pc;
    int numToDo = ((legendHeight==1) ? numClasses : (numClasses/2));
    for (int i=0;i<numToDo;i++) {
     
      gx.setColor((Color)m_colorList.elementAt(i));
      // can we fit the full label or will each need to be trimmed?
      if ((numToDo * maxLabelLen) > (w-(m_HorizontalPad*2))) {
	String val;
	val = m_Instances.attribute(m_cIndex).value(i);

	int sw = m_labelMetrics.stringWidth(val);
	int rm=0;
	// truncate string if necessary
	if (sw > ((w-(m_HorizontalPad*2)) / (numToDo))) {
	  int incr = (sw / val.length());
	  rm = (sw -  ((w-(m_HorizontalPad*2)) / numToDo)) / incr;
	  if (rm <= 0) {
	    rm = 0;
	  }
	  if (rm >= val.length()) {
	    rm = val.length() - 1;
	  }
	  val = val.substring(0,val.length()-rm);
	  sw = m_labelMetrics.stringWidth(val);
	}
	NomLabel jj = new NomLabel(val, i);
	jj.setFont(gx.getFont());

	jj.setSize(m_labelMetrics.stringWidth(jj.getText()),
		   m_labelMetrics.getAscent() + 4);
	this.add(jj);
	jj.setLocation(x, y);
	jj.setForeground((Color)m_colorList.
			 elementAt(i % m_colorList.size()));

	x += sw + 2;
      } else {
	
	NomLabel jj;
	jj = new NomLabel(m_Instances.attribute(m_cIndex).value(i), i);

	jj.setFont(gx.getFont());

	jj.setSize(m_labelMetrics.stringWidth(jj.getText()),
		   m_labelMetrics.getAscent() + 4);
	this.add(jj);
	jj.setLocation(x, y);
	jj.setForeground((Color)m_colorList.
			 elementAt(i % m_colorList.size()));

  

	x += ((w-(m_HorizontalPad*2)) / numToDo);
      }	  
    }

    x = m_HorizontalPad;
    y = 1+ hf + 5 +hf;
    for (int i=numToDo;i<numClasses;i++) {
      
      gx.setColor((Color)m_colorList.elementAt(i));
      if (((numClasses-numToDo+1) * maxLabelLen) > 
	  (w - (m_HorizontalPad*2))) {
	String val;
	val = m_Instances.attribute(m_cIndex).value(i);

	int sw = m_labelMetrics.stringWidth(val);
	int rm=0;
	// truncate string if necessary
	if (sw > ((w-(m_HorizontalPad*2)) / (numClasses-numToDo+1))) {
	  int incr = (sw / val.length());
	  rm = (sw -  ((w-(m_HorizontalPad*2)) / (numClasses-numToDo))) 
	    / incr;
	  if (rm <= 0) {
	    rm = 0;
	  }
	  if (rm >= val.length()) {
	    rm = val.length() - 1;
	  }
	  val = val.substring(0,val.length()-rm);
	  sw = m_labelMetrics.stringWidth(val);
	}
	//this is the clipped string
	NomLabel jj = new NomLabel(val, i);
	jj.setFont(gx.getFont());

	jj.setSize(m_labelMetrics.stringWidth(jj.getText()),
		   m_labelMetrics.getAscent() + 4);

	this.add(jj);
	jj.setLocation(x, y);
	jj.setForeground((Color)m_colorList.
			 elementAt(i % m_colorList.size()));
	
	x += sw +2;
      } else {
	//this is the full string
	NomLabel jj;
	jj = new NomLabel(m_Instances.attribute(m_cIndex).value(i), i);

	jj.setFont(gx.getFont());

	jj.setSize(m_labelMetrics.stringWidth(jj.getText()),
		   m_labelMetrics.getAscent() + 4);
	this.add(jj);
	jj.setLocation(x, y);
	jj.setForeground((Color)m_colorList.
			 elementAt(i % m_colorList.size()));

	x += ((w - (m_HorizontalPad*2)) / (numClasses-numToDo));
      }	  
    }

  }

  /**
   * Renders the legend for a numeric colouring attribute
   * @param gx the graphics context
   */
  protected void paintNumeric(Graphics gx) {

    setFonts(gx);
    if (m_HorizontalPad == 0) {
      setCindex(m_cIndex);
    }

    int w = this.getWidth();
    double rs = 15;
    double incr = 240.0 / (double)(w-(m_HorizontalPad*2));
    int hf = m_labelMetrics.getAscent();
      
    for (int i=m_HorizontalPad;i<
	   (w-m_HorizontalPad);i++) {
      Color c = new Color((int)rs,150,(int)(255-rs));
      gx.setColor(c);
      gx.drawLine(i,0,
		  i,0+m_spectrumHeight);
      rs += incr;
    }

    int whole = (int)Math.abs(m_maxC);
    double decimal = Math.abs(m_maxC) - whole;
    int nondecimal;
    nondecimal = (whole > 0) 
      ? (int)(Math.log(whole) / Math.log(10))
      : 1;
    
    m_precisionC = (decimal > 0) 
      ? (int)Math.abs(((Math.log(Math.abs(m_maxC)) / 
			Math.log(10))))+2
      : 1;
    if (m_precisionC > VisualizeUtils.MAX_PRECISION) {
      m_precisionC = 1;
    }

    String maxStringC = Utils.doubleToString(m_maxC,
					     nondecimal+1+m_precisionC
					     ,m_precisionC);

	
    int mswc = m_labelMetrics.stringWidth(maxStringC);
    int tmsc = mswc;
    if (w > (2 * tmsc)) {
      gx.setColor(Color.black);
      gx.drawLine(m_HorizontalPad,
		  (m_spectrumHeight+5),
		  w-m_HorizontalPad,
		  (m_spectrumHeight+5));

      gx.drawLine(w-m_HorizontalPad,
		  (m_spectrumHeight+5),
		  w-m_HorizontalPad,
		  (m_spectrumHeight+5+m_tickSize));

      gx.drawString(maxStringC, 
		    (w-m_HorizontalPad)-(mswc/2),
		    (m_spectrumHeight+5+m_tickSize+hf));

      gx.drawLine(m_HorizontalPad,
		  (m_spectrumHeight+5),
		  m_HorizontalPad,
		  (m_spectrumHeight+5+m_tickSize));

      whole = (int)Math.abs(m_minC);
      decimal = Math.abs(m_minC) - whole;
      nondecimal = (whole > 0) 
	? (int)(Math.log(whole) / Math.log(10))
	: 1;
      
      m_precisionC = (decimal > 0) 
	? (int)Math.abs(((Math.log(Math.abs(m_minC)) / 
			  Math.log(10))))+2
	: 1;

      if (m_precisionC > VisualizeUtils.MAX_PRECISION) {
	m_precisionC = 1;
      }
      
      maxStringC = Utils.doubleToString(m_minC,
					nondecimal+1+m_precisionC
					,m_precisionC);

      mswc = m_labelMetrics.stringWidth(maxStringC);
      gx.drawString(maxStringC, 
		    m_HorizontalPad-(mswc/2),
		    (m_spectrumHeight+5+m_tickSize+hf));

      // draw the middle value if there is space
      if (w > (3 * tmsc)) {
	double mid = m_minC+((m_maxC-m_minC)/2.0);
	gx.drawLine(m_HorizontalPad+((w-(2*m_HorizontalPad))/2),
		    (m_spectrumHeight+5),
		    m_HorizontalPad+((w-(2*m_HorizontalPad))/2),
		    (m_spectrumHeight+5+m_tickSize));

	whole = (int)Math.abs(mid);
	decimal = Math.abs(mid) - whole;
	nondecimal = (whole > 0) 
	  ? (int)(Math.log(whole) / Math.log(10))
	  : 1;
	
	m_precisionC = (decimal > 0) 
	  ? (int)Math.abs(((Math.log(Math.abs(mid)) / 
			    Math.log(10))))+2
	  : 1;
	if (m_precisionC > VisualizeUtils.MAX_PRECISION) {
	  m_precisionC = 1;
	}
	
	maxStringC = Utils.doubleToString(mid,
					  nondecimal+1+m_precisionC
					  ,m_precisionC);

	mswc = m_labelMetrics.stringWidth(maxStringC);
	gx.drawString(maxStringC,
		      m_HorizontalPad+((w-(2*m_HorizontalPad))/2)-(mswc/2),
		      (m_spectrumHeight+5+m_tickSize+hf));
      }
    }
  }

  /**
   * Renders this component
   * @param gx the graphics context
   */
  public void paintComponent(Graphics gx) {
    super.paintComponent(gx);
    if (m_isEnabled) {
      if (m_isNumeric) {
	m_oldWidth = -9000;   //done so that if change back to nom, it will
	//work
	this.removeAll();
	paintNumeric(gx);
      } else {
	if (m_Instances != null && 
	    m_Instances.numInstances() > 0 && 
	    m_Instances.numAttributes() > 0) {
	  if (m_oldWidth != this.getWidth()) {
	    this.removeAll();
	    m_oldWidth = this.getWidth();
	    paintNominal(gx);
	  }
	}
      }
    }
  }

  /**
   * Main method for testing this class.
   * @param args first argument must specify an arff file. Second can
   * specify an optional index to colour labels on
   */
  public static void main(String [] args) {
    try {
      if (args.length < 1) {
	System.err.println("Usage : weka.gui.visualize.ClassPanel <dataset> "
			   +"[class col]");
	System.exit(1);
      }
      final javax.swing.JFrame jf = 
	new javax.swing.JFrame("Weka Explorer: Class");
      jf.setSize(500,100);
      jf.getContentPane().setLayout(new BorderLayout());
      final ClassPanel p2 = new ClassPanel();
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
