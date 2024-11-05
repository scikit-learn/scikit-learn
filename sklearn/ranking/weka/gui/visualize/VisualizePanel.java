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
 *    VisualizePanel.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.visualize;

import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.labelranking.PreferenceAttribute;
import weka.gui.ExtensionFileFilter;
import weka.gui.Logger;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionAdapter;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.Writer;
import java.util.Random;

import javax.swing.BorderFactory;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.SwingConstants;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.filechooser.FileFilter;

/** 
 * This panel allows the user to visualize a dataset (and if provided) a
 * classifier's/clusterer's predictions in two dimensions.
 *
 * If the user selects a nominal attribute as the colouring attribute then
 * each point is drawn in a colour that corresponds to the discrete value
 * of that attribute for the instance. If the user selects a numeric
 * attribute to colour on, then the points are coloured using a spectrum
 * ranging from blue to red (low values to high).
 *
 * When a classifier's predictions are supplied they are plotted in one
 * of two ways (depending on whether the class is nominal or numeric).<br>
 * For nominal class: an error made by a classifier is plotted as a square
 * in the colour corresponding to the class it predicted.<br>
 * For numeric class: predictions are plotted as varying sized x's, where
 * the size of the x is related to the magnitude of the error.
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @author Malcolm Ware (mfw4@cs.waikato.ac.nz)
 * @version $Revision: 5747 $
 */
public class VisualizePanel
  extends PrintablePanel {

  /** for serialization */
  private static final long serialVersionUID = 240108358588153943L;

  /** Inner class to handle plotting */
  protected class PlotPanel
    extends PrintablePanel
    implements Plot2DCompanion {

    /** for serialization */
    private static final long serialVersionUID = -4823674171136494204L;

    /** The actual generic plotting panel */
    protected Plot2D m_plot2D = new Plot2D();

    /** The instances from the master plot */
    protected Instances m_plotInstances=null;

    /** The master plot */
    protected PlotData2D m_originalPlot=null;
    
    /** Indexes of the attributes to go on the x and y axis and the attribute
	to use for colouring and the current shape for drawing */
    protected int m_xIndex=0;
    protected int m_yIndex=0;
    protected int m_cIndex=0;
    protected int m_sIndex=0;

    /**the offsets of the axes once label metrics are calculated */
    private int m_XaxisStart=0;
    private int m_YaxisStart=0;
    private int m_XaxisEnd=0;
    private int m_YaxisEnd=0;

    /** True if the user is currently dragging a box. */
    private boolean m_createShape;
    
    /** contains all the shapes that have been drawn for these attribs */
    private FastVector m_shapes;

    /** contains the points of the shape currently being drawn. */
    private FastVector m_shapePoints;

    /** contains the position of the mouse (used for rubberbanding). */
    private Dimension m_newMousePos;

    /** Constructor */
    public PlotPanel() {
      this.setBackground(m_plot2D.getBackground());
      this.setLayout(new BorderLayout());
      this.add(m_plot2D, BorderLayout.CENTER);
      m_plot2D.setPlotCompanion(this);

      m_createShape = false;        
      m_shapes = null;////
      m_shapePoints = null;
      m_newMousePos = new Dimension();

      this.addMouseListener(new MouseAdapter() {
	  ///////      
	  public void mousePressed(MouseEvent e) {
	    if ((e.getModifiers() & MouseEvent.BUTTON1_MASK) == MouseEvent.BUTTON1_MASK) {
	      //
	      if (m_sIndex == 0) {
		//do nothing it will get dealt to in the clicked method
	      }
	      else if (m_sIndex == 1) {
		m_createShape = true;
		m_shapePoints = new FastVector(5);
		m_shapePoints.addElement(new Double(m_sIndex));
		m_shapePoints.addElement(new Double(e.getX()));
		m_shapePoints.addElement(new Double(e.getY()));
		m_shapePoints.addElement(new Double(e.getX()));
		m_shapePoints.addElement(new Double(e.getY()));
		//		Graphics g = PlotPanel.this.getGraphics();
		Graphics g = m_plot2D.getGraphics();
		g.setColor(Color.black);
		g.setXORMode(Color.white);
		g.drawRect(((Double)m_shapePoints.elementAt(1)).intValue(),
			   ((Double)m_shapePoints.elementAt(2)).intValue(),
			   ((Double)m_shapePoints.elementAt(3)).intValue() -
			   ((Double)m_shapePoints.elementAt(1)).intValue(), 
			   ((Double)m_shapePoints.elementAt(4)).intValue() -
			   ((Double)m_shapePoints.elementAt(2)).intValue());
		g.dispose();
	      }
	      //System.out.println("clicked");
	    }
	    //System.out.println("clicked");
	  }
	  //////
	  public void mouseClicked(MouseEvent e) {
	    
	    if ((m_sIndex == 2 || m_sIndex == 3) && 
		(m_createShape || 
		 (e.getModifiers() & MouseEvent.BUTTON1_MASK) == MouseEvent.BUTTON1_MASK)) {
	      if (m_createShape) {
		//then it has been started already.

		Graphics g = m_plot2D.getGraphics();
		g.setColor(Color.black);
		g.setXORMode(Color.white);
		if ((e.getModifiers() & MouseEvent.BUTTON1_MASK) == MouseEvent.BUTTON1_MASK &&
                    !e.isAltDown()) {
		  m_shapePoints.addElement(new 
		    Double(m_plot2D.convertToAttribX(e.getX())));
		  
		  m_shapePoints.addElement(new 
		    Double(m_plot2D.convertToAttribY(e.getY())));
		  
		  m_newMousePos.width = e.getX();
		  m_newMousePos.height = e.getY();
		  g.drawLine((int)Math.ceil
			     (m_plot2D.convertToPanelX
			      (((Double)m_shapePoints.
				elementAt(m_shapePoints.size() - 2)).
			       doubleValue())),
			     
			     (int)Math.ceil
			     (m_plot2D.convertToPanelY
			      (((Double)m_shapePoints.
				elementAt(m_shapePoints.size() - 1)).
			       doubleValue())),
			     m_newMousePos.width, m_newMousePos.height);
		  
		}
		else if (m_sIndex == 3) {
		  //then extend the lines to infinity 
		  //(100000 or so should be enough).
		  //the area is selected by where the user right clicks 
		  //the mouse button
		  
		  m_createShape = false;
		  if (m_shapePoints.size() >= 5) {
		    double cx = Math.ceil
		      (m_plot2D.convertToPanelX
		       (((Double)m_shapePoints.elementAt
			 (m_shapePoints.size() - 4)).doubleValue()));
		    
		    double cx2 = Math.ceil
		      (m_plot2D.convertToPanelX
		       (((Double)m_shapePoints.elementAt
			 (m_shapePoints.size() - 2)).doubleValue())) - 
		      cx;
		    
		    cx2 *= 50000;
		    
		    double cy = Math.ceil
		      (m_plot2D.
		       convertToPanelY(((Double)m_shapePoints.
					elementAt(m_shapePoints.size() - 3)).
				       doubleValue()));
		    double cy2 = Math.ceil
		      (m_plot2D.convertToPanelY(((Double)m_shapePoints.
					  elementAt(m_shapePoints.size() - 1)).
					  doubleValue())) - cy;
		    cy2 *= 50000;
			    
		    
		    double cxa = Math.ceil(m_plot2D.convertToPanelX
					   (((Double)m_shapePoints.
					     elementAt(3)).
					    doubleValue()));
		    double cxa2 = Math.ceil(m_plot2D.convertToPanelX
					    (((Double)m_shapePoints.
					      elementAt(1)).
					     doubleValue())) - cxa;
		    cxa2 *= 50000;
		    
		    
		    double cya = Math.ceil
		      (m_plot2D.convertToPanelY
		       (((Double)m_shapePoints.elementAt(4)).
			doubleValue()));
		    double cya2 = Math.ceil
		      (m_plot2D.convertToPanelY
		       (((Double)m_shapePoints.elementAt(2)).
			doubleValue())) - cya;
		    
		    cya2 *= 50000;
		    
		    m_shapePoints.setElementAt
		      (new Double(m_plot2D.convertToAttribX(cxa2 + cxa)), 1);
		    
		    m_shapePoints.setElementAt
		      (new Double(m_plot2D.convertToAttribY(cy2 + cy)), 
		       m_shapePoints.size() - 1);
		    
		    m_shapePoints.setElementAt
		      (new Double(m_plot2D.convertToAttribX(cx2 + cx)), 
		       m_shapePoints.size() - 2);
		    
		    m_shapePoints.setElementAt
		      (new Double(m_plot2D.convertToAttribY(cya2 + cya)), 2);
		    
		    
		    //determine how infinity line should be built
		    
		    cy = Double.POSITIVE_INFINITY;
		    cy2 = Double.NEGATIVE_INFINITY;
		    if (((Double)m_shapePoints.elementAt(1)).
			doubleValue() > 
			((Double)m_shapePoints.elementAt(3)).
			doubleValue()) {
		      if (((Double)m_shapePoints.elementAt(2)).
			  doubleValue() == 
			  ((Double)m_shapePoints.elementAt(4)).
			  doubleValue()) {
			cy = ((Double)m_shapePoints.elementAt(2)).
			  doubleValue();
		      }
		    }
		    if (((Double)m_shapePoints.elementAt
			 (m_shapePoints.size() - 2)).doubleValue() > 
			((Double)m_shapePoints.elementAt
			 (m_shapePoints.size() - 4)).doubleValue()) {
		      if (((Double)m_shapePoints.elementAt
			   (m_shapePoints.size() - 3)).
			  doubleValue() == 
			  ((Double)m_shapePoints.elementAt
			   (m_shapePoints.size() - 1)).doubleValue()) {
			cy2 = ((Double)m_shapePoints.lastElement()).
			  doubleValue();
		      }
		    }
		    m_shapePoints.addElement(new Double(cy));
		    m_shapePoints.addElement(new Double(cy2));
		    
		    if (!inPolyline(m_shapePoints, m_plot2D.convertToAttribX
				    (e.getX()), 
				    m_plot2D.convertToAttribY(e.getY()))) {
		      Double tmp = (Double)m_shapePoints.
			elementAt(m_shapePoints.size() - 2);
		      m_shapePoints.setElementAt
			(m_shapePoints.lastElement(), 
			 m_shapePoints.size() - 2);
		      m_shapePoints.setElementAt
			(tmp, m_shapePoints.size() - 1);
		    }
		    
		    if (m_shapes == null) {
		      m_shapes = new FastVector(4);
		    }
		    m_shapes.addElement(m_shapePoints);

		    m_submit.setText("Submit");
		    m_submit.setActionCommand("Submit");
		    
		    m_submit.setEnabled(true);
		  }
		  
		  m_shapePoints = null;
		  PlotPanel.this.repaint();
		  
		}
		else {
		  //then close the shape
		  m_createShape = false;
		  if (m_shapePoints.size() >= 7) {
		    m_shapePoints.addElement(m_shapePoints.elementAt(1));
		    m_shapePoints.addElement(m_shapePoints.elementAt(2));
		    if (m_shapes == null) {
		      m_shapes = new FastVector(4);
		    }
		    m_shapes.addElement(m_shapePoints);
			   
		    m_submit.setText("Submit");
		    m_submit.setActionCommand("Submit");
		    
		    m_submit.setEnabled(true);
		  }
		  m_shapePoints = null;
		  PlotPanel.this.repaint();
		}
		g.dispose();
		//repaint();
	      }
	      else if ((e.getModifiers() & MouseEvent.BUTTON1_MASK) == MouseEvent.BUTTON1_MASK) {
		//then this is the first point
		m_createShape = true;
		m_shapePoints = new FastVector(17);
		m_shapePoints.addElement(new Double(m_sIndex));
		m_shapePoints.addElement(new 
		  Double(m_plot2D.convertToAttribX(e.getX()))); //the new point
		m_shapePoints.addElement(new 
		  Double(m_plot2D.convertToAttribY(e.getY())));
		m_newMousePos.width = e.getX();      //the temp mouse point
		m_newMousePos.height = e.getY();

		Graphics g = m_plot2D.getGraphics();
		g.setColor(Color.black);
		g.setXORMode(Color.white);
		g.drawLine((int)Math.ceil
			   (m_plot2D.convertToPanelX(((Double)m_shapePoints.
					     elementAt(1)).doubleValue())),
			   (int)Math.ceil
			   (m_plot2D.convertToPanelY(((Double)m_shapePoints.
					     elementAt(2)).doubleValue())),
			   m_newMousePos.width, m_newMousePos.height);
		g.dispose();
	      }
	    }
	    else {
	      if ((e.getModifiers() & InputEvent.BUTTON1_MASK) == 
		  InputEvent.BUTTON1_MASK) {
		
		m_plot2D.searchPoints(e.getX(),e.getY(), false);
	      } else {
		m_plot2D.searchPoints(e.getX(), e.getY(), true);
	      }
	    }
	  }
	  
	  /////////             
	  public void mouseReleased(MouseEvent e) {

	    if (m_createShape) {
	      if (((Double)m_shapePoints.elementAt(0)).intValue() == 1) {
		m_createShape = false;
		Graphics g = m_plot2D.getGraphics();
		g.setColor(Color.black);
		g.setXORMode(Color.white);
		g.drawRect(((Double)m_shapePoints.elementAt(1)).
			   intValue(), 
			   ((Double)m_shapePoints.elementAt(2)).intValue(),
			   ((Double)m_shapePoints.elementAt(3)).intValue() -
			   ((Double)m_shapePoints.elementAt(1)).intValue(), 
			   ((Double)m_shapePoints.elementAt(4)).intValue() -
			   ((Double)m_shapePoints.elementAt(2)).intValue());
		
		g.dispose();
		if (checkPoints(((Double)m_shapePoints.elementAt(1)).
				doubleValue(), 
				((Double)m_shapePoints.elementAt(2)).
				doubleValue()) &&
		    checkPoints(((Double)m_shapePoints.elementAt(3)).
				doubleValue(), 
				((Double)m_shapePoints.elementAt(4)).
				doubleValue())) {
		  //then the points all land on the screen
		  //now do special check for the rectangle
		  if (((Double)m_shapePoints.elementAt(1)).doubleValue() <
		      ((Double)m_shapePoints.elementAt(3)).doubleValue() 
		      &&
		      ((Double)m_shapePoints.elementAt(2)).doubleValue() <
		      ((Double)m_shapePoints.elementAt(4)).doubleValue()) {
		    //then the rectangle is valid
		    if (m_shapes == null) {
		      m_shapes = new FastVector(2);
		    }
		    m_shapePoints.setElementAt(new 
		      Double(m_plot2D.convertToAttribX(((Double)m_shapePoints.
					       elementAt(1)).
					      doubleValue())), 1);
		    m_shapePoints.setElementAt(new 
		      Double(m_plot2D.convertToAttribY(((Double)m_shapePoints.
					       elementAt(2)).
					      doubleValue())), 2);
		    m_shapePoints.setElementAt(new 
		      Double(m_plot2D.convertToAttribX(((Double)m_shapePoints.
					       elementAt(3)).
					      doubleValue())), 3);
		    m_shapePoints.setElementAt(new 
		      Double(m_plot2D.convertToAttribY(((Double)m_shapePoints.
					       elementAt(4)).
					      doubleValue())), 4);
		    
		    m_shapes.addElement(m_shapePoints);
		    
		    m_submit.setText("Submit");
		    m_submit.setActionCommand("Submit");
		    
		    m_submit.setEnabled(true);

		    PlotPanel.this.repaint();
		  }
		}
		m_shapePoints = null;
	      }
	    }
	  }
	});
      
      this.addMouseMotionListener(new MouseMotionAdapter() {
	  public void mouseDragged(MouseEvent e) {
	    //check if the user is dragging a box
	    if (m_createShape) {
	      if (((Double)m_shapePoints.elementAt(0)).intValue() == 1) {
		Graphics g = m_plot2D.getGraphics();
		g.setColor(Color.black);
		g.setXORMode(Color.white);
		g.drawRect(((Double)m_shapePoints.elementAt(1)).intValue(), 
			   ((Double)m_shapePoints.elementAt(2)).intValue(),
			   ((Double)m_shapePoints.elementAt(3)).intValue() -
			   ((Double)m_shapePoints.elementAt(1)).intValue(), 
			   ((Double)m_shapePoints.elementAt(4)).intValue() -
			   ((Double)m_shapePoints.elementAt(2)).intValue());
		
		m_shapePoints.setElementAt(new Double(e.getX()), 3);
		m_shapePoints.setElementAt(new Double(e.getY()), 4);
		
		g.drawRect(((Double)m_shapePoints.elementAt(1)).intValue(), 
			   ((Double)m_shapePoints.elementAt(2)).intValue(),
			   ((Double)m_shapePoints.elementAt(3)).intValue() -
			   ((Double)m_shapePoints.elementAt(1)).intValue(), 
			   ((Double)m_shapePoints.elementAt(4)).intValue() -
			   ((Double)m_shapePoints.elementAt(2)).intValue());
		g.dispose();
	      }
	    }
	  }
	  
	  public void mouseMoved(MouseEvent e) {
	    if (m_createShape) {
	      if (((Double)m_shapePoints.elementAt(0)).intValue() == 2 || 
		  ((Double)m_shapePoints.elementAt(0)).intValue() == 3) {
		Graphics g = m_plot2D.getGraphics();
		g.setColor(Color.black);
		g.setXORMode(Color.white);
		g.drawLine((int)Math.ceil(m_plot2D.convertToPanelX
					  (((Double)m_shapePoints.elementAt
					    (m_shapePoints.size() - 2)).
					   doubleValue())),
			   (int)Math.ceil(m_plot2D.convertToPanelY
					  (((Double)m_shapePoints.elementAt
					    (m_shapePoints.size() - 1)).
					   doubleValue())),
			   m_newMousePos.width, m_newMousePos.height);
		
		m_newMousePos.width = e.getX();
		m_newMousePos.height = e.getY();
		
		g.drawLine((int)Math.ceil(m_plot2D.convertToPanelX
					  (((Double)m_shapePoints.elementAt
					    (m_shapePoints.size() - 2)).
					   doubleValue())),
			   (int)Math.ceil(m_plot2D.convertToPanelY
					  (((Double)m_shapePoints.elementAt
					    (m_shapePoints.size() - 1)).
					   doubleValue())),
			   m_newMousePos.width, m_newMousePos.height);
		g.dispose();
	      }
	    }
	  }
	});
      
      m_submit.addActionListener(new ActionListener() {
	  public void actionPerformed(ActionEvent e) {
	 
	    if (e.getActionCommand().equals("Submit")) {
	      if (m_splitListener != null && m_shapes != null) {
		//then send the split to the listener
		Instances sub_set1 = new Instances(m_plot2D.getMasterPlot().
						   m_plotInstances, 500);
		Instances sub_set2 = new Instances(m_plot2D.getMasterPlot().
						   m_plotInstances, 500);
		
		if (m_plot2D.getMasterPlot().
		    m_plotInstances != null) {
		  
		  for (int noa = 0 ; noa < m_plot2D.getMasterPlot().
			 m_plotInstances.numInstances(); noa++) {
		    if (!m_plot2D.getMasterPlot().
			m_plotInstances.instance(noa).isMissing(m_xIndex) &&
			!m_plot2D.getMasterPlot().
			m_plotInstances.instance(noa).isMissing(m_yIndex)){
		      
		      if (inSplit(m_plot2D.getMasterPlot().
				  m_plotInstances.instance(noa))) {
			sub_set1.add(m_plot2D.getMasterPlot().
				     m_plotInstances.instance(noa));
		      }
		      else {
			sub_set2.add(m_plot2D.getMasterPlot().
				     m_plotInstances.instance(noa));
		      }
		    }
		  }
		  FastVector tmp = m_shapes;
		  cancelShapes();
		  m_splitListener.userDataEvent(new 
		    VisualizePanelEvent(tmp, sub_set1, sub_set2, m_xIndex, 
					m_yIndex));
		}
	      }
	      else if (m_shapes != null && 
		       m_plot2D.getMasterPlot().m_plotInstances != null) { 
		Instances sub_set1 = new Instances(m_plot2D.getMasterPlot().
						   m_plotInstances, 500);
		int count = 0;
		for (int noa = 0 ; noa < m_plot2D.getMasterPlot().
		       m_plotInstances.numInstances(); noa++) {
		  if (inSplit(m_plot2D.getMasterPlot().
			      m_plotInstances.instance(noa))) {
		    sub_set1.add(m_plot2D.getMasterPlot().
				 m_plotInstances.instance(noa));
		    count++;
		  }
		  
		}

		int [] nSizes = null;
		int [] nTypes = null;
		int x = m_xIndex;
		int y = m_yIndex;

		if (m_originalPlot == null) {
		  //this sets these instances as the instances 
		  //to go back to.
		  m_originalPlot = m_plot2D.getMasterPlot();
		}

		if (count > 0) {
		  nTypes = new int[count];
		  nSizes = new int[count];
		  count = 0;
		  for (int noa = 0; noa < m_plot2D.getMasterPlot().
			 m_plotInstances.numInstances(); 
		       noa++) {
		    if (inSplit(m_plot2D.getMasterPlot().
				m_plotInstances.instance(noa))) {

		      nTypes[count] = m_plot2D.getMasterPlot().
			m_shapeType[noa];
		      nSizes[count] = m_plot2D.getMasterPlot().
			m_shapeSize[noa];
		      count++;
		    }
		  }
		}
		cancelShapes();

		PlotData2D newPlot = new PlotData2D(sub_set1);

		try {
		  newPlot.setShapeSize(nSizes);
		  newPlot.setShapeType(nTypes);
		
		  m_plot2D.removeAllPlots();
		  
		  VisualizePanel.this.addPlot(newPlot);
		} catch (Exception ex) {
		  System.err.println(ex);
		  ex.printStackTrace();
		}

		try {
		  VisualizePanel.this.setXIndex(x);
		  VisualizePanel.this.setYIndex(y);
		} catch(Exception er) {
		  System.out.println("Error : " + er);
		  //  System.out.println("Part of user input so had to" +
		  //		 " catch here");
		}
	      }
	    }
	    else if (e.getActionCommand().equals("Reset")) {
	      int x = m_xIndex;
	      int y = m_yIndex;

	      m_plot2D.removeAllPlots();
	      try {
		VisualizePanel.this.addPlot(m_originalPlot);
	      } catch (Exception ex) {
		System.err.println(ex);
		ex.printStackTrace();
	      }

	      try {
		VisualizePanel.this.setXIndex(x);
		VisualizePanel.this.setYIndex(y);
	      } catch(Exception er) {
		System.out.println("Error : " + er);
	      }
	    }
	  }  
	});

      m_cancel.addActionListener(new ActionListener() {
	  public void actionPerformed(ActionEvent e) {
	    cancelShapes();
	    PlotPanel.this.repaint();
	  }
	});
      ////////////
    }

    /**
     * Removes all the plots.
     */
    public void removeAllPlots() {
      m_plot2D.removeAllPlots();
      m_legendPanel.setPlotList(m_plot2D.getPlots());
    }
    
    /**
     * @return The FastVector containing all the shapes.
     */
    public FastVector getShapes() {
      
      return m_shapes;
    }
    
    /**
     * Sets the list of shapes to empty and also cancels
     * the current shape being drawn (if applicable).
     */
    public void cancelShapes() {
       
      if (m_splitListener == null) {
	m_submit.setText("Reset");
	m_submit.setActionCommand("Reset");

	if (m_originalPlot == null || 
	    m_originalPlot.m_plotInstances == m_plotInstances) {
	  m_submit.setEnabled(false);
	}
	else {
	  m_submit.setEnabled(true);
	}
      }
      else {
	m_submit.setEnabled(false);
      }
      
      m_createShape = false;
      m_shapePoints = null;
      m_shapes = null;
      this.repaint();
    }

    /**
     * This can be used to set the shapes that should appear.
     * @param v The list of shapes.
     */
    public void setShapes(FastVector v) {
      //note that this method should be fine for doubles,
      //but anything else that uses something other than doubles 
      //(or uneditable objects) could have unsafe copies.
      if (v != null) {
	FastVector temp;
	m_shapes = new FastVector(v.size());
	for (int noa = 0; noa < v.size(); noa++) {
	  temp = new FastVector(((FastVector)v.elementAt(noa)).size());
	  m_shapes.addElement(temp);
	  for (int nob = 0; nob < ((FastVector)v.elementAt(noa)).size()
		 ; nob++) {
	    
	    temp.addElement(((FastVector)v.elementAt(noa)).elementAt(nob));
	    
	  }
	}
      }
      else {
	m_shapes = null;
      }
      this.repaint();
    }
    
    /** 
     * This will check the values of the screen points passed and make sure 
     * that they land on the screen
     * @param x1 The x coord.
     * @param y1 The y coord.
     * @return true if the point would land on the screen
     */
    private boolean checkPoints(double x1, double y1) {
      if (x1 < 0 || x1 > this.getSize().width || y1 < 0 
	  || y1 > this.getSize().height) {
	return false;
      }
      return true;
    }
    
    /**
     * This will check if an instance is inside or outside of the current
     * shapes.
     * @param i The instance to check.
     * @return True if 'i' falls inside the shapes, false otherwise.
     */
    public boolean inSplit(Instance i) {
      //this will check if the instance lies inside the shapes or not
      
      if (m_shapes != null) {
	FastVector stmp;
	double x1, y1, x2, y2;
	for (int noa = 0; noa < m_shapes.size(); noa++) {
	  stmp = (FastVector)m_shapes.elementAt(noa);
	  if (((Double)stmp.elementAt(0)).intValue() == 1) {
	    //then rectangle
	    x1 = ((Double)stmp.elementAt(1)).doubleValue();
	    y1 = ((Double)stmp.elementAt(2)).doubleValue();
	    x2 = ((Double)stmp.elementAt(3)).doubleValue();
	    y2 = ((Double)stmp.elementAt(4)).doubleValue();
	    if (i.value(m_xIndex) >= x1 && i.value(m_xIndex) <= x2 &&
		i.value(m_yIndex) <= y1 && i.value(m_yIndex) >= y2) {
	      //then is inside split so return true;
	      return true;
	    }
	  }
	  else if (((Double)stmp.elementAt(0)).intValue() == 2) {
	    //then polygon
	    if (inPoly(stmp, i.value(m_xIndex), i.value(m_yIndex))) {
	      return true;
	    }
	  }
	  else if (((Double)stmp.elementAt(0)).intValue() == 3) {
	    //then polyline
	    if (inPolyline(stmp, i.value(m_xIndex), i.value(m_yIndex))) {
	      return true;
	    }
	  }
	}
      }
      return false;
    }
    
    /**
     * Checks to see if the coordinate passed is inside the ployline
     * passed, Note that this is done using attribute values and not there
     * respective screen values.
     * @param ob The polyline.
     * @param x The x coord.
     * @param y The y coord.
     * @return True if it falls inside the polyline, false otherwise.
     */
    private boolean inPolyline(FastVector ob, double x, double y) {
      //this works similar to the inPoly below except that
      //the first and last lines are treated as extending infinite in one 
      //direction and 
      //then infinitly in the x dirction their is a line that will 
      //normaly be infinite but
      //can be finite in one or both directions
      
      int countx = 0;
      double vecx, vecy;
      double change;
      double x1, y1, x2, y2;
      
      for (int noa = 1; noa < ob.size() - 4; noa+= 2) {
	y1 = ((Double)ob.elementAt(noa+1)).doubleValue();
	y2 = ((Double)ob.elementAt(noa+3)).doubleValue();
	x1 = ((Double)ob.elementAt(noa)).doubleValue();
	x2 = ((Double)ob.elementAt(noa+2)).doubleValue();
	
	//System.err.println(y1 + " " + y2 + " " + x1 + " " + x2);
	vecy = y2 - y1;
	vecx = x2 - x1;
	if (noa == 1 && noa == ob.size() - 6) {
	  //then do special test first and last edge
	  if (vecy != 0) {
	    change = (y - y1) / vecy;
	    if (vecx * change + x1 >= x) {
	      //then intersection
	      countx++;
	    }
	  }
	}
	else if (noa == 1) {
	  if ((y < y2 && vecy > 0) || (y > y2 && vecy < 0)) {
	    //now just determine intersection or not
	    change = (y - y1) / vecy;
	    if (vecx * change + x1 >= x) {
	      //then intersection on horiz
	      countx++;
	    }
	  }
	}
	else if (noa == ob.size() - 6) {
	  //then do special test on last edge
	  if ((y <= y1 && vecy < 0) || (y >= y1 && vecy > 0)) {
	    change = (y - y1) / vecy;
	    if (vecx * change + x1 >= x) {
	      countx++;
	    }
	  }
	}
	else if ((y1 <= y && y < y2) || (y2 < y && y <= y1)) {
	  //then continue tests.
	  if (vecy == 0) {
	    //then lines are parallel stop tests in 
	    //ofcourse it should never make it this far
	  }
	  else {
	    change = (y - y1) / vecy;
	    if (vecx * change + x1 >= x) {
	      //then intersects on horiz
	      countx++;
	    }
	  }
	}
      }
      
      //now check for intersection with the infinity line
      y1 = ((Double)ob.elementAt(ob.size() - 2)).doubleValue();
      y2 = ((Double)ob.elementAt(ob.size() - 1)).doubleValue();
      
      if (y1 > y2) {
	//then normal line
	if (y1 >= y && y > y2) {
	  countx++;
	}
      }
      else {
	//then the line segment is inverted
	if (y1 >= y || y > y2) {
	  countx++;
	}
      }
      
      if ((countx % 2) == 1) {
	return true;
      }
      else {
	return false;
      }
    }


    /**
     * This checks to see if The coordinate passed is inside
     * the polygon that was passed.
     * @param ob The polygon.
     * @param x The x coord.
     * @param y The y coord.
     * @return True if the coordinate is in the polygon, false otherwise.
     */
    private boolean inPoly(FastVector ob, double x, double y) {
      //brief on how this works
      //it draws a line horizontally from the point to the right (infinitly)
      //it then sees how many lines of the polygon intersect this, 
      //if it is even then the point is
      // outside the polygon if it's odd then it's inside the polygon
      int count = 0;
      double vecx, vecy;
      double change;
      double x1, y1, x2, y2;
      for (int noa = 1; noa < ob.size() - 2; noa += 2) {
	y1 = ((Double)ob.elementAt(noa+1)).doubleValue();
	y2 = ((Double)ob.elementAt(noa+3)).doubleValue();
	if ((y1 <= y && y < y2) || (y2 < y && y <= y1)) {
	  //then continue tests.
	  vecy = y2 - y1;
	  if (vecy == 0) {
	    //then lines are parallel stop tests for this line
	  }
	  else {
	    x1 = ((Double)ob.elementAt(noa)).doubleValue();
	    x2 = ((Double)ob.elementAt(noa+2)).doubleValue();
	    vecx = x2 - x1;
	    change = (y - y1) / vecy;
	    if (vecx * change + x1 >= x) {
	      //then add to count as an intersected line
	      count++;
	    }
	  }
	}
      }
      if ((count % 2) == 1) {
	//then lies inside polygon
	//System.out.println("in");
	return true;
      }
      else {
	//System.out.println("out");
	return false;
      }
      //System.out.println("WHAT?!?!?!?!!?!??!?!");
      //return false;
    }

    /**
     * Set level of jitter and repaint the plot using the new jitter value
     * @param j the level of jitter
     */
    public void setJitter(int j) {
      m_plot2D.setJitter(j);
    }

    /**
     * Set the index of the attribute to go on the x axis
     * @param x the index of the attribute to use on the x axis
     */
    public void setXindex(int x) {

      // this just ensures that the shapes get disposed of 
      //if the attribs change
      if (x != m_xIndex) {
	cancelShapes();
      }
      m_xIndex = x;
      m_plot2D.setXindex(x);
      if (m_showAttBars) {
	m_attrib.setX(x);
      }
      //      this.repaint();
    }
    
    /**
     * Set the index of the attribute to go on the y axis
     * @param y the index of the attribute to use on the y axis
     */
    public void setYindex(int y) {
    
      // this just ensures that the shapes get disposed of 
      //if the attribs change
      if (y != m_yIndex) {
	cancelShapes();
      }
      m_yIndex = y;
      m_plot2D.setYindex(y);
      if (m_showAttBars) {
	m_attrib.setY(y);
      }
      //      this.repaint();
    }

    /**
     * Set the index of the attribute to use for colouring
     * @param c the index of the attribute to use for colouring
     */
    public void setCindex(int c) {
      m_cIndex = c;
      m_plot2D.setCindex(c);
      if (m_showAttBars) {
	m_attrib.setCindex(c, m_plot2D.getMaxC(), m_plot2D.getMinC());
      }
      m_classPanel.setCindex(c);
      this.repaint();
    }

    /**
     * Set the index of the attribute to use for the shape.
     * @param s the index of the attribute to use for the shape
     */
    public void setSindex(int s) {
      if (s != m_sIndex) {
	m_shapePoints = null;
	m_createShape = false;
      }
      m_sIndex = s;
      this.repaint();
    }
    
    /**
     * Clears all existing plots and sets a new master plot
     * @param newPlot the new master plot
     * @exception Exception if plot could not be added
     */
    public void setMasterPlot(PlotData2D newPlot) throws Exception {
      m_plot2D.removeAllPlots();
      this.addPlot(newPlot);
    }

    /**
     * Adds a plot. If there are no plots so far this plot becomes
     * the master plot and, if it has a custom colour defined then
     * the class panel is removed.
     * @param newPlot the plot to add.
     * @exception Exception if plot could not be added
     */
    public void addPlot(PlotData2D newPlot) throws Exception {
      if (m_plot2D.getPlots().size() == 0) {
	m_plot2D.addPlot(newPlot);
	if (m_plotSurround.getComponentCount() > 1 && 
	    m_plotSurround.getComponent(1) == m_attrib &&
	    m_showAttBars) {
	  try {
	    m_attrib.setInstances(newPlot.m_plotInstances);
	    m_attrib.setCindex(0);m_attrib.setX(0); m_attrib.setY(0);
	  } catch (Exception ex) {
	    // more attributes than the panel can handle?
	    // Due to hard coded constraints in GridBagLayout
	    m_plotSurround.remove(m_attrib);
	    System.err.println("Warning : data contains more attributes "
			       +"than can be displayed as attribute bars.");
	    if (m_Log != null) {
	      m_Log.logMessage("Warning : data contains more attributes "
			       +"than can be displayed as attribute bars.");
	    }
	  }
	} else if (m_showAttBars) {
	  try {
	    m_attrib.setInstances(newPlot.m_plotInstances);
	    m_attrib.setCindex(0);m_attrib.setX(0); m_attrib.setY(0);
	    GridBagConstraints constraints = new GridBagConstraints();
	    constraints.fill = GridBagConstraints.BOTH;
	    constraints.insets = new Insets(0, 0, 0, 0);
	    constraints.gridx=4;constraints.gridy=0;constraints.weightx=1;
	    constraints.gridwidth=1;constraints.gridheight=1;
	    constraints.weighty=5;
	    m_plotSurround.add(m_attrib, constraints);
	  } catch (Exception ex) {
	    System.err.println("Warning : data contains more attributes "
			       +"than can be displayed as attribute bars.");
	    if (m_Log != null) {
	      m_Log.logMessage("Warning : data contains more attributes "
			       +"than can be displayed as attribute bars.");
	    }
	  }
	}
	m_classPanel.setInstances(newPlot.m_plotInstances);

	plotReset(newPlot.m_plotInstances, newPlot.getCindex());
	if (newPlot.m_useCustomColour && m_showClassPanel) {
	  VisualizePanel.this.remove(m_classSurround);
	  switchToLegend();
	  m_legendPanel.setPlotList(m_plot2D.getPlots());
	  m_ColourCombo.setEnabled(false);
	}
      } else  {
	if (!newPlot.m_useCustomColour && m_showClassPanel) {
	  VisualizePanel.this.add(m_classSurround, BorderLayout.SOUTH);
	  m_ColourCombo.setEnabled(true);
	}
	if (m_plot2D.getPlots().size() == 1) {
	  switchToLegend();
	}
	m_plot2D.addPlot(newPlot);
	m_legendPanel.setPlotList(m_plot2D.getPlots());
      }
    }

    /**
     * Remove the attibute panel and replace it with the legend panel
     */
    protected void switchToLegend() {

      if (m_plotSurround.getComponentCount() > 1 && 
	  m_plotSurround.getComponent(1) == m_attrib) {
	m_plotSurround.remove(m_attrib);
      }
	
      if (m_plotSurround.getComponentCount() > 1 &&
	  m_plotSurround.getComponent(1) == m_legendPanel) {
	return;
      }

      GridBagConstraints constraints = new GridBagConstraints();
      constraints.fill = GridBagConstraints.BOTH;
      constraints.insets = new Insets(0, 0, 0, 0);
      constraints.gridx=4;constraints.gridy=0;constraints.weightx=1;
      constraints.gridwidth=1;constraints.gridheight=1;
      constraints.weighty=5;
      m_plotSurround.add(m_legendPanel, constraints);
      setSindex(0);
      m_ShapeCombo.setEnabled(false);
    }
    
    protected void switchToBars() {
      if (m_plotSurround.getComponentCount() > 1 && 
          m_plotSurround.getComponent(1) == m_legendPanel) {
        m_plotSurround.remove(m_legendPanel);
      }
      
      if (m_plotSurround.getComponentCount() > 1 &&
          m_plotSurround.getComponent(1) == m_attrib) {
        return;
      }
      
      if (m_showAttBars) {
        try {
          m_attrib.setInstances(m_plot2D.getMasterPlot().m_plotInstances);
          m_attrib.setCindex(0);m_attrib.setX(0); m_attrib.setY(0);
          GridBagConstraints constraints = new GridBagConstraints();
          constraints.fill = GridBagConstraints.BOTH;
          constraints.insets = new Insets(0, 0, 0, 0);
          constraints.gridx=4;constraints.gridy=0;constraints.weightx=1;
          constraints.gridwidth=1;constraints.gridheight=1;
          constraints.weighty=5;
          m_plotSurround.add(m_attrib, constraints);
        } catch (Exception ex) {
          System.err.println("Warning : data contains more attributes "
                             +"than can be displayed as attribute bars.");
          if (m_Log != null) {
            m_Log.logMessage("Warning : data contains more attributes "
                             +"than can be displayed as attribute bars.");
          }
        }
      }
    }

    /**
     * Reset the visualize panel's buttons and the plot panels instances
     * 
     * @param inst	the data
     * @param cIndex	the color index
     */
    private void plotReset(Instances inst, int cIndex) {
      if (m_splitListener == null) {
	m_submit.setText("Reset");
	m_submit.setActionCommand("Reset");
	//if (m_origInstances == null || m_origInstances == inst) {
	if (m_originalPlot == null || m_originalPlot.m_plotInstances == inst) {
	  m_submit.setEnabled(false);
	}
	else {
	  m_submit.setEnabled(true);
	}
      } 
      else {
	m_submit.setEnabled(false);
      }

      m_plotInstances = inst;
      if (m_splitListener != null) {
	m_plotInstances.randomize(new Random());
      }
      m_xIndex=0;
      m_yIndex=0;
      m_cIndex=cIndex;
      cancelShapes();
    }

    /**
     * Set a list of colours to use for plotting points
     * @param cols a list of java.awt.Colors
     */
    public void setColours(FastVector cols) {
      m_plot2D.setColours(cols);
      m_colorList = cols;
    }
    
    /**
     * This will draw the shapes created onto the panel.
     * For best visual, this should be the first thing to be drawn
     * (and it currently is).
     * @param gx The graphics context.
     */
    private void drawShapes(Graphics gx) {
      //FastVector tmp = m_plot.getShapes();
      
      if (m_shapes != null) {
	FastVector stmp;
	int x1, y1, x2, y2;
	for (int noa = 0; noa < m_shapes.size(); noa++) {
	  stmp = (FastVector)m_shapes.elementAt(noa);
	  if (((Double)stmp.elementAt(0)).intValue() == 1) {
	    //then rectangle
	    x1 = (int)m_plot2D.convertToPanelX(((Double)stmp.elementAt(1)).
				      doubleValue());
	    y1 = (int)m_plot2D.convertToPanelY(((Double)stmp.elementAt(2)).
				      doubleValue());
	    x2 = (int)m_plot2D.convertToPanelX(((Double)stmp.elementAt(3)).
				      doubleValue());
	    y2 = (int)m_plot2D.convertToPanelY(((Double)stmp.elementAt(4)).
				      doubleValue());
	    
	    gx.setColor(Color.gray);
	    gx.fillRect(x1, y1, x2 - x1, y2 - y1);
	    gx.setColor(Color.black);
	    gx.drawRect(x1, y1, x2 - x1, y2 - y1);
	    
	  }
	  else if (((Double)stmp.elementAt(0)).intValue() == 2) {
	    //then polygon
	    int[] ar1, ar2;
	    ar1 = getXCoords(stmp);
	    ar2 = getYCoords(stmp);
	    gx.setColor(Color.gray);
	    gx.fillPolygon(ar1, ar2, (stmp.size() - 1) / 2); 
	    gx.setColor(Color.black);
	    gx.drawPolyline(ar1, ar2, (stmp.size() - 1) / 2);
	  }
	  else if (((Double)stmp.elementAt(0)).intValue() == 3) {
	    //then polyline
	    int[] ar1, ar2;
	    FastVector tmp = makePolygon(stmp);
	    ar1 = getXCoords(tmp);
	    ar2 = getYCoords(tmp);
	    
	    gx.setColor(Color.gray);
	    gx.fillPolygon(ar1, ar2, (tmp.size() - 1) / 2);
	    gx.setColor(Color.black);
	    gx.drawPolyline(ar1, ar2, (tmp.size() - 1) / 2);
	  }
	}
      }
      
      if (m_shapePoints != null) {
	//then the current image needs to be refreshed
	if (((Double)m_shapePoints.elementAt(0)).intValue() == 2 ||
	    ((Double)m_shapePoints.elementAt(0)).intValue() == 3) {
	  gx.setColor(Color.black);
	  gx.setXORMode(Color.white);
	  int[] ar1, ar2;
	  ar1 = getXCoords(m_shapePoints);
	  ar2 = getYCoords(m_shapePoints);
	  gx.drawPolyline(ar1, ar2, (m_shapePoints.size() - 1) / 2);
	  m_newMousePos.width = (int)Math.ceil
	    (m_plot2D.convertToPanelX(((Double)m_shapePoints.elementAt
			      (m_shapePoints.size() - 2)).doubleValue()));
	  
	  m_newMousePos.height = (int)Math.ceil
	    (m_plot2D.convertToPanelY(((Double)m_shapePoints.elementAt
			      (m_shapePoints.size() - 1)).doubleValue()));
	  
	  gx.drawLine((int)Math.ceil
		     (m_plot2D.convertToPanelX(((Double)m_shapePoints.elementAt
						(m_shapePoints.size() - 2)).
					       doubleValue())),
		      (int)Math.ceil(m_plot2D.convertToPanelY
				     (((Double)m_shapePoints.elementAt
				       (m_shapePoints.size() - 1)).
				      doubleValue())),
		      m_newMousePos.width, m_newMousePos.height);
	  gx.setPaintMode();
	}
      }
    }
    
    /**
     * This is called for polylines to see where there two lines that
     * extend to infinity cut the border of the view.
     * @param x1 an x point along the line
     * @param y1 the accompanying y point.
     * @param x2 The x coord of the end point of the line.
     * @param y2 The y coord of the end point of the line.
     * @param x 0 or the width of the border line if it has one.
     * @param y 0 or the height of the border line if it has one.
     * @param offset The offset for the border line (either for x or y
     * dependant on which one doesn't change).
     * @return double array that contains the coordinate for the point 
     * that the polyline cuts the border (which ever side that may be).
     */
    private double[] lineIntersect(double x1, double y1, double x2, double y2, 
				   double x, double y, double offset) {
      //the first 4 params are thestart and end points of a line
      //the next param is either 0 for no change in x or change in x, 
      //the next param is the same for y
      //the final 1 is the offset for either x or y (which ever has no change)
      double xval;
      double yval;
      double xn = -100, yn = -100;
      double change;
      if (x == 0) {
	if ((x1 <= offset && offset < x2) || (x1 >= offset && offset > x2)) {
	  //then continue
	  xval = x1 - x2;
	  change = (offset - x2) / xval;
	  yn = (y1 - y2) * change + y2;
	  if (0 <= yn && yn <= y) {
	    //then good
	    xn = offset;
	  }
	  else {
	    //no intersect
	    xn = -100;
	  }
	}
      }
      else if (y == 0) {
	if ((y1 <= offset && offset < y2) || (y1 >= offset && offset > y2)) {
	  //the continue
	  yval = (y1 - y2);
	  change = (offset - y2) / yval;
	  xn = (x1 - x2) * change + x2;
	  if (0 <= xn && xn <= x) {
	    //then good
	    yn = offset;
	  }
	  else {
	    xn = -100;
	  }
	}
      }
      double[] ret = new double[2];
      ret[0] = xn;
      ret[1] = yn;
      return ret;
    }


    /**
     * This will convert a polyline to a polygon for drawing purposes
     * So that I can simply use the polygon drawing function.
     * @param v The polyline to convert.
     * @return A FastVector containing the polygon.
     */
    private FastVector makePolygon(FastVector v) {
      FastVector building = new FastVector(v.size() + 10);
      double x1, y1, x2, y2;
      int edge1 = 0, edge2 = 0;
      for (int noa = 0; noa < v.size() - 2; noa++) {
	building.addElement(new Double(((Double)v.elementAt(noa)).
				       doubleValue()));
      }
      
      //now clip the lines
      double[] new_coords;
      //note lineIntersect , expects the values to have been converted to 
      //screen coords
      //note the first point passed is the one that gets shifted.
      x1 = m_plot2D.convertToPanelX(((Double)v.elementAt(1)).doubleValue());
      y1 = m_plot2D.convertToPanelY(((Double)v.elementAt(2)).doubleValue());
      x2 = m_plot2D.convertToPanelX(((Double)v.elementAt(3)).doubleValue());
      y2 = m_plot2D.convertToPanelY(((Double)v.elementAt(4)).doubleValue());

      if (x1 < 0) {
	//test left
	new_coords = lineIntersect(x1, y1, x2, y2, 0, this.getHeight(), 0);
	edge1 = 0;
	if (new_coords[0] < 0) {
	  //then not left
	  if (y1 < 0) {
	    //test top
	    new_coords = lineIntersect(x1, y1, x2, y2, this.getWidth(), 0, 0);
	    edge1 = 1;
	  }
	  else {
	    //test bottom
	    new_coords = lineIntersect(x1, y1, x2, y2, this.getWidth(), 0, 
				       this.getHeight());
	    edge1 = 3;
	  }
	}
      }
      else if (x1 > this.getWidth()) {
	//test right
	new_coords = lineIntersect(x1, y1, x2, y2, 0, this.getHeight(), 
				   this.getWidth());
	edge1 = 2;
	if (new_coords[0] < 0) {
	  //then not right
	  if (y1 < 0) {
	    //test top
	    new_coords = lineIntersect(x1, y1, x2, y2, this.getWidth(), 0, 0);
	    edge1 = 1;
	  }
	  else {
	    //test bottom
	    new_coords = lineIntersect(x1, y1, x2, y2, this.getWidth(), 0, 
				       this.getHeight());
	    edge1 = 3;
	  }
	}
      }
      else if (y1 < 0) {
	//test top
	new_coords = lineIntersect(x1, y1, x2, y2, this.getWidth(), 0, 0);
	edge1 = 1;
      }
      else {
	//test bottom
	new_coords = lineIntersect(x1, y1, x2, y2, this.getWidth(), 0, 
				   this.getHeight());
	edge1 = 3;
      }
      
      building.setElementAt(new 
	Double(m_plot2D.convertToAttribX(new_coords[0])), 1);
      building.setElementAt(new 
	Double(m_plot2D.convertToAttribY(new_coords[1])), 2);

      x1 = m_plot2D.convertToPanelX(((Double)v.elementAt(v.size() - 4)).
				    doubleValue());
      y1 = m_plot2D.convertToPanelY(((Double)v.elementAt(v.size() - 3)).
				    doubleValue());
      x2 = m_plot2D.convertToPanelX(((Double)v.elementAt(v.size() - 6)).
				    doubleValue());
      y2 = m_plot2D.convertToPanelY(((Double)v.elementAt(v.size() - 5)).
				    doubleValue());
      
      if (x1 < 0) {
	//test left
	new_coords = lineIntersect(x1, y1, x2, y2, 0, this.getHeight(), 0);
	edge2 = 0;
	if (new_coords[0] < 0) {
	  //then not left
	  if (y1 < 0) {
	    //test top
	    new_coords = lineIntersect(x1, y1, x2, y2, this.getWidth(), 0, 0);
	    edge2 = 1;
	  }
	  else {
	    //test bottom
	    new_coords = lineIntersect(x1, y1, x2, y2, this.getWidth(), 0, 
				       this.getHeight());
	    edge2 = 3;
	  }
	}
      }
      else if (x1 > this.getWidth()) {
	//test right
	new_coords = lineIntersect(x1, y1, x2, y2, 0, this.getHeight(), 
				   this.getWidth());
	edge2 = 2;
	if (new_coords[0] < 0) {
	  //then not right
	  if (y1 < 0) {
	    //test top
	    new_coords = lineIntersect(x1, y1, x2, y2, this.getWidth(), 0, 0);
	    edge2 = 1;
	  }
	  else {
	    //test bottom
	    new_coords = lineIntersect(x1, y1, x2, y2, this.getWidth(), 0, 
				       this.getHeight());
	    edge2 = 3;
	  }
	}
      }
      else if (y1 < 0) {
	//test top
	new_coords = lineIntersect(x1, y1, x2, y2, this.getWidth(), 0, 0);
	edge2 = 1;
      }
      else {
	//test bottom
	new_coords = lineIntersect(x1, y1, x2, y2, this.getWidth(), 0, 
				   this.getHeight());
	edge2 = 3;
      }
      
      building.setElementAt(new 
	Double(m_plot2D.convertToAttribX(new_coords[0])), building.size() - 2);
      building.setElementAt(new 
	Double(m_plot2D.convertToAttribY(new_coords[1])), building.size() - 1);
      

      //trust me this complicated piece of code will
      //determine what points on the boundary of the view to add to the polygon
      int xp, yp;

      xp = this.getWidth() * ((edge2 & 1) ^ ((edge2 & 2) / 2));
      yp = this.getHeight() * ((edge2 & 2) / 2);
      //System.out.println(((-1 + 4) % 4) + " hoi");
      
      if (inPolyline(v, m_plot2D.convertToAttribX(xp), 
		     m_plot2D.convertToAttribY(yp))) {
	//then add points in a clockwise direction
	building.addElement(new Double(m_plot2D.convertToAttribX(xp)));
	building.addElement(new Double(m_plot2D.convertToAttribY(yp)));
	for (int noa = (edge2 + 1) % 4; noa != edge1; noa = (noa + 1) % 4) {
	  xp = this.getWidth() * ((noa & 1) ^ ((noa & 2) / 2));
	  yp = this.getHeight() * ((noa & 2) / 2);
	  building.addElement(new Double(m_plot2D.convertToAttribX(xp)));
	  building.addElement(new Double(m_plot2D.convertToAttribY(yp)));
	}
      }
      else {
	xp = this.getWidth() * ((edge2 & 2) / 2);
	yp = this.getHeight() * (1 & ~((edge2 & 1) ^ ((edge2 & 2) / 2)));
	if (inPolyline(v, m_plot2D.convertToAttribX(xp), 
		       m_plot2D.convertToAttribY(yp))) {
	  //then add points in anticlockwise direction
	  building.addElement(new Double(m_plot2D.convertToAttribX(xp)));
	  building.addElement(new Double(m_plot2D.convertToAttribY(yp)));
	  for (int noa = (edge2 + 3) % 4; noa != edge1; noa = (noa + 3) % 4) {
	    xp = this.getWidth() * ((noa & 2) / 2);
	    yp = this.getHeight() * (1 & ~((noa & 1) ^ ((noa & 2) / 2)));
	    building.addElement(new Double(m_plot2D.convertToAttribX(xp)));
	    building.addElement(new Double(m_plot2D.convertToAttribY(yp)));
	  }
	}
      }
      return building;
    }

    /**
     * This will extract from a polygon shape its x coodrdinates
     * so that an awt.Polygon can be created.
     * @param v The polygon shape.
     * @return an int array containing the screen x coords for the polygon.
     */
    private int[] getXCoords(FastVector v) {
      int cach = (v.size() - 1) / 2;
      int[] ar = new int[cach];
      for (int noa = 0; noa < cach; noa ++) {
	ar[noa] = (int)m_plot2D.convertToPanelX(((Double)v.elementAt(noa * 2 +
						1)).doubleValue());
      }
      return ar;
    }

    /**
     * This will extract from a polygon shape its y coordinates
     * so that an awt.Polygon can be created.
     * @param v The polygon shape.
     * @return an int array containing the screen y coords for the polygon.
     */
    private int[] getYCoords(FastVector v) {
      int cach = (v.size() - 1) / 2;
      int[] ar = new int[cach];
      for (int noa = 0; noa < cach; noa ++) {
	ar[noa] = (int)m_plot2D.
	  convertToPanelY(((Double)v.elementAt(noa * 2 + 2)).
			  doubleValue());
      }
      return ar;
    }
    
    /**
     * Renders the polygons if necessary
     * @param gx the graphics context
     */
    public void prePlot(Graphics gx) {
      super.paintComponent(gx);
      if (m_plotInstances != null) {
	drawShapes(gx); // will be in paintComponent of ShapePlot2D
      }
    }
  }



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
  
  /** Lets the user select the attribute for the x axis */
  protected JComboBox m_XCombo = new JComboBox();

  /** Lets the user select the attribute for the y axis */
  protected JComboBox m_YCombo = new JComboBox();

  /** Lets the user select the attribute to use for colouring */
  protected JComboBox m_ColourCombo = new JComboBox();
  
  /** Lets the user select the shape they want to create for instance 
   * selection. */
  protected JComboBox m_ShapeCombo = new JComboBox();

  /** Button for the user to enter the splits. */
  protected JButton m_submit = new JButton("Submit");
  
  /** Button for the user to remove all splits. */
  protected JButton m_cancel = new JButton("Clear");

  /** Button for the user to open the visualized set of instances */
  protected JButton m_openBut = new JButton("Open");

  /** Button for the user to save the visualized set of instances */
  protected JButton m_saveBut = new JButton("Save");

  /** Stop the combos from growing out of control */
  private Dimension COMBO_SIZE = new Dimension(250, m_saveBut
					       .getPreferredSize().height);

  /** file chooser for saving instances */
  protected JFileChooser m_FileChooser 
    = new JFileChooser(new File(System.getProperty("user.dir")));

  /** Filter to ensure only arff files are selected */  
  protected FileFilter m_ArffFilter =
    new ExtensionFileFilter(Instances.FILE_EXTENSION, "Arff data files");

  /** Label for the jitter slider */
  protected JLabel m_JitterLab= new JLabel("Jitter",SwingConstants.RIGHT);

  /** The jitter slider */
  protected JSlider m_Jitter = new JSlider(0,50,0);

  /** The panel that displays the plot */
  protected PlotPanel m_plot = new PlotPanel();

  /** The panel that displays the attributes , using color to represent 
   * another attribute. */
  protected AttributePanel m_attrib = 
    new AttributePanel(m_plot.m_plot2D.getBackground());

  /** The panel that displays legend info if there is more than one plot */
  protected LegendPanel m_legendPanel = new LegendPanel();

  /** Panel that surrounds the plot panel with a titled border */
  protected JPanel m_plotSurround = new JPanel();

  /** Panel that surrounds the class panel with a titled border */
  protected JPanel m_classSurround = new JPanel();

  /** An optional listener that we will inform when ComboBox selections
      change */
  protected ActionListener listener = null;

  /** An optional listener that we will inform when the user creates a 
   * split to seperate instances. */
  protected VisualizePanelListener m_splitListener = null;

  /** The name of the plot (not currently displayed, but can be used
      in the containing Frame or Panel) */
  protected String m_plotName = "";

  /** The panel that displays the legend for the colouring attribute */
  protected ClassPanel m_classPanel = 
    new ClassPanel(m_plot.m_plot2D.getBackground());
  
  /** The list of the colors used */
  protected FastVector m_colorList;

  /** These hold the names of preferred columns to visualize on---if the
      user has defined them in the Visualize.props file */
  protected String m_preferredXDimension = null;
  protected String m_preferredYDimension = null;
  protected String m_preferredColourDimension = null;

  /** Show the attribute bar panel */
  protected boolean m_showAttBars = true;
  
  /** Show the class panel **/
  protected boolean m_showClassPanel = true;

  /** the logger */
  protected Logger m_Log;
  
  /**
   * Sets the Logger to receive informational messages
   *
   * @param newLog the Logger that will now get info messages
   */
  public void setLog(Logger newLog) {
    m_Log = newLog;
  }
  
  /**
   * Set whether the attribute bars should be shown or not.
   * If turned off via this method then any setting in the
   * properties file (if exists) is ignored.
   * 
   * @param sab false if attribute bars are not to be displayed.
   */
  public void setShowAttBars(boolean sab) {
    if (!sab && m_showAttBars) {
      m_plotSurround.remove(m_attrib);
    } else if (sab && !m_showAttBars) {
      GridBagConstraints constraints = new GridBagConstraints();
      constraints.insets = new Insets(0, 0, 0, 0);
      constraints.gridx=4;constraints.gridy=0;constraints.weightx=1;
      constraints.gridwidth=1;constraints.gridheight=1;constraints.weighty=5;
      m_plotSurround.add(m_attrib, constraints);
    }
    m_showAttBars = sab;
    repaint();
  }
  
  /**
   * Gets whether or not attribute bars are being displayed.
   * 
   * @return true if attribute bars are being displayed.
   */
  public boolean getShowAttBars() {
    return m_showAttBars;
  }
  
  /**
   * Set whether the class panel should be shown or not.
   * 
   * @param scp false if class panel is not to be displayed
   */
  public void setShowClassPanel(boolean scp) {
    if (!scp && m_showClassPanel) {
      remove(m_classSurround);
    } else if (scp && !m_showClassPanel) {
      add(m_classSurround, BorderLayout.SOUTH);
    }
    m_showClassPanel = scp;
    repaint();
  }
  
  /**
   * Gets whether or not the class panel is being displayed.
   * 
   * @return true if the class panel is being displayed.
   */
  public boolean getShowClassPanel() {
    return m_showClassPanel;
  }

  /** 
   * This constructor allows a VisualizePanelListener to be set. 
   * 
   * @param ls		the listener to use
   */
  public VisualizePanel(VisualizePanelListener ls) {
    this();
    m_splitListener = ls;
  }

  /**
   * Set the properties for the VisualizePanel
   * 
   * @param relationName	the name of the relation, can be null
   */
  private void setProperties(String relationName) {
    if (VisualizeUtils.VISUALIZE_PROPERTIES != null) {
      String thisClass = this.getClass().getName();
      if (relationName == null) {
	
	String showAttBars = thisClass+".displayAttributeBars";

	String val = VisualizeUtils.VISUALIZE_PROPERTIES.
	  getProperty(showAttBars);
	if (val == null) {
	  //System.err.println("Displaying attribute bars ");
//	  m_showAttBars = true;
	} else {
	  // only check if this hasn't been turned off programatically 
	  if (m_showAttBars) {
	    if (val.compareTo("true") == 0 || val.compareTo("on") == 0) {
	      //System.err.println("Displaying attribute bars ");
	      m_showAttBars = true;
	    } else {
	      m_showAttBars = false;
	    }
	  }
	}
      } else {
	/*
	System.err.println("Looking for preferred visualization dimensions for "
			   +relationName);
	*/
	String xcolKey = thisClass+"."+relationName+".XDimension";
	String ycolKey = thisClass+"."+relationName+".YDimension";
	String ccolKey = thisClass+"."+relationName+".ColourDimension";
      
	m_preferredXDimension = VisualizeUtils.VISUALIZE_PROPERTIES.
	  getProperty(xcolKey);
	/*
	if (m_preferredXDimension == null) {
	  System.err.println("No preferred X dimension found in "
			     +VisualizeUtils.PROPERTY_FILE
			     +" for "+xcolKey);
	} else {
	  System.err.println("Setting preferred X dimension to "
			     +m_preferredXDimension);
			     }*/
	m_preferredYDimension = VisualizeUtils.VISUALIZE_PROPERTIES.
	  getProperty(ycolKey);
	/*
	if (m_preferredYDimension == null) {
	  System.err.println("No preferred Y dimension found in "
			     +VisualizeUtils.PROPERTY_FILE
			     +" for "+ycolKey);
	} else {
	  System.err.println("Setting preferred dimension Y to "
			     +m_preferredYDimension);
			     }*/
	m_preferredColourDimension = VisualizeUtils.VISUALIZE_PROPERTIES.
	  getProperty(ccolKey);
	/*
	if (m_preferredColourDimension == null) {
	  System.err.println("No preferred Colour dimension found in "
			     +VisualizeUtils.PROPERTY_FILE
			     +" for "+ycolKey);
	} else {
	  System.err.println("Setting preferred Colour dimension to "
			     +m_preferredColourDimension);
			     }*/
      }
    }
  }

  /**
   * Constructor
   */
  public VisualizePanel() {
    super();
    setProperties(null);
    m_FileChooser.setFileFilter(m_ArffFilter);
    m_FileChooser.setFileSelectionMode(JFileChooser.FILES_ONLY);

    m_XCombo.setToolTipText("Select the attribute for the x axis");
    m_YCombo.setToolTipText("Select the attribute for the y axis");
    m_ColourCombo.setToolTipText("Select the attribute to colour on");
    m_ShapeCombo.setToolTipText("Select the shape to use for data selection"); 

    m_XCombo.setPreferredSize(COMBO_SIZE);
    m_YCombo.setPreferredSize(COMBO_SIZE);
    m_ColourCombo.setPreferredSize(COMBO_SIZE);
    m_ShapeCombo.setPreferredSize(COMBO_SIZE);

    m_XCombo.setMaximumSize(COMBO_SIZE);
    m_YCombo.setMaximumSize(COMBO_SIZE);
    m_ColourCombo.setMaximumSize(COMBO_SIZE);
    m_ShapeCombo.setMaximumSize(COMBO_SIZE);
    
    m_XCombo.setMinimumSize(COMBO_SIZE);
    m_YCombo.setMinimumSize(COMBO_SIZE);
    m_ColourCombo.setMinimumSize(COMBO_SIZE);
    m_ShapeCombo.setMinimumSize(COMBO_SIZE);
    //////////
    m_XCombo.setEnabled(false);
    m_YCombo.setEnabled(false);
    m_ColourCombo.setEnabled(false);
    m_ShapeCombo.setEnabled(false);

    // tell the class panel and the legend panel that we want to know when 
    // colours change
    m_classPanel.addRepaintNotify(this);
    m_legendPanel.addRepaintNotify(this);
    
    // Check the default colours against the background colour of the
    // plot panel. If any are equal to the background colour then
    // change them (so they are visible :-)
    for (int i = 0; i < m_DefaultColors.length; i++) {
      Color c = m_DefaultColors[i];
      if (c.equals(m_plot.m_plot2D.getBackground())) {
        int red = c.getRed();
        int blue = c.getBlue();
        int green = c.getGreen();
        red += (red < 128) ? (255 - red) / 2 : -(red / 2);
        blue += (blue < 128) ? (blue - red) / 2 : -(blue / 2);
        green += (green< 128) ? (255 - green) / 2 : -(green / 2);
        m_DefaultColors[i] = new Color(red, green, blue);
      }
    }
    m_classPanel.setDefaultColourList(m_DefaultColors);
    m_attrib.setDefaultColourList(m_DefaultColors);

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
    m_plot.setColours(m_colorList);
    m_classPanel.setColours(m_colorList);
    m_attrib.setColours(m_colorList);
    m_attrib.addAttributePanelListener(new AttributePanelListener() {
	public void attributeSelectionChange(AttributePanelEvent e) {
	  if (e.m_xChange) {
	    m_XCombo.setSelectedIndex(e.m_indexVal);
	  } else {
	    m_YCombo.setSelectedIndex(e.m_indexVal);
	  }
	}
      });
    
    m_XCombo.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e) {
	  int selected = m_XCombo.getSelectedIndex();
	  if (selected < 0) {
	    selected = 0;
	  }
	  m_plot.setXindex(selected);
	 
	  // try sending on the event if anyone is listening
	  if (listener != null) {
	    listener.actionPerformed(e);
	  }
	}
      });

    m_YCombo.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e) {
	  int selected = m_YCombo.getSelectedIndex();
	  if (selected < 0) {
	    selected = 0;
	  }
	  m_plot.setYindex(selected);
	 
	  // try sending on the event if anyone is listening
	  if (listener != null) {
	    listener.actionPerformed(e);
	  }
	}
      });

    m_ColourCombo.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e) {
	  int selected = m_ColourCombo.getSelectedIndex();
	  if (selected < 0) {
	    selected = 0;
	  }
	  m_plot.setCindex(selected);

	  if (listener != null) {
	    listener.actionPerformed(e);
	  }
	}
      });
    
    ///////
    m_ShapeCombo.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e) {
	  int selected = m_ShapeCombo.getSelectedIndex();
	  if (selected < 0) {
	    selected = 0;
	  }
	  m_plot.setSindex(selected);
	  // try sending on the event if anyone is listening
	  if (listener != null) {
	    listener.actionPerformed(e);
	  }
	}
      });


    ///////////////////////////////////////

    m_Jitter.addChangeListener(new ChangeListener() {
	public void stateChanged(ChangeEvent e) {
	  m_plot.setJitter(m_Jitter.getValue());
	}
      });

    m_openBut.setToolTipText("Loads previously saved instances from a file");
    m_openBut.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e) {
	  openVisibleInstances();
	}
      });
    
    m_saveBut.setEnabled(false);
    m_saveBut.setToolTipText("Save the visible instances to a file");
    m_saveBut.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e) {
	  saveVisibleInstances();
	}
      });
    
    JPanel combos = new JPanel();
    GridBagLayout gb = new GridBagLayout();
    GridBagConstraints constraints = new GridBagConstraints();


    m_XCombo.setLightWeightPopupEnabled(false);
    m_YCombo.setLightWeightPopupEnabled(false);
    m_ColourCombo.setLightWeightPopupEnabled(false);
    m_ShapeCombo.setLightWeightPopupEnabled(false);
    combos.setBorder(BorderFactory.createEmptyBorder(10, 5, 10, 5));

    combos.setLayout(gb);
    constraints.gridx=0;constraints.gridy=0;constraints.weightx=5;
    constraints.fill = GridBagConstraints.HORIZONTAL;
    constraints.gridwidth=2;constraints.gridheight=1;
    constraints.insets = new Insets(0,2,0,2);
    combos.add(m_XCombo,constraints);
    constraints.gridx=2;constraints.gridy=0;constraints.weightx=5;
    constraints.gridwidth=2;constraints.gridheight=1;
    combos.add(m_YCombo,constraints);
    constraints.gridx=0;constraints.gridy=1;constraints.weightx=5;
    constraints.gridwidth=2;constraints.gridheight=1;
    combos.add(m_ColourCombo,constraints);
    //
    constraints.gridx=2;constraints.gridy=1;constraints.weightx=5;
    constraints.gridwidth=2;constraints.gridheight=1;
    combos.add(m_ShapeCombo,constraints);

   
    JPanel mbts = new JPanel();
    mbts.setLayout(new GridLayout(1,4));
    mbts.add(m_submit); mbts.add(m_cancel); mbts.add(m_openBut); mbts.add(m_saveBut);

    constraints.gridx=0;constraints.gridy=2;constraints.weightx=5;
    constraints.gridwidth=2;constraints.gridheight=1;
    combos.add(mbts, constraints);

    ////////////////////////////////
    constraints.gridx=2;constraints.gridy=2;constraints.weightx=5;
    constraints.gridwidth=1;constraints.gridheight=1;
    constraints.insets = new Insets(10,0,0,5);
    combos.add(m_JitterLab,constraints);
    constraints.gridx=3;constraints.gridy=2;
    constraints.weightx=5;
    constraints.insets = new Insets(10,0,0,0);
    combos.add(m_Jitter,constraints);

    m_classSurround = new JPanel();
    m_classSurround.
      setBorder(BorderFactory.createTitledBorder("Class colour")); 
    m_classSurround.setLayout(new BorderLayout());

    m_classPanel.setBorder(BorderFactory.createEmptyBorder(15,10,10,10));
    m_classSurround.add(m_classPanel, BorderLayout.CENTER);

    GridBagLayout gb2 = new GridBagLayout();
    m_plotSurround.setBorder(BorderFactory.createTitledBorder("Plot"));
    m_plotSurround.setLayout(gb2);

    constraints.fill = GridBagConstraints.BOTH;
    constraints.insets = new Insets(0, 0, 0, 10);
    constraints.gridx=0;constraints.gridy=0;constraints.weightx=3;
    constraints.gridwidth=4;constraints.gridheight=1;constraints.weighty=5;
    m_plotSurround.add(m_plot, constraints);
    
    if (m_showAttBars) {
      constraints.insets = new Insets(0, 0, 0, 0);
      constraints.gridx=4;constraints.gridy=0;constraints.weightx=1;
      constraints.gridwidth=1;constraints.gridheight=1;constraints.weighty=5;
      m_plotSurround.add(m_attrib, constraints);
    }

    setLayout(new BorderLayout());
    add(combos, BorderLayout.NORTH);
    add(m_plotSurround, BorderLayout.CENTER);
    add(m_classSurround, BorderLayout.SOUTH);
    
    String [] SNames = new String [4];
    SNames[0] = "Select Instance";
    SNames[1] = "Rectangle";
    SNames[2] = "Polygon";
    SNames[3] = "Polyline";

    m_ShapeCombo.setModel(new DefaultComboBoxModel(SNames));
    m_ShapeCombo.setEnabled(true);
  }

  /**
   * displays the previously saved instances
   * 
   * @param insts	the instances to display
   * @throws Exception	if display is not possible
   */
  protected void openVisibleInstances(Instances insts) throws Exception {
    PlotData2D tempd = new PlotData2D(insts);
    tempd.setPlotName(insts.relationName());
    tempd.addInstanceNumberAttribute();
    m_plot.m_plot2D.removeAllPlots();
    addPlot(tempd);
    
    // modify title
    Component parent = getParent();
    while (parent != null) {
      if (parent instanceof JFrame) {
	((JFrame) parent).setTitle(
	    "Weka Classifier Visualize: " 
	    + insts.relationName() 
	    + " (display only)");
	break;
      }
      else {
	parent = parent.getParent();
      }
    }
  }
  
  /**
   * Loads previously saved instances from a file
   */
  protected void openVisibleInstances() {
    try {
      int returnVal = m_FileChooser.showOpenDialog(this);
      if (returnVal == JFileChooser.APPROVE_OPTION) {
	File sFile = m_FileChooser.getSelectedFile();
	if (!sFile.getName().toLowerCase().
	    endsWith(Instances.FILE_EXTENSION)) {
	  sFile = new File(sFile.getParent(), sFile.getName() + Instances.FILE_EXTENSION);
	}
	File selected = sFile;
	Instances insts = new Instances(new BufferedReader(new FileReader(selected)));
	openVisibleInstances(insts);
      }
    } catch (Exception ex) {
      ex.printStackTrace();
      m_plot.m_plot2D.removeAllPlots();
      JOptionPane.showMessageDialog(
	  this, 
	  ex.getMessage(), 
	  "Error loading file...", 
	  JOptionPane.ERROR_MESSAGE);
    }
  }

  /**
   * Save the currently visible set of instances to a file
   */
  private void saveVisibleInstances() {
    FastVector plots = m_plot.m_plot2D.getPlots();
    if (plots != null) {
      PlotData2D master = (PlotData2D)plots.elementAt(0);
      Instances saveInsts = new Instances(master.getPlotInstances());
      for (int i = 1; i < plots.size(); i++) {
	PlotData2D temp = (PlotData2D)plots.elementAt(i);
	Instances addInsts = temp.getPlotInstances();
	for (int j = 0; j < addInsts.numInstances(); j++) {
	  saveInsts.add(addInsts.instance(j));
	}
      }
      try {
	int returnVal = m_FileChooser.showSaveDialog(this);
	if (returnVal == JFileChooser.APPROVE_OPTION) {
	  File sFile = m_FileChooser.getSelectedFile();
	  if (!sFile.getName().toLowerCase().
	      endsWith(Instances.FILE_EXTENSION)) {
	    sFile = new File(sFile.getParent(), sFile.getName() 
			     + Instances.FILE_EXTENSION);
	  }
	  File selected = sFile;
	  Writer w = new BufferedWriter(new FileWriter(selected));
	  w.write(saveInsts.toString());
	  w.close();
	}
      } catch (Exception ex) {
	ex.printStackTrace();
      }
    }
  }


  /**
   * Sets the index used for colouring. If this method is called then
   * the supplied index is used and the combo box for selecting colouring
   * attribute is disabled.
   * @param index the index of the attribute to use for colouring
   */
  public void setColourIndex(int index) {
    if (index >= 0) {
      m_ColourCombo.setSelectedIndex(index);
    } else {
      m_ColourCombo.setSelectedIndex(0);
    }
    m_ColourCombo.setEnabled(false);
  }
  
  
  /**
   * Set the index of the attribute for the x axis 
   * @param index the index for the x axis
   * @exception Exception if index is out of range.
   */
  public void setXIndex(int index) throws Exception {
    if (index >= 0 && index < m_XCombo.getItemCount()) {
      m_XCombo.setSelectedIndex(index);
    } else {
      throw new Exception("x index is out of range!");
    }
  }

  /**
   * Get the index of the attribute on the x axis
   * @return the index of the attribute on the x axis
   */
  public int getXIndex() {
    return m_XCombo.getSelectedIndex();
  }

  /**
   * Set the index of the attribute for the y axis 
   * @param index the index for the y axis
   * @exception Exception if index is out of range.
   */
  public void setYIndex(int index) throws Exception {
    if (index >= 0 && index < m_YCombo.getItemCount()) {
      m_YCombo.setSelectedIndex(index);
    } else {
      throw new Exception("y index is out of range!");
    }
  }
  
  /**
   * Get the index of the attribute on the y axis
   * @return the index of the attribute on the x axis
   */
  public int getYIndex() {
    return m_YCombo.getSelectedIndex();
  }

  /**
   * Get the index of the attribute selected for coloring
   * @return the index of the attribute on the x axis
   */
  public int getCIndex() {
    return m_ColourCombo.getSelectedIndex();
  }

  /**
   * Get the index of the shape selected for creating splits.
   * @return The index of the shape.
   */
  public int getSIndex() {
    return m_ShapeCombo.getSelectedIndex();
  }
  
  /** 
   * Set the shape for creating splits.
   * @param index The index of the shape.
   * @exception Exception if index is out of range.
   */
  public void setSIndex(int index) throws Exception {
    if (index >= 0 && index < m_ShapeCombo.getItemCount()) {
      m_ShapeCombo.setSelectedIndex(index);
    }
    else {
      throw new Exception("s index is out of range!");
    }
  }

  /**
   * Add a listener for this visualize panel
   * @param act an ActionListener
   */
  public void addActionListener(ActionListener act) {
    listener = act;
  }

  /**
   * Set a name for this plot
   * @param plotName the name for the plot
   */
  public void setName(String plotName) {
    m_plotName = plotName;
  }

  /**
   * Returns the name associated with this plot. "" is returned if no
   * name is set.
   * @return the name of the plot
   */
  public String getName() {
    return m_plotName;
  }

  /**
   * Get the master plot's instances
   * @return the master plot's instances
   */
  public Instances getInstances() {
    return m_plot.m_plotInstances;
  }

  /**
   * Sets the Colors in use for a different attrib
   * if it is not a nominal attrib and or does not have
   * more possible values then this will do nothing.
   * otherwise it will add default colors to see that
   * there is a color for the attrib to begin with.
   * @param a The index of the attribute to color.
   * @param i The instances object that contains the attribute.
   */
  protected void newColorAttribute(int a, Instances i) {
    if (i.attribute(a).isNominal() || i.attribute(a).isRanking()) {
      for (int noa = m_colorList.size(); noa < i.attribute(a).numValues();
	   noa++) {
	Color pc = m_DefaultColors[noa % 10];
	int ija =  noa / 10;
	ija *= 2; 
	for (int j=0;j<ija;j++) {
	  pc = pc.brighter();
	}
	
	m_colorList.addElement(pc);
      }
      m_plot.setColours(m_colorList);
      m_attrib.setColours(m_colorList);
      m_classPanel.setColours(m_colorList);
    }
  }


  /**
   * This will set the shapes for the instances.
   * @param l A list of the shapes, providing that
   * the objects in the lists are non editable the data will be
   * kept intact.
   */
  public void setShapes(FastVector l) {
    m_plot.setShapes(l);
  }

  /**
   * Tells the panel to use a new set of instances.
   * @param inst a set of Instances
   */
  public void setInstances(Instances inst) {
    if (inst.numAttributes() > 0 && inst.numInstances() > 0) {
      newColorAttribute(inst.numAttributes()-1, inst);
    }

    PlotData2D temp = new PlotData2D(inst);
    temp.setPlotName(inst.relationName());
    
    try {
      setMasterPlot(temp);
    } catch (Exception ex) {
      System.err.println(ex);
      ex.printStackTrace();
    } 
  }

  /**
   * initializes the comboboxes based on the data
   * 
   * @param inst	the data to base the combobox-setup on
   */
  public void setUpComboBoxes(Instances inst) {
    setProperties(inst.relationName());
    int prefX = -1;
    int prefY = -1;
    if (inst.numAttributes() > 1) {
      prefY = 1;
    }
    int prefC = -1;
    String [] XNames = new String [inst.numAttributes()];
    String [] YNames = new String [inst.numAttributes()];
    String [] CNames = new String [inst.numAttributes()];
    for (int i = 0; i < XNames.length; i++) {
      String type = "";
      switch (inst.attribute(i).type()) {
      case Attribute.NOMINAL:
	type = " (Nom)";
	break;
      case Attribute.NUMERIC:
	type = " (Num)";
	break;
      case Attribute.STRING:
	type = " (Str)";
	break;
      case Attribute.DATE:
	type = " (Dat)";
	break;
      case Attribute.RELATIONAL:
	type = " (Rel)";
	break;
      case PreferenceAttribute.RANKING:
    type = " (Rnk)";
    break;
      default:
	type = " (???)";
      }
      XNames[i] = "X: "+ inst.attribute(i).name()+type;
      YNames[i] = "Y: "+ inst.attribute(i).name()+type;
      CNames[i] = "Colour: "+ inst.attribute(i).name()+type;
      if (m_preferredXDimension != null) {
	if (m_preferredXDimension.compareTo(inst.attribute(i).name()) == 0) {
	  prefX = i;
	  //System.err.println("Found preferred X dimension");
	}
      }
      if (m_preferredYDimension != null) {
	if (m_preferredYDimension.compareTo(inst.attribute(i).name()) == 0) {
	  prefY = i;
	  //System.err.println("Found preferred Y dimension");
	}
      }
      if (m_preferredColourDimension != null) {
	if (m_preferredColourDimension.
	    compareTo(inst.attribute(i).name()) == 0) {
	  prefC = i;
	  //System.err.println("Found preferred Colour dimension");
	}
      }
    }
    m_XCombo.setModel(new DefaultComboBoxModel(XNames));
    m_YCombo.setModel(new DefaultComboBoxModel(YNames));

    m_ColourCombo.setModel(new DefaultComboBoxModel(CNames));
    //m_ShapeCombo.setModel(new DefaultComboBoxModel(SNames));
    //m_ShapeCombo.setEnabled(true);
    m_XCombo.setEnabled(true);
    m_YCombo.setEnabled(true);
    
    if (m_splitListener == null) {
      m_ColourCombo.setEnabled(true);
      m_ColourCombo.setSelectedIndex(inst.numAttributes()-1);
    }
    m_plotSurround.setBorder((BorderFactory.createTitledBorder("Plot: "
			      +inst.relationName())));
    try {
      if (prefX != -1) {
	setXIndex(prefX);
      }
      if (prefY != -1) {
	setYIndex(prefY);
      }
      if (prefC != -1) {
	m_ColourCombo.setSelectedIndex(prefC);
      }
    } catch (Exception ex) {
      System.err.println("Problem setting preferred Visualization dimensions");
    }
  }

  /**
   * Removes all the plots.
   */
  public void removeAllPlots() {
    m_plot.removeAllPlots();
  }
  
  /**
   * Set the master plot for the visualize panel
   * @param newPlot the new master plot
   * @exception Exception if the master plot could not be set
   */
  public void setMasterPlot(PlotData2D newPlot) throws Exception {
    m_plot.setMasterPlot(newPlot);
    setUpComboBoxes(newPlot.m_plotInstances);
    m_saveBut.setEnabled(true);
    repaint();
  }

  /**
   * Set a new plot to the visualize panel
   * @param newPlot the new plot to add
   * @exception Exception if the plot could not be added
   */
  public void addPlot(PlotData2D newPlot) throws Exception {
    m_plot.addPlot(newPlot);
    if (m_plot.m_plot2D.getMasterPlot() != null) {
      setUpComboBoxes(newPlot.m_plotInstances);
    }
    m_saveBut.setEnabled(true);
    repaint();
  }
  
  /**
   * Returns the underlying plot panel.
   * 
   * @return		the plot panel
   */
  public PlotPanel getPlotPanel() {
    return m_plot;
  }

  /**
   * Main method for testing this class
   * 
   * @param args	the commandline parameters
   */
  public static void main(String [] args) {
    try {
      if (args.length < 1) {
	System.err.println("Usage : weka.gui.visualize.VisualizePanel "
			   +"<dataset> [<dataset> <dataset>...]");
	System.exit(1);
      }

      weka.core.logging.Logger.log(weka.core.logging.Logger.Level.INFO, "Logging started");
      final javax.swing.JFrame jf = 
	new javax.swing.JFrame("Weka Explorer: Visualize");
      jf.setSize(500,400);
      jf.getContentPane().setLayout(new BorderLayout());
      final VisualizePanel sp = new VisualizePanel();
      
      jf.getContentPane().add(sp, BorderLayout.CENTER);
      jf.addWindowListener(new java.awt.event.WindowAdapter() {
	public void windowClosing(java.awt.event.WindowEvent e) {
	  jf.dispose();
	  System.exit(0);
	}
      });

      jf.setVisible(true);
      if (args.length >= 1) {
	for (int j = 0; j < args.length; j++) {
	  System.err.println("Loading instances from " + args[j]);
	  java.io.Reader r = new java.io.BufferedReader(
			     new java.io.FileReader(args[j]));
	  Instances i = new Instances(r);
	  i.setClassIndex(i.numAttributes()-1);
	  PlotData2D pd1 = new PlotData2D(i);
	  
	  if (j == 0) {
	    pd1.setPlotName("Master plot");
	    sp.setMasterPlot(pd1);
	  } else {
	    pd1.setPlotName("Plot "+(j+1));
	    pd1.m_useCustomColour = true;
	    pd1.m_customColour = (j % 2 == 0) ? Color.red : Color.blue; 
	    sp.addPlot(pd1);
	  }
	}
      }
    } catch (Exception ex) {
      ex.printStackTrace();
      System.err.println(ex.getMessage());
    }
  }
}
