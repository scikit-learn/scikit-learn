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
 *    AttributeVisualizationPanel.java
 *    Copyright (C) 2003 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui;

import weka.core.Attribute;
import weka.core.AttributeStats;
import weka.core.FastVector;
import weka.core.Instances;
import weka.core.SparseInstance;
import weka.core.Utils;
import weka.core.labelranking.PreferenceAttribute;
import weka.core.labelranking.PreferenceDenseInstance;
import weka.core.labelranking.RankUtilities;
import weka.gui.visualize.PrintableComponent;
import weka.gui.visualize.PrintablePanel;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.MouseEvent;
import java.io.FileReader;

import javax.swing.JComboBox;
import javax.swing.JFrame;

/**
 * Creates a panel that shows a visualization of an
 * attribute in a dataset. For nominal attribute it
 * shows a bar plot, with each bar corresponding to
 * each nominal value of the attribute with its height
 * equal to the frequecy that value appears in the
 * dataset. For numeric attributes, it displays a
 * histogram. The width of an interval in the
 * histogram is calculated using Scott's(1979)
 * method: <br>
 *    intervalWidth = Max(1, 3.49*Std.Dev*numInstances^(1/3))
 * Then the number of intervals is calculated by: <br>
 *   intervals = max(1, Math.round(Range/intervalWidth);
 *
 * @author Ashraf M. Kibriya (amk14@cs.waikato.ac.nz)
 * @version $Revision: 5721 $
 */
public class AttributeVisualizationPanel
  extends PrintablePanel {

  /** for serialization */
  private static final long serialVersionUID = -8650490488825371193L;
  
  /** This holds the current set of instances */
  protected Instances m_data;
  
  //RANKING BEGIN
  protected int[][] backupPairRankWeights;
  //RANKING END
  
  /**
   * This holds the attribute stats of the current attribute on display. It is
   * calculated in setAttribute(int idx) when it is called to set a new
   * attribute index.
   */
  protected AttributeStats m_as;
  
  
  /** This holds the index of the current attribute on display and should be
   *  set through setAttribute(int idx).
   */
  protected int m_attribIndex;
  
  /**
   * This holds the max value of the current attribute. In case of nominal
   * attribute it is the highest count that a nominal value has in the
   * attribute (given by m_as.nominalWeights[i]), otherwise in case of numeric
   * attribute it is simply the maximum value present in the attribute (given by
   * m_as.numericStats.max). It is used to calculate the ratio of the height of
   * the bars with respect to the height of the display area.
   */
  protected double m_maxValue;
  
  /**
   * This array holds the count (or height) for the each of the bars in a
   * barplot or a histogram. In case of barplots (and current attribute being
   * nominal) its length (and the number of bars) is equal to the number of
   * nominal values in the current attribute, with each field of the array being
   * equal to the count of each nominal that it represents (the count of ith
   * nominal value of an attribute is given by m_as.nominalWeights[i]). Whereas,
   * in case of histograms (and current attribute being numeric) the width of
   * its intervals is calculated by Scott's(1979) method: <br>
   *    intervalWidth = Max(1, 3.49*Std.Dev*numInstances^(1/3))
   * And the number of intervals by: <br>
   *   intervals = max(1, Math.round(Range/intervalWidth);
   * Then each field of this array contains the number of values of the current
   * attribute that fall in the histogram interval that it represents. <br>
   * NOTE: The values of this array are only calculated if the class attribute
   * is not set or if it is numeric.
   */
  protected double[] m_histBarCounts;
  
  /**
   * This array holds the per class count (or per class height) of the each of
   * the bars in a barplot or a histogram.
   * For nominal attributes the format is: <br>
   *    m_histBarClassCounts[nominalValue][classValue+1].
   * For numeric attributes the format is: <br>
   *    m_histBarClassCounts[interval][classValues+1], <br>
   *      where the number of intervals is calculated by the Scott's method as
   *            mentioned above.
   * The array is initialized to have 1+numClasses to accomodate for instances
   * with missing class value. The ones with missing class value are displayed
   * as a black sub par in a histogram or a barplot.
   *
   * NOTE: The values of this array are only calculated if the class attribute
   * is set and it is nominal.
   */
  SparseInstance m_histBarClassCounts[];
  
  /**
   * Contains the range of each bar in a histogram. It is used to work out the
   * range of bar the mouse pointer is on in getToolTipText().
   */
  protected double m_barRange;
  
  /** Contains the current class index. */
  protected int m_classIndex;
  
  /** This stores the BarCalc or HistCalc thread while a new barplot or
   * histogram is being calculated. */
  private Thread m_hc;
  
  /** True if the thread m_hc above is running. */
  private boolean m_threadRun=false;
  
  private boolean m_doneCurrentAttribute = false;
  private boolean m_displayCurrentAttribute = false;
  
  /** This stores and lets the user select a class attribute. It also has
   * an entry "No Class" if the user does not want to set a class attribute
   * for colouring.
   */
  protected JComboBox m_colorAttrib;
  
  /**
   * Fontmetrics used to get the font size which is required for calculating
   * displayable area size, bar height ratio and width of strings that are
   * displayed on top of bars indicating their count.
   */
  private FontMetrics m_fm;
  
  /**
   * Lock variable to synchronize the different threads running currently in
   * this class. There are two to three threads in this class, AWT paint thread
   * which is handled differently in paintComponent() which checks on
   * m_threadRun to determine if it can perform full paint or not, the second
   * thread is the main execution thread and the third is the one represented by
   * m_hc which we start when we want to calculate the internal fields for a bar
   * plot or a histogram.
   */
  private Integer m_locker = new Integer(1);
  
  //Image img;
  
  /** Contains discrete colours for colouring of subbars of histograms and
   * bar plots when the class attribute is set and is nominal
   */
  private FastVector m_colorList = new FastVector();
  
  /** default colour list */
  private static final Color [] m_defaultColors = {Color.blue,
  Color.red,
  Color.cyan,
  new Color(75, 123, 130),
  Color.pink,
  Color.green,
  Color.orange,
  new Color(255, 0, 255),
  new Color(255, 0, 0),
  new Color(0, 255, 0),
  };
   
  
  /**
   * Constructor - If used then the class will not show the class selection
   * combo box.
   */
  public AttributeVisualizationPanel() {
    this(false);
  }
  
  /**
   * Constructor.
   * @param showColouringOption - should be true if the class selection combo
   * box is to be displayed with the histogram/barplot, or false otherwise.
   * P.S: the combo box is always created it just won't be shown if
   * showColouringOption is false.
   */
  public AttributeVisualizationPanel(boolean showColouringOption) {
    this.setFont( new Font("Default", Font.PLAIN, 9) );
    m_fm = this.getFontMetrics( this.getFont() );
    this.setToolTipText("");
    FlowLayout fl= new FlowLayout(FlowLayout.LEFT);
    this.setLayout(fl);
    this.addComponentListener( new ComponentAdapter() {
      public void componentResized(ComponentEvent ce) {
        if(m_data!=null) {
//          calcGraph();
        }
      }
    });
    
    m_colorAttrib = new JComboBox();
    m_colorAttrib.addItemListener( new ItemListener() {
      public void itemStateChanged(ItemEvent ie) {
        if(ie.getStateChange()==ItemEvent.SELECTED) {
          m_classIndex = m_colorAttrib.getSelectedIndex() - 1;
          if (m_as != null) {
            setAttribute(m_attribIndex);
          }
        }
      }
    });
    
    if(showColouringOption) {
      //m_colorAttrib.setVisible(false);
      this.add(m_colorAttrib);
      validate();
    }
  }
  
  /**
   * Sets the instances for use
   *
   * @param newins a set of Instances
   */
  public void setInstances(Instances newins) {
    m_attribIndex = 0;
    m_as = null;
    m_data = new Instances(newins);
    if(m_colorAttrib!=null) {
      m_colorAttrib.removeAllItems();
      m_colorAttrib.addItem("No class");
      for(int i=0; i<m_data.numAttributes(); i++) {
	String type = "";
	switch (m_data.attribute(i).type()) {
	
	  case PreferenceAttribute.RANKING: 
		type =" (Rnk) ";
		break;
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
	  default:
	    type = "(???) ";
	}
        m_colorAttrib.addItem(new String("Class: "+m_data.attribute(i).name()+
        " " + type));
      }
      if (m_data.classIndex() >= 0) {
        m_colorAttrib.setSelectedIndex(m_data.classIndex() + 1);
      } else {
        m_colorAttrib.setSelectedIndex(m_data.numAttributes());
      }
      //if (m_data.classIndex() >= 0) {
      //    m_colorAttrib.setSelectedIndex(m_data.classIndex());
      //}
    }
    if (m_data.classIndex() >= 0) {
      m_classIndex = m_data.classIndex();
    } else {
      m_classIndex = m_data.numAttributes()-1;
    }
    
  }
  
  /**
   * Returns the class selection combo box if the parent component wants to
   * place it in itself or in some component other than this component.
   */
  public JComboBox getColorBox() {
    return m_colorAttrib;
  }
  
  /**
   * Get the coloring (class) index for the plot
   *
   * @return an <code>int</code> value
   */
  public int getColoringIndex() {
    return m_classIndex; //m_colorAttrib.getSelectedIndex();
  }
  
  /**
   * Set the coloring (class) index for the plot
   *
   * @param ci an <code>int</code> value
   */
  public void setColoringIndex(int ci) {
    m_classIndex = ci;
    if(m_colorAttrib!=null)
      m_colorAttrib.setSelectedIndex(ci + 1);
    else
      setAttribute(m_attribIndex);
  }
  
  /**
   * Tells the panel which attribute to visualize.
   *
   * @param index The index of the attribute
   */
  public void setAttribute(int index) {
    
    synchronized (m_locker) {
      //m_threadRun = true;
      m_threadRun = false;
      m_doneCurrentAttribute = false;
      m_displayCurrentAttribute = true;
      //if(m_hc!=null && m_hc.isAlive()) m_hc.stop();
      m_attribIndex = index;
      m_as = m_data.attributeStats(m_attribIndex);
      //m_classIndex = m_colorAttrib.getSelectedIndex();
    }
    this.repaint();
    // calcGraph();
  }
  
  /**
   * Recalculates the barplot or histogram to display, required usually when the
   * attribute is changed or the component is resized.
   */
  public void calcGraph(int panelWidth, int panelHeight) {
    
    synchronized (m_locker) {
      m_threadRun = true;
      if(m_as.nominalWeights!=null) {
        m_hc = new BarCalc(panelWidth, panelHeight);
        m_hc.setPriority(m_hc.MIN_PRIORITY);
        m_hc.start();
      }
      else if(m_as.numericStats!=null) {
        m_hc = new HistCalc();
        m_hc.setPriority(m_hc.MIN_PRIORITY);
        m_hc.start();
      } else {
        m_histBarCounts = null;
        m_histBarClassCounts = null;
        m_doneCurrentAttribute = true;
        m_threadRun = false;
        this.repaint();
      }
    }
  }
  
  /**
   * Internal class that calculates the barplot to display, in a separate
   * thread. In particular it initializes some of the crucial internal fields
   * required by paintComponent() to display the histogram for the current
   * attribute. These include: m_histBarCounts or m_histBarClassCounts,
   * m_maxValue and m_colorList.
   */
  private class BarCalc extends Thread {
    private int m_panelWidth;
    private int m_panelHeight;
    
    public BarCalc(int panelWidth, int panelHeight) {
      m_panelWidth = panelWidth;
      m_panelHeight = panelHeight;
    }
    
    public void run() {
      synchronized (m_locker) {
        // there is no use doing/displaying anything if the resolution
        // of the panel is less than the number of values for this attribute
        if (m_data.attribute(m_attribIndex).numValues() > m_panelWidth) {
          m_histBarClassCounts = null;
          m_threadRun = false;
          m_doneCurrentAttribute = true;
          m_displayCurrentAttribute = false;
          AttributeVisualizationPanel.this.repaint();
          return;
        }
        
        
        if((m_classIndex >= 0) &&
//                (m_data.attribute(m_classIndex).isNominal())) {
        		(m_data.attribute(m_attribIndex).isNominal())) {
                  SparseInstance histClassCounts[];
                  histClassCounts = new SparseInstance[m_data.attribute(m_attribIndex).numValues()];
                                          //[m_data.attribute(m_classIndex).numValues()+1];
                  
                  if (m_as.nominalWeights.length > 0) {
                    m_maxValue = m_as.nominalWeights[0];
                    for(int i=0; i<m_data.attribute(m_attribIndex).numValues(); i++) {
                      if(m_as.nominalWeights[i]>m_maxValue)
                	m_maxValue = m_as.nominalWeights[i];
                    }
                  }
                  else {
                    m_maxValue = 0;
                  }
                  
                  if(m_colorList.size()==0)
                    m_colorList.addElement(Color.black);
                  for(int i=m_colorList.size();
                  i < m_data.attribute(m_classIndex).numValues()+1; i++) {
                    Color pc = m_defaultColors[(i-1) % 10];
                    int ija =  (i-1) / 10;
                    ija *= 2;
                    
                    for (int j=0;j<ija;j++) {
                      pc = pc.darker();
                    }
                    
                    m_colorList.addElement(pc);
                  }
                  
                  // first sort data on attribute values
                  m_data.sort(m_attribIndex);
                  double[] tempClassCounts = null;
                  int tempAttValueIndex = -1;
                  
                  for(int k=0; k<m_data.numInstances(); k++) {
                    //System.out.println("attrib: "+
                    //                   m_data.instance(k).value(m_attribIndex)+
                    //                   " class: "+
                    //                   m_data.instance(k).value(m_classIndex));
                    if(!m_data.instance(k).isMissing(m_attribIndex)) {
                      // check to see if we need to allocate some space here
                      if (m_data.instance(k).value(m_attribIndex) != tempAttValueIndex) {
                        if (tempClassCounts != null) {
                          // set up the sparse instance for the previous bar (if any)
                          int numNonZero = 0;
                          for (int z = 0; z < tempClassCounts.length; z++) {
                            if (tempClassCounts[z] > 0) {
                              numNonZero++;
                            }
                          }
                          double[] nonZeroVals = new double[numNonZero];
                          int[] nonZeroIndices = new int[numNonZero];
                          int count = 0;
                          for (int z = 0; z < tempClassCounts.length; z++) {
                            if (tempClassCounts[z] > 0) {
                              nonZeroVals[count] = tempClassCounts[z];
                              nonZeroIndices[count++] = z;
                            }
                          }
                          SparseInstance tempS = 
                            new SparseInstance(1.0, nonZeroVals, nonZeroIndices, tempClassCounts.length);
                          histClassCounts[tempAttValueIndex] = tempS;
                        }
                        
                        tempClassCounts = new double[m_data.attribute(m_classIndex).numValues() + 1];
                        tempAttValueIndex = (int)m_data.instance(k).value(m_attribIndex);
                        
                        /* histClassCounts[(int)m_data.instance(k).value(m_attribIndex)] = 
                          new double[m_data.attribute(m_classIndex).numValues()+1]; */ 
                      }
                      if(m_data.instance(k).isMissing(m_classIndex)) {
                        /* histClassCounts[(int)m_data.instance(k).value(m_attribIndex)]
                                       [0] += m_data.instance(k).weight(); */
                        tempClassCounts[0] += m_data.instance(k).weight();
                      } else {
                        tempClassCounts[(int)m_data.instance(k).value(m_classIndex)+1] 
                                        += m_data.instance(k).weight();
                        
                        /*histClassCounts[(int)m_data.instance(k).value(m_attribIndex)]
                                      [(int)m_data.instance(k).value(m_classIndex)+1] += m_data.instance(k).weight();*/
                      }
                    }
                  }
                  
                  // set up sparse instance for last bar?
                  if (tempClassCounts != null) {
                    // set up the sparse instance for the previous bar (if any)
                    int numNonZero = 0;
                    for (int z = 0; z < tempClassCounts.length; z++) {
                      if (tempClassCounts[z] > 0) {
                        numNonZero++;
                      }
                    }
                    double[] nonZeroVals = new double[numNonZero];
                    int[] nonZeroIndices = new int[numNonZero];
                    int count = 0;
                    for (int z = 0; z < tempClassCounts.length; z++) {
                      if (tempClassCounts[z] > 0) {
                        nonZeroVals[count] = tempClassCounts[z];
                        nonZeroIndices[count++] = z;
                      }
                    }
                    SparseInstance tempS = 
                      new SparseInstance(1.0, nonZeroVals, nonZeroIndices, tempClassCounts.length);
                    histClassCounts[tempAttValueIndex] = tempS;
                  }
                  
                  //for(int i=0; i<histClassCounts.length; i++) {
                  //int sum=0;
                  //for(int j=0; j<histClassCounts[i].length; j++) {
                  //    sum = sum+histClassCounts[i][j];
                  //}
                  //System.out.println("histCount: "+sum+" Actual: "+
                  //                   m_as.nominalWeights[i]);
                  //}
                  
                  m_threadRun=false;
                  m_doneCurrentAttribute = true;
                  m_displayCurrentAttribute = true;
                  m_histBarClassCounts = histClassCounts;
                  //Image tmpImg = new BufferedImage(getWidth(), getHeight(),
                  //                                 BufferedImage.TYPE_INT_RGB);
                  //drawGraph( tmpImg.getGraphics() );
                  //img = tmpImg;
                  AttributeVisualizationPanel.this.repaint();
                }
        
        //RANKING BEGIN
        //paint attribute evaluation for ranking attribute.
        if((m_classIndex >= 0) &&
        (/*m_data.attribute(m_classIndex).isNominal()||*/m_data.attribute(m_attribIndex).isRanking() && m_data.attribute(m_classIndex).isRanking() )) {
          SparseInstance histClassCounts[];
          histClassCounts = new SparseInstance[m_data.attribute(m_attribIndex).numValues()];
                                  //[m_data.attribute(m_classIndex).numValues()+1];
          
          if (m_as.nominalWeights.length > 0) {
            m_maxValue = m_as.nominalWeights[0];
            for(int i=0; i<m_data.attribute(m_attribIndex).numValues(); i++) {
              if(m_as.nominalWeights[i]>m_maxValue)
        	m_maxValue = m_as.nominalWeights[i];
            }
          }
          else {
            m_maxValue = 0;
          }
          
          if(m_colorList.size()==0)
            m_colorList.addElement(Color.black);
          for(int i=m_colorList.size();
          i < m_data.attribute(m_classIndex).numValues()+1; i++) {
            Color pc = m_defaultColors[(i-1) % 10];
            int ija =  (i-1) / 10;
            ija *= 2;
            
            for (int j=0;j<ija;j++) {
              pc = pc.darker();
            }
            
            m_colorList.addElement(pc);
          }
          
          // first sort data on attribute values
          m_data.sort(m_attribIndex);
          double[] tempClassCounts = null;
          int tempAttValueIndex = -1;
          
          for(int k=0; k<m_data.numInstances(); k++) {
            //System.out.println("attrib: "+
            //                   m_data.instance(k).value(m_attribIndex)+
            //                   " class: "+
            //                   m_data.instance(k).value(m_classIndex));
            if(!m_data.instance(k).isMissing(m_attribIndex)) {
              // check to see if we need to allocate some space here
              if (m_data.instance(k).value(m_attribIndex) != tempAttValueIndex) {
                if (tempClassCounts != null) {
                  // set up the sparse instance for the previous bar (if any)
                  int numNonZero = 0;
                  for (int z = 0; z < tempClassCounts.length; z++) {
                    if (tempClassCounts[z] > 0) {
                      numNonZero++;
                    }
                  }
                  double[] nonZeroVals = new double[numNonZero];
                  int[] nonZeroIndices = new int[numNonZero];
                  int count = 0;
                  for (int z = 0; z < tempClassCounts.length; z++) {
                    if (tempClassCounts[z] > 0) {
                      nonZeroVals[count] = tempClassCounts[z];
                      nonZeroIndices[count++] = z;
                    }
                  }
                  SparseInstance tempS = 
                    new SparseInstance(1.0, nonZeroVals, nonZeroIndices, tempClassCounts.length);
                  histClassCounts[tempAttValueIndex] = tempS;
                }
                
                tempClassCounts = new double[m_data.attribute(m_classIndex).numValues() + 1];
                tempAttValueIndex = (int)m_data.instance(k).value(m_attribIndex);
                
                /* histClassCounts[(int)m_data.instance(k).value(m_attribIndex)] = 
                  new double[m_data.attribute(m_classIndex).numValues()+1]; */
              }
              if(m_data.instance(k).isMissing(m_classIndex)) {
                /* histClassCounts[(int)m_data.instance(k).value(m_attribIndex)]
                               [0] += m_data.instance(k).weight(); */
                tempClassCounts[0] += m_data.instance(k).weight();
              } else {
                tempClassCounts[(int)m_data.instance(k).value(m_classIndex)+1] 
                                += m_data.instance(k).weight();
                
                /*histClassCounts[(int)m_data.instance(k).value(m_attribIndex)]
                              [(int)m_data.instance(k).value(m_classIndex)+1] += m_data.instance(k).weight();*/
              }
            }
          }
          
          // set up sparse instance for last bar?
          if (tempClassCounts != null) {
            // set up the sparse instance for the previous bar (if any)
            int numNonZero = 0;
            for (int z = 0; z < tempClassCounts.length; z++) {
              if (tempClassCounts[z] > 0) {
                numNonZero++;
              }
            }
            double[] nonZeroVals = new double[numNonZero];
            int[] nonZeroIndices = new int[numNonZero];
            int count = 0;
            for (int z = 0; z < tempClassCounts.length; z++) {
              if (tempClassCounts[z] > 0) {
                nonZeroVals[count] = tempClassCounts[z];
                nonZeroIndices[count++] = z;
              }
            }
            
            
           
            if(m_data.attribute(m_classIndex)instanceof PreferenceAttribute || m_data.instance(0) instanceof PreferenceDenseInstance ){
            	
            	//Making a backup. If filters are applied, this backup is used.
            	if(m_as.pairRankWeights!=null){
            		backupPairRankWeights = m_as.pairRankWeights;
            	}
            	else{
            		m_as.pairRankWeights = backupPairRankWeights;
            	}
            	if(m_data.getLabels().size()==0)
            			m_data.setLabels(RankUtilities.labels);

            	if(m_as.pairRankWeights.length==0)m_as = m_data.attributeStats(m_data.classIndex());
            	           	           	

            	if(histClassCounts.length != Math.pow(m_as.pairRankWeights.length,2)){
            		histClassCounts = new SparseInstance[(int)Math.pow(m_as.pairRankWeights.length,2)];
            	}
            	
            	try{
            		m_histBarCounts = new double[(int)Math.pow(m_as.pairRankWeights.length,2)];
            	}
            	catch(Exception e){
            		System.out.println("Exception in AttributeVisualization.");
            	}
            	
            	double[] nzero = new double[(int)Math.pow(m_as.pairRankWeights.length,2)];
            	int[] nzeroIndex = new int[(int)Math.pow(m_as.pairRankWeights.length,2)];
            	
          	  int cnt=0;
          	  for(int a=0; a<m_as.pairRankWeights.length; a++){
          		  for(int b=0; b<m_as.pairRankWeights.length; b++){
          			  nzero[cnt] = m_as.pairRankWeights[a][b];
          			  nzeroIndex[cnt] = cnt;
              		  SparseInstance si = new SparseInstance(1.0, nzero, nzeroIndex, (int)Math.pow(m_as.pairRankWeights.length,2));
              		  histClassCounts[cnt] = si;      
                  	nzero = new double[(int)Math.pow(m_as.pairRankWeights.length,2)];
                	nzeroIndex = new int[(int)Math.pow(m_as.pairRankWeights.length,2)];
          			  cnt++;
          		  }
          	  }
            }
            else{
            	SparseInstance tempS = 
            		new SparseInstance(1.0, nonZeroVals, nonZeroIndices, tempClassCounts.length);
            	histClassCounts[tempAttValueIndex] = tempS;
            }
            //RANKING END
            AttributeVisualizationPanel.this.repaint();
            
          }
          
          //for(int i=0; i<histClassCounts.length; i++) {
          //int sum=0;
          //for(int j=0; j<histClassCounts[i].length; j++) {
          //    sum = sum+histClassCounts[i][j];
          //}
          //System.out.println("histCount: "+sum+" Actual: "+
          //                   m_as.nominalWeights[i]);
          //}
          
          
          
          m_threadRun=false;
          m_doneCurrentAttribute = true;
          m_displayCurrentAttribute = true;
          m_histBarClassCounts = histClassCounts;
          //Image tmpImg = new BufferedImage(getWidth(), getHeight(),
          //                                 BufferedImage.TYPE_INT_RGB);
          //drawGraph( tmpImg.getGraphics() );
          //img = tmpImg;
          AttributeVisualizationPanel.this.repaint();
        }
        else {
          double histCounts[];
          histCounts  = new double[m_data.attribute(m_attribIndex).numValues()];
          
        	  if (m_as.nominalWeights.length > 0) {
        		  m_maxValue = m_as.nominalWeights[0];
        		  for(int i=0; i<m_data.attribute(m_attribIndex).numValues(); i++) {
        			  if(m_as.nominalWeights[i]>m_maxValue)
        				  m_maxValue = m_as.nominalWeights[i];
        		  }
        	  }
        	  else {
        		  m_maxValue = 0;
        	  }
          
        	  for(int k=0; k<m_data.numInstances(); k++) {
        		  if(!m_data.instance(k).isMissing(m_attribIndex))
        			  histCounts[(int)m_data.instance(k).value(m_attribIndex)] += 
        				  m_data.instance(k).weight();
        	  }
          m_threadRun=false;
          m_displayCurrentAttribute = true;
          m_doneCurrentAttribute = true;
          m_histBarCounts = histCounts;
          //Image tmpImg = new BufferedImage(getWidth(), getHeight(),
          //                                 BufferedImage.TYPE_INT_RGB);
          //drawGraph( tmpImg.getGraphics() );
          //img = tmpImg;
          AttributeVisualizationPanel.this.repaint();
        }
      } //end synchronized
      //end run()
    }
  }
  
  /**
   * Internal class that calculates the histogram to display, in a separate
   * thread. In particular it initializes some of the crucial internal fields
   * required by paintComponent() to display the histogram for the current
   * attribute. These include: m_histBarCounts or m_histBarClassCounts,
   * m_maxValue and m_colorList.
   */
  private class HistCalc extends Thread {
    public void run() {
      synchronized (m_locker) {
    	  //RANKING BEGIN
       //modified for handling ranking attributes, too.
    	 if((m_classIndex >= 0) &&
           (/*m_data.attribute(m_classIndex).isNominal() ||*/ m_data.attribute(m_classIndex).isRanking())) {
          
          int intervals; double intervalWidth=0.0;
          
          //This uses the M.P.Wand's method to calculate the histogram's
          //interval width. See "Data-Based Choice of Histogram Bin Width", in
          //The American Statistician, Vol. 51, No. 1, Feb., 1997, pp. 59-64.
          //intervalWidth = Math.pow(6D/( -psi(2, g21())*m_data.numInstances()),
          //                          1/3D );
          
          //This uses the Scott's method to calculate the histogram's interval
          //width. See "On optimal and data-based histograms".
          // See Biometrika, 66, 605-610 OR see the same paper mentioned above.
          intervalWidth =  3.49 * m_as.numericStats.stdDev *
                           Math.pow(m_data.numInstances(), -1/3D);
          //The Math.max is introduced to remove the possibility of
          //intervals=0 and =NAN that can happen if respectively all the numeric
          //values are the same or the interval width is evaluated to zero.
          intervals = Math.max(1,
          (int)Math.round( (m_as.numericStats.max - m_as.numericStats.min) /
                           intervalWidth) );
          
          //System.out.println("Max: "+m_as.numericStats.max+
          //                   " Min: "+m_as.numericStats.min+
          //                   " stdDev: "+m_as.numericStats.stdDev+
          //                   "intervalWidth: "+intervalWidth);
          
          //The number 4 below actually represents a padding of 3 pixels on
          //each side of the histogram, and is also reflected in other parts of 
          //the code in the shape of numerical constants like "6" here.
          if(intervals > AttributeVisualizationPanel.this.getWidth()) {
            intervals = AttributeVisualizationPanel.this.getWidth()-6;
            if(intervals<1)//if width is too small then use 1 and forget padding
              intervals = 1;
          }
          double histClassCounts[][]  =
                          new double[intervals]
                                 [m_data.attribute(m_classIndex).numValues()+1];
          
          double barRange   = (m_as.numericStats.max - m_as.numericStats.min) /
                              (double)histClassCounts.length;
                 
          m_maxValue = 0;
          
          if(m_colorList.size()==0)
            m_colorList.addElement(Color.black);
          for(int i = m_colorList.size();
          i < m_data.attribute(m_classIndex).numValues()+1; i++) {
            Color pc = m_defaultColors[(i-1) % 10];
            int ija =  (i-1) / 10;
            ija *= 2;
            for (int j=0;j<ija;j++) {
              pc = pc.darker();
            }
                       
            m_colorList.addElement(pc);
          }
          
          for(int k=0; k<m_data.numInstances(); k++) {
            int t=0; //This holds the interval that the attibute value of the
                     //new instance belongs to.
            try {
              if(!m_data.instance(k).isMissing(m_attribIndex)) {
                //1. see footnote at the end of this file
                t = (int)Math.ceil( (float)(
                (m_data.instance(k).value(m_attribIndex)-m_as.numericStats.min)
                / barRange) );
                if(t==0) {
                  if(m_data.instance(k).isMissing(m_classIndex))
                    histClassCounts[t][0] += m_data.instance(k).weight();
                  else
                    histClassCounts[t][(int)m_data.instance(k).value(m_classIndex)+1] +=
                      m_data.instance(k).weight();
                  //if(histCounts[t]>m_maxValue)
                  //  m_maxValue = histCounts[t];
                }
                else {
                  if(m_data.instance(k).isMissing(m_classIndex))
                    histClassCounts[t-1][0] += m_data.instance(k).weight();
                  else
                    histClassCounts[t-1][(int)m_data.instance(k).value(m_classIndex)+1] +=
                      m_data.instance(k).weight();
                  //if(histCounts[t-1]>m_maxValue)
                  //  m_maxValue = histCounts[t-1];
                }
              }
            }
            catch(ArrayIndexOutOfBoundsException ae) {
              System.out.println("t:"+(t)+
              " barRange:"+barRange+
              " histLength:"+histClassCounts.length+
              " value:"+m_data.instance(k).value(m_attribIndex)+
              " min:"+m_as.numericStats.min+
              " sumResult:"+
              (m_data.instance(k).value(m_attribIndex) -
              m_as.numericStats.min)+
              " divideResult:"+
              (float)((m_data.instance(k).value(m_attribIndex) -
              m_as.numericStats.min) / barRange)+
              " finalResult:"+
              Math.ceil((float)((m_data.instance(k).value(m_attribIndex)-
              m_as.numericStats.min) / barRange)) );
            }
          }
          for(int i=0; i<histClassCounts.length; i++) {
            double sum=0;
            for(int j=0; j<histClassCounts[i].length; j++)
              sum = sum+histClassCounts[i][j];
            if(m_maxValue<sum)
              m_maxValue = sum;
          }
          
          // convert to sparse instances
          SparseInstance[] histClassCountsSparse = 
            new SparseInstance[histClassCounts.length];
          
          for (int i = 0; i < histClassCounts.length; i++) {
            int numSparseValues = 0;
            for (int j = 0; j < histClassCounts[i].length; j++) {
              if (histClassCounts[i][j] > 0) {
                numSparseValues++;
              }
            }
            double[] sparseValues = new double[numSparseValues];
            int[] sparseIndices = new int[numSparseValues];
            int count = 0;
            for (int j = 0; j < histClassCounts[i].length; j++) {
              if (histClassCounts[i][j] > 0) {
                sparseValues[count] = histClassCounts[i][j];
                sparseIndices[count++] = j;
              }
            }
            
            SparseInstance tempS = 
              new SparseInstance(1.0, sparseValues, sparseIndices, 
                  histClassCounts[i].length);
            histClassCountsSparse[i] = tempS;
            
          }
          
          m_histBarClassCounts = histClassCountsSparse;
          m_barRange =  barRange;
          
        }
    	//RANKING END
    	 if((m_classIndex >= 0) &&
    	           (m_data.attribute(m_attribIndex).isNominal())) {
    	          
    	          int intervals; double intervalWidth=0.0;
    	          
    	          //This uses the M.P.Wand's method to calculate the histogram's
    	          //interval width. See "Data-Based Choice of Histogram Bin Width", in
    	          //The American Statistician, Vol. 51, No. 1, Feb., 1997, pp. 59-64.
    	          //intervalWidth = Math.pow(6D/( -psi(2, g21())*m_data.numInstances()),
    	          //                          1/3D );
    	          
    	          //This uses the Scott's method to calculate the histogram's interval
    	          //width. See "On optimal and data-based histograms".
    	          // See Biometrika, 66, 605-610 OR see the same paper mentioned above.
    	          intervalWidth =  3.49 * m_as.numericStats.stdDev *
    	                           Math.pow(m_data.numInstances(), -1/3D);
    	          //The Math.max is introduced to remove the possibility of
    	          //intervals=0 and =NAN that can happen if respectively all the numeric
    	          //values are the same or the interval width is evaluated to zero.
    	          intervals = Math.max(1,
    	          (int)Math.round( (m_as.numericStats.max - m_as.numericStats.min) /
    	                           intervalWidth) );
    	          
    	          //System.out.println("Max: "+m_as.numericStats.max+
    	          //                   " Min: "+m_as.numericStats.min+
    	          //                   " stdDev: "+m_as.numericStats.stdDev+
    	          //                   "intervalWidth: "+intervalWidth);
    	          
    	          //The number 4 below actually represents a padding of 3 pixels on
    	          //each side of the histogram, and is also reflected in other parts of 
    	          //the code in the shape of numerical constants like "6" here.
    	          if(intervals > AttributeVisualizationPanel.this.getWidth()) {
    	            intervals = AttributeVisualizationPanel.this.getWidth()-6;
    	            if(intervals<1)//if width is too small then use 1 and forget padding
    	              intervals = 1;
    	          }
    	          double histClassCounts[][]  =
    	                          new double[intervals]
    	                                 [m_data.attribute(m_classIndex).numValues()+1];
    	          
    	          double barRange   = (m_as.numericStats.max - m_as.numericStats.min) /
    	                              (double)histClassCounts.length;
    	          
    	          m_maxValue = 0;
    	          
    	          if(m_colorList.size()==0)
    	            m_colorList.addElement(Color.black);
    	          for(int i = m_colorList.size();
    	          i < m_data.attribute(m_classIndex).numValues()+1; i++) {
    	            Color pc = m_defaultColors[(i-1) % 10];
    	            int ija =  (i-1) / 10;
    	            ija *= 2;
    	            for (int j=0;j<ija;j++) {
    	              pc = pc.darker();
    	            }
    	            m_colorList.addElement(pc);
    	          }
    	          
    	          for(int k=0; k<m_data.numInstances(); k++) {
    	            int t=0; //This holds the interval that the attibute value of the
    	                     //new instance belongs to.
    	            try {
    	              if(!m_data.instance(k).isMissing(m_attribIndex)) {
    	                //1. see footnote at the end of this file
    	                t = (int)Math.ceil( (float)(
    	                (m_data.instance(k).value(m_attribIndex)-m_as.numericStats.min)
    	                / barRange) );
    	                if(t==0) {
    	                  if(m_data.instance(k).isMissing(m_classIndex))
    	                    histClassCounts[t][0] += m_data.instance(k).weight();
    	                  else
    	                    histClassCounts[t][(int)m_data.instance(k).value(m_classIndex)+1] +=
    	                      m_data.instance(k).weight();
    	                  //if(histCounts[t]>m_maxValue)
    	                  //  m_maxValue = histCounts[t];
    	                }
    	                else {
    	                  if(m_data.instance(k).isMissing(m_classIndex))
    	                    histClassCounts[t-1][0] += m_data.instance(k).weight();
    	                  else
    	                    histClassCounts[t-1][(int)m_data.instance(k).value(m_classIndex)+1] +=
    	                      m_data.instance(k).weight();
    	                  //if(histCounts[t-1]>m_maxValue)
    	                  //  m_maxValue = histCounts[t-1];
    	                }
    	              }
    	            }
    	            catch(ArrayIndexOutOfBoundsException ae) {
    	              System.out.println("t:"+(t)+
    	              " barRange:"+barRange+
    	              " histLength:"+histClassCounts.length+
    	              " value:"+m_data.instance(k).value(m_attribIndex)+
    	              " min:"+m_as.numericStats.min+
    	              " sumResult:"+
    	              (m_data.instance(k).value(m_attribIndex) -
    	              m_as.numericStats.min)+
    	              " divideResult:"+
    	              (float)((m_data.instance(k).value(m_attribIndex) -
    	              m_as.numericStats.min) / barRange)+
    	              " finalResult:"+
    	              Math.ceil((float)((m_data.instance(k).value(m_attribIndex)-
    	              m_as.numericStats.min) / barRange)) );
    	            }
    	          }
    	          for(int i=0; i<histClassCounts.length; i++) {
    	            double sum=0;
    	            for(int j=0; j<histClassCounts[i].length; j++)
    	              sum = sum+histClassCounts[i][j];
    	            if(m_maxValue<sum)
    	              m_maxValue = sum;
    	          }
    	          
    	          // convert to sparse instances
    	          SparseInstance[] histClassCountsSparse = 
    	            new SparseInstance[histClassCounts.length];
    	          
    	          for (int i = 0; i < histClassCounts.length; i++) {
    	            int numSparseValues = 0;
    	            for (int j = 0; j < histClassCounts[i].length; j++) {
    	              if (histClassCounts[i][j] > 0) {
    	                numSparseValues++;
    	              }
    	            }
    	            double[] sparseValues = new double[numSparseValues];
    	            int[] sparseIndices = new int[numSparseValues];
    	            int count = 0;
    	            for (int j = 0; j < histClassCounts[i].length; j++) {
    	              if (histClassCounts[i][j] > 0) {
    	                sparseValues[count] = histClassCounts[i][j];
    	                sparseIndices[count++] = j;
    	              }
    	            }
    	            
    	            SparseInstance tempS = 
    	              new SparseInstance(1.0, sparseValues, sparseIndices, 
    	                  histClassCounts[i].length);
    	            histClassCountsSparse[i] = tempS;
    	            
    	          }
    	          
    	          m_histBarClassCounts = histClassCountsSparse;
    	          m_barRange =  barRange;
    	          
    	        }
    	 
    	 
    	 
        else { //else if the class attribute is numeric or the class is not set
          
          int intervals; double intervalWidth;
          //At the time of this coding the
          //possibility of datasets with zero instances
          //was being dealt with in the
          //PreProcessPanel of weka Explorer.
          
          //old method of calculating number of intervals
          //intervals =  m_as.totalCount>10 ?
          //                  (int)(m_as.totalCount*0.1):(int)m_as.totalCount;
          
          //This uses the M.P.Wand's method to calculate the histogram's
          //interval width. See "Data-Based Choice of Histogram Bin Width", in
          //The American Statistician, Vol. 51, No. 1, Feb., 1997, pp. 59-64.
          //intervalWidth = Math.pow(6D/(-psi(2, g21())*m_data.numInstances() ),
          //                          1/3D );
          
          //This uses the Scott's method to calculate the histogram's interval
          //width. See "On optimal and data-based histograms".
          // See Biometrika, 66, 605-610 OR see the same paper mentioned above.
          intervalWidth =  3.49 * m_as.numericStats.stdDev *
                           Math.pow(m_data.numInstances(), -1/3D);
          //The Math.max is introduced to remove the possibility of
          //intervals=0 and =NAN that can happen if respectively all the numeric
          //values are the same or the interval width is evaluated to zero.
          intervals = Math.max(1,
          (int)Math.round( (m_as.numericStats.max - m_as.numericStats.min) /
                           intervalWidth) );
          
          //The number 4 below actually represents a padding of 3 pixels on
          //each side of the histogram, and is also reflected in other parts of 
          //the code in the shape of numerical constants like "6" here.
          if(intervals > AttributeVisualizationPanel.this.getWidth()) {
            intervals = AttributeVisualizationPanel.this.getWidth()-6;
            if(intervals<1)
              intervals = 1;
          }
          
          double[] histCounts  = new double[intervals];
          double barRange   = (m_as.numericStats.max - m_as.numericStats.min) /
                              (double)histCounts.length;
          
          m_maxValue = 0;
          
          for(int k=0; k<m_data.numInstances(); k++) {
            int t=0; //This holds the interval to which the current attribute's
                    //value of this particular instance k belongs to.
            
            if(m_data.instance(k).isMissing(m_attribIndex)) 
              continue; //ignore missing values
            
            try {
              //1. see footnote at the end of this file
              t =(int) Math.ceil((
              (m_data.instance(k).value(m_attribIndex)-m_as.numericStats.min)
              / barRange));
              if(t==0) {
                histCounts[t] += m_data.instance(k).weight();
                if(histCounts[t]>m_maxValue)
                  m_maxValue = histCounts[t];
              }
              else {
                histCounts[t-1] += m_data.instance(k).weight();
                if(histCounts[t-1]>m_maxValue)
                  m_maxValue = histCounts[t-1];
              }
            }
            catch(ArrayIndexOutOfBoundsException ae) {
              ae.printStackTrace();
              System.out.println("t:"+(t)+
              " barRange:"+barRange+
              " histLength:"+histCounts.length+
              " value:"+m_data.instance(k).value(m_attribIndex)+
              " min:"+m_as.numericStats.min+
              " sumResult:"+
              (m_data.instance(k).value(m_attribIndex)-m_as.numericStats.min)+
              " divideResult:"+
              (float)((m_data.instance(k).value(m_attribIndex) -
              m_as.numericStats.min)/barRange)+
              " finalResult:"+
              Math.ceil( (float)((m_data.instance(k).value(m_attribIndex) -
              m_as.numericStats.min) / barRange)) );
            }
          }
          m_histBarCounts = histCounts;
          m_barRange =  barRange;
        }
        
        m_threadRun=false;
        m_displayCurrentAttribute = true;
        m_doneCurrentAttribute = true;
        //Image tmpImg = new BufferedImage(getWidth(), getHeight(),
        //                                 BufferedImage.TYPE_INT_RGB);
        //drawGraph( tmpImg.getGraphics() );
        //img = tmpImg;
        AttributeVisualizationPanel.this.repaint();
      }
    }
    
    /****Code for M.P.Wand's method of histogram bin width selection.
     *   There is some problem with it. It always comes up -ve value
     *   which is raised to the power 1/3 and gives an NAN.
     * private static final int M=400;
     * private double psi(int r, double g) {
     * double val;
     *
     * double sum=0.0;
     * for(int i=0; i<M; i++) {
     * double valCjKj=0.0;
     * for(int j=0; j<M; j++) {
     * valCjKj += c(j) * k(r, j-i, g);
     * }
     * sum += valCjKj*c(i);
     * }
     *
     * val = Math.pow(m_data.numInstances(), -2) * sum;
     * //System.out.println("psi returns: "+val);
     * return val;
     * }
     * private double g21() {
     * double val;
     *
     * val = Math.pow(2 / ( Math.sqrt(2D*Math.PI)*psi(4, g22()) * 
     *                      m_data.numInstances() ), 1/5D)
     *       * Math.sqrt(2) * m_as.numericStats.stdDev;
     * //System.out.println("g21 returns: "+val);
     * return val;
     * }
     * private double g22() {
     * double val;
     *
     * val = Math.pow( 2D/(5*m_data.numInstances()), 1/7D) * 
     *       Math.sqrt(2) * m_as.numericStats.stdDev;
     * //System.out.println("g22 returns: "+val);
     * return val;
     * }
     * private double c(int j) {
     * double val=0.0;
     * double sigma = (m_as.numericStats.max - m_as.numericStats.min)/(M-1);
     *
     * //System.out.println("In c before doing the sum we have");
     * //System.out.println("max: " +m_as.numericStats.max+" min: "+
     * //                   m_as.numericStats.min+" sigma: "+sigma);
     *
     * for(int i=0; i<m_data.numInstances(); i++) {
     * if(!m_data.instance(i).isMissing(m_attribIndex))
     * val += Math.max( 0,
     * ( 1 - Math.abs( Math.pow(sigma, -1)*(m_data.instance(i).value(m_attribIndex) - j) ) )
     * );
     * }
     * //System.out.println("c returns: "+val);
     * return val;
     * }
     * private double k(int r, int j, double g) {
     * double val;
     * double sigma = (m_as.numericStats.max - m_as.numericStats.min)/(M-1);
     * //System.out.println("Before calling L we have");
     * //System.out.println("Max: "+m_as.numericStats.max+" Min: "+m_as.numericStats.min+"\n"+
     * //			 "r: "+r+" j: "+j+" g: "+g);
     * val = Math.pow( g, -r-1) * L(sigma*j/g);
     * //System.out.println("k returns: "+val);
     * return val;
     * }
     * private double L(double x) {
     * double val;
     *
     * val = Math.pow( 2*Math.PI, -1/2D ) * Math.exp( -(x*x)/2D );
     * //System.out.println("L returns: "+val);
     * return val;
     * }
     *******End of Wand's method
     */
  }
  
  
  /**
   * Returns "&lt;nominal value&gt; [&lt;nominal value count&gt;]"
   * if displaying a bar plot and mouse is on some bar.
   * If displaying histogram then it
   *     <li>returns "count &lt;br&gt; [&lt;bars Range&gt;]" if mouse is
   *     on the first bar. </li>
   *     <li>returns "count &lt;br&gt; (&lt;bar's Range&gt;]" if mouse is
   *     on some bar other than the first one. </li>
   * Otherwise it returns ""
   *
   * @param ev The mouse event
   */
  public String getToolTipText(MouseEvent ev) { 
    if(m_as!=null && m_as.nominalWeights!=null) { //if current attrib is nominal
      
      float intervalWidth = this.getWidth() / (float)m_as.nominalWeights.length;
      double heightRatio;      
      int barWidth, x=0, y=0;
      
      //if intervalWidth is at least six then bar width is 80% of intervalwidth
      if(intervalWidth>5) //the rest is padding
        barWidth = (int)Math.floor(intervalWidth*0.8F);
      else
        barWidth = 1;  //Otherwise barwidth is 1 & padding would be at least 1.
      
      //initializing x to maximum of 1 or 10% of interval width (i.e. half of 
      //the padding which is 20% of interval width, as there is 10% on each 
      //side of the bar) so that it points to the start of the 1st bar
      x = x + (int)( (Math.floor(intervalWidth*0.1F))<1 ? 
                     1:(Math.floor(intervalWidth*0.1F)) );

      //Adding to x the appropriate value so that it points to the 1st bar of 
      //our "centered" barplot. If subtracting barplots width from panel width 
      //gives <=2 then the barplot is not centered.
      if(this.getWidth() - (m_as.nominalWeights.length*barWidth+
                           (int)( (Math.floor(intervalWidth*0.2F))<1 ? 
                                   1:(Math.floor(intervalWidth*0.2F)) ) * 
                           m_as.nominalWeights.length) > 2 ) {
        
        //The following amounts to adding to x the half of the area left after
        //subtracting from the components width the width of the whole barplot
        //(i.e. width of all the bars plus the width of all the bar paddings, 
        //thereby equaling to the whole barplot), since our barplot is centered.
        x += ( this.getWidth() - (m_as.nominalWeights.length*barWidth + 
                                 (int)( (Math.floor(intervalWidth*0.2F))<1 ? 
                                        1:(Math.floor(intervalWidth*0.2F)) ) * 
                                 m_as.nominalWeights.length) ) / 2;
      }
        
      for(int i=0; i<m_as.nominalWeights.length; i++) {        
        heightRatio = (this.getHeight()-(double)m_fm.getHeight())/m_maxValue;
        //initializing y to point to (from top) the start of the bar
        y = (int) (this.getHeight()-Math.round(m_as.nominalWeights[i]*heightRatio));
        
        //if our mouse is on a bar then return the count of this bar in our
        //barplot 
        if(ev.getX() >= x && ev.getX()<=x+barWidth && 
           ev.getY() >= this.getHeight() - 
                        Math.round(m_as.nominalWeights[i]*heightRatio) )
          return(m_data.attribute(m_attribIndex).value(i)+
                 " ["+Utils.doubleToString(m_as.nominalWeights[i], 3)+"]");
        //otherwise advance x to next bar and check that. Add barwidth to x 
        //and padding which is max(1, 20% of interval width)
        x = x+barWidth+(int)( (Math.floor(intervalWidth*0.2F))<1 ? 
                              1:(Math.floor(intervalWidth*0.2F)) );
      }
    }
    else if(m_threadRun==false &&     //if attrib is numeric
            (m_histBarCounts!=null || m_histBarClassCounts!=null)) {

      double heightRatio, intervalWidth;
      int x=0, y=0,  barWidth;
      double bar = m_as.numericStats.min;
      //if the class attribute is set and it is nominal or a ranking.
      if((m_classIndex >= 0) && (m_data.attribute(m_classIndex).isNominal() || m_data.attribute(m_classIndex).isRanking())) {
        //there is 3 pixels of padding on each side of the histogram
        //the barwidth is 1 if after removing the padding its width is less
        //then the displayable width
        barWidth = ((this.getWidth()-6)/m_histBarClassCounts.length)<1 ? 
                   1:((this.getWidth()-6)/m_histBarClassCounts.length);
        
        //initializing x to 3 adding appropriate value to make it point to the 
        //start of the 1st bar of our "centered" histogram.
        x = 3;
        if( (this.getWidth() - (x + m_histBarClassCounts.length*barWidth)) > 5 )
          x += (this.getWidth() - (x + m_histBarClassCounts.length*barWidth))/2;
        
        heightRatio = (this.getHeight()-(double)m_fm.getHeight())/m_maxValue;
        
        if( ev.getX()-x >= 0) {
          //The temp holds the index of the current interval that we are looking
          //at
          int temp = (int)((ev.getX()-x)/(barWidth+0.0000000001));
          if(temp == 0){  //handle the special case temp==0. see footnote 1
            double sum=0;
            for(int k=0; k<m_histBarClassCounts[0].numValues(); k++)
              sum += m_histBarClassCounts[0].valueSparse(k);
            //return the count of the interval mouse is pointing to plus 
            //the range of values that fall into this interval
            return ("<html><center><font face=Dialog size=-1>" + 
                Utils.doubleToString(sum, 3) + "<br>"+
                    "["+Utils.doubleToString(bar+m_barRange*temp,3)+
                    ", "+Utils.doubleToString((bar+m_barRange*(temp+1)),3)+
                    "]"+"</font></center></html>");
          }
          else if( temp < m_histBarClassCounts.length ) { //handle case temp!=0
            double sum=0;
            for(int k=0; k<m_histBarClassCounts[temp].numValues(); k++)
              sum+=m_histBarClassCounts[temp].valueSparse(k);
            //return the count of the interval mouse is pointing to plus 
            //the range of values that fall into this interval
            return ("<html><center><font face=Dialog size=-1>" + 
                Utils.doubleToString(sum, 3) + "<br>("+
                    Utils.doubleToString(bar+m_barRange*temp,3)+", "+
                    Utils.doubleToString((bar+m_barRange*(temp+1)),3)+
                    "]</font></center></html>");
          }
        }
      }
      else {  //else if the class attribute is not set or is numeric
        barWidth = ((this.getWidth()-6)/m_histBarCounts.length) < 1 ? 
                   1 : ((this.getWidth()-6)/m_histBarCounts.length);
        
        //initializing x to 3 adding appropriate value to make it point to the 
        //start of the 1st bar of our "centered" histogram.
        x = 3;
        if( (this.getWidth() - (x + m_histBarCounts.length*barWidth)) > 5 )
          x += (this.getWidth() - (x + m_histBarCounts.length*barWidth))/2;
        
        heightRatio = (this.getHeight()-(float)m_fm.getHeight())/m_maxValue;
        
        if( ev.getX()-x >= 0) {
          //Temp holds the index of the current bar we are looking at.
          int temp = (int)((ev.getX()-x)/(barWidth+0.0000000001));
          
          //return interval count as well as its range
          if(temp == 0) //handle special case temp==0. see footnote 1.
            return ("<html><center><font face=Dialog size=-1>"+
                    m_histBarCounts[0]+"<br>"+
                    "["+Utils.doubleToString(bar+m_barRange*temp,3)+", "+
                    Utils.doubleToString((bar+m_barRange*(temp+1)),3)+
                    "]"+
                    "</font></center></html>");
          else if(temp < m_histBarCounts.length) //handle case temp!=0
            return ("<html><center><font face=Dialog size=-1>"+
                    m_histBarCounts[temp]+"<br>"+
                    "("+Utils.doubleToString(bar+m_barRange*temp,3)+", "+
                    Utils.doubleToString((bar+m_barRange*(temp+1)),3)+
                    "]"+
                    "</font></center></html>");
        }
      }
    }
    return PrintableComponent.getToolTipText(m_Printer);
  }
  
  
  /**
   * Paints this component
   *
   * @param g The graphics object for this component
   */
  public void paintComponent(Graphics g) {
    g.clearRect(0,0,this.getWidth(), this.getHeight());
    
    if(m_as!=null) {    //If calculations have been done and histogram/barplot
      if (!m_doneCurrentAttribute && !m_threadRun) {
        calcGraph(this.getWidth(), this.getHeight());
      }
      if(m_threadRun==false && m_displayCurrentAttribute) {  //calculation thread is not running
        int buttonHeight=0;
        
        if(m_colorAttrib!=null)
          buttonHeight =m_colorAttrib.getHeight()+m_colorAttrib.getLocation().y;
     
        //if current attribute is nominal then draw barplot.
        if(m_as.nominalWeights != null && 
           (m_histBarClassCounts!=null || m_histBarCounts!=null) ) {
          double heightRatio, intervalWidth;
          int x=0, y=0, barHeight, barWidth;
          
          
          if((m_classIndex >= 0) && 
                  (m_data.attribute(m_attribIndex).isNominal())) {
                    
                 intervalWidth=(this.getWidth()/(float)m_histBarClassCounts.length);
                 
                 //Barwidth is 80% of interval width.The remaining 20% is padding,
                 //10% on each side of the bar. If interval width is less then 5 the
                 //20% of that value is less than 1, in that case we use bar width of
                 //1 and padding of 1 pixel on each side of the bar.
                 if(intervalWidth>5)
                   barWidth = (int)Math.floor(intervalWidth*0.8F);
                 else
                   barWidth = 1;

                 //initializing x to 10% of interval width or to 1 if 10% is <1. This
                 //is essentially the LHS padding of the 1st bar.
                 x = x + (int)( (Math.floor(intervalWidth*0.1F))<1 ?
                                1 : (Math.floor(intervalWidth*0.1F)) );
                 
                 //Add appropriate value to x so that it starts at the 1st bar of 
                 //a "centered" barplot.
                 if(this.getWidth() - (m_histBarClassCounts.length*barWidth + 
                                      (int)( (Math.floor(intervalWidth*0.2F))<1 ? 
                                             1 :(Math.floor(intervalWidth*0.2F))
                                           ) * m_histBarClassCounts.length) > 2 ) {
                   //We take the width of all the bars and all the paddings (20%
                   //of interval width), and subtract it from the width of the panel
                   //to get the extra space that would be left after drawing. We 
                   //divide that space by 2 to get its mid-point and add that to our
                   //x, thus making the whole bar plot drawn centered in our 
                   //component.
                   x += (this.getWidth()-(m_histBarClassCounts.length*barWidth+
                                         (int)( (Math.floor(intervalWidth*0.2F))<1 ?
                                                1 : (Math.floor(intervalWidth*0.2F))
                                              ) * m_histBarClassCounts.length))/2;
                 }
                 
                 //this holds the count of the bar and will be calculated by adding
                 //up the counts of individual subbars. It is displayed at the top
                 //of each bar.
                 double sum=0;
                 for(int i=0; i<m_histBarClassCounts.length; i++) {

                   //calculating the proportion of the components height compared to 
                   //the maxvalue in our attribute, also taking into account the 
                   //height of font to display bars count and the height of the class 
                   //ComboBox.
                   heightRatio = ( this.getHeight()-(double)m_fm.getHeight() - 
                       buttonHeight ) / m_maxValue;              
                   y=this.getHeight();
                   if (m_histBarClassCounts[i] != null) {
                     for(int j=0; j<m_histBarClassCounts[i].numAttributes(); j++) {
                       sum = sum + m_histBarClassCounts[i].value(j);
                       y = (int) (y-Math.round(m_histBarClassCounts[i].value(j) * heightRatio));
                       //selecting the colour corresponding to the current class.
                       if(m_data.attribute(m_classIndex).isRanking()){
                    	   g.setColor(Color.BLUE);
                       }
                       else
                    	   g.setColor( (Color)m_colorList.elementAt(j) );
                       g.fillRect(x, y, barWidth, 
                           (int) Math.round(m_histBarClassCounts[i].value(j) * heightRatio));
                       g.setColor(Color.black);
                     }
                   }
                   //drawing the bar count at the top of the bar if it is less than
                   //interval width. draw it 1px up to avoid touching the bar.
                   if(m_fm.stringWidth(Utils.doubleToString(sum, 1))<intervalWidth)
                     g.drawString(Utils.doubleToString(sum, 1), x, y-1);
                   //advancing x to the next bar by adding bar width and padding
                   //of both the bars (i.e. RHS padding of the bar just drawn and LHS
                   //padding of the new bar).
                   x = x+barWidth+(int)( (Math.floor(intervalWidth*0.2F))<1 ? 
                       1:(Math.floor(intervalWidth*0.2F)) );
                   //reseting sum for the next bar.
                   sum=0;

                 }
               }
          
        
          //RANKING BEGIN
          else if((m_classIndex >= 0) && 
             (m_data.attribute(m_attribIndex).isRanking())) {
              
        	  
            intervalWidth=(this.getWidth()/(float)m_histBarClassCounts.length);
            
            //Barwidth is 80% of interval width.The remaining 20% is padding,
            //10% on each side of the bar. If interval width is less then 5 the
            //20% of that value is less than 1, in that case we use bar width of
            //1 and padding of 1 pixel on each side of the bar.
            if(intervalWidth>5)
              barWidth = (int)Math.floor(intervalWidth*0.8F);
            else
              barWidth = 1;

            //initializing x to 10% of interval width or to 1 if 10% is <1. This
            //is essentially the LHS padding of the 1st bar.
            x = x + (int)( (Math.floor(intervalWidth*0.1F))<1 ?
                           1 : (Math.floor(intervalWidth*0.1F)) );
            
            //Add appropriate value to x so that it starts at the 1st bar of 
            //a "centered" barplot.
            if(this.getWidth() - (m_histBarClassCounts.length*barWidth + 
                                 (int)( (Math.floor(intervalWidth*0.2F))<1 ? 
                                        1 :(Math.floor(intervalWidth*0.2F))
                                      ) * m_histBarClassCounts.length) > 2 ) {
              //We take the width of all the bars and all the paddings (20%
              //of interval width), and subtract it from the width of the panel
              //to get the extra space that would be left after drawing. We 
              //divide that space by 2 to get its mid-point and add that to our
              //x, thus making the whole bar plot drawn centered in our 
              //component.
              x += (this.getWidth()-(m_histBarClassCounts.length*barWidth+
                                    (int)( (Math.floor(intervalWidth*0.2F))<1 ?
                                           1 : (Math.floor(intervalWidth*0.2F))
                                         ) * m_histBarClassCounts.length))/2;
            }
            
            //this holds the count of the bar and will be calculated by adding
            //up the counts of individual subbars. It is displayed at the top
            //of each bar.
            
        
            double sum=0;
            boolean ranking = false;
            Graphics2D g2d = (Graphics2D)g;
            for(int i=0; i<m_histBarClassCounts.length; i++) {
            	if(m_data.attribute(m_classIndex) instanceof PreferenceAttribute || m_data.instance(0) instanceof PreferenceDenseInstance) ranking = true;
              //calculating the proportion of the components height compared to 
              //the maxvalue in our attribute, also taking into account the 
              //height of font to display bars count and the height of the class 
              //ComboBox.
            	
                       
            	if(m_maxValue==0.0){
            		for(int a=0; a<m_as.pairRankWeights.length; a++){
            			for(int b=0; b<m_as.pairRankWeights.length; b++){
            				if(m_as.pairRankWeights[a][b]>m_maxValue){
            					m_maxValue = m_as.pairRankWeights[a][b];
            				}
            			}
            		}
            	}
            	
              heightRatio = ( this.getHeight()-(double)m_fm.getHeight() - 
                  buttonHeight ) / m_maxValue;              
              y=this.getHeight();
              
              	  
            	  if (m_histBarClassCounts[i] != null) {
            		  if(!ranking){
            			  try{
            				  for(int j=0; j<m_histBarClassCounts[i].numAttributes(); j++) {
            					  sum = sum + m_histBarClassCounts[i].value(j);
            					  y = (int) (y-Math.round(m_histBarClassCounts[i].value(j) * heightRatio));
            					  //selecting the colour corresponding to the current class.
            					  g2d.setColor( (Color)m_colorList.elementAt(j) );
            					  g2d.fillRect(x, y, barWidth, 
            			  			(int) Math.round(m_histBarClassCounts[i].value(j) * heightRatio));
            					  g2d.setColor(Color.black);
            				  }	
            			  }
            			  catch(Exception e){
            				  ranking =true;
            			  }
            		  }
            		  
            		  if(ranking){              
            				for(int k=0; k<m_as.pairRankWeights.length;k++){
            					for(int j=0; j<m_as.pairRankWeights.length; j++){ 
            						
            						double num = m_data.numInstances();
            						double act = m_as.pairRankWeights[k][j];
            						
            						int grey =(int) ((act/num)*255);
            						if(grey>255)grey=255;
            						g.setColor(new Color(255-grey,255-grey,255-grey));
            						
           						
                       		/*		if(m_as.pairRankWeights.length>9){
                       					int len = m_as.pairRankWeights.length;
                       					//int len = this.getWidth();
                       					
                       					int tmpX = 30+(50-2*len)*j+(50-2*len)/4;
                       						
                       						
                       						g2d.fill3DRect(30+(50-2*len)*j,(k+2)*(15-len/3),50-2*len, 15-len/3,true); 
                       						g2d.drawRect(30,2*(15-len/3),m_as.pairRankWeights.length*(50-2*len),m_as.pairRankWeights.length*(15-len/3));
                       						g2d.setColor(Color.black);
                       						g2d.drawString(">", 10, 15-len/3);
                       						g2d.drawString(m_data.labels.get(j)[0],tmpX,(15-len/3));
                       						g2d.drawString(m_data.labels.get(k)[0], 10, (k+3)*(15-len/3)-2);

                       				}*/
                       			//	else{   
                       					int marginLeftRight=30;
                       					int marginUpDown=30;
                       					int len = this.getWidth()-2*marginLeftRight; //subtracting margins.
                       				    int heigth = this.getHeight()-marginUpDown;
                       					int rectLen = len/m_as.pairRankWeights.length;
                       					int rectHeight = heigth/m_as.pairRankWeights.length;
                       					
                       					
                       					g2d.fill3DRect(marginLeftRight+rectLen*j, k*rectHeight+marginUpDown,rectLen,rectHeight,true);
                       			
                       					g2d.setColor(Color.black);
                       					g2d.drawRect(marginLeftRight,marginUpDown,rectLen*(m_as.pairRankWeights.length), (m_as.pairRankWeights.length)*rectHeight);
                       					
                       					int fontSize = 25-m_as.pairRankWeights.length;
                       					if(fontSize<=0)fontSize=5;
                       					Font visFont =new Font("Arial", Font.ITALIC|Font.PLAIN, fontSize);
                       					g2d.setFont(visFont);
                       					g2d.drawString(">", 0,  marginUpDown-fontSize/2);  
//                       					int labelSize=m_data.labels.get(j)[0].length();
                       					//if(labelSize<=2)
                       						g2d.drawString(m_data.labels.get(j)[0],marginLeftRight+rectLen*(j+1)-rectLen/2-rectLen/3, marginUpDown-fontSize/2);
                       					//else
                       						//g2d.drawString(m_data.labels.get(j)[0],marginLeftRight+rectLen*(j+1)-labelSize-rectLen/2, marginUpDown-fontSize/2);
                       					g2d.drawString(m_data.labels.get(k)[0], 0, marginUpDown+rectHeight*(k+1)-rectHeight/5);
                       					
                       					/*
                       					g2d.fill3DRect(30+40*j,(k+2)*15,40, 15,true);                     					
                       					g2d.drawRect(30, 30, 40*(m_as.pairRankWeights.length), (m_as.pairRankWeights.length)*15);
                       					g2d.setColor(Color.black);
                       					g2d.drawString(">", 5, 20);
                       					g2d.drawString(m_data.labels.get(j)[0],50+40*j , 20);
                       					g2d.drawString(m_data.labels.get(k)[0], 5, (k+2)*15+10);
                       					*/
                       			//	}
            					}
            				}            				
            		  }
            	  
            	  }             
              //drawing the bar count at the top of the bar if it is less than
              //interval width. draw it 1px up to avoid touching the bar.
              if(!ranking){
            	  if(m_fm.stringWidth(Utils.doubleToString(sum, 1))<intervalWidth)
            		  g2d.drawString(Utils.doubleToString(sum, 1), x, y-1);
              //advancing x to the next bar by adding bar width and padding
              //of both the bars (i.e. RHS padding of the bar just drawn and LHS
              //padding of the new bar).
            	  x = x+barWidth+(int)( (Math.floor(intervalWidth*0.2F))<1 ? 
            			  1:(Math.floor(intervalWidth*0.2F)) );
              //reseting sum for the next bar.
            	  sum=0;
              }

            }
          }
          //RANKING END
          
          //else if class attribute is numeric or not set then draw black bars.
          else {
            intervalWidth =  (this.getWidth()/(float)m_histBarCounts.length);
            
            //same as in the case of nominal class (see inside of if stmt 
            //corresponding to the current else above).
            if(intervalWidth>5)
              barWidth = (int)Math.floor(intervalWidth*0.8F);
            else
              barWidth = 1;
            
            //same as in the case of nominal class (see inside of if stmt 
            //corresponding to the current else above).
            x = x + (int)( (Math.floor(intervalWidth*0.1F))<1 ? 
                           1:(Math.floor(intervalWidth*0.1F)) );
            
            //same as in the case of nominal class
            if( this.getWidth() - (m_histBarCounts.length*barWidth+
                                  (int)( (Math.floor(intervalWidth*0.2F))<1 ? 
                                         1:(Math.floor(intervalWidth*0.2F)) ) * 
                                  m_histBarCounts.length) > 2 ) {
              x += (this.getWidth() -(m_histBarCounts.length*barWidth + 
                                     (int)((Math.floor(intervalWidth*0.2F))<1 ? 
                                           1:(Math.floor(intervalWidth*0.2F)))*
                                     m_histBarCounts.length))/2;
            }
            
            for(int i=0; i<m_histBarCounts.length; i++) {
              //calculating the proportion of the height of the component 
              //compared to the maxValue in our attribute.
              heightRatio = (this.getHeight()-(float)m_fm.getHeight() - 
                             buttonHeight) / m_maxValue;
              y = (int) (this.getHeight()-Math.round(m_histBarCounts[i]*heightRatio));
              g.fillRect(x, y, barWidth, 
                         (int) Math.round(m_histBarCounts[i]*heightRatio));
              //draw the bar count if it's width is smaller than intervalWidth.
              //draw it 1px above to avoid touching the bar.
              if(m_fm.stringWidth(Utils.doubleToString(m_histBarCounts[i], 1)) < 
                                    intervalWidth)
                g.drawString(Utils.doubleToString(m_histBarCounts[i], 1), x, y-1);
              //Advance x to the next bar by adding bar-width and padding
              //of the bars (RHS padding of current bar & LHS padding of next 
              //bar).
              x = x+barWidth+(int)( (Math.floor(intervalWidth*0.2F))<1 ? 
                                     1:(Math.floor(intervalWidth*0.2F)) );
            }
          }
          
        } //<--end if m_as.nominalCount!=null
        //if the current attribute is numeric then draw a histogram.
        else if(m_as.numericStats != null && 
                (m_histBarClassCounts!=null || m_histBarCounts!=null)) {

        	double heightRatio;
          float intervalWidth;
          int x=0, y=0,  barWidth;
          
          //If the class attribute is set and is not numeric then draw coloured 
          //subbars for the histogram bars
          
          if((m_classIndex >=0) && 
                  (m_data.attribute(m_attribIndex).isNominal())) {
                 
                 //There is a padding of 3px on each side of the histogram.
                 barWidth = ((this.getWidth()-6)/m_histBarClassCounts.length)<1 ? 
                            1 : ((this.getWidth()-6)/m_histBarClassCounts.length);
                 
                 //initializing x to start at the start of the 1st bar after padding.
                 x = 3;
                 //Adding appropriate value to x to account for a "centered" 
                 //histogram
                 if( (this.getWidth() - 
                     (x + m_histBarClassCounts.length*barWidth)) > 5 ) {
                   //We take the current value of x (histogram's RHS padding) and add
                   //the barWidths of all the bars to it to us the size of 
                   //our histogram. We subtract that from the width of the panel 
                   //giving us the extra space that would be left if the histogram is
                   //drawn and divide that by 2 to get the midpoint of that extra
                   //space. That space is then added to our x, hence making the 
                   //histogram centered.
                   x += ( this.getWidth() - 
                         (x + m_histBarClassCounts.length*barWidth) ) / 2;
                 }
                 
                 for(int i=0; i<m_histBarClassCounts.length; i++) {
                   if (m_histBarClassCounts[i] != null) {
                     //Calculating height ratio. Leave space of 19 for an axis line at 
                     //the bottom
                     heightRatio = (this.getHeight()-(float)m_fm.getHeight() - 
                         buttonHeight-19) / m_maxValue;
                     y = this.getHeight()-19;
                     //This would hold the count of the bar (sum of sub-bars).
                     double sum = 0;
                     for(int j=0; j<m_histBarClassCounts[i].numValues(); j++) {
                       y = (int) (y-Math.round(m_histBarClassCounts[i].valueSparse(j) * heightRatio));
                       //System.out.println("Filling x:"+x+" y:"+y+" width:"+barWidth+
                       //                   " height:"+
                       //                   (m_histBarClassCounts[i][j]*heightRatio));
                       //selecting the color corresponding to our class
                       g.setColor( (Color)m_colorList.elementAt(m_histBarClassCounts[i].index(j)) );
                       
                       //RANKING BEGIN
                       if(m_data.attribute(m_classIndex).isRanking()){
                    	   g.setColor(Color.BLUE);
                       }
                       //RANKING END
                       
                       //drawing the bar if its width is greater than 1
                       if(barWidth>1)
                         g.fillRect(x, y, 
                             barWidth, 
                             (int) Math.round(m_histBarClassCounts[i].valueSparse(j)*heightRatio));
                       //otherwise drawing a line
                       else if((m_histBarClassCounts[i].valueSparse(j) * heightRatio)>0)
                         g.drawLine(x, y, x, 
                             (int) (y+Math.round(m_histBarClassCounts[i].valueSparse(j)*heightRatio)));
                       g.setColor(Color.black);
                       
                       //RANKING BEGIN
                       if(m_data.attribute(m_classIndex).isRanking()){
                    	   g.setColor(Color.BLUE);
                       }
                       //RANKING END
                       
                       sum = sum + m_histBarClassCounts[i].valueSparse(j);
                     }
                     //Drawing bar count on the top of the bar if it is < barWidth
                     if(m_fm.stringWidth(" "+Utils.doubleToString(sum, 1))<barWidth)
                       g.drawString(" "+Utils.doubleToString(sum, 1), x, y-1);
                     //Moving x to the next bar
                     x = x+barWidth;
                   }
                 }
                 
                 //Now drawing the axis line at the bottom of the histogram
                 //initializing x again to the start of the plot
                 x = 3;
                 if( (this.getWidth() - 
                     (x + m_histBarClassCounts.length*barWidth)) > 5 )
                   x += (this.getWidth() - 
                        (x + m_histBarClassCounts.length*barWidth))/2;
                 
                 g.drawLine(x, this.getHeight()-17,
                            (barWidth==1)?x+barWidth*m_histBarClassCounts.length-1 : 
                                          x+barWidth*m_histBarClassCounts.length,
                            this.getHeight()-17); //axis line -- see footnote 2.
                 g.drawLine(x, this.getHeight()-16, 
                            x, this.getHeight()-12); //minimum line
                 g.drawString(Utils.doubleToString(m_as.numericStats.min, 2),
                              x,
                              this.getHeight()-12+m_fm.getHeight()); //minimum value
                 g.drawLine(x+(barWidth*m_histBarClassCounts.length)/2,
                            this.getHeight()-16,
                            x+(barWidth*m_histBarClassCounts.length)/2,
                            this.getHeight()-12); //median line
                 //Drawing median value. X position for drawing the value is: from 
                 //start of the plot take the mid point and subtract from it half
                 //of the width of the value to draw.
                 g.drawString(Utils.doubleToString(m_as.numericStats.max/2+m_as.numericStats.min/2, 2),
                              x+(barWidth*m_histBarClassCounts.length)/2 - 
                                m_fm.stringWidth(Utils.doubleToString(m_as.numericStats.max/2+m_as.numericStats.min/2, 2))/2,
                              this.getHeight()-12+m_fm.getHeight()); //median value
                 g.drawLine((barWidth==1) ? x+barWidth*m_histBarClassCounts.length-1:
                                            x+barWidth*m_histBarClassCounts.length,
                            this.getHeight()-16,
                            (barWidth==1) ? x+barWidth*m_histBarClassCounts.length-1:
                                            x+barWidth*m_histBarClassCounts.length,
                            this.getHeight()-12); //maximum line
                 g.drawString(Utils.doubleToString(m_as.numericStats.max, 2),
                              (barWidth==1) ?
                   x+barWidth*m_histBarClassCounts.length-m_fm.stringWidth(Utils.doubleToString(m_as.numericStats.max, 2))-1:
                   x+barWidth*m_histBarClassCounts.length-m_fm.stringWidth(Utils.doubleToString(m_as.numericStats.max, 2)),
                   this.getHeight()-12+m_fm.getHeight()); //maximum value -- see 2.
               }
          
          //RANKING BEGIN
          //modified for usage with ranking attributes.
          if((m_classIndex >=0) && 
             (/*m_data.attribute(m_classIndex).isNominal() ||*/ m_data.attribute(m_classIndex).isRanking())) {
            
            //There is a padding of 3px on each side of the histogram.
            barWidth = ((this.getWidth()-6)/m_histBarClassCounts.length)<1 ? 
                       1 : ((this.getWidth()-6)/m_histBarClassCounts.length);
            
            //initializing x to start at the start of the 1st bar after padding.
            x = 3;
            //Adding appropriate value to x to account for a "centered" 
            //histogram
            if( (this.getWidth() - 
                (x + m_histBarClassCounts.length*barWidth)) > 5 ) {
              //We take the current value of x (histogram's RHS padding) and add
              //the barWidths of all the bars to it to us the size of 
              //our histogram. We subtract that from the width of the panel 
              //giving us the extra space that would be left if the histogram is
              //drawn and divide that by 2 to get the midpoint of that extra
              //space. That space is then added to our x, hence making the 
              //histogram centered.
              x += ( this.getWidth() - 
                    (x + m_histBarClassCounts.length*barWidth) ) / 2;
            }
            
            for(int i=0; i<m_histBarClassCounts.length; i++) {
              if (m_histBarClassCounts[i] != null) {
                //Calculating height ratio. Leave space of 19 for an axis line at 
                //the bottom
                heightRatio = (this.getHeight()-(float)m_fm.getHeight() - 
                    buttonHeight-19) / m_maxValue;
                y = this.getHeight()-19;
                //This would hold the count of the bar (sum of sub-bars).
                double sum = 0;
                for(int j=0; j<m_histBarClassCounts[i].numValues(); j++) {
                  y = (int) (y-Math.round(m_histBarClassCounts[i].valueSparse(j) * heightRatio));
                  //System.out.println("Filling x:"+x+" y:"+y+" width:"+barWidth+
                  //                   " height:"+
                  //                   (m_histBarClassCounts[i][j]*heightRatio));
                  //selecting the color corresponding to our class
     
                  if(m_data.attribute(m_classIndex) instanceof PreferenceAttribute || m_data.get(0) instanceof PreferenceDenseInstance){
                	  g.setColor(Color.blue);
                  }
                  else{
                	  g.setColor( (Color)m_colorList.elementAt(m_histBarClassCounts[i].index(j)) );
                  }
        
                  //drawing the bar if its width is greater than 1
                  if(barWidth>1)
                    g.fillRect(x, y, 
                        barWidth, 
                        (int) Math.round(m_histBarClassCounts[i].valueSparse(j)*heightRatio));
                  //otherwise drawing a line
                  else if((m_histBarClassCounts[i].valueSparse(j) * heightRatio)>0)
                    g.drawLine(x, y, x, 
                        (int) (y+Math.round(m_histBarClassCounts[i].valueSparse(j)*heightRatio)));
                  g.setColor(Color.black);
                  sum = sum + m_histBarClassCounts[i].valueSparse(j);
                }
                //Drawing bar count on the top of the bar if it is < barWidth
                if(m_fm.stringWidth(" "+Utils.doubleToString(sum, 1))<barWidth)
                  g.drawString(" "+Utils.doubleToString(sum, 1), x, y-1);
                //Moving x to the next bar
                x = x+barWidth;
              }
            }
            
            //Now drawing the axis line at the bottom of the histogram
            //initializing x again to the start of the plot
            x = 3;
            if( (this.getWidth() - 
                (x + m_histBarClassCounts.length*barWidth)) > 5 )
              x += (this.getWidth() - 
                   (x + m_histBarClassCounts.length*barWidth))/2;
            
            g.drawLine(x, this.getHeight()-17,
                       (barWidth==1)?x+barWidth*m_histBarClassCounts.length-1 : 
                                     x+barWidth*m_histBarClassCounts.length,
                       this.getHeight()-17); //axis line -- see footnote 2.
            g.drawLine(x, this.getHeight()-16, 
                       x, this.getHeight()-12); //minimum line
            g.drawString(Utils.doubleToString(m_as.numericStats.min, 2),
                         x,
                         this.getHeight()-12+m_fm.getHeight()); //minimum value
            g.drawLine(x+(barWidth*m_histBarClassCounts.length)/2,
                       this.getHeight()-16,
                       x+(barWidth*m_histBarClassCounts.length)/2,
                       this.getHeight()-12); //median line
            //Drawing median value. X position for drawing the value is: from 
            //start of the plot take the mid point and subtract from it half
            //of the width of the value to draw.
            g.drawString(Utils.doubleToString(m_as.numericStats.max/2+m_as.numericStats.min/2, 2),
                         x+(barWidth*m_histBarClassCounts.length)/2 - 
                           m_fm.stringWidth(Utils.doubleToString(m_as.numericStats.max/2+m_as.numericStats.min/2, 2))/2,
                         this.getHeight()-12+m_fm.getHeight()); //median value
            g.drawLine((barWidth==1) ? x+barWidth*m_histBarClassCounts.length-1:
                                       x+barWidth*m_histBarClassCounts.length,
                       this.getHeight()-16,
                       (barWidth==1) ? x+barWidth*m_histBarClassCounts.length-1:
                                       x+barWidth*m_histBarClassCounts.length,
                       this.getHeight()-12); //maximum line
            g.drawString(Utils.doubleToString(m_as.numericStats.max, 2),
                         (barWidth==1) ?
              x+barWidth*m_histBarClassCounts.length-m_fm.stringWidth(Utils.doubleToString(m_as.numericStats.max, 2))-1:
              x+barWidth*m_histBarClassCounts.length-m_fm.stringWidth(Utils.doubleToString(m_as.numericStats.max, 2)),
              this.getHeight()-12+m_fm.getHeight()); //maximum value -- see 2.
          }//RANKING END
          else {  //if class attribute is numeric
            //There is a padding of 3px on each side of the histogram.
            barWidth = ((this.getWidth()-6)/m_histBarCounts.length) < 1 ? 
                        1:((this.getWidth()-6)/m_histBarCounts.length);

            //Same as above. Pls inside of the if stmt.
            x = 3;
            if( (this.getWidth() - (x + m_histBarCounts.length*barWidth)) > 5 )
              x += (this.getWidth() - (x + m_histBarCounts.length*barWidth))/2;
            
            //Same as above
            for(int i=0; i<m_histBarCounts.length; i++) {
              //calculating the ration of the component's height compared to 
              //the maxValue in our current attribute. Leaving 19 pixels to
              //draw the axis at the bottom of the histogram.
              heightRatio = (this.getHeight()-(float)m_fm.getHeight() - 
                             buttonHeight-19) / m_maxValue;
              y = (int) (this.getHeight() - 
                  Math.round(m_histBarCounts[i]*heightRatio)-19);
              //System.out.println("Filling x:"+x+" y:"+y+" width:"+barWidth+
              //                   " height:"+(m_histBarCounts[i]*heightRatio));
              //same as in the if stmt above
              if(barWidth>1)
                g.drawRect(x, y, barWidth, 
                           (int) Math.round(m_histBarCounts[i]*heightRatio));
              else if((m_histBarCounts[i]*heightRatio)>0)
                g.drawLine(x, y, 
                           x, (int) (y+Math.round(m_histBarCounts[i]*heightRatio)));
              if(m_fm.stringWidth(" "+Utils.doubleToString(m_histBarCounts[i], 1)) < 
                    barWidth)
                g.drawString(" "+Utils.doubleToString(m_histBarCounts[i], 1), x, y-1);
              
              x = x+barWidth;
            }
            
            //Now drawing the axis at the bottom of the histogram
            x = 3;
            if( (this.getWidth() - (x + m_histBarCounts.length*barWidth)) > 5 )
              x += (this.getWidth() - (x + m_histBarCounts.length*barWidth))/2;
            
            //This is exact the same as in the if stmt above. See the inside of
            //the stmt for details
            g.drawLine(x, this.getHeight()-17,
                       (barWidth==1) ? x+barWidth*m_histBarCounts.length-1 : 
                                       x+barWidth*m_histBarCounts.length,
                       this.getHeight()-17); //axis line
            g.drawLine(x, this.getHeight()-16, 
                       x, this.getHeight()-12); //minimum line
            g.drawString(Utils.doubleToString(m_as.numericStats.min, 2),
                         x,
                         this.getHeight()-12+m_fm.getHeight()); //minimum value
            g.drawLine(x+(barWidth*m_histBarCounts.length)/2,
                       this.getHeight()-16,
                       x+(barWidth*m_histBarCounts.length)/2,
                       this.getHeight()-12); //median line
            g.drawString(Utils.doubleToString(m_as.numericStats.max/2+m_as.numericStats.min/2, 2),
                         x+(barWidth*m_histBarCounts.length)/2 - 
                           m_fm.stringWidth(Utils.doubleToString(m_as.numericStats.max/2+m_as.numericStats.min/2, 2))/2,
                         this.getHeight()-12+m_fm.getHeight()); //median value
            g.drawLine((barWidth==1) ? x+barWidth*m_histBarCounts.length-1 : 
                                        x+barWidth*m_histBarCounts.length,
                       this.getHeight()-16,
                       (barWidth==1) ? x+barWidth*m_histBarCounts.length-1 : 
                                       x+barWidth*m_histBarCounts.length,
                       this.getHeight()-12); //maximum line
            g.drawString(Utils.doubleToString(m_as.numericStats.max, 2),
                         (barWidth==1) ? 
              x+barWidth*m_histBarCounts.length-m_fm.stringWidth(Utils.doubleToString(m_as.numericStats.max, 2))-1 : 
              x+barWidth*m_histBarCounts.length-m_fm.stringWidth(Utils.doubleToString(m_as.numericStats.max, 2)),
              this.getHeight()-12+m_fm.getHeight()); //maximum value
          }
          //System.out.println("barWidth:"+barWidth+
          //                   " histBarCount:"+m_histBarCounts.length);
          
        } else {
          g.clearRect(0, 0, this.getWidth(), this.getHeight());
          g.drawString("Attribute is neither numeric nor nominal.",
          this.getWidth()/2 - m_fm.
          stringWidth("Attribute is neither numeric nor nominal.")/2,
          this.getHeight()/2-m_fm.getHeight()/2);
        }
      } //<--end if of calculation thread
      else if (m_displayCurrentAttribute) {   //if still calculation thread is running plot
        g.clearRect(0, 0, this.getWidth(), this.getHeight());
        g.drawString("Calculating. Please Wait...",
        this.getWidth()/2 - m_fm.stringWidth("Calculating. Please Wait...")/2,
        this.getHeight()/2-m_fm.getHeight()/2);
      } else if (!m_displayCurrentAttribute) {
        g.clearRect(0, 0, this.getWidth(), this.getHeight());
        g.drawString("Too many values to display.",
        this.getWidth()/2 - m_fm.stringWidth("Too many values to display.")/2,
        this.getHeight()/2-m_fm.getHeight()/2);
      }
    } //<--end if(m_as==null) this means 
  }
  
  
  /**
   * Main method to test this class from command line
   *
   * @param args The arff file and the index of the attribute to use
   */
  public static void main(String [] args) {
    if(args.length!=3) {
      final JFrame jf = new JFrame("AttribVisualization");
      AttributeVisualizationPanel ap = new AttributeVisualizationPanel();
      try {
        Instances ins = new Instances( new FileReader(args[0]) );
        ap.setInstances(ins);
        System.out.println("Loaded: "+args[0]+
                           "\nRelation: "+ap.m_data.relationName()+
                           "\nAttributes: "+ap.m_data.numAttributes());
        ap.setAttribute( Integer.parseInt(args[1]) );
      }
      catch(Exception ex) { ex.printStackTrace(); System.exit(-1); }
      System.out.println("The attributes are: ");
      for(int i=0; i<ap.m_data.numAttributes(); i++)
        System.out.println(ap.m_data.attribute(i).name());
      
      jf.setSize(500, 300);
      jf.getContentPane().setLayout( new BorderLayout() );
      jf.getContentPane().add(ap, BorderLayout.CENTER );
      jf.setDefaultCloseOperation( jf.EXIT_ON_CLOSE );
      jf.setVisible(true);
    }
    else
      System.out.println("Usage: java AttributeVisualizationPanel"+
                         " [arff file] [index of attribute]");
  }
}


/*
 * t =(int) Math.ceil((float)(
 *              (m_data.instance(k).value(m_attribIndex)-m_as.numericStats.min)
 *                           / barRange));
 * 1. 
 * This equation gives a value between (i-1)+smallfraction and i if the 
 * attribute m_attribIndex for the current instances lies in the ith
 * interval. We then increment the value of our i-1th field of our 
 * histogram/barplot array. 
 * If, for example, barRange=3 then, apart from the 1st 
 * interval, each interval i has values in the range 
 * (minValue+3*i-1, minValue+3*i]. The 1st interval has range 
 * [minValue, minValue+i]. Hence it can be seen in the code we specifically 
 * handle t=0 separately.
 *
 */


/**
 * (barWidth==1)?x+barWidth*m_histBarClassCounts.length-1 : 
 *                                    x+barWidth*m_histBarClassCounts.length
 * 2. 
 * In the case barWidth==1 we subtract 1 otherwise the line become one pixel
 * longer than the actual size of the histogram
 */