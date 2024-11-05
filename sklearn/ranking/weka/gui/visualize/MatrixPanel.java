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
 *    MatrixPanel.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.visualize;

import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instances;
import weka.core.labelranking.PreferenceAttribute;
import weka.gui.ExtensionFileFilter;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Image;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Random;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSlider;
import javax.swing.JSplitPane;
import javax.swing.JTextField;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

/** 
 * This panel displays a plot matrix of the user selected attributes
 * of a given data set. 
 * 
 * The datapoints are coloured using a discrete colouring set if the 
 * user has selected a nominal attribute for colouring. If the user
 * has selected a numeric attribute then the datapoints are coloured
 * using a colour spectrum ranging from blue to red (low values to
 * high). Datapoints missing a class value are displayed in black.
 * 
 * @author Ashraf M. Kibriya (amk14@cs.waikato.ac.nz)
 * @version $Revision: 6613 $
 */
public class MatrixPanel
  extends JPanel{

  /** for serialization */
  private static final long serialVersionUID = -1232642719869188740L;

  /** The that panel contains the actual matrix */
  private final Plot m_plotsPanel;

  /** The panel that displays the legend of the colouring attribute */
  protected final ClassPanel m_cp = new ClassPanel();

  /** The panel that contains all the buttons and tools, i.e. resize, jitter bars and sub-sampling buttons etc
      on the bottom of the panel */
  protected JPanel optionsPanel;

  /** Split pane for splitting the matrix and the buttons and bars */
  protected JSplitPane jp;
  /** The button that updates the display to reflect the changes made by the user. 
      E.g. changed attribute set for the matrix    */
  protected JButton m_updateBt = new JButton("Update");

  /** The button to display a window to select attributes */
  protected JButton m_selAttrib = new JButton("Select Attributes");

  /** The dataset for which this panel will display the plot matrix for  */
  protected Instances m_data=null;

  /** The list for selecting the attributes to display the plot matrix */
  protected JList m_attribList = new JList();

  /** The scroll pane to scrolling the matrix */
  protected final JScrollPane m_js = new JScrollPane();

  /** The combo box to allow user to select the colouring attribute */
  protected JComboBox m_classAttrib = new JComboBox();

  /** The slider to adjust the size of the cells in the matrix  */  
  protected JSlider m_plotSize = new JSlider(50, 200, 100);

  /** The slider to adjust the size of the datapoints  */  
  protected JSlider m_pointSize = new JSlider(1, 10, 1);

  /** The slider to add jitter to the plots */  
  protected JSlider m_jitter = new JSlider(0, 20, 0); 

  /** For adding random jitter */
  private Random rnd = new Random();
    
  /** Array containing precalculated jitter values */
  private int jitterVals[][];
 
  /** This stores the size of the datapoint */
  private int datapointSize=1;

  /** The text area for percentage to resample data */
  protected JTextField m_resamplePercent = new JTextField(5);

  /** The label for resample percentage */
  protected JButton m_resampleBt =  new JButton("SubSample % :");

  /** Random seed for random subsample */
  protected JTextField m_rseed = new JTextField(5);
 
  /** Displays the current size beside the slider bar for cell size */
  private final JLabel m_plotSizeLb = new JLabel("PlotSize: [100]");

  /** Displays the current size beside the slider bar for point size */
  private final JLabel m_pointSizeLb = new JLabel("PointSize: [10]");

  /** This array contains the indices of the attributes currently selected  */
  private int [] m_selectedAttribs;

  /** This contains the index of the currently selected colouring attribute  */
  private int m_classIndex;

  /** This is a local array cache for all the instance values for faster rendering */
  private int [][] m_points;

  /** This is an array cache for the colour of each of the instances depending on the 
      colouring attribute. If the colouring attribute is nominal then it contains the 
      index of the colour in our colour list. Otherwise, for numeric colouring attribute,
      it contains the precalculated red component for each instance's colour */
  private int [] m_pointColors;

  /** Contains true for each attribute value (only the selected attributes+class attribute) 
      that is  missing, for each instance.
      m_missing[i][j] == true if m_selectedAttribs[j] is missing in instance i. 
      m_missing[i][m_missing[].length-1] == true  if class value is missing in instance i. */ 
  private boolean [][] m_missing;

  /** This array contains for the classAttribute: <br>
      m_type[0] = [type of attribute, nominal, string or numeric]<br>
      m_type[1] = [number of discrete values of nominal or string attribute <br>
      or same as m_type[0] for numeric attribute] */
  private int [] m_type;

  /** Stores the maximum size for PlotSize label to keep it's size constant */
  private Dimension m_plotLBSizeD;

  /** Stores the maximum size for PointSize label to keep it's size constant */
  private Dimension m_pointLBSizeD;

  /** Contains discrete colours for colouring for nominal attributes */
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
						   Color.black};

  /** color for the font used in column and row names */
  private final Color fontColor = new Color(98, 101, 156);

  /** font used in column and row names */
  private final java.awt.Font f = new java.awt.Font("Dialog", java.awt.Font.BOLD, 11);

  protected transient Image m_osi = null;
  protected boolean[][] m_plottedCells;
  protected boolean m_regenerateOSI = true;
  protected boolean m_clearOSIPlottedCells;
  protected double m_previousPercent = -1;
  
  protected JCheckBox m_fastScroll = 
    new JCheckBox("Fast scrolling (uses more memory)");

  /** 
   * Constructor
   */
  public MatrixPanel() {
    m_rseed.setText("1");

    /** Setting up GUI **/
    m_selAttrib.addActionListener( new ActionListener() {
	public void actionPerformed(ActionEvent ae) {
	  final JDialog jd = new JDialog((JFrame) MatrixPanel.this.getTopLevelAncestor(), 
					 "Attribute Selection Panel",
					 true);

	  JPanel jp = new JPanel();
	  JScrollPane js = new JScrollPane(m_attribList);
	  JButton okBt = new JButton("OK");
	  JButton cancelBt = new JButton("Cancel");
	  final int [] savedSelection = m_attribList.getSelectedIndices();
					
	  okBt.addActionListener( new ActionListener() {	
	      public void actionPerformed(ActionEvent e) {
		jd.dispose(); }
	    } );

	  cancelBt.addActionListener( new ActionListener() {
	      public void actionPerformed(ActionEvent e) {
		m_attribList.setSelectedIndices(savedSelection);
		jd.dispose();}
	    });
	  jd.addWindowListener( new WindowAdapter() {
	      public void windowClosing(WindowEvent e) {
		m_attribList.setSelectedIndices(savedSelection);
		jd.dispose();}
	    });
	  jp.add(okBt);
	  jp.add(cancelBt);

	  jd.getContentPane().add(js, BorderLayout.CENTER); 
	  jd.getContentPane().add(jp, BorderLayout.SOUTH);

	  if(js.getPreferredSize().width < 200)
	    jd.setSize( 250, 250 );
	  else
	    jd.setSize( (int) js.getPreferredSize().width+10, 250);
					
	  jd.setLocation( m_selAttrib.getLocationOnScreen().x,
			  m_selAttrib.getLocationOnScreen().y-jd.getHeight() );
	  jd.setVisible(true);
	}
      });
      
    m_updateBt.addActionListener( new ActionListener() {
	public void actionPerformed(ActionEvent e) {
	    //m_selectedAttribs = m_attribList.getSelectedIndices();
	  initInternalFields();
					
	  Plot a = m_plotsPanel;
	  a.setCellSize( m_plotSize.getValue() );					
	  Dimension d = new Dimension((m_selectedAttribs.length)*(a.cellSize+a.extpad)+2, 
				      (m_selectedAttribs.length)*(a.cellSize+a.extpad)+2
				     );
	  //System.out.println("Size: "+a.cellSize+" Extpad: "+
	  //		   a.extpad+" selected: "+
	  //		   m_selectedAttribs.length+' '+d); 
	  a.setPreferredSize(d);
	  a.setSize( a.getPreferredSize() );
	  a.setJitter( m_jitter.getValue() );
	  
	  if (m_fastScroll.isSelected() && m_clearOSIPlottedCells) {
	    m_plottedCells = new boolean[m_selectedAttribs.length][m_selectedAttribs.length];
	    m_clearOSIPlottedCells = false;
	  }
	  
	  if (m_regenerateOSI) {
	    m_osi = null;
	  }
	  m_js.revalidate();
	  m_cp.setColours(m_colorList);
	  m_cp.setCindex(m_classIndex);
	  m_regenerateOSI = false;
					
	  repaint();
	}
      });
    m_updateBt.setPreferredSize( m_selAttrib.getPreferredSize() );
    
    m_jitter.addChangeListener(new ChangeListener() {
      public void stateChanged(ChangeEvent ce) {
        if (m_fastScroll.isSelected()) {
          m_clearOSIPlottedCells = true;
        }
      }
    });
      
    m_plotSize.addChangeListener( new ChangeListener() {
	public void stateChanged(ChangeEvent ce) {
	  m_plotSizeLb.setText("PlotSize: ["+m_plotSize.getValue()+"]");
	  m_plotSizeLb.setPreferredSize( m_plotLBSizeD );
	  m_jitter.setMaximum( m_plotSize.getValue()/5 ); //20% of cell Size
	  m_regenerateOSI = true;
	}
      });
 
    m_pointSize.addChangeListener( new ChangeListener() {
	public void stateChanged(ChangeEvent ce) {
	  m_pointSizeLb.setText("PointSize: ["+m_pointSize.getValue()+"]");
	  m_pointSizeLb.setPreferredSize( m_pointLBSizeD );
	  datapointSize = m_pointSize.getValue();
	  if (m_fastScroll.isSelected()) {
	    m_clearOSIPlottedCells = true;
	  }
	}
      });
 
    m_resampleBt.addActionListener( new ActionListener() { 
	public void actionPerformed(ActionEvent e) {	  	  
	  JLabel rseedLb = new JLabel("Random Seed: ");
	  JTextField rseedTxt = m_rseed;
	  JLabel percentLb = new JLabel("Subsample as");
	  JLabel percent2Lb = new JLabel("% of input: ");
	  final JTextField percentTxt = new JTextField(5);
	  percentTxt.setText( m_resamplePercent.getText() );
	  JButton doneBt = new JButton("Done");

	  final JDialog jd = new JDialog((JFrame) MatrixPanel.this.getTopLevelAncestor(), 
					 "Subsample % Panel",
					 true) {
	      private static final long serialVersionUID = -269823533147146296L;
	      
	      public void dispose() { 
		m_resamplePercent.setText(percentTxt.getText());
		super.dispose();
	      } 
	    };
	  jd.setDefaultCloseOperation( JDialog.DISPOSE_ON_CLOSE );
			       
	  doneBt.addActionListener( new ActionListener(){ 
	      public void actionPerformed(ActionEvent ae) {
		jd.dispose();  
	      }
	    });
	  GridBagLayout gbl = new GridBagLayout();
	  GridBagConstraints gbc = new GridBagConstraints();
	  JPanel p1 = new JPanel( gbl );		
	  gbc.anchor = GridBagConstraints.WEST; gbc.fill = GridBagConstraints.HORIZONTAL;
	  gbc.insets = new Insets(0,2,2,2);
	  gbc.gridwidth = GridBagConstraints.RELATIVE;
	  p1.add(rseedLb, gbc); gbc.weightx = 0;
	  gbc.gridwidth = GridBagConstraints.REMAINDER; gbc.weightx=1;
	  p1.add(rseedTxt, gbc);
	  gbc.insets = new Insets(8,2,0,2); gbc.weightx=0;
	  p1.add(percentLb, gbc);
	  gbc.insets = new Insets(0,2,2,2); gbc.gridwidth = GridBagConstraints.RELATIVE;
	  p1.add(percent2Lb, gbc);
	  gbc.gridwidth = GridBagConstraints.REMAINDER; gbc.weightx=1;
	  p1.add(percentTxt, gbc);
	  gbc.insets = new Insets(8,2,2,2);

	  JPanel p3 = new JPanel( gbl );
	  gbc.fill = GridBagConstraints.HORIZONTAL; gbc.gridwidth = GridBagConstraints.REMAINDER;
	  gbc.weightx = 1;  gbc.weighty = 0;
	  p3.add(p1, gbc);
	  gbc.insets = new Insets(8,4,8,4);
	  p3.add(doneBt, gbc);
					   
	  jd.getContentPane().setLayout( new BorderLayout() );
	  jd.getContentPane().add(p3, BorderLayout.NORTH);
	  jd.pack();
	  jd.setLocation( m_resampleBt.getLocationOnScreen().x,
			  m_resampleBt.getLocationOnScreen().y-jd.getHeight() );
	  jd.setVisible(true);
	}		
      });

    optionsPanel = new JPanel( new GridBagLayout() ); //all the rest of the panels are in here.
    final JPanel p2 = new JPanel( new BorderLayout() );  //this has class colour panel
    final JPanel p3 = new JPanel( new GridBagLayout() ); //this has update and select buttons
    final JPanel p4 = new JPanel( new GridBagLayout() ); //this has the slider bars and combobox
    GridBagConstraints gbc = new GridBagConstraints();
     
    m_plotLBSizeD = m_plotSizeLb.getPreferredSize();
    m_pointLBSizeD = m_pointSizeLb.getPreferredSize();
    m_pointSizeLb.setText("PointSize: [1]");
    m_pointSizeLb.setPreferredSize( m_pointLBSizeD );
    m_resampleBt.setPreferredSize( m_selAttrib.getPreferredSize() );

    gbc.fill = GridBagConstraints.HORIZONTAL;
    gbc.anchor = GridBagConstraints.NORTHWEST;
    gbc.insets = new Insets(2,2,2,2);
    p4.add(m_plotSizeLb, gbc);
    gbc.weightx=1; gbc.gridwidth = GridBagConstraints.REMAINDER;
    p4.add(m_plotSize, gbc);
    gbc.weightx=0; gbc.gridwidth = GridBagConstraints.RELATIVE;
    p4.add(m_pointSizeLb, gbc);
    gbc.weightx=1; gbc.gridwidth = GridBagConstraints.REMAINDER;
    p4.add(m_pointSize, gbc);
    gbc.weightx=0; gbc.gridwidth = GridBagConstraints.RELATIVE;
    p4.add( new JLabel("Jitter: "), gbc);
    gbc.weightx=1; gbc.gridwidth = GridBagConstraints.REMAINDER;
    p4.add(m_jitter, gbc);
    p4.add(m_classAttrib, gbc);
      
    gbc.gridwidth = GridBagConstraints.REMAINDER;
    gbc.weightx=1;
    gbc.fill = GridBagConstraints.NONE;
    p3.add(m_fastScroll, gbc);
    p3.add(m_updateBt, gbc);
    p3.add(m_selAttrib, gbc);
    gbc.gridwidth = GridBagConstraints.RELATIVE;
    gbc.weightx = 0;
    gbc.fill = GridBagConstraints.VERTICAL;
    gbc.anchor = GridBagConstraints.WEST;
    p3.add(m_resampleBt, gbc);
    gbc.gridwidth = GridBagConstraints.REMAINDER;
    p3.add(m_resamplePercent, gbc);
    
    p2.setBorder(BorderFactory.createTitledBorder("Class Colour"));
    p2.add(m_cp, BorderLayout.SOUTH);

    gbc.insets = new Insets(8,5,2,5);
    gbc.anchor = GridBagConstraints.SOUTHWEST; gbc.fill = GridBagConstraints.HORIZONTAL; gbc.weightx=1;
    gbc.gridwidth = GridBagConstraints.RELATIVE;
    optionsPanel.add(p4, gbc);
    gbc.gridwidth = GridBagConstraints.REMAINDER;
    optionsPanel.add(p3, gbc);
    optionsPanel.add(p2, gbc);
    
    m_fastScroll.setSelected(false);
    m_fastScroll.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        if (!m_fastScroll.isSelected()) {
          m_osi = null;
        } else {
          m_plottedCells = new boolean[m_selectedAttribs.length][m_selectedAttribs.length];
        }
        MatrixPanel.this.invalidate();
        MatrixPanel.this.repaint();
      }
    });

    this.addComponentListener( new ComponentAdapter() {
	public void componentResized(ComponentEvent cv) {
	  m_js.setMinimumSize( new Dimension(MatrixPanel.this.getWidth(),
					     MatrixPanel.this.getHeight()
					     -optionsPanel.getPreferredSize().height-10));
	  jp.setDividerLocation( MatrixPanel.this.getHeight()-optionsPanel.getPreferredSize().height-10 );
	}
      });

    optionsPanel.setMinimumSize( new Dimension(0,0) );
    jp = new JSplitPane(JSplitPane.VERTICAL_SPLIT, m_js, optionsPanel);
    jp.setOneTouchExpandable(true);
    jp.setResizeWeight(1);
    this.setLayout( new BorderLayout() );
    this.add(jp, BorderLayout.CENTER);

    /** Setting up the initial color list **/
    for(int i=0; i<m_defaultColors.length-1; i++)
      m_colorList.addElement(m_defaultColors[i]);
      
    /** Initializing internal fields and components **/
    m_selectedAttribs = m_attribList.getSelectedIndices();
    m_plotsPanel = new Plot();
    m_plotsPanel.setLayout(null);
    m_js.getHorizontalScrollBar().setUnitIncrement( 10 );
    m_js.getVerticalScrollBar().setUnitIncrement( 10 ); 
    m_js.setViewportView( m_plotsPanel );
    m_js.setColumnHeaderView( m_plotsPanel.getColHeader() );
    m_js.setRowHeaderView( m_plotsPanel.getRowHeader() );
    final JLabel lb = new JLabel(" Plot Matrix");
    lb.setFont(f); lb.setForeground(fontColor);
    lb.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
    m_js.setCorner(JScrollPane.UPPER_LEFT_CORNER, lb);
    m_cp.setInstances(m_data);
    m_cp.setBorder(BorderFactory.createEmptyBorder(15,10,10,10));
    m_cp.addRepaintNotify(m_plotsPanel);
    //m_updateBt.doClick(); //not until setting up the instances
  }



  /** Initializes internal data fields, i.e. data values, type, missing and color cache arrays 
   */
  public void initInternalFields() {
    Instances inst = m_data;
    m_classIndex = m_classAttrib.getSelectedIndex();
    m_selectedAttribs = m_attribList.getSelectedIndices();
    double minC=0, maxC=0;

    /** Resampling  **/
    double currentPercent = Double.parseDouble(m_resamplePercent.getText());
    if(currentPercent <= 100) {
      if (currentPercent != m_previousPercent) {
        m_clearOSIPlottedCells = true;
      }
        inst = new Instances(m_data, 0, m_data.numInstances());
        inst.randomize( new Random(Integer.parseInt(m_rseed.getText())) );
        
        //System.err.println("gettingPercent: " +
        //                   Math.round(
        //                     Double.parseDouble(m_resamplePercent.getText())
        //                     / 100D * m_data.numInstances()
        //                             )
        //                  );
        
        inst = new Instances(inst, 
                 0,
                 (int)Math.round(currentPercent
                 / 100D*inst.numInstances())
                            );
        m_previousPercent = currentPercent;
    }
    m_points = new int[inst.numInstances()][m_selectedAttribs.length]; //changed
    m_pointColors = new int[inst.numInstances()];
    m_missing = new boolean[inst.numInstances()][m_selectedAttribs.length+1]; //changed
    m_type = new int[2]; //[m_selectedAttribs.length]; //changed
    jitterVals = new int[inst.numInstances()][2];
      
    /** Setting up the color list for non-numeric attribute as well as jittervals**/
    if(!(inst.attribute(m_classIndex).isNumeric())) {
	  
      for(int i=m_colorList.size(); i<inst.attribute(m_classIndex).numValues()+1; i++) {
	Color pc = m_defaultColors[i % 10];
	int ija =  i / 10;
	ija *= 2; 
	for (int j=0;j<ija;j++) {
	    pc = pc.darker();
	}
	m_colorList.addElement(pc);
      }
	  
      for(int i=0; i<inst.numInstances(); i++) {
	//set to black for missing class value which is last colour is default list
	if(inst.instance(i).isMissing(m_classIndex))
	  m_pointColors[i] =  m_defaultColors.length-1;
	else
	  m_pointColors[i] = (int) inst.instance(i).value(m_classIndex);

	jitterVals[i][0] = rnd.nextInt(m_jitter.getValue()+1)
	  - m_jitter.getValue()/2;
	jitterVals[i][1] = rnd.nextInt(m_jitter.getValue()+1)
	  - m_jitter.getValue()/2;
	      
      }
    }
    /** Setting up color variations for numeric attribute as well as jittervals **/
    else {
      for(int i=0; i<inst.numInstances(); i++) {
	if(!(inst.instance(i).isMissing(m_classIndex))) {
	  minC = maxC = inst.instance(i).value(m_classIndex);
	  break;
	}
      }
	  
      for(int i=1; i<inst.numInstances(); i++) {
	if(!(inst.instance(i).isMissing(m_classIndex))) {
	  if(minC > inst.instance(i).value(m_classIndex))
	    minC = inst.instance(i).value(m_classIndex);
	  if(maxC < inst.instance(i).value(m_classIndex))
	    maxC = inst.instance(i).value(m_classIndex);
	}
      }
	  
      for(int i=0; i<inst.numInstances(); i++) {
	double r = (inst.instance(i).value(m_classIndex) - minC) / (maxC - minC);
	r = (r * 240) + 15;
	m_pointColors[i] = (int)r;

	jitterVals[i][0] = rnd.nextInt(m_jitter.getValue()+1)
	  - m_jitter.getValue()/2;
	jitterVals[i][1] = rnd.nextInt(m_jitter.getValue()+1)
	  - m_jitter.getValue()/2;
      }
    }

    /** Creating local cache of the data values **/
    double min[]=new double[m_selectedAttribs.length], max=0;  //changed
    double ratio[] = new double[m_selectedAttribs.length];     //changed
    double cellSize = m_plotSize.getValue(), temp1=0, temp2=0;

    for(int j=0; j<m_selectedAttribs.length; j++) {
      int i;
      for(i=0; i<inst.numInstances(); i++) {
	min[j] = max = 0;
	if(!(inst.instance(i).isMissing(m_selectedAttribs[j]))) {
	  min[j] = max = inst.instance(i).value(m_selectedAttribs[j]);
	  break;
	}
      }
      for( i=i; i<inst.numInstances(); i++ ) {
	if(!(inst.instance(i).isMissing(m_selectedAttribs[j]))) {
	  if(inst.instance(i).value(m_selectedAttribs[j]) < min[j])
	    min[j] = inst.instance(i).value(m_selectedAttribs[j]);
	  if(inst.instance(i).value(m_selectedAttribs[j]) > max)
	    max = inst.instance(i).value(m_selectedAttribs[j]);
	}
      }
      ratio[j] =  cellSize / (max - min[j]);
    }

    boolean classIndexProcessed=false;
    for(int j=0; j<m_selectedAttribs.length; j++) {
      if(inst.attribute(m_selectedAttribs[j]).isNominal() || inst.attribute(m_selectedAttribs[j]).isRanking() || inst.attribute(m_selectedAttribs[j]).isString()) {
	  //m_type[0][j] = 1;  m_type[1][j] = inst.attribute(m_selectedAttribs[j]).numValues();

	temp1 = cellSize/(double)inst.attribute(m_selectedAttribs[j]).numValues(); //m_type[1][j];
	temp2 = temp1/2;
	for(int i=0; i<inst.numInstances(); i++) {
	  m_points[i][j] = (int) Math.round(temp2+temp1*inst.instance(i).value(m_selectedAttribs[j]));
	  if(inst.instance(i).isMissing(m_selectedAttribs[j])) {
	    m_missing[i][j] = true;    //represents missing value
	    if(m_selectedAttribs[j]==m_classIndex) {
		m_missing[i][m_missing[0].length-1] = true;
		classIndexProcessed = true;
	    }
	  }
	}
      }
      else {
	  //m_type[0][j] = m_type[1][j] = 0;
	for(int i=0; i<inst.numInstances(); i++) {
	  m_points[i][j] = (int) Math.round((inst.instance(i).value(m_selectedAttribs[j])
					     -min[j])*ratio[j]);	
	  if(inst.instance(i).isMissing(m_selectedAttribs[j])) {
	    m_missing[i][j] = true;    //represents missing value
	    if(m_selectedAttribs[j]==m_classIndex) {
		m_missing[i][m_missing[0].length-1] = true;
		classIndexProcessed = true;
	    }
	  }
	}
      }
    }

    if(inst.attribute(m_classIndex).isNominal()|| inst.attribute(m_classIndex).isRanking() || inst.attribute(m_classIndex).isString()) {
	m_type[0] = 1; m_type[1] = inst.attribute(m_classIndex).numValues();
    }
    else
	m_type[0] = m_type[1] = 0;

    if(classIndexProcessed==false) {  //class Index has not been processed as class index is not among the selected attribs
	for(int i=0; i<inst.numInstances(); i++) {
	    if(inst.instance(i).isMissing(m_classIndex))
		m_missing[i][m_missing[0].length-1] = true;
	}
    }

    m_cp.setColours(m_colorList);
  }

  /** Sets up the UI's attributes lists 
   */  
  public void setupAttribLists() {
    String [] tempAttribNames = new String[m_data.numAttributes()];
    String type;

    m_classAttrib.removeAllItems();
    for(int i=0; i<tempAttribNames.length; i++) {
      switch (m_data.attribute(i).type()) {
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
      default:
	type = " (???)";
      }
      tempAttribNames[i] = new String("Colour: "+m_data.attribute(i).name()+" "+type);
      m_classAttrib.addItem(tempAttribNames[i]);
    }
    if (m_data.classIndex() == -1)
      m_classAttrib.setSelectedIndex(tempAttribNames.length - 1);
    else
      m_classAttrib.setSelectedIndex(m_data.classIndex());
    m_attribList.setListData(tempAttribNames);
    m_attribList.setSelectionInterval(0, tempAttribNames.length-1);
  }

  /** Calculates the percentage to resample 
   */
  public void setPercent() {
    if(m_data.numInstances() > 700) {
      double percnt = 500D/m_data.numInstances()*100;     
      percnt *= 100;
      percnt = Math.round(percnt);
      percnt /= 100;

      m_resamplePercent.setText(""+percnt);
    }
    else
      m_resamplePercent.setText("100");
  }


  /** This method changes the Instances object of this class to a new one. It also does all the necessary
      initializations for displaying the panel. This must be called before trying to display the panel.
      @param newInst The new set of Instances
  */
  public void setInstances(Instances newInst) {

    m_osi = null;
    m_fastScroll.setSelected(false);
    m_data = newInst;
    setPercent();
    setupAttribLists();
    m_rseed.setText("1");
    initInternalFields();
    m_cp.setInstances(m_data);
    m_cp.setCindex(m_classIndex);
    m_updateBt.doClick();
  }


  /**
     Main method for testing this class
  */
  public static void main(String [] args)  {
    final JFrame jf = new JFrame("Weka Explorer: MatrixPanel");
    final JButton setBt = new JButton("Set Instances");
    Instances data = null;
    try {
      if(args.length==1)
	data = new Instances( new BufferedReader( new FileReader(args[0])) ); 
      else {
	System.out.println("Usage: MatrixPanel <arff file>"); 
	System.exit(-1);
      }
    } catch(IOException ex) { ex.printStackTrace(); System.exit(-1); }
     
    final MatrixPanel mp = new MatrixPanel();
    mp.setInstances(data);
    setBt.addActionListener( new ActionListener() {
	public void actionPerformed(ActionEvent e) {
	  JFileChooser chooser = new JFileChooser(new java.io.File(System.getProperty("user.dir")));
	  ExtensionFileFilter myfilter = new ExtensionFileFilter("arff", "Arff data files");
	  chooser.setFileFilter(myfilter);
	  int returnVal = chooser.showOpenDialog(jf);
		  
	  if(returnVal == JFileChooser.APPROVE_OPTION)
	    {
	      try{
		System.out.println("You chose to open this file: " +chooser.getSelectedFile().getName());
		Instances in = new Instances ( new FileReader(chooser.getSelectedFile().getAbsolutePath()) );
		mp.setInstances(in);
	      }
	      catch(Exception ex) { ex.printStackTrace(); }
	    }
	}
      });
    //System.out.println("Loaded: "+args[0]+"\nRelation: "+data.relationName()+"\nAttributes: "+data.numAttributes());
    //System.out.println("The attributes are: ");
    //for(int i=0; i<data.numAttributes(); i++)
    //  System.out.println(data.attribute(i).name());

    //RepaintManager.currentManager(jf.getRootPane()).setDoubleBufferingEnabled(false);
    jf.getContentPane().setLayout( new BorderLayout() );
    jf.getContentPane().add(mp, BorderLayout.CENTER);
    jf.getContentPane().add(setBt, BorderLayout.SOUTH);
    jf.getContentPane().setFont( new java.awt.Font( "SansSerif", java.awt.Font.PLAIN, 11) );
    jf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    jf.setSize(800, 600);
    jf.setVisible(true);
    jf.repaint();
  }

  
  /**
     Internal class responsible for displaying the actual matrix
     Requires the internal data fields of the parent class to be properly initialized
     before being created
  */
  private class Plot
    extends JPanel
    implements MouseMotionListener, MouseListener {

    /** for serialization */
    private static final long serialVersionUID = -1721245738439420882L;    

    int extpad=3, intpad=4, cellSize=100, cellRange=100, lastx=0, lasty=0, jitter=0;
    java.awt.Rectangle r;
    java.awt.FontMetrics fm;
    int lastxpos, lastypos;
    JPanel jPlColHeader, jPlRowHeader;


    /** Constructor 
     */
    public Plot() {
      super();
      this.setToolTipText("blah");
      this.addMouseMotionListener( this );
      this.addMouseListener( this );
      initialize();
    }

    /** Initializes the internal fields */
    public void initialize() {
      lastxpos = lastypos = 0;	  
      cellRange = cellSize; cellSize = cellRange + 2*intpad;

      jPlColHeader = new JPanel() {
	private static final long serialVersionUID = -9098547751937467506L;
        java.awt.Rectangle r;
        public void paint(Graphics g) {
          r = g.getClipBounds();
          g.setColor(this.getBackground());
          g.fillRect(r.x, r.y, r.width, r.height);
          g.setFont( f );
          fm = g.getFontMetrics();
          int xpos = 0, ypos = 0, attribWidth=0;
          
          g.setColor(fontColor);
          xpos = extpad;
          ypos=extpad+fm.getHeight();
          
          for(int i=0; i<m_selectedAttribs.length; i++) {
            if( xpos+cellSize < r.x)
            { xpos += cellSize+extpad; continue; }
            else if(xpos > r.x+r.width)
            { break; }
            else {
              attribWidth = fm.stringWidth(m_data.attribute(m_selectedAttribs[i]).name());
              g.drawString(m_data.attribute(m_selectedAttribs[i]).name(),
              (attribWidth<cellSize) ? (xpos + (cellSize/2 - attribWidth/2)):xpos,
              ypos);
            }
            xpos += cellSize+extpad;
          }
          fm = null; r=null;
        }
        
        public Dimension getPreferredSize() {
          fm = this.getFontMetrics(this.getFont());
          return new Dimension( m_selectedAttribs.length*(cellSize+extpad),
          2*extpad + fm.getHeight() );
        }
      };

      jPlRowHeader = new JPanel() {
	private static final long serialVersionUID = 8474957069309552844L;
	
        java.awt.Rectangle r;
        public void paint(Graphics g) {
          r = g.getClipBounds();
          g.setColor(this.getBackground());
          g.fillRect(r.x, r.y, r.width, r.height);
          g.setFont( f );
          fm = g.getFontMetrics();
          int xpos = 0, ypos = 0;
          
          g.setColor(fontColor);
          xpos = extpad;
          ypos=extpad;
          
          for(int j=m_selectedAttribs.length-1; j>=0; j--) {
            if( ypos+cellSize < r.y )
            { ypos += cellSize+extpad;  continue; }
            else if( ypos > r.y+r.height )
              break;
            else {
              g.drawString(m_data.attribute(m_selectedAttribs[j]).name(), xpos+extpad, ypos+cellSize/2);
            }
            xpos = extpad;
            ypos += cellSize+extpad;
          }
          r=null;
        }
        
        public Dimension getPreferredSize() {
          return new Dimension( 100+extpad,
          m_selectedAttribs.length*(cellSize+extpad)
          );
        }
      };
      jPlColHeader.setFont(f);
      jPlRowHeader.setFont(f);
      this.setFont(f);
    }      

    public JPanel getRowHeader() {
	  return jPlRowHeader;
    }

    public JPanel getColHeader() {
	return jPlColHeader;
    }

    public void mouseMoved(MouseEvent e) {
      Graphics g = this.getGraphics();
      int xpos=extpad, ypos=extpad;

      for(int j=m_selectedAttribs.length-1; j>=0; j--) {
	for(int i=0; i<m_selectedAttribs.length; i++) {
	  if(e.getX()>=xpos && e.getX()<=xpos+cellSize+extpad)
	    if(e.getY()>=ypos && e.getY()<=ypos+cellSize+extpad) {
	      if(xpos!=lastxpos || ypos!=lastypos) {
		g.setColor( Color.red );
		g.drawRect(xpos-1, ypos-1, cellSize+1, cellSize+1);
		if(lastxpos!=0 && lastypos!=0) {
		  g.setColor( this.getBackground().darker() );
		  g.drawRect(lastxpos-1, lastypos-1, cellSize+1, cellSize+1); }
		lastxpos = xpos; lastypos = ypos;
	      }
	      return;
	    }
	  xpos+=cellSize+extpad;
	}
	xpos=extpad;
	ypos+=cellSize+extpad;
      }
      if(lastxpos!=0 && lastypos!=0) {
	g.setColor( this.getBackground().darker() );
	g.drawRect(lastxpos-1, lastypos-1, cellSize+1, cellSize+1); }
      lastxpos=lastypos=0;
    }

    public void mouseDragged(MouseEvent e){ }

    public void mouseClicked(MouseEvent e) {
      int i=0, j=0, found=0;
	  
      int xpos=extpad, ypos=extpad;
      for(j=m_selectedAttribs.length-1; j>=0; j--) {
	for(i=0; i<m_selectedAttribs.length; i++) {
	  if(e.getX()>=xpos && e.getX()<=xpos+cellSize+extpad)
	    if(e.getY()>=ypos && e.getY()<=ypos+cellSize+extpad) {
	      found=1; break;
	    }
	  xpos+=cellSize+extpad;
	}
	if(found==1)
	  break;
	xpos=extpad;
	ypos+=cellSize+extpad;
      }
      if(found==0)
	return;

      JFrame jf = new JFrame("Weka Explorer: Visualizing "+m_data.relationName() );
      VisualizePanel vp = new VisualizePanel();
      try {
	PlotData2D pd = new PlotData2D(m_data);
	pd.setPlotName("Master Plot");
	vp.setMasterPlot(pd);
	//System.out.println("x: "+i+" y: "+j);
	vp.setXIndex(m_selectedAttribs[i]);
	vp.setYIndex(m_selectedAttribs[j]);
	vp.m_ColourCombo.setSelectedIndex( m_classIndex );
      }
      catch(Exception ex) { ex.printStackTrace(); }
      jf.getContentPane().add(vp);
      jf.setSize(800,600);
      jf.setVisible(true);
    } 

    public void mouseEntered(MouseEvent e){ }
    public void mouseExited(MouseEvent e){ }
    public void mousePressed(MouseEvent e){ }
    public void mouseReleased(MouseEvent e){ }

    /** sets the new jitter value for the plots
     */
    public void setJitter(int newjitter) {
      jitter = newjitter;
    }
      
    /** sets the new size for the plots
     */
    public void setCellSize(int newCellSize) {
      cellSize = newCellSize;
      initialize();
    }

    /** Returns the X and Y attributes of the plot the mouse is currently
	on
    */
    public String getToolTipText(MouseEvent event) {
      int xpos=extpad, ypos=extpad;
	  
      for(int j=m_selectedAttribs.length-1; j>=0; j--) {
	for(int i=0; i<m_selectedAttribs.length; i++) {
	  if(event.getX()>=xpos && event.getX()<=xpos+cellSize+extpad)
	    if(event.getY()>=ypos && event.getY()<=ypos+cellSize+extpad)
	      return("X: "+m_data.attribute(m_selectedAttribs[i]).name()+
		     " Y: "+m_data.attribute(m_selectedAttribs[j]).name()+
		     " (click to enlarge)");
	  xpos+=cellSize+extpad;
	}
	xpos=extpad;
	ypos+=cellSize+extpad;
      }
      return ("Matrix Panel");
    }
        

    /**  Paints a single Plot at xpos, ypos. and xattrib and yattrib on X and
	 Y axes
    */
    public void paintGraph(Graphics g, int xattrib, int yattrib, int xpos, int ypos) {
      int x, y;
      g.setColor( this.getBackground().darker().darker() );
      g.drawRect(xpos-1, ypos-1, cellSize+1, cellSize+1);
      g.setColor(Color.white);
      g.fillRect(xpos, ypos, cellSize, cellSize);
      for(int i=0; i<m_points.length; i++) {
        
        if( !(m_missing[i][yattrib] || m_missing[i][xattrib]) ) {
          
          if(m_type[0]==0)
            if(m_missing[i][m_missing[0].length-1])
              g.setColor(m_defaultColors[m_defaultColors.length-1]);
            else
              g.setColor( new Color(m_pointColors[i],150,(255-m_pointColors[i])) );
          else
            g.setColor((Color)m_colorList.elementAt(m_pointColors[i]));
          
          if(m_points[i][xattrib]+jitterVals[i][0]<0 || m_points[i][xattrib]+jitterVals[i][0]>cellRange)
            if(cellRange-m_points[i][yattrib]+jitterVals[i][1]<0 || cellRange-m_points[i][yattrib]+jitterVals[i][1]>cellRange) {
              //both x and y out of range don't add jitter
              x=intpad+m_points[i][xattrib];
              y=intpad+(cellRange - m_points[i][yattrib]);
            }
            else {
              //only x out of range
              x=intpad+m_points[i][xattrib];
              y=intpad+(cellRange - m_points[i][yattrib])+jitterVals[i][1];
            }
          else if(cellRange-m_points[i][yattrib]+jitterVals[i][1]<0 || cellRange-m_points[i][yattrib]+jitterVals[i][1]>cellRange) {
            //only y out of range
            x=intpad+m_points[i][xattrib]+jitterVals[i][0];
            y=intpad+(cellRange - m_points[i][yattrib]);
          }
          else {
            //none out of range
            x=intpad+m_points[i][xattrib]+jitterVals[i][0];
            y=intpad+(cellRange - m_points[i][yattrib])+jitterVals[i][1];
          }
          if(datapointSize==1)
            g.drawLine(x+xpos, y+ypos, x+xpos, y+ypos);
          else
            g.drawOval(x+xpos-datapointSize/2, y+ypos-datapointSize/2, datapointSize, datapointSize);
        }
      }
      g.setColor( fontColor );
    }
    
    private void createOSI() {
      int iwidth = this.getWidth();
      int iheight = this.getHeight();
      m_osi = this.createImage(iwidth, iheight);
      clearOSI();
    }
    
    private void clearOSI() {
      if (m_osi == null) {
        return;
      }
      
      int iwidth = this.getWidth();
      int iheight = this.getHeight();
      Graphics m = m_osi.getGraphics();
      m.setColor(this.getBackground().darker().darker());
      m.fillRect(0, 0, iwidth, iheight);
    }
    

    /**
       Paints the matrix of plots in the current visible region
    */
    public void paintME(Graphics g) {
      Graphics g2 = g;
      if (m_osi == null && m_fastScroll.isSelected()) {
        createOSI();
      }
      if (m_osi != null && m_fastScroll.isSelected()) {
        g2 = m_osi.getGraphics();
      }
      r = g.getClipBounds();
      
      g.setColor( this.getBackground() );
      g.fillRect(r.x, r.y, r.width, r.height);
      g.setColor( fontColor );
      
      int xpos = 0, ypos = 0;
      
      xpos = extpad;
      ypos=extpad;
      
      
      for(int j=m_selectedAttribs.length-1; j>=0; j--) {
        if( ypos+cellSize < r.y )
        { ypos += cellSize+extpad;  continue; }
        else if( ypos > r.y+r.height )
          break;
        else {
          for(int i=0; i<m_selectedAttribs.length; i++) {
            if( xpos+cellSize < r.x) {
              xpos += cellSize+extpad; continue; }
            else if(xpos > r.x+r.width)
              break;
            else if (m_fastScroll.isSelected()) {
              if (!m_plottedCells[i][j]) {
                paintGraph(g2, i, j, xpos, ypos); //m_selectedAttribs[i], m_selectedAttribs[j], xpos, ypos);
                m_plottedCells[i][j] = true;
              }
            } else {
              paintGraph(g2, i, j, xpos, ypos);
            }
            xpos += cellSize+extpad;
          }
        }
        xpos = extpad;
        ypos += cellSize+extpad;
      }
    }
      
    /** paints this JPanel (PlotsPanel)
     */
    public void paintComponent(Graphics g) {
      paintME(g);
      if (m_osi != null && m_fastScroll.isSelected()) {
        g.drawImage(m_osi, 0, 0, this);
      }
    }
  }
}
