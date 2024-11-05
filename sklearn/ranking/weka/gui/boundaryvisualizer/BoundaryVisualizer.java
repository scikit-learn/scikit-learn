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
 *   BoundaryVisualizer.java
 *   Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.boundaryvisualizer;

import weka.classifiers.Classifier;
import weka.classifiers.AbstractClassifier;
import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instances;
import weka.core.DenseInstance;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformationHandler;
import weka.core.Utils;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;
import weka.core.labelranking.PreferenceAttribute;
import weka.gui.ExtensionFileFilter;
import weka.gui.GenericObjectEditor;
import weka.gui.PropertyPanel;
import weka.gui.visualize.ClassPanel;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.io.FileOutputStream;
import java.io.ObjectOutputStream;
import java.util.Vector;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTextField;

/**
 * BoundaryVisualizer. Allows the visualization of classifier decision
 * boundaries in two dimensions. A supplied classifier is first
 * trained on supplied training data, then a data generator (currently
 * using kernels) is used to generate new instances at points fixed in
 * the two visualization dimensions but random in the other
 * dimensions. These instances are classified by the classifier and
 * plotted as points with colour corresponding to the probability
 * distribution predicted by the classifier. At present, 2 * 2^(#
 * non-fixed dimensions) points are generated from each kernel per
 * pixel in the display. In practice, fewer points than this are
 * actually classified because kernels are weighted (on a per-pixel
 * basis) according to the fixexd dimensions and kernels corresponding
 * to the lowest 1% of the weight mass are discarded. Predicted
 * probability distributions are weighted (acording to the fixed
 * visualization dimensions) and averaged to produce an RGB value for
 * the pixel. For more information, see<p>
 * 
 * Eibe Frank and Mark Hall (2003). Visualizing Class Probability
 * Estimators. Working Paper 02/03, Department of Computer Science,
 * University of Waikato.
 *
 * @author <a href="mailto:mhall@cs.waikato.ac.nz">Mark Hall</a>
 * @version $Revision: 6482 $
 * @since 1.0
 * @see JPanel 
 */
public class BoundaryVisualizer
  extends JPanel implements TechnicalInformationHandler {

  /** for serialization */
  private static final long serialVersionUID = 3933877580074013208L;

  /**
   * Inner class to handle rendering the axis
   *
   * @author <a href="mailto:mhall@cs.waikato.ac.nz">Mark Hall</a>
   * @version $Revision: 6482 $
   * @since 1.0
   * @see JPanel
   */
  private class AxisPanel
    extends JPanel {

    /** for serialization */
    private static final long serialVersionUID = -7421022416674492712L;
    
    private static final int MAX_PRECISION = 10;
    private boolean m_vertical = false;
    private final int PAD = 5;
    private FontMetrics m_fontMetrics;
    private int m_fontHeight;
    
    public AxisPanel(boolean vertical) {
      m_vertical = vertical;
      this.setBackground(Color.black);
      //      Graphics g = this.getGraphics();
      String fontFamily = this.getFont().getFamily();
      Font newFont = new Font(fontFamily, Font.PLAIN, 10);
      this.setFont(newFont);
    }

    public Dimension getPreferredSize() {
      if (m_fontMetrics == null) {
	Graphics g = this.getGraphics();
	m_fontMetrics = g.getFontMetrics();
	m_fontHeight = m_fontMetrics.getHeight();
      }
      if (!m_vertical) {
	return new Dimension(this.getSize().width, PAD+2+m_fontHeight);
      }
      return new Dimension(50, this.getSize().height);
    }

    public void paintComponent(Graphics g) {
      super.paintComponent(g);
      this.setBackground(Color.black);
      if (m_fontMetrics == null) {
	m_fontMetrics = g.getFontMetrics();
	m_fontHeight = m_fontMetrics.getHeight();
      }

      Dimension d = this.getSize();
      Dimension d2 = m_boundaryPanel.getSize();
      g.setColor(Color.gray);
      int hf = m_fontMetrics.getAscent();
      if (!m_vertical) {
	g.drawLine(d.width, PAD, d.width-d2.width, PAD);
	// try and draw some scale values
	if (getInstances() != null) {
	  int precisionXmax = 1;
	  int precisionXmin = 1;
	  int whole = (int)Math.abs(m_maxX);
	  double decimal = Math.abs(m_maxX) - whole;
	  int nondecimal;
	  nondecimal = (whole > 0) 
	    ? (int)(Math.log(whole) / Math.log(10))
	    : 1;
	  
	  precisionXmax = (decimal > 0) 
	    ? (int)Math.abs(((Math.log(Math.abs(m_maxX)) / 
			      Math.log(10))))+2
	    : 1;
	  if (precisionXmax > MAX_PRECISION) {
	    precisionXmax = 1;
	  }
	  String maxStringX = Utils.doubleToString(m_maxX,
						   nondecimal+1+precisionXmax
						   ,precisionXmax);
	  
	  whole = (int)Math.abs(m_minX);
	  decimal = Math.abs(m_minX) - whole;
	  nondecimal = (whole > 0) 
	    ? (int)(Math.log(whole) / Math.log(10))
	    : 1;
	  precisionXmin = (decimal > 0) 
	    ? (int)Math.abs(((Math.log(Math.abs(m_minX)) / 
			      Math.log(10))))+2
	    : 1;
	  if (precisionXmin > MAX_PRECISION) {
	    precisionXmin = 1;
	  }
	  
	  String minStringX = Utils.doubleToString(m_minX,
						   nondecimal+1+precisionXmin,
						   precisionXmin);
	  g.drawString(minStringX,  d.width-d2.width, PAD+hf+2);
	  int maxWidth = m_fontMetrics.stringWidth(maxStringX);
	  g.drawString(maxStringX, d.width-maxWidth, PAD+hf+2);
	}
      } else {
	g.drawLine(d.width-PAD, 0, d.width-PAD, d2.height);
	// try and draw some scale values
	if (getInstances() != null) {
	  int precisionYmax = 1;
	  int precisionYmin = 1;
	  int whole = (int)Math.abs(m_maxY);
	  double decimal = Math.abs(m_maxY) - whole;
	  int nondecimal;
	  nondecimal = (whole > 0) 
	    ? (int)(Math.log(whole) / Math.log(10))
	    : 1;
	  
	  precisionYmax = (decimal > 0) 
	    ? (int)Math.abs(((Math.log(Math.abs(m_maxY)) / 
			      Math.log(10))))+2
	    : 1;
	  if (precisionYmax > MAX_PRECISION) {
	    precisionYmax = 1;
	  }
	  String maxStringY = Utils.doubleToString(m_maxY,
						   nondecimal+1+precisionYmax
						   ,precisionYmax);
	  
	  whole = (int)Math.abs(m_minY);
	  decimal = Math.abs(m_minY) - whole;
	  nondecimal = (whole > 0) 
	    ? (int)(Math.log(whole) / Math.log(10))
	    : 1;
	  precisionYmin = (decimal > 0) 
	    ? (int)Math.abs(((Math.log(Math.abs(m_minY)) / 
			      Math.log(10))))+2
	    : 1;
	  if (precisionYmin > MAX_PRECISION) {
	    precisionYmin = 1;
	  }
	  
	  String minStringY = Utils.doubleToString(m_minY,
						   nondecimal+1+precisionYmin,
						   precisionYmin);
	  int maxWidth = m_fontMetrics.stringWidth(minStringY);
	  g.drawString(minStringY,  d.width-PAD-maxWidth-2, d2.height);
	  maxWidth = m_fontMetrics.stringWidth(maxStringY);
	  g.drawString(maxStringY, d.width-PAD-maxWidth-2, hf);
	}
      }
    }
  }
  
  /** the number of visualizer windows we have open. */
  protected static int m_WindowCount = 0; 
  
  /** whether the exit if there are no more windows open */
  protected static boolean m_ExitIfNoWindowsOpen = true;

  /** the training instances */
  private Instances m_trainingInstances;

  /** the classifier to use */
  private Classifier m_classifier;

  // plot area dimensions
  protected int m_plotAreaWidth = 384;
  //protected int m_plotAreaHeight = 384;
  protected int m_plotAreaHeight = 384;

  /** the plotting panel */
  protected BoundaryPanel m_boundaryPanel;

  // combo boxes for selecting the class attribute, class values (for
  // colouring pixels), and visualization attributes
  protected JComboBox m_classAttBox = new JComboBox();
  protected JComboBox m_xAttBox = new JComboBox();
  protected JComboBox m_yAttBox = new JComboBox();

  protected Dimension COMBO_SIZE = 
    new Dimension((int)(m_plotAreaWidth * 0.75),
		  m_classAttBox.getPreferredSize().height);

  protected JButton m_startBut = new JButton("Start");

  protected JCheckBox m_plotTrainingData = new JCheckBox("Plot training data");

  protected JPanel m_controlPanel;

  protected ClassPanel m_classPanel = new ClassPanel();

  // separate panels for rendering axis information
  private AxisPanel m_xAxisPanel;
  private AxisPanel m_yAxisPanel;

  // min and max values for visualization dimensions
  private double m_maxX;
  private double m_maxY;
  private double m_minX;
  private double m_minY;

  private int m_xIndex;
  private int m_yIndex;

  /* Kernel density estimator/generator */
  private KDDataGenerator m_dataGenerator;

  /* number of samples per pixel (fixed dimensions only) */
  private int m_numberOfSamplesFromEachRegion;

  /** base for sampling in the non-fixed dimensions */
  private int m_generatorSamplesBase;

  /** Set the kernel bandwidth to cover this many nearest neighbours */
  private int m_kernelBandwidth;
  
  private JTextField m_regionSamplesText = 
    new JTextField(""+0);

  private JTextField m_generatorSamplesText = 
    new JTextField(""+0);

  private JTextField m_kernelBandwidthText = 
    new JTextField(""+3+"  ");

  //jimmy
  protected GenericObjectEditor m_classifierEditor = new GenericObjectEditor(); //the widget to select the classifier
  protected PropertyPanel m_ClassifierPanel = new PropertyPanel(m_classifierEditor);
  /** The file chooser for selecting arff files */
  protected JFileChooser m_FileChooser 
    = new JFileChooser(new File(System.getProperty("user.dir")));
  protected ExtensionFileFilter m_arffFileFilter = 
    new ExtensionFileFilter(Instances.FILE_EXTENSION,
			    "Arff data files");
  protected JLabel dataFileLabel = new JLabel(); //stores the name of the data file (currently stores relation name rather than filename)
  protected JPanel m_addRemovePointsPanel = new JPanel(); //a panel which contains the controls to add and remove points
  protected JComboBox m_classValueSelector = new JComboBox(); //a widget to select the class attribute.
  protected JRadioButton m_addPointsButton = new JRadioButton(); //when this is selected, clicking on the BoundaryPanel will add points.
  protected JRadioButton m_removePointsButton = new JRadioButton(); //when this is selected, clicking on the BoundaryPanel will remove points.
  protected ButtonGroup m_addRemovePointsButtonGroup = new ButtonGroup();
  protected JButton removeAllButton = new JButton ("Remove all"); //button to remove all points
  protected JButton chooseButton = new JButton("Open File"); //button to choose a data file
  
  /* Register the property editors we need */
  static {
    GenericObjectEditor.registerEditors();
  }
  
  /**
   * Returns a string describing this tool
   * @return a description of the tool suitable for
   * displaying in various Weka GUIs
   */
  public String globalInfo() {
    return "Class for visualizing class probability estimates.\n\n"
    + "For more information, see\n\n"
    + getTechnicalInformation().toString();
  }
  
  /**
   * Returns an instance of a TechnicalInformation object, containing 
   * detailed information about the technical background of this class,
   * e.g., paper reference or book this class is based on.
   * 
   * @return the technical information about this class
   */
  public TechnicalInformation getTechnicalInformation() {
    TechnicalInformation        result;

    result = new TechnicalInformation(Type.INPROCEEDINGS);
    result.setValue(Field.AUTHOR, "Eibe Frank and Mark Hall");
    result.setValue(Field.TITLE, "Visualizing class probability estimators");
    result.setValue(Field.BOOKTITLE, "European Conference on Principles and Practice of " +
    		"Knowledge Discovery in Databases");
    result.setValue(Field.YEAR, "2003");
    result.setValue(Field.PAGES, "168-169");
    result.setValue(Field.PUBLISHER, "Springer-Verlag");
    result.setValue(Field.ADDRESS, "Cavtat-Dubrovnik");

    return result;
  }


  /**
   * Creates a new <code>BoundaryVisualizer</code> instance.
   */
  public BoundaryVisualizer() {
    
    setLayout(new BorderLayout());
    m_classAttBox.setMinimumSize(COMBO_SIZE);
    m_classAttBox.setPreferredSize(COMBO_SIZE);
    m_classAttBox.setMaximumSize(COMBO_SIZE);
    m_classAttBox.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        if (m_classAttBox.getItemCount() != 0)
	{
		try { 
			m_classPanel.setCindex(m_classAttBox.getSelectedIndex());
			plotTrainingData();
			System.err.println("Here in class att box listener");
		} catch (Exception ex) {ex.printStackTrace();}
		
		//set up the add points selector combo box. -jimmy
		setUpClassValueSelectorCB();
	}
      }
      });
	    

    m_xAttBox.setMinimumSize(COMBO_SIZE);
    m_xAttBox.setPreferredSize(COMBO_SIZE);
    m_xAttBox.setMaximumSize(COMBO_SIZE);

    m_yAttBox.setMinimumSize(COMBO_SIZE);
    m_yAttBox.setPreferredSize(COMBO_SIZE);
    m_yAttBox.setMaximumSize(COMBO_SIZE);

    m_classPanel.setMinimumSize(new 
      Dimension((int)COMBO_SIZE.getWidth()*2, 
		(int)COMBO_SIZE.getHeight()*2));
    m_classPanel.setPreferredSize(new 
      Dimension((int)COMBO_SIZE.getWidth()*2, 
		(int)COMBO_SIZE.getHeight()*2));


    m_controlPanel = new JPanel();
    m_controlPanel.setLayout(new BorderLayout());
    
    //jimmy
    JPanel dataChooseHolder = new JPanel(new BorderLayout());
    dataChooseHolder.setBorder(BorderFactory.createTitledBorder("Dataset"));
    dataChooseHolder.add(dataFileLabel, BorderLayout.WEST);
    
    m_FileChooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
    m_FileChooser.addChoosableFileFilter(m_arffFileFilter);
    chooseButton.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        try {
		setInstancesFromFileQ();
		int classIndex = m_classAttBox.getSelectedIndex();
		if (m_trainingInstances != null && m_classifier != null && ((m_trainingInstances.attribute(classIndex).isNominal() || (m_trainingInstances.attribute(classIndex).isRanking())) )) {
			m_startBut.setEnabled(true);
//			plotTrainingData();
		}
		
		
	} catch (Exception ex) {
		ex.printStackTrace(System.out);
		System.err.println("exception");
	}
	
      }
    });
    dataChooseHolder.add(chooseButton, BorderLayout.EAST);
    
    JPanel classifierHolder = new JPanel();
    classifierHolder.setBorder(BorderFactory.createTitledBorder("Classifier"));
    classifierHolder.setLayout(new BorderLayout());
    m_classifierEditor.setClassType(weka.classifiers.Classifier.class);
    
    m_classifierEditor.addPropertyChangeListener(new PropertyChangeListener() {
    	public void propertyChange(PropertyChangeEvent evt) {
		m_classifier = (Classifier)m_classifierEditor.getValue();
		try {
			int classIndex = m_classAttBox.getSelectedIndex();
			if (m_trainingInstances != null && m_classifier != null && (m_trainingInstances.attribute(classIndex).isNominal() || m_trainingInstances.attribute(classIndex).isRanking())) {
				m_startBut.setEnabled(true);
			}
		} catch (Exception ex) {};
	}
    });
    classifierHolder.add(m_ClassifierPanel, BorderLayout.CENTER);
    
        

    JPanel cHolder = new JPanel();
    cHolder.setBorder(BorderFactory.createTitledBorder("Class Attribute"));
    cHolder.add(m_classAttBox);

    JPanel vAttHolder = new JPanel();
    vAttHolder.setLayout(new GridLayout(2,1));
    vAttHolder.setBorder(BorderFactory.
			 createTitledBorder("Visualization Attributes"));
    vAttHolder.add(m_xAttBox);
    vAttHolder.add(m_yAttBox);

    JPanel colOne = new JPanel();
    colOne.setLayout(new BorderLayout());
    colOne.add(dataChooseHolder, BorderLayout.NORTH); //jimmy
    colOne.add(cHolder, BorderLayout.CENTER);
    //colOne.add(vAttHolder, BorderLayout.SOUTH);

    JPanel tempPanel = new JPanel();
    tempPanel.setBorder(BorderFactory.
			createTitledBorder("Sampling control"));
    tempPanel.setLayout(new GridLayout(3,1));

    JPanel colTwo = new JPanel();
    colTwo.setLayout(new BorderLayout());
    JPanel gsP = new JPanel(); gsP.setLayout(new BorderLayout());
    gsP.add(new JLabel(" Base for sampling (r)"), BorderLayout.CENTER);
    gsP.add(m_generatorSamplesText, BorderLayout.WEST);
    tempPanel.add(gsP);

    JPanel rsP = new JPanel(); rsP.setLayout(new BorderLayout());
    rsP.add(new JLabel(" Num. locations per pixel"), BorderLayout.CENTER);
    rsP.add(m_regionSamplesText, BorderLayout.WEST);
    tempPanel.add(rsP);

    JPanel ksP = new JPanel(); ksP.setLayout(new BorderLayout());
    ksP.add(new JLabel(" Kernel bandwidth (k)"), BorderLayout.CENTER);
    ksP.add(m_kernelBandwidthText, BorderLayout.WEST);
    tempPanel.add(ksP);
    
    colTwo.add(classifierHolder,BorderLayout.NORTH);//jimmy
    //colTwo.add(tempPanel, BorderLayout.CENTER);
    colTwo.add(vAttHolder, BorderLayout.CENTER);

    JPanel startPanel = new JPanel();
    startPanel.setBorder(BorderFactory.
			 createTitledBorder("Plotting"));
    startPanel.setLayout(new BorderLayout());
    startPanel.add(m_startBut, BorderLayout.CENTER);
    startPanel.add(m_plotTrainingData, BorderLayout.WEST);

    //colTwo.add(startPanel, BorderLayout.SOUTH);

    m_controlPanel.add(colOne, BorderLayout.WEST);
    m_controlPanel.add(colTwo, BorderLayout.CENTER);
    JPanel classHolder = new JPanel();
    classHolder.setLayout(new BorderLayout()); //jimmy
    classHolder.setBorder(BorderFactory.createTitledBorder("Class color"));
    classHolder.add(m_classPanel, BorderLayout.CENTER);
    m_controlPanel.add(classHolder, BorderLayout.SOUTH);
    
    JPanel aboutAndControlP = new JPanel();
    aboutAndControlP.setLayout(new BorderLayout());
    aboutAndControlP.add(m_controlPanel, BorderLayout.SOUTH);
    
    weka.gui.PropertySheetPanel psp = new weka.gui.PropertySheetPanel();
    psp.setTarget(BoundaryVisualizer.this);
    JPanel aboutPanel = psp.getAboutPanel();
    
    aboutAndControlP.add(aboutPanel, BorderLayout.NORTH);

    add(aboutAndControlP, BorderLayout.NORTH);
    
    //classHolder.add(newWindowButton, BorderLayout.EAST);
   
    // set up the add-remove points widgets
    m_addRemovePointsPanel.setBorder(BorderFactory.createTitledBorder("Add / remove data points"));
    m_addRemovePointsPanel.setLayout(new GridBagLayout());
    GridBagConstraints constraints = new GridBagConstraints();
    constraints.weightx = 1.0;
    constraints.weighty = 1.0;
    constraints.gridx = 0;
    constraints.gridy = 0;
    constraints.fill = GridBagConstraints.BOTH;
    m_addRemovePointsPanel.add(m_addPointsButton);
    constraints.gridx = 1;
    m_addRemovePointsPanel.add(new JLabel("Add points"), constraints);
    constraints.gridx = 2;
    m_addRemovePointsPanel.add(m_classValueSelector);
    constraints.gridx = 0;
    constraints.gridy = 1;
    m_addRemovePointsPanel.add(m_removePointsButton, constraints);
    constraints.gridx = 1;
    m_addRemovePointsPanel.add(new JLabel("Remove points"),constraints);
    
    
    	removeAllButton.addActionListener(new ActionListener() {
		public void actionPerformed(ActionEvent e) {
			if (m_trainingInstances != null)
			{
				if (m_startBut.getText().equals("Stop")) //we are plotting
					return;
				m_boundaryPanel.removeAllInstances();
				computeBounds();
				m_xAxisPanel.repaint(0,0,0,m_xAxisPanel.getWidth(), m_xAxisPanel.getHeight());
				m_yAxisPanel.repaint(0,0,0,m_yAxisPanel.getWidth(), m_yAxisPanel.getHeight());
	
				try {m_boundaryPanel.plotTrainingData(); } catch (Exception ex) {}
			}
		}	
	});
    constraints.gridx = 2;
    m_addRemovePointsPanel.add(removeAllButton, constraints);
    
//     m_addRemovePointsPanel.add(addPointsFrame, BorderLayout.NORTH);
//     m_addRemovePointsPanel.add(removePointsFrame, BorderLayout.CENTER);
    //m_addRemovePointsPanel.add(removeAllButton, BorderLayout.SOUTH);
    
    
    m_addRemovePointsButtonGroup.add(m_addPointsButton);
    m_addRemovePointsButtonGroup.add(m_removePointsButton);
    m_addPointsButton.setSelected(true);
        
    //classHolder.add(m_addRemovePointsPanel, BorderLayout.SOUTH);
    

    m_boundaryPanel = new BoundaryPanel(m_plotAreaWidth, m_plotAreaHeight);
    m_numberOfSamplesFromEachRegion = m_boundaryPanel.getNumSamplesPerRegion();
    m_regionSamplesText.setText(""+m_numberOfSamplesFromEachRegion+"  ");
    m_generatorSamplesBase = (int)m_boundaryPanel.getGeneratorSamplesBase();
    m_generatorSamplesText.setText(""+m_generatorSamplesBase+"  ");

    m_dataGenerator = new KDDataGenerator();
    m_kernelBandwidth = m_dataGenerator.getKernelBandwidth();
    m_kernelBandwidthText.setText(""+m_kernelBandwidth+"  ");
    m_boundaryPanel.setDataGenerator(m_dataGenerator);
    
     
    JPanel gfxPanel = new JPanel();
    gfxPanel.setLayout(new BorderLayout());
    gfxPanel.setBorder(BorderFactory.createEtchedBorder());
    //add(gfxPanel, BorderLayout.CENTER);
        
   // gfxPanel.add(m_addRemovePointsPanel, BorderLayout.NORTH);
    gfxPanel.add(m_boundaryPanel, BorderLayout.CENTER);
    m_xAxisPanel = new AxisPanel(false);
    gfxPanel.add(m_xAxisPanel, BorderLayout.SOUTH);
    m_yAxisPanel = new AxisPanel(true);
    gfxPanel.add(m_yAxisPanel, BorderLayout.WEST);
    
    JPanel containerPanel = new JPanel();
    containerPanel.setLayout(new BorderLayout());
    containerPanel.add(gfxPanel, BorderLayout.CENTER);
    add(containerPanel, BorderLayout.WEST);
    
    JPanel rightHandToolsPanel = new JPanel(); //this panel contains the widgets to the right of the BoundaryPanel.
    rightHandToolsPanel.setLayout(new BoxLayout(rightHandToolsPanel, BoxLayout.PAGE_AXIS));

    rightHandToolsPanel.add(m_addRemovePointsPanel);
    
    JButton newWindowButton = new JButton("Open a new window"); //the button for spawning a new window for the program.
    //newWindowButton.setMaximumSize(new Dimension(100, 100));
    //newWindowButton.setPreferredSize(new Dimension(120, m_addRemovePointsPanel.getHeight()));
    newWindowButton.addActionListener(new ActionListener() {
    	public void actionPerformed(ActionEvent e) {
		try {
			Instances newTrainingData = null;
			Classifier newClassifier = null;
			if (m_trainingInstances != null)
				newTrainingData = new Instances(m_trainingInstances);
			if (m_classifier != null)
				newClassifier = AbstractClassifier.makeCopy(m_classifier);
			createNewVisualizerWindow(newClassifier, newTrainingData);
		} catch (Exception ex) {  ex.printStackTrace();}
	}
    });
    JPanel newWindowHolder = new JPanel();
    newWindowHolder.add(newWindowButton);
    rightHandToolsPanel.add(newWindowHolder);
    rightHandToolsPanel.add(tempPanel);
    rightHandToolsPanel.add(startPanel);
    
    containerPanel.add(rightHandToolsPanel, BorderLayout.EAST);
        
    /*add(m_boundaryPanel, BorderLayout.CENTER);

    m_xAxisPanel = new AxisPanel(false);
    add(m_xAxisPanel, BorderLayout.SOUTH);
    m_yAxisPanel = new AxisPanel(true);
    add(m_yAxisPanel, BorderLayout.WEST);*/

    m_startBut.setEnabled(false);
    m_startBut.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e) {
	  if (m_startBut.getText().equals("Start")) {
	    if (m_trainingInstances != null && m_classifier != null) {
		try {
			
			int BPSuccessCode = setUpBoundaryPanel(); //set up the boundary panel, find out if it was successful or not.
			
			if (BPSuccessCode == 1)
				JOptionPane.showMessageDialog(null,"Error: Kernel Bandwidth can't be less than zero!");
			else if (BPSuccessCode == 2) {
				JOptionPane.showMessageDialog(null,"Error: Kernel Bandwidth must be less than the number of training instances!");
			} else {
				m_boundaryPanel.start();
				m_startBut.setText("Stop");
				setControlEnabledStatus(false);
			}
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	    }
	  } else {
	    m_boundaryPanel.stopPlotting();
	    m_startBut.setText("Start");
	    setControlEnabledStatus(true);
	  }
	}
      });

    m_boundaryPanel.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e) {
	  m_startBut.setText("Start");
	  setControlEnabledStatus(true);
	}
      });

    m_classPanel.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e) {

	  try {
	    // save color vector to a file
	    FastVector colors = m_boundaryPanel.getColors();
	    FileOutputStream fos = new FileOutputStream("colors.ser");
	    ObjectOutputStream oos = new ObjectOutputStream(fos);
	    oos.writeObject(colors);
	    oos.flush();
	    oos.close();
	  } catch (Exception ex) {}

	  m_boundaryPanel.replot();
	  
	}
      });
      
    //set up a mouse listener for the boundary panel.
    m_boundaryPanel.addMouseListener(new MouseAdapter() {
    	public void mouseClicked(MouseEvent e) {
// 		System.err.println("boundary panel mouseClick " + e.getX() + " " + e.getY());
		if (m_trainingInstances != null) {
			if (m_startBut.getText().equals("Stop")) //we are plotting
				return;
		
			if (m_addPointsButton.isSelected()) {//we are in add mode
				double classVal = 0;
				boolean validInput = true;
				if (m_trainingInstances.attribute(m_classAttBox.getSelectedIndex()).isNominal() || m_trainingInstances.attribute(m_classAttBox.getSelectedIndex()).isRanking()) //class is nominal
					classVal = (double)m_classValueSelector.getSelectedIndex();
				else {
					String indexStr = "";
					try {					
						indexStr = (String)m_classValueSelector.getSelectedItem();
						classVal = Double.parseDouble(indexStr);
					} catch (Exception ex) {
						if (indexStr == null) indexStr = "";
						JOptionPane.showMessageDialog(null,"Error adding a point: \"" + indexStr + "\""
							+ " is not a valid class value.");
						validInput = false;
					}
				}
				//System.err.println("classVal is " + classVal);
				if (validInput)
					m_boundaryPanel.addTrainingInstanceFromMouseLocation(e.getX(), e.getY(), m_classAttBox.getSelectedIndex(), classVal);
			}
			else { //remove mode
				m_boundaryPanel.removeTrainingInstanceFromMouseLocation(e.getX(), e.getY());
			}
			try{ plotTrainingData(); } catch (Exception ex) {} //jimmy
			m_xAxisPanel.repaint(0,0,0,m_xAxisPanel.getWidth(), m_xAxisPanel.getHeight());
    			m_yAxisPanel.repaint(0,0,0,m_yAxisPanel.getWidth(), m_yAxisPanel.getHeight());
		}
	}
    });
  }
    
  /**
   * Set the enabled status of the controls
   *
   * @param status a <code>boolean</code> value
   */
  private void setControlEnabledStatus(boolean status) {
    m_classAttBox.setEnabled(status);
    m_xAttBox.setEnabled(status);
    m_yAttBox.setEnabled(status);
    m_regionSamplesText.setEnabled(status);
    m_generatorSamplesText.setEnabled(status);
    m_kernelBandwidthText.setEnabled(status);
    m_plotTrainingData.setEnabled(status);
    removeAllButton.setEnabled(status);
    m_classValueSelector.setEnabled(status);
    m_addPointsButton.setEnabled(status);
    m_removePointsButton.setEnabled(status);
    m_FileChooser.setEnabled(status);
    chooseButton.setEnabled(status);
  }

  /**
   * Set a classifier to use
   *
   * @param newClassifier the classifier to use
   * @exception Exception if an error occurs
   */
  public void setClassifier(Classifier newClassifier) throws Exception {

    m_classifier = newClassifier;
    
    try {
	int classIndex = m_classAttBox.getSelectedIndex();
	
	if ((m_classifier != null) && (m_trainingInstances != null) &&
		(m_trainingInstances.attribute(classIndex).isNominal() || m_trainingInstances.attribute(classIndex).isRanking())) {
		m_startBut.setEnabled(true);
	}
	else
		m_startBut.setEnabled(false);
    } catch (Exception e) {}
    
  }
  
  /** Sets up the bounds on our x and y axes to fit the dataset.
      Also repaints the x and y axes.
  */
  private void computeBounds() {
  
    m_boundaryPanel.computeMinMaxAtts(); //delegate to the BoundaryPanel
  
    String xName = (String)m_xAttBox.getSelectedItem();
    if (xName == null) {
      return;
    }
    xName = Utils.removeSubstring(xName, "X: ");
    xName = Utils.removeSubstring(xName, " (Num)");
    String yName = (String)m_yAttBox.getSelectedItem();
    yName = Utils.removeSubstring(yName, "Y: ");
    yName = Utils.removeSubstring(yName, " (Num)");

    m_xIndex = -1;
    m_yIndex = -1;
    for (int i = 0; i < m_trainingInstances.numAttributes(); i++) {
      if (m_trainingInstances.attribute(i).name().equals(xName)) {
	m_xIndex = i;
      } 
      if (m_trainingInstances.attribute(i).name().equals(yName)) {
	m_yIndex = i;
      }
    }
    
    m_minX = m_boundaryPanel.getMinXBound();
    m_minY = m_boundaryPanel.getMinYBound();
    m_maxX = m_boundaryPanel.getMaxXBound();
    m_maxY = m_boundaryPanel.getMaxYBound();
    //System.err.println("setting bounds to " + m_minX + " " + m_minY + " " + m_maxX + " " + m_maxY);
    m_xAxisPanel.repaint(0,0,0,m_xAxisPanel.getWidth(), m_xAxisPanel.getHeight());
    m_yAxisPanel.repaint(0,0,0,m_yAxisPanel.getWidth(), m_yAxisPanel.getHeight());
  }

  /**
   * Get the training instances
   *
   * @return the training instances
   */
  public Instances getInstances() {
    return m_trainingInstances;
  }

  /**
   * Set the training instances
   *
   * @param inst the instances to use
   */
  public void setInstances(Instances inst) throws Exception {
    if (inst == null) {
    	m_trainingInstances = inst;
    	m_classPanel.setInstances(m_trainingInstances);
	return;
    }
    
    // count the number of numeric attributes
    int numCount = 0;
    for (int i = 0; i < inst.numAttributes(); i++) {
      if (inst.attribute(i).isNumeric()) {
        numCount++;
      }
    }
    
    if (numCount < 2) {
      JOptionPane.showMessageDialog(null,"We need at least two numeric " +
      		"attributes in order to visualize!");
      return;
    }
        
    m_trainingInstances = inst;
    m_classPanel.setInstances(m_trainingInstances);
    // setup combo boxes
    String [] classAttNames = new String [m_trainingInstances.numAttributes()];
    final Vector xAttNames = new Vector();
    Vector yAttNames = new Vector();

    for (int i = 0; i < m_trainingInstances.numAttributes(); i++) {
      classAttNames[i] = m_trainingInstances.attribute(i).name();
      String type = "";
      switch (m_trainingInstances.attribute(i).type()) {
    case PreferenceAttribute.RANKING:
      type = " (Rnk)";
      break;
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
	default:
	  type = " (???)";
      }
      classAttNames[i] += type;
      if (m_trainingInstances.attribute(i).isNumeric()) {
	xAttNames.addElement("X: "+classAttNames[i]);
	yAttNames.addElement("Y: "+classAttNames[i]);
      }
    }

    m_classAttBox.setModel(new DefaultComboBoxModel(classAttNames));
    m_xAttBox.setModel(new DefaultComboBoxModel(xAttNames));
    m_yAttBox.setModel(new DefaultComboBoxModel(yAttNames));
    if (xAttNames.size() > 1) {
      m_yAttBox.setSelectedIndex(1);
    }
    
    m_classAttBox.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e) {
	  configureForClassAttribute();
	}
      });

    m_xAttBox.addItemListener(new ItemListener() {
	public void itemStateChanged(ItemEvent e) {
	  if (e.getStateChange() == ItemEvent.SELECTED) {
/*	    if (xAttNames.size() > 1) {
	      if (m_xAttBox.getSelectedIndex() == 
		  m_yAttBox.getSelectedIndex()) {
		m_xAttBox.setSelectedIndex((m_xAttBox.getSelectedIndex() + 1) %
					   xAttNames.size());
	      }
	    } */
	    computeBounds();
	    repaint();
	    try{ plotTrainingData(); } catch (Exception ex) {ex.printStackTrace();} //jimmy	    
	  }
	}
      });

    m_yAttBox.addItemListener(new ItemListener() {
	public void itemStateChanged(ItemEvent e) {
	  if (e.getStateChange() == ItemEvent.SELECTED) {
/*	    if (xAttNames.size() > 1) {
	      if (m_yAttBox.getSelectedIndex() == 
		  m_xAttBox.getSelectedIndex()) {
		m_yAttBox.setSelectedIndex((m_yAttBox.getSelectedIndex() + 1) %
					   xAttNames.size());
	      }
	    } */
	    computeBounds();
	    repaint();
	    try{ plotTrainingData(); } catch (Exception ex) {ex.printStackTrace();}
	  }
	}
      });
    
    if (classAttNames.length > 0)
      m_classAttBox.setSelectedIndex(classAttNames.length - 1); //select last attribute as class by default.  -jimmy
      
    //set up the add points selector combo box
    setUpClassValueSelectorCB();
    
    configureForClassAttribute();
    
    m_classPanel.setCindex(m_classAttBox.getSelectedIndex());
    plotTrainingData();      
    computeBounds();
    revalidate();
    repaint();
    
    if (getTopLevelAncestor() instanceof java.awt.Window) {
      ((java.awt.Window)getTopLevelAncestor()).pack();
    }
  }
  
  /** Set up the combo box that chooses which class values to use when adding data points.
  */
  private void setUpClassValueSelectorCB() {
    m_classValueSelector.removeAllItems();
    int classAttribute = m_classAttBox.getSelectedIndex();
    //System.err.println(m_trainingInstances.numClasses() + " classes");
    m_trainingInstances.setClassIndex(classAttribute);
    if (m_trainingInstances.attribute(classAttribute).isNominal() || m_trainingInstances.attribute(classAttribute).isRanking()) {
    	m_classValueSelector.setEditable(false);
    	for (int i = 0; i < /*m_trainingInstances.numDistinctValues(classAttribute)*/m_trainingInstances.numClasses(); i++)
    		m_classValueSelector.insertItemAt(m_trainingInstances.attribute(classAttribute).value(i) , i);
	m_classValueSelector.setSelectedIndex(0);
    }
    else {
    	m_classValueSelector.setEditable(true);
    }
  }
  
  /**
   * Set up the class values combo boxes
   */
  private void configureForClassAttribute() {
    int classIndex = m_classAttBox.getSelectedIndex();
    if (classIndex >= 0) {
      // see if this is a nominal attribute
      if ((!m_trainingInstances.attribute(classIndex).isNominal() && !m_trainingInstances.attribute(classIndex).isRanking())  || m_classifier == null) {
	m_startBut.setEnabled(false);
      } else {
	m_startBut.setEnabled(true);
      }
      // set up class colours
	FastVector colors = new FastVector();
	if ((!m_trainingInstances.attribute(m_classAttBox.getSelectedIndex()).isNominal() && !m_trainingInstances.attribute(m_classAttBox.getSelectedIndex()).isRanking())) //this if by jimmy
	{
		for (int i = 0; i < BoundaryPanel.DEFAULT_COLORS.length; i++)
			colors.addElement(BoundaryPanel.DEFAULT_COLORS[i]);
	}
	else {
		for (int i = 0; i < 
		m_trainingInstances.attribute(classIndex).numValues(); i++) {
			colors.addElement(BoundaryPanel.
				DEFAULT_COLORS[i % BoundaryPanel.DEFAULT_COLORS.length]);
// 			m_classPanel.setColours(colors);	  
// 			m_boundaryPanel.setColors(colors);
		}
	}
	m_classPanel.setColours(colors); //jimmy
	m_boundaryPanel.setColors(colors);
   }
  }
  
    
  /**
   * Queries the user for a file to load instances from, then loads the
   * instances in a background process. This is done in the IO
   * thread, and an error message is popped up if the IO thread is busy.
   */
  public void setInstancesFromFileQ() {
    
//     if (m_IOThread == null) {
      int returnVal = m_FileChooser.showOpenDialog(this);
      if (returnVal == JFileChooser.APPROVE_OPTION) {
	File selected = m_FileChooser.getSelectedFile();
	
	try
	{
	java.io.Reader r = new java.io.BufferedReader(
				new java.io.FileReader(selected));
	Instances i = new Instances(r);
	setInstances(i);
	
	//dataFileLabel.setText(selected.getName());
	String relationName = i.relationName();
	String truncatedN = relationName;
	if (relationName.length() > 25) {
	  truncatedN = relationName.substring(0, 25) + "...";
	}
	dataFileLabel.setText(truncatedN);
	dataFileLabel.setToolTipText(relationName);
	} catch (Exception e)
	{
		JOptionPane.showMessageDialog(this,"Can't load at this time,\n"
				    + "currently busy with other IO",
				    "Load Instances",
				    JOptionPane.WARNING_MESSAGE);
		    e.printStackTrace();
	
	}
      }
  }
  
  /** Sets up the BoundaryPanel object so that it is ready for plotting.
   * @return an error code:<br/>
   *		0 - SUCCESS<br/>
   *		1 - ERROR - Kernel bandwidth < 0<br/>
   *		2 - ERROR - Kernel bandwidth >= number of training instances.
   */
  public int setUpBoundaryPanel() throws Exception {
  	int returner = 0; //OK code.
  	int tempSamples = m_numberOfSamplesFromEachRegion;
		try {
		  tempSamples = 
		    Integer.parseInt(m_regionSamplesText.getText().trim());
		} catch (Exception ex) {
		  m_regionSamplesText.setText(""+tempSamples);
		}
		m_numberOfSamplesFromEachRegion = tempSamples;
		m_boundaryPanel.
		  setNumSamplesPerRegion(tempSamples);

		tempSamples = m_generatorSamplesBase;
		try {
		  tempSamples = 
		    Integer.parseInt(m_generatorSamplesText.getText().trim());
		} catch (Exception ex) {
		  m_generatorSamplesText.setText(""+tempSamples);
		}
		m_generatorSamplesBase = tempSamples;
		m_boundaryPanel.setGeneratorSamplesBase((double)tempSamples);

		tempSamples = m_kernelBandwidth;
		try {
		  tempSamples = 
		    Integer.parseInt(m_kernelBandwidthText.getText().trim());
		} catch (Exception ex) {
		  m_kernelBandwidthText.setText(""+tempSamples);
		}
		m_kernelBandwidth = tempSamples;
		m_dataGenerator.setKernelBandwidth(tempSamples);
		
		if (m_kernelBandwidth < 0)	returner = 1;
		if (m_kernelBandwidth >= m_trainingInstances.numInstances())	returner = 2;

		m_trainingInstances.
		  setClassIndex(m_classAttBox.getSelectedIndex());
		m_boundaryPanel.setClassifier(m_classifier);
		m_boundaryPanel.setTrainingData(m_trainingInstances);
		m_boundaryPanel.setXAttribute(m_xIndex);
		m_boundaryPanel.setYAttribute(m_yIndex);
		m_boundaryPanel.
		  setPlotTrainingData(m_plotTrainingData.isSelected());

	return returner;
  }
  
  /** Plots the training data on-screen.  Also does all of the setup required 
   *  for this to work.
  */
  public void plotTrainingData() throws Exception {
  	m_boundaryPanel.initialize();
 	setUpBoundaryPanel();
	computeBounds();
	m_boundaryPanel.plotTrainingData();
  }
  
  /** Stops the plotting thread.
  */
  public void stopPlotting() {
  	m_boundaryPanel.stopPlotting();
  }
  
  /**
   * Sets whether System.exit gets called when no more windows are open.
   * 
   * @param value	if TRUE then a System.exit call is ossued after the 
   * 			last window gets closed.
   */
  public static void setExitIfNoWindowsOpen(boolean value) {
    m_ExitIfNoWindowsOpen = value;
  }
  
  /**
   * Gets whether System.exit gets called after the last window gets closed
   * 
   * @return		TRUE if System.exit gets called after last window
   * 			got closed.
   */
  public static boolean getExitIfNoWindowsOpen() {
    return m_ExitIfNoWindowsOpen;
  }
  
  /** Creates a new GUI window with all of the BoundaryVisualizer trappings,
   *  @param classifier The classifier to use in the new window.  May be null.
   *  @param instances  The dataset to visualize on in the new window.  May be null.
   */
  public static void createNewVisualizerWindow(Classifier classifier, Instances instances) throws Exception {
      m_WindowCount++;
  
      final javax.swing.JFrame jf = 
	new javax.swing.JFrame("Weka classification boundary visualizer");
      jf.getContentPane().setLayout(new BorderLayout());
      final BoundaryVisualizer bv = new BoundaryVisualizer();
      jf.getContentPane().add(bv, BorderLayout.CENTER);
      jf.setSize(bv.getMinimumSize());
      jf.addWindowListener(new java.awt.event.WindowAdapter() {
	  public void windowClosing(java.awt.event.WindowEvent e) {
	    m_WindowCount--;
	    bv.stopPlotting();
	    jf.dispose();
	    if ((m_WindowCount == 0) && m_ExitIfNoWindowsOpen) {
		System.exit(0);
	    }
	  }
	});

      jf.pack();
      jf.setVisible(true);
      jf.setResizable(false);
      
      if (classifier == null)
      	bv.setClassifier(null);
      else {
	bv.setClassifier(classifier);
	bv.m_classifierEditor.setValue(classifier);
      }
      
      if (instances == null)
      	bv.setInstances(null);
      else
      {
	bv.setInstances(instances);
	
	try{
		bv.dataFileLabel.setText(instances.relationName());
		bv.plotTrainingData();
		bv.m_classPanel.setCindex(bv.m_classAttBox.getSelectedIndex());
		bv.repaint(0,0,0,bv.getWidth(), bv.getHeight());
	} catch (Exception ex) {}
      }
  
  }

  /**
   * Main method for testing this class
   *
   * @param args a <code>String[]</code> value
   */
  public static void main(String [] args) {
    weka.core.logging.Logger.log(weka.core.logging.Logger.Level.INFO, "Logging started");
    try {
    	if (args.length < 2) {
		createNewVisualizerWindow(null, null);
	}
	else {
		String [] argsR = null;
		if (args.length > 2) {
			argsR = new String [args.length-2];
			for (int j = 2; j < args.length; j++) {
			argsR[j-2] = args[j];
			}
		}
		Classifier c = AbstractClassifier.forName(args[1], argsR);
		
		System.err.println("Loading instances from : "+args[0]);
		java.io.Reader r = new java.io.BufferedReader(
					new java.io.FileReader(args[0]));
		Instances i = new Instances(r);
	
		createNewVisualizerWindow(c, i);
	}
      
    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }
}


