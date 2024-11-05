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
 *    SerializedModelSaverCustomizer.java
 *    Copyright (C) 2008 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

import weka.gui.GenericObjectEditor;
import weka.gui.PropertySheetPanel;
import weka.core.Environment;
import weka.core.EnvironmentHandler;
import weka.core.Tag;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.beans.Customizer;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeSupport;
import java.io.File;

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.BorderFactory;

import javax.swing.JTextField;
import javax.swing.SwingConstants;
import javax.swing.filechooser.FileFilter;

/**
 * GUI Customizer for the SerializedModelSaver bean
 *
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}org
 * @version $Revision: 5563 $
 */
public class SerializedModelSaverCustomizer
  extends JPanel
  implements Customizer, CustomizerCloseRequester, EnvironmentHandler {

  /** for serialization */
  private static final long serialVersionUID = -4874208115942078471L;

  static {
     GenericObjectEditor.registerEditors();
  }

  private PropertyChangeSupport m_pcSupport = 
    new PropertyChangeSupport(this);

  private weka.gui.beans.SerializedModelSaver m_smSaver;

  private PropertySheetPanel m_SaverEditor = 
    new PropertySheetPanel();

  private JFileChooser m_fileChooser 
    = new JFileChooser(new File(System.getProperty("user.dir")));
  

  private JFrame m_parentFrame;
  
  private JFrame m_fileChooserFrame;
  
  //private JTextField m_prefixText;
  private EnvironmentField m_prefixText;

  private JComboBox m_fileFormatBox;

  private JCheckBox m_relativeFilePath;
  
  private JCheckBox m_includeRelationName;
  
  private Environment m_env = Environment.getSystemWide();
  
  private EnvironmentField m_directoryText;
  

  /** Constructor */  
  public SerializedModelSaverCustomizer() {

    try {
      m_SaverEditor.addPropertyChangeListener(
	  new PropertyChangeListener() {
	      public void propertyChange(PropertyChangeEvent e) {
		repaint();
		if (m_smSaver != null) {
		  System.err.println("Property change!!");
                  //		  m_smSaver.setSaver(m_smSaver.getSaver());
		}
	      }
	    });
      repaint();
    } catch (Exception ex) {
      ex.printStackTrace();
    }
    setLayout(new BorderLayout());

    m_fileChooser.setDialogType(JFileChooser.SAVE_DIALOG);
    m_fileChooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
    m_fileChooser.setApproveButtonText("Select directory and prefix");

    m_fileChooser.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e) {
	  if (e.getActionCommand().equals(JFileChooser.APPROVE_SELECTION)) {
	    try {
              m_smSaver.setPrefix(m_prefixText.getText());
//              m_smSaver.setDirectory(m_fileChooser.getSelectedFile());
              
              File selectedFile = m_fileChooser.getSelectedFile();
              m_directoryText.setText(selectedFile.toString());
              
	    } catch (Exception ex) {
	      ex.printStackTrace();
	    }
	  }
	  // closing
	  if (m_parentFrame != null) {
	    m_fileChooserFrame.dispose();
	  }
	}
      });   
  }

  public void setParentFrame(JFrame parent) {
    m_parentFrame = parent;
  }
  
  private void setUpOther() {
    removeAll();
    add(m_SaverEditor, BorderLayout.CENTER);
    validate();
    repaint();
  }
  
  /** Sets up dialog for saving models to a file */  
  public void setUpFile() {
    removeAll();
    m_fileChooser.setFileFilter(new FileFilter() { 
      public boolean accept(File f) { 
        return f.isDirectory();
      }
      public String getDescription() {
        return "Directory";
      }
    });

    m_fileChooser.setAcceptAllFileFilterUsed(false);

    try{
      if (!m_smSaver.getDirectory().getPath().equals("")) {
       // File tmp = m_smSaver.getDirectory();
        String dirStr = m_smSaver.getDirectory().toString();
        if (Environment.containsEnvVariables(dirStr)) {
          try {
            dirStr = m_env.substitute(dirStr);
          } catch (Exception ex) {
            // ignore
          }
        }
        File tmp = new File(dirStr);;
        tmp = new File (tmp.getAbsolutePath());
        m_fileChooser.setCurrentDirectory(tmp);
      }
    } catch(Exception ex) {
      System.out.println(ex);
    }

    JPanel innerPanel = new JPanel();
    innerPanel.setLayout(new BorderLayout());
    
    JPanel alignedP = new JPanel();
    GridBagLayout gbLayout = new GridBagLayout();
    alignedP.setLayout(gbLayout);
    
    JLabel prefixLab = new JLabel("Prefix for file name", SwingConstants.RIGHT);
    prefixLab.setBorder(BorderFactory.createEmptyBorder(0, 5, 0, 0));
    GridBagConstraints gbConstraints = new GridBagConstraints();
    gbConstraints.anchor = GridBagConstraints.EAST;
    gbConstraints.fill = GridBagConstraints.HORIZONTAL;
    gbConstraints.gridy = 0; gbConstraints.gridx = 0;
    gbLayout.setConstraints(prefixLab, gbConstraints);
    alignedP.add(prefixLab);
    
//    m_prefixText = new JTextField(m_smSaver.getPrefix(), 25);
    m_prefixText = new EnvironmentField();
    m_prefixText.setEnvironment(m_env);
/*    int width = m_prefixText.getPreferredSize().width;
    int height = m_prefixText.getPreferredSize().height;
    m_prefixText.setMinimumSize(new Dimension(width * 2, height));
    m_prefixText.setPreferredSize(new Dimension(width * 2, height)); */
    m_prefixText.setText(m_smSaver.getPrefix());
    gbConstraints = new GridBagConstraints();
    gbConstraints.anchor = GridBagConstraints.EAST;
    gbConstraints.fill = GridBagConstraints.HORIZONTAL;
    gbConstraints.gridy = 0; gbConstraints.gridx = 1;
    gbLayout.setConstraints(m_prefixText, gbConstraints);
    alignedP.add(m_prefixText);
    
    JLabel ffLab = new JLabel("File format", SwingConstants.RIGHT);
    ffLab.setBorder(BorderFactory.createEmptyBorder(0, 5, 0, 0));
    gbConstraints = new GridBagConstraints();
    gbConstraints.anchor = GridBagConstraints.EAST;
    gbConstraints.fill = GridBagConstraints.HORIZONTAL;
    gbConstraints.gridy = 1; gbConstraints.gridx = 0;
    gbLayout.setConstraints(ffLab, gbConstraints);
    alignedP.add(ffLab);
    
    setUpFileFormatComboBox();
    m_fileFormatBox.setBorder(BorderFactory.createEmptyBorder(0, 5, 0, 5));
    gbConstraints = new GridBagConstraints();
    gbConstraints.anchor = GridBagConstraints.EAST;
    gbConstraints.fill = GridBagConstraints.HORIZONTAL;
    gbConstraints.gridy = 1; gbConstraints.gridx = 1;
    gbLayout.setConstraints(m_fileFormatBox, gbConstraints);
    alignedP.add(m_fileFormatBox);

    JPanel about = m_SaverEditor.getAboutPanel();
    if (about != null) {
      innerPanel.add(about, BorderLayout.NORTH);
    }
    add(innerPanel, BorderLayout.NORTH);
    
    JLabel directoryLab = new JLabel("Directory", SwingConstants.RIGHT);
    directoryLab.setBorder(BorderFactory.createEmptyBorder(0, 5, 0, 0));
    gbConstraints = new GridBagConstraints();
    gbConstraints.anchor = GridBagConstraints.EAST;
    gbConstraints.fill = GridBagConstraints.HORIZONTAL;
    gbConstraints.gridy = 2; gbConstraints.gridx = 0;
    gbLayout.setConstraints(directoryLab, gbConstraints);
    alignedP.add(directoryLab);
    
    m_directoryText = new EnvironmentField();
    m_directoryText.setEnvironment(m_env);  
/*    width = m_directoryText.getPreferredSize().width;
    height = m_directoryText.getPreferredSize().height;
    m_directoryText.setMinimumSize(new Dimension(width * 2, height));
    m_directoryText.setPreferredSize(new Dimension(width * 2, height)); */
    
    m_directoryText.setText(m_smSaver.getDirectory().toString());
    
    JButton browseBut = new JButton("Browse...");
    browseBut.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        try {
          final JFrame jf = new JFrame("Choose directory");
          jf.getContentPane().setLayout(new BorderLayout());
          jf.getContentPane().add(m_fileChooser, BorderLayout.CENTER);
          jf.pack();
          jf.setVisible(true);
          m_fileChooserFrame = jf;
        } catch (Exception ex) {
          ex.printStackTrace();
        }
      }
    });
    
    JPanel efHolder = new JPanel();
    efHolder.setLayout(new BorderLayout());
    JPanel bP = new JPanel(); bP.setLayout(new BorderLayout());
    bP.setBorder(BorderFactory.createEmptyBorder(5,0,5,5));
    bP.add(browseBut, BorderLayout.CENTER);
    efHolder.add(bP, BorderLayout.EAST);
    efHolder.add(m_directoryText, BorderLayout.CENTER);
    gbConstraints = new GridBagConstraints();
    gbConstraints.anchor = GridBagConstraints.EAST;
    gbConstraints.fill = GridBagConstraints.HORIZONTAL;
    gbConstraints.gridy = 2; gbConstraints.gridx = 1;
    gbConstraints.weightx = 5; // make sure that extra horizontal space gets allocated to this column
    gbLayout.setConstraints(efHolder, gbConstraints);
    alignedP.add(efHolder);

    
    JLabel relativeLab = new JLabel("Use relative file paths", SwingConstants.RIGHT);
    relativeLab.setBorder(BorderFactory.createEmptyBorder(0, 5, 0, 0));
    gbConstraints = new GridBagConstraints();
    gbConstraints.anchor = GridBagConstraints.EAST;
    gbConstraints.fill = GridBagConstraints.HORIZONTAL;
    gbConstraints.gridy = 3; gbConstraints.gridx = 0;
    gbLayout.setConstraints(relativeLab, gbConstraints);
    alignedP.add(relativeLab);
    
    m_relativeFilePath = new JCheckBox();
    m_relativeFilePath.
      setSelected(m_smSaver.getUseRelativePath());

    m_relativeFilePath.addActionListener(new ActionListener() {
        public void actionPerformed(ActionEvent e) {
          m_smSaver.setUseRelativePath(m_relativeFilePath.isSelected());
        }
      });
    gbConstraints = new GridBagConstraints();
    gbConstraints.anchor = GridBagConstraints.EAST;
    gbConstraints.fill = GridBagConstraints.HORIZONTAL;
    gbConstraints.gridy = 3; gbConstraints.gridx = 1;
    gbLayout.setConstraints(m_relativeFilePath, gbConstraints);
    alignedP.add(m_relativeFilePath);
    
    
    JLabel relationLab = new JLabel("Include relation name in file name", SwingConstants.RIGHT);
    relationLab.setBorder(BorderFactory.createEmptyBorder(0, 5, 0, 0));
    gbConstraints = new GridBagConstraints();
    gbConstraints.anchor = GridBagConstraints.EAST;
    gbConstraints.fill = GridBagConstraints.HORIZONTAL;
    gbConstraints.gridy = 4; gbConstraints.gridx = 0;
    gbLayout.setConstraints(relationLab, gbConstraints);
    alignedP.add(relationLab);
    
    m_includeRelationName = new JCheckBox();
    m_includeRelationName.setToolTipText("Include the relation name of the training data used "
        + "to create the model in the file name.");
    m_includeRelationName.setSelected(m_smSaver.getIncludeRelationName());
        
    gbConstraints = new GridBagConstraints();
    gbConstraints.anchor = GridBagConstraints.EAST;
    gbConstraints.fill = GridBagConstraints.HORIZONTAL;
    gbConstraints.gridy = 4; gbConstraints.gridx = 1;
    gbLayout.setConstraints(m_includeRelationName, gbConstraints);
    alignedP.add(m_includeRelationName);
    
    JButton OKBut = new JButton("OK");
    OKBut.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        try {          
          m_smSaver.setPrefix(m_prefixText.getText());
          m_smSaver.setDirectory(new File(m_directoryText.getText()));
          m_smSaver.
            setIncludeRelationName(m_includeRelationName.isSelected());
        } catch (Exception ex) {
          ex.printStackTrace();
        }
        
        m_parentFrame.dispose();
      }
    });

    JButton CancelBut = new JButton("Cancel");
    CancelBut.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        m_parentFrame.dispose();
      }
    });
    
    JPanel butHolder = new JPanel();
    butHolder.setLayout(new FlowLayout());
    butHolder.add(OKBut);
    butHolder.add(CancelBut);
  
    JPanel holderPanel = new JPanel();
    holderPanel.setLayout(new BorderLayout());
    holderPanel.add(alignedP, BorderLayout.NORTH);
    holderPanel.add(butHolder, BorderLayout.SOUTH);
    add(holderPanel, BorderLayout.SOUTH);
  }

  /**
   * Set the model saver to be customized
   *
   * @param object a weka.gui.beans.SerializedModelSaver
   */
  public void setObject(Object object) {
    m_smSaver = (weka.gui.beans.SerializedModelSaver)object;
    m_SaverEditor.setTarget(m_smSaver);

    setUpFile();    
  }

  private void setUpFileFormatComboBox() {
    m_fileFormatBox = new JComboBox();
    for (int i = 0; i < SerializedModelSaver.s_fileFormatsAvailable.size(); i++) {
      Tag temp = SerializedModelSaver.s_fileFormatsAvailable.get(i);
      m_fileFormatBox.addItem(temp);
    }

    Tag result = m_smSaver.validateFileFormat(m_smSaver.getFileFormat());
    if (result == null) {
      m_fileFormatBox.setSelectedIndex(0);
    } else {
      m_fileFormatBox.setSelectedItem(result);
    }

    m_fileFormatBox.addActionListener(new ActionListener() {
        public void actionPerformed(ActionEvent e) {
          Tag selected = (Tag)m_fileFormatBox.getSelectedItem();
          if (selected != null) {
            m_smSaver.setFileFormat(selected);
          }
        }
      });
  }

  /**
   * Add a property change listener
   *
   * @param pcl a <code>PropertyChangeListener</code> value
   */
  public void addPropertyChangeListener(PropertyChangeListener pcl) {
    m_pcSupport.addPropertyChangeListener(pcl);
  }

  /**
   * Remove a property change listener
   *
   * @param pcl a <code>PropertyChangeListener</code> value
   */
  public void removePropertyChangeListener(PropertyChangeListener pcl) {
    m_pcSupport.removePropertyChangeListener(pcl);
  }

  /* (non-Javadoc)
   * @see weka.core.EnvironmentHandler#setEnvironment(weka.core.Environment)
   */
  public void setEnvironment(Environment env) {
    m_env = env;
  }
}
