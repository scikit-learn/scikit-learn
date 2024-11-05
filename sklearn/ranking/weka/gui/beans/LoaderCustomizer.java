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
 *    LoaderCustomizer.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

import weka.core.Environment;
import weka.core.EnvironmentHandler;
import weka.core.converters.DatabaseConverter;
import weka.core.converters.DatabaseLoader;
import weka.core.converters.FileSourcedConverter;
import weka.gui.ExtensionFileFilter;
import weka.gui.GenericObjectEditor;
import weka.gui.PropertySheetPanel;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.Customizer;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeSupport;
import java.io.File;
import java.util.ArrayList;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JPasswordField;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
import javax.swing.JCheckBox;

/**
 * GUI Customizer for the loader bean
 *
 * @author <a href="mailto:mhall@cs.waikato.ac.nz">Mark Hall</a>
 * @version $Revision: 5406 $
 */
public class LoaderCustomizer
  extends JPanel
  implements Customizer, CustomizerCloseRequester, EnvironmentHandler {

  /** for serialization */
  private static final long serialVersionUID = 6990446313118930298L;

  static {
     GenericObjectEditor.registerEditors();
  }

  private PropertyChangeSupport m_pcSupport = 
    new PropertyChangeSupport(this);

  private weka.gui.beans.Loader m_dsLoader;

  private PropertySheetPanel m_LoaderEditor = 
    new PropertySheetPanel();

  private JFileChooser m_fileChooser 
    = new JFileChooser(new File(System.getProperty("user.dir")));
  /*  private JDialog m_chooserDialog = 
    new JDialog((JFrame)getTopLevelAncestor(),
    true); */

  private JFrame m_parentFrame;
  private JFrame m_fileChooserFrame;
  
  private EnvironmentField m_dbaseURLText;
  
  private EnvironmentField m_userNameText;
  
  private EnvironmentField m_queryText;
   
  private EnvironmentField m_keyText;
  
  private JPasswordField m_passwordText;

  private JCheckBox m_relativeFilePath;
  
  private EnvironmentField m_fileText;
  
  private Environment m_env = Environment.getSystemWide();

  public LoaderCustomizer() {
    /*    m_fileEditor.addPropertyChangeListener(new PropertyChangeListener() {
	public void propertyChange(PropertyChangeEvent e) {
	  if (m_dsLoader != null) {
	    m_dsLoader.setDataSetFile((File)m_fileEditor.getValue());
	  }
	}
	}); */

    try {
      /*      m_LoaderEditor.setClassType(weka.core.converters.Loader.class);
	      m_LoaderEditor.setValue(new weka.core.converters.ArffLoader()); */
      m_LoaderEditor.addPropertyChangeListener(
	  new PropertyChangeListener() {
	      public void propertyChange(PropertyChangeEvent e) {
		repaint();
		if (m_dsLoader != null) {
		  System.err.println("Property change!!");
		  m_dsLoader.setLoader(m_dsLoader.getLoader());
		}
	      }
	    });
      repaint();
    } catch (Exception ex) {
      ex.printStackTrace();
    }

    setLayout(new BorderLayout());
    //    add(m_fileEditor.getCustomEditor(), BorderLayout.CENTER);
    //    add(m_LoaderEditor, BorderLayout.CENTER);
    m_fileChooser.setDialogType(JFileChooser.OPEN_DIALOG);
    m_fileChooser.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e) {
	  if (e.getActionCommand().equals(JFileChooser.APPROVE_SELECTION)) {
	    try {
              File selectedFile = m_fileChooser.getSelectedFile();
/*              EnvironmentField ef = m_environmentFields.get(0);
              ef.setText(selectedFile.toString()); */
              m_fileText.setText(selectedFile.toString());
              
/*	      ((FileSourcedConverter)m_dsLoader.getLoader()).
		setFile(selectedFile);
	      // tell the loader that a new file has been selected so
	      // that it can attempt to load the header
	      //m_dsLoader.setLoader(m_dsLoader.getLoader());
	      m_dsLoader.newFileSelected(); */
	    } catch (Exception ex) {
	      ex.printStackTrace();
	    }
	  }
	  // closing
	  if (m_fileChooserFrame != null) {
	    // m_parentFrame.dispose();
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
    add(m_LoaderEditor, BorderLayout.CENTER);
    validate();
    repaint();
  }
  
  
  /** Sets up a customizer window for a Database Connection*/
  private void setUpDatabase() {

    removeAll();

    JPanel db = new JPanel();
    GridBagLayout gbLayout = new GridBagLayout();
    //db.setLayout(new GridLayout(6, 1));
    db.setLayout(gbLayout);

    JLabel urlLab = new JLabel("Database URL", SwingConstants.RIGHT);
    urlLab.setBorder(BorderFactory.createEmptyBorder(0, 5, 0, 0));
    GridBagConstraints gbConstraints = new GridBagConstraints();
    gbConstraints.anchor = GridBagConstraints.EAST;
    gbConstraints.fill = GridBagConstraints.HORIZONTAL;
    gbConstraints.gridy = 0; gbConstraints.gridx = 0;
    gbLayout.setConstraints(urlLab, gbConstraints);
    db.add(urlLab);

    m_dbaseURLText = new EnvironmentField();
    m_dbaseURLText.setEnvironment(m_env);
/*    int width = m_dbaseURLText.getPreferredSize().width;
    int height = m_dbaseURLText.getPreferredSize().height;
    m_dbaseURLText.setMinimumSize(new Dimension(width * 2, height));
    m_dbaseURLText.setPreferredSize(new Dimension(width * 2, height)); */
    m_dbaseURLText.setText(((DatabaseConverter)m_dsLoader.getLoader()).getUrl());
    gbConstraints = new GridBagConstraints();
    gbConstraints.anchor = GridBagConstraints.EAST;
    gbConstraints.fill = GridBagConstraints.HORIZONTAL;
    gbConstraints.gridy = 0; gbConstraints.gridx = 1;
    gbConstraints.weightx = 5;
    gbLayout.setConstraints(m_dbaseURLText, gbConstraints);
    db.add(m_dbaseURLText);
    
    JLabel userLab = new JLabel("Username", SwingConstants.RIGHT);
    userLab.setBorder(BorderFactory.createEmptyBorder(0, 5, 0, 0));
    gbConstraints = new GridBagConstraints();
    gbConstraints.anchor = GridBagConstraints.EAST;
    gbConstraints.fill = GridBagConstraints.HORIZONTAL;
    gbConstraints.gridy = 1; gbConstraints.gridx = 0;
    gbLayout.setConstraints(userLab, gbConstraints);
    db.add(userLab);
    
    m_userNameText = new EnvironmentField();
    m_userNameText.setEnvironment(m_env);
/*    m_userNameText.setMinimumSize(new Dimension(width * 2, height));
    m_userNameText.setPreferredSize(new Dimension(width * 2, height)); */
    m_userNameText.setText(((DatabaseConverter)m_dsLoader.getLoader()).getUser()); 
    gbConstraints = new GridBagConstraints();
    gbConstraints.anchor = GridBagConstraints.EAST;
    gbConstraints.fill = GridBagConstraints.HORIZONTAL;
    gbConstraints.gridy = 1; gbConstraints.gridx = 1;
    gbLayout.setConstraints(m_userNameText, gbConstraints);
    db.add(m_userNameText);

    JLabel passwordLab = new JLabel("Password ", SwingConstants.RIGHT);
    passwordLab.setBorder(BorderFactory.createEmptyBorder(0, 5, 0, 0));
    gbConstraints = new GridBagConstraints();
    gbConstraints.anchor = GridBagConstraints.EAST;
    gbConstraints.fill = GridBagConstraints.HORIZONTAL;
    gbConstraints.gridy = 2; gbConstraints.gridx = 0;
    gbLayout.setConstraints(passwordLab, gbConstraints);
    db.add(passwordLab);
    
    m_passwordText = new JPasswordField();     
    JPanel passwordHolder = new JPanel();
    passwordHolder.setLayout(new BorderLayout());
    passwordHolder.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
//    passwordHolder.add(passwordLab, BorderLayout.WEST);
    passwordHolder.add(m_passwordText, BorderLayout.CENTER);
/*    passwordHolder.setMinimumSize(new Dimension(width * 2, height));
    passwordHolder.setPreferredSize(new Dimension(width * 2, height)); */
    gbConstraints = new GridBagConstraints();
    gbConstraints.anchor = GridBagConstraints.EAST;
    gbConstraints.fill = GridBagConstraints.HORIZONTAL;
    gbConstraints.gridy = 2; gbConstraints.gridx = 1;
    gbLayout.setConstraints(passwordHolder, gbConstraints);
    db.add(passwordHolder);

    JLabel queryLab = new JLabel("Query", SwingConstants.RIGHT);
    queryLab.setBorder(BorderFactory.createEmptyBorder(0, 5, 0, 0));
    gbConstraints = new GridBagConstraints();
    gbConstraints.anchor = GridBagConstraints.EAST;
    gbConstraints.fill = GridBagConstraints.HORIZONTAL;
    gbConstraints.gridy = 3; gbConstraints.gridx = 0;
    gbLayout.setConstraints(queryLab, gbConstraints);
    db.add(queryLab);
    
    m_queryText = new EnvironmentField();
    m_queryText.setEnvironment(m_env);
/*    m_queryText.setMinimumSize(new Dimension(width * 2, height));
    m_queryText.setPreferredSize(new Dimension(width * 2, height)); */
    m_queryText.setText(((DatabaseLoader)m_dsLoader.getLoader()).getQuery());
    gbConstraints = new GridBagConstraints();
    gbConstraints.anchor = GridBagConstraints.EAST;
    gbConstraints.fill = GridBagConstraints.HORIZONTAL;
    gbConstraints.gridy = 3; gbConstraints.gridx = 1;
    gbLayout.setConstraints(m_queryText, gbConstraints);
    db.add(m_queryText);
    
    JLabel keyLab = new JLabel("Key columns", SwingConstants.RIGHT);
    keyLab.setBorder(BorderFactory.createEmptyBorder(0, 5, 0, 0));
    gbConstraints = new GridBagConstraints();
    gbConstraints.anchor = GridBagConstraints.EAST;
    gbConstraints.fill = GridBagConstraints.HORIZONTAL;
    gbConstraints.gridy = 4; gbConstraints.gridx = 0;
    gbLayout.setConstraints(keyLab, gbConstraints);
    db.add(keyLab);
    
    m_keyText = new EnvironmentField();
    m_keyText.setEnvironment(m_env);
    /*m_keyText.setMinimumSize(new Dimension(width * 2, height));
    m_keyText.setPreferredSize(new Dimension(width * 2, height)); */
    m_keyText.setText(((DatabaseLoader)m_dsLoader.getLoader()).getKeys());
    gbConstraints = new GridBagConstraints();
    gbConstraints.anchor = GridBagConstraints.EAST;
    gbConstraints.fill = GridBagConstraints.HORIZONTAL;
    gbConstraints.gridy = 4; gbConstraints.gridx = 1;
    gbLayout.setConstraints(m_keyText, gbConstraints);
    db.add(m_keyText);

    JPanel buttonsP = new JPanel();
    buttonsP.setLayout(new FlowLayout());
    JButton ok,cancel;
    buttonsP.add(ok = new JButton("OK"));
    buttonsP.add(cancel=new JButton("Cancel"));
    ok.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent evt){
        ((DatabaseLoader)m_dsLoader.getLoader()).resetStructure();  
        ((DatabaseConverter)m_dsLoader.getLoader()).setUrl(m_dbaseURLText.getText());
        ((DatabaseConverter)m_dsLoader.getLoader()).setUser(m_userNameText.getText());
        ((DatabaseConverter)m_dsLoader.getLoader()).setPassword(new String(m_passwordText.getPassword()));
        ((DatabaseLoader)m_dsLoader.getLoader()).setQuery(m_queryText.getText());
        ((DatabaseLoader)m_dsLoader.getLoader()).setKeys(m_keyText.getText());
        try{
          m_dsLoader.notifyStructureAvailable(((DatabaseLoader)m_dsLoader.getLoader()).getStructure());
          //database connection has been configured
          m_dsLoader.setDB(true);
        }catch (Exception ex){
        }
        if (m_parentFrame != null) {
          m_parentFrame.dispose();
        }
      }
    });
    cancel.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent evt){
        if (m_parentFrame != null) {
          m_parentFrame.dispose();
        }
      }
    });
    
    JPanel holderP = new JPanel();
    holderP.setLayout(new BorderLayout());
    holderP.add(db, BorderLayout.NORTH);
    holderP.add(buttonsP, BorderLayout.SOUTH);
    
    //db.add(buttonsP);
    JPanel about = m_LoaderEditor.getAboutPanel();
    if (about != null) {
      add(about, BorderLayout.NORTH);
    }
    add(holderP, BorderLayout.SOUTH);
  }

  public void setUpFile() {
    removeAll();

    boolean currentFileIsDir = false;
    File tmp = ((FileSourcedConverter)m_dsLoader.getLoader()).retrieveFile();
    String tmpString = tmp.toString();
    if (Environment.containsEnvVariables(tmpString)) {
      try {
        tmpString = m_env.substitute(tmpString);
      } catch (Exception ex) {
        // ignore
      }
    }
    File tmp2 = new File((new File(tmpString)).getAbsolutePath());

    if (tmp2.isDirectory()) {
      m_fileChooser.setCurrentDirectory(tmp2);
      currentFileIsDir = true;
    } else {
      m_fileChooser.setSelectedFile(tmp2);
    }
    
    FileSourcedConverter loader = (FileSourcedConverter) m_dsLoader.getLoader();
    String[] ext = loader.getFileExtensions();
    ExtensionFileFilter firstFilter = null;
    for (int i = 0; i < ext.length; i++) {
      ExtensionFileFilter ff =
	new ExtensionFileFilter(
	    ext[i], loader.getFileDescription() + " (*" + ext[i] + ")");
      if (i == 0)
	firstFilter = ff;
      m_fileChooser.addChoosableFileFilter(ff);
    }
    if (firstFilter != null)
      m_fileChooser.setFileFilter(firstFilter);
    JPanel about = m_LoaderEditor.getAboutPanel();
    JPanel northPanel = new JPanel();
    northPanel.setLayout(new BorderLayout());
    if (about != null) {
      northPanel.add(about, BorderLayout.NORTH);
    }
    add(northPanel, BorderLayout.NORTH);
    final EnvironmentField ef = new EnvironmentField();
    JPanel efHolder = new JPanel();
    efHolder.setLayout(new BorderLayout());

    ef.setEnvironment(m_env);
    /*int width = ef.getPreferredSize().width;
    int height = ef.getPreferredSize().height;
//    ef.setMinimumSize(new Dimension(width * 2, height));
    ef.setPreferredSize(new Dimension(width * 2, height)); */
    m_fileText = ef;
    
    // only set the text on the EnvironmentField if the current file is not a directory
    if (!currentFileIsDir) {
      ef.setText(tmp.toString());
    }
    
    efHolder.add(ef, BorderLayout.CENTER);
    JButton browseBut = new JButton("Browse...");
    browseBut.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        try {
          final JFrame jf = new JFrame("Choose file");
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
    
    JPanel bP = new JPanel(); bP.setLayout(new BorderLayout());
    bP.setBorder(BorderFactory.createEmptyBorder(5,0,5,5));
    bP.add(browseBut, BorderLayout.CENTER);
    efHolder.add(bP, BorderLayout.EAST);
    JPanel alignedP = new JPanel();
    GridBagLayout gbLayout = new GridBagLayout();
    alignedP.setLayout(gbLayout);
    JLabel efLab = new JLabel("Filename", SwingConstants.RIGHT);
    efLab.setBorder(BorderFactory.createEmptyBorder(0, 5, 0, 0));
    GridBagConstraints gbConstraints = new GridBagConstraints();
    gbConstraints.anchor = GridBagConstraints.EAST;
    gbConstraints.fill = GridBagConstraints.HORIZONTAL;
    gbConstraints.gridy = 0; gbConstraints.gridx = 0;
    gbLayout.setConstraints(efLab, gbConstraints);    
    alignedP.add(efLab);
    
    gbConstraints = new GridBagConstraints();
    gbConstraints.anchor = GridBagConstraints.EAST;
    gbConstraints.fill = GridBagConstraints.HORIZONTAL;
    gbConstraints.gridy = 0; gbConstraints.gridx = 1;
    gbConstraints.weightx = 5;
    gbLayout.setConstraints(efHolder, gbConstraints);
    alignedP.add(efHolder);
    
    
    northPanel.add(alignedP, BorderLayout.SOUTH);
    
    // add(m_fileChooser, BorderLayout.CENTER);

    // m_relativeFilePath = new JCheckBox("Use relative file paths");
    m_relativeFilePath = new JCheckBox();
    m_relativeFilePath.
      setSelected(((FileSourcedConverter)m_dsLoader.getLoader()).getUseRelativePath());

    m_relativeFilePath.addActionListener(new ActionListener() {
        public void actionPerformed(ActionEvent e) {
          ((FileSourcedConverter)m_dsLoader.getLoader()).
            setUseRelativePath(m_relativeFilePath.isSelected());
        }
      });
    
    JLabel relativeLab = new JLabel("Use relative file paths", SwingConstants.RIGHT);
    relativeLab.setBorder(BorderFactory.createEmptyBorder(0, 5, 0, 0));
    gbConstraints = new GridBagConstraints();
    gbConstraints.anchor = GridBagConstraints.EAST;
    gbConstraints.fill = GridBagConstraints.HORIZONTAL;
    gbConstraints.gridy = 1; gbConstraints.gridx = 0;
    gbLayout.setConstraints(relativeLab, gbConstraints);
    alignedP.add(relativeLab);
    
    
    gbConstraints = new GridBagConstraints();
    gbConstraints.anchor = GridBagConstraints.EAST;
    gbConstraints.fill = GridBagConstraints.HORIZONTAL;
    gbConstraints.gridy = 1; gbConstraints.gridx = 1;
    gbLayout.setConstraints(m_relativeFilePath, gbConstraints);
    alignedP.add(m_relativeFilePath);
        
    JPanel butHolder = new JPanel();
    //butHolder.setLayout(new GridLayout(1,2));
    butHolder.setLayout(new FlowLayout());
    JButton OKBut = new JButton("OK");
    OKBut.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        try {
          ((FileSourcedConverter)m_dsLoader.getLoader()).
          setFile(new File(ef.getText()));
          // tell the loader that a new file has been selected so
          // that it can attempt to load the header
          m_dsLoader.setLoader(m_dsLoader.getLoader());
          m_dsLoader.newFileSelected();
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
    
    butHolder.add(OKBut);
    butHolder.add(CancelBut);
    
    add(butHolder, BorderLayout.SOUTH);
  }

  /**
   * Set the loader to be customized
   *
   * @param object a weka.gui.beans.Loader
   */
  public void setObject(Object object) {
    m_dsLoader = (weka.gui.beans.Loader)object;
    m_LoaderEditor.setTarget(m_dsLoader.getLoader());
    //    m_fileEditor.setValue(m_dsLoader.getDataSetFile());
    if (m_dsLoader.getLoader() instanceof FileSourcedConverter) {
      setUpFile();
    } else{ 
        if(m_dsLoader.getLoader() instanceof DatabaseConverter) {
            setUpDatabase();
        }
        else
      setUpOther();
    }
  }
  
  public void setEnvironment(Environment env) {
    m_env = env;
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
}
