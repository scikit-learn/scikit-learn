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
 * FileScriptingPanel.java
 * Copyright (C) 2009 University of Waikato, Hamilton, New Zealand
 * Copyright (c) 1995 - 2008 Sun Microsystems, Inc.  
 */

package weka.gui.scripting;

import weka.core.Utils;
import weka.gui.ComponentHelper;
import weka.gui.DocumentPrinting;
import weka.gui.ExtensionFileFilter;
import weka.gui.PropertyDialog;
import weka.gui.scripting.event.ScriptExecutionEvent;
import weka.gui.scripting.event.ScriptExecutionListener;
import weka.gui.scripting.event.TitleUpdatedEvent;
import weka.gui.scripting.event.ScriptExecutionEvent.Type;

import java.awt.BorderLayout;
import java.awt.Dialog;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.io.File;
import java.util.HashMap;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.BorderFactory;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JInternalFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JTextPane;
import javax.swing.KeyStroke;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.event.UndoableEditEvent;
import javax.swing.event.UndoableEditListener;
import javax.swing.text.DefaultEditorKit;
import javax.swing.text.Document;
import javax.swing.text.JTextComponent;
import javax.swing.undo.CannotRedoException;
import javax.swing.undo.CannotUndoException;
import javax.swing.undo.UndoManager;

/**
 * Supports loading/saving of files.
 * 
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @author  Sun Microsystems Inc (see <a href="http://java.sun.com/docs/books/tutorial/uiswing/examples/components/TextComponentDemoProject/src/components/TextComponentDemo.java">TextComponentDemo.java</a>)
 * @version $Revision: 5144 $
 */
public abstract class FileScriptingPanel
  extends ScriptingPanel 
  implements ScriptExecutionListener {

  /** for serialization. */
  private static final long serialVersionUID = 1583670545010241816L;

  /**
   * A slightly extended action class.
   * 
   * @author  fracpete (fracpete at waikato dot ac dot nz)
   * @version $Revision: 5144 $
   */
  public abstract class BasicAction
    extends AbstractAction {
    
    /** for serialization. */
    private static final long serialVersionUID = 2821117985661550385L;

    /**
     * Constructor for setting up an action.
     * 
     * @param name	the name of the action (to be displayed in menu/button)
     * @param icon	the icon name (no path required if in weka/gui/images), can be null
     * @param accel	the accelerator command, e.g., "ctrl N", can be null
     * @param mnemonic	the mnemonic character
     */
    public BasicAction(String name, String icon, String accel, Character mnemonic) {
      super(name);
      
      if ((icon != null) && (icon.length() > 0))
	putValue(Action.SMALL_ICON, ComponentHelper.getImageIcon(icon));
      if ((accel != null) && (accel.length() > 0))
	putValue(Action.ACCELERATOR_KEY, KeyStroke.getKeyStroke(accel));
      if (mnemonic != null)
	putValue(Action.MNEMONIC_KEY, new Integer(mnemonic.charValue()));
    }
  }
  
  /**
   * The New action.
   * 
   * @author  fracpete (fracpete at waikato dot ac dot nz)
   * @version $Revision: 5144 $
   */
  protected class NewAction
    extends BasicAction {
    
    /** for serialization. */
    private static final long serialVersionUID = -8665722554539726090L;

    /**
     * Initializes the action.
     */
    public NewAction() {
      super("New", "new.gif", "ctrl N", 'N');
      setEnabled(true);
    }

    /**
     * Fired when action got executed.
     * 
     * @param e		the event
     */
    public void actionPerformed(ActionEvent e) {
      m_Script.empty();
      notifyTitleUpdatedListeners(new TitleUpdatedEvent(FileScriptingPanel.this));
    }
  }

  /**
   * The Open action.
   * 
   * @author  fracpete (fracpete at waikato dot ac dot nz)
   * @version $Revision: 5144 $
   */
  protected class OpenAction
    extends BasicAction {
    
    /** for serialization. */
    private static final long serialVersionUID = -4496148485267789162L;

    /**
     * Initializes the action.
     */
    public OpenAction() {
      super("Open...", "open.gif", "ctrl O", 'O');
      setEnabled(true);
    }

    /**
     * Fired when action got executed.
     * 
     * @param e		the event
     */
    public void actionPerformed(ActionEvent e) {
      boolean	ok;
      int	retVal;
      
      if (!checkModified())
	return;
      
      retVal = m_FileChooser.showOpenDialog(FileScriptingPanel.this);
      if (retVal != JFileChooser.APPROVE_OPTION)
        return;
      
      ok = m_Script.open(m_FileChooser.getSelectedFile());
      m_TextCode.setCaretPosition(0);
      if (!ok)
	JOptionPane.showMessageDialog(
	    FileScriptingPanel.this, 
	    "Couldn't open file '" + m_FileChooser.getSelectedFile() + "'!");
      
      notifyTitleUpdatedListeners(new TitleUpdatedEvent(FileScriptingPanel.this));
    }
  }

  /**
   * The Save action.
   * 
   * @author  fracpete (fracpete at waikato dot ac dot nz)
   * @version $Revision: 5144 $
   */
  protected class SaveAction
    extends BasicAction {
    
    /** for serialization. */
    private static final long serialVersionUID = -74651145892063975L;
    
    /** whether to bring up the save dialog all the time. */
    protected boolean m_ShowDialog;
    
    /**
     * Initializes the action.
     * 
     * @param name		the name of the action
     * @param showDialog	whether to always show the dialog
     */
    public SaveAction(String name, boolean showDialog) {
      super(name, (showDialog ? "" : "save.gif"), (showDialog ? "ctrl shift S" : "ctrl S"), (showDialog ? 'a' : 'S'));
      m_ShowDialog = showDialog;
      setEnabled(true);
    }

    /**
     * Fired when action got executed.
     * 
     * @param e		the event
     */
    public void actionPerformed(ActionEvent e) {
      boolean	ok;
      int	retVal;
      
      if (m_ShowDialog || (m_Script.getFilename() == null)) {
	retVal = m_FileChooser.showSaveDialog(FileScriptingPanel.this);
	if (retVal != JFileChooser.APPROVE_OPTION)
	  return;
	ok = m_Script.saveAs(m_FileChooser.getSelectedFile());
      }
      else {
	ok = m_Script.save();
      }
      
      if (!ok) {
	if (m_Script.getFilename() != null)
	  JOptionPane.showMessageDialog(
	      FileScriptingPanel.this, 
	      "Failed to save file '" + m_FileChooser.getSelectedFile() + "'!");
	else
	  JOptionPane.showMessageDialog(
	      FileScriptingPanel.this, 
	      "Failed to save file!");
      }
      else {
	m_SaveAction.setEnabled(false);
      }

      notifyTitleUpdatedListeners(new TitleUpdatedEvent(FileScriptingPanel.this));
    }
  }

  /**
   * The Print action.
   * 
   * @author  fracpete (fracpete at waikato dot ac dot nz)
   * @version $Revision: 5144 $
   */
  protected class PrintAction
    extends BasicAction {
    
    /** for serialization. */
    private static final long serialVersionUID = -6246539539545724632L;

    /**
     * Initializes the action.
     */
    public PrintAction() {
      super("Print...", "print.gif", "ctrl P", 'P');
      setEnabled(true);
    }

    /**
     * Fired when action got executed.
     * 
     * @param e		the event
     */
    public void actionPerformed(ActionEvent e) {
      JTextPane		pane;
      DocumentPrinting 	doc;
      
      pane = newCodePane();
      pane.setText(m_TextCode.getText());
      doc = new DocumentPrinting();
      doc.print(pane);
    }
  }

  /**
   * The Clear output action.
   * 
   * @author  fracpete (fracpete at waikato dot ac dot nz)
   * @version $Revision: 5144 $
   */
  protected class ClearOutputAction
    extends BasicAction {
    
    /** for serialization. */
    private static final long serialVersionUID = 47986890456997211L;

    /**
     * Initializes the action.
     */
    public ClearOutputAction() {
      super("Clear output", "", "F2", 'C');
      setEnabled(true);
    }

    /**
     * Fired when action got executed.
     * 
     * @param e		the event
     */
    public void actionPerformed(ActionEvent e) {
      m_TextOutput.setText("");
    }
  }

  /**
   * The Exit action. Sends out a WindowEvent/WINDOW_CLOSED to all 
   * WindowListener objects of a jframe.
   * 
   * @author  fracpete (fracpete at waikato dot ac dot nz)
   * @version $Revision: 5144 $
   */
  protected class ExitAction
    extends BasicAction {
    
    /** for serialization. */
    private static final long serialVersionUID = -5884709836238884180L;

    /**
     * Initializes the action.
     */
    public ExitAction() {
      super("Exit", "", "", 'x');
      setEnabled(true);
    }

    /**
     * Fired when action got executed.
     * 
     * @param e		the event
     */
    public void actionPerformed(ActionEvent e) {
      Dialog		dialog;
      Frame		frame;
      JFrame		jframe;
      JInternalFrame	jintframe;
      int		i;
      WindowListener[] 	listeners;
      WindowEvent	event;
      
      if (!checkModified())
	return;
      
      if (PropertyDialog.getParentDialog(FileScriptingPanel.this) != null) {
	dialog = PropertyDialog.getParentDialog(FileScriptingPanel.this);
	dialog.setVisible(false);
      }
      else if (PropertyDialog.getParentFrame(FileScriptingPanel.this) != null) {
	jintframe = PropertyDialog.getParentInternalFrame(FileScriptingPanel.this);
	if (jintframe != null) {
	  jintframe.doDefaultCloseAction();
	}
	else {
	  frame = PropertyDialog.getParentFrame(FileScriptingPanel.this);
	  if (frame instanceof JFrame) {
	    jframe = (JFrame) frame;
	    if (jframe.getDefaultCloseOperation() == JFrame.HIDE_ON_CLOSE)
	      jframe.setVisible(false);
	    else if (jframe.getDefaultCloseOperation() == JFrame.DISPOSE_ON_CLOSE)
	      jframe.dispose();
	    else if (jframe.getDefaultCloseOperation() == JFrame.EXIT_ON_CLOSE)
	      System.exit(0);

	    // notify listeners
	    listeners = jframe.getWindowListeners();
	    event     = new WindowEvent(jframe, WindowEvent.WINDOW_CLOSED);
	    for (i = 0; i < listeners.length; i++)
	      listeners[i].windowClosed(event);
	  }
	  else {
	    frame.dispose();
	}
	}
      }
    }
  }

  /**
   * The Undo action.
   * 
   * @author  Sun Microsystems Inc (see <a href="http://java.sun.com/docs/books/tutorial/uiswing/examples/components/TextComponentDemoProject/src/components/TextComponentDemo.java">TextComponentDemo.java</a>)
   * @version $Revision: 5144 $
   */
  protected class UndoAction 
    extends BasicAction {
    
    /** for serialization. */
    private static final long serialVersionUID = 4298096648424808522L;

    /**
     * Initializes the action.
     */
    public UndoAction() {
      super("Undo", "undo.gif", "ctrl Z", 'U');
      setEnabled(false);
    }

    /**
     * Fired when action got executed.
     * 
     * @param e		the event
     */
    public void actionPerformed(ActionEvent e) {
      try {
	m_Undo.undo();
      }
      catch (CannotUndoException ex) {
	System.out.println("Unable to undo: " + ex);
	ex.printStackTrace();
      }
      updateUndoState();
      m_RedoAction.updateRedoState();
    }

    /**
     * Updates the redo state.
     */
    protected void updateUndoState() {
      if (m_Undo.canUndo()) {
	setEnabled(true);
	putValue(Action.NAME, m_Undo.getUndoPresentationName());
      }
      else {
	setEnabled(false);
	putValue(Action.NAME, "Undo");
      }
    }
  }

  /**
   * The Redo action.
   * 
   * @author  Sun Microsystems Inc (see <a href="http://java.sun.com/docs/books/tutorial/uiswing/examples/components/TextComponentDemoProject/src/components/TextComponentDemo.java">TextComponentDemo.java</a>)
   * @version $Revision: 5144 $
   */
  protected class RedoAction
    extends BasicAction {
    
    /** for serialization. */
    private static final long serialVersionUID = 4533966901523279350L;

    /**
     * Initializes the action.
     */
    public RedoAction() {
      super("Redo", "redo.gif", "ctrl Y", 'R');
      setEnabled(false);
    }

    /**
     * Fired when action got executed.
     * 
     * @param e		the event
     */
    public void actionPerformed(ActionEvent e) {
      try {
	m_Undo.redo();
      }
      catch (CannotRedoException ex) {
	System.out.println("Unable to redo: " + ex);
	ex.printStackTrace();
      }
      updateRedoState();
      m_UndoAction.updateUndoState();
    }

    /**
     * Updates the redo state.
     */
    protected void updateRedoState() {
      if (m_Undo.canRedo()) {
	setEnabled(true);
	putValue(Action.NAME, m_Undo.getRedoPresentationName());
      }
      else {
	setEnabled(false);
	putValue(Action.NAME, "Redo");
      }
    }
  }

  /**
   * The Commandline args action.
   * 
   * @author  fracpete (fracpete at waikato dot ac dot nz)
   * @version $Revision: 5144 $
   */
  protected class CommandlineArgsAction
    extends BasicAction {
    
    /** for serialization. */
    private static final long serialVersionUID = -3183470039010826204L;

    /**
     * Initializes the action.
     */
    public CommandlineArgsAction() {
      super("Arguments...", "properties.gif", "", 'g');
      setEnabled(true);
    }

    /**
     * Fired when action got executed.
     * 
     * @param e		the event
     */
    public void actionPerformed(ActionEvent e) {
      String	retVal;
      
      retVal = JOptionPane.showInputDialog(
	  FileScriptingPanel.this, 
	  "Please enter the command-line arguments", 
	  Utils.joinOptions(m_Args));
      if (retVal == null)
	return;
      
      try {
	m_Args = Utils.splitOptions(retVal);
      }
      catch (Exception ex) {
	m_Args = new String[0];
	ex.printStackTrace();
	JOptionPane.showMessageDialog(
	    FileScriptingPanel.this, 
	    "Error setting command-line arguments:\n" + ex, 
	    "Error", 
	    JOptionPane.ERROR_MESSAGE);
      }
    }
  }

  /**
   * The Start action.
   * 
   * @author  fracpete (fracpete at waikato dot ac dot nz)
   * @version $Revision: 5144 $
   */
  protected class StartAction
    extends BasicAction {
    
    /** for serialization. */
    private static final long serialVersionUID = -7936456072955996220L;

    /**
     * Initializes the action.
     */
    public StartAction() {
      super((m_Script.canExecuteScripts() ? "Start" : "Start (missing classes?)"), "run.gif", "ctrl R", 'S');
      setEnabled(false);
    }

    /**
     * Fired when action got executed.
     * 
     * @param e		the event
     */
    public void actionPerformed(ActionEvent e) {
      if (!checkModified())
	return;
      
      if (m_Script.getFilename() == null)
	return;
      
      try {
	m_Script.start(m_Args);
      }
      catch (Exception ex) {
	ex.printStackTrace();
	JOptionPane.showMessageDialog(
	    FileScriptingPanel.this, 
	    "Error running script:\n" + ex, 
	    "Error", 
	    JOptionPane.ERROR_MESSAGE);
      }
    }
  }

  /**
   * The Stop action.
   * 
   * @author  fracpete (fracpete at waikato dot ac dot nz)
   * @version $Revision: 5144 $
   */
  protected class StopAction
    extends BasicAction {
    
    /** for serialization. */
    private static final long serialVersionUID = 8764023289575718872L;

    /**
     * Initializes the action.
     */
    public StopAction() {
      super("Stop", "stop.gif", "ctrl shift R", 'o');
      setEnabled(false);
    }

    /**
     * Fired when action got executed.
     * 
     * @param e		the event
     */
    public void actionPerformed(ActionEvent e) {
      try {
	m_Script.stop();
      }
      catch (Exception ex) {
	// ignored
      }
    }
  }

  /**
   * The About action.
   * 
   * @author  fracpete (fracpete at waikato dot ac dot nz)
   * @version $Revision: 5144 $
   */
  protected class AboutAction
    extends BasicAction {
    
    /** for serialization. */
    private static final long serialVersionUID = -6420463480569171227L;

    /**
     * Initializes the action.
     */
    public AboutAction() {
      super("About...", "", "F1", 'A');
      setEnabled(true);
    }

    /**
     * Fired when action got executed.
     * 
     * @param e		the event
     */
    public void actionPerformed(ActionEvent e) {
      JDialog	dialog;
      
      if (PropertyDialog.getParentDialog(FileScriptingPanel.this) != null)
	dialog = new JDialog(PropertyDialog.getParentDialog(FileScriptingPanel.this), getName());
      else
	dialog = new JDialog(PropertyDialog.getParentFrame(FileScriptingPanel.this), getName());
      dialog.setTitle((String) getValue(Action.NAME));
      dialog.getContentPane().setLayout(new BorderLayout());
      dialog.getContentPane().add(getAboutPanel());
      dialog.pack();
      dialog.setLocationRelativeTo(FileScriptingPanel.this);
      dialog.setVisible(true);
    }
  }
  
  /**
   * This listener class listens for edits that can be undone.
   * 
   * @author  Sun Microsystems Inc (see <a href="http://java.sun.com/docs/books/tutorial/uiswing/examples/components/TextComponentDemoProject/src/components/TextComponentDemo.java">TextComponentDemo.java</a>)
   * @version $Revision: 5144 $
   */
  protected class ScriptUndoableEditListener
    implements UndoableEditListener {
    
    /**
     * Gets called when an undoable event gets triggered.
     * 
     * @param e		the event
     */
    public void undoableEditHappened(UndoableEditEvent e) {
      //Remember the edit and update the menus.
      m_Undo.addEdit(e.getEdit());
      m_UndoAction.updateUndoState();
      m_RedoAction.updateRedoState();
    }
  }
  
  /** the directory with the scripting-specific images. */
  public final static String IMAGES_DIR = "weka/gui/scripting/images";
  
  /** for loading/saving file. */
  protected JFileChooser m_FileChooser;
  
  /** the script. */
  protected Script m_Script;

  /** the script area. */
  protected JTextArea m_ScriptArea;

  /** the output area. */
  protected JTextArea m_OutputArea;
  
  /** for informing the user. */
  protected JLabel m_LabelInfo;
  
  /** for storing the actions under their name. */
  protected HashMap<Object, Action> m_Actions;

  /** the new action. */
  protected NewAction m_NewAction;
  /** the open action. */
  protected OpenAction m_OpenAction;
  /** the Save action. */
  protected SaveAction m_SaveAction;
  /** the Save as action. */
  protected SaveAction m_SaveAsAction;
  /** the Print action. */
  protected PrintAction m_PrintAction;
  /** the clear output action. */
  protected ClearOutputAction m_ClearOutputAction;
  /** the exit action. */
  protected ExitAction m_ExitAction;
  /** the undo action. */
  protected UndoAction m_UndoAction;
  /** the redo action. */
  protected RedoAction m_RedoAction;
  /** the cut action. */
  protected Action m_CutAction;
  /** the copy action. */
  protected Action m_CopyAction;
  /** the paste action. */
  protected Action m_PasteAction;
  /** the start action. */
  protected StartAction m_StartAction;
  /** the stop action. */
  protected StopAction m_StopAction;
  /** the arguments action. */
  protected CommandlineArgsAction m_ArgsAction;
  /** the about action. */
  protected AboutAction m_AboutAction;

  /** the undo manager. */
  protected UndoManager m_Undo;
  
  /** the text pane with the code. */
  protected JTextPane m_TextCode;
  
  /** the text pane for the output. */
  protected JTextPane m_TextOutput;
  
  /** the commandline arguments to use. */
  protected String[] m_Args;
  
  /**
   * For initializing member variables.
   */
  protected void initialize() {
    super.initialize();
    
    m_FileChooser = new JFileChooser();
    m_FileChooser.setAcceptAllFileFilterUsed(true);
    m_FileChooser.setMultiSelectionEnabled(false);
    
    m_Undo = new UndoManager();
    m_Args = new String[0];
  }
  
  /**
   * Sets up the GUI after initializing the members.
   */
  protected void initGUI() {
    JPanel	panel;
    
    super.initGUI();

    setLayout(new BorderLayout(0, 5));
    
    m_TextCode = newCodePane();
    m_TextCode.setFont(new Font("monospaced", Font.PLAIN, 12));
    m_TextCode.getDocument().addUndoableEditListener(new ScriptUndoableEditListener());
    m_TextCode.getDocument().addDocumentListener(new DocumentListener() {
      public void changedUpdate(DocumentEvent e) {
	update();
      }
      public void insertUpdate(DocumentEvent e) {
	update();
      }
      public void removeUpdate(DocumentEvent e) {
	update();
      }
      protected void update() {
	Document doc = m_TextCode.getDocument();
	m_StartAction.setEnabled((doc.getLength() > 0) && m_Script.canExecuteScripts());
	m_SaveAction.setEnabled(true);
	notifyTitleUpdatedListeners(new TitleUpdatedEvent(FileScriptingPanel.this));
      }
    });
    add(new JScrollPane(m_TextCode), BorderLayout.CENTER);
    
    panel = new JPanel(new BorderLayout(0, 5));
    panel.setPreferredSize(new Dimension(50, 200));
    add(panel, BorderLayout.SOUTH);
    
    m_TextOutput = new JTextPane();
    panel.add(new JScrollPane(m_TextOutput), BorderLayout.CENTER);
    
    m_LabelInfo = new JLabel(" ");
    m_LabelInfo.setBorder(BorderFactory.createLoweredBevelBorder());
    panel.add(m_LabelInfo, BorderLayout.SOUTH);
  }
  
  /**
   * Finishes up after initializing members and setting up the GUI.
   */
  protected void initFinish() {
    ExtensionFileFilter[]	filters;
    int				i;
    
    super.initFinish();
    
    m_Script = newScript(m_TextCode.getDocument());
    m_Script.addScriptFinishedListener(this);
    filters  = m_Script.getFilters();
    for (i = filters.length - 1; i >= 0; i--)
      m_FileChooser.addChoosableFileFilter(filters[i]);
    
    m_Actions = createActionTable(m_TextCode);

    // file
    m_NewAction         = new NewAction();
    m_OpenAction        = new OpenAction();
    m_SaveAction        = new SaveAction("Save", false);
    m_SaveAsAction      = new SaveAction("Save As...", true);
    m_PrintAction       = new PrintAction();
    m_ClearOutputAction = new ClearOutputAction();
    m_ExitAction        = new ExitAction();
    
    // edit
    m_UndoAction        = new UndoAction();
    m_RedoAction        = new RedoAction();
    m_CutAction         = updateAction(m_Actions.get(DefaultEditorKit.cutAction), "Cut", "cut.gif", "ctrl X", 'C');
    m_CopyAction        = updateAction(m_Actions.get(DefaultEditorKit.copyAction), "Copy", "copy.gif", "ctrl C", 'o');
    m_PasteAction       = updateAction(m_Actions.get(DefaultEditorKit.pasteAction), "Paste", "paste.gif", "ctrl V", 'P');
    
    // script
    m_StartAction       = new StartAction();
    m_StopAction        = new StopAction();
    m_ArgsAction        = new CommandlineArgsAction();
    
    // help
    m_AboutAction       = new AboutAction();
  }
  
  /**
   * Updates the action and returns it.
   * 
   * @param action	the action to update
   * @param name	the name to be used as display, can be null
   * @param icon	the icon to use (if located in weka/gui/images, not path required), can be null
   * @param accel	the accelerator command to use (e.g., "ctrl N"), can be null
   * @param mnemonic	the mnemonic character to use, can be null
   * @return		the updated action
   */
  protected Action updateAction(Action action, String name, String icon, String accel, Character mnemonic) {
    Action	result;
    
    // did we already update that action for another component?
    if (action == null) {
      result = m_Actions.get(name);
      return result;
    }
    
    result = action;
    
    if ((name != null) && (name.length() > 0))
      result.putValue(Action.NAME, name);
    if ((icon != null) && (icon.length() > 0))
      result.putValue(Action.SMALL_ICON, ComponentHelper.getImageIcon(icon));
    if ((accel != null) && (accel.length() > 0))
      result.putValue(Action.ACCELERATOR_KEY, KeyStroke.getKeyStroke(accel));
    if (mnemonic != null)
      result.putValue(Action.MNEMONIC_KEY, new Integer(mnemonic.charValue()));
    
    return result;
  }
  
  /**
   * Creates a new JTextPane for the code.
   * 
   * @return		the text pane
   */
  protected abstract JTextPane newCodePane();

  /**
   * Returns an initialized script object.
   * 
   * @param doc		the document to use as basis
   * @return		the initialized script
   */
  protected abstract Script newScript(Document doc);

  /**
   * Gets sent when a script finishes execution.
   * 
   * @param e		the event
   */
  public void scriptFinished(ScriptExecutionEvent e) {
    if (e.getType() == Type.FINISHED)
      showInfo("Script execution finished");
    else if (e.getType() == Type.STOPPED)
      showInfo("Script execution stopped by user");
    else if (e.getType() == Type.ERROR)
      showInfo("Script execution failed" + (e.hasAdditional() ? (": " + e.getAdditional()) : ""));
    
    if (e.getType() != Type.STARTED) {
      m_NewAction.setEnabled(true);
      m_OpenAction.setEnabled(true);
      m_SaveAction.setEnabled(true);
      m_SaveAsAction.setEnabled(true);
      m_CutAction.setEnabled(true);
      m_CopyAction.setEnabled(true);
      m_PasteAction.setEnabled(true);
      m_StartAction.setEnabled(true);
      m_StopAction.setEnabled(false);
    }
    else {
      m_NewAction.setEnabled(false);
      m_OpenAction.setEnabled(false);
      m_SaveAction.setEnabled(false);
      m_SaveAsAction.setEnabled(false);
      m_CutAction.setEnabled(false);
      m_CopyAction.setEnabled(false);
      m_PasteAction.setEnabled(false);
      m_StartAction.setEnabled(false);
      m_StopAction.setEnabled(true);
    }
  }
  
  /**
   * The following two methods allow us to find an
   * action provided by the editor kit by its name.
   * 
   * @param comp	the component to get the actions from
   * @return		the relation
   */
  protected HashMap<Object, Action> createActionTable(JTextComponent comp) {
    HashMap<Object, Action> 	result;
    int				i;
    Action[] 			actions;
    Action 			action;
    
    result  = new HashMap<Object, Action>();
    actions = comp.getActions();
    for (i = 0; i < actions.length; i++) {
      action = actions[i];
      result.put(action.getValue(Action.NAME), action);
    }
    
    return result;
  }
  
  /**
   * Returns a panel to be displayed with the AboutAction.
   * 
   * @return		the panel with some information on the scripting panel
   */
  protected abstract JPanel getAboutPanel();
  
  /**
   * Returns the title (without the filename).
   * 
   * @return		the plain title
   */
  public abstract String getPlainTitle();
  
  /**
   * Returns the current title for the frame/dialog.
   * 
   * @return		the title
   */
  public String getTitle() {
    String	result;
    
    result = getPlainTitle();
    
    if (m_Script.isModified())
      result = "*" + result;
    
    if (m_Script.getFilename() != null)
      result += " [" + m_Script.getFilename() + "]";
    
    return result;
  }
  
  /**
   * Returns the text area that is used for displaying output on stdout
   * and stderr.
   * 
   * @return		the JTextArea
   */
  public JTextPane getOutput() {
    return m_TextOutput;
  }
  
  /**
   * Returns the menu bar to to be displayed in the frame.
   * 
   * @return		the menu bar, null if not applicable
   */
  public JMenuBar getMenuBar() {
    JMenuBar	result;
    JMenu	menu;
    JMenuItem	menuitem;
    
    result = new JMenuBar();
    
    // File
    menu = new JMenu("File");
    menu.setMnemonic('F');
    result.add(menu);
    
    // File/New
    menuitem = new JMenuItem(m_NewAction);
    menu.add(menuitem);
    
    // File/Open
    menuitem = new JMenuItem(m_OpenAction);
    menu.addSeparator();
    menu.add(menuitem);
    
    // File/Save
    menuitem = new JMenuItem(m_SaveAction);
    menu.add(menuitem);
    
    // File/SaveAs
    menuitem = new JMenuItem(m_SaveAsAction);
    menu.add(menuitem);
    
    // File/Print
    menuitem = new JMenuItem(m_PrintAction);
    menu.addSeparator();
    menu.add(menuitem);
    
    // File/Clear output
    menuitem = new JMenuItem(m_ClearOutputAction);
    menu.addSeparator();
    menu.add(menuitem);
    
    // File/Exit
    menuitem = new JMenuItem(m_ExitAction);
    menu.addSeparator();
    menu.add(menuitem);
    
    // Edit
    menu = new JMenu("Edit");
    menu.setMnemonic('E');
    result.add(menu);
    
    // Edit/Undo
    menuitem = new JMenuItem(m_UndoAction);
    menu.add(menuitem);
    
    // Edit/Redo
    menuitem = new JMenuItem(m_RedoAction);
    menu.add(menuitem);
    
    // Edit/Cut
    menuitem = new JMenuItem(m_CutAction);
    menu.addSeparator();
    menu.add(menuitem);
    
    // Edit/Copy
    menuitem = new JMenuItem(m_CopyAction);
    menu.add(menuitem);
    
    // Edit/Paste
    menuitem = new JMenuItem(m_PasteAction);
    menu.add(menuitem);
    
    // Script
    menu = new JMenu("Script");
    menu.setMnemonic('S');
    result.add(menu);
    
    // Script/Start
    menuitem = new JMenuItem(m_StartAction);
    menu.add(menuitem);
    
    // Script/Stop
    menuitem = new JMenuItem(m_StopAction);
    menu.add(menuitem);
    
    // Script/Arguments
    menuitem = new JMenuItem(m_ArgsAction);
    menu.add(menuitem);
    
    // Help
    menu = new JMenu("Help");
    menu.setMnemonic('H');
    result.add(menu);
    
    // Help/About
    menuitem = new JMenuItem(m_AboutAction);
    menu.add(menuitem);
    
    return result;
  }
  
  /**
   * Updates the info shown in the bottom panel.
   * 
   * @param msg		the message to display
   */
  protected void showInfo(String msg) {
    if (msg == null)
      msg = " ";
    m_LabelInfo.setText(msg);
  }
  
  /**
   * Opens the specified file.
   * 
   * @param file	the file to open
   */
  public void open(File file) {
    m_Script.open(file);
  }
  
  /**
   * Checks whether the script is modified and asks the user to save it or not.
   * If everything is fine and one can ignore the modified state, true is 
   * returned.
   * 
   * @return		true if one can proceed
   */
  protected boolean checkModified() {
    boolean	result;
    int		retVal;
    
    result = true;
    
    if (m_Script.isModified()) {
      retVal = JOptionPane.showConfirmDialog(
	  FileScriptingPanel.this, 
	  "Script not saved - save it now?", 
	  "Confirm", 
	  JOptionPane.YES_NO_CANCEL_OPTION);
      
      if (retVal == JOptionPane.YES_OPTION) {
	if (m_Script.getFilename() != null)
	  m_Script.save();
	else
	  m_SaveAsAction.actionPerformed(null);
	result = !m_Script.isModified();
      }
      else if (retVal == JOptionPane.CANCEL_OPTION) {
	result = false;
      }
    }
    
    return result;
  }
}
