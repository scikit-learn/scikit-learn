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
 * LogWindow.java
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui;

import weka.core.Tee;
import weka.core.Utils;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.PrintStream;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JTextPane;
import javax.swing.SpinnerNumberModel;
import javax.swing.event.CaretEvent;
import javax.swing.event.CaretListener;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.text.Style;
import javax.swing.text.StyleConstants;
import javax.swing.text.StyleContext;
import javax.swing.text.StyledDocument;

/** 
 * Frame that shows the output from stdout and stderr.
 *
 * @author  FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 4973 $
 */
public class LogWindow 
  extends JFrame
  implements CaretListener, ChangeListener {

  /** for serialization */
  private static final long serialVersionUID = 5650947361381061112L;

  /** the name of the style for stdout */
  public final static String STYLE_STDOUT = "stdout";

  /** the name of the style for stderr */
  public final static String STYLE_STDERR = "stderr";

  /** the color of the style for stdout */
  public final static Color COLOR_STDOUT = Color.BLACK;

  /** the Color of the style for stderr */
  public final static Color COLOR_STDERR = Color.RED;

  /** whether we're debugging - enables output on stdout */
  public final static boolean DEBUG = false;

  /** whether the JTextPane has wordwrap or not */
  public boolean m_UseWordwrap = true;
  
  /** the output */
  protected JTextPane m_Output = new JTextPane();

  /** the clear button */
  protected JButton m_ButtonClear = new JButton("Clear");

  /** the close button */
  protected JButton m_ButtonClose = new JButton("Close");

  /** the current size */
  protected JLabel m_LabelCurrentSize = new JLabel("currently: 0");

  /** the spinner for the max number of chars */
  protected JSpinner m_SpinnerMaxSize = new JSpinner();

  /** whether to allow wordwrap or not */
  protected JCheckBox m_CheckBoxWordwrap = new JCheckBox("Use wordwrap");

  /** for redirecting stdout */
  protected static Tee m_TeeOut = null;

  /** for redirecting stderr */
  protected static Tee m_TeeErr = null;

  /** inner class for printing to the window, is used instead of standard
   * System.out and System.err */
  protected class LogWindowPrintStream extends PrintStream {
    /** the parent */
    protected LogWindow m_Parent = null;

    /** the style of the printstream */
    protected String m_Style = null;
    
    /**
     * the constructor
     * @param parent      the parent frame
     * @param stream      the stream (used for constructor of superclass)
     * @param style       the style name associated with this output
     */
    public LogWindowPrintStream( LogWindow parent, 
                                 PrintStream stream, 
                                 String style ) {
      super(stream);

      m_Parent = parent;
      m_Style  = style;
    }
    
    /**
     * flushes the printstream
     */
    public synchronized void flush() {
      // ignored
    }

    /**
     * prints the given int
     */
    public synchronized void print(int x) {
      print(new Integer(x).toString());
    }

    /**
     * prints the given boolean 
     */
    public synchronized void print(boolean x) {
      print(new Boolean(x).toString());
    }

    /**
     * prints the given string 
     */
    public synchronized void print(String x) {
      StyledDocument      doc;
      int                 size;
      int                 maxSize;
      int                 pos;
      
      doc = m_Parent.m_Output.getStyledDocument();

      try {
        // insert text
        doc.insertString(doc.getLength(), x, doc.getStyle(m_Style));
        
        // move cursor to end
        m_Parent.m_Output.setCaretPosition(doc.getLength());

        // trim size if necessary
        m_Parent.trim();
      }
      catch (Exception e) {
        e.printStackTrace();
      }
    }

    /**
     * prints the given object
     */
    public synchronized void print(Object x) {
      String                  line;
      Throwable               t;
      StackTraceElement[]     trace;
      int                     i;

      if (x instanceof Throwable) {
        t     = (Throwable) x;
        trace = t.getStackTrace();
        line  = t.getMessage() + "\n";
        for (i = 0; i < trace.length; i++)
          line += "\t" + trace[i].toString() + "\n";
        x = line;
      }

      if (x == null)
	print("null");
      else
	print(x.toString());
    }

    /**
     * prints a new line
     */
    public synchronized void println() {
      print("\n");
    }

    /**
     * prints the given int
     */
    public synchronized void println(int x) {
      print(x);
      println();
    }

    /**
     * prints the given boolean
     */
    public synchronized void println(boolean x) {
      print(x);
      println();
    }

    /**
     * prints the given string
     */
    public synchronized void println(String x) {
      print(x);
      println();
    }

    /**
     * prints the given object (for Throwables we print the stack trace)
     */
    public synchronized void println(Object x) {
      print(x);
      println();
    }
  }
  
  /**
   * creates the frame
   */
  public LogWindow() {
    super("Weka - Log");

    createFrame();

    // styles
    StyledDocument      doc;
    Style               style;
    boolean             teeDone;

    doc   = m_Output.getStyledDocument();
    style = StyleContext.getDefaultStyleContext()
                        .getStyle(StyleContext.DEFAULT_STYLE);
    style = doc.addStyle(STYLE_STDOUT, style);
    StyleConstants.setFontFamily(style, "monospaced");
    StyleConstants.setForeground(style, COLOR_STDOUT);
    
    style = StyleContext.getDefaultStyleContext()
                        .getStyle(StyleContext.DEFAULT_STYLE);
    style = doc.addStyle(STYLE_STDERR, style);
    StyleConstants.setFontFamily(style, "monospaced");
    StyleConstants.setForeground(style, COLOR_STDERR);

    // print streams (instantiate only once!)
    teeDone = !((m_TeeOut == null) && (m_TeeErr == null));
    if (!DEBUG) {
      if (!teeDone) {
        m_TeeOut = new Tee(System.out);
        System.setOut(m_TeeOut);
      }
      m_TeeOut.add(
          new LogWindowPrintStream(this, m_TeeOut.getDefault(), STYLE_STDOUT));
    }

    if (!teeDone) {
      m_TeeErr = new Tee(System.err);
      System.setErr(m_TeeErr);
    }
    m_TeeErr.add(
        new LogWindowPrintStream(this, m_TeeErr.getDefault(), STYLE_STDERR));
  }

  /**
   * creates the frame and all its components
   */
  protected void createFrame() {
    JPanel                panel;
    JPanel                panel2;
    JPanel                panel3;
    JPanel                panel4;
    SpinnerNumberModel    model;
    int                   width;
    JLabel                label;

    // set layout
    setSize(600, 400);
    width = getBounds().width;
    setLocation(
        getGraphicsConfiguration().getBounds().width - width, getLocation().y);
    getContentPane().setLayout(new BorderLayout());
    
    // output 
    getContentPane().add(new JScrollPane(m_Output), BorderLayout.CENTER);
    setWordwrap(m_UseWordwrap);
    
    // button(s)
    panel = new JPanel(new BorderLayout());
    getContentPane().add(panel, BorderLayout.SOUTH);
    panel3 = new JPanel(new BorderLayout());
    panel.add(panel3, BorderLayout.SOUTH);
    panel2 = new JPanel(new FlowLayout(FlowLayout.RIGHT));
    panel3.add(panel2, BorderLayout.EAST);

    m_ButtonClear.setMnemonic('C');
    m_ButtonClear.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e) {
	  clear();
	}
      });
    panel2.add(m_ButtonClear);

    m_ButtonClose.setMnemonic('l');
    m_ButtonClose.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e) {
	  close();
	}
      });
    panel2.add(m_ButtonClose);

    // size + current size + wordwrap
    panel2 = new JPanel(new GridLayout(1, 3));
    panel3.add(panel2, BorderLayout.WEST);
    
    // size
    panel4 = new JPanel(new FlowLayout());
    panel2.add(panel4);
    model = (SpinnerNumberModel) m_SpinnerMaxSize.getModel();
    model.setMinimum(new Integer(1));
    model.setStepSize(new Integer(1000));
    model.setValue(new Integer(100000));
    model.addChangeListener(this);

    label = new JLabel("max. Size");
    label.setDisplayedMnemonic('m');
    label.setLabelFor(m_SpinnerMaxSize);

    panel4.add(label);
    panel4.add(m_SpinnerMaxSize);

    // current size
    panel4 = new JPanel(new FlowLayout());
    panel2.add(panel4);
    panel4.add(m_LabelCurrentSize);

    // wordwrap
    panel4 = new JPanel(new FlowLayout());
    panel2.add(panel4);
    m_CheckBoxWordwrap.setSelected(m_UseWordwrap);
    m_CheckBoxWordwrap.addItemListener(new ItemListener() {
	public void itemStateChanged(ItemEvent e) {
	  setWordwrap(m_CheckBoxWordwrap.isSelected());
	}
      });
    panel4.add(m_CheckBoxWordwrap);
  }

  /**
   * clears the output
   */
  public void clear() {
    m_Output.setText("");
  }

  /**
   * closes the frame
   */
  public void close() {
    setVisible(false);
  }

  /**
   * trims the JTextPane, if too big
   */
  public void trim() {
    StyledDocument      doc;
    int                 size;
    int                 maxSize;
    int                 pos;
    
    doc = m_Output.getStyledDocument();

    // too large?
    size    = doc.getLength();
    maxSize = ((Integer) m_SpinnerMaxSize.getValue()).intValue();
    if (size > maxSize) {
      try {
        // determine EOL after which to cut
        pos = size - maxSize;
        while (!doc.getText(pos, 1).equals("\n"))
          pos++;
        while (doc.getText(pos, 1).equals("\n")) 
          pos++;
        // delete text
        doc.remove(0, pos);
      }
      catch (Exception ex) {
        // don't print it, otherwise we get an endless loop!
        if (DEBUG)
          System.out.println(ex);
      }
    }
    
    // move cursor to end
    m_Output.setCaretPosition(doc.getLength());
  }

  /**
   * returns a string representation (#RGB) of the given color
   */
  protected String colorToString(Color c) {
    String      result;
    
    result = "#" + Utils.padLeft(Integer.toHexString(c.getRed()),   2)
                 + Utils.padLeft(Integer.toHexString(c.getGreen()), 2)
                 + Utils.padLeft(Integer.toHexString(c.getBlue()),  2);

    result = result.replaceAll("\\ ", "0").toUpperCase();
    
    return result;
  }

  /**
   * toggles the wordwrap<br/>
   * override wordwrap from: 
   * http://forum.java.sun.com/thread.jspa?threadID=498535&messageID=2356174
   */
  public void setWordwrap(boolean wrap) {
    Container   parent;
    JTextPane   outputOld;
    
    m_UseWordwrap = wrap;
    if (m_CheckBoxWordwrap.isSelected() != m_UseWordwrap)
      m_CheckBoxWordwrap.setSelected(m_UseWordwrap);

    // create new JTextPane
    parent    = m_Output.getParent();
    outputOld = m_Output;
    if (m_UseWordwrap)
      m_Output = new JTextPane();
    else
      m_Output = new JTextPane(){
        private static final long serialVersionUID = -8275856175921425981L;
        public void setSize(Dimension d) {    
          if (d.width < getGraphicsConfiguration().getBounds().width) 
            d.width = getGraphicsConfiguration().getBounds().width; 
          super.setSize(d);
        }

        public boolean getScrollableTracksViewportWidth() { 
          return false; 
        }
      };
    m_Output.setEditable(false);
    m_Output.addCaretListener(this);
    m_Output.setDocument(outputOld.getDocument());
    m_Output.setCaretPosition(m_Output.getDocument().getLength());
    //m_Output.setToolTipText(
    //      "stdout = " + colorToString(COLOR_STDOUT) + ", "
    //    + "stderr = " + colorToString(COLOR_STDERR));
    parent.add(m_Output);
    parent.remove(outputOld);
  }

  /**
   * Called when the caret position is updated.
   */
  public void caretUpdate(CaretEvent e) {
    m_LabelCurrentSize.setText(
        "currently: " + m_Output.getStyledDocument().getLength());

    if (DEBUG)
      System.out.println(e);
  }

  /**
   * Invoked when the target of the listener has changed its state.
   */
  public void stateChanged(ChangeEvent e) {
    // check max size if Spinner is changed
    if (e.getSource() == m_SpinnerMaxSize.getModel()) {
      trim();
      validate();
      caretUpdate(null);
    }
  }

  /**
   * for testing only
   */
  public static void main(String[] args) {
    LogWindow       log;

    LookAndFeel.setLookAndFeel();
    
    log = new LogWindow();
    log.setVisible(true);
    log.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

    // test output
    System.out.print("a");
    System.err.print("a");
    System.out.print("a");
    System.out.println();
    System.err.println(new java.util.Date());
  }
}
