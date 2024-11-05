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
 * GroovyPanel.java
 * Copyright (C) 2009 University of Waikato, Hamilton, New Zealand
 */

package weka.gui.scripting;

import weka.core.Utils;
import weka.gui.BrowserHelper;
import weka.gui.ComponentHelper;
import weka.gui.visualize.VisualizeUtils;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Font;
import java.awt.GridLayout;
import java.util.Properties;

import javax.swing.BorderFactory;
import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextPane;
import javax.swing.text.Document;

/**
 * A scripting panel for <a href="http://groovy.codehaus.org/" target="_blank">Groovy</a>.
 * 
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5837 $
 */
public class GroovyPanel
  extends FileScriptingPanel {
  
  /** for serialization. */
  private static final long serialVersionUID = -3475707604414854111L;
  
  /** the Groovy setup. */
  public final static String PROPERTIES_FILE = "weka/gui/scripting/Groovy.props";
  
  /**
   * Creates a new JTextPane for the code.
   * 
   * @return		the text pane
   */
  protected JTextPane newCodePane() {
    JTextPane		result;
    SyntaxDocument	doc;
    Properties		props;
    
    try {
      props = Utils.readProperties(PROPERTIES_FILE);
    }
    catch (Exception e) {
      e.printStackTrace();
      props = new Properties();
    }
    
    result = new JTextPane();
    if (props.getProperty("Syntax", "false").equals("true")) {
      doc    = new SyntaxDocument(props);
      result.setDocument(doc);
      result.setBackground(doc.getBackgroundColor());
    }
    else {
      result.setForeground(VisualizeUtils.processColour(props.getProperty("ForegroundColor", "black"), Color.BLACK));
      result.setBackground(VisualizeUtils.processColour(props.getProperty("BackgroundColor", "white"), Color.WHITE));
      result.setFont(new Font(props.getProperty("FontName", "monospaced"), Font.PLAIN, Integer.parseInt(props.getProperty("FontSize", "12"))));
    }
    
    return result;
  }
  
  /**
   * Returns an icon to be used in a frame.
   * 
   * @return		the icon
   */
  public ImageIcon getIcon() {
    return ComponentHelper.getImageIcon(IMAGES_DIR + "/groovy_small.png");
  }
  
  /**
   * Returns a panel to be displayed with the AboutAction.
   * 
   * @return		the panel with some information on the scripting panel
   */
  protected JPanel getAboutPanel() {
    JPanel	result;
    JPanel	panel;
    
    result = new JPanel(new BorderLayout());
    result.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));
    
    // image
    result.add(new JLabel(ComponentHelper.getImageIcon(IMAGES_DIR + "/groovy_medium.png")), BorderLayout.CENTER);
    
    // links
    panel = new JPanel(new GridLayout(5, 1));
    panel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));
    result.add(panel, BorderLayout.SOUTH);
    
    panel.add(new JLabel("Groovy homepage"));
    panel.add(BrowserHelper.createLink("http://groovy.codehaus.org/", null));
    panel.add(new JLabel(" "));
    panel.add(new JLabel("Weka and Groovy"));
    panel.add(BrowserHelper.createLink("http://weka.wikispaces.com/Using+Weka+from+Groovy", null));
    
    return result;
  }

  /**
   * Returns the title (without the filename).
   * 
   * @return		the plain title
   */
  public String getPlainTitle() {
    return "Groovy Console";
  }

  /**
   * Returns an initialized script object.
   * 
   * @param doc		the document to use as basis
   * @return		the initialized script
   */
  protected Script newScript(Document doc) {
    return new GroovyScript(doc);
  }
  
  /**
   * Displays the panel in a frame.
   * 
   * @param args	can take a file as first argument
   */
  public static void main(String[] args) {
    showPanel(new GroovyPanel(), args);
  }
}
