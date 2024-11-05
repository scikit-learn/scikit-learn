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
 *    SimpleCLI.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui;

import weka.gui.scripting.ScriptingPanel;

import java.awt.BorderLayout;

import javax.swing.JFrame;

/**
 * Creates a very simple command line for invoking the main method of
 * classes. System.out and System.err are redirected to an output area.
 * Features a simple command history -- use up and down arrows to move
 * through previous commmands. This gui uses only AWT (i.e. no Swing).
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5142 $
 */
public class SimpleCLI
  extends JFrame {
  
  /** for serialization. */
  static final long serialVersionUID = -50661410800566036L;
  
  /**
   * Constructor.
   *
   * @throws Exception if an error occurs
   */
  public SimpleCLI() throws Exception {
    SimpleCLIPanel	panel;

    panel = new SimpleCLIPanel();
    
    setLayout(new BorderLayout());
    setTitle(panel.getTitle());
    setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
    setIconImage(panel.getIcon().getImage());
    add(panel);
    pack();
    setSize(600, 500);
    setLocationRelativeTo(null);
    setVisible(true);
  }

  /**
   * Method to start up the simple cli.
   *
   * @param args 	Not used.
   */
  public static void main(String[] args) {
    ScriptingPanel.showPanel(new SimpleCLIPanel(), args, 600, 500);
  }
}
