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
 * ComponentHelper.java
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui;

import java.awt.Component;
import java.awt.Image;
import java.net.URL;

import javax.swing.ImageIcon;
import javax.swing.JOptionPane;


/**
 * A helper class for some common tasks with Dialogs, Icons, etc.
 *
 *
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 1.3 $ 
 */

public class ComponentHelper {
  /** the default directories for images */
  public final static String[] IMAGES = {"weka/gui/", "weka/gui/images/"};
  
  /**
   * returns the ImageIcon for a given filename and directory, NULL if not successful
   * 
   * @param dir           the directory to look in for the file
   * @param filename      the file to retrieve
   * @return              the imageicon if found, otherwise null
   */
  public static ImageIcon getImageIcon(String dir, String filename) {
    URL            url;
    ImageIcon      result;
    int            i;
    
    result = null;
    url    = Loader.getURL(dir, filename);
    
    // try to find the image at default locations
    if (url == null) {
      for (i = 0; i < IMAGES.length; i++) {
        url = Loader.getURL(IMAGES[i], filename);
        if (url != null)
          break;
      }
    }
    
    if (url != null)
      result = new ImageIcon(url);
    
    return result;
  }
  
  /**
   * returns the ImageIcon for a given filename, NULL if not successful
   *
   * @param filename      the file to retrieve
   * @return              the imageicon if found, otherwise null
   */
  public static ImageIcon getImageIcon(String filename) {
    return getImageIcon("", filename);
  }
  
  /**
   * returns the Image for a given directory and filename, NULL if not successful
   * 
   * @param dir           the directory to look in for the file
   * @param filename      the file to retrieve
   * @return              the image if found, otherwise null
   */
  public static Image getImage(String dir, String filename) {
    ImageIcon      img;
    Image          result;
    
    result = null;
    img    = getImageIcon(dir, filename);
    
    if (img != null)
      result = img.getImage();
    
    return result;
  }
  
  /**
   * returns the Image for a given filename, NULL if not successful
   * 
   * @param filename      the file to retrieve
   * @return              the image if found, otherwise null
   */
  public static Image getImage(String filename) {
    ImageIcon      img;
    Image          result;
    
    result = null;
    img    = getImageIcon(filename);
    
    if (img != null)
      result = img.getImage();
    
    return result;
  }
  
  /**
   * displays a message box with the given title, message, buttons and icon
   * ant the dimension.
   * it returns the pressed button.
   * @param parent         the parent component 
   * @param title          the title of the message box
   * @param msg            the text to display
   * @param buttons        the captions of the buttons to display
   * @param messageType    the type of message like defined in <code>JOptionPane</code> (the icon is determined on this basis)
   * @return               the button that was pressed
   * @see JOptionPane
   */
  public static int showMessageBox(Component parent, String title, String msg, int buttons, int messageType) { 
    String        icon;
    
    switch (messageType) {
      case JOptionPane.ERROR_MESSAGE:
      icon = "weka/gui/images/error.gif";
      break;
    case JOptionPane.INFORMATION_MESSAGE:
      icon = "weka/gui/images/information.gif";
      break;
    case JOptionPane.WARNING_MESSAGE:
      icon = "weka/gui/images/information.gif";
      break;
    case JOptionPane.QUESTION_MESSAGE:
      icon = "weka/gui/images/question.gif";
      break;
    default:
      icon = "weka/gui/images/information.gif";
      break;
    }
    
    return JOptionPane.showConfirmDialog(parent, msg, title, buttons, messageType, getImageIcon(icon));
  }
  
  /**
   * pops up an input dialog
   * 
   * @param parent        the parent of this dialog, can be <code>null</code>
   * @param title         the title to display, can be <code>null</code>
   * @param msg           the message to display
   * @param initialValue  the initial value to display as input
   * @return              the entered value, or if cancelled <code>null</code>  
   */
  public static String showInputBox(Component parent, String title, String msg, Object initialValue) {
    Object        result;
    
    if (title == null)
      title = "Input...";
    
    result = JOptionPane.showInputDialog(
             parent, msg, title, JOptionPane.QUESTION_MESSAGE, getImageIcon("question.gif"), null, initialValue);
    
    if (result != null)
      return result.toString();
    else
      return null;
  }
}

