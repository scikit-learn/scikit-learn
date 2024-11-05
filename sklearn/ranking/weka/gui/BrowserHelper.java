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
 * BrowserHelper.java
 * Copyright (C) 2006 University of Waikato, Hamilton, New Zealand
 */

package weka.gui;

import java.awt.Color;
import java.awt.Component;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.lang.reflect.Method;

import javax.swing.JLabel;
import javax.swing.JOptionPane;


/**
 * A little helper class for browser related stuff. <p/>
 * 
 * The <code>openURL</code> method is based on 
 * <a href="http://www.centerkey.com/java/browser/" target="_blank">Bare Bones Browser Launch</a>,
 * which is placed in the public domain.
 * 
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 6314 $
 */
public class BrowserHelper {

  /** Linux/Unix binaries to look for */
  public final static String[] LINUX_BROWSERS = 
    {"firefox", "google-chrome", "opera", "konqueror", "epiphany", "mozilla", "netscape"};
  
  /**
   * opens the URL in a browser.
   * 
   * @param url		the URL to open
   */
  public static void openURL(String url) {
    openURL(null, url);
  }

  /**
   * opens the URL in a browser.
   * 
   * @param parent	the parent component
   * @param url		the URL to open
   */
  public static void openURL(Component parent, String url) {
    openURL(parent, url, true);
  }

  /**
   * opens the URL in a browser.
   * 
   * @param parent	the parent component
   * @param url		the URL to open
   * @param showDialog	whether to display a dialog in case of an error or
   * 			just print the error to the console
   */
  public static void openURL(Component parent, String url, boolean showDialog) {
    String osName = System.getProperty("os.name"); 
    try { 
      // Mac OS
      if (osName.startsWith("Mac OS")) { 
	Class fileMgr = Class.forName("com.apple.eio.FileManager"); 
	Method openURL = fileMgr.getDeclaredMethod("openURL", new Class[] {String.class}); 
	openURL.invoke(null, new Object[] {url}); 
      } 
      // Windows
      else if (osName.startsWith("Windows")) {
	Runtime.getRuntime().exec("rundll32 url.dll,FileProtocolHandler " + url); 
      }
      // assume Unix or Linux
      else {  
	String browser = null; 
	for (int count = 0; count < LINUX_BROWSERS.length && browser == null; count++) {
	  // look for binaries and take first that's available
	  if (Runtime.getRuntime().exec(new String[] {"which", LINUX_BROWSERS[count]}).waitFor() == 0) {
	    browser = LINUX_BROWSERS[count];
	    break;
	  }
	}
	if (browser == null) 
	  throw new Exception("Could not find web browser");
	else
	  Runtime.getRuntime().exec(new String[] {browser, url});
      }
    }
    catch (Exception e) {
      String errMsg = "Error attempting to launch web browser:\n" + e.getMessage();
      
      if (showDialog)
	JOptionPane.showMessageDialog(
	    parent, errMsg);
      else
	System.err.println(errMsg);
    }
  } 
  
  /**
   * Generates a label with a clickable link.
   * 
   * @param url		the url of the link
   * @param text	the text to display instead of URL. if null or of 
   * 			length 0 then the URL is used
   * @return		the generated label
   */
  public static JLabel createLink(String url, String text) {
    final String urlF = url;
    final JLabel result = new JLabel();
    result.setText((text == null) || (text.length() == 0) ? url : text);
    result.setToolTipText("Click to open link in browser");
    result.setForeground(Color.BLUE);
    result.addMouseListener(new MouseAdapter() {
      public void mouseClicked(MouseEvent e) {
	if (e.getButton() == MouseEvent.BUTTON1) {
	  BrowserHelper.openURL(urlF);
	}
	else {
	  super.mouseClicked(e);
	}
      }
      public void mouseEntered(MouseEvent e) {
	result.setForeground(Color.RED);
      }
      public void mouseExited(MouseEvent e) {
	result.setForeground(Color.BLUE);
      }
    });
    
    return result;
  }
}
