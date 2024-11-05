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
 *    KnowledgeFlow.java
 *    Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

/**
 * Startup class for the KnowledgeFlow. Displays a splash screen.
 *
 * @author Mark Hall
 * @version  $Revision: 4717 $
 */
public class KnowledgeFlow {

  /**
   * Static method that can be called from a running program
   * to launch the KnowledgeFlow
   */
  public static void startApp() {
    KnowledgeFlowApp.addStartupListener(new StartUpListener() {
        public void startUpComplete() {
          weka.gui.SplashWindow.disposeSplash();
        }
      });
                                        
    weka.gui.SplashWindow.splash(ClassLoader.
                                 getSystemResource("weka/gui/beans/icons/splash.jpg"));

    Thread nt = new Thread() {
        public void run() {
          weka.gui.SplashWindow.invokeMethod("weka.gui.beans.KnowledgeFlowApp", 
                                             "createSingleton", null);
        }};
      nt.start();
  }

    /**
     * Shows the splash screen, launches the application and then disposes
     * the splash screen.
     * @param args the command line arguments
     */
    public static void main(String[] args) {
      weka.core.logging.Logger.log(weka.core.logging.Logger.Level.INFO, "Logging started");
      weka.gui.SplashWindow.splash(ClassLoader.
                                   getSystemResource("weka/gui/beans/icons/splash.jpg"));
      weka.gui.SplashWindow.invokeMain("weka.gui.beans.KnowledgeFlowApp", args);
      weka.gui.SplashWindow.disposeSplash();
    }
  
}
