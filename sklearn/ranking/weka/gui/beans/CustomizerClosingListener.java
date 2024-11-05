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
 *    CustomizerClosingListener.java
 *    Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 *
 */
package weka.gui.beans;

/**
 * @author Mark Hall
 * @version $Revision: 1.4 $
 */
public interface CustomizerClosingListener {

  /**
   * Customizer classes that want to know when they are being
   * disposed of can implement this method.
   */
  void customizerClosing();

}
