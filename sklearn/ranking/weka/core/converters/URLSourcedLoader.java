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
 *    URLSourcedLoader.java
 *    Copyright (C) 2006 University of Waikato, Hamilton, New Zealand
 *
 */
package weka.core.converters;

import java.net.URL;

/**
 * Interface to a loader that can load from a http url
 *
 * @author Mark Hall
 * @version $Revision 1.0 $
 */
public interface URLSourcedLoader {

  /**
   * Set the url to load from
   *
   * @param url the url to load from
   * @exception Exception if the url can't be set.
   */
  void setURL(String url) throws Exception;

  /**
   * Return the current url
   *
   * @return the current url
   */
  String retrieveURL();
}
