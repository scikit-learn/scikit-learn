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

package weka.gui.beans;

/**
 * Listener interface that customizer classes that are interested
 * in data format changes can implement.
 *
 * @author <a href="mailto:mhall@cs.waikato.ac.nz">Mark Hall</a>
 * @version $Revision: 1.3 $
 */
public interface DataFormatListener {
  
  /**
   * Recieve a DataSetEvent that encapsulates a new data format. The
   * DataSetEvent may contain null for the encapsulated format. This indicates
   * that there is no data format available (ie. user may have disconnected
   * an input source of data in the KnowledgeFlow).
   *
   * @param e a <code>DataSetEvent</code> value
   */
  void newDataFormat(DataSetEvent e);
}
