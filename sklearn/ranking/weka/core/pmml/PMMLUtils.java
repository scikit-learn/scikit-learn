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
 *    PMMLUtils.java
 *    Copyright (C) 2008 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core.pmml;

/**
 * Utility routines.
 * 
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com)
 * @version $Revision: 5987 $
 */
public class PMMLUtils {

  /**
   * Utility method to left or right pad strings with arbitrary characters.
   * 
   * @param source the source string
   * @param padChar the character to pad with
   * @param length the length of the resulting string
   * @param leftPad pad to the left instead of the right
   * @return a padded string
   */
  public static String pad(String source, String padChar, 
                            int length, boolean leftPad) {
    StringBuffer temp = new StringBuffer();

    if (leftPad) {
      for (int i = 0; i< length; i++) {
        temp.append(padChar);
      }
      temp.append(source);
    } else {
      temp.append(source);
      for (int i = 0; i< length; i++) {
        temp.append(padChar);
      }
    }
    return temp.toString();
  }
}
