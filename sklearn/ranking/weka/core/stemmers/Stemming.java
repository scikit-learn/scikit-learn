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
 * Stemming.java
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core.stemmers;

import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.core.Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.Reader;
import java.io.Writer;
import java.util.Enumeration;
import java.util.Vector;

/**
 * A helper class for using the stemmers. Run with option '-h' to list
 * all the available options.
 *
 * @author  FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5953 $
 */
public class Stemming
  implements RevisionHandler {

  /**
   * lists all the options on the command line
   *
   * @param stemmer     the stemmer to list the parameters for
   * @return 		the option string
   */
  protected static String makeOptionsString(Stemmer stemmer) {
    Vector<Option>          options;
    Enumeration     enm;
    StringBuffer    result;

    options = new Vector<Option>();
    
    // general options
    options.add(
        new Option(
          "\tDisplays this help.",
          "h", 0, "-h"));

    options.add(
        new Option(
          "\tThe file to process.",
          "i", 1, "-i <input-file>"));

    options.add(
        new Option(
          "\tThe file to output the processed data to (default stdout).",
          "o", 1, "-o <output-file>"));

    options.add(
        new Option(
          "\tUses lowercase strings.",
          "l", 0, "-l"));

    // stemmer options?
    if (stemmer instanceof OptionHandler) {
      enm = ((OptionHandler) stemmer).listOptions();
      while (enm.hasMoreElements())
        options.add((Option)enm.nextElement());
    }

    // print options
    result = new StringBuffer();
    result.append("\nStemmer options:\n\n");
    enm = options.elements();
    while (enm.hasMoreElements()) {
      Option option = (Option) enm.nextElement();
      result.append(option.synopsis() + "\n");
      result.append(option.description() + "\n");
    }

    return result.toString();
  }
  
  /**
   * Applies the given stemmer according to the given options. '-h' lists
   * all the available options for the given stemmer.
   *
   * @param stemmer     the stemmer to use
   * @param options     the options for the stemmer
   * @throws Exception	if something goes wrong
   */
  public static void useStemmer(Stemmer stemmer, String[] options) 
    throws Exception {

    Reader          reader;
    StringBuffer    input;
    Writer          output;
    String          tmpStr;
    boolean         lowerCase;
    
    // help?
    if (Utils.getFlag('h', options)) {
      System.out.println(makeOptionsString(stemmer));
      return;
    }

    // input file
    tmpStr = Utils.getOption('i', options);
    if (tmpStr.length() == 0)
      throw new IllegalArgumentException(
          "No input file defined!" + makeOptionsString(stemmer));
    else
      reader = new BufferedReader(
                  new InputStreamReader(new FileInputStream(tmpStr)));

    input = new StringBuffer();

    // output file?
    tmpStr = Utils.getOption('o', options);
    if (tmpStr.length() == 0)
      output = new BufferedWriter(
                  new OutputStreamWriter(System.out));
    else
      output = new BufferedWriter(
                  new OutputStreamWriter(new FileOutputStream(tmpStr)));

    // lowercase?
    lowerCase = Utils.getFlag('l', options);

    // stemmer options
    if (stemmer instanceof OptionHandler)
      ((OptionHandler) stemmer).setOptions(options);
    
    // unknown options?
    try {
      Utils.checkForRemainingOptions(options);
    }
    catch (Exception e) {
      System.out.println(e.getMessage());
      System.out.println(makeOptionsString(stemmer));
      return;
    }
    
    // process file
    int character;
    while ((character = reader.read()) != -1) {
      char ch = (char) character;
      if (Character.isWhitespace((char) ch)) {
        if (input.length() > 0) {
          output.write(stemmer.stem(input.toString()));
          input = new StringBuffer();
        }
        output.write(ch);
      } 
      else {
        if (lowerCase)
          input.append(Character.toLowerCase(ch));
        else
          input.append(ch);
      }
    }
    output.flush();
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5953 $");
  }
}
