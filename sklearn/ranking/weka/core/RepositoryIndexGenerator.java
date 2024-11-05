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
 *    RepositoryIndexGenerator.java
 *    Copyright (C) 2000-2010 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Iterator;
import java.util.Map;
import java.util.Properties;
import java.util.Set;
import java.util.TreeMap;

/**
 * Class for generating html index files and supporting text files
 * for a Weka package meta data repository. To Run<br><br>
 * 
 * <code>java weka.core.RepositoryIndexGenerator <path to repository></code>
 * 
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com)
 * @version $Revision: 6776 $
 */
public class RepositoryIndexGenerator {
  
  public static String HEADER = "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0 Transitional//EN\">\n" +
  		"<html>\n<head>\n<title>Waikato Environment for Knowledge Analysis (WEKA)</title>\n" +
  		"<!-- CSS Stylesheet -->\n<style>body\n{\nbackground: #ededed;\ncolor: #666666;\n" +
  		"font: 14px Tahoma, Helvetica, sans-serif;;\nmargin: 5px 10px 5px 10px;\npadding: 0px;\n" +
  		"}\n</style>\n\n</head>\n<body bgcolor=\"#ededed\" text=\"#666666\">\n";
  				
  public static String BIRD_IMAGE1 = "<img src=\"Title-Bird-Header.gif\">\n";
  public static String BIRD_IMAGE2 = "<img src=\"../Title-Bird-Header.gif\">\n";
  public static String PENTAHO_IMAGE1 = "<img src=\"pentaho_logo_rgb_sm.png\">\n\n";
  public static String PENTAHO_IMAGE2 = "<img src=\"../pentaho_logo_rgb_sm.png\">\n\n";
  
  private static int[] parseVersion(String version) {
    int major = 0;
    int minor = 0;
    int revision = 0;
    int[] majMinRev = new int[3];

    try {
      String tmpStr = version;
      tmpStr = tmpStr.replace('-', '.');
      if (tmpStr.indexOf(".") > -1) {
        major  = Integer.parseInt(tmpStr.substring(0, tmpStr.indexOf(".")));
        tmpStr = tmpStr.substring(tmpStr.indexOf(".") + 1);
        if (tmpStr.indexOf(".") > -1) {
          minor  = Integer.parseInt(tmpStr.substring(0, tmpStr.indexOf(".")));
          tmpStr = tmpStr.substring(tmpStr.indexOf(".") + 1);
          if (!tmpStr.equals("")) {
            revision = Integer.parseInt(tmpStr);
          } else {
            revision = 0;
          }
        } else {
          if (!tmpStr.equals("")) {
            minor = Integer.parseInt(tmpStr);
          } else {
            minor = 0;
          }
        }
      } else {
        if (!tmpStr.equals("")) {
          major = Integer.parseInt(tmpStr);
        } else {
          major = 0;
        }
      }
    } catch (Exception e) {
      e.printStackTrace();
      major    = -1;
      minor    = -1;
      revision = -1;
    } finally {
      majMinRev[0] = major;
      majMinRev[1] = minor;
      majMinRev[2] = revision;
    }
    
    return majMinRev;
  }
  
protected static int compare(String version1, String version2) {
    
    // parse both of the versions
    int[] majMinRev1 = parseVersion(version1);
    int[] majMinRev2 = parseVersion(version2);
    
    int result;
    
    if (majMinRev1[0] < majMinRev2[0]) {
      result = -1;
    } else if (majMinRev1[0] == majMinRev2[0]) {
      if (majMinRev1[1] < majMinRev2[1]) {
        result = -1;
      } else if (majMinRev1[1] == majMinRev2[1]) {
        if (majMinRev1[2] < majMinRev2[2]) {
          result = -1;
        } else if (majMinRev1[2] == majMinRev2[2]) {
          result = 0;
        } else {
          result = 1;
        }
      } else {
        result = 1;
      }
    } else {
      result = 1;
    }
    
    return result;
  }
  
  private static String[] processPackage(File packageDirectory) throws Exception {
    System.err.println("Processing " + packageDirectory);
    File[] contents = packageDirectory.listFiles();
    File latest = null;
    ArrayList<File> propsFiles = new ArrayList<File>();
    StringBuffer versionsTextBuffer = new StringBuffer();
    
    for (int i = 0; i < contents.length; i++) {
      if (contents[i].isFile() && contents[i].getName().endsWith(".props")) {
        propsFiles.add(contents[i]);
        if (contents[i].getName().equals("Latest.props")) {
          latest = contents[i];
        } /*else {
          String versionNumber = 
            contents[i].getName().substring(0, contents[i].getName().indexOf(".props"));
          versionsTextBuffer.append(versionNumber + "\n");
        } */
      }
    }        
        
    File[] sortedPropsFiles = propsFiles.toArray(new File[0]);
    Arrays.sort(sortedPropsFiles, new Comparator<File>() {
      public int compare (File first, File second) {
        String firstV = first.getName().substring(0, first.getName().indexOf(".props"));
        String secondV = second.getName().substring(0, second.getName().indexOf(".props"));
        if (firstV.equalsIgnoreCase("Latest")) {
          return -1;
        } else if (secondV.equalsIgnoreCase("Latest")) {
          return 1;
        }
        return -RepositoryIndexGenerator.compare(firstV, secondV);
      }
    });
    
    StringBuffer indexBuff = new StringBuffer();
    indexBuff.append(HEADER + BIRD_IMAGE2);
    indexBuff.append(PENTAHO_IMAGE2);
    Properties latestProps = new Properties();
    latestProps.load(new BufferedReader(new FileReader(latest)));
    String packageName = latestProps.getProperty("PackageName") + ": ";
    String packageTitle = latestProps.getProperty("Title");
    String packageCategory = latestProps.getProperty("Category");
    if (packageCategory == null) {
      packageCategory = "";
    }
    indexBuff.append("<h2>" + packageName + packageTitle + "</h2>\n\n");
    
    String description = latestProps.getProperty("Description");
    indexBuff.append("<p>" + description.replace("\n", "<br/>") + "</p>\n\n");
    
    indexBuff.append("Versions:<br>\n");
    for (int i = 0; i < sortedPropsFiles.length; i++) {
      if (i > 0) {
        String versionNumber = 
          sortedPropsFiles[i].getName().substring(0, sortedPropsFiles[i].getName().indexOf(".props"));
        versionsTextBuffer.append(versionNumber + "\n");
        System.err.println(versionNumber);
      }
      String name = sortedPropsFiles[i].getName();
      name = name.substring(0, name.indexOf(".props"));
      indexBuff.append("<a href=\"" + name + ".html" + "\">" + name + "</a><br>\n");

      StringBuffer version = new StringBuffer();
      version.append(HEADER + BIRD_IMAGE2);
      version.append(PENTAHO_IMAGE2);
      version.append("<table summary=\"Package " + packageName + " summary\">\n");
      Properties versionProps = new Properties();
      versionProps.load(new BufferedReader(new FileReader(sortedPropsFiles[i])));

      Set<Object> keys = versionProps.keySet();
      String[] sortedKeys = keys.toArray(new String[0]);
      Arrays.sort(sortedKeys);
      // Iterator<Object> keyI = keys.iterator();
      //while (keyI.hasNext()) {
      for (String key : sortedKeys) {
        // String key = (String)keyI.next();
        if (key.equalsIgnoreCase("PackageName") || 
            key.equalsIgnoreCase("Title") || 
            key.equalsIgnoreCase("Description") ||
            key.equalsIgnoreCase("DoNotLoadIfFileNotPresentMessage") || 
            key.equalsIgnoreCase("DoNotLoadIfClassNotPresentMessage")) {

        } else {
          version.append("<tr><td valign=top>"+ key + ":</td><td width=50></td>");
          String propValue = versionProps.getProperty(key);
          
          propValue = propValue.replace("<", "&#60;");
          propValue = propValue.replace(">", "&#62;");
          propValue = propValue.replace("@", "{[at]}");
          propValue = propValue.replace("\n", "<br/>");

          /*if (key.equals("Author") || key.equals("Maintainer")) {
            propValue = propValue.replace(".", "{[dot]}");
          } */
          
          if (key.equals("PackageURL") || key.equals("URL")) {
            propValue = "<a href=\"" + propValue +"\">" + propValue + "</a>";
          }

          version.append("<td>" + propValue + "</td></tr>\n");
        }
      }
      
      version.append("</table>\n</body>\n</html>\n");
      String versionHTMLFileName = 
        packageDirectory.getAbsolutePath() + File.separator + name + ".html";
      BufferedWriter br = new BufferedWriter(new FileWriter(versionHTMLFileName));
      br.write(version.toString());
      br.flush();
      br.close();
    }
        
    indexBuff.append("</body>\n</html>\n");
    String packageIndexName = 
      packageDirectory.getAbsolutePath() + File.separator + "index.html";
    BufferedWriter br = new BufferedWriter(new FileWriter(packageIndexName));
    br.write(indexBuff.toString());
    br.flush();
    br.close();
    
    // write the versions file to the directory
    String packageVersionsName = 
      packageDirectory.getAbsolutePath() + File.separator + "versions.txt";
    br = new BufferedWriter(new FileWriter(packageVersionsName));
    br.write(versionsTextBuffer.toString());
    br.flush();
    br.close();
    
    // return indexBuff.toString();
    String[] returnInfo = new String[2];
    returnInfo[0] = packageTitle; returnInfo[1] = packageCategory;
    return returnInfo;
  }
  
  private static void writeMainIndex(Map<String, String[]> packages, 
      File repositoryHome) throws Exception {
    StringBuffer indexBuf = new StringBuffer();
    StringBuffer packageList = new StringBuffer();
    
    indexBuf.append(HEADER + BIRD_IMAGE1);
    indexBuf.append(PENTAHO_IMAGE1);
    indexBuf.append("<h1>WEKA Packages </h1>\n\n");
    indexBuf.append(/*"<h3>Download WekaLite</h3>\n<a href=\"wekaLite.jar\">wekaLite.jar</a>" + */
    		"<p><b>IMPORTANT: make sure there are no old versions of Weka (<3.7.2) in " +
    		"your CLASSPATH before starting Weka</b>\n\n");
    
    indexBuf.append("<h3>Installation of Packages</h3>\n");
    indexBuf.append("A GUI package manager is available from the \"Tools\" menu of" +
    		" the GUIChooser<br><br><code>java -jar weka.jar</code><p>\n\n");
    
    indexBuf.append("For a command line package manager type" +
    		":<br><br<code>java weka.core.WekaPackageManager -h" +
    		"</code><br><br>\n");
    indexBuf.append("<hr/>\n");
    
    indexBuf.append("<h3>Running packaged algorithms from the command line</h3>"
        + "<code>java weka.Run [algorithm name]</code><p>"
        + "Substring matching is also supported. E.g. try:<br><br>"
        + "<code>java weka.Run Bayes</code><hr/>");

    
    Set<String> names = packages.keySet();
    indexBuf.append("<h3>Available Packages (" 
        + packages.keySet().size() +")</h3>\n\n");
    
    indexBuf.append("<table>\n");
    Iterator<String> i = names.iterator();
    while (i.hasNext()) {
      String packageName = i.next();
      String[] info  = packages.get(packageName);
      String packageTitle = info[0];
      String packageCategory = info[1];
      String href = "<a href=\"./" + packageName + "/index.html\">" + packageName + "</a>";
      
      indexBuf.append("<tr valign=\"top\">\n");
      indexBuf.append("<td>" + href +
          "</td><td width=50></td><td>" + packageCategory +
          "</td><td width=50></td><td>" + packageTitle + "</td></tr>\n");
      
      // append to the package list
      packageList.append(packageName + "\n");
    }
    
    indexBuf.append("</table>\n<hr/>\n</body></html>\n");
    String indexName = 
      repositoryHome.getAbsolutePath() + File.separator + "index.html";
    BufferedWriter br = new BufferedWriter(new FileWriter(indexName));
    br.write(indexBuf.toString());
    br.flush();
    br.close();
    
    String packageListName = 
      repositoryHome.getAbsolutePath() + File.separator + "packageList.txt";
    br = new BufferedWriter(new FileWriter(packageListName));
    br.write(packageList.toString());
    br.flush();
    br.close();
    
    String numPackagesName = 
      repositoryHome.getAbsolutePath() + File.separator + "numPackages.txt";
    br = new BufferedWriter(new FileWriter(numPackagesName));
    br.write(packages.keySet().size() + "\n");
    br.flush();
    br.close();
  }
  
  /**
   * Main method for running the RepositoryIndexGenerator
   * 
   * @param args first argument needs to be the path the the
   * repository.
   */
  public static void main(String[] args) {
    try {
      
      if (args.length < 1) {
        System.err.println("Usage:\n\n\tRepositoryIndexGenerator <path to repository>");
        System.exit(1);
      }
      
      StringBuffer mainIndex = new StringBuffer();
      File repositoryHome = new File(args[0]);
      TreeMap<String, String[]> packages = new TreeMap<String, String[]>();

      // ArrayList<File> packages = new ArrayList<File>();
      File[] contents = repositoryHome.listFiles();

      for (int i = 0; i < contents.length; i++) {
        if (contents[i].isDirectory()) {
         // packages.add(contents[i]);

          // process this package and create its index
          String[] packageInfo = processPackage(contents[i]);
          packages.put(contents[i].getName(), packageInfo);
        }
      }
      
      // write the main index file
      writeMainIndex(packages, repositoryHome);
    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }
}
