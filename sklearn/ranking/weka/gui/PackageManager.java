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
 *    PackageManager.java
 *    Copyright (C) 2009 University of Waikato, Hamilton, New Zealand
 */

package weka.gui;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridLayout;
import java.awt.Image;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.DefaultCellEditor;
import javax.swing.DefaultComboBoxModel;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JEditorPane;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTable;
import javax.swing.JTextArea;
import javax.swing.JToolBar;
import javax.swing.ListSelectionModel;
import javax.swing.SwingWorker;
import javax.swing.event.HyperlinkEvent;
import javax.swing.event.HyperlinkListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.JTableHeader;
import javax.swing.table.TableColumnModel;

import org.pentaho.packageManagement.Dependency;
import org.pentaho.packageManagement.Package;
import org.pentaho.packageManagement.PackageConstraint;

import weka.core.Environment;
import weka.core.Utils;
import weka.core.WekaPackageManager;

/**
 * A GUI interface the the package management system.
 * 
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com)
 * @version $Revision: 6807 $
 */
public class PackageManager extends JPanel {
  
  /** For serialization */
  private static final long serialVersionUID = -7463821313750352385L;
  
  protected static final String PACKAGE_COLUMN = "Package";
  protected static final String CATEGORY_COLUMN = "Category";
  protected static final String INSTALLED_COLUMN = "Installed version";
  protected static final String REPOSITORY_COLUMN = "Repository version";
  protected static final String LOADED_COLUMN = "Loaded";
  
  /** The JTable for displaying the package names and version numbers */
  protected JTable m_table = new ETable();
  
  protected JSplitPane m_splitP;
    
  // protected JTextArea m_packageDescription;
  
  /** An editor pane to display package information */
  protected JEditorPane m_infoPane;
  
  /** Installed radio button */
  protected JRadioButton m_installedBut = new JRadioButton("Installed");
  
  /** Available radio button */
  protected JRadioButton m_availableBut = new JRadioButton("Available");
  
  /** All radio button */
  protected JRadioButton m_allBut = new JRadioButton("All");
  
  /** Button for installing the selected package */
  protected JButton m_installBut = new JButton("Install");
  protected JCheckBox m_forceBut = new JCheckBox("Ignore dependencies/conflicts");
  
  /** Button for uninstalling the selected package */
  protected JButton m_uninstallBut = new JButton("Uninstall");
  
  protected JButton m_refreshCacheBut = new JButton("Refresh repository cache");
  
  protected JProgressBar m_progress = new JProgressBar(0, 100);
  protected JLabel m_detailLabel = new JLabel();
  
  protected JButton m_backB;
  protected LinkedList<URL> m_browserHistory = new LinkedList<URL>();
  protected static final String BROWSER_HOME = "http://www.cs.waikato.ac.nz/ml/weka";
  protected JButton m_homeB;
  
  protected DefaultTableModel m_model;

  protected Map<String, List<Object>> m_packageLookupInfo;

  protected List<Package> m_allPackages;
  protected List<Package> m_installedPackages;
  protected List<Package> m_availablePackages;
  
  /** The column in the table to sort the entries by */
  protected int m_sortColumn = 0;
  
  /** Reverse the sort order if the user clicks the same column header twice */
  protected boolean m_reverseSort = false;
  
  protected Comparator<Package> m_packageComparator = new Comparator<Package>() {

    @Override
    public int compare(Package o1, Package o2) {
      String meta1 = "";
      String meta2 = "";
      if (m_sortColumn == 0) {
        meta1 = o1.getName();
        meta2 = o2.getName();
      } else {
        if (o1.getPackageMetaDataElement("Category") != null) {
          meta1 = o1.getPackageMetaDataElement("Category").toString();
        }
        
        if (o2.getPackageMetaDataElement("Category") != null) {
          meta2 = o2.getPackageMetaDataElement("Category").toString();
        }
      }
      
      
      int result = meta1.compareTo(meta2);
      if (m_reverseSort) {
        result = -result;
      }
      return result;
    }    
  };
  
  protected boolean m_installing = false;    
  
  class ProgressPrintStream extends PrintStream {
    
    private Progressable m_listener;
    
    public ProgressPrintStream(Progressable listener) {
      // have to invoke a super class constructor
      super(System.out);
      m_listener = listener;      
    }
    
    public void println(String string) {
      System.out.println(string); // make sure the log picks it up
      m_listener.makeProgress(string);
    }
    
    public void println(Object obj) {
      println(obj.toString());
    }
    
    public void print(String string) {
      System.out.print(string); // make sure the log picks it up
      m_listener.makeProgress(string);
    }
    
    public void print(Object obj) {
      print (obj.toString());
    }    
  }
  
  interface Progressable {
    void makeProgress(String progressMessage);
  }
  
  class EstablishCache extends SwingWorker<Void, Void> implements Progressable {
    private int m_progressCount = 0;
    private Exception m_error = null;
    
    private javax.swing.ProgressMonitor m_progress;
    
    public void makeProgress(String progressMessage) {
      m_progress.setNote(progressMessage);
      m_progressCount++;
      m_progress.setProgress(m_progressCount);
    }
    
    public Void doInBackground() {
      int numPackages = WekaPackageManager.numRepositoryPackages();
      if (numPackages < 0) {
        // there was some problem getting the file that holds this
        // information from the repository server - try to continue
        // anyway with a max value of 100 for the number of packages
        // (since all we use this for is setting the upper bound on
        // the progress bar).
        numPackages = 100;
      }
      m_progress = new javax.swing.ProgressMonitor(PackageManager.this, "Establising cache...", 
          "", 0, numPackages);
      ProgressPrintStream pps = new ProgressPrintStream(this);
      m_error = WekaPackageManager.establishCacheIfNeeded(pps);

      m_cacheEstablished = true;
      return null;
    }
    
    public void done() {
      m_progress.close();
      if (m_error != null) {
        displayErrorDialog("There was a problem establishing the package\n" +
        		"meta data cache. We'll try to use the repository" +
        		"directly.", m_error);
      }
    }
  }
  
  class RefreshCache extends SwingWorker<Void, Void> implements Progressable {
    private int m_progressCount = 0;
    private Exception m_error = null;
    
    public void makeProgress(String progressMessage) {
      m_detailLabel.setText(progressMessage);
      m_progressCount++;
      m_progress.setValue(m_progressCount);
    }
    
    public Void doInBackground() {
      m_cacheRefreshInProgress = true;
      int numPackages = WekaPackageManager.numRepositoryPackages();
      if (numPackages < 0) {
        // there was some problem getting the file that holds this
        // information from the repository server - try to continue
        // anyway with a max value of 100 for the number of packages
        // (since all we use this for is setting the upper bound on
        // the progress bar).
        numPackages = 100;
      }
      m_progress.setMaximum(numPackages);
      m_refreshCacheBut.setEnabled(false);
      m_installBut.setEnabled(false);
      m_installedBut.setEnabled(false);
      m_availableBut.setEnabled(false);
      m_allBut.setEnabled(false);      
      ProgressPrintStream pps = new ProgressPrintStream(this);
      m_error = WekaPackageManager.refreshCache(pps);
      getAllPackages();
      return null;
    }
    
    public void done() {
      m_progress.setValue(m_progress.getMinimum());
      if (m_error != null) {
        displayErrorDialog("There was a problem refreshing the package\n" +
            "meta data cache. We'll try to use the repository" +
            "directly.", m_error);
        m_detailLabel.setText("");
      } else {
        m_detailLabel.setText("Cache refresh completed");
      }
      
      m_installBut.setEnabled(true);
      m_refreshCacheBut.setEnabled(true);
      m_installedBut.setEnabled(true);
      m_availableBut.setEnabled(true);
      m_allBut.setEnabled(true);
      
      updateTable();
      m_cacheRefreshInProgress = false;
    }
  }
  
  private void pleaseCloseAppWindowsPopUp() {
    if (!Utils.getDontShowDialog("weka.gui.PackageManager.PleaseCloseApplicationWindows")) {
      JCheckBox dontShow = new JCheckBox("Do not show this message again");
      Object[] stuff = new Object[2];
      stuff[0] = "Please close any open Weka application windows\n" +
        "(Explorer, Experimenter, KnowledgeFlow, SimpleCLI)\n" + 
        "before proceeding.\n";
      stuff[1] = dontShow;

      JOptionPane.showMessageDialog(PackageManager.this, stuff, 
          "Weka Package Manager", JOptionPane.OK_OPTION);

      if (dontShow.isSelected()) {
        try {
          Utils.setDontShowDialog("weka.gui.PackageManager.PleaseCloseApplicationWindows");
        } catch (Exception ex) {
          // quietly ignore
        }
      }
    }
  }
  
  class UninstallTask extends SwingWorker<Void, Void> implements Progressable {
    
    private List<String> m_packageNamesToUninstall;
//    private String m_packageName;
  //   private boolean m_successfulUninstall = false;
    private List<String> m_unsuccessfulUninstalls = new ArrayList<String>();
    
    private int m_progressCount = 0;
    
    public void setPackages(List<String> packageNames) {
      m_packageNamesToUninstall = packageNames;
    }
    
    public void makeProgress(String progressMessage) {
      m_detailLabel.setText(progressMessage);
      m_progressCount++;
      m_progress.setValue(m_progressCount);
      if (m_progressCount == m_progress.getMaximum()) {
        m_progress.setMaximum(m_progressCount + 5);
      }
    }
    
    public Void doInBackground() {
      m_installing = true;
      m_installBut.setEnabled(false);
      m_uninstallBut.setEnabled(false);      
      m_refreshCacheBut.setEnabled(false);
      m_availableBut.setEnabled(false);
      m_allBut.setEnabled(false);
      m_installedBut.setEnabled(false);
      
      ProgressPrintStream pps = new ProgressPrintStream(this);
      m_progress.setMaximum(m_packageNamesToUninstall.size() * 5);

      for (int zz = 0; zz < m_packageNamesToUninstall.size(); zz++) {

        String packageName = m_packageNamesToUninstall.get(zz);

        boolean explorerPropertiesExist = 
          WekaPackageManager.installedPackageResourceExists(packageName, "Explorer.props");

        if (!m_forceBut.isSelected()) {
          List<Package> compromised = new ArrayList<Package>();

          // Now check to see which other installed packages depend on this one
          List<Package> installedPackages;
          try {
            installedPackages = WekaPackageManager.getInstalledPackages();
          } catch (Exception e) {
            e.printStackTrace();
            displayErrorDialog("Can't determine which packages are installed!", e);
            // return null; // can't proceed
            m_unsuccessfulUninstalls.add(packageName);
            continue;
          }
          for (Package p : installedPackages) {
            List<Dependency> tempDeps;
            try {
              tempDeps = p.getDependencies();
            } catch (Exception e) {
              e.printStackTrace();
              displayErrorDialog("Problem determining dependencies for package : " 
                  + p.getName(), e);
              //return null; // can't proceed
              m_unsuccessfulUninstalls.add(packageName);
              continue;
            }

            for (Dependency d : tempDeps) {
              if (d.getTarget().getPackage().getName().equals(packageName)) {
                // add this installed package to the list
                compromised.add(p);
                break;
              }
            }
          }

          if (compromised.size() > 0) {
            StringBuffer message = new StringBuffer();
            message.append("The following installed packages depend on " 
                + packageName + " :\n\n");
            for (Package p : compromised) {
              message.append("\t" + p.getName() + "\n");
            }

            message.append("\nDo you wish to proceed?");
            int result = JOptionPane.showConfirmDialog(PackageManager.this, message.toString(), 
                "Weka Package Manager", JOptionPane.YES_NO_OPTION);

            if (result == JOptionPane.NO_OPTION) {
              // bail out here
              //return null;
              continue;
            }          
          }
        }

//        m_progress.setMaximum(10);
        try {
          if (explorerPropertiesExist) {
            // need to remove any set Explorer properties first
            WekaPackageManager.removeExplorerProps(packageName);
          }
          WekaPackageManager.uninstallPackage(packageName, true, pps);

        } catch (Exception e) {
          e.printStackTrace();
          displayErrorDialog("Unable to uninstall package: " + packageName, e);
          //return null;
          m_unsuccessfulUninstalls.add(packageName);
          continue;
        }
      }
      
      WekaPackageManager.refreshGOEProperties();
      // m_successfulUninstall = true;
      
      return null;
    }
    
    public void done() {
      m_progress.setValue(m_progress.getMinimum());      
      if (m_unsuccessfulUninstalls.size() == 0) {
        m_detailLabel.setText("Packages removed successfully.");
        
        if (!Utils.getDontShowDialog("weka.gui.PackageManager.RestartAfterUninstall")) {
          JCheckBox dontShow = new JCheckBox("Do not show this message again");
          Object[] stuff = new Object[2];
          stuff[0] = "Weka might need to be restarted for\n" +
            "the changes to come into effect.\n";
          stuff[1] = dontShow;

          JOptionPane.showMessageDialog(PackageManager.this, stuff, 
              "Weka Package Manager", JOptionPane.OK_OPTION);

          if (dontShow.isSelected()) {
            try {
              Utils.setDontShowDialog("weka.gui.PackageManager.RestartAfterUninstall");
            } catch (Exception ex) {
              // quietly ignore
            }
          }
        }        
      } else {
        StringBuffer failedPackageNames = new StringBuffer();
        for (String p : m_unsuccessfulUninstalls) {
          failedPackageNames.append(p + "\n");
        }
        displayErrorDialog("The following package(s) could not be uninstalled\n"
            + "for some reason (check the log)\n" + failedPackageNames.toString()
            , "");
        m_detailLabel.setText("Finished uninstalling.");
      }
      
      m_refreshCacheBut.setEnabled(true);
      m_availableBut.setEnabled(true);
      m_allBut.setEnabled(true);
      m_installedBut.setEnabled(true);
      
      // force refresh of installed and available packages
      m_installedPackages = null;
      m_availablePackages = null;
//      m_installBut.setEnabled(true);
      m_installing = false;
      updateTable();
      if (m_table.getSelectedRow() >= 0) {
        // mainly to update the install/uninstall button status
        //displayPackageInfo(m_table.getSelectedRow());
        updateInstallUninstallButtonEnablement();
      }
    }
  }
   
  class InstallTask extends SwingWorker<Void, Void> implements Progressable {

    private List<String> m_packageNamesToInstall;
    private List<Object> m_versionsToInstall;
    
//    private boolean m_successfulInstall = false;
    private List<Package> m_unsuccessfulInstalls = new ArrayList<Package>();

    private int m_progressCount = 0;

    public void setPackages(List<String> packagesToInstall) {
      m_packageNamesToInstall = packagesToInstall;
    }

    public void setVersions(List<Object> versionsToInstall) {
      m_versionsToInstall = versionsToInstall;
    }

    public void makeProgress(String progressMessage) {
      m_detailLabel.setText(progressMessage);
      m_progressCount++;
      m_progress.setValue(m_progressCount);
      if (m_progressCount == m_progress.getMaximum()) {
        m_progress.setMaximum(m_progressCount + 5);
      }
    }

    /*
     * Main task. Executed in background thread.
     */
    @Override
    public Void doInBackground() {
      m_installing = true;
      m_installBut.setEnabled(false);
      m_uninstallBut.setEnabled(false);
      m_refreshCacheBut.setEnabled(false);
      m_availableBut.setEnabled(false);
      m_allBut.setEnabled(false);
      m_installedBut.setEnabled(false);
      ProgressPrintStream pps = new ProgressPrintStream(this);
      m_progress.setMaximum(m_packageNamesToInstall.size() * 30);

      for (int zz = 0; zz < m_packageNamesToInstall.size(); zz++) {
        Package packageToInstall = null;
        String packageName = m_packageNamesToInstall.get(zz);
        Object versionToInstall = m_versionsToInstall.get(zz);
        try {
          packageToInstall = WekaPackageManager.getRepositoryPackageInfo(packageName, 
             versionToInstall.toString());
        } catch (Exception e) {
          e.printStackTrace();
          displayErrorDialog("Unable to obtain package info for package: " 
              + packageName, e);
//          return null; // bail out here
          m_unsuccessfulInstalls.add(packageToInstall);
          continue;
        }

        // check for any special installation instructions
        Object specialInstallMessage = 
          packageToInstall.getPackageMetaDataElement("MessageToDisplayOnInstallation");
        if (specialInstallMessage != null && 
            specialInstallMessage.toString().length() > 0) {
          String siM = specialInstallMessage.toString();
          try {
            siM = Environment.getSystemWide().substitute(siM);
          } catch (Exception ex) {
            // quietly ignore
          }
          JOptionPane.showMessageDialog(PackageManager.this, 
              packageToInstall + "\n\n" + siM, 
              "Weka Package Manager", JOptionPane.OK_OPTION);
        }

        if (!m_forceBut.isSelected()) {
          try {
            if (!packageToInstall.isCompatibleBaseSystem()) {
              List<Dependency> baseSysDep = packageToInstall.getBaseSystemDependency();
              StringBuffer depList = new StringBuffer();
              for (Dependency bd : baseSysDep) {
                depList.append(bd.getTarget().toString() + " ");
              }

              JOptionPane.showMessageDialog(PackageManager.this, "Unable to install package " +
                  "\n" + packageName + " because it requires" +
                  "\n" + depList.toString(), 
                  "Weka Package Manager", JOptionPane.ERROR_MESSAGE);
              // bail out here
              //return null;
              m_unsuccessfulInstalls.add(packageToInstall);
              continue;
            }                    
          } catch (Exception e) {
            e.printStackTrace();
            displayErrorDialog("Problem determining dependency on base system" +
                " for package: " + packageName, e);
            //return null; // can't proceed
            m_unsuccessfulInstalls.add(packageToInstall);
            continue;
          }

          // check to see if package is already installed
          boolean upOrDowngrading = false;
          if (packageToInstall.isInstalled()) {
            Package installedVersion = null;
            try {
              installedVersion = WekaPackageManager.getInstalledPackageInfo(packageName);
            } catch (Exception e) {
              e.printStackTrace();
              displayErrorDialog("Problem obtaining package info for package: " 
                  + packageName, e);
              //return null; // can't proceed
              m_unsuccessfulInstalls.add(packageToInstall);
              continue;
            }

            if (!packageToInstall.equals(installedVersion)) {
              int result = JOptionPane.showConfirmDialog(PackageManager.this, "Package " + 
                  installedVersion + " is already installed. Replace with " + 
                  packageToInstall + "?", "Weka Package Manager", JOptionPane.YES_NO_OPTION);
              if (result == JOptionPane.NO_OPTION) {
                // bail out here
                //return null;
                m_unsuccessfulInstalls.add(packageToInstall);
                continue;
              }

              if (!Utils.getDontShowDialog("weka.gui.PackageManager.RestartAfterUpgrade")) {
                JCheckBox dontShow = new JCheckBox("Do not show this message again");
                Object[] stuff = new Object[2];
                stuff[0] = "Weka will need to be restared after installation for\n" +
                "the changes to come into effect.\n";
                stuff[1] = dontShow;

                JOptionPane.showMessageDialog(PackageManager.this, stuff, 
                    "Weka Package Manager", JOptionPane.OK_OPTION);

                if (dontShow.isSelected()) {
                  try {
                    Utils.setDontShowDialog("weka.gui.PackageManager.RestartAfterUpgrade");
                  } catch (Exception ex) {
                    // quietly ignore
                  }
                }
              }
            } else {
              int result = JOptionPane.showConfirmDialog(PackageManager.this, "Package " + 
                  installedVersion + " is already installed. Install again?",
                  "Weka Package Manager", JOptionPane.YES_NO_OPTION);
              if (result == JOptionPane.NO_OPTION) {
                // bail out here
                //return null;
                m_unsuccessfulInstalls.add(packageToInstall);
                continue;
              }
            }                    
          }


          // Now get a full list of dependencies for this package and
          // check for any conflicts
          Map<String, List<Dependency>> conflicts = new HashMap<String, List<Dependency>>();
          List<Dependency> dependencies = null;
          try {          
            dependencies = 
              WekaPackageManager.getAllDependenciesForPackage(packageToInstall, conflicts);
          } catch (Exception e) {
            e.printStackTrace();
            displayErrorDialog("Problem determinining dependencies for package: "
                + packageToInstall.getName(), e);
            //return null; // can't proceed
            m_unsuccessfulInstalls.add(packageToInstall);
            continue;
          }

          if (conflicts.size() > 0) {
            StringBuffer message = new StringBuffer();
            message.append("Package " + packageName + " requires the following packages:\n\n");
            Iterator<Dependency> depI = dependencies.iterator();
            while (depI.hasNext()) {
              Dependency d = depI.next();
              message.append("\t" + d +"\n");
            }

            message.append("\nThere are conflicting dependencies:\n\n");
            Set<String> pNames = conflicts.keySet();
            Iterator<String> pNameI = pNames.iterator();
            while (pNameI.hasNext()) {
              String pName = pNameI.next();
              message.append("Conflicts for " + pName + "\n");
              List<Dependency> confsForPackage = conflicts.get(pName);
              Iterator<Dependency> confs = confsForPackage.iterator();
              while (confs.hasNext()) {
                Dependency problem = confs.next();
                message.append("\t" + problem + "\n");
              }
            }

            JOptionPane.showConfirmDialog(PackageManager.this, message.toString(), 
                "Weka Package Manager", JOptionPane.OK_OPTION);

            // bail out here
            //return null;
            m_unsuccessfulInstalls.add(packageToInstall);
            continue;
          }

          // Next check all dependencies against what is installed and
          // inform the user about which installed packages will be altered. Also
          // build the list of only those packages that need to be installed or
          // upgraded (excluding those that are already installed and are OK).
          List<PackageConstraint> needsUpgrade = new ArrayList<PackageConstraint>();
          List<Package> finalListToInstall = new ArrayList<Package>();

          Iterator<Dependency> depI = dependencies.iterator();
          boolean depsOk = true;
          while (depI.hasNext()) {
            Dependency toCheck = depI.next();
            if (toCheck.getTarget().getPackage().isInstalled()) {
              String toCheckName = 
                toCheck.getTarget().getPackage().
                getPackageMetaDataElement("PackageName").toString();
              try {
                Package installedVersion = WekaPackageManager.getInstalledPackageInfo(toCheckName);
                if (!toCheck.getTarget().checkConstraint(installedVersion)) {
                  needsUpgrade.add(toCheck.getTarget());
                  Package mostRecent = toCheck.getTarget().getPackage();
                  if (toCheck.getTarget() instanceof 
                      org.pentaho.packageManagement.VersionPackageConstraint) {
                    mostRecent = 
                      WekaPackageManager.mostRecentVersionWithRespectToConstraint(toCheck.getTarget());
                  }
                  finalListToInstall.add(mostRecent);
                }
              } catch (Exception ex) {
                ex.printStackTrace();
                displayErrorDialog("An error has occurred while checking " +
                    "package dependencies", ex);
                // bail out here
                //return null;
                depsOk = false;
                break;
              }
            } else {
              try {
                Package mostRecent = toCheck.getTarget().getPackage();
                if (toCheck.getTarget() instanceof 
                    org.pentaho.packageManagement.VersionPackageConstraint) {
                  mostRecent = 
                    WekaPackageManager.mostRecentVersionWithRespectToConstraint(toCheck.getTarget());
                }
                finalListToInstall.add(mostRecent);
              } catch (Exception ex) {
                ex.printStackTrace();
                displayErrorDialog("An error has occurred while checking " +
                    "package dependencies", ex);
                // bail out here
                //return null;
                depsOk = false;
                break;
              }
            }
          }
          
          if (!depsOk) {
            // bail out on this package
            m_unsuccessfulInstalls.add(packageToInstall);
            continue;
          }

          if (needsUpgrade.size() > 0) {
            StringBuffer temp = new StringBuffer();
            for (PackageConstraint pc : needsUpgrade) {
              temp.append(pc + "\n");
            }
            int result = JOptionPane.showConfirmDialog(PackageManager.this, 
                "The following packages will be upgraded in order to install:\n\n" + temp.toString(),
                "Weka Package Manager", JOptionPane.YES_NO_OPTION);

            if (result == JOptionPane.NO_OPTION) {
              // bail out here
              //return null;
              m_unsuccessfulInstalls.add(packageToInstall);
              continue;
            }

            // now take a look at the other installed packages and see if
            // any would have a problem when these ones are upgraded
            boolean conflictsAfterUpgrade = false;
            List<Package> installed = null;
            try {
              installed = WekaPackageManager.getInstalledPackages();
            } catch (Exception e) {
              e.printStackTrace();
              displayErrorDialog("Unable to determine what packages are installed!", e);
              //return null; // can't proceed
              m_unsuccessfulInstalls.add(packageToInstall);
              continue;
            }
            List<Package> toUpgrade = new ArrayList<Package>();
            for (PackageConstraint pc : needsUpgrade) {
              toUpgrade.add(pc.getPackage());
            }


            // add the actual package the user is wanting to install if it
            // is going to be an up/downgrade rather than a first install since
            // other installed packages may depend on the currently installed version
            // and thus could be affected after the up/downgrade
            toUpgrade.add(packageToInstall);

            StringBuffer tempM = new StringBuffer();
            depsOk = true;
            for (int i = 0; i < installed.size(); i++) {
              Package tempP = installed.get(i);
              String tempPName = tempP.getName();
              boolean checkIt = true;
              for (int j = 0; j < needsUpgrade.size(); j++) {
                if (tempPName.equals(needsUpgrade.get(j).getPackage().getName())) {
                  checkIt = false;
                  break;
                }
              }

              if (checkIt) {
                List<Dependency> problem = null;
                try {
                  problem = tempP.getIncompatibleDependencies(toUpgrade);
                } catch (Exception e) {
                  e.printStackTrace();
                  displayErrorDialog("An error has occurred while checking " +
                      "package dependencies", e);
                  // return null; // can't continue
                  depsOk = false;
                  break;
                }
                if (problem.size() > 0) {
                  conflictsAfterUpgrade = true;

                  tempM.append("Package " + tempP.getName() + " will have a compatibility" +
                  "problem with the following packages after upgrading them:\n");
                  Iterator<Dependency> dI = problem.iterator();
                  while (dI.hasNext()) {
                    tempM.append("\t" + dI.next().getTarget().getPackage() + "\n");
                  }
                }
              }
            }
            
            if (!depsOk) {
              m_unsuccessfulInstalls.add(packageToInstall);
              continue;
            }

            if (conflictsAfterUpgrade) {
              JOptionPane.showConfirmDialog(PackageManager.this, tempM.toString() + "\n"
                  + "Unable to continue with installation.", 
                  "Weka Package Manager", JOptionPane.OK_OPTION);

             // return null; //bail out here
              m_unsuccessfulInstalls.add(packageToInstall);
              continue;
            }
          }

          if (finalListToInstall.size() > 0) {
            StringBuffer message = new StringBuffer();
            message.append("To install " + packageName + " the following packages will" +
            " be installed/upgraded:\n\n");
            for (Package p : finalListToInstall) {
              message.append("\t" + p + "\n");
            }

            int result = JOptionPane.showConfirmDialog(PackageManager.this, 
                message.toString(),
                "Weka Package Manager", JOptionPane.YES_NO_OPTION);

            if (result == JOptionPane.NO_OPTION) {
              // bail out here
              //return null;
              m_unsuccessfulInstalls.add(packageToInstall);
              continue;
            }
            m_progress.setMaximum(m_progress.getMaximum() + (finalListToInstall.size() * 30));
          }

          // OK, now we can download and install everything

          // first install the final list of dependencies
          try {
            WekaPackageManager.installPackages(finalListToInstall, pps);
          } catch (Exception e) {
            e.printStackTrace();
            displayErrorDialog("An error has occurred while installing " +
                "dependent packages", e);
            //return null;
            m_unsuccessfulInstalls.add(packageToInstall);
            continue;
          }

          // Now install the package itself
          //m_progress.setMaximum(finalListToInstall.size() * 10 + 10);
          try {
            WekaPackageManager.installPackageFromRepository(packageName, 
                versionToInstall.toString(), pps);
          } catch (Exception e) {
            e.printStackTrace();
            displayErrorDialog("Problem installing package: " + packageName,
                e);
            // return null;
            m_unsuccessfulInstalls.add(packageToInstall);
            continue;
          }
        } else {
          //m_progress.setMaximum(10);
          // just install this package without checking/downloading dependencies etc.
          try {
            WekaPackageManager.installPackageFromRepository(packageName, 
                versionToInstall.toString(), pps);
          } catch (Exception e) {
            e.printStackTrace();
            displayErrorDialog("Problem installing package: " + packageName,
                e);
            //return null;
            m_unsuccessfulInstalls.add(packageToInstall);
            continue;
          }
        }
      }

//      m_successfulInstall = true;
      
      // Make sure that the new stuff is available to all GUIs
      WekaPackageManager.refreshGOEProperties();
      return null;
    }

    public void done() {
      m_progress.setValue(m_progress.getMinimum());
      if (m_unsuccessfulInstalls.size() == 0) {
//      if (m_successfulInstall) {
        m_detailLabel.setText("Package(s) installed successfully.");        
      } else {
        StringBuffer failedPackageNames = new StringBuffer();
        for (Package p : m_unsuccessfulInstalls) {
          failedPackageNames.append(p.getName() + "\n");
        }
        displayErrorDialog("The following package(s) could not be installed\n"
            + "for some reason (check the log)\n" + failedPackageNames.toString()
            , "");
        m_detailLabel.setText("Install complete.");
      }
      
      m_refreshCacheBut.setEnabled(true);
      m_availableBut.setEnabled(true);
      m_allBut.setEnabled(true);
      m_installedBut.setEnabled(true);
      
      // force refresh of installed and available packages
      m_installedPackages = null;
      m_availablePackages = null;
      
//      m_installBut.setEnabled(true);
      m_installing = false;
      updateTable();
      if (m_table.getSelectedRow() >= 0) {
        // mainly to update the install/uninstall button status
        //displayPackageInfo(m_table.getSelectedRow());
        updateInstallUninstallButtonEnablement();
      }
    }
  }
  
  /*public class ComboBoxRenderer extends JComboBox implements TableCellRenderer {
    public ComboBoxRenderer(String[] items) {
      super(items);
    }

    public Component getTableCellRendererComponent(JTable table, Object value,
        boolean isSelected, boolean hasFocus, int row, int column) {
      if (isSelected) {
        setForeground(table.getSelectionForeground());
        super.setBackground(table.getSelectionBackground());
      } else {
        setForeground(table.getForeground());
        setBackground(table.getBackground());
      }

      // Select the current value
      setSelectedItem(value);
      return this;
    }
  }*/
  
  protected class ComboBoxEditor extends DefaultCellEditor {
    public ComboBoxEditor() {
      super(new JComboBox(new String[] {"one", "two"}));
    }
    
    public Component getTableCellEditorComponent(JTable table, Object value,
        boolean isSelected, int row, int column) {
      String packageName = m_table.getValueAt(row, 
          getColumnIndex(PACKAGE_COLUMN)).toString();
      List<Object> catAndVers = m_packageLookupInfo.get(packageName);
      List<Object> repVersions = (List<Object>)catAndVers.get(1);
      
      String[] versions = repVersions.toArray(new String[1]);
      Component combo = getComponent();
      if (combo instanceof JComboBox) {
        ((JComboBox)combo).setModel(new DefaultComboBoxModel(versions));
        ((JComboBox)combo).setSelectedItem(value);
      } else {
        System.err.println("Uh oh!!!!!");
      }
      return combo;
    }
  }
  
  protected boolean m_cacheEstablished = false;
  protected boolean m_cacheRefreshInProgress = false;
  public static String PAGE_HEADER = "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0 Transitional//EN\">\n" +
    "<html>\n<head>\n<title>Waikato Environment for Knowledge Analysis (WEKA)</title>\n" +
    "<!-- CSS Stylesheet -->\n<style>body\n{\nbackground: #ededed;\ncolor: #666666;\n" +
    "font: 14px Tahoma, Helvetica, sans-serif;;\nmargin: 5px 10px 5px 10px;\npadding: 0px;\n" +
    "}\n</style>\n\n</head>\n<body bgcolor=\"#ededed\" text=\"#666666\">\n";
  
  private static String initialPage() {
    StringBuffer initialPage = new StringBuffer();
    initialPage.append(PAGE_HEADER);
    initialPage.append("<h1>WEKA Package Manager</h1>\n\n</body></html>\n");
    return initialPage.toString();
  }
  
  protected class HomePageThread extends Thread {
    public void run() {
      try {
        m_homeB.setEnabled(false);
        m_backB.setEnabled(false);
        URLConnection conn = null;
        URL homeURL = new URL(BROWSER_HOME);
        org.pentaho.packageManagement.PackageManager pm = 
          WekaPackageManager.getUnderlyingPackageManager();
        if (pm.setProxyAuthentication()) {
          conn = homeURL.openConnection(pm.getProxy());
        } else {
          conn = homeURL.openConnection();
        }
        
        // read the html for the home page - all we want to do here is make
        // sure that the web server is responding, so that we don't tie
        // up the JEditorPane indefinitely, since there seems to be no
        // way to set a timeout in JEditorPane
        conn.setConnectTimeout(10000); // 10 seconds
        BufferedReader bi = 
          new BufferedReader(new InputStreamReader(conn.getInputStream()));
        while (bi.readLine() != null) {
          //
        }
        
        m_infoPane.setPage(BROWSER_HOME);
      } catch (Exception ex) {
        // don't make a fuss
      } finally {
        m_homeB.setEnabled(true);
        m_backB.setEnabled(true);
      }
    }
  }
  
  private int getColumnIndex(String columnName) {
    return m_table.getColumn(columnName).getModelIndex();
  }
  
  public PackageManager() {

    EstablishCache ec = new EstablishCache();
    ec.execute();


       
    while (!m_cacheEstablished) {
      try {
        Thread.sleep(1000);
      } catch (InterruptedException e1) {
        e1.printStackTrace();
      }
    }
    
    // first try and get the full list of packages
    getAllPackages();
                        
    setLayout(new BorderLayout());
    
    ButtonGroup bGroup = new ButtonGroup();
    bGroup.add(m_installedBut);
    bGroup.add(m_availableBut);
    bGroup.add(m_allBut);
    
    JPanel butPanel = new JPanel();
    butPanel.setLayout(new BorderLayout());
    
    JPanel packageDisplayP = new JPanel();
    packageDisplayP.setLayout(new BorderLayout());    
    JPanel packageDHolder = new JPanel();
    packageDHolder.setLayout(new FlowLayout());
    packageDHolder.add(m_installedBut);
    packageDHolder.add(m_availableBut);
    packageDHolder.add(m_allBut);
    packageDisplayP.add(packageDHolder, BorderLayout.SOUTH);
    packageDisplayP.add(m_refreshCacheBut, BorderLayout.NORTH);
    butPanel.add(packageDisplayP, BorderLayout.WEST);
    
    m_refreshCacheBut.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        RefreshCache r = new RefreshCache();
        r.execute();
      }
    });

    JPanel installP = new JPanel();
    JPanel buttP = new JPanel();
    buttP.setLayout(new GridLayout(1,2));
    installP.setLayout(new BorderLayout());
    buttP.add(m_installBut);
    buttP.add(m_uninstallBut);
    m_installBut.setEnabled(false);
    m_uninstallBut.setEnabled(false);
    installP.add(buttP, BorderLayout.NORTH);
    installP.add(m_forceBut, BorderLayout.SOUTH);
    butPanel.add(installP, BorderLayout.EAST);
    m_installBut.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        int[] selectedRows = m_table.getSelectedRows();
        
        if (selectedRows.length > 0) {

          //int selected = m_table.getSelectedRow();
          //if (selected != -1) {
          List<String> packageNames = new ArrayList<String>();
          List<Object> versions = new ArrayList<Object>();
          StringBuffer confirmList = new StringBuffer();
          for (int i = 0; i < selectedRows.length; i++) {
            String packageName = m_table.getValueAt(selectedRows[i], 
                getColumnIndex(PACKAGE_COLUMN)).toString();
            packageNames.add(packageName);
            Object packageVersion = m_table.getValueAt(selectedRows[i], 
                getColumnIndex(REPOSITORY_COLUMN));
            versions.add(packageVersion);
            confirmList.append(packageName + " " + packageVersion.toString() 
                + "\n");
          }
          
          JTextArea jt = new JTextArea("The following packages will be " +
          		"installed/upgraded:\n\n" + confirmList.toString(), 10, 40);
          int result = JOptionPane.showConfirmDialog(PackageManager.this, 
              new JScrollPane(jt), 
              "Weka Package Manager", JOptionPane.YES_NO_OPTION);
          
          if (result == JOptionPane.YES_OPTION) {
            pleaseCloseAppWindowsPopUp();

            InstallTask task = new InstallTask();
            task.setPackages(packageNames);
            task.setVersions(versions);
            task.execute();
          }
        }
      }
    });
    
    m_uninstallBut.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        // int selected = m_table.getSelectedRow();
        
        int[] selectedRows = m_table.getSelectedRows();
        
        if (selectedRows.length > 0) {
          List<String> packageNames = new ArrayList<String>();
          StringBuffer confirmList = new StringBuffer();
          
          for (int i = 0; i < selectedRows.length; i++) {
            String packageName = m_table.getValueAt(selectedRows[i], 
                getColumnIndex(PACKAGE_COLUMN)).toString();
            Package p = null;
            try {
              p = WekaPackageManager.getRepositoryPackageInfo(packageName);
            } catch (Exception e1) {         
//              e1.printStackTrace();
  //            continue;
              // see if we can get installed package info
              try {
              p = WekaPackageManager.getInstalledPackageInfo(packageName);
              } catch (Exception e2) {
                e2.printStackTrace();
                continue;
              }
            }
            
            if (p.isInstalled()) {
              packageNames.add(packageName);
              confirmList.append(packageName + "\n");
            }
          }
          
          if (packageNames.size() > 0) {
            JTextArea jt = new JTextArea("The following packages will be " +
            		"uninstalled:\n" + confirmList.toString(), 10, 40);
            int result = JOptionPane.showConfirmDialog(PackageManager.this, 
                 new JScrollPane(jt), 
                "Weka Package Manager", JOptionPane.YES_NO_OPTION);

            if (result == JOptionPane.YES_OPTION) {
              pleaseCloseAppWindowsPopUp();
              UninstallTask task = new UninstallTask();
              task.setPackages(packageNames);
              task.execute();
            }
          }
        }
        
/*        if (selected != -1) {
          String packageName = m_table.getValueAt(selected, 
              getColumnIndex(PACKAGE_COLUMN)).toString();
          
          pleaseCloseAppWindowsPopUp();
          UninstallTask task = new UninstallTask();
          task.setPackage(packageName);
          task.execute();
        } */
      }
    });
    
    JPanel progressP = new JPanel();
    progressP.setLayout(new BorderLayout());
    progressP.add(m_progress, BorderLayout.NORTH);
    progressP.add(m_detailLabel, BorderLayout.SOUTH);
    butPanel.add(progressP, BorderLayout.CENTER);
    
    JPanel topPanel = new JPanel();
    topPanel.setLayout(new BorderLayout());
    topPanel.setBorder(BorderFactory.createTitledBorder("Packages"));
    topPanel.add(butPanel, BorderLayout.NORTH);
    m_allBut.setSelected(true);
    
    m_allBut.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        m_table.clearSelection();
        updateTable();
        updateInstallUninstallButtonEnablement();
      }
    });
    
    m_availableBut.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        m_table.clearSelection();
        updateTable();
        updateInstallUninstallButtonEnablement();
      }
    });
    
    m_installedBut.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        m_table.clearSelection();
        updateTable();
        updateInstallUninstallButtonEnablement();
      }
    });
    
    m_model = 
      new DefaultTableModel(new String[] {PACKAGE_COLUMN, CATEGORY_COLUMN, 
          INSTALLED_COLUMN, REPOSITORY_COLUMN, LOADED_COLUMN}, 15) {
      
        private static final long serialVersionUID = -2886328542412471039L;

      public boolean isCellEditable(int row, int col) {
        if (col != 3) {
          return false;
        } else {
          return true;
        }
      }
    }; 
    
    m_table.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
    m_table.setColumnSelectionAllowed(false);
    m_table.setPreferredScrollableViewportSize(new Dimension(550, 200));
    m_table.setModel(m_model);
    if (System.getProperty("os.name").contains("Mac")) {
      m_table.setShowVerticalLines(true);
    } else {
      m_table.setShowVerticalLines(false);
    }
    m_table.setShowHorizontalLines(false);    
    m_table.getColumn("Repository version").setCellEditor(new ComboBoxEditor());
    m_table.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
      public void valueChanged(ListSelectionEvent e) {
        if (!e.getValueIsAdjusting() && !m_cacheRefreshInProgress) {         
          ListSelectionModel lm = (ListSelectionModel) e.getSource();
          boolean infoDisplayed = false;
          for (int i = e.getFirstIndex(); i <= e.getLastIndex(); i++) {
            if (lm.isSelectedIndex(i)) {
              if (!infoDisplayed) {
                // display package info for the first one in the list
                displayPackageInfo(i);
                infoDisplayed = true;
                break;
              }              
            }
          }
          updateInstallUninstallButtonEnablement();
        }
      }
    });
    
    JTableHeader header = m_table.getTableHeader();
    header.addMouseListener(new MouseAdapter() {
      public void mouseClicked(MouseEvent evt) {
        TableColumnModel colModel = m_table.getColumnModel();

        // The index of the column whose header was clicked
        int vColIndex = colModel.getColumnIndexAtX(evt.getX());
        
        // Return if not clicked on any column header or
        // clicked on the version number cols
        if (vColIndex == -1 || vColIndex > 1) {
            return;
        }
        
        if (vColIndex == m_sortColumn) {
          // toggle the sort order
          m_reverseSort = !m_reverseSort;
        } else {
          m_reverseSort = false;
        }
        m_sortColumn = vColIndex;
        updateTable();        
      }
    });
    
    topPanel.add(new JScrollPane(m_table), BorderLayout.CENTER);
        
//    add(topPanel, BorderLayout.NORTH);
    
/*    m_packageDescription = new JTextArea(10,10);
    m_packageDescription.setLineWrap(true); */
    
    try {
      //m_infoPane = new JEditorPane(BROWSER_HOME);
      String initialPage = initialPage();
      m_infoPane = new JEditorPane("text/html", initialPage);
    } catch (Exception ex) {
      m_infoPane = new JEditorPane();
    }
        
    m_infoPane.setEditable(false);
    m_infoPane.addHyperlinkListener(new HyperlinkListener() {
      public void hyperlinkUpdate(HyperlinkEvent event) {
        if (event.getEventType() == HyperlinkEvent.EventType.ACTIVATED) {
          try {
            if (event.getURL().toExternalForm().endsWith(".zip") ||
                event.getURL().toExternalForm().endsWith(".jar")) {
              // don't render archives!
            } else {
              if (m_browserHistory.size() == 0) {
                m_backB.setEnabled(true);
              }
              m_browserHistory.add(m_infoPane.getPage());
              m_infoPane.setPage(event.getURL());
            }
          } catch(IOException ioe) {
            
          }
        }
      }
    });
    
    //JScrollPane sp = new JScrollPane(m_packageDescription);
    //JScrollPane sp = new JScrollPane(m_infoPane);
    //sp.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_AS_NEEDED);
    //sp.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_AS_NEEDED);
    JPanel browserP = new JPanel();
    browserP.setLayout(new BorderLayout());
    m_backB = new JButton(new ImageIcon(loadImage("weka/gui/images/back.gif")));
    m_backB.setToolTipText("Back");
    m_backB.setEnabled(false);
    m_backB.setBorder(BorderFactory.createEmptyBorder(0, 4, 0, 4));
    m_homeB = new JButton(new ImageIcon(loadImage("weka/gui/images/home.gif")));
    m_homeB.setBorder(BorderFactory.createEmptyBorder(0, 4, 0, 4));
    m_homeB.setToolTipText("Home");
    JToolBar browserTools = new JToolBar();
    browserTools.add(m_backB);
    browserTools.add(m_homeB);
    browserTools.setFloatable(false);
    
    // Start loading the home page
    Thread homePageThread = new HomePageThread();
    
    homePageThread.setPriority(Thread.MIN_PRIORITY);
    homePageThread.start();
    
    m_backB.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        URL previous = m_browserHistory.removeLast();
        try {
          m_infoPane.setPage(previous);
          if (m_browserHistory.size() == 0) {
            m_backB.setEnabled(false);
          }
        } catch (IOException ex) {
          //
        }
      }
    });
    
    m_homeB.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        try {
          URL back = m_infoPane.getPage();
          if (back != null) {
            m_browserHistory.add(back);
          }
          
          String initialPage = initialPage();
          m_infoPane.setContentType("text/html");
          m_infoPane.setText(initialPage);          
          HomePageThread hp = new HomePageThread();
          hp.setPriority(Thread.MIN_PRIORITY);
          hp.start();
        } catch (Exception ex) {
          // don't make a fuss
        }
      }
    });
    
    browserP.add(browserTools, BorderLayout.NORTH);
    browserP.add(new JScrollPane(m_infoPane), BorderLayout.CENTER);
  //  add(browserP, BorderLayout.CENTER);
    
    m_splitP = new JSplitPane(JSplitPane.VERTICAL_SPLIT, topPanel, browserP);
    m_splitP.setOneTouchExpandable(true);

    add(m_splitP, BorderLayout.CENTER);

    updateTable();
  }
  
  private void updateInstallUninstallButtonEnablement() {
    boolean enableInstall = false;
    boolean enableUninstall = false;
    
    if (!m_installing) {
      int[] selectedRows = m_table.getSelectedRows(); 
      // check the package to see whether we should enable the
      // install button or uninstall button. Once we've determined
      // that the list contains at least one package to be installed
      // and uninstalled we don't have to check any further

      for (int i = 0; i < selectedRows.length; i++) {
        if (!enableInstall || !enableUninstall) {
          enableInstall = true; // we should always be able to install an already installed package
          String packageName = m_table.getValueAt(selectedRows[i], 
              getColumnIndex(PACKAGE_COLUMN)).toString();
          try {
            Package p = WekaPackageManager.getRepositoryPackageInfo(packageName);
            if (!enableUninstall) {
              enableUninstall = p.isInstalled(); 
            }

            /*if (!enableInstall) {
              enableInstall = !p.isInstalled();
            } */
          } catch (Exception e1) {
            // not a repository package - just enable the uninstall button
            enableUninstall = true;
            enableInstall = false;
          }
        }
      }
    }
    
    // now set the button enablement
    m_installBut.setEnabled(enableInstall);
    m_uninstallBut.setEnabled(enableUninstall);
  }
  
  private Image loadImage(String path) {
    Image pic = null;
    URL imageURL = this.getClass().getClassLoader().getResource(path);
    if (imageURL == null) {
      // ignore
    } else {
      pic = Toolkit.getDefaultToolkit().getImage(imageURL);
    }
    
    return pic;
  }
  
  /*private String getRepVersions(String packageName, String packageVersion) {
    List<Object> repVersions = m_packageVersionsLookup.get(packageName);
    StringBuffer repString = new StringBuffer();
    
    if (repVersions.size() > 1) {
      repString.append("(");
      for (int i = 0; i < repVersions.size(); i++) {
        if (!repVersions.get(i).equals(packageVersion)) {
          repString.append(repVersions.get(i).toString());
          if (i < repVersions.size() - 1) {
            repString.append(", ");
          }
        }
      }
      repString.append(")");
    } else {
      return "";
    }
      
    return repString.toString();
  } */ 
  
  private void updateTable() {
    
    if (m_installedPackages == null || m_availablePackages == null) {
      // update the loaded status
      for (Package p : m_allPackages) {
        String loadStatus = "";
        if (p.isInstalled()) {
          File packageRoot = new File(WekaPackageManager.getPackageHome().toString()
              + File.separator + p.getName());
          boolean loaded = WekaPackageManager.loadCheck(p, packageRoot);
          loadStatus = (loaded) ? "Yes" : "No - check log";
        }
        List<Object> catAndVers = m_packageLookupInfo.get(p.getName());
        catAndVers.set(2, loadStatus);
      }      
    }
    
    if (m_allBut.isSelected()) {
      m_model.setRowCount(m_allPackages.size());
      
      Collections.sort(m_allPackages, m_packageComparator);
      int row = 0;
      for (Package p : m_allPackages) {
        m_model.setValueAt(p.getName(), row, getColumnIndex(PACKAGE_COLUMN));
        
        String category = "";
        if (p.getPackageMetaDataElement("Category") != null) {
          category = (String)p.getPackageMetaDataElement("Category");
        }
        m_model.setValueAt(category, row, getColumnIndex(CATEGORY_COLUMN));
        
        String installedV = "";
        Object repositoryV = p.getPackageMetaDataElement("Version");
        // String repString = getRepVersions(p.getName(), repositoryV);
  //      String[] repVersions = getRepVersions2(p.getName(), repositoryV);
//        repositoryV = repositoryV + " " + repString;
        
        if (p.isInstalled()) {        
          try {
            Package installed = WekaPackageManager.getInstalledPackageInfo(p.getName());
            installedV = installed.getPackageMetaDataElement("Version").toString();            
          } catch (Exception ex) {
            ex.printStackTrace();
            displayErrorDialog("An error has occurred while trying to obtain" +
            		" installed package info", ex);
          }
        }
        m_model.setValueAt(installedV, row, getColumnIndex(INSTALLED_COLUMN));
        m_model.setValueAt(repositoryV, row, getColumnIndex(REPOSITORY_COLUMN));
        List<Object> catAndVers = m_packageLookupInfo.get(p.getName());
        String loadStatus = (String)catAndVers.get(2);
        m_model.setValueAt(loadStatus, row, getColumnIndex(LOADED_COLUMN));
        row++;
      }
      
      m_table.revalidate();
      m_table.repaint();
    } else if (m_installedBut.isSelected()) {
      try {
        if (m_installedPackages == null) {
          m_installedPackages = WekaPackageManager.getInstalledPackages();
        }
        
        m_model.setRowCount(m_installedPackages.size());

        int row = 0;
        for (Package p : m_installedPackages) {
          m_model.setValueAt(p.getName(), row, getColumnIndex(PACKAGE_COLUMN));

          String installedV = p.getPackageMetaDataElement("Version").toString();
          String category = "";
          if (p.getPackageMetaDataElement("Category") != null) {
            category = p.getPackageMetaDataElement("Category").toString();
          }                    
          
          List<Object> catAndVers = m_packageLookupInfo.get(p.getName());
          Object repositoryV = "-----";
          if (catAndVers != null) {
            // handle non-repository packages
            List<Object> repVersions = (List<Object>) catAndVers.get(1);
            repositoryV = repVersions.get(0);
          }
//          String repString = getRepVersions(p.getName(), repositoryV);
//          repositoryV = repositoryV + " " + repString;          

          m_model.setValueAt(category, row, getColumnIndex(CATEGORY_COLUMN));
          m_model.setValueAt(installedV, row, getColumnIndex(INSTALLED_COLUMN));
          m_model.setValueAt(repositoryV, row, getColumnIndex(REPOSITORY_COLUMN));
          if (catAndVers != null) {
            String loadStatus = (String)catAndVers.get(2);
            m_model.setValueAt(loadStatus, row, getColumnIndex(LOADED_COLUMN));
          } else {
            // handle non-repository packages
            File packageRoot = new File(WekaPackageManager.getPackageHome().toString()
                + File.separator + p.getName());
            boolean loaded = WekaPackageManager.loadCheck(p, packageRoot);
            String loadStatus = (loaded) ? "Yes" : "No - check log";
            m_model.setValueAt(loadStatus, row, getColumnIndex(LOADED_COLUMN));
          }
          row++;
        }

      } catch (Exception ex) {
        ex.printStackTrace();
      }
    } else {
      try {
        if (m_availablePackages == null) {
          m_availablePackages = WekaPackageManager.getAvailablePackages();
        }
        
        m_model.setRowCount(m_availablePackages.size());
        
        int row = 0;
        for (Package p : m_availablePackages) {
          m_model.setValueAt(p.getName(), row, getColumnIndex(PACKAGE_COLUMN));          
          String category = "";
          if (p.getPackageMetaDataElement("Category") != null) {
            category = p.getPackageMetaDataElement("Category").toString();
          }
          
          String installedV = "";
          List<Object> catAndVers = m_packageLookupInfo.get(p.getName());
          List<Object> repVersions = (List<Object>) catAndVers.get(1);          
          Object repositoryV = repVersions.get(0);
//          String repString = getRepVersions(p.getName(), repositoryV);
  //        repositoryV = repositoryV + " " + repString;

          m_model.setValueAt(category, row, getColumnIndex(CATEGORY_COLUMN));
          m_model.setValueAt(installedV, row, getColumnIndex(INSTALLED_COLUMN));
          m_model.setValueAt(repositoryV, row, getColumnIndex(REPOSITORY_COLUMN));
          String loadStatus = (String)catAndVers.get(2);
          m_model.setValueAt(loadStatus, row, getColumnIndex(LOADED_COLUMN));
          row++; 
        }
      } catch (Exception ex) {
        ex.printStackTrace();
      }
    }
  }
  
  private void displayPackageInfo(int i) {
    String packageName = m_table.getValueAt(i, 
        getColumnIndex(PACKAGE_COLUMN)).toString();

    boolean repositoryPackage = true;
    try {
      Package repP = WekaPackageManager.getRepositoryPackageInfo(packageName);
    } catch (Exception ex) {
      repositoryPackage = false;
    }
    String versionURL = WekaPackageManager.getPackageRepositoryURL().toString() 
      + "/" + packageName + "/index.html";
    
    try {
      URL back = m_infoPane.getPage();
      if (m_browserHistory.size() == 0 && back != null) {
        m_backB.setEnabled(true);
      }
      if (back != null) {
        m_browserHistory.add(back);
      }
      
      if (repositoryPackage) {
        m_infoPane.setPage(new URL(versionURL));
      } else {
        // try and display something on this non-official package
        try {
          Package p = WekaPackageManager.getInstalledPackageInfo(packageName);
          Map<?, ?> meta = p.getPackageMetaData();
          Set<?> keys = meta.keySet();
          StringBuffer sb = new StringBuffer();
          sb.append(weka.core.RepositoryIndexGenerator.HEADER);
          sb.append("<H1>" + packageName + " (Unofficial) </H1>");
          for (Object k : keys) {
            if (!k.toString().equals("PackageName")) {
              Object value = meta.get(k);
              sb.append(k + " : " + value + "<p>");
            }
          }
          sb.append("</html>\n");
          m_infoPane.setText(sb.toString());
        } catch (Exception e) {
          // ignore
        }
      }
    } catch (Exception ex) {
      ex.printStackTrace();
    }
    
    updateInstallUninstallButtonEnablement();
    if (m_availableBut.isSelected()) {
      m_uninstallBut.setEnabled(false);
    }
    
/*    if (m_installing) {
      m_installBut.setEnabled(false);
      m_uninstallBut.setEnabled(false);
    } else {
      m_installBut.setEnabled(true);
      if (m_availableBut.isSelected()) {
        m_uninstallBut.setEnabled(false);
      } else {
        try {
          Package p = WekaPackageManager.getRepositoryPackageInfo(packageName);
          m_uninstallBut.setEnabled(p.isInstalled());
        } catch (Exception ex) {
          m_uninstallBut.setEnabled(false);
        }
      }
    } */
  }
  
  private void getPackagesAndEstablishLookup() throws Exception {
    m_allPackages = WekaPackageManager.getAllPackages();

    // now fill the lookup map
    m_packageLookupInfo = new TreeMap<String, List<Object>>();
    //Iterator<Package> i = allP.iterator();
    
    for (Package p : m_allPackages) {
      // Package p = i.next();
      String packageName = p.getName();
      String category = "";
      if (p.getPackageMetaDataElement("Category") != null) {
        category = p.getPackageMetaDataElement("Category").toString();
      }
      
      // check the load status of this package (if installed)
      String loadStatus = "";
      if (p.isInstalled()) {
        File packageRoot = new File(WekaPackageManager.getPackageHome().toString());
        boolean loaded = WekaPackageManager.loadCheck(p, packageRoot);          
        loadStatus = (loaded) ? "Yes" : "No - check log";
      }

      List<Object> versions = 
        WekaPackageManager.getRepositoryPackageVersions(packageName);
      List<Object> catAndVers = new ArrayList<Object>();
      catAndVers.add(category); catAndVers.add(versions); catAndVers.add(loadStatus);
      m_packageLookupInfo.put(packageName, catAndVers);
    }
  }
  
  private void getAllPackages() {
    try {
      getPackagesAndEstablishLookup();
    } catch (Exception ex) {
      // warn the user that we were unable to get the list of packages
      // from the repository
      ex.printStackTrace();
      System.err.println("A problem has occurred whilst trying to get all " +
      		"package information. Trying a cache refresh...");
      WekaPackageManager.refreshCache(System.out);
      try {
        // try again
        getPackagesAndEstablishLookup();
      } catch (Exception e) {
        e.printStackTrace();
      }
    }
  }
  
  private void displayErrorDialog(String message, Exception e) {
    java.io.StringWriter sw = new java.io.StringWriter();
    e.printStackTrace(new java.io.PrintWriter(sw));
    
    String result = sw.toString();
    displayErrorDialog(message, result);
  }
  
  private void displayErrorDialog(String message, String stackTrace) {
    Object[] options = null;
    
    if (stackTrace != null && stackTrace.length() > 0) {
      options = new Object[2];
      options[0] = "OK"; options[1] = "Show error";
    } else {
      options = new Object[1];
      options[0] = "OK";
    }
    int result = JOptionPane.showOptionDialog(this,
        message,
        "Weka Package Manager",
        JOptionPane.YES_NO_OPTION,
        JOptionPane.ERROR_MESSAGE,
        null,
        options,
        options[0]);

    if (result == 1) {
      JTextArea jt = new JTextArea(stackTrace, 10, 40);
      JOptionPane.showMessageDialog(PackageManager.this, new JScrollPane(jt), 
          "Weka Package Manager", JOptionPane.OK_OPTION);
    }
  }
  
  /**
   * Setting the initial placement of the divider line on a JSplitPane
   * is problematic. Most of the time it positions itself just fine based
   * on the preferred and minimum sizes of the two things it divides. However,
   * sometimes it seems to set itself such that the top component is not visible
   * without manually setting the position. This method can be called (after
   * the containing frame is visible) to set the divider location to 40% of the
   * way down the window.
   */
  public void setInitialSplitPaneDividerLocation() {
    m_splitP.setDividerLocation(0.4);
  }
  
  public static void main(String[] args) {
    weka.core.logging.Logger.log(weka.core.logging.Logger.Level.INFO, "Logging started");
    LookAndFeel.setLookAndFeel();
    
    PackageManager pm = new PackageManager();    
    
    final javax.swing.JFrame jf =
      new javax.swing.JFrame("Weka Package Manager");
    jf.getContentPane().setLayout(new BorderLayout());
    jf.getContentPane().add(pm, BorderLayout.CENTER);
    jf.addWindowListener(new java.awt.event.WindowAdapter() {
      public void windowClosing(java.awt.event.WindowEvent e) {
        jf.dispose();
        System.exit(0);
      }
    });
    Dimension screenSize = jf.getToolkit().getScreenSize();
    int width = screenSize.width * 8 / 10;
    int height = screenSize.height * 8 / 10;
    jf.setBounds(width/8, height/8, width, height);
    jf.setVisible(true);
    pm.setInitialSplitPaneDividerLocation();
  }
  
}
